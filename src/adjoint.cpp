/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   adjoint.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    General Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 ** 		Georg Reuber
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 **     The current framework is developed by Georg Reuber (JGU Mainz)
 **     
 **     If you think it is helpful, please cite the following paper:
 **     Georg S. Reuber, Anton A. Popov, Boris J.P. Kaus, (2018) Deriving scaling laws in geodynamics using adjoint gradients,
 **      Tectonophysics, Vol. 746, p. 352-363. doi:10.1016/j.tecto.2017.07.017
 **  
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

// FRAMEWORK CODE FOR LaMEM TO USE ADJOINT GRADIENT (INVERSION)
//---------------------------------------------------------------------------
// COMPUTATION OF ADJOINT INVERSION
//---------------------------------------------------------------------------
// RECIPE:
// Objective function    F(x,x(p)) = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)]     // p = parameter ; x = converged solution ; xini = comparison solution (same size as jr->gsol) ; P = Projection vector containing the proportions of solution influence
// Derivative I          dF/dx     = P*x-P*x_ini
// Adjoint operation     psi       = (J^T)^-1 * dF/dx                       // J = converged Jacobian matrix
// Derivative II         dr/dp     = [r(p+h) - r(p)]/h                      // finite difference approximation of derivative of residual r vs. parameter
// Gradients             dF/dp     = -psi^T * dr/dp
//                      
// ------------------------------------------------------------
// **** TO BE MODIFIED TO ACCOUNT FOR NOMENCLATURE CHANGES **:  
//
// EXAMPLE IN INPUT FILE:
//  # General Adjoint parameters:
//  Adjoint_mode  = AdjointGradients      # options: [None; AdjointGradients, GradientDescent; Inversion]
//  Inv_Ap        = 1
//  Inv_OFdef     = 1

//  # Parameters that are mores specific to inversion cases
//  Inv_Tao       = 1
//  # Parameters
//  <AdjointParameterStart>
//		ID  		= 1		# Phase of the parameter
//		Type 		= n   	# Type of parameter   
//		InitGuess 	= 2     # Initial guess of the parameter value
//		LowerBound  = 1		# [optional] lower bound of parameter (used with TAO)
//		UpperBound  = 3		# [optional] upper bound of parameter (used with TAO)
//	<AdjointParameterEnd>
//
//  # Define the coordinates of observation points
//	<AdjointObservationPointStart>
//		Coordinate 			= 0.5 0.5 0.5					
//		VelocityComponent 	= z
//		Value  				= -0.04248
//	<AdjointObservationPointEnd>
//
// ------------------------------------------------------------
// LINEAR SOLVER:
// You can control the behaviour of the KSP object for the adjoint with the prefix "as_" (the same way as "js_")
//
// ------------------------------------------------------------
// FULL INVERSION REMARKS:
// 1) In case you want to perform the full adjoint inversion (use = 3)  and you do not want to specify comparison points by hand (OFdef = 1) make sure that you have a
//    comparison file with a petsc vector the same size as jr->gsol and called Forward_Solution.bin. Most likely you want to run a forward simulation and then change
//    the parameters to do so just run your forward model with use = 4 which automatically saves this file. Perturb the values within the input script and solve again
//    with use = 3.
//
// 2) You can control the behaviour of the TAO object for the adjoint with the prefix "tao_" (example: '-tao_type lmvm' ; '-tao_fatol 1e-15' ; '-tao_converged_reason')
//
// 3) In case you want to use powerlaw viscosities make sure the reference strainrate is the same as the one that you use in the viscosity computation!!
//
// ------------------------------------------------------------
// IMPORTANT REMARKS:
// 1) Since the Adjoint needs the Jacobian matrix for computing the gradients it's crucial to make sure that
//    you compute the Jacobian matrix in the timesteps where you want to compute the gradients
//    (f.e. a linear problem would need a low value for -snes_atol [1e-20] and a low max iteration count -snes_max_it [2] to
//    guarantee the computation of the Jacobian + the option '-snes_type ksponly')
// 2) This code does not actually solve the system with the transposed Jacobian but uses the original Jacobian as approximation. So you should make
//    sure that your Jacobian is symmteric (e.g. use lithostatic pressure in plasticity, etc..)
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "adjoint.h"
#include "phase.h"
#include "tools.h"
#include "fdstag.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "objFunct.h"
#include "constEq.h"
#include "parsing.h"
#include "gravity.h"
#include <petscsys.h> 
//-----------------------------------------------------------------------------
// A bit stupid that this has to be twice declared, but the original function is only in AVD.cpp and not in a header file anymore ...
PetscInt FindPointInCellAdjoint(
	PetscScalar *px, // node coordinates
	PetscInt     L,  // index of the leftmost node
	PetscInt     R,  // index of the rightmost node
	PetscScalar  x)  // point coordinate
{
	// find ID of the cell containing point (call this function for local point only!)
	if(x < px[L] || x > px[R])
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Non-local marker");
	}
	// get initial guess assuming uniform grid
	PetscInt M = L + (PetscInt)((x-px[L])/((px[R]-px[L])/(PetscScalar)(R-L)));

	if(M == R) return R-1;

	if(px[M]   <= x) L=M;
	if(px[M+1] >= x) R=M+1;

	while((R-L) > 1)
	{
		M = (L+R)/2;
		if(px[M] <= x) L=M;
		if(px[M] >= x) R=M;

	}
	return(L);
}
//---------------------------------------------------------------------------
void swapStruct(struct Material_t *A, struct Material_t *B){
    struct Material_t temp = *A;
    *A = *B;
    *B = temp;
}
//---------------------------------------------------------------------------
void AddParamToList(PetscInt ID, PetscScalar value, const char par_str[_str_len_], PetscInt iP, 
		char type_name[][_str_len_],
		PetscInt 	*phsar,
		PetscScalar *Par,
		PetscInt    *FDgrad,
		PetscScalar *FDeps)
{
	PetscScalar val;
	PetscBool 	found;	
	char     	*dbkey;

	strcpy(type_name[iP], par_str);
	phsar[iP] 		=  	ID;
	
	// Check if there is a command-line option & use that instead
	asprintf(&dbkey, "-%s[%i]", par_str,ID);
	PetscOptionsGetScalar(NULL, NULL, dbkey, &val, &found);
	if (found){
		value = val;	// found a command-line option
	}


	Par[iP] 		=  	value;
	Parameter_SetFDgrad_Option(&FDgrad[iP], type_name[iP]);	
	FDeps[iP] 		=	0.0;
}

//---------------------------------------------------------------------------
/* 
	This (optionally) reads in all material parameters of the ParamFile for 
	which we will compute gradients
*/
#undef __FUNCT__
#define __FUNCT__ "Adjoint_ScanForMaterialParameters"
PetscErrorCode Adjoint_ScanForMaterialParameters(FB *fb, Scaling *scal, PetscInt *iP, 
		char type_name[][_str_len_],
		PetscInt 	*phsar,
		PetscScalar *Par,
		PetscInt    *FDgrad,
		PetscScalar *FDeps)
{
	PetscFunctionBegin;
	PetscErrorCode 	ierr;
	PetscInt 		i, jj, ID;
	PetscBool 		ReadAllMatParams=PETSC_FALSE, AddParamToGradient;
	char 			par_str[_str_len_], adjointstr[_str_len_];
	char        	ndiff[_str_len_], ndisl[_str_len_], npeir[_str_len_];
	Material_t 		m;

	
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<AdjointParameterStart>", "<AdjointParameterEnd>"); CHKERRQ(ierr);

	// error checking
	if(fb->nblocks > _MAX_PAR_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many adjoint parameters specified! Max allowed: %lld", (LLD)_MAX_PAR_);
	}
	if(!fb->nblocks)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define adjoint Parameters (mdN) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

	// Loop over all AdjointParameterBlocks to scan for the keyword 
	for(i = 0; i < fb->nblocks; i++)
	{
		ierr = getStringParam(fb, _OPTIONAL_, "Type", par_str, NULL); CHKERRQ(ierr);
		if (!strcmp(par_str,"AllMaterialParameters")){ ReadAllMatParams=PETSC_TRUE;	}
	}
	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	if (!ReadAllMatParams){ PetscFunctionReturn(0); }		// only continue if the keyword is present

	// Start reading all material parameters if requested
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| Adjoint: using all listed material parameters to compute gradients:   \n");

	// setup block access mode
	ierr = FBFindBlocks(fb, _REQUIRED_, "<MaterialStart>", "<MaterialEnd>"); CHKERRQ(ierr);

	// read each individual phase
	for(jj = 0; jj < fb->nblocks; jj++)
	{
		char     	*ptr, *line, **lines, *par_str;
		PetscInt  	i, lnbeg, lnend;
		PetscScalar value;
			
		ierr 	= getIntParam(fb, _REQUIRED_, "ID", &ID, 1, _max_num_phases_); CHKERRQ(ierr); // phase ID

		// get line buffer & pointers
		line  	= fb->lbuf;
		lines 	= FBGetLineRanges(fb, &lnbeg, &lnend);
		for(i = lnbeg; i < lnend; i++)
		{
			strcpy(line, lines[i]);					// copy line for parsing
			ptr 	= strtok(line, " ");
			par_str = ptr;							// name of the parameter
			if (par_str){							// if not empty
				if (strcmp(par_str,"ID")){
				
					// In case there are some parameters within the MaterialStructure that should NOT be taken into account, check that here
					AddParamToGradient = PETSC_TRUE;
					if 		(!(strcmp(par_str,"visID"))	)		{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"rho_ph")))		{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"rho_ph_dir")))	{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"disl_prof")))	{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"diff_prof")))	{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"peir_prof")))	{AddParamToGradient = PETSC_FALSE;	}
					else if (!(strcmp(par_str,"Name")))			{AddParamToGradient = PETSC_FALSE;	}
						
					// And some parameters, in particular creep laws, actually set a few other parameters (such as powerlaw exponent). 
					// It would be good to automatically add th
					if (AddParamToGradient){
						ptr		= 	strtok(NULL, " ");	// Space before equal sign
						ptr   	= 	strtok(NULL, " ");	// Space after equal sign
						
						// retrieve values after equal sign [NOTE: we can only deal with a single value & assume it to be a scalar]
						value 	= 	(PetscScalar)strtod(ptr, NULL);
						
						// ADD them to the database
						AddParamToList(ID, value, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps);
						
						*iP             =  *iP + 1;

					}
				}
			}
		}


		// We need to treat creep profiles in a different manner, as they set several material parameters at once
		ierr = GetProfileName(fb, scal, ndisl, "disl_prof"); CHKERRQ(ierr);
		if(strlen(ndisl)){
			ierr = SetDislProfile(&m, ndisl);  CHKERRQ(ierr);
			ierr = PetscMalloc((size_t)_str_len_*sizeof(char), par_str); CHKERRQ(ierr);

			// Set parameters 
			strcpy(par_str,"Bn"); AddParamToList(ID, m.Bn, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"n");  AddParamToList(ID, m.n,  par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"En"); AddParamToList(ID, m.En, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"Vn"); AddParamToList(ID, m.Vn, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
		}

		ierr = GetProfileName(fb, scal, ndiff, "diff_prof"); 	CHKERRQ(ierr);
		if(strlen(ndiff)){
			ierr = SetDiffProfile(&m, ndiff);  					CHKERRQ(ierr);
			ierr = PetscMalloc((size_t)_str_len_*sizeof(char), par_str); CHKERRQ(ierr);

			// Set parameters 	
			strcpy(par_str,"Bd"); AddParamToList(ID, m.Bd, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"Ed"); AddParamToList(ID, m.Ed, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"Vd"); AddParamToList(ID, m.Vd, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
		}

		ierr = GetProfileName(fb, scal, npeir, "peir_prof"); 	CHKERRQ(ierr);
		if(strlen(npeir)){
			ierr = SetPeirProfile(&m, npeir);         			CHKERRQ(ierr);
			ierr = PetscMalloc((size_t)_str_len_*sizeof(char), &par_str); CHKERRQ(ierr);

			// Set parameters 	
			strcpy(par_str,"Bp"); 		AddParamToList(ID, m.Bp,    par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"Ep"); 		AddParamToList(ID, m.Ep, 	par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"Vp"); 		AddParamToList(ID, m.Vp, 	par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"taup"); 	AddParamToList(ID, m.taup, 	par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"gamma"); 	AddParamToList(ID, m.gamma, par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
			strcpy(par_str,"q"); 		AddParamToList(ID, m.q, 	par_str, *iP, type_name, phsar, Par, FDgrad, FDeps); ++*iP;
		}

		fb->blockID++;

	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// Print overview & indicate which parameters are not specified
	for(jj = 0; jj < *iP; jj++){
		strcpy(par_str, type_name[jj]);

		if (FDgrad[jj]==0){strcpy(adjointstr, "adjoint"); }
		else{strcpy(adjointstr, "FD     "); }
			
		// Print overview & indicate which parameters are not specified
		if (FDgrad[jj]){
			PetscPrintf(PETSC_COMM_WORLD, "|   %-2i: %5s       %+6s[%-2i] = %-9.4g   \n",jj+1,adjointstr, par_str,phsar[jj],Par[jj]);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD, "|   %-2i: %5s       %+6s[%-2i] = %-9.4g   \n",jj+1,adjointstr, par_str,phsar[jj],Par[jj]);
		}
	}

	if (*iP>_MAX_PAR_){
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many inverse parameters specified! Max allowed: %lld", (LLD)  _MAX_PAR_);
	}
	


	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
/* This reads the input parameters from the file/command-line & sets default 
values. Also performs error-checking
*/
#undef __FUNCT__
#define __FUNCT__ "LaMEMAdjointReadInputSetDefaults"
PetscErrorCode LaMEMAdjointReadInputSetDefaults(ModParam *IOparam, Adjoint_Vecs *Adjoint_Vectors)
{
	PetscFunctionBegin;
	PetscErrorCode 	ierr;
	FB 				*fb;
	PetscScalar     *gradar, *Ubar, *Lbar, ts, *Par, mean, var;
	PetscInt         i, j, ti, ID, iStart, p, ct1, ct2;
	char             str[_str_len_], par_str[_str_len_], Vel_comp[_str_len_], ParType[_str_len_];
	Scaling          scal;

	fb 					=	IOparam->fb;	// filebuffer

	// Some defaults
	IOparam->FS         		= 0;
	IOparam->MfitType           = 0;
	IOparam->Gr         		= 1;
	IOparam->SCF        		= 0;
	IOparam->mdI        		= 0;
	IOparam->Ab         		= 0;
	IOparam->Ap         		= 1;
	IOparam->Adv        		= 0;
	IOparam->OFdef      		= 1;
	IOparam->Tao        		= 1;
	IOparam->tol        		= 1e-10;
	IOparam->facLS      		= 2;
	IOparam->facB       		= 0.5;
	IOparam->maxfac     		= 100;
	IOparam->Scale_Grad 		= 0.1;
	IOparam->maxit     	 		= 50;
	IOparam->maxitLS    		= 20;
	IOparam->ScalLaws   		= 0;
	IOparam->ReferenceDensity 	= 0;
	IOparam->SCF 				= 0; 		
	
    // Create scaling object
	ierr = ScalingCreate(&scal, fb, PETSC_FALSE); CHKERRQ(ierr);

	// Some general Adjoint Gradient parameters:
    ierr = getStringParam(fb, _OPTIONAL_, "Adjoint_GradientCalculation", str, NULL); CHKERRQ(ierr);  // must have component
    if     	(!strcmp(str, "CostFunction"))      IOparam->Gr=0;
	else if (!strcmp(str, "Solution"))          IOparam->Gr=1;
	else{	SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [Solution; CostFunction] as parameter for Adjoint_GradientCalculation, not %s",str);} 

	ierr = getStringParam(fb, _OPTIONAL_, "Adjoint_ScaleCostFunction", str, NULL); CHKERRQ(ierr);  // must have component
    if     	(!strcmp(str, "None"))      IOparam->SCF=0;
	else if (!strcmp(str, "Mean"))      IOparam->SCF=1;
	else if (!strcmp(str, "Var"))       IOparam->SCF=2;
	else{	SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [None; Mean; Var] as parameter for Adjoint_ScaleCostFunction, not %s",str);} 

	ierr = getIntParam   (fb, _OPTIONAL_, "Adjoint_FieldSensitivity"         , &IOparam->FS,        1, 1        ); CHKERRQ(ierr);  // Do a field sensitivity test? -> Will do the test for the first InverseParStart that is given!
	ierr = getIntParam   (fb, _OPTIONAL_, "Adjoint_ObservationPoints"        , &IOparam->Ap,        1, 3        ); CHKERRQ(ierr);  // 1 = several indices ; 2 = the whole domain ; 3 = surface
	ierr = getIntParam   (fb, _OPTIONAL_, "Adjoint_AdvectPoint"              , &IOparam->Adv,       1, 1        ); CHKERRQ(ierr);  // 1 = advect the point
	ierr = getIntParam   (fb, _OPTIONAL_, "Adjoint_ObjectiveFunctionDef"     , &IOparam->OFdef,     1, 1        ); CHKERRQ(ierr);  // Objective function defined by hand?
	ierr = getIntParam   (fb, _OPTIONAL_, "Adjoint_PrintScalingLaws"     	 , &IOparam->ScalLaws,  1, 1        ); CHKERRQ(ierr);  // Print scaling laws (combined with AdjointGradients)
	ierr = getStringParam(fb, _OPTIONAL_, "Adjoint_ScalingLawFilename"     	 , str,  "ScalingLaw.dat"  ); 		   CHKERRQ(ierr);  // Scaling law filename
	ierr = getScalarParam(fb, _OPTIONAL_, "Adjoint_DII_ref"       			 , &IOparam->DII_ref,   1, 1        ); CHKERRQ(ierr);  // Reference strainrate needed for direct FD for pointwise kernels for powerlaw viscosity (very unflexible so far)
	ierr  = PetscMemcpy(IOparam->ScalLawFilename, 	str,   (size_t)_str_len_*sizeof(char) ); 		  	 		   CHKERRQ(ierr); 
   
  	ierr = getScalarParam(fb, _OPTIONAL_, "Adjoint_ReferenceDensity"       	 , &IOparam->ReferenceDensity, 1, 1 ); CHKERRQ(ierr);  // Reference density (density parameters are computed w.r.t. this)
	
	// If we do inversion, additional parameters can be specified:
	ierr = getIntParam   (fb, _OPTIONAL_, "Inversion_maxit"     			, &IOparam->maxit,     1, 1500     ); CHKERRQ(ierr);  // maximum number of inverse iterations
	ierr = getIntParam   (fb, _OPTIONAL_, "Inversion_maxit_linesearch"   	, &IOparam->maxitLS,   1, 1500     ); CHKERRQ(ierr);  // maximum number of backtracking	
	ierr = getIntParam   (fb, _OPTIONAL_, "Inversion_ApplyBounds"        	, &IOparam->Ab,        1, 1        ); CHKERRQ(ierr);  // Apply bounds?
	ierr = getIntParam   (fb, _OPTIONAL_, "Inversion_EmployTAO"       		, &IOparam->Tao,       1, 1        ); CHKERRQ(ierr);  // Use TAO?
	ierr = getScalarParam(fb, _OPTIONAL_, "Inversion_rtol"       			, &IOparam->tol,       1, 1        ); CHKERRQ(ierr);  // tolerance for F/Fini after which code has converged
	ierr = getScalarParam(fb, _OPTIONAL_, "Inversion_factor_linesearch"     , &IOparam->facLS,     1, 1        ); CHKERRQ(ierr);  // factor in the line search that multiplies current line search parameter if GD update was successful (increases convergence speed)
	ierr = getScalarParam(fb, _OPTIONAL_, "Inversion_facB"      			, &IOparam->facB,      1, 1        ); CHKERRQ(ierr);  // backtrack factor that multiplies current line search parameter if GD update was not successful
	ierr = getScalarParam(fb, _OPTIONAL_, "Inversion_maxfac"    			, &IOparam->maxfac,    1, 1        ); CHKERRQ(ierr);  // limit on the factor (only used without tao)
	ierr = getScalarParam(fb, _OPTIONAL_, "Inversion_Scale_Grad"			, &IOparam->Scale_Grad,1, 1        ); CHKERRQ(ierr);  // Magnitude of initial parameter update (factor_ini = Scale_Grad/Grad)

	PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"|                                      LaMEM                                \n");
	PetscPrintf(PETSC_COMM_WORLD,"|                        Adjoint Gradient Framework Active                  \n");
	PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");

    
    PetscPrintf(PETSC_COMM_WORLD,"| Adjoint parameters:  \n");
	
	if(IOparam->use == _adjointgradients_ ) 
	{
		PetscPrintf(PETSC_COMM_WORLD, "|    Adjoint mode                             : AdjointGradients  \n");
		if (IOparam->Gr==0){ PetscPrintf(PETSC_COMM_WORLD, "|    Gradients are computed w.r.t.            : CostFunction \n", IOparam->Gr); }
		else               { PetscPrintf(PETSC_COMM_WORLD, "|    Gradients are computed w.r.t.            : Solution     \n", IOparam->Gr); }
		PetscPrintf(PETSC_COMM_WORLD, "|    Field-based gradient evaluation          : %d    \n", IOparam->FS);		

		if 		(IOparam->Ap == 1){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : several observation points \n"); }
		else if (IOparam->Ap == 2){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : whole domain \n"); }
		else if (IOparam->Ap == 3){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : surface      \n"); }
		
		PetscPrintf(PETSC_COMM_WORLD, "|    Advect evaluation points with flow       : %d    \n", IOparam->Adv);

		PetscPrintf(PETSC_COMM_WORLD, "|    Objective function type                  : %d    \n", IOparam->MfitType);
		
		PetscPrintf(PETSC_COMM_WORLD, "|    Objective function defined in input      : %d    \n", IOparam->OFdef);

		if ((IOparam->Gr==0) & (IOparam->ScalLaws==1) ){
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "| If you want scaling laws, you need to have Adjoint_GradientCalculation=Solution rather than CostFunction \n",IOparam->use);
		}
	}
	else if(IOparam->use == _gradientdescent_) 
	{
		PetscPrintf(PETSC_COMM_WORLD, "|    Adjoint mode                             : Gradient descent (or Quasi-Newton) inversion  \n");
		PetscPrintf(PETSC_COMM_WORLD, "|    Use Tao BLMVM (or LaMEM steepest descent): %d    \n", IOparam->Tao);
		if 		(IOparam->Ap == 1){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : several observation points  \n"); }
		else if (IOparam->Ap == 2){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : whole domain   \n"); }
		else if (IOparam->Ap == 3){PetscPrintf(PETSC_COMM_WORLD, "|    Gradient evaluation points               : surface       \n"); }
		PetscPrintf(PETSC_COMM_WORLD, "|    Advect evaluation points with flow       : %d    \n", IOparam->Adv);

		PetscPrintf(PETSC_COMM_WORLD, "|    Objective function type                  : %d    \n", IOparam->MfitType);

		PetscPrintf(PETSC_COMM_WORLD, "|    Objective function defined in input      : %d    \n", IOparam->OFdef);

		PetscPrintf(PETSC_COMM_WORLD, "|    Maximum gradient descent iterations      : %d    \n", IOparam->maxit);
		PetscPrintf(PETSC_COMM_WORLD, "|    Maximum linesearch iterations            : %d    \n", IOparam->maxitLS);
		PetscPrintf(PETSC_COMM_WORLD, "|    Apply bounds                             : %d    \n", IOparam->Ab);
		PetscPrintf(PETSC_COMM_WORLD, "|    Tolerance (F/Fini)                       : %.5e  \n", IOparam->tol);
		if (IOparam->Tao == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "|    Not employing TAO, but instead our build-in gradient algorithm, with the following parameters: \n", IOparam->facLS);
			PetscPrintf(PETSC_COMM_WORLD, "|     Linesearch factor (successful update)    : %.5e  \n", IOparam->facLS);
			PetscPrintf(PETSC_COMM_WORLD, "|     Linesearch factor (overstep)             : %.5e  \n", IOparam->facB);
			PetscPrintf(PETSC_COMM_WORLD, "|     Maximum linesearch factor                : %.5e  \n", IOparam->maxfac);
			PetscPrintf(PETSC_COMM_WORLD, "|     Scale for initial parameter update       : %.5e  \n", IOparam->Scale_Grad);
		}
	} 
	else if (IOparam->use == _syntheticforwardrun_) 
	{
		PetscPrintf(PETSC_COMM_WORLD, "|    Adjoint mode                             : SyntheticForwardRun  (saving forward run for debugging purposes)  \n");
	}
	else
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "\n| Use = %d not known; should be within [0-4]\n",IOparam->use);
	}

	ierr = AdjointVectorsCreate(Adjoint_Vectors, IOparam);      CHKERRQ(ierr);
	//ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->maxit, PETSC_DETERMINE, &IOparam->fcconv);   CHKERRQ(ierr);   // 1500 is the maximum inversion iterations that are accepted

	// TEMPORARY VARIABLES
	PetscInt		phsar[_MAX_PAR_];
	PetscInt		FDgrad[_MAX_PAR_];
	PetscScalar 	FDeps[_MAX_PAR_];
	PetscInt 		vec_log10[_MAX_PAR_];
	char 			type_name[_MAX_PAR_][_str_len_];
	
	PetscScalar     Ax[  _MAX_OBS_];
	PetscScalar     Ay[  _MAX_OBS_];
	PetscScalar     Az[  _MAX_OBS_];
	PetscInt        Av[  _MAX_OBS_];
	PetscScalar     Ae[  _MAX_OBS_];
	char 			ObsName[_MAX_OBS_][5];
	VecGetArray(Adjoint_Vectors->P,&Par);
	
    // Read material input parameters from file
    ierr = FBFindBlocks(fb, _REQUIRED_, "<MaterialStart>", "<MaterialEnd>"); CHKERRQ(ierr);
    
	// initialize ID for consistency checks
	//for(jj = 0; jj < _max_num_phases_; jj++) dbm.phases[jj].ID = -1;

	// error checking
	if(fb->nblocks > _max_num_phases_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)_max_num_phases_);
	}

	// PARAMETERS
	

	// 1) Check whether we want to take ALL material parameters into account
	iStart 	= 0;
	ierr 	= Adjoint_ScanForMaterialParameters(fb, &scal, &iStart, type_name, phsar, Par, FDgrad, FDeps); CHKERRQ(ierr);
	for(j = 0; j < iStart; j++){vec_log10[j]=0; }	// initialize (no log10 in LaMEM material parameters)

	// 2) Check the AdjointParameter blocks for additional parameters
	// Get parameter / typ / etc.
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<AdjointParameterStart>", "<AdjointParameterEnd>"); CHKERRQ(ierr);

	// error checking
	if(fb->nblocks > _MAX_PAR_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many adjoint parameters specified! Max allowed: %lld", (LLD)_MAX_PAR_);
	}
	if(!fb->nblocks)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define adjoint Parameters (mdN) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

    // temporary strings
    char        *lb_str, *ub_str, *val_str, logstr[_str_len_], adjointstr[_str_len_];
    PetscScalar par_val;
	PetscInt 	grad;


	// Read parameters which we will use in the inversion, or for which we will compute gradients
    VecGetArray(Adjoint_Vectors->Ub,&Ubar);
	VecGetArray(Adjoint_Vectors->Lb,&Lbar);
	VecGetArray(Adjoint_Vectors->grad,&gradar);
	fb->blockID = 0;
	i 			= iStart;
	//PetscPrintf(PETSC_COMM_WORLD, "| \n|    # of adjoint parameter blocks found      : %i   \n", fb->nblocks);
	for(j = 0; j < fb->nblocks; j++)
	{
		// Retrieve name of the parameter
		ierr = getStringParam(fb, _OPTIONAL_, "Type", par_str, NULL); CHKERRQ(ierr);
		
		if (strcmp(par_str,"AllMaterialParameters")){
			// If it is not the keyword mentioned above
			strcpy(type_name[i], par_str);

			ierr 		= getIntParam   (fb, _REQUIRED_, "ID" , &ID, 1, _max_num_phases_); CHKERRQ(ierr);		// phase at which it applies
			phsar[i]  	= ID;                    // PHASE

			// Parameter value
			par_val     = 0;
			ierr        = getScalarParam(fb, _OPTIONAL_, "InitGuess", &par_val, 1, 1); CHKERRQ(ierr);
			Par[i]      = par_val;                   

			// Upper bound    
			ts = 0;
			ierr = getScalarParam(fb, _OPTIONAL_, "UpperBound", &ts, 1, 1 ); CHKERRQ(ierr);
			if (ts){
				asprintf(&ub_str, "%-9.4g", ts);
				Ubar[i]   = ts;  
				IOparam->Ab = 1;
			} 
			else{
				asprintf(&ub_str, "%-9s", "-");
				Ubar[i]   = par_val;  
			}

			// Lower bound 
			ts = 0;
			ierr = getScalarParam(fb, _OPTIONAL_, "LowerBound", &ts, 1, 1 ); CHKERRQ(ierr);
			if (ts){
				asprintf(&lb_str, "%-9.4g", ts);
				Lbar[i]   	= ts;  
				IOparam->Ab = 1; // tell code that we employ bounds
			} 
			else{
				asprintf(&lb_str, "%-9s", "-");
				Lbar[i]   = par_val;  
			}
			
			// Compute gradient by finite differences (optional)
			grad = -1;
			ierr = getIntParam(fb, _OPTIONAL_, "FD_gradient", &grad, 1,1); CHKERRQ(ierr);  // must have component
			if (grad<0){
				// Retrieve scaling for the parameter
				ierr = Parameter_SetFDgrad_Option(&FDgrad[i], par_str);	CHKERRQ(ierr);
			}
			else{
				FDgrad[i] = grad;
			}

			// do we cpmpute with the log10(param) internally (implying that the parameter value )
			p 		=	0;
			ierr 	= getIntParam(fb, _OPTIONAL_, "log10", &p, 1, 1 ); CHKERRQ(ierr);	// eps for brute force FD gradients
			vec_log10[i] = p;

			ts 		= 0;
			ierr 	= getScalarParam(fb, _OPTIONAL_, "FD_eps", &ts, 1, 1 ); CHKERRQ(ierr);	// eps for brute force FD gradients
			FDeps[i] = ts;

			gradar[i]    = 0.0;                     // GRADIENTS

			// PARAMETER VALUES
			if (par_val){
				// Add option to options database
				ierr = AddMaterialParameterToCommandLineOptions(par_str, ID, par_val); CHKERRQ(ierr);
				
				asprintf(&val_str, "%-9.4g", par_val); 
				//PetscPrintf(PETSC_COMM_WORLD,"Adding parameter to Options Database: %s \n", option);
			}
			else{
				asprintf(&val_str, "%-9s", "-"); 
			}

			if (vec_log10[i]==1){strcpy(logstr, "log10"); }
			else{strcpy(logstr, "     "); }
			if (FDgrad[i]==0){strcpy(adjointstr, "adjoint"); }
			else{strcpy(adjointstr, "FD     "); }
			


			
			// Print overview & indicate which parameters are not specified
			PetscPrintf(PETSC_COMM_WORLD, "|   %-2i: %s %s %+6s[%-2i] = %s; bnd=[%s; %s]   \n",i+1,adjointstr,logstr,par_str, ID,val_str,lb_str,ub_str);
			
			i = i+1;
		}

		fb->blockID++;
	}
	
    ierr  = PetscMemcpy(IOparam->grd,       	gradar,      (size_t)_MAX_PAR_*sizeof(PetscScalar) ); 				CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->type_name, 	type_name,   (size_t)_MAX_PAR_*(size_t)_str_len_*sizeof(char) ); 	CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->phs,       	phsar,       (size_t)_MAX_PAR_*sizeof(PetscInt)     ); 				CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->FD_gradient, 	FDgrad,      (size_t)_MAX_PAR_*sizeof(PetscInt)    ); 				CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->FD_eps,       	FDeps,       (size_t)_MAX_PAR_*sizeof(PetscScalar) ); 				CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->par_log10,     vec_log10,   (size_t)_MAX_PAR_*sizeof(PetscInt) ); 					CHKERRQ(ierr);


	VecRestoreArray(Adjoint_Vectors->P,&Par);
	VecRestoreArray(Adjoint_Vectors->Ub,&Ubar);
	VecRestoreArray(Adjoint_Vectors->Lb,&Lbar);
	VecRestoreArray(Adjoint_Vectors->grad,&gradar);
	IOparam->mdN = i;

	// Count # of FD gradients vs adjoint gradients and display it
	PetscInt numFD=0, numAdjoint=0;
	for(j = 0; j < IOparam->mdN; j++)
	{
		if  (IOparam->FD_gradient[j]==1){ numFD++; }
		else                            { numAdjoint++;}
	}
	PetscPrintf(PETSC_COMM_WORLD, "|   Total number of adjoint gradients      : %i   \n",numAdjoint);
	PetscPrintf(PETSC_COMM_WORLD, "|   Total number of FD gradients           : %i   \n",numFD);
	
	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// LOCATIONS
	// Get location / value / etc.
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<AdjointObservationPointStart>", "<AdjointObservationPointEnd>"); CHKERRQ(ierr);

	// error checking
	if(fb->nblocks >   _MAX_OBS_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many inverse indices specified! Max allowed: %lld", (LLD)  _MAX_OBS_);
	}
	if(!fb->nblocks && IOparam->Ap == 1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define sample points (mdi) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

	// Catch that if we have scaling none the index values don't go to inf:
	if(scal.utype == _NONE_)
	{
		scal.length   = 1;
		scal.velocity = 1;
	}

	// read each individual index
	if ( (fb->nblocks>0) & (IOparam->Ap==1)){
		PetscPrintf(PETSC_COMM_WORLD, "| \n|   Total number of observation points     : %i   \n", fb->nblocks);
	}
    else
    {
      	PetscPrintf(PETSC_COMM_WORLD, "| \n ");
    }

	ct1 = 0;
	ct2 = 0;   // used to check that only one observation type is used    	
    for(i = 0; i < fb->nblocks; i++)
	{
		// retrieve the coordinates of the sample points	
		ierr 	= getScalarParam(fb, _OPTIONAL_, "Coordinate",  IOparam->Coord, 3, 1);        CHKERRQ(ierr);       // not required if we compute 
		Ax[i] 	= (IOparam->Coord[0])/scal.length;
		Ay[i] 	= (IOparam->Coord[1])/scal.length;
		Az[i] 	= (IOparam->Coord[2])/scal.length;

		// Determine what cost function type is used
		ierr 	= getStringParam(fb, _REQUIRED_, "Parameter", ParType, NULL); CHKERRQ(ierr);  // must have component
		if     	(!strcmp(ParType, "Vx") | !strcmp(ParType, "Vy") | !strcmp(ParType, "Vz"))
		{    
			IOparam->MfitType = 0;  ct1++;	// these parameters are compute by projection 
		}
		else if	(!strcmp(ParType, "PSD") | !strcmp(ParType, "Exx") | !strcmp(ParType, "Eyy") | !strcmp(ParType, "Ezz" ) | 
				 !strcmp(ParType, "Exy") | !strcmp(ParType, "Eyz") | !strcmp(ParType, "Exz") | !strcmp(ParType, "E2nd"))
		{    
			IOparam->MfitType = 1;  ct2++;	// these parameters are compyted @ the center for the FDSTAG point
		}
		else
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [Vx,Vy,Vz,PSD,Exx,Eyy,Ezz,Exy,Exz,Eyz,E2nd] as Parameter\n");
		} 
		if (ct1 > 0 && ct2 > 0)  SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Only one Parameter for observation is allowed [Vx,Vy,Vz] OR [PSD,Exx,Eyy,Ezz,Exy,Eyz,Exz,E2nd]\n");
	
		strcpy(ObsName[i], ParType);		// store

		if (IOparam->MfitType == 0)
		{
			if     	(!strcmp(ParType, "Vx"))    ti=1;
			else if (!strcmp(ParType, "Vy"))    ti=2;
			else 								ti=3;
			Av[i] 	= ti;                     // VELOCITY COMPONENT
		
			if (IOparam->Gr==0)
			{
		    	ierr = getScalarParam(fb, _REQUIRED_, "Value", &ts, 1, 1 ); CHKERRQ(ierr);
			}
			else
			{
				ierr = getScalarParam(fb, _OPTIONAL_, "Value", &ts, 1, 1 ); CHKERRQ(ierr);
			}
			Ae[i] = ts /scal.velocity;     // VELOCITY VALUE

			if ((fb->nblocks<6) & (IOparam->Ap==1))
			{
				// Print overview 
				if (IOparam->Gr==0)
				{
					PetscPrintf(PETSC_COMM_WORLD, "|       [%f,%f,%f] has target velocity %s=%7.5f\n", IOparam->Coord[0],IOparam->Coord[1],IOparam->Coord[2], ParType, ts);  // cost function
				}
				else
				{
					PetscPrintf(PETSC_COMM_WORLD, "|       [%f,%f,%f] will compute gradient w.r.t. %s\n", IOparam->Coord[0],IOparam->Coord[1],IOparam->Coord[2], ParType);  // w.r.t. solution
				}
			}

		}
		else if(IOparam->MfitType == 1)	// parmeters defined at center
		{
			if (IOparam->Gr==0 )
			{
				ierr = getScalarParam(fb, _REQUIRED_, "Value", &ts, 1, 1 ); CHKERRQ(ierr);
			}
			else
			{
				ierr = getScalarParam(fb, _OPTIONAL_, "Value", &ts, 1, 1 ); CHKERRQ(ierr);
			}
			if   (!strcmp(ParType, "PSD")){ 
				Ae[i] = ts * 0.01745329251;     // PSD VALUE is given in degree; transfer to radians
			}
			else {
				Ae[i] = ts ;     				// value
			}
			Av[i] = 0;      // Placeholder

			if ((fb->nblocks<6) & (IOparam->Ap==1))
			{
				// Print overview 
				if (IOparam->Gr==0)
				{
					PetscPrintf(PETSC_COMM_WORLD, "|       [%f,%f,%f] has target %s=%7.5f\n", IOparam->Coord[0],IOparam->Coord[1],IOparam->Coord[2], ParType, ts);  // cost function
				}
				else
				{
					PetscPrintf(PETSC_COMM_WORLD, "|       [%f,%f,%f] will compute gradient w.r.t. %s \n", IOparam->Coord[0],IOparam->Coord[1],IOparam->Coord[2], ParType);  // w.r.t. solution
				}
			}
		}
		else
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [0; 1] as parameter for Adjoint_CostFunction, not %d",IOparam->MfitType);
		}

        if (IOparam->Ap>1)
		{
            const char *str_vec;
            if   (IOparam->Ap==2){str_vec="everywhere";}
            else {str_vec="at the internal free surface";}

            PetscPrintf(PETSC_COMM_WORLD, "|    We will compute the gradient %s w.r.t. V%s\n", str_vec, Vel_comp);  // inform where the gradient will be computed, and w.r.t. which component 
        }

		fb->blockID++;
	}
	PetscPrintf(PETSC_COMM_WORLD, "| \n");

    ierr  = PetscMemcpy(IOparam->Ax, 		Ax, 		(size_t)_MAX_OBS_*sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->Ay, 		Ay, 		(size_t)_MAX_OBS_*sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->Az, 		Az,			(size_t)_MAX_OBS_*sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->Av, 		Av, 		(size_t)_MAX_OBS_*sizeof(PetscInt));    CHKERRQ(ierr);
    ierr  = PetscMemcpy(IOparam->Ae, 		Ae, 		(size_t)_MAX_OBS_*sizeof(PetscScalar)); CHKERRQ(ierr);
	ierr  = PetscMemcpy(IOparam->ObsName, 	ObsName,   	(size_t)_MAX_OBS_*5*sizeof(char) ); 	CHKERRQ(ierr);

	IOparam->mdI = i;

	// Create the scaling for the cost function with the statistics of the observations
	IOparam->vel_scale = 1;   // Default 1
	if (IOparam->SCF == 1)
	{
		mean = 0;
		for(i=0;i<IOparam->mdI;i++)
		{
				mean += IOparam->Ae[i];
		}
		mean = abs(mean)/IOparam->mdI;
		IOparam->vel_scale = mean;
	}
	else if(IOparam->SCF == 2)
	{
		mean = 0;
		for(i=0;i<IOparam->mdI;i++)
		{
				mean += IOparam->Ae[i];
		}
		mean = abs(mean)/IOparam->mdI;
		var  = 0; 
		for(i=0;i<IOparam->mdI;i++)
		{
			var = pow(IOparam->Ae[i] - mean,2);
		}
		IOparam->vel_scale = var/(IOparam->mdI-1);
	}


	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// Error checking
	if (IOparam->use != 0 && !IOparam->Ap)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define observation points (Ap) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "LaMEMAdjointMain"
PetscErrorCode LaMEMAdjointMain(ModParam *IOparam)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscScalar    	F, *fcconvar, *Par;
	PetscInt        i;
	Adjoint_Vecs    Adjoint_Vectors;
	PetscLogDouble  cputime_start, cputime_end;

	PetscTime(&cputime_start);

	IOparam->count = 1;  // iteration counter for the initial cost function
	F              = 1e100;

	// Read input parameters
    ierr = LaMEMAdjointReadInputSetDefaults(IOparam, &Adjoint_Vectors); 			CHKERRQ(ierr);

	// Copy adjoint parameters to command-line	options
	VecCopy(Adjoint_Vectors.P,IOparam->P);		// Copy 
	VecGetArray(IOparam->P,&Par);
	for(i=0;i<IOparam->mdN;i++)
	{
		ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[i], i);			CHKERRQ(ierr);
    }
	VecRestoreArray(IOparam->P,&Par);


	//===========================================================
	// SOLVE ADJOINT (by calling LaMEMLibMain)
	//===========================================================
	// only compute the adjoint gradients or simply forward code
	if(IOparam->use == _adjointgradients_)
 	{
		// Compute Gradients 
		ierr = ComputeGradientsAndObjectiveFunction(Adjoint_Vectors.P, &F, Adjoint_Vectors.grad, IOparam);	CHKERRQ(ierr);

		// Print scaling laws
		ierr = PrintScalingLaws(IOparam);	CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint-based gradient inversion
 	else if(IOparam->use == _gradientdescent_)
 	{
 		if(IOparam->Tao == 1)
 		{
            // if TAO is employed, try the LMVM/BLMVM algorithms    
 	 		Tao tao;
			TaoLineSearch ls;

 	 		ierr = TaoCreate(PETSC_COMM_WORLD,&tao); CHKERRQ(ierr);

			PetscPrintf(PETSC_COMM_WORLD,"| TAO: IOparam->Ab bounds = %i \n",IOparam->Ab);	

 	 	 	// 1.Check if bounds are available and sets them
 	 		if (IOparam->Ab == 1)
 	 		{
 	 	 	 	ierr = TaoSetVariableBounds(tao,Adjoint_Vectors.Lb,Adjoint_Vectors.Ub);	 								CHKERRQ(ierr);
 	 	 	 	ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr);
 	 		}
 	 		else
 	 		{
 	 			ierr = TaoSetType(tao,TAOLMVM);CHKERRQ(ierr);                    // TAOLMVM, TAOBLMVM, TAOBMRM (bad), TAOCG (all 4), TAOTRON (might crash but is fast)
 	 		}

 	 		// 2. Set up Tao
 	 	 	ierr = TaoSetObjectiveAndGradientRoutine(tao, AdjointOptimisationTAO, IOparam);	 	CHKERRQ(ierr);  // sets the forward routine as well
 	 	 	ierr = TaoSetInitialVector(tao,Adjoint_Vectors.P);	 							    CHKERRQ(ierr);
 	 	 	ierr = TaoSetTolerances(tao,1e-30,1e-30,1e-30);	                                    CHKERRQ(ierr);
 	 	 	ierr = TaoSetFunctionLowerBound(tao,1e-10);                                          CHKERRQ(ierr);
 	 	 	ierr = TaoSetFromOptions(tao);	 										            CHKERRQ(ierr);
			
			// Line-Search
			ierr = TaoGetLineSearch(tao, &ls);                                                  CHKERRQ(ierr);
			ierr = TaoLineSearchSetFromOptions(ls);                                             CHKERRQ(ierr);

 	 	 	// 3. Solve Tao & view result
 	 	 	ierr = TaoSolve(tao);	 												            CHKERRQ(ierr);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");
	
 	 	 	TaoView(tao,PETSC_VIEWER_STDOUT_WORLD);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");
	
 	 	 	// 4. Clean
 	 	 	ierr = TaoDestroy(&tao);
 		}
 		else
 		{
            // Without TAO try our own implementation of a line search tuned gradient descent
			ierr = AdjointOptimisation(Adjoint_Vectors.P, F, Adjoint_Vectors.grad, IOparam);
 		}

 		PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");
		PetscPrintf(PETSC_COMM_WORLD,"| *                         INVERSION RESULT SUMMARY                      * \n");
		PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n");
		PetscPrintf(PETSC_COMM_WORLD,"| Number of inversion iterations: %d\n",IOparam->count);
 		PetscPrintf(PETSC_COMM_WORLD,"| F/Fini:\n");
 		VecGetArray(IOparam->fcconv,&fcconvar);
 		for(i=1;i<IOparam->count;i++)
 		{
 			PetscPrintf(PETSC_COMM_WORLD,"| %.5e\n",fcconvar[i]);
 		}
 		VecRestoreArray(IOparam->fcconv,&fcconvar);
 		PetscPrintf(PETSC_COMM_WORLD,"| \n| Final cost function:\n");
 		PetscPrintf(PETSC_COMM_WORLD,"| %.5e\n",IOparam->mfit);
 		PetscPrintf(PETSC_COMM_WORLD,"| \n| Final Parameters: \n");
		VecGetArray(IOparam->P,&Par);
		for(i=0;i<IOparam->mdN;i++)
		{
			
			PetscInt 		CurPhase;
			char 			CurName[_str_len_];
	
			CurPhase = 	IOparam->phs[i];			// phase of the parameter
    		strcpy(CurName, IOparam->type_name[i]);	// name

			PetscPrintf(PETSC_COMM_WORLD,"| %s[%i] = %.5e\n",CurName, CurPhase, Par[i]);
		}
		VecRestoreArray(IOparam->P,&Par);
 		PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------------------------------------- \n\n");
		
 	}
 	// this is a forward simulation that we want to save as comparison solution
 	else if(IOparam->use == _syntheticforwardrun_)
 	{
 		// call LaMEM main library function
 		ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);

 		// Save output
 		PetscViewer     viewerVel;
 		PetscViewerCreate(PETSC_COMM_WORLD,&viewerVel);
 		PetscViewerSetType(viewerVel,PETSCVIEWERBINARY);
 		PetscViewerFileSetMode(viewerVel,FILE_MODE_WRITE);
 		PetscViewerFileSetName(viewerVel,"Forward_Solution_Vel.bin");
 	 	VecView(IOparam->xini,viewerVel);
 	 	PetscViewerDestroy(&viewerVel);

 	 	PetscPrintf(PETSC_COMM_WORLD,"| ------------------------------------------\n|         Forward Solution successfully saved\n| ------------------------------------------\n");
 	}

	ierr = AdjointVectorsDestroy(&Adjoint_Vectors, IOparam); CHKERRQ(ierr);

	PetscTime(&cputime_end);
    PetscPrintf(PETSC_COMM_WORLD,"| Adjoint computation was successful & took %g s                         	 \n",cputime_end - cputime_start);
	PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************ \n");


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointVectorsCreate"
PetscErrorCode AdjointVectorsCreate(Adjoint_Vecs *Adjoint_Vectors, ModParam *IOparam)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// VECTORS
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Adjoint_Vectors->Lb);      CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Adjoint_Vectors->Ub);      CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Adjoint_Vectors->val);     CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Adjoint_Vectors->P);       CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Adjoint_Vectors->grad);    CHKERRQ(ierr);

	ierr = VecDuplicate(Adjoint_Vectors->P,&IOparam->P);											CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->maxit, PETSC_DETERMINE, &IOparam->fcconv);  	 	CHKERRQ(ierr); 


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointVectorsDestroy"
PetscErrorCode AdjointVectorsDestroy(Adjoint_Vecs *Adjoint_Vectors, ModParam *IOparam)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecDestroy(&Adjoint_Vectors->Lb);      	CHKERRQ(ierr);
	ierr = VecDestroy(&Adjoint_Vectors->Ub);      	CHKERRQ(ierr);
	ierr = VecDestroy(&Adjoint_Vectors->val);     	CHKERRQ(ierr);
	ierr = VecDestroy(&Adjoint_Vectors->P);       	CHKERRQ(ierr);
	ierr = VecDestroy(&Adjoint_Vectors->grad);    	CHKERRQ(ierr);

	ierr = VecDestroy(&IOparam->P);					CHKERRQ(ierr);
	ierr = VecDestroy(&IOparam->fcconv);			CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointCreate"
PetscErrorCode AdjointCreate(AdjGrad *aop, JacRes *jr, ModParam *IOparam)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Create everything
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vx);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vy);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vz);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->sty); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (jr->fs->DA_CEN, &aop->gradfield);      CHKERRQ(ierr);

	ierr = VecDuplicate(jr->gsol, &aop->dPardu);          CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &aop->dF);              CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &aop->pro);             CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &IOparam->xini);  	  CHKERRQ(ierr);  // create a new one

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointDestroy"
PetscErrorCode AdjointDestroy(AdjGrad *aop, ModParam *IOparam)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	
	ierr = VecDestroy(&aop->vx);         CHKERRQ(ierr);
	ierr = VecDestroy(&aop->vy);         CHKERRQ(ierr);
	ierr = VecDestroy(&aop->vz);         CHKERRQ(ierr);
	ierr = VecDestroy(&aop->sty);        CHKERRQ(ierr);
	ierr = VecDestroy(&aop->gradfield);  CHKERRQ(ierr);

	ierr = VecDestroy(&aop->dPardu);     CHKERRQ(ierr);
	ierr = VecDestroy(&aop->dF);         CHKERRQ(ierr);
	ierr = VecDestroy(&aop->pro);        CHKERRQ(ierr);
	ierr = VecDestroy(&IOparam->xini); 	 CHKERRQ(ierr); 

	// Destroy the Adjoint gradients structures
	// ierr = PetscMemzero(aop, sizeof(AdjGrad)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ComputeGradientsAndObjectiveFunction"
PetscErrorCode ComputeGradientsAndObjectiveFunction(Vec Parameters, PetscScalar *ObjectiveValue, Vec Gradient, ModParam *IOparam)
{
	/* 
		This computes the objective & gradients (either adjoint or finite difference 
		or a combination of both)
	*/
	PetscErrorCode 	ierr;
	PetscScalar 	*Par, *Grad;
	PetscInt 		i;

	PetscFunctionBegin;

	// Copy adjoint parameters to command-line	options
	VecCopy(Parameters,IOparam->P);		// Copy 
	VecGetArray(IOparam->P,&Par);
	for(i=0;i<IOparam->mdN;i++)
	{
		ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[i], i);			CHKERRQ(ierr);
    }
	VecRestoreArray(IOparam->P,&Par);

	/* 	FD gradients: call the FD gradient routine (& LaMEM many times) for those parameters for which we request it 
		Note: should be computed first, before computing adjoint to ensure that the cost function is computed for the standard (and nt perturbed) parameters
	*/
	IOparam->BruteForce_FD = PETSC_TRUE;
	ierr = AdjointFiniteDifferenceGradients(IOparam);	CHKERRQ(ierr);

	// Adjoint gradients: Call LaMEM main library function once (computes gradients @ the end)
 	IOparam->BruteForce_FD = PETSC_FALSE;
	ierr = LaMEMLibMain(IOparam); 															CHKERRQ(ierr);

	// Print overview of cost function & gradients 
	ierr = PrintCostFunction(IOparam);					CHKERRQ(ierr);
	ierr = PrintGradientsAndObservationPoints(IOparam); CHKERRQ(ierr);

	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr); 

	// Copy back gradients from IOparam   
	// Gradient info is thus in IOparam & as vectors
	VecGetArray(Gradient,&Grad);
	for(i=0;i<IOparam->mdN;i++){
		Grad[i] = IOparam->grd[i];		// Gradient
	}
	VecRestoreArray(Gradient,&Grad);

	*ObjectiveValue = IOparam->mfit;		// Objective function

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointOptimisation"
PetscErrorCode AdjointOptimisation(Vec P, PetscScalar F, Vec grad, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	PetscInt 		i, j, LScount;
	PetscScalar 	*Par, *Paroldar, *gradar, *gradoldar, *dPtemp, *fcconvar;
	PetscScalar   	Fold;
	ModParam    	*IOparam;
	IOparam     	= (ModParam*)ctx;
	Vec         	dP,dgrad,Pold,gradold,r;

	// get parameter values
	VecDuplicate(IOparam->P,&dP);
	VecDuplicate(IOparam->P,&Pold);
	VecDuplicate(grad,&gradold);
	VecDuplicate(grad,&dgrad);
	VecDuplicate(grad,&r);
	VecCopy(P,IOparam->P);
	VecCopy(P,Pold);

	// Initialize cost functions
	F 		= 1e100;
	Fold 	= 1e100;

	for(j = 0; j < IOparam->mdN; j++)
	{
		IOparam->factor2array[j] = 1;
	}

	while(F>IOparam->tol)
	{
		
		// Give the updated values to the code
		VecCopy(P,IOparam->P);		
		//VecGetArray(P,&Par);
		//for(i=0;i<IOparam->mdN;i++)
		//{
		//	ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[i], i);			CHKERRQ(ierr);
    	//}
		//VecRestoreArray(P,&Par);

		// Reset line search counter
		LScount = 1;

		// call LaMEM main library function
		//ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);
		ierr = ComputeGradientsAndObjectiveFunction(P, &F, grad, IOparam);	CHKERRQ(ierr);
				
		// Save initial cost function & create initial Hessian
		if(IOparam->count==1)
		{
			IOparam->mfitini = IOparam->mfit;
		}

		// Save cost function
		F = IOparam->mfit;

		VecGetArray(P,&Par);
		VecGetArray(grad,&gradar);
		for(j = 0; j < IOparam->mdN; j++)
		{
			gradar[j] = IOparam->grd[j];
			if (IOparam->count==1)
			{
				IOparam->factor2array[j] = fabs((IOparam->Scale_Grad*Par[j])/fabs(gradar[j]));
			}
		}
		VecRestoreArray(grad,&gradar);
		VecRestoreArray(P,&Par);

		PetscPrintf(PETSC_COMM_WORLD,"| AdjointOptimisation: Gradients. [0]=%e, [1]=%e \n", IOparam->grd[0], IOparam->grd[1]);
		
		// If cost function in this timestep is larger then before perform bisection line search (striclty speacking we only require to evaluate the cost function during LS..)
		while(F>Fold)
		{

			// Reduce line search alpha
			for(j = 0; j < IOparam->mdN; j++)
			{
				IOparam->factor2array[j] *= IOparam->facB;
			}

			PetscPrintf(PETSC_COMM_WORLD,"\n| - - - - - - - - - - - - - - - - - - - - - - - - - - - \n");
			PetscPrintf(PETSC_COMM_WORLD,"|               LINE SEARCH IT %d                       \n",LScount);

			VecGetArray(P,&Par);
			VecGetArray(Pold,&Paroldar);
			VecGetArray(gradold,&gradoldar);

			/*VecGetArray(dP,&dPtemp);
			for(i=0;i<IOparam->mdN;i++)
			{
				for(j=0;j<IOparam->mdN;j++)
				{
					dPtemp[i] = dPtemp[i]*LSalpha;
				}
			}
			*/

			// Update parameter
			for(i=0;i<IOparam->mdN;i++)
			{
				Par[i] 	= 	Paroldar[i] - gradoldar[i] * IOparam->factor2array[i] ;
			//	ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[i], i);			CHKERRQ(ierr);
    		}
			VecRestoreArray(P,&Par);
			VecRestoreArray(Pold,&Paroldar);
			VecRestoreArray(gradold,&gradoldar);
			//VecRestoreArray(dP,&dPtemp);

			// Store the updated values (now done inside the routine below)
			//VecCopy(P,IOparam->P);  

			// call LaMEM main library function
			//ierr 	= LaMEMLibMain(IOparam); CHKERRQ(ierr);
			ierr = ComputeGradientsAndObjectiveFunction(P, &F, grad, IOparam);	CHKERRQ(ierr);
			//F 		= IOparam->mfit;

			LScount+=1;
			if(LScount>IOparam->maxitLS)
			{
				PetscPrintf(PETSC_COMM_WORLD,"| ******************************************************\n");
				PetscPrintf(PETSC_COMM_WORLD,"| *              SOLUTION DIVERGED                     *\n");
				PetscPrintf(PETSC_COMM_WORLD,"| ******************************************************\n\n");

				// Return parameters for final output
				VecCopy(P,IOparam->P);

				VecDestroy(&dP);
				VecDestroy(&Pold);
				VecDestroy(&gradold);
				VecDestroy(&dgrad);
				VecDestroy(&r);

				PetscFunctionReturn(0);
			}

			PetscPrintf(PETSC_COMM_WORLD,"|    F = %10.6e,  Fold = %10.6e                      \n",F,Fold);
			PetscPrintf(PETSC_COMM_WORLD,"\n| - - - - - - - - - - - - - - - - - - - - - - - - - - - \n");

		}

		PetscPrintf(PETSC_COMM_WORLD,"\n| ------------------------------------------------------------------------ \n");
		PetscPrintf(PETSC_COMM_WORLD,"| %d. IT INVERSION RESULT: line search its = %d ; F / FINI = %.5e\n| \n",IOparam->count,LScount-1,IOparam->mfit/IOparam->mfitini);
		PetscPrintf(PETSC_COMM_WORLD,"| Fold = %.5e \n|    F = %.5e\n| \n",Fold,F);

		VecGetArray(grad,&gradar);
		VecGetArray(P,&Par);
		// VecGetArray(dP,&dPtemp);

		// Store old values in case we overshoot next time
		VecCopy(P,Pold);
		VecCopy(grad,gradold);
		Fold = F;

		// Display the current state of the parameters
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"| %D. Diff parameter value = %.5e\n",j+1,-gradar[j] * IOparam->factor2array[j]);
		}

		PetscPrintf(PETSC_COMM_WORLD,"| \n");


		// Update parameter
		for(i=0;i<IOparam->mdN;i++)
		{
			if(F>IOparam->tol)
			{
				Par[i] 	= 	Par[i] -gradar[i] * IOparam->factor2array[i];
				ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[i], i);			CHKERRQ(ierr);
    		}
		}
		VecRestoreArray(grad,&gradar);
		VecRestoreArray(P,&Par);
		// VecRestoreArray(dP,&dPtemp);

		// Display the current state of the parameters
		VecGetArray(P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"| %D. Parameter value = %.5e\n",j+1,Par[j]);
		}
		VecRestoreArray(P,&Par);


		PetscPrintf(PETSC_COMM_WORLD,"| -------------------------------------------------------------------------\n\n");
		

		// Give the updated values to the code  (actually unfortunately necessary here and at the top of this function - need to rearrange that)
		VecCopy(P,IOparam->P);

		VecGetArray(IOparam->fcconv,&fcconvar);
		fcconvar[IOparam->count] = IOparam->mfit/IOparam->mfitini;
		VecRestoreArray(IOparam->fcconv,&fcconvar);

		// Increase line search alpha after succesful iteration
		for(i=0;	i<IOparam->mdN;	i++)
		{
			// dPtemp[i] = - (dPtemp[i] + (gradar[i])) * IOparam->factor2array[i];
			IOparam->factor2array[i] *= IOparam->facLS;
			if(IOparam->factor2array[i]>IOparam->maxfac)
			{
				IOparam->factor2array[i] = IOparam->maxfac;
			}
			PetscPrintf(PETSC_COMM_WORLD,"| LS factor for %d.Parameter = %.5e\n",i+1,IOparam->factor2array[i]);
		}	
		PetscPrintf(PETSC_COMM_WORLD,"| \n");

		// count
		IOparam->count += 1;
		if(IOparam->count>IOparam->maxit)
		{
			PetscPrintf(PETSC_COMM_WORLD,"\n\n| Maximum number of inverse iterations reached\n\n");
			break;
		}
	}

	// clean up vectors
	VecDestroy(&dP);
	VecDestroy(&Pold);
	VecDestroy(&gradold);
	VecDestroy(&dgrad);
	VecDestroy(&r);


	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointOptimisationTAO"
PetscErrorCode AdjointOptimisationTAO(Tao tao, Vec P, PetscReal *F, Vec grad, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	PetscInt j;
	PetscScalar *Par, *gradar, *fcconvar, Fmisfit;
	ModParam    *IOparam;
	char		CurName[_str_len_];
	
	IOparam    =  (ModParam*)ctx;

	// get parameter values
	VecCopy(P,IOparam->P);

	// Print parameters
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| *************************************************************************\n");
	PetscPrintf(PETSC_COMM_WORLD,"| TAO start iteration %i: \n", IOparam->count);
	PetscPrintf(PETSC_COMM_WORLD,"| Currently employed parameters: \n");
	VecGetArray(IOparam->P,&Par);
	for(j = 0; j < IOparam->mdN; j++){
	
		//ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[j], j);					CHKERRQ(ierr);
    
	    strcpy(CurName, IOparam->type_name[j]);	// name
		PetscPrintf(PETSC_COMM_WORLD,"|  %s[%i]=%10.10e \n",CurName,IOparam->phs[j],Par[j]);
	}
	VecRestoreArray(IOparam->P,&Par);
	PetscPrintf(PETSC_COMM_WORLD,"| \n");

	// Compute gradient & objective function for the parameters P
	ierr = ComputeGradientsAndObjectiveFunction(P, &Fmisfit, grad, IOparam);	CHKERRQ(ierr);
		
	//ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);

	// Print the gradient & misfit
	for(j = 0; j < IOparam->mdN; j++)
	{
	//	gradar[j] = IOparam->grd[j];
		
	}

	*F = IOparam->mfit;  // objective function

	// Save initial cost function
	if(IOparam->count==1){IOparam->mfitini = IOparam->mfit;}

	// Display the current state of the parameters and gradient
	PetscPrintf(PETSC_COMM_WORLD,"| *************************************************************************\n");
	PetscPrintf(PETSC_COMM_WORLD,"| TAO results of iteration %i: \n", IOparam->count);
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| Parameter values: \n");
	
	VecGetArray(IOparam->P,&Par);
	VecGetArray(grad,&gradar);
	for(j = 0; j < IOparam->mdN; j++)
	{
	//	ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[j], j);					CHKERRQ(ierr);

		strcpy(CurName, IOparam->type_name[j]);	// name

		PetscPrintf(PETSC_COMM_WORLD,"|   %D. %s[%i] = %- 10.5e, gradient=%- 10.5e\n",j+1,CurName,IOparam->phs[j],Par[j],gradar[j]);
		
	}
	VecRestoreArray(grad,&gradar);
	VecRestoreArray(IOparam->P,&Par);

	// Print cost function 
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| misfit           = %2.8e \n",IOparam->mfit);
	PetscPrintf(PETSC_COMM_WORLD,"| misfit / misfit0 = %2.8e\n| ------------------------------------------\n\n",IOparam->mfit/IOparam->mfitini);

	VecGetArray(IOparam->fcconv,&fcconvar);
	fcconvar[IOparam->count] = IOparam->mfit/IOparam->mfitini;
	VecRestoreArray(IOparam->fcconv,&fcconvar);

	// count
	IOparam->count += 1;
	if(IOparam->count>1500)
	{
		PetscPrintf(PETSC_COMM_WORLD,"\n\n\n| EXCEEDED 1500 FUNCTION EVALUATIONS (consider changing inversion options)\n\n\n");
		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
 //---------------------------------------------------------------------------
 #undef __FUNCT__
 #define __FUNCT__ "AdjointObjectiveAndGradientFunction"
 PetscErrorCode AdjointObjectiveAndGradientFunction(AdjGrad *aop, JacRes *jr, NLSol *nl, ModParam *IOparam, SNES snes, FreeSurf *surf)
 {
	// This computes the objective function and adjoint gradients (not the 'brute-force' FD gradients)

 	PetscErrorCode ierr;
 	PetscFunctionBegin;

	//========================================
	// COMPUTE OBJECTIVE FUNCTION & GRADIENT
	//========================================

	// Get the cost function and derivative
	ierr = AdjointObjectiveFunction(aop, jr, IOparam, surf); CHKERRQ(ierr);

	if (IOparam->BruteForce_FD){PetscFunctionReturn(0); }	// in case we compute Brute force FD, we can return at this stage

 	// Get the adjoint gradients
 	ierr = AdjointComputeGradients(jr, aop, nl, snes, IOparam, surf);        CHKERRQ(ierr);

 	PetscFunctionReturn(0);
 }
 //---------------------------------------------------------------------------
 #undef __FUNCT__
 #define __FUNCT__ "AdjointObjectiveFunction"
 PetscErrorCode AdjointObjectiveFunction(AdjGrad *aop, JacRes *jr, ModParam *IOparam, FreeSurf *surf)
 {
	// This computes the objective function

	Scaling             *scal;
	Vec                  xini, sqrtpro;
	PetscScalar          Ad; 

 	PetscErrorCode ierr;
 	PetscFunctionBegin;

 	scal = jr->scal;

	// Create projection vector
	ierr = VecDuplicate(jr->gsol, &xini);           CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sqrtpro);           CHKERRQ(ierr);

	//========================================
	// COMPUTE OBJECTIVE FUNCTION
	//========================================
	// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
	ierr = AdjointPointInPro(jr, aop, IOparam, surf); 	CHKERRQ(ierr);

	ierr = VecCopy(IOparam->xini,xini);             	CHKERRQ(ierr);

	if (IOparam->Gr == 1)		// Gradients w.r.t. Solution
	{
		PetscScalar value;

		if(IOparam->MfitType == 0)
		{
			// Velocity

			// -------- Only get gradients with respect to the solution --------
			ierr = VecCopy(aop->pro,aop->dF); 				CHKERRQ(ierr); // dF/dx = P

			ierr = VecDot(aop->pro, jr->gsol,&value);     	CHKERRQ(ierr);
			IOparam->mfit 	   = value*scal->velocity;    
		}
		else if(IOparam->MfitType == 1)
		{
			// principal stress direction
			ierr 				= AdjointGet_F_dFdu_Center(jr, aop, IOparam);                       CHKERRQ(ierr);

			if (!strcmp(IOparam->ObsName[0],"PSD")){
				IOparam->mfit 	  	= 	IOparam->mfitCenter; 						// stress angle is nondimensional
			}
			else{
				IOparam->mfit 	  	= 	IOparam->mfitCenter*(1.0/scal->time_si); 	// strain rate is dimensional [SI]
			}
			PetscPrintf(PETSC_COMM_WORLD,"| IOparam->mfit = %e \n",IOparam->mfit);
		

		}
	}
	else if(IOparam->Gr == 0)	// Gradients w.r.t. CostFunction
	{

		if(IOparam->MfitType == 0)
		{
			ierr = VecCopy(aop->pro,sqrtpro);             	CHKERRQ(ierr);

			// Incorporate projection vector (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
			ierr = VecAYPX(xini,-1,jr->gsol);                                       CHKERRQ(ierr);
			ierr = VecScale(xini,sqrt(1/IOparam->vel_scale));                   	CHKERRQ(ierr);

			ierr =  VecSqrtAbs(sqrtpro);  CHKERRQ(ierr);

			ierr = VecPointwiseMult(xini, xini,sqrtpro);                           CHKERRQ(ierr);
			
			// Compute objective function value (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
			ierr 	           = VecDot(xini,xini,&Ad);
			Ad 		          /= 2;
			IOparam->mfit 	   = Ad*pow(scal->velocity,2); // Dimensional misfit function

			ierr = VecCopy(IOparam->xini,xini);             	CHKERRQ(ierr);

			// Incorporate projection vector (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
			ierr = VecAYPX(xini,-1,jr->gsol);                                       CHKERRQ(ierr);
			ierr = VecScale(xini,1/IOparam->vel_scale);                             CHKERRQ(ierr);
			ierr = VecPointwiseMult(xini, xini,aop->pro);                           CHKERRQ(ierr);

			// Compute it's derivative (dF/dx = P*x-P*x_ini)
			ierr = VecCopy(xini,aop->dF); 	
		}
		if(IOparam->MfitType == 1)
		{
			ierr 				= 	AdjointGet_F_dFdu_Center(jr, aop, IOparam);                       CHKERRQ(ierr);
			if (!strcmp(IOparam->ObsName[0],"PSD")){
				IOparam->mfit 	  	= 	IOparam->mfitCenter; 								// stress angle is nondimensional
			}
			else{
				IOparam->mfit 	  	= 	IOparam->mfitCenter*pow(1.0/scal->time_si,2.0); 	// strain rate is dimensional [SI]
			}
			PetscPrintf(PETSC_COMM_WORLD,"| IOparam->mfit = %e \n",IOparam->mfit);
		}
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"| ERROR choose Inv_Gr = 0 or = 1\n");
	}

	// Destroy
	ierr = VecDestroy(&xini);
	ierr = VecDestroy(&sqrtpro);


 	PetscFunctionReturn(0);
 }
//---------------------------------------------------------------------------
/* 
	Brute Force finite difference computation of the gradient, by calling LaMEM twice & perturbing the parameter 
*/
#undef __FUNCT__
#define __FUNCT__ "AdjointFiniteDifferenceGradients"
PetscErrorCode AdjointFiniteDifferenceGradients(ModParam *IOparam)
{
	PetscErrorCode 	ierr;
	PetscInt 		j;
	PetscScalar 	*Par, Perturb, Misfit_ref, Misfit_pert, FD_gradients_eps=1e-6, FD_eps, Grad; 
	char 			CurName[_str_len_];
	PetscBool 		flg, FD_Adjoint = PETSC_FALSE;


	// 0) Retrieve (optional) command-line parameters
	ierr = PetscOptionsGetScalar(NULL, NULL,"-FD_gradients_eps",&FD_gradients_eps,&flg); CHKERRQ(ierr);
    if (flg){
		PetscPrintf(PETSC_COMM_WORLD,"| Updated eps used for computing finite difference gradients to: %2.5e  \n",FD_gradients_eps);
    }

	// 1) Compute 'reference' state using LaMEM & the current set of parameters
	// Set parameters as command-line options
	VecGetArray(IOparam->P,&Par);
	for(j = 0; j < IOparam->mdN; j++){
		ierr	=	CopyParameterToLaMEMCommandLine(IOparam,  Par[j], j);					CHKERRQ(ierr);
    }
	VecRestoreArray(IOparam->P,&Par);

	// check if we actually need to compute FD gradients
	for(j = 0; j < IOparam->mdN; j++){
		if (IOparam->FD_gradient[j]>0){	// only if we want to compute a FD gradient for this paramater
			FD_Adjoint = PETSC_TRUE;
		}
	}

	if (FD_Adjoint){

		// Call LaMEM
		ierr 		= 	LaMEMLibMain(IOparam); CHKERRQ(ierr);		// call LaMEM
		Misfit_ref	=	IOparam->mfit;
		
		PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************ \n");
		PetscPrintf(PETSC_COMM_WORLD,"|                       FINITE DIFFERENCE GRADIENTS                        \n");
		PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************ \n");
		PetscPrintf(PETSC_COMM_WORLD,"| Reference objective function: %- 2.6e \n",IOparam->mfit);


		// 2) Loop over all parameters & compute a solution  with the perturbed parameters
		
		// Set parameters as command-line options
		VecGetArray(IOparam->P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscInt 	CurPhase;
			PetscScalar CurVal;

			if (IOparam->FD_gradient[j]>0)
			{	// only if we want to compute a FD gradient for this paramater
					
				CurPhase 		= 	IOparam->phs[j];
				CurVal 	 		= 	Par[j];
				FD_eps 			=	IOparam->FD_eps[j];
				if (FD_eps==0.0)
				{
					FD_eps 		=	FD_gradients_eps;	// use default value
				}
				strcpy(CurName, IOparam->type_name[j]);	// name

				// Only execute this for those parameters for which we want to know the FD gradient (to be added here)
				Perturb 	= 	CurVal*FD_eps;

				ierr		=	CopyParameterToLaMEMCommandLine(IOparam, CurVal + Perturb, j);		CHKERRQ(ierr);
			
				// Compute solution with updated parameter
				ierr 		= 	LaMEMLibMain(IOparam); 												CHKERRQ(ierr);
				Misfit_pert = 	IOparam->mfit;

				// FD gradient
				Grad 			=	(Misfit_pert-Misfit_ref)/Perturb;
				IOparam->grd[j] = 	Grad;									// store gradient

				// Set back parameter
				ierr			=	CopyParameterToLaMEMCommandLine(IOparam,  CurVal, j);		CHKERRQ(ierr);
			
				PetscPrintf(PETSC_COMM_WORLD,"|  Perturbed Misfit value     : %- 2.6e \n", Misfit_pert);
				PetscPrintf(PETSC_COMM_WORLD,"|  Brute force FD gradient %+5s[%2i] = %e, with eps=%1.4e \n", CurName, CurPhase, Grad, FD_eps);
			}
		}
		VecRestoreArray(IOparam->P,&Par);
		IOparam->mfit = Misfit_ref; // reset value, just in case it is overwritten by a perturbed value
	}


 	PetscFunctionReturn(0);
 }


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointComputeGradients"
PetscErrorCode AdjointComputeGradients(JacRes *jr, AdjGrad *aop, NLSol *nl, SNES snes, ModParam *IOparam, FreeSurf *surf)
{
	
	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG              *fs;
	KSP                 ksp_as;
	KSPConvergedReason  reason;
	PetscInt            i, j, lrank, grank;
	PetscScalar         dt, grd, Perturb, coord_local[3], *vx, *vy, *vz, *Par, *sty, CurVal;
	Vec 				res_pert, sol, psi, psiPar, drdp, res;
	PC                  ipc_as;
	Scaling             *scal;
    PetscBool           flg;
    char                CurName[_str_len_];
	PetscMPIInt    		rank;
	PetscScalar    		*rbuf1=PETSC_NULL;

	fs = jr->fs;
	dt = jr->ts->dt;

	scal = jr->scal;

	// Create all needed vectors in the same size as the solution vector
	ierr = VecDuplicate(jr->gsol, &psi);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &psiPar);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &res_pert);    CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &res);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sol); 		 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &drdp);	 	 CHKERRQ(ierr);
	ierr = VecCopy(jr->gsol,sol); 				 CHKERRQ(ierr);
	ierr = VecCopy(jr->gres,res); 				 CHKERRQ(ierr);;

	//========
	// SOLVE
	//========
	// Solve the adjoint equation (psi = J^-T * dF/dx)
	// (A side note that I figured out, ksp still sometimes results in a > 0 gradient even if cost function is zero.. possibly really bad condition number?)
	if(IOparam->MfitType == 0)
	{
		ierr = SNESGetKSP(snes, &ksp_as);         		CHKERRQ(ierr);
		ierr = KSPSetOptionsPrefix(ksp_as,"as_"); 		CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp_as);         		CHKERRQ(ierr);
		ierr = KSPGetPC(ksp_as, &ipc_as);           	CHKERRQ(ierr);
		ierr = PCSetType(ipc_as, PCMAT);          		CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp_as,nl->J,nl->P);		CHKERRQ(ierr);
		ierr = KSPSolve(ksp_as,aop->dF,psi);			CHKERRQ(ierr);
		ierr = KSPGetConvergedReason(ksp_as,&reason);	CHKERRQ(ierr);
	}
	else if(IOparam->MfitType == 1)
 	{
 		ierr = SNESGetKSP(snes, &ksp_as);         		CHKERRQ(ierr);
 		ierr = KSPSetOptionsPrefix(ksp_as,"as_"); 		CHKERRQ(ierr);
 		ierr = KSPSetFromOptions(ksp_as);         		CHKERRQ(ierr);
 		ierr = KSPGetPC(ksp_as, &ipc_as);            	CHKERRQ(ierr);
 		ierr = PCSetType(ipc_as, PCMAT);          		CHKERRQ(ierr);
 		ierr = KSPSetOperators(ksp_as,nl->J,nl->P);	CHKERRQ(ierr);
 		ierr = KSPSolve(ksp_as,aop->dPardu,psiPar);	CHKERRQ(ierr);
 		ierr = KSPGetConvergedReason(ksp_as,&reason);	CHKERRQ(ierr);
 	}

    // Set the FD step-size for computing dr/dp (or override it with a command-line option, which is more for advanced users/testing)
    aop->FD_epsilon = 1e-6;
    ierr = PetscOptionsGetScalar(NULL, NULL,"-FD_epsilon_adjoint",&aop->FD_epsilon,&flg); CHKERRQ(ierr);
    if (flg)
    {
        PetscPrintf(PETSC_COMM_WORLD,"|     Finite difference step size for Adjoint dr/dp calculation set to %e \n", aop->FD_epsilon);
    }


	// Field sensitivity (aka geodynamic sensitivity kernels) or 'classic' phase based gradients?
 	if(IOparam->FS == 1)
	{
        /*  This plots the sensitivity of the solution to incremental changes in the parameter (density here) at 
             every point in the computational domain. This allows a visual inspection of where changes in density have the largest impact 
             on say, surface velocity
        */

		// Compute residual with the converged Jacobian analytically
		
    	strcpy(CurName, IOparam->type_name[0]);	// name of the parameter
		if (IOparam->mdN>1)
		{
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"| Field-based gradients can currently only be computed for a single parameter \n");		// error check
		}
		if (!strcmp(CurName,"rho") | !strcmp(CurName,"eta") | !strcmp(CurName,"n"))		// we need some way to ensure that we do this for one field at a time only
		{
			PetscPrintf(PETSC_COMM_WORLD,"| Starting computation of Field-based gradients.  \n");	
			// Compute the gradient
			if(IOparam->MfitType == 0)
			{
				aop->CurScal   = (scal->velocity)/(1);
				ierr = AdjointFormResidualFieldFD(snes, sol, psi, nl, aop, IOparam);          CHKERRQ(ierr);
			}
			else if(IOparam->MfitType == 1)
			{
				if (!strcmp(IOparam->ObsName[0],"PSD")){
					// PSD is dimensionless
					aop->CurScal    =   1.0;
				}
				else{
					// if NOT PSD, it should be strain rate and have units of 1/s. 
					aop->CurScal =    (1.0/scal->time_si);	
				}

				ierr = AdjointFormResidualFieldFD(snes, sol, psiPar, nl, aop, IOparam);          CHKERRQ(ierr);
			}
			PetscPrintf(PETSC_COMM_WORLD,"| Finished gradient computation & added it to VTK\n");
			PetscPrintf(PETSC_COMM_WORLD,"| Add '-out_gradient = 1' to your parameter file. \n");
		}
		else 
		{
			PetscPrintf(PETSC_COMM_WORLD,"| Field based gradient only for [rho,eta,n] programmed not for %s! \n",CurName);
		}
	}
	else // Phase based gradients
	{
 		//=================
		// PARAMETER LOOP
		//=================
		VecGetArray(IOparam->P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			if (!IOparam->FD_gradient[j]){	// only if we want to compute an adjoint gradient for this parameter

				// Get the initial residual since it is overwritten in VecAYPX
				ierr = VecCopy(jr->gres,res); 			CHKERRQ(ierr);

				// Get current phase and parameter which is being perturbed
				//CurPhase 		= 	IOparam->phs[j];
				CurVal 	 		= 	Par[j];
				strcpy(CurName, IOparam->type_name[j]);	// name

				// Perturb parameter
				Perturb 		= 	aop->FD_epsilon*CurVal;

				//PetscPrintf(PETSC_COMM_WORLD,"| *** dr/dp: Perturbing parameter %s[%i]=%e to value %e  \n",CurName,CurVal, CurVal + Perturb);

				// Set as command-line option & create updated material database
				ierr 			=	CopyParameterToLaMEMCommandLine(IOparam,  CurVal + Perturb, j);			CHKERRQ(ierr);
				ierr 			= 	CreateModifiedMaterialDatabase(IOparam);     							CHKERRQ(ierr);			// update LaMEM material DB (to call directly call the LaMEM residual routine)

				// Swap material structure of phase with that of LaMEM Material DB
				for (i=0; i < nl->pc->pm->jr->dbm->numPhases; i++)
				{
					ierr =   PetscMemzero(&nl->pc->pm->jr->dbm->phases[i],  sizeof(Material_t));   CHKERRQ(ierr);
					swapStruct(&nl->pc->pm->jr->dbm->phases[i], &IOparam->dbm_modified.phases[i]);  
				}

				ierr 			= 	FormResidual(snes, sol, res_pert, nl);         							CHKERRQ(ierr);        // compute the residual with the perturbed parameter
				ierr 			=	VecAYPX(res,-1.0,res_pert);                      							CHKERRQ(ierr);        // res = (res_perturbed-res)
				ierr 			=	VecScale(res,1.0/Perturb);                      					CHKERRQ(ierr);        // res = (res_perturbed-res)/Perturb
				
				ierr 			=	VecCopy(res,drdp);														CHKERRQ(ierr);        
				
				// Reset parameter again
				ierr 			=	CopyParameterToLaMEMCommandLine(IOparam,  CurVal, j);					CHKERRQ(ierr);
				ierr 			= 	CreateModifiedMaterialDatabase(IOparam);     							CHKERRQ(ierr);			// update LaMEM material DB (to call directly call the LaMEM residual routine)

				// Swap material structure of phase with that of LaMEM Material DB back
				for (i=0; i < nl->pc->pm->jr->dbm->numPhases; i++)
				{
					ierr =   PetscMemzero(&nl->pc->pm->jr->dbm->phases[i],  sizeof(Material_t));   CHKERRQ(ierr);
					swapStruct(&nl->pc->pm->jr->dbm->phases[i], &IOparam->dbm_modified.phases[i]);  
				}
				
				// Compute the gradient (dF/dp = -psi^T * dr/dp) & Save gradient
				if (IOparam->MfitType == 0)
				{
					ierr          	=   VecDot(drdp,psi,&grd);                       CHKERRQ(ierr);
					if (IOparam->Gr == 1)	
					{
						aop->CurScal = scal->velocity;
					}
					else if (IOparam->Gr == 0)	
					{
						aop->CurScal = pow(scal->velocity,2);
					}
				}
				else if (IOparam->MfitType == 1)		// PSD or strainrate 
				{
					ierr          	=   VecDot(drdp,psiPar,&grd);                       CHKERRQ(ierr);
				
					if (!strcmp(IOparam->ObsName[0],"PSD")){
						// PSD is dimensionless
						aop->CurScal    =   1.0;
					}
					else{
						// if NOT PSD, it should be strain rate and have units of 1/s. 
						// if we later add other variables at the center points (e.g., stress, this needs to be expanded)
						if	 	(IOparam->Gr == 1){		aop->CurScal =    (1.0/scal->time_si);		}
						else if (IOparam->Gr == 0){		aop->CurScal = pow(1.0/scal->time_si,2.0);	}
					}
					PetscPrintf(PETSC_COMM_WORLD,"| grad=%e, aop->CurScal=%e vel-scale =%e\n",grd, aop->CurScal, scal->length);
				}
					
				IOparam->grd[j]	=   -grd*aop->CurScal;						// gradient

			}
		}
		VecRestoreArray(IOparam->P,&Par);
	}

	if(IOparam->mdI<_MAX_OBS_ && IOparam->Ap == 1)
	{

		// get the current velocities at the observation point
		ierr = AdjointPointInPro(jr, aop, IOparam, surf);    CHKERRQ(ierr);

		VecGetArray(aop->vx,&vx);
		VecGetArray(aop->vy,&vy);
		VecGetArray(aop->vz,&vz);

		// Print the solution variable at the user defined index (if there are sufficiently few)
		for (i=0; i<IOparam->mdI; i++)
		{

			coord_local[0] 				= 	IOparam->Ax[i];
			coord_local[1] 				= 	IOparam->Ay[i];
			coord_local[2] 				= 	IOparam->Az[i];
			IOparam->Apoint_on_proc[i]	=	PETSC_FALSE;
				

			// get global & local ranks of a marker
			ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

			IOparam->Avel_num[i]	= 0.0;
			
			// If lrank is not 13 the point is not on this processor
			if(lrank == 13)
			{
                char vel_com[20];
                PetscScalar vel=0;
	
	
                if (IOparam->Av[i] == 1){strcpy(vel_com, "Vx"); vel = vx[i]*scal->velocity;}
                if (IOparam->Av[i] == 2){strcpy(vel_com, "Vy"); vel = vy[i]*scal->velocity;}
                if (IOparam->Av[i] == 3){strcpy(vel_com, "Vz"); vel = vz[i]*scal->velocity;}
               IOparam->Apoint_on_proc[i]=PETSC_TRUE;

				if(IOparam->MfitType == 0)
				{
					IOparam->Avel_num[i]  = vel;
				}
				else if (IOparam->MfitType == 1)
				{
					VecGetArray(aop->sty,&sty);
					IOparam->Avel_num[i]  = sty[i];
					VecRestoreArray(aop->sty,&sty);
				}

				if (IOparam->Adv == 1)     // advect the point? - Note that this should be done only once per timestep; need to ensure that this is still the case if 
				{
					IOparam->Ax[i] += vx[i]*dt;
					IOparam->Ay[i] += vy[i]*dt;
					IOparam->Az[i] += vz[i]*dt;
				}
			}
		}

		VecRestoreArray(aop->vx,&vx);
		VecRestoreArray(aop->vy,&vy);
		VecRestoreArray(aop->vz,&vz);

		// Send velocity array to rank 0 (sum over all arrays), to deal with points residing on different processors
		{

			MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

			if ( rank == 0) 
			{ 
       			rbuf1 = (PetscScalar *)malloc(_MAX_PAR_*sizeof(PetscScalar)); 
       		} 

			// send from all processors -> rank 0
			MPI_Reduce( IOparam->Avel_num, rbuf1, _MAX_PAR_, MPIU_SCALAR, MPI_SUM, 0, PETSC_COMM_WORLD);
   			
			if ( rank == 0)
			{ 
				ierr  = PetscMemcpy(IOparam->Avel_num,   rbuf1,  (size_t)_MAX_PAR_*sizeof(PetscScalar) ); CHKERRQ(ierr);		// copy array to correct point
			}

			// send from rank 0 to all other processors
			MPI_Bcast(IOparam->Avel_num, _MAX_PAR_, MPIU_SCALAR, 0, PETSC_COMM_WORLD);	

			if ( rank == 0)
			{ 
				free(rbuf1);
			}
		}
	}
 

 	// Print overview of gradients (if requested @ this stage)
	//ierr = PrintGradientsAndObservationPoints(IOparam); CHKERRQ(ierr);

	// Clean
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&psiPar);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&drdp);
	ierr = VecDestroy(&res_pert);
	ierr = VecDestroy(&res);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PrintCostFunction"
/*
    This prints the current value of the objective function
*/
PetscErrorCode PrintCostFunction(ModParam *IOparam)
{
	PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"|                       COMPUTATION OF THE COST FUNCTION                    \n");
    PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************\n");

	PetscPrintf(PETSC_COMM_WORLD,"| Current Cost function = %2.10e\n",IOparam->mfit);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PrintGradientsAndObservationPoints"
/*
    This prints an overview of the gradients (adjoint & fd) and the velocity values
*/
PetscErrorCode PrintGradientsAndObservationPoints(ModParam *IOparam)
{
	PetscErrorCode 	ierr;
	char 			CurName[_str_len_], logstr[_str_len_];
	PetscInt 		j, CurPhase;
	PetscScalar 	*Par;
	Scaling	 		scal;

	// retrieve some of the required data
	ierr = ScalingCreate(&scal, IOparam->fb, PETSC_FALSE); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************ \n");
    PetscPrintf(PETSC_COMM_WORLD,"|                       COMPUTATION OF THE GRADIENTS                       \n");
    PetscPrintf(PETSC_COMM_WORLD,"| ************************************************************************ \n| ");


   	PetscPrintf(PETSC_COMM_WORLD,"\n| Gradients: \n");
   	
	if (IOparam->FS==0){   
		PetscPrintf(PETSC_COMM_WORLD,"|                    Parameter             |  Gradient (dimensional)  \n");    
		PetscPrintf(PETSC_COMM_WORLD,"|                  -----------------------   ------------------------ \n");    

		VecGetArray(IOparam->P,&Par);
		for(j = 0; j < IOparam->mdN; j++){
			// Get current phase and parameter which is being perturbed
			CurPhase 		= 	IOparam->phs[j];
			strcpy(CurName, IOparam->type_name[j]);	// name
			if (IOparam->par_log10[j]==1){strcpy(logstr, "log10"); }
			else{strcpy(logstr, "     "); }
		

			// Print result
			if (IOparam->FD_gradient[j]>0){
				PetscPrintf(PETSC_COMM_WORLD,"|       FD %5d:   %+5s%+5s[%2i]           %- 1.6e \n",j+1, logstr, CurName, CurPhase, IOparam->grd[j]);
			}
			else{
				PetscPrintf(PETSC_COMM_WORLD,"|  adjoint %5d:   %+5s%+5s[%2i]           %- 1.6e \n",j+1, logstr, CurName, CurPhase, IOparam->grd[j]);
			}

		}
		VecRestoreArray(IOparam->P,&Par);
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD,"|    Computed field-based gradients \n");    
	}
	PetscPrintf(PETSC_COMM_WORLD,"| \n| ");
	

	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr); // because of PETSC_COMM_SELF below
	
	if(IOparam->mdI<_MAX_OBS_ && IOparam->Ap == 1)
	{

        PetscPrintf(PETSC_COMM_WORLD,"\n| Observation points: \n");
		if (IOparam->MfitType == 0) {PetscPrintf(PETSC_COMM_WORLD,"|                                                   Velocity %s       \n",scal.lbl_velocity);    }
        else if (IOparam->MfitType == 1) {PetscPrintf(PETSC_COMM_WORLD,"|                                                    Center values         \n");    }
        PetscPrintf(PETSC_COMM_WORLD,"|                       Location         |      Target         Value     \n");    
        PetscPrintf(PETSC_COMM_WORLD,"|       --------------------------------- --- ------------- ------------- \n");    

		// Print the solution variable at the user defined index (if there are sufficiently few)
		for (j=0; j<IOparam->mdI; j++)
		{

			if(IOparam->Apoint_on_proc[j])
			{
                char vel_com[20];
                PetscScalar x,y,z,CostFunc;

				if (IOparam->MfitType == 0)
				{
					if (IOparam->Av[j] == 1){strcpy(vel_com, "Vx"); }
					if (IOparam->Av[j] == 2){strcpy(vel_com, "Vy"); }
					if (IOparam->Av[j] == 3){strcpy(vel_com, "Vz"); }
					CostFunc = IOparam->Ae[j]*scal.velocity;

					x = IOparam->Ax[j]*scal.length;
					y = IOparam->Ay[j]*scal.length;
					z = IOparam->Az[j]*scal.length;
			
					PetscPrintf(PETSC_COMM_SELF,"| %-4d: [%9.3f; %9.3f; %9.3f]  %s % 8.5e  % 8.5e \n",j+1,x,y,z, vel_com, CostFunc, IOparam->Avel_num[j]);

					// PetscScalar misfit;
					// misfit = 0.5*pow(IOparam->Avel_num[j]-CostFunc, 2);
					//PetscPrintf(PETSC_COMM_SELF,"| Debugging: Velocity misfit =  % 15.8e \n",misfit);
				}
				else if (IOparam->MfitType == 1)
				{
					CostFunc = IOparam->Ae[j];
					
					x = IOparam->Ax[j]*scal.length;
					y = IOparam->Ay[j]*scal.length;
					z = IOparam->Az[j]*scal.length;
					if (!strcmp(IOparam->ObsName[j],"PSD"))
					{
						CostFunc = 	CostFunc*(180/3.14159265359);		// transfer to degrees
					}

					PetscPrintf(PETSC_COMM_SELF,"| %-4d: [%9.3f; %9.3f; %9.3f] %s  % 8.5e  % 8.5e\n",j+1,x,y,z, IOparam->ObsName[j], CostFunc, IOparam->Avel_num[j]);
				}
			}

		}
	  	PetscPrintf(PETSC_COMM_WORLD,"| \n");
 
	}
	
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr); // because of PETSC_COMM_SELF above
	
	PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPointInPro"
/*
    AdjointPointInPro creates a projection vector which projects the requested velocity component(s) 
        from the global solution vector to the observation point(s), assuming linear interpolation
*/
PetscErrorCode AdjointPointInPro(JacRes *jr, AdjGrad *aop, ModParam *IOparam, FreeSurf *surf)
{
	PetscErrorCode      ierr;
	FDSTAG              *fs;
	Vec                 lproX, lproY, lproZ, gproX, gproY, gproZ, pro, xini, lxiniX, lxiniY, lxiniZ, gxiniX, gxiniY, gxiniZ;
	PetscScalar         coord_local[3], *temppro, ***llproX, ***llproY, ***llproZ, *dggproX, *dggproY, *dggproZ, *tempxini, ***llxiniX, ***llxiniY, ***llxiniZ, *dggxiniX, *dggxiniY, *dggxiniZ;
	PetscScalar         *vx, *vy, *vz;
	PetscScalar         f1,f2,f3,f4,f5,f6,f7,f8;
	PetscInt            j, i, ii, sx, sy, sz, nx, ny, nz, I, J, K, II, JJ, KK, lrank, grank, level;
	PetscScalar         w, z, xb, yb, zb, xe, ye, ze, xc, yc, zc, *iter, *ncx, *ncy, *ncz, *ccx, *ccy, *ccz, ***lvx, ***lvy, ***lvz, ***vgrid, ***topo, ***vsurf;
	Discret1D           *dsz;
	InterpFlags         iflag;

	PetscFunctionBegin;

	fs = jr->fs;

	// create vectors with correct layout (doesn't copy values!)
 	ierr = VecDuplicate(jr->gsol, &pro);             CHKERRQ(ierr);
 	ierr = VecDuplicate(jr->gsol, &xini);            CHKERRQ(ierr);
	VecZeroEntries(pro);
	VecZeroEntries(xini);
	

	// Access the local velocities
	ierr = DMDAVecGetArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_X, &gproX);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gproY);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gproZ);   CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lproX);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lproY);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lproZ);   CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_X, &gxiniX);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gxiniY);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gxiniZ);   CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lxiniX);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lxiniY);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lxiniZ);   CHKERRQ(ierr);

	// Zero out entries
	VecZeroEntries(gproX);
	VecZeroEntries(gproY);
	VecZeroEntries(gproZ);
	VecZeroEntries(lproX);
	VecZeroEntries(lproY);
	VecZeroEntries(lproZ);

	VecZeroEntries(gxiniX);
	VecZeroEntries(gxiniY);
	VecZeroEntries(gxiniZ);
	VecZeroEntries(lxiniX);
	VecZeroEntries(lxiniY);
	VecZeroEntries(lxiniZ);

	// If we want only a few observation points we need to interpolate velocity from the grid to the observation point
	if(IOparam->Ap == 1)
	{
		ierr = VecGetArray(aop->vx,&vx); CHKERRQ(ierr);
		ierr = VecGetArray(aop->vy,&vy); CHKERRQ(ierr);
		ierr = VecGetArray(aop->vz,&vz); CHKERRQ(ierr);

		//============================================
		// INDEXING LOOP OVER ALL OBSERVATION POINTS
		//============================================
		for(ii = 0; ii < IOparam->mdI; ii++)
		{
			// Create coordinate vector
			coord_local[0] = IOparam->Ax[ii];
			coord_local[1] = IOparam->Ay[ii];
			coord_local[2] = IOparam->Az[ii];

			// get global & local ranks (processor) of the point
			ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

			// If lrank is not 13 the point is not on this processor
			if(lrank == 13)
			{
				// starting indices & number of cells
				sx = fs->dsx.pstart; nx = fs->dsx.ncels;
				sy = fs->dsy.pstart; ny = fs->dsy.ncels;
				sz = fs->dsz.pstart; nz = fs->dsz.ncels;

				// node & cell coordinates
				ncx = fs->dsx.ncoor; ccx = fs->dsx.ccoor;
				ncy = fs->dsy.ncoor; ccy = fs->dsy.ccoor;
				ncz = fs->dsz.ncoor; ccz = fs->dsz.ccoor;

				// find I, J, K indices by bisection algorithm
				I = FindPointInCellAdjoint(ncx, 0, nx, coord_local[0]);
				J = FindPointInCellAdjoint(ncy, 0, ny, coord_local[1]);
				K = FindPointInCellAdjoint(ncz, 0, nz, coord_local[2]);

				// get coordinates of cell center
				xc = ccx[I];
				yc = ccy[J];
				zc = ccz[K];

				// map marker on the cells of X, Y, Z & center grids
				if(coord_local[0] > xc) { II = I; } else { II = I-1; }
				if(coord_local[1] > yc) { JJ = J; } else { JJ = J-1; }
				if(coord_local[2] > zc) { KK = K; } else { KK = K-1; }

				ierr = DMDAVecGetArray(fs->DA_X, lproX, &llproX);      CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Y, lproY, &llproY);      CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Z, lproZ, &llproZ);      CHKERRQ(ierr);

				if(IOparam->Av[ii] == 1) // Vx velocity
				{
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					xe = (coord_local[0] - ncx[I])/(ncx[I+1] - ncx[I]); xb = 1.0 - xe;
					ye = (coord_local[1] - ccy[JJ])/(ccy[JJ+1] - ccy[JJ]); yb = 1.0 - ye;
					ze = (coord_local[2] - ccz[KK])/(ccz[KK+1] - ccz[KK]); zb = 1.0 - ze;
	
					llproX[sz+KK  ][sy+JJ  ][sx+I  ] = (xb*yb*zb);    
					llproX[sz+KK  ][sy+JJ  ][sx+I+1] = (xe*yb*zb);  
					llproX[sz+KK  ][sy+JJ+1][sx+I  ] = (xb*ye*zb);  
					llproX[sz+KK  ][sy+JJ+1][sx+I+1] = (xe*ye*zb);  
					llproX[sz+KK+1][sy+JJ  ][sx+I  ] = (xb*yb*ze);  
					llproX[sz+KK+1][sy+JJ  ][sx+I+1] = (xe*yb*ze);  
					llproX[sz+KK+1][sy+JJ+1][sx+I  ] = (xb*ye*ze);   
					llproX[sz+KK+1][sy+JJ+1][sx+I+1] = (xe*ye*ze); 

					if(IOparam->OFdef == 1 ) 
					{
						ierr = DMDAVecGetArray(fs->DA_X, lxiniX, &llxiniX);      CHKERRQ(ierr);

						f1 = vx[ii]/lvx[sz+KK  ][sy+JJ  ][sx+I  ];
						f2 = vx[ii]/lvx[sz+KK  ][sy+JJ  ][sx+I+1];
						f3 = vx[ii]/lvx[sz+KK  ][sy+JJ+1][sx+I  ];
						f4 = vx[ii]/lvx[sz+KK  ][sy+JJ+1][sx+I+1];
						f5 = vx[ii]/lvx[sz+KK+1][sy+JJ  ][sx+I  ];
						f6 = vx[ii]/lvx[sz+KK+1][sy+JJ  ][sx+I+1];
						f7 = vx[ii]/lvx[sz+KK+1][sy+JJ+1][sx+I  ];
						f8 = vx[ii]/lvx[sz+KK+1][sy+JJ+1][sx+I+1];

						llxiniX[sz+KK  ][sy+JJ  ][sx+I  ] = IOparam->Ae[ii]/f1;
						llxiniX[sz+KK  ][sy+JJ  ][sx+I+1] = IOparam->Ae[ii]/f2;
						llxiniX[sz+KK  ][sy+JJ+1][sx+I  ] = IOparam->Ae[ii]/f3;
						llxiniX[sz+KK  ][sy+JJ+1][sx+I+1] = IOparam->Ae[ii]/f4;
						llxiniX[sz+KK+1][sy+JJ  ][sx+I  ] = IOparam->Ae[ii]/f5;
						llxiniX[sz+KK+1][sy+JJ  ][sx+I+1] = IOparam->Ae[ii]/f6;
						llxiniX[sz+KK+1][sy+JJ+1][sx+I  ] = IOparam->Ae[ii]/f7;
						llxiniX[sz+KK+1][sy+JJ+1][sx+I+1] = IOparam->Ae[ii]/f8;

						ierr = DMDAVecRestoreArray(fs->DA_X, lxiniX, &llxiniX);      CHKERRQ(ierr);
					}
				}
				else if(IOparam->Av[ii] == 2) // Vy velocity
				{
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					xe = (coord_local[0] - ccx[II])/(ccx[II+1] - ccx[II]); xb = 1.0 - xe;
					ye = (coord_local[1] - ncy[J])/(ncy[J+1] - ncy[J]); yb = 1.0 - ye;
					ze = (coord_local[2] - ccz[KK])/(ccz[KK+1] - ccz[KK]); zb = 1.0 - ze;
	
					llproY[sz+KK  ][sy+J  ][sx+II  ] = (xb*yb*zb);    
					llproY[sz+KK  ][sy+J  ][sx+II+1] = (xe*yb*zb);    
					llproY[sz+KK  ][sy+J+1][sx+II  ] = (xb*ye*zb);  
					llproY[sz+KK  ][sy+J+1][sx+II+1] = (xe*ye*zb);     
					llproY[sz+KK+1][sy+J  ][sx+II  ] = (xb*yb*ze);     
					llproY[sz+KK+1][sy+J  ][sx+II+1] = (xe*yb*ze);  
					llproY[sz+KK+1][sy+J+1][sx+II  ] = (xb*ye*ze);   
					llproY[sz+KK+1][sy+J+1][sx+II+1] = (xe*ye*ze); 

					if(IOparam->OFdef == 1 )
					{
						ierr = DMDAVecGetArray(fs->DA_Y, lxiniY, &llxiniY);      CHKERRQ(ierr);

						f1 = vy[ii]/lvy[sz+KK  ][sy+J  ][sx+II  ];
						f2 = vy[ii]/lvy[sz+KK  ][sy+J  ][sx+II+1];
						f3 = vy[ii]/lvy[sz+KK  ][sy+J+1][sx+II  ];
						f4 = vy[ii]/lvy[sz+KK  ][sy+J+1][sx+II+1];
						f5 = vy[ii]/lvy[sz+KK+1][sy+J  ][sx+II  ];
						f6 = vy[ii]/lvy[sz+KK+1][sy+J  ][sx+II+1];
						f7 = vy[ii]/lvy[sz+KK+1][sy+J+1][sx+II  ];
						f8 = vy[ii]/lvy[sz+KK+1][sy+J+1][sx+II+1];

						llxiniY[sz+KK  ][sy+J  ][sx+II  ] = IOparam->Ae[ii]/f1;
						llxiniY[sz+KK  ][sy+J  ][sx+II+1] = IOparam->Ae[ii]/f2;
						llxiniY[sz+KK  ][sy+J+1][sx+II  ] = IOparam->Ae[ii]/f3;
						llxiniY[sz+KK  ][sy+J+1][sx+II+1] = IOparam->Ae[ii]/f4;
						llxiniY[sz+KK+1][sy+J  ][sx+II  ] = IOparam->Ae[ii]/f5;
						llxiniY[sz+KK+1][sy+J  ][sx+II+1] = IOparam->Ae[ii]/f6;
						llxiniY[sz+KK+1][sy+J+1][sx+II  ] = IOparam->Ae[ii]/f7;
						llxiniY[sz+KK+1][sy+J+1][sx+II+1] = IOparam->Ae[ii]/f8;

						ierr = DMDAVecRestoreArray(fs->DA_Y, lxiniY, &llxiniY);      CHKERRQ(ierr);
					}
				}
				else if(IOparam->Av[ii] == 3) // Vz velocity
				{
					
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					xe = (coord_local[0] - ccx[II])/(ccx[II+1] - ccx[II]); xb = 1.0 - xe;
					ye = (coord_local[1] - ccy[JJ])/(ccy[JJ+1] - ccy[JJ]); yb = 1.0 - ye;
					ze = (coord_local[2] - ncz[K])/(ncz[K+1] - ncz[K]); zb = 1.0 - ze;

					llproZ[sz+K  ][sy+JJ  ][sx+II  ] = (xb*yb*zb);    
					llproZ[sz+K  ][sy+JJ  ][sx+II+1] = (xe*yb*zb);    
					llproZ[sz+K  ][sy+JJ+1][sx+II  ] = (xb*ye*zb);    
					llproZ[sz+K  ][sy+JJ+1][sx+II+1] = (xe*ye*zb);     
					llproZ[sz+K+1][sy+JJ  ][sx+II  ] = (xb*yb*ze);     
					llproZ[sz+K+1][sy+JJ  ][sx+II+1] = (xe*yb*ze);  
					llproZ[sz+K+1][sy+JJ+1][sx+II  ] = (xb*ye*ze);   
					llproZ[sz+K+1][sy+JJ+1][sx+II+1] = (xe*ye*ze); 

					if(IOparam->OFdef == 1 )
					{
						ierr = DMDAVecGetArray(fs->DA_Z, lxiniZ, &llxiniZ);      CHKERRQ(ierr);

						f1 = vz[ii]/lvz[sz+K  ][sy+JJ  ][sx+II  ];
						f2 = vz[ii]/lvz[sz+K  ][sy+JJ  ][sx+II+1];
						f3 = vz[ii]/lvz[sz+K  ][sy+JJ+1][sx+II  ];
						f4 = vz[ii]/lvz[sz+K  ][sy+JJ+1][sx+II+1];
						f5 = vz[ii]/lvz[sz+K+1][sy+JJ  ][sx+II  ];
						f6 = vz[ii]/lvz[sz+K+1][sy+JJ  ][sx+II+1];
						f7 = vz[ii]/lvz[sz+K+1][sy+JJ+1][sx+II  ];
						f8 = vz[ii]/lvz[sz+K+1][sy+JJ+1][sx+II+1];

						llxiniZ[sz+K  ][sy+JJ  ][sx+II  ] = IOparam->Ae[ii]/f1;
						llxiniZ[sz+K  ][sy+JJ  ][sx+II+1] = IOparam->Ae[ii]/f2;
						llxiniZ[sz+K  ][sy+JJ+1][sx+II  ] = IOparam->Ae[ii]/f3;
						llxiniZ[sz+K  ][sy+JJ+1][sx+II+1] = IOparam->Ae[ii]/f4;
						llxiniZ[sz+K+1][sy+JJ  ][sx+II  ] = IOparam->Ae[ii]/f5;
						llxiniZ[sz+K+1][sy+JJ  ][sx+II+1] = IOparam->Ae[ii]/f6;
						llxiniZ[sz+K+1][sy+JJ+1][sx+II  ] = IOparam->Ae[ii]/f7;
						llxiniZ[sz+K+1][sy+JJ+1][sx+II+1] = IOparam->Ae[ii]/f8;

						ierr = DMDAVecRestoreArray(fs->DA_Z, lxiniZ, &llxiniZ);      CHKERRQ(ierr);
					}
				}
				/*
				else
				{
					SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Unknown velocity component ");
				}
				*/

				ierr = DMDAVecRestoreArray(fs->DA_X, lproX, &llproX);      CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Y, lproY, &llproY);      CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Z, lproZ, &llproZ);      CHKERRQ(ierr);
			}

		}

		ierr = VecRestoreArray(aop->vx,&vx);        CHKERRQ(ierr);
		ierr = VecRestoreArray(aop->vy,&vy);        CHKERRQ(ierr);
		ierr = VecRestoreArray(aop->vz,&vz);        CHKERRQ(ierr);
	}

	// We want the whole domain as comparison
	else if (IOparam->Ap == 2)
	{
		for(ii = 0; ii < 3; ii++)
		{
			if(IOparam->Av[ii] == 1)        // vx velocity
			{
				ierr = VecSet(lproX,1);
			}
			else if(IOparam->Av[ii] == 2)  // vy velocity
			{
				ierr = VecSet(lproY,1);
			}
			else if(IOparam->Av[ii] == 3)  // Vz velocity
			{
				ierr = VecSet(lproZ,1);
			}
		}
	}
	else if (IOparam->Ap == 3)     // take the topography velocity as comparison
	{
		for(ii = 0; ii < 3; ii++)
		{
			if (IOparam->Av[ii] == 1)
			{
				dsz   = &fs->dsz;
				level = dsz->rank;
			
				// create column communicator
				ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
			
				// set interpolation flags
				iflag.update    = PETSC_FALSE;
				iflag.use_bound = PETSC_TRUE;
			
				ierr = DMDAVecRestoreArray(fs->DA_X, jr->lvx, &lvx);   CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_X, lproX, &llproX);      CHKERRQ(ierr);
				
				// interpolate velocity component from grid faces to corners
				ierr = InterpXFaceCorner(fs, jr->lvx, jr->lbcor, iflag); CHKERRQ(ierr);
			
				// load ghost values
				LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)
			
				// access topograpy, grid and surface velocity
				ierr = DMDAVecGetArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
			
				// scan all free surface local points
				ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);
			
				START_PLANE_LOOP
				{
					
					// get topography
					z = topo[level][j][i];
			
					// check whether point belongs to domain
					if(z >= dsz->gcrdbeg && z < dsz->gcrdend)
					{
						// find containing cell
						K = FindPointInCellAdjoint(dsz->ncoor, 0, dsz->ncels, z);
			
						// get interpolation weight
						w = (z - dsz->ncoor[K])/(dsz->ncoor[K+1] - dsz->ncoor[K]);
						
						llproX[sz+K][j][i]   = 1.0 - w;
						llproX[sz+K+1][j][i] = w;
					}
				}
				END_PLANE_LOOP
	
				// restore access
				ierr = DMDAVecRestoreArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_X, lproX, &llproX);            CHKERRQ(ierr);
				
			}
			else if (IOparam->Av[ii] == 2)
			{
				dsz   = &fs->dsz;
				level = dsz->rank;
			
				// create column communicator
				ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
			
				// set interpolation flags
				iflag.update    = PETSC_FALSE;
				iflag.use_bound = PETSC_TRUE;
			
				ierr = DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &lvy);   CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Y, lproY, &llproY);      CHKERRQ(ierr);
				
				// interpolate velocity component from grid faces to corners
				ierr = InterpYFaceCorner(fs, jr->lvy, jr->lbcor, iflag); CHKERRQ(ierr);
			
				// load ghost values
				LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)
			
				// access topograpy, grid and surface velocity
				ierr = DMDAVecGetArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
			
				// scan all free surface local points
				ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);
			
				START_PLANE_LOOP
				{
					// get topography
					z = topo[level][j][i];
			
					// check whether point belongs to domain
					if(z >= dsz->gcrdbeg && z < dsz->gcrdend)
					{
						// find containing cell
						K = FindPointInCellAdjoint(dsz->ncoor, 0, dsz->ncels, z);
			
						// get interpolation weight
						w = (z - dsz->ncoor[K])/(dsz->ncoor[K+1] - dsz->ncoor[K]);
						
						llproY[sz+K][j][i]   = 1.0 - w;
						llproY[sz+K+1][j][i] = w;
					}
				}
				END_PLANE_LOOP
	
				// restore access
				ierr = DMDAVecRestoreArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Y, lproY, &llproY);            CHKERRQ(ierr);
			}
			else if (IOparam->Av[ii] == 3)
			{
				dsz   = &fs->dsz;
				level = dsz->rank;
			
				// create column communicator
				ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
			
				// set interpolation flags
				iflag.update    = PETSC_FALSE;
				iflag.use_bound = PETSC_TRUE;
				
				ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &lvz);   CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Z, lproZ, &llproZ);      CHKERRQ(ierr);
				
				// interpolate velocity component from grid faces to corners
				ierr = InterpZFaceCorner(fs, jr->lvz, jr->lbcor, iflag); CHKERRQ(ierr);
			
				// load ghost values
				LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)
			
				// access topograpy, grid and surface velocity
				ierr = DMDAVecGetArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
			
				// scan all free surface local points
				ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);
			
				START_PLANE_LOOP
				{
					// get topography
					z = topo[level][j][i];
			
					// check whether point belongs to domain
					if(z >= dsz->gcrdbeg && z < dsz->gcrdend)
					{
						// find containing cell
						K = FindPointInCellAdjoint(dsz->ncoor, 0, dsz->ncels, z);
			
						// get interpolation weight
						w = (z - dsz->ncoor[K])/(dsz->ncoor[K+1] - dsz->ncoor[K]);
						
						llproZ[sz+K][j][i]   = 1.0 - w;
						llproZ[sz+K+1][j][i] = w;
					}
				}
				END_PLANE_LOOP
	
				// restore access
				ierr = DMDAVecRestoreArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Z, lproZ, &llproZ);            CHKERRQ(ierr);
			}
		}
		
	}

	LOCAL_TO_GLOBAL(fs->DA_X, lproX, gproX);
	LOCAL_TO_GLOBAL(fs->DA_Y, lproY, gproY);
	LOCAL_TO_GLOBAL(fs->DA_Z, lproZ, gproZ);

	LOCAL_TO_GLOBAL(fs->DA_X, lxiniX, gxiniX);
	LOCAL_TO_GLOBAL(fs->DA_Y, lxiniY, gxiniY);
	LOCAL_TO_GLOBAL(fs->DA_Z, lxiniZ, gxiniZ);

	ierr = VecGetArray(gproX, &dggproX);      CHKERRQ(ierr);
	ierr = VecGetArray(gproY, &dggproY);      CHKERRQ(ierr);
	ierr = VecGetArray(gproZ, &dggproZ);      CHKERRQ(ierr);

	ierr = VecGetArray(gxiniX, &dggxiniX);      CHKERRQ(ierr);
	ierr = VecGetArray(gxiniY, &dggxiniY);      CHKERRQ(ierr);
	ierr = VecGetArray(gxiniZ, &dggxiniZ);      CHKERRQ(ierr);

	// Put the proportion into the Projection vector where the user defined the computation coordinates (P)
	ierr = VecGetArray(pro, &temppro);        CHKERRQ(ierr);
	iter = temppro;

	ierr  = PetscMemcpy(iter, dggproX, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, dggproY, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, dggproZ, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	// restore & destroy
	ierr = VecRestoreArray(pro,  &temppro);         CHKERRQ(ierr);

	// Put the proportion into the comparison solution vector where the user defined the computation coordinates (P)
	ierr = VecGetArray(xini, &tempxini);        CHKERRQ(ierr);
	iter = tempxini;

	ierr  = PetscMemcpy(iter, dggxiniX, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, dggxiniY, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, dggxiniZ, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	// restore & destroy
	ierr = VecRestoreArray(xini, &tempxini);         CHKERRQ(ierr);

	ierr = VecRestoreArray(gproX, &dggproX);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gproY, &dggproY);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gproZ, &dggproZ);       CHKERRQ(ierr);

	ierr = VecRestoreArray(gxiniX, &dggxiniX);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gxiniY, &dggxiniY);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gxiniZ, &dggxiniZ);       CHKERRQ(ierr);

	ierr = VecCopy(pro, aop->pro);

	// If it is 0 we want to keep the initially loaded solution and just ignore copying
	if(IOparam->OFdef == 1)
	{
		ierr = VecCopy(xini,IOparam->xini);       CHKERRQ(ierr);        // This is a projection vector, whch projects 
	}

	ierr = VecDestroy(&lproX);
	ierr = VecDestroy(&lproY);
	ierr = VecDestroy(&lproZ);
	ierr = VecDestroy(&gproX);
	ierr = VecDestroy(&gproY);
	ierr = VecDestroy(&gproZ);
	ierr = VecDestroy(&pro);

	ierr = VecDestroy(&lxiniX);
	ierr = VecDestroy(&lxiniY);
	ierr = VecDestroy(&lxiniZ);
	ierr = VecDestroy(&gxiniX);
	ierr = VecDestroy(&gxiniY);
	ierr = VecDestroy(&gxiniZ);
	ierr = VecDestroy(&xini);

	ierr = DMDAVecRestoreArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "AdjointGradientResetParameter"
PetscErrorCode AdjointGradientResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop)
{
	PetscScalar  ini;
	Material_t  *phases;

	PetscFunctionBegin;

	// Get initial value of currently perturbed parameter
	ini    = aop->Ini;
	phases = nl->pc->pm->jr->dbm->phases;

	// Set the current parameter back to its original value
	if (CurPar==_RHO0_)        		{phases[CurPhase].rho    = ini;
	}else if (CurPar==_RHON_)  		{phases[CurPhase].rho_n  = ini;
	}else if (CurPar==_RHOC_)  		{phases[CurPhase].rho_c  = ini;
	}else if (CurPar==_BULK_)  		{phases[CurPhase].Kb     = ini;
	}else if (CurPar==_KP_)  		{phases[CurPhase].Kp     = ini;
	}else if (CurPar==_SHEAR_)  	{phases[CurPhase].G      = ini;
	}else if (CurPar==_ETA_)  		{phases[CurPhase].Bd     = ini;
	}else if (CurPar==_ED_)  		{phases[CurPhase].Ed     = ini;
	}else if (CurPar==_VD_) 		{phases[CurPhase].Vd     = ini;
	}else if (CurPar==_ETA0_) 		{phases[CurPhase].Bn     = ini;
	}else if (CurPar==_N_) 			{phases[CurPhase].n      = ini;   phases[CurPhase].Bn = aop->Ini2;
	}else if (CurPar==_EN_) 		{phases[CurPhase].En     = ini;
	}else if (CurPar==_VN_) 		{phases[CurPhase].Vn     = ini;
	}else if (CurPar==_TAUP_) 		{phases[CurPhase].taup   = ini;
	}else if (CurPar==_GAMMA_) 		{phases[CurPhase].gamma  = ini;
	}else if (CurPar==_Q_) 			{phases[CurPhase].q      = ini;
	}else if (CurPar==_FRICTION_) 	{phases[CurPhase].fr 	 = ini;
	}else if (CurPar==_COHESION_) 	{phases[CurPhase].ch 	 = ini;
	}else if (CurPar==_CP_) 	    {phases[CurPhase].Cp     = ini;
	}else if (CurPar==_A_) 			{phases[CurPhase].A      = ini;}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointFormResidualFieldFD"
PetscErrorCode AdjointFormResidualFieldFD(SNES snes, Vec x, Vec psi, NLSol *nl, AdjGrad *aop, ModParam *IOparam  )
{
    // "geodynamic sensitivity kernels" 
	// ONLY PROGRAMMED FOR DENSITY!!
	// -> Just a debug function for the field gradient ...
	// -> This thing produces the negative of the gradient (multiply with minus; or compare abs value)

	ConstEqCtx  ctx;
	JacRes     *jr;
	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	PetscInt    iter, temprank;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, mcx, mcy, mcz;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz, dx, dy, dz, Le;
	PetscScalar gx, gy, gz, tx, ty, tz, sxx, syy, szz, sxy, sxz, syz, gres;
	PetscScalar J2Inv, DII, z, rho, Tc, pc, pc_lith, pc_pore, dt, fssa, *grav;
	PetscScalar grdt;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc, ***bcp, ***llgradfield;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***T, ***p_lith, ***p_pore;
	Vec         rpl, res;
	char 		CurName[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr = nl->pc->pm->jr;

	// Create stuff
	ierr = VecDuplicate(jr->gres, &res);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &rpl);		 CHKERRQ(ierr);
	
	ierr = VecZeroEntries(jr->lgradfield);	 	 CHKERRQ(ierr);

	strcpy(CurName, IOparam->type_name[0]);	// name

	/*// apply pressure limit at the first visco-plastic timestep and iteration
    if(jr->ts->istep == 1 && jr->ctrl->pLimPlast == PETSC_TRUE)
    {
    	jr->matLim.presLimFlg = PETSC_TRUE;
	}*/

	// access context
	fs  = jr->fs;
	bc  = jr->bc;

	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	mx  = fs->dsx.tnods - 1;
	my  = fs->dsy.tnods - 1;
	mz  = fs->dsz.tnods - 1;

	// access residual context variables
	fssa   =  jr->ctrl.FSSA; // density gradient penalty parameter
	grav   =  jr->ctrl.grav; // gravity acceleration
	dt     =  jr->ts->dt;    // time step

	// Recompute correct strainrates (necessary!!)
	// ierr =  JacResGetEffStrainRate(jr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lgradfield,&llgradfield);      CHKERRQ(ierr);

	for(PetscInt kk=0;kk<fs->dsx.tcels;kk++)
	{
		for(PetscInt jk=0;jk<fs->dsy.tcels;jk++)
		{
			for(PetscInt ik=0;ik<fs->dsz.tcels;ik++)
			{

				// setup constitutive equation evaluation context parameters
				ierr = setUpConstEq(&ctx, jr); CHKERRQ(ierr);

				// clear local residual vectors
				ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
				ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
				ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);
				ierr = VecZeroEntries(jr->gc);  CHKERRQ(ierr);
				temprank = 100;
				aop->Perturb = 1;

				// access work vectors
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->gc,      &gc);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,      &T);      CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx,    &dxx);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy,    &dyy);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz,    &dzz);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy,    &dxy);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz,    &dxz);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz,    &dyz);    CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_X,   jr->lfx,     &fx);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Y,   jr->lfy,     &fy);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Z,   jr->lfz,     &fz);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,     &vx);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,     &vy);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,     &vz);     CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &p_pore); CHKERRQ(ierr);
				ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,     &bcp);    CHKERRQ(ierr);

				//-------------------------------
				// central points
				//-------------------------------
				iter = 0;
				GET_CELL_RANGE(nx, sx, fs->dsx)
				GET_CELL_RANGE(ny, sy, fs->dsy)
				GET_CELL_RANGE(nz, sz, fs->dsz)

				START_STD_LOOP
				{

					// access solution variables
					svCell = &jr->svCell[iter++];

					//=================
					// SECOND INVARIANT
					//=================

					// access strain rates
					XX = dxx[k][j][i];
					YY = dyy[k][j][i];
					ZZ = dzz[k][j][i];

					// x-y plane, i-j indices
					XY1 = dxy[k][j][i];
					XY2 = dxy[k][j+1][i];
					XY3 = dxy[k][j][i+1];
					XY4 = dxy[k][j+1][i+1];

					// x-z plane, i-k indices
					XZ1 = dxz[k][j][i];
					XZ2 = dxz[k+1][j][i];
					XZ3 = dxz[k][j][i+1];
					XZ4 = dxz[k+1][j][i+1];

					// y-z plane, j-k indices
					YZ1 = dyz[k][j][i];
					YZ2 = dyz[k+1][j][i];
					YZ3 = dyz[k][j+1][i];
					YZ4 = dyz[k+1][j+1][i];

					// compute second invariant
					J2Inv = 0.5*(XX*XX + YY*YY + ZZ*ZZ) +
					0.25*(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
					0.25*(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
					0.25*(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

					DII = sqrt(J2Inv);

					//=======================
					// CONSTITUTIVE EQUATIONS
					//=======================

					// access current pressure
					pc = p[k][j][i];

					// current temperature
					Tc = T[k][j][i];

					// access current lithostatic pressure
					pc_lith = p_lith[k][j][i];

					// access current pore pressure (zero if deactivated)
					pc_pore = p_pore[k][j][i];

					// z-coordinate of control volume
					z = COORD_CELL(k, sz, fs->dsz);

					// get characteristic element size
					dx = SIZE_CELL(i, sx, fs->dsx);
					dy = SIZE_CELL(j, sy, fs->dsy);
					dz = SIZE_CELL(k, sz, fs->dsz);
					Le = sqrt(dx*dx + dy*dy + dz*dz);

					// setup control volume parameters
					ierr = setUpCtrlVol(&ctx, svCell->phRat, &svCell->svDev, &svCell->svBulk, pc, pc_lith, pc_pore, Tc, DII, z, Le); CHKERRQ(ierr);

					// evaluate constitutive equations on the cell
					ierr = cellConstEqFD(&ctx, svCell, XX, YY, ZZ, sxx, syy, szz, gres, rho, aop, IOparam,  i,  j,  k,  ik,  jk,  kk); CHKERRQ(ierr);

					// Set perturbation paramter for the finite differences
					if (!strcmp(CurName,"rho"))
					{
						if ((i)==ik && (j)==jk && (k)==kk)
						{
							aop->Perturb = rho*aop->FD_epsilon;
							rho += aop->Perturb;
						}
					}
					
					// compute gravity terms 
					gx = rho*grav[0];
					gy = rho*grav[1];
					gz = rho*grav[2];

					// compute stabilization terms (lumped approximation)
					tx = -fssa*dt*gx;
					ty = -fssa*dt*gy;
					tz = -fssa*dt*gz;

					//=========
					// RESIDUAL
					//=========

					// get mesh steps for the backward and forward derivatives
					bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
					bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
					bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

					// momentum
					fx[k][j][i] -= (sxx + vx[k][j][i]*tx)/bdx + gx/2.0;   fx[k][j][i+1] += (sxx + vx[k][j][i+1]*tx)/fdx - gx/2.0;
					fy[k][j][i] -= (syy + vy[k][j][i]*ty)/bdy + gy/2.0;   fy[k][j+1][i] += (syy + vy[k][j+1][i]*ty)/fdy - gy/2.0;
					fz[k][j][i] -= (szz + vz[k][j][i]*tz)/bdz + gz/2.0;   fz[k+1][j][i] += (szz + vz[k+1][j][i]*tz)/fdz - gz/2.0;

					// pressure boundary constraints
					if(i == 0   && bcp[k][j][i-1] != DBL_MAX) fx[k][j][i]   += -p[k][j][i-1]/bdx;
					if(i == mcx && bcp[k][j][i+1] != DBL_MAX) fx[k][j][i+1] -= -p[k][j][i+1]/fdx;
					if(j == 0   && bcp[k][j-1][i] != DBL_MAX) fy[k][j][i]   += -p[k][j-1][i]/bdy;
					if(j == mcy && bcp[k][j+1][i] != DBL_MAX) fy[k][j+1][i] -= -p[k][j+1][i]/fdy;
					if(k == 0   && bcp[k-1][j][i] != DBL_MAX) fz[k][j][i]   += -p[k-1][j][i]/bdz;
					if(k == mcz && bcp[k+1][j][i] != DBL_MAX) fz[k+1][j][i] -= -p[k+1][j][i]/fdz;
					

					// mass (volume)
					gc[k][j][i] = gres;
					
				}
				END_STD_LOOP

				//-------------------------------
				// xy edge points
				//-------------------------------
				iter = 0;
				GET_NODE_RANGE(nx, sx, fs->dsx)
				GET_NODE_RANGE(ny, sy, fs->dsy)
				GET_CELL_RANGE(nz, sz, fs->dsz)

				START_STD_LOOP
				{
					// access solution variables
					svEdge = &jr->svXYEdge[iter++];

					//=================
					// SECOND INVARIANT
					//=================

					// check index bounds
					I1 = i;   if(I1 == mx) I1--;
					I2 = i-1; if(I2 == -1) I2++;
					J1 = j;   if(J1 == my) J1--;
					J2 = j-1; if(J2 == -1) J2++;

					// access strain rates
					XY = dxy[k][j][i];

					// x-y plane, i-j indices (i & j - bounded)
					XX1 = dxx[k][J1][I1];
					XX2 = dxx[k][J1][I2];
					XX3 = dxx[k][J2][I1];
					XX4 = dxx[k][J2][I2];

					// x-y plane, i-j indices (i & j - bounded)
					YY1 = dyy[k][J1][I1];
					YY2 = dyy[k][J1][I2];
					YY3 = dyy[k][J2][I1];
					YY4 = dyy[k][J2][I2];

					// x-y plane, i-j indices (i & j - bounded)
					ZZ1 = dzz[k][J1][I1];
					ZZ2 = dzz[k][J1][I2];
					ZZ3 = dzz[k][J2][I1];
					ZZ4 = dzz[k][J2][I2];

					// y-z plane j-k indices (j - bounded)
					XZ1 = dxz[k][J1][i];
					XZ2 = dxz[k+1][J1][i];
					XZ3 = dxz[k][J2][i];
					XZ4 = dxz[k+1][J2][i];

					// x-z plane i-k indices (i - bounded)
					YZ1 = dyz[k][j][I1];
					YZ2 = dyz[k+1][j][I1];
					YZ3 = dyz[k][j][I2];
					YZ4 = dyz[k+1][j][I2];

					// compute second invariant
					J2Inv = XY*XY +
					0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
					0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
					0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
					0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
					0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

					DII = sqrt(J2Inv);

					//=======================
					// CONSTITUTIVE EQUATIONS
					//=======================

					// access current pressure (x-y plane, i-j indices)
					pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k][j-1][i] + p[k][j-1][i-1]);

					// current temperature (x-y plane, i-j indices)
					Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k][j-1][i] + T[k][j-1][i-1]);

					// access current lithostatic pressure (x-y plane, i-j indices)
					pc_lith = 0.25*(p_lith[k][j][i] + p_lith[k][j][i-1] + p_lith[k][j-1][i] + p_lith[k][j-1][i-1]);

					// access current pore pressure (x-y plane, i-j indices)
					pc_pore = 0.25*(p_pore[k][j][i] + p_pore[k][j][i-1] + p_pore[k][j-1][i] + p_pore[k][j-1][i-1]);

					// get characteristic element size
					dx = SIZE_NODE(i, sx, fs->dsx);
					dy = SIZE_NODE(j, sy, fs->dsy);
					dz = SIZE_CELL(k, sz, fs->dsz);
					Le = sqrt(dx*dx + dy*dy + dz*dz);

					// setup control volume parameters
					ierr = setUpCtrlVol(&ctx, svEdge->phRat, &svEdge->svDev, NULL, pc, pc_lith, pc_pore, Tc, DII, DBL_MAX, Le); CHKERRQ(ierr);


					// evaluate constitutive equations on the edge
					ierr = edgeConstEqFD(&ctx, svEdge, XY, sxy, aop, IOparam,  i,  j,  k,  ik,  jk,  kk); CHKERRQ(ierr);
					// ierr = edgeConstEq(&ctx, svEdge, XY, sxy); CHKERRQ(ierr);

					//=========
					// RESIDUAL
					//=========

					// get mesh steps for the backward and forward derivatives
					bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
					bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

					// momentum
					fx[k][j-1][i] -= sxy/bdy;   fx[k][j][i] += sxy/fdy;
					fy[k][j][i-1] -= sxy/bdx;   fy[k][j][i] += sxy/fdx;

				}
				END_STD_LOOP

				//-------------------------------
				// xz edge points
				//-------------------------------
				iter = 0;
				GET_NODE_RANGE(nx, sx, fs->dsx)
				GET_CELL_RANGE(ny, sy, fs->dsy)
				GET_NODE_RANGE(nz, sz, fs->dsz)

				START_STD_LOOP
				{
					// access solution variables
					svEdge = &jr->svXZEdge[iter++];

					//=================
					// SECOND INVARIANT
					//=================

					// check index bounds
					I1 = i;   if(I1 == mx) I1--;
					I2 = i-1; if(I2 == -1) I2++;
					K1 = k;   if(K1 == mz) K1--;
					K2 = k-1; if(K2 == -1) K2++;

					// access strain rates
					XZ = dxz[k][j][i];

					// x-z plane, i-k indices (i & k - bounded)
					XX1 = dxx[K1][j][I1];
					XX2 = dxx[K1][j][I2];
					XX3 = dxx[K2][j][I1];
					XX4 = dxx[K2][j][I2];

					// x-z plane, i-k indices (i & k - bounded)
					YY1 = dyy[K1][j][I1];
					YY2 = dyy[K1][j][I2];
					YY3 = dyy[K2][j][I1];
					YY4 = dyy[K2][j][I2];

					// x-z plane, i-k indices (i & k - bounded)
					ZZ1 = dzz[K1][j][I1];
					ZZ2 = dzz[K1][j][I2];
					ZZ3 = dzz[K2][j][I1];
					ZZ4 = dzz[K2][j][I2];

					// y-z plane, j-k indices (k - bounded)
					XY1 = dxy[K1][j][i];
					XY2 = dxy[K1][j+1][i];
					XY3 = dxy[K2][j][i];
					XY4 = dxy[K2][j+1][i];

					// xy plane, i-j indices (i - bounded)
					YZ1 = dyz[k][j][I1];
					YZ2 = dyz[k][j+1][I1];
					YZ3 = dyz[k][j][I2];
					YZ4 = dyz[k][j+1][I2];

					// compute second invariant
					J2Inv = XZ*XZ +
					0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
					0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
					0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
					0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
					0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

					DII = sqrt(J2Inv);

					//=======================
					// CONSTITUTIVE EQUATIONS
					//=======================

					// access current pressure (x-z plane, i-k indices)
					pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k-1][j][i] + p[k-1][j][i-1]);

					// current temperature (x-z plane, i-k indices)
					Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k-1][j][i] + T[k-1][j][i-1]);

					// access current lithostatic pressure (x-z plane, i-k indices)
					pc_lith = 0.25*(p_lith[k][j][i] + p_lith[k][j][i-1] + p_lith[k-1][j][i] + p_lith[k-1][j][i-1]);

					// access current pore pressure (x-z plane, i-k indices)
					pc_pore = 0.25*(p_pore[k][j][i] + p_pore[k][j][i-1] + p_pore[k-1][j][i] + p_pore[k-1][j][i-1]);

					// get characteristic element size
					dx = SIZE_NODE(i, sx, fs->dsx);
					dy = SIZE_CELL(j, sy, fs->dsy);
					dz = SIZE_NODE(k, sz, fs->dsz);
					Le = sqrt(dx*dx + dy*dy + dz*dz);

					// setup control volume parameters
					ierr = setUpCtrlVol(&ctx, svEdge->phRat, &svEdge->svDev, NULL, pc, pc_lith, pc_pore, Tc, DII, DBL_MAX, Le); CHKERRQ(ierr);

					// evaluate constitutive equations on the edge
					ierr = edgeConstEqFD(&ctx, svEdge, XZ, sxz, aop, IOparam,  i,  j,  k,  ik,  jk,  kk); CHKERRQ(ierr);
					// ierr = edgeConstEq(&ctx, svEdge, XZ, sxz); CHKERRQ(ierr);

					//=========
					// RESIDUAL
					//=========

					// get mesh steps for the backward and forward derivatives
					bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
					bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

					// momentum
					fx[k-1][j][i] -= sxz/bdz;   fx[k][j][i] += sxz/fdz;
					fz[k][j][i-1] -= sxz/bdx;   fz[k][j][i] += sxz/fdx;

				}
				END_STD_LOOP

				//-------------------------------
				// yz edge points
				//-------------------------------
				iter = 0;
				GET_CELL_RANGE(nx, sx, fs->dsx)
				GET_NODE_RANGE(ny, sy, fs->dsy)
				GET_NODE_RANGE(nz, sz, fs->dsz)

				START_STD_LOOP
				{
					// access solution variables
					svEdge = &jr->svYZEdge[iter++];

					//=================
					// SECOND INVARIANT
					//=================

					// check index bounds
					J1 = j;   if(J1 == my) J1--;
					J2 = j-1; if(J2 == -1) J2++;
					K1 = k;   if(K1 == mz) K1--;
					K2 = k-1; if(K2 == -1) K2++;

					// access strain rates
					YZ = dyz[k][j][i];

					// y-z plane, j-k indices (j & k - bounded)
					XX1 = dxx[K1][J1][i];
					XX2 = dxx[K1][J2][i];
					XX3 = dxx[K2][J1][i];
					XX4 = dxx[K2][J2][i];

					// y-z plane, j-k indices (j & k - bounded)
					YY1 = dyy[K1][J1][i];
					YY2 = dyy[K1][J2][i];
					YY3 = dyy[K2][J1][i];
					YY4 = dyy[K2][J2][i];

					// y-z plane, j-k indices (j & k - bounded)
					ZZ1 = dzz[K1][J1][i];
					ZZ2 = dzz[K1][J2][i];
					ZZ3 = dzz[K2][J1][i];
					ZZ4 = dzz[K2][J2][i];

					// x-z plane, i-k indices (k -bounded)
					XY1 = dxy[K1][j][i];
					XY2 = dxy[K1][j][i+1];
					XY3 = dxy[K2][j][i];
					XY4 = dxy[K2][j][i+1];

					// x-y plane, i-j indices (j - bounded)
					XZ1 = dxz[k][J1][i];
					XZ2 = dxz[k][J1][i+1];
					XZ3 = dxz[k][J2][i];
					XZ4 = dxz[k][J2][i+1];

					// compute second invariant
					J2Inv = YZ*YZ +
					0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
					0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
					0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
					0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
					0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4);

					DII = sqrt(J2Inv);

					//=======================
					// CONSTITUTIVE EQUATIONS
					//=======================

					// access current pressure (y-z plane, j-k indices)
					pc = 0.25*(p[k][j][i] + p[k][j-1][i] + p[k-1][j][i] + p[k-1][j-1][i]);

					// current temperature (y-z plane, j-k indices)
					Tc = 0.25*(T[k][j][i] + T[k][j-1][i] + T[k-1][j][i] + T[k-1][j-1][i]);

					// access current lithostatic pressure (y-z plane, j-k indices)
					pc_lith = 0.25*(p_lith[k][j][i] + p_lith[k][j-1][i] + p_lith[k-1][j][i] + p_lith[k-1][j-1][i]);

					// access current pore pressure (y-z plane, j-k indices)
					pc_pore = 0.25*(p_pore[k][j][i] + p_pore[k][j-1][i] + p_pore[k-1][j][i] + p_pore[k-1][j-1][i]);

					// get characteristic element size
					dx = SIZE_CELL(i, sx, fs->dsx);
					dy = SIZE_NODE(j, sy, fs->dsy);
					dz = SIZE_NODE(k, sz, fs->dsz);
					Le = sqrt(dx*dx + dy*dy + dz*dz);

					// setup control volume parameters
					ierr = setUpCtrlVol(&ctx, svEdge->phRat, &svEdge->svDev, NULL, pc, pc_lith, pc_pore, Tc, DII, DBL_MAX, Le); CHKERRQ(ierr);

					// evaluate constitutive equations on the edge
					ierr = edgeConstEqFD(&ctx, svEdge, YZ, syz, aop, IOparam,  i,  j,  k,  ik,  jk,  kk); CHKERRQ(ierr);
					// ierr = edgeConstEq(&ctx, svEdge, YZ, syz); CHKERRQ(ierr);

					//=========
					// RESIDUAL
					//=========

					// get mesh steps for the backward and forward derivatives
					bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
					bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

					// update momentum residuals
					fy[k-1][j][i] -= syz/bdz;   fy[k][j][i] += syz/fdz;
					fz[k][j-1][i] -= syz/bdy;   fz[k][j][i] += syz/fdy;

				}
				END_STD_LOOP

				// restore vectors
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,      &gc);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,      &T);      CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx,    &dxx);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy,    &dyy);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz,    &dzz);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy,    &dxy);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz,    &dxz);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz,    &dyz);    CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,     &fx);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,     &fy);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,     &fz);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,     &vx);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,     &vy);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,     &vz);     CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &p_pore); CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,     &bcp);    CHKERRQ(ierr);

				// assemble global residuals from local contributions
				LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
				LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
				LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

				// check convergence of constitutive equations
				ierr = checkConvConstEq(&ctx); CHKERRQ(ierr);

				// copy residuals to global vector
				ierr = JacResCopyRes(jr, rpl); CHKERRQ(ierr);

				ierr = FormResidual(snes, x, res, nl);              CHKERRQ(ierr);

				ierr = VecAYPX(res,-1,rpl);                           CHKERRQ(ierr);
				ierr = VecScale(res,1/aop->Perturb);   CHKERRQ(ierr);

				// Compute the gradient
				ierr = VecDot(res,psi,&grdt);    CHKERRQ(ierr);

				GET_CELL_RANGE(nx, sx, fs->dsx)
				GET_CELL_RANGE(ny, sy, fs->dsy)
				GET_CELL_RANGE(nz, sz, fs->dsz)

				if (ik >= sx && ik < sx+nx && jk >= sy && jk < sy+ny && kk >= sz && kk < sz+nz)
				{
					temprank = 13;
				}
				else
				{
					temprank = 100;
				}

				// If temprank is not 13 the point is not on this processor
				if(temprank == 13)
				{
					llgradfield[kk][jk][ik] = -grdt*aop->CurScal;
				}
			}
		}
	}

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lgradfield,&llgradfield);CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lgradfield);

	// give it back to the adjoint context
	ierr = VecCopy(jr->lgradfield,aop->gradfield);

	// deactivate pressure limit after it has been activated
	// jr->matLim.presLimFlg = PETSC_FALSE;

	ierr = VecDestroy(&rpl);   CHKERRQ(ierr);
	ierr = VecDestroy(&res);   CHKERRQ(ierr);


	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AddMaterialParameterToCommandLineOptions"
PetscErrorCode AddMaterialParameterToCommandLineOptions(char *name, PetscInt ID, PetscScalar val)
{
    PetscErrorCode  ierr;
	char            *option, *option_value;
    PetscBool       PrintOutput=PETSC_FALSE;
    
    PetscFunctionBegin;
    if (ID<0){	asprintf(&option, "-%s ", name); }
	else{ 		asprintf(&option, "-%s[%i]", name, ID); }
	asprintf(&option_value, "%10.20e", val);

	ierr = PetscOptionsSetValue(NULL, option, option_value);    CHKERRQ(ierr);   // this
    
    //PrintOutput = PETSC_TRUE;
    if (PrintOutput){
        PetscPrintf(PETSC_COMM_WORLD,"| **** Added option %s=%s to the database. **** \n",option,option_value);
        PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
    }

    PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------
// Copies a parameter to the LaMEM command-line database
#undef __FUNCT__
#define __FUNCT__ "CopyParameterToLaMEMCommandLine"
PetscErrorCode CopyParameterToLaMEMCommandLine(ModParam *IOparam, PetscScalar CurVal, PetscInt j)
{
    PetscErrorCode  ierr;
    PetscBool       PrintOutput=PETSC_FALSE;
	PetscInt 		CurPhase;
	PetscScalar 	val;
	char 			CurName[_str_len_];
	
    PetscFunctionBegin;

	CurPhase = 	IOparam->phs[j];			// phase of the parameter
    strcpy(CurName, IOparam->type_name[j]);	// name

	ierr = DeleteMaterialParameterFromCommandLineOptions(CurName, CurPhase); 		CHKERRQ(ierr);

	if (IOparam->par_log10[j]==1){
		val = pow(10,CurVal);
		ierr = AddMaterialParameterToCommandLineOptions(CurName, CurPhase, val); 	CHKERRQ(ierr);
	}
	else{
		ierr = AddMaterialParameterToCommandLineOptions(CurName, CurPhase, CurVal); 	CHKERRQ(ierr);
	}
	
	//PrintOutput=PETSC_TRUE;
	if (PrintOutput){
		PetscPrintf(PETSC_COMM_WORLD,"| *** Added parameter %s[%i]=%10.20e to the LaMEM database \n",CurName,CurPhase,CurVal);
	}

    PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DeleteMaterialParameterFromCommandLineOptions"
PetscErrorCode DeleteMaterialParameterFromCommandLineOptions(char *name, PetscInt ID)
{
    PetscErrorCode  ierr;
	char            *option;
    PetscBool       PrintOutput=PETSC_FALSE;
    
    PetscFunctionBegin;
    
    asprintf(&option, "-%s[%i]", name, ID); 
    ierr = PetscOptionsClearValue(NULL, option);    CHKERRQ(ierr);  
    
  //  PrintOutput = PETSC_TRUE;
    if (PrintOutput){
    	PetscPrintf(PETSC_COMM_WORLD,"| **** Deleted option %s from the database. **** \n",option);
       	PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
    }

    PetscFunctionReturn(0);
}


/*---------------------------------------------------------------------------
	This will re-compute the Material options DB, based on the latest command-line 
	parameters and store the result in IOparam->dbm_modified, which is passed to 
	LaMEM @ can be used there 
*/
#undef __FUNCT__
#define __FUNCT__ "CreateModifiedMaterialDatabase"
PetscErrorCode CreateModifiedMaterialDatabase(ModParam *IOparam)
{
    PetscErrorCode  ierr;
	PetscBool       PrintOutput=PETSC_FALSE;
    Scaling         scal;
    FB             *fb;
    PetscFunctionBegin;

    fb          = IOparam->fb;

    // Create scaling object
	ierr = ScalingCreate(&scal, fb, PETSC_FALSE); CHKERRQ(ierr);
	IOparam->dbm_modified.scal    = &scal;

    // Call material database with modified parameters
    ierr = DBMatCreate(&IOparam->dbm_modified, fb, PETSC_FALSE); 	CHKERRQ(ierr);  

    //PrintOutput = PETSC_TRUE;
    if (PrintOutput){
        PetscPrintf(PETSC_COMM_WORLD,"| **** Created Modified material database **** \n");
    }

    PetscFunctionReturn(0);
}

 
 
/*---------------------------------------------------------------------------
	This computes the scaling and sets a default FD method (adjoi nt or not) 
	for a material parameter
*/
#undef __FUNCT__
#define __FUNCT__ "Parameter_SetFDgrad_Option"
PetscErrorCode Parameter_SetFDgrad_Option(PetscInt *FD_grad, char *name)
{
 	PetscFunctionBegin;
	PetscBool PrintOutput=PETSC_FALSE, found=PETSC_FALSE;

	*FD_grad 	= 1;

	// density
	if		 (!strcmp("rho",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("rho_c",name))	{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("beta",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	
	else if  (!strcmp("eta",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("eta0",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("e0",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	
	// diffusion creep
	else if  (!strcmp("Bd",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("Vd",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	else if  (!strcmp("Ed",name))		{ found=PETSC_TRUE;	*FD_grad=0; }
	
	// dislocation creep (power-law)
	else if  (!strcmp("n",name))		{ found=PETSC_TRUE; *FD_grad=0; } 	
	else if  (!strcmp("Bn",name))		{ found=PETSC_TRUE; *FD_grad=0; } 	
	else if  (!strcmp("Vn",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("En",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	
	// Peierls creep
	else if  (!strcmp("Bp",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("Ep",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("Vp",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("taup",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("gamma",name))	{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("q",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	

	// dc-creep
	else if  (!strcmp("Bdc",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("mu",name))		{ found=PETSC_TRUE; *FD_grad=0; }

	// ps-creep
	else if  (!strcmp("Bps",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("d",name))		{ found=PETSC_TRUE; *FD_grad=0; }

	// elasticity
	else if  (!strcmp("G",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("Kb",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	else if  (!strcmp("nu",name))		{ found=PETSC_TRUE; *FD_grad=0; }
	
	// plasticity
	else if  (!strcmp("ch",name))		{ found=PETSC_TRUE; *FD_grad=1; }
	else if  (!strcmp("fr",name))		{ found=PETSC_TRUE; *FD_grad=1; }
	else if  (!strcmp("eta_st",name))	{ found=PETSC_TRUE; *FD_grad=1; }
	
	// temperature
	else if  (!strcmp("alpha",name))	{ found=PETSC_TRUE; *FD_grad=1; }
	else if  (!strcmp("Cp",name))		{ found=PETSC_TRUE; *FD_grad=1; }
	else if  (!strcmp("k",name))		{ found=PETSC_TRUE; *FD_grad=1; }
	else if  (!strcmp("A",name))		{ found=PETSC_TRUE; *FD_grad=1; }
	
	if (!found){
		PetscPrintf(PETSC_COMM_WORLD,"| WARNING: Unknown Adjoint parameter = %s; I therefore use brute-force FD to compute gradients; Please expand Parameter_SetFDgrad_Option in adjoint.cpp \n", name);
	}

	//PrintOutput=PETSC_TRUE;
	if (PrintOutput){
		if (*FD_grad){ 	PetscPrintf(PETSC_COMM_WORLD,"| Parameter %s computed with finite differences \n", name);}
		else{			PetscPrintf(PETSC_COMM_WORLD,"| Parameter %s computed as adjoint \n", name);	}
	}

	
	PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------------
	This prints scaling laws
*/
#undef __FUNCT__
#define __FUNCT__ "PrintScalingLaws"
PetscErrorCode PrintScalingLaws(ModParam *IOparam)
{
 	PetscFunctionBegin;
	PetscErrorCode 	ierr;
	FILE        	*db;
	PetscInt 		j, k=0, CurPhase, idx[IOparam->mdN], maxNum=10;
	PetscScalar 	Exponent[IOparam->mdN], ExpMag[IOparam->mdN], P, grad, *Par, F, A, Vel_check, b;
	char 			CurName[_str_len_], PhaseDescription[_str_len_], logstr[_str_len_], adjointstr[_str_len_], comp_str[_str_len_];
	PetscBool 		isRhoParam=PETSC_FALSE;

	if (!IOparam->ScalLaws){ PetscFunctionReturn(0);}  // do we want to print them?

	ierr 			= 	CreateModifiedMaterialDatabase(IOparam); CHKERRQ(ierr);			// to retrieve phase descriptions
	// Should be adjoint Gradients
	if (!(IOparam->use == _adjointgradients_ )) { PetscFunctionReturn(0);}  // do we want to print them?

	// Should be vs. Cost function
	if (IOparam->Gr==0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, " Scaling laws require: Adjoint_GradientCalculation = Solution");
	}

	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| -------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"| Scaling laws: \n");
	PetscPrintf(PETSC_COMM_WORLD,"| \n");
	PetscPrintf(PETSC_COMM_WORLD,"|   Assumption: \n");
	PetscPrintf(PETSC_COMM_WORLD,"|              Vel = A * p[0]^b[0] * p[1]^b[1] * p[2]^b[2] * ... \n");
	PetscPrintf(PETSC_COMM_WORLD,"|                where \n");
	PetscPrintf(PETSC_COMM_WORLD,"|                      Vel - velocity;   p[] - parameter \n");
	PetscPrintf(PETSC_COMM_WORLD,"|                      b[] - exponent;   A   - prefactor (computed for p>0) \n");
	PetscPrintf(PETSC_COMM_WORLD,"|       \n");

	if (IOparam->mdN<10){	PetscPrintf(PETSC_COMM_WORLD,"|   Results: \n"); maxNum = IOparam->mdN;	}
	else{					PetscPrintf(PETSC_COMM_WORLD,"|   Results (10 most important): \n"); 	}
	/* 	Compute exponents of scaling laws
	 	We assume that velocity is given by
			Vz = A * p[0]^b[0] * p[1]^b[1] * p[2]^b[2] * ..
		with 
				p[0] - parameters	
				b[0] - exponent
				A    - prefactor of scaling law
	*/
	F = IOparam->mfit;		// velocity value
	A = F;
	VecGetArray(IOparam->P,&Par);
	for(j = 0; j < IOparam->mdN; j++)
	{
		grad 		= 	IOparam->grd[j];					// gradient
		strcpy(CurName, IOparam->type_name[j]);	// name
		P    		= 	Par[j];								// parameter value	
		if (!strcmp("rho",CurName))
		{
			P = P - IOparam->ReferenceDensity;				// Compute with density difference
			isRhoParam=PETSC_TRUE;
		}	

		// Compute exponent
		Exponent[j] = 	grad*P/F;							// b value
		
		if (PetscIsInfOrNanScalar(Exponent[j]))
		{
			ExpMag[j] 	=	0; 			
		}
		else
		{
			ExpMag[j] 	=	-PetscAbs(Exponent[j]); 			// magnitude of exponent (for sorting later)
		}

		idx[j] 		=	j;
		if (P>0)
		{ // only non-zero positive parameters contribute
			if (!isnan(grad))
			{
				A 	=   A*1.0/(PetscPowScalar(P,Exponent[j]));	// prefactor
			}
			else
			{
				PetscPrintf(PETSC_COMM_WORLD,"Detected nan for parameter %s A=%f\n", CurName,A);
			}
		}
	}
	VecRestoreArray(IOparam->P,&Par);

	// Sort according to magnitude
	PetscSortRealWithPermutation(IOparam->mdN,ExpMag,idx);

	
	PetscPrintf(PETSC_COMM_WORLD,"|            Parameter      |    Exponent b[]  |  Phase Description    \n");
	PetscPrintf(PETSC_COMM_WORLD,"|     ----------------------  -----------------  ----------------------- \n");
	for(j = 0; j < maxNum; j++)
	{
		k 				= idx[j];
		CurPhase 		= 	IOparam->phs[k];
		strcpy(CurName, IOparam->type_name[k]);	// name
		if (!strcmp("rho",CurName) & (IOparam->ReferenceDensity!=0.0))
		{
			char *Name;	
			asprintf(&Name, "delta(%s)", CurName);	// w compute w.r.t. Reference Density
			strcpy(CurName, Name);	// name
		}
		if (IOparam->par_log10[k]==1){strcpy(logstr, "log10"); }
		else{strcpy(logstr, "     "); }

		strcpy(PhaseDescription, IOparam->dbm_modified.phases[CurPhase].Name);	// name
		if (!strlen(PhaseDescription)){strcpy(PhaseDescription, "-");} 			// if no name is indicated in input file	
		if (CurPhase<0)
		{
			PetscPrintf(PETSC_COMM_WORLD,"|      %-5s%10s       %- 8.3f          %s\n",logstr, CurName, Exponent[k],PhaseDescription);		
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD,"|      %-5s%10s[%3i]      %- 8.3f         %s\n",logstr, CurName, CurPhase, Exponent[k],PhaseDescription);
		}
		
	}
	PetscPrintf(PETSC_COMM_WORLD,"|       \n");
	PetscPrintf(PETSC_COMM_WORLD,"|   Prefactor A               : %- 2.8e \n",A);
	if (isRhoParam)
	{
		PetscPrintf(PETSC_COMM_WORLD,"|   Reference Density         : %- 2.2f  \n",IOparam->ReferenceDensity);
	}

	// Check that velocity is indeed computed
	Vel_check = A;
	VecGetArray(IOparam->P,&Par);
	for(j = 0; j < IOparam->mdN; j++)
	{
		P    		= 	Par[j];	
		strcpy(CurName, IOparam->type_name[j]);	// name
		if (!strcmp("rho",CurName) )
		{
			P = P - IOparam->ReferenceDensity;				// Compute with density difference
		}							
		b 			=	Exponent[j];
		if (P>0)
		{ 
			Vel_check	*=  PetscPowScalar(P,b);
		}
	}
	VecRestoreArray(IOparam->P,&Par);
	PetscPrintf(PETSC_COMM_WORLD,"|   Velocity check            : %- 2.8e \n",Vel_check);
		

	// Save output to file
	
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		char filename[_str_len_];
		strcpy(filename, "ScalingLaw.dat");	// name
		if (IOparam->ScalLawFilename)
		{
			PetscMemcpy(filename, IOparam->ScalLawFilename,   (size_t)_str_len_*sizeof(char) ); 		
		}

		db = fopen(filename, "w");

		fprintf(db,"# Scaling Law, computed on %s %s  \n",__DATE__,__TIME__);
		fprintf(db,"#  \n");
		fprintf(db,"#   Vel = A * p[0]^b[0] * p[1]^b[1] * p[2]^b[2] * ... \n");
		fprintf(db,"#  \n");
		fprintf(db,"# Prefactor A       : %- 10.8e  \n",A);
		fprintf(db,"# Reference Density : %- 10.8f  \n",IOparam->ReferenceDensity);
		fprintf(db,"#  \n");
		fprintf(db,"# Observation points:  \n");
		
		fprintf(db,"#     x              y              z               Component   Measured value   \n");
		fprintf(db,"# --- -------------- -------------- --------------  ----------  --------------- \n");
		
		// Print observation points info
		for(j = 0; j < IOparam->mdI; j++)
		{
			if 		(IOparam->Av[j]==1){strcpy(comp_str, "Vx");	}
			else if (IOparam->Av[j]==2){strcpy(comp_str, "Vy");	}
			else if (IOparam->Av[j]==3){strcpy(comp_str, "Vz");	}
		
			fprintf(db,"# %3i %- 14.5f %- 14.5f %- 14.5f   %s         %- 14.5e\n",j+1, IOparam->Ax[j], IOparam->Ay[j], IOparam->Az[j], comp_str, IOparam->Avel_num[j]);
		}
		fprintf(db,"#  \n");
		fprintf(db,"#  \n");
		


		// Print Scaling laws info
		fprintf(db,"# Scaling law parameters:\n");
		fprintf(db,"# Parameter             Phase    Exponent b[]       Value p[]          Type     Phase Description   \n");
		fprintf(db,"# --------------------  -------  -----------------  -----------------  -------  --------------------\n");
		VecGetArray(IOparam->P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			k 				= 	idx[j];	
			CurPhase 		= 	IOparam->phs[k];
			strcpy(CurName, IOparam->type_name[k]);	// name
			strcpy(PhaseDescription, IOparam->dbm_modified.phases[CurPhase].Name);	// name

			P = Par[k];
			if (!strcmp("rho",CurName) & (IOparam->ReferenceDensity!=0.0))
			{
				char *Name;	
				P = P - IOparam->ReferenceDensity;				// Compute with density difference
				asprintf(&Name, "delta(%s)", CurName);			// clarify that we compute w.r.t. Reference Density
				strcpy(CurName, Name);							// name
			}	

			if (IOparam->par_log10[k]==1){strcpy(logstr, "log10"); }
			else{strcpy(logstr, "     "); }
			if (IOparam->FD_gradient[k]==0){strcpy(adjointstr, "adjoint"); }
			else{strcpy(adjointstr, "FD     "); }
		

			fprintf(db,"  %s %13s    %3i     %- 18.9e %- 18.9e %s  %s \n",logstr, CurName, CurPhase, Exponent[k],P, adjointstr, PhaseDescription);

		}
		VecRestoreArray(IOparam->P,&Par);
		fclose(db);
		PetscPrintf(PETSC_COMM_WORLD,"|   Scaling law data saved to :  %s \n",IOparam->ScalLawFilename);
		
	}
	PetscPrintf(PETSC_COMM_WORLD,"|       \n");
	PetscPrintf(PETSC_COMM_WORLD,"| -------------------------------------------------------------------------\n");



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointGet_F_dFdu_Center"
PetscErrorCode AdjointGet_F_dFdu_Center(JacRes *jr, AdjGrad *aop, ModParam *IOparam)
{
	// Compute derivative of stress objective function with respect to the solution (dF/du) (dF/dst = (P*st-P*st_ini) * dphi/de * de/du)       
	// dphi/de = (1/(2*pow(exx-eyy,2)) * (-2*exy,2exy,exx-eyy)); e = deviatoric strainrate

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    ii, i, j, k, nx, ny, nz, sx, sy, sz, iterat, lrank, grank;
	PetscInt    I, J, K, ID;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar XX, YY, ZZ, XY, XZ, YZ, XY2, XZ2, YZ2, E2;
	PetscScalar bdx, bdy, bdz, fdx, fdy, fdz, dx, dy, dz;
	PetscScalar phival, Parameter, Param_local, mfitParam;
	PetscScalar *tempPar,  *tempdPardu;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***vx,  ***vy,  ***vz;
	Vec         gxPar, gyPar, gzPar, gxdPardu, gydPardu, gzdPardu;
	Vec         lxPar, lyPar, lzPar, lxdPardu, lydPardu, lzdPardu;
	PetscScalar *dggxPar, *dggyPar, *dggzPar, *dggxdPardu, *dggydPardu, *dggzdPardu, *iter, *sty;
	PetscScalar ***xPar, ***yPar, ***zPar, ***xdPardu, ***ydPardu, ***zdPardu;
	PetscScalar dPardu_local;
	PetscScalar coord_local[3],  Cons;
	PetscInt    As_Ind[IOparam->mdI+1];
	Scaling 	*scal;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs 		= jr->fs;
	scal 	= jr->scal;

	// Initialize vector to store observations:
	ierr = VecZeroEntries(aop->sty);         CHKERRQ(ierr);

	// Initialize the cost function
	IOparam->mfitCenter = 	0.0;
	mfitParam 			= 	0.0;
	Param_local 		= 	0.0;
	Parameter 			=	0.0;

	/* 
		For the cost function,  determine in which FDSTAG cell the observation are made.
		Note: we compute the parameters only at the center of the FDSTAG cell & interpolate
		values such as Exy from edges->center
	*/
	for(ii = 0; ii < IOparam->mdI; ii++)
	{
		// Create coordinate vector
		coord_local[0] = IOparam->Ax[ii];
		coord_local[1] = IOparam->Ay[ii];
		coord_local[2] = IOparam->Az[ii];

		ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

		if(lrank == 13) // point is located on the current processor
		{
			// starting indices & number of cells
			nx = fs->dsx.ncels;
			ny = fs->dsy.ncels;
			nz = fs->dsz.ncels;

			// node & cell coordinates
			ncx = fs->dsx.ncoor;
			ncy = fs->dsy.ncoor;
			ncz = fs->dsz.ncoor;

			// find I, J, K indices by bisection algorithm
			I = FindPointInCellAdjoint(ncx, 0, nx, coord_local[0]);
			J = FindPointInCellAdjoint(ncy, 0, ny, coord_local[1]);
			K = FindPointInCellAdjoint(ncz, 0, nz, coord_local[2]);

			GET_CELL_ID(ID, I, J, K, nx, ny);

			As_Ind[ii] = ID;

			/*
			// Compute the local finite differences for dF/du (debugging)
			PetscScalar h;
			// h = 1e-6;
			h = 0.0;
			ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);
			vx[K][J][I] +=h;
			PetscPrintf(PETSC_COMM_SELF,"%d %d %d\n\n",I,J,K);
			*/

			// PetscPrintf(PETSC_COMM_WORLD,"DEBUG 1a %.10f %.10f %.10f %d, %d\n\n",I,J,K,As_Ind[ii],ii,ncx[I],ncy[J],ncz[K] );

		}
	}

	// clear local residual vectors
	ierr = VecZeroEntries(jr->phi);         CHKERRQ(ierr);
	ierr = VecZeroEntries(aop->dPardu); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_X, &gxPar);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gyPar);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gzPar);   CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lxPar);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lyPar);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lzPar);   CHKERRQ(ierr);

	VecZeroEntries(gxPar);
	VecZeroEntries(gyPar);
	VecZeroEntries(gzPar);

	VecZeroEntries(lxPar);
	VecZeroEntries(lyPar);
	VecZeroEntries(lzPar);

	ierr = DMCreateGlobalVector(fs->DA_X, &gxdPardu);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gydPardu);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gzdPardu);   CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lxdPardu);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lydPardu);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lzdPardu);   CHKERRQ(ierr);

	VecZeroEntries(gxdPardu);
	VecZeroEntries(gydPardu);
	VecZeroEntries(gzdPardu);

	VecZeroEntries(lxdPardu);
	VecZeroEntries(lydPardu);
	VecZeroEntries(lzdPardu);

	// Recompute correct strainrates 
	// ierr =  JacResGetEffStrainRate(jr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx,    &dxx);    	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy,    &dyy);    	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz,    &dzz);    	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy,    &dxy);    	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz,    &dxz);    	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz,    &dyz);    	CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,     &vx); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,     &vy); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,     &vz); 		CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_X, lxPar,    &xPar);    		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, lyPar,    &yPar);    		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, lzPar,    &zPar);    		CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_X, lxdPardu, &xdPardu); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, lydPardu, &ydPardu); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, lzdPardu, &zdPardu); 		CHKERRQ(ierr);

	ierr = VecGetArray(aop->sty,&sty); CHKERRQ(ierr);

	iterat = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);
		
		svCell = &jr->svCell[iterat++];
		svBulk = &svCell->svBulk;

		// access strain rates for the current cell
		XX = dxx[k][j][i];
		YY = dyy[k][j][i];
		ZZ = dzz[k][j][i];
		
		// x-y plane, i-j indices
		XY = (dxy[k][j][i] + dxy[k][j][i+1] + dxy[k][j+1][i] + dxy[k][j+1][i+1])/4.0;

		// x-z plane
		XZ = (dxz[k][j][i] + dxz[k][j][i+1] + dxz[k+1][j][i] + dxz[k+1][j][i+1])/4.0;

		// y-z plane
		YZ = (dyz[k][j][i] + dyz[k][j+1][i] + dyz[k+1][j][i] + dyz[k+1][j+1][i])/4.0;

		// 2nd invariant: first square shear components, before projecting them to the center!
		XY2 = (pow(dxy[k][j][i],2.0) + pow(dxy[k][j][i+1],2.0) + pow(dxy[k][j+1][i],2.0) + pow(dxy[k][j+1][i+1],2.0))/4.0;
		XZ2 = (pow(dxz[k][j][i],2.0) + pow(dxz[k][j][i+1],2.0) + pow(dxz[k+1][j][i],2.0) + pow(dxz[k+1][j][i+1],2.0))/4.0;
		YZ2 = (pow(dyz[k][j][i],2.0) + pow(dyz[k][j+1][i],2.0) + pow(dyz[k+1][j][i],2.0) + pow(dyz[k+1][j+1][i],2.0))/4.0;

		E2  =  pow( 0.5*(XX*XX + YY*YY + ZZ*ZZ) + XY2 + XZ2 + YZ2, 0.5);		// second invariant @ center
		
		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// if PSD is chosen we first have to calculate it at every point  TAKE CARE HERE THIS NEEDS TO BE CHANGED TO BE MORE FLEXIBLE -> ONLY CHECKS FIRST OBS
		if (!strcmp(IOparam->ObsName[0],"PSD"))
		{			
			// Perform the adjoint computation for PSD:
			phival 	= 0.5*atan(((2*XY)/(XX-YY)));
			// PetscPrintf(PETSC_COMM_WORLD,"%.20f %.20f %.20f\n   ",XY,XX,YY);
			if(phival<0) 
			{
				phival = -phival;
			} 
			else
			{
				phival = ((3.14159265359/4 - phival) + 3.14159265359/4);    // minu could be set to -1 here to mimic the discontiunuity in the ambuigity of phival in the gradient
			} 
			if(XY>0) phival = phival + 3.14159265359/2;
			svBulk->phi  	= phival * (180/3.14159265359);
		}

		// Loop over observations
		for(ii = 0; ii < IOparam->mdI; ii++)
		{
			if (As_Ind[ii] == iterat-1)
			{
				// Create coordinate vector
				coord_local[0] = IOparam->Ax[ii];
				coord_local[1] = IOparam->Ay[ii];
				coord_local[2] = IOparam->Az[ii];

				ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);
				if(lrank == 13)
				{
					/* Retrieve parameter */
					if (!strcmp(IOparam->ObsName[ii],"PSD")) {Parameter       = phival;}
					// Perform adjoint computation for strainrate components
					else if (!strcmp(IOparam->ObsName[ii],"Exx")) {	Parameter =	XX;		}
					else if (!strcmp(IOparam->ObsName[ii],"Eyy")) {	Parameter =	YY;		}
					else if (!strcmp(IOparam->ObsName[ii],"Ezz")) {	Parameter =	ZZ;		}
					else if (!strcmp(IOparam->ObsName[ii],"Exy")) {	Parameter =	XY;		}
					else if (!strcmp(IOparam->ObsName[ii],"Exz")) {	Parameter =	XZ;		}
					else if (!strcmp(IOparam->ObsName[ii],"Eyz")) {	Parameter =	YZ;		}
					else if (!strcmp(IOparam->ObsName[ii],"E2nd")){	Parameter =	E2;		}


					if (!strcmp(IOparam->ObsName[ii],"PSD"))
					{					
						/* If we perform the computation for PSD (Principal Stress Direction) */	
						// Compute objective function derivative (dFdu = P*(st-st_ini))
						dPardu_local = 0;
						if (IOparam->Gr == 0)
						{
							// dphidu_local = svBulk->phi-IOparam->Ae[ii];
							dPardu_local = phival-IOparam->Ae[ii]; 

							// Compute objective function value (F += (1/2)*[P*(st-st_ini)' * P*(st-st_ini)])
							mfitParam += (phival-IOparam->Ae[ii])*(phival-IOparam->Ae[ii]);
						}
						else if (IOparam->Gr == 1)
						{
							dPardu_local = 1;    // -> use this for sens Kernel

							// Just give back norm of solution at the points
							mfitParam += (phival);
						}

						Cons = -1 * (1/(pow(XX-YY,2)+4*pow(XY,2)));   // See up there that if phival < 0 -> dphi = - & if phival > 0 -> dphi = - as well
						xdPardu[k][j  ][i  ] += dPardu_local * (Cons*(-XY)*(-1.0/dx) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*( 1.0/bdy-1.0/fdy)*(1.0/8.0));	// 1
						xdPardu[k][j  ][i+1] += dPardu_local * (Cons*(-XY)*( 1.0/dx) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*( 1.0/bdy-1.0/fdy)*(1.0/8.0));	// 2
						xdPardu[k][j-1][i  ] += dPardu_local * (Cons*(-XY)*( 0     ) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*(-1.0/bdy      )*(1.0/8.0));   // 5
						xdPardu[k][j-1][i+1] += dPardu_local * (Cons*(-XY)*( 0     ) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*(-1.0/bdy      )*(1.0/8.0));	// 6
						xdPardu[k][j+1][i  ] += dPardu_local * (Cons*(-XY)*( 0     ) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*( 1.0/fdy      )*(1.0/8.0));	// 7
						xdPardu[k][j+1][i+1] += dPardu_local * (Cons*(-XY)*( 0     ) + Cons*( XY)*( 0   ) + Cons*(XX-YY)*( 1.0/fdy      )*(1.0/8.0));	// 8

						ydPardu[k][j  ][i  ] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*(-1.0/dy) + Cons*(XX-YY)*( 1.0/bdx-1.0/fdx)*(1.0/8.0));	// 3
						ydPardu[k][j+1][i  ] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*( 1.0/dy) + Cons*(XX-YY)*( 1.0/bdx-1.0/fdx)*(1.0/8.0));	// 4
						ydPardu[k][j  ][i-1] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*( 0     ) + Cons*(XX-YY)*(-1.0/bdx        )*(1.0/8.0));	// 9
						ydPardu[k][j+1][i-1] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*( 0     ) + Cons*(XX-YY)*(-1.0/bdx        )*(1.0/8.0));	// 10
						ydPardu[k][j  ][i+1] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*( 0     ) + Cons*(XX-YY)*( 1.0/fdx        )*(1.0/8.0));	// 11
						ydPardu[k][j+1][i+1] += dPardu_local * (Cons*(-XY)*( 0   ) + Cons*( XY)*( 0     ) + Cons*(XX-YY)*( 1.0/fdx        )*(1.0/8.0));	// 12

						// Store observation
						sty[ii] = phival * (180/3.14159265359);
					}
					else if (!strcmp(IOparam->ObsName[ii],"Exx") | !strcmp(IOparam->ObsName[ii],"Eyy") | !strcmp(IOparam->ObsName[ii],"Ezz") | 
							 !strcmp(IOparam->ObsName[ii],"Exy") | !strcmp(IOparam->ObsName[ii],"Eyz") | !strcmp(IOparam->ObsName[ii],"Exz") | !strcmp(IOparam->ObsName[ii],"E2nd"))
					{		

						// Perform computation for strainrate tensor components 
						dPardu_local = 0;
						if (IOparam->Gr == 0)
						{
							// dphidu_local = svBulk->phi-IOparam->Ae[ii];
							Param_local = Parameter-IOparam->Ae[ii]; 

							// Compute objective function value (F += (1/2)*[P*(st-st_ini)' * P*(st-st_ini)])
							mfitParam += (Parameter-IOparam->Ae[ii])*(Parameter-IOparam->Ae[ii]);
						}
						else if (IOparam->Gr == 1)
						{
							Param_local = 1;    // -> use this for sens Kernel

							// Just give back norm of solution at the points
							mfitParam += (Parameter);
						}

						if 		(!strcmp(IOparam->ObsName[ii],"Exx"))
						{
							// derivative of Exx vs Vx:
							xdPardu[k][j  ][i  ] += Param_local * (-1.0/dx);	// 1
							xdPardu[k][j  ][i+1] += Param_local * ( 1.0/dx);	// 2
						}
						else if (!strcmp(IOparam->ObsName[ii],"Eyy"))
						{
							// derivative of Eyy vs Vy:
							ydPardu[k][j  ][i  ] += Param_local * (-1.0/dy) ;	// 3
							ydPardu[k][j+1][i  ] += Param_local * ( 1.0/dy) ;	// 4
						}
						else if (!strcmp(IOparam->ObsName[ii],"Ezz"))
						{
							// derivative of Ezz vs Vz:
							zdPardu[k  ][j  ][i  ] += Param_local * (-1.0/dz) ;	
							zdPardu[k+1][j  ][i  ] += Param_local * ( 1.0/dz) ;	
						}
						else if (!strcmp(IOparam->ObsName[ii],"Exy"))
						{
							// derivative of Exy (@ lower-left corner) vs Vx & Vy:
							xdPardu[k][j  ][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdy);	// 1
							xdPardu[k][j-1][i  ] += Param_local * (1.0/8.0) * (-1.0/bdy);	// 5
							ydPardu[k][j  ][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdx);	// 3
							ydPardu[k][j  ][i-1] += Param_local * (1.0/8.0) * (-1.0/bdx);	// 9

							// derivative of Exy (@ lower-right corner) vs Vx & Vy:
							xdPardu[k][j  ][i+1] += Param_local * (1.0/8.0) * ( 1.0/bdy);	// 2
							xdPardu[k][j-1][i+1] += Param_local * (1.0/8.0) * (-1.0/bdy);	// 6
							ydPardu[k][j  ][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdx);	// 11
							ydPardu[k][j  ][i  ] += Param_local * (1.0/8.0) * (-1.0/fdx);	// 3

							// derivative of Exy (@ upper-right corner) vs Vx & Vy:
							xdPardu[k][j+1][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdy);	// 8
							xdPardu[k][j  ][i+1] += Param_local * (1.0/8.0) * (-1.0/fdy);	// 2
							ydPardu[k][j+1][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdx);	// 12
							ydPardu[k][j+1][i  ] += Param_local * (1.0/8.0) * (-1.0/fdx);	// 3

							// derivative of Exy (@ upper-left corner) vs Vx & Vy:
							xdPardu[k][j+1][i  ] += Param_local * (1.0/8.0) * ( 1.0/fdy);	// 7
							xdPardu[k][j  ][i  ] += Param_local * (1.0/8.0) * (-1.0/fdy);	// 1
							ydPardu[k][j+1][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdx);	// 4
							ydPardu[k][j+1][i-1] += Param_local * (1.0/8.0) * (-1.0/bdx);	// 10
						}
						else if (!strcmp(IOparam->ObsName[ii],"Exz"))
						{
							// derivative of Exz (@ bottom-left corner) vs Vx & Vz:
							xdPardu[k  ][j][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdz);	// 
							xdPardu[k-1][j][i  ] += Param_local * (1.0/8.0) * (-1.0/bdz);	// 
							zdPardu[k  ][j][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdx);	// 
							zdPardu[k  ][j][i-1] += Param_local * (1.0/8.0) * (-1.0/bdx);	// 

							// derivative of Exz (@ bottom-right corner) vs Vx & Vz:
							xdPardu[k  ][j][i+1] += Param_local * (1.0/8.0) * ( 1.0/bdz);	// 
							xdPardu[k-1][j][i+1] += Param_local * (1.0/8.0) * (-1.0/bdz);	// 
							zdPardu[k  ][j][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdx);	// 
							zdPardu[k  ][j][i  ] += Param_local * (1.0/8.0) * (-1.0/fdx);	// 

							// derivative of Exz (@ upper-left corner) vs Vx & Vz:
						  	xdPardu[k+1][j][i  ] += Param_local * (1.0/8.0) * ( 1.0/fdz);	// 
							xdPardu[k  ][j][i  ] += Param_local * (1.0/8.0) * (-1.0/fdz);	// 
							zdPardu[k+1][j][i  ] += Param_local * (1.0/8.0) * ( 1.0/bdx);	// 
							zdPardu[k+1][j][i-1] += Param_local * (1.0/8.0) * (-1.0/bdx);	// 

							// derivative of Exz (@ upper-right corner) vs Vx & Vz:
							xdPardu[k+1][j][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdz);	// 
							xdPardu[k  ][j][i+1] += Param_local * (1.0/8.0) * (-1.0/fdz);	// 
							zdPardu[k+1][j][i+1] += Param_local * (1.0/8.0) * ( 1.0/fdx);	// 
							zdPardu[k+1][j][i  ] += Param_local * (1.0/8.0) * (-1.0/fdx);	// 
						}
						else if (!strcmp(IOparam->ObsName[ii],"Eyz"))
						{
							// derivative of Eyz (@ bottom-front corner) vs Vy & Vz:
							ydPardu[k  ][j  ][i] += Param_local * (1.0/8.0) * ( 1.0/bdz);	 
							ydPardu[k-1][j  ][i] += Param_local * (1.0/8.0) * (-1.0/bdz);	 
							zdPardu[k  ][j  ][i] += Param_local * (1.0/8.0) * ( 1.0/bdy);	 
							zdPardu[k  ][j-1][i] += Param_local * (1.0/8.0) * (-1.0/bdy);	 

							// derivative of Eyz (@ bottom-back corner) vs Vy & Vz:
							ydPardu[k  ][j+1][i] += Param_local * (1.0/8.0) * ( 1.0/bdz);	 
							ydPardu[k-1][j+1][i] += Param_local * (1.0/8.0) * (-1.0/bdz);	 
							zdPardu[k  ][j+1][i] += Param_local * (1.0/8.0) * ( 1.0/fdy);	 
							zdPardu[k  ][j  ][i] += Param_local * (1.0/8.0) * (-1.0/fdy);	 

							// derivative of Eyz (@ top-front corner) vs Vy & Vz:
							ydPardu[k+1][j  ][i] += Param_local * (1.0/8.0) * ( 1.0/fdz);	 
							ydPardu[k  ][j  ][i] += Param_local * (1.0/8.0) * (-1.0/fdz);	 
							zdPardu[k+1][j  ][i] += Param_local * (1.0/8.0) * ( 1.0/bdy);	 
							zdPardu[k+1][j-1][i] += Param_local * (1.0/8.0) * (-1.0/bdy);	 

							// derivative of Eyz (@ top-back corner) vs Vy & Vz:
							ydPardu[k+1][j+1][i] += Param_local * (1.0/8.0) * ( 1.0/fdz);	 
							ydPardu[k  ][j+1][i] += Param_local * (1.0/8.0) * (-1.0/fdz);	 
							zdPardu[k+1][j+1][i] += Param_local * (1.0/8.0) * ( 1.0/fdy);	 
							zdPardu[k+1][j  ][i] += Param_local * (1.0/8.0) * (-1.0/fdy);	
						}
						else if (!strcmp(IOparam->ObsName[ii],"E2nd"))
						{
							// not yet implemented!
						}
					
						// Store observation
						sty[ii] = Parameter*(1.0/scal->time_si);
						
					}




					// PetscPrintf(PETSC_COMM_SELF,"DEBUGDEBUGDEBUGDEBUG %.10f, %.10f, %.10f, %.10f; %.10f %d %d %d %d\n\n",svBulk->phi,IOparam->Ae[ii],dphidu_local,mfitCenter,xdphidu[k][j  ][i  ],i,j,k,ii);
				}
			}
		}
	}
	END_STD_LOOP

	ierr = VecRestoreArray(aop->sty,&sty); CHKERRQ(ierr);
	ierr = PetscBarrier((PetscObject)aop->sty); CHKERRQ(ierr);
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&mfitParam, &IOparam->mfitCenter, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		IOparam->mfitCenter = mfitParam;
	}
	if (IOparam->Gr == 0)
	{
		IOparam->mfitCenter /= 2.0;
	}

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx,    &dxx);    	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy,    &dyy);    	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz,    &dzz);    	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy,    &dxy);    	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz,    &dxz);    	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz,    &dyz);    	CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,     &vx);     	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,     &vy);     	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,     &vz);     	CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X, lxPar,    &xPar);    		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, lyPar,    &yPar);    		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, lzPar,    &zPar);    		CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X, lxdPardu, &xdPardu); 		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, lydPardu, &ydPardu); 		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, lzdPardu, &zdPardu); 		CHKERRQ(ierr);

	LOCAL_TO_GLOBAL(fs->DA_X, lxPar, gxPar);
	LOCAL_TO_GLOBAL(fs->DA_Y, lyPar, gyPar);
	LOCAL_TO_GLOBAL(fs->DA_Z, lzPar, gzPar);

	ierr = VecGetArray(gxPar, &dggxPar);      CHKERRQ(ierr);
	ierr = VecGetArray(gyPar, &dggyPar);      CHKERRQ(ierr);
	ierr = VecGetArray(gzPar, &dggzPar);      CHKERRQ(ierr);

	ierr = VecGetArray(jr->phi, &tempPar);        CHKERRQ(ierr);
	iter = tempPar;

	ierr  = PetscMemcpy(iter, dggxPar, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, dggyPar, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, dggzPar, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr = VecRestoreArray(jr->phi, &tempPar);         CHKERRQ(ierr);

	ierr = VecRestoreArray(gxPar, &dggxPar);      CHKERRQ(ierr);
	ierr = VecRestoreArray(gyPar, &dggyPar);      CHKERRQ(ierr);
	ierr = VecRestoreArray(gzPar, &dggzPar);      CHKERRQ(ierr);

	LOCAL_TO_GLOBAL(fs->DA_X, lxdPardu, gxdPardu);
	LOCAL_TO_GLOBAL(fs->DA_Y, lydPardu, gydPardu);
	LOCAL_TO_GLOBAL(fs->DA_Z, lzdPardu, gzdPardu);

	ierr = VecGetArray(gxdPardu, &dggxdPardu);      CHKERRQ(ierr);
	ierr = VecGetArray(gydPardu, &dggydPardu);      CHKERRQ(ierr);
	ierr = VecGetArray(gzdPardu, &dggzdPardu);      CHKERRQ(ierr);

	ierr = VecGetArray(aop->dPardu, &tempdPardu);        CHKERRQ(ierr);
	iter = tempdPardu;

	ierr  = PetscMemcpy(iter, dggxdPardu, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, dggydPardu, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, dggzdPardu, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr = VecRestoreArray(aop->dPardu, &tempdPardu);         CHKERRQ(ierr);

	ierr = VecRestoreArray(gxdPardu, &dggxdPardu);      CHKERRQ(ierr);
	ierr = VecRestoreArray(gydPardu, &dggydPardu);      CHKERRQ(ierr);
	ierr = VecRestoreArray(gzdPardu, &dggzdPardu);      CHKERRQ(ierr);

	// Destroy
	ierr = VecDestroy(&gxPar);   CHKERRQ(ierr);
	ierr = VecDestroy(&gyPar);   CHKERRQ(ierr);
	ierr = VecDestroy(&gzPar);   CHKERRQ(ierr);

	ierr = VecDestroy(&lxPar);   CHKERRQ(ierr);
	ierr = VecDestroy(&lyPar);   CHKERRQ(ierr);
	ierr = VecDestroy(&lzPar);   CHKERRQ(ierr);

	ierr = VecDestroy(&gxdPardu);   CHKERRQ(ierr);
	ierr = VecDestroy(&gydPardu);   CHKERRQ(ierr);
	ierr = VecDestroy(&gzdPardu);   CHKERRQ(ierr);

	ierr = VecDestroy(&lxdPardu);   CHKERRQ(ierr);
	ierr = VecDestroy(&lydPardu);   CHKERRQ(ierr);
	ierr = VecDestroy(&lzdPardu);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "cellConstEqFD"
PetscErrorCode cellConstEqFD(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarCell  *svCell, // solution variables
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz,    // ...
		PetscScalar &sxx,    // Cauchy stress components
		PetscScalar &syy,    // ...
		PetscScalar &szz,    // ...
		PetscScalar &gres,   // volumetric residual
		PetscScalar &rho,    // effective density
		AdjGrad *aop,
		ModParam *IOparam,
		PetscInt ii, 
		PetscInt jj, 
		PetscInt k, 
		PetscInt ik, 
		PetscInt jk, 
		PetscInt kk)
{
	// evaluate constitutive equations on the cell

	SolVarDev   *svDev;
	SolVarBulk  *svBulk;
	Controls    *ctrl;
	PetscScalar  eta_st, ptotal, txx, tyy, tzz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	svDev  = ctx->svDev;
	svBulk = ctx->svBulk;
	ctrl   = ctx->ctrl;

	// evaluate deviatoric constitutive equation
	ierr = devConstEqFD(ctx, aop, IOparam,  ii,  jj,  k,  ik,  jk,  kk); CHKERRQ(ierr);

	// evaluate volumetric constitutive equation
	ierr = volConstEq(ctx); CHKERRQ(ierr);

	// get stabilization viscosity
	if(ctrl->initGuess) eta_st = 0.0;
	else                eta_st = svDev->eta_st;

	// compute stabilization stresses
	sxx = 2.0*eta_st*svCell->dxx;
	syy = 2.0*eta_st*svCell->dyy;
	szz = 2.0*eta_st*svCell->dzz;

	// compute history shear stress
	svCell->sxx = 2.0*ctx->eta*dxx;
	svCell->syy = 2.0*ctx->eta*dyy;
	svCell->szz = 2.0*ctx->eta*dzz;

	// compute plastic strain-rate components
	txx = ctx->DIIpl*dxx;
	tyy = ctx->DIIpl*dyy;
	tzz = ctx->DIIpl*dzz;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = 0.5*(txx*txx + tyy*tyy + tzz*tzz);

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	txx = svCell->dxx - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tyy = svCell->dyy - svDev->I2Gdt*(svCell->syy - svCell->hyy);
	tzz = svCell->dzz - svDev->I2Gdt*(svCell->szz - svCell->hzz);

	// compute shear heating term contribution
	svDev->Hr =
		txx*svCell->sxx + tyy*svCell->syy + tzz*svCell->szz +
		sxx*svCell->dxx + syy*svCell->dyy + szz*svCell->dzz;

	// compute total viscosity
	svDev->eta = ctx->eta + eta_st;

	// get total pressure (effective pressure + pore pressure)
	ptotal = ctx->p + ctrl->biot*ctx->p_pore;

	// compute total Cauchy stresses
	sxx += svCell->sxx - ptotal;
	syy += svCell->syy - ptotal;
	szz += svCell->szz - ptotal;

	// save output variables
	svCell->eta_cr = ctx->eta_cr; // creep viscosity
	svCell->DIIdif = ctx->DIIdif; // relative diffusion creep strain rate
	svCell->DIIdis = ctx->DIIdis; // relative dislocation creep strain rate
	svCell->DIIprl = ctx->DIIprl; // relative Peierls creep strain rate
	svCell->yield  = ctx->yield;  // average yield stress in control volume

	// compute volumetric residual
	if(ctrl->actExp)
	{
		gres = -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta + svBulk->alpha*(ctx->T - svBulk->Tn)/ctx->dt;
	}
	else
	{
		gres = -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta;
	}

	// store effective density
	rho = svBulk->rho;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "setUpPhaseFD"
PetscErrorCode setUpPhaseFD(ConstEqCtx *ctx, PetscInt ID, AdjGrad *aop, ModParam *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk)
{
	// setup phase parameters for deviatoric constitutive equation
	// evaluate dependence on constant parameters (pressure, temperature)

	Material_t  *mat;
	Soft_t      *soft;
	Controls    *ctrl;
	PData       *Pd;
	PetscScalar  APS, Le, dt, p, p_lith, p_pore, T, mf, mfd, mfn;
	PetscScalar  Q, RT, ch, fr, p_visc, p_upper, p_lower, dP, p_total;
	PetscScalar  Inin,Inieta,ViscTemp;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	mat    = ctx->phases + ID;
	soft   = ctx->soft;
	ctrl   = ctx->ctrl;
	Pd     = ctx->Pd;
	APS    = ctx->svDev->APS;
	Le     = ctx->Le;
	dt     = ctx->dt;
	p      = ctx->p;
	p_lith = ctx->p_lith;
	p_pore = ctx->p_pore;
	T      = ctx->T;
	mf     = 0.0;

	if(mat->pdAct == 1)
	{
		// compute melt fraction from phase diagram
		ierr = setDataPhaseDiagram(Pd, p, T, mat->pdn); CHKERRQ(ierr);

		// store melt fraction
		mf = Pd->mf;
	}

	// set RT
	RT         =  ctrl->Rugc*T;
	if(!RT) RT = -1.0;

	// initialize phase parameters
	ctx->A_els = 0.0; // elasticity constant
	ctx->A_dif = 0.0; // diffusion constant
	ctx->A_max = 0.0; // upper bound constant
	ctx->A_dis = 0.0; // dislocation constant
	ctx->N_dis = 1.0; // dislocation exponent
	ctx->A_prl = 0.0; // Peierls constant
	ctx->N_prl = 1.0; // Peierls exponent
	ctx->taupl = 0.0; // plastic yield stress

	// MELT FRACTION
	mfd = 1.0;
	mfn = 1.0;

	if(!strcmp(IOparam->type_name[0],"eta"))
	{
		if ((ii)==ik && (jj)==jk && (k)==kk)
		{
			ViscTemp = (pow((mat->Bn * pow(2,mat->n) * pow(IOparam->DII_ref, mat->n-1)) , -1/mat->n));
			Inieta = ViscTemp;
			aop->Perturb = ViscTemp*aop->FD_epsilon;
			ViscTemp += aop->Perturb;
			mat->Bn = pow (2.0*ViscTemp, -mat->n) * pow(IOparam->DII_ref, 1.0 - mat->n);
		}
	}

	if(mf)
	{
		// limit melt fraction
		if(mf > ctrl->mfmax) mf = ctrl->mfmax;

		// compute corrections factors for diffusion & dislocation creep
		mfd = exp(mat->mfc*mf);
		mfn = exp(mat->mfc*mf*mat->n);
	}

	// PRESSURE

	// pore pressure
	if(ctrl->gwType == _GW_NONE_) p_pore = 0.0;

	// total pressure
	p_total = p + ctrl->biot*p_pore;

	// assign pressure for viscous laws
	if(ctrl->pLithoVisc)  p_visc = p_lith;
	else                  p_visc = p_total;

	// ELASTICITY
	if(mat->G)
	{
		// Elasticity correction can only DECREASE the viscosity.
		// eta/G << dt (viscous regime)  eta*(dt/(dt + eta/G)) -> eta
		// eta/G >> dt (elastic regime)  eta*(dt/(dt + eta/G)) -> G*dt < eta
		// Elasticity doesn't normally interact with the bottom viscosity limit,
		// instead it rather acts as a smooth limiter for maximum viscosity.

		ctx->A_els = 1.0/(mat->G*dt)/2.0;
	}

	// LINEAR DIFFUSION CREEP (NEWTONIAN)
	if(mat->Bd)
	{
		Q          = (mat->Ed + p_visc*mat->Vd)/RT;
		ctx->A_dif =  mat->Bd*exp(-Q)*mfd;
	}

	// PS-CREEP
	else if(mat->Bps && T)
	{
		Q          = mat->Eps/RT;
		ctx->A_dif = mat->Bps*exp(-Q)/T/pow(mat->d, 3.0);
	}

	// UPPER BOUND CREEP
	if(ctrl->eta_max)
	{
		ctx->A_max = 1.0/(ctrl->eta_max)/2.0;
	}

	// DISLOCATION CREEP (POWER LAW)
	if(mat->Bn)
	{
		Q          = (mat->En + p_visc*mat->Vn)/RT;
		if(!strcmp(IOparam->type_name[0],"n"))
		{
			if ((ii)==ik && (jj)==jk && (k)==kk)
			{
				ViscTemp = (pow((mat->Bn * pow(2,mat->n) * pow(IOparam->DII_ref, mat->n-1)) , -1/mat->n));
				Inin = mat->n;
				aop->Perturb = mat->n*aop->FD_epsilon;
				mat->n += aop->Perturb;
				mat->Bn = pow (2.0*ViscTemp, -mat->n) * pow(IOparam->DII_ref, 1.0 - mat->n);
			}
		}
		ctx->N_dis =  mat->n;
		ctx->A_dis =  mat->Bn*exp(-Q)*mfn;
	}

	// DC-CREEP
	else if(mat->Bdc && T)
	{
		Q          = mat->Edc/RT;
		ctx->N_dis = Q;
		ctx->A_dis = mat->Bdc*exp(-Q*log(mat->Rdc))*pow(mat->mu, -Q);
	}

	// PEIERLS CREEP (LOW TEMPERATURE RATE-DEPENDENT PLASTICITY, POWER-LAW APPROXIMATION)
	if(mat->Bp && T)
	{
		Q           = (mat->Ep + p_visc*mat->Vp)/RT;
		ctx->N_prl =  Q*pow(1.0-mat->gamma, mat->q-1.0)*mat->q*mat->gamma;
		ctx->A_prl =  mat->Bp/pow(mat->gamma*mat->taup, ctx->N_prl)*exp(-Q*pow(1.0-mat->gamma, mat->q));
	}

	if(!strcmp(IOparam->type_name[0],"n"))
	{
		if ((ii)==ik && (jj)==jk && (k)==kk)
		{
			mat->n = Inin;
			mat->Bn = pow (2.0*ViscTemp, -mat->n) * pow(IOparam->DII_ref, 1.0 - mat->n);
		}
	}
	if(!strcmp(IOparam->type_name[0],"eta"))
	{
		if ((ii)==ik && (jj)==jk && (k)==kk)
		{
			mat->Bn = pow (2.0*Inieta, -mat->n) * pow(IOparam->DII_ref, 1.0 - mat->n);
		}
	}

	if(PetscIsInfOrNanScalar(ctx->A_dif)) ctx->A_dif = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_dis)) ctx->A_dis = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_prl)) ctx->A_prl = 0.0;

	// PLASTICITY
	if(!mat->ch && !mat->fr)
	{
		PetscFunctionReturn(0); // return if no plasticity is set
	}

	// apply strain softening to friction and cohesion
	ch = applyStrainSoft(soft, mat->chSoftID, APS, Le, mat->ch);
	fr = applyStrainSoft(soft, mat->frSoftID, APS, Le, mat->fr);

	// fit to limits
	if(ch < ctrl->minCh) ch = ctrl->minCh;
	if(fr < ctrl->minFr) fr = ctrl->minFr;

	// override total pressure with lithostatic if requested
	if(ctrl->pLithoPlast)
	{
		// Use lithostatic, rather than dynamic pressure to evaluate yield stress
		// This converges better, but does not result in localization of deformation & shear banding,
		// so only apply it for large-scale simulations where plasticity does not matter much

		p_total = p_lith;
	}
	else if(ctrl->pLimPlast)
	{
		// apply pressure limits

		// yielding surface: (S1-S3)/2 = (S1+S3)/2*sin(phi) + C*cos(phi)
		// pressure can be written as: P = (S1+S2+S3)/3 and P~=S2,then P=(S1+S3)/2
		// so the yield surface can be rewritten as:
		// P-S3=P*sin(phi) + C*cos(phi)   --> compression
		// S1-P=P*sin(phi) + C*cos(phi)   --> extension
		// under pure shear compression, S3=P_Lithos and S1=P_Lithos when extension
		// P = -( S3+C*cos(phi))/(sin(phi)-1)  --> compression
		// P = -(-S1+C*cos(phi))/(sin(phi)+1)  --> extension

		p_upper = -( p_lith + ch * cos(fr))/(sin(fr) - 1.0); // compression
		p_lower = -(-p_lith + ch * cos(fr))/(sin(fr) + 1.0); // extension

		if(p_total > p_upper) p_total = p_upper;
		if(p_total < p_lower) p_total = p_lower;
	}

	// compute cohesion and friction coefficient
	ch = cos(fr)*ch;
	fr = sin(fr);

	// compute effective mean stress
	dP = (p_total - p_pore);

	// compute yield stress
	if(dP < 0.0) ctx->taupl =         ch; // Von-Mises model for extension
	else         ctx->taupl = dP*fr + ch; // Drucker-Prager model for compression

	// correct for ultimate yield stress (if defined)
	if(ctrl->tauUlt) { if(ctx->taupl > ctrl->tauUlt) ctx->taupl = ctrl->tauUlt; }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "devConstEqFD"
PetscErrorCode devConstEqFD(ConstEqCtx *ctx, AdjGrad *aop, ModParam *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk)
{
	// evaluate deviatoric constitutive equations in control volume

	Controls    *ctrl;
	PetscScalar *phRat;
	SolVarDev   *svDev;
	Material_t  *phases;
	PetscInt     i, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ctrl      = ctx->ctrl;
	numPhases = ctx->numPhases;
	phRat     = ctx->phRat;
	svDev     = ctx->svDev;
	phases    = ctx->phases;

	// zero out results
	ctx->eta    = 0.0; // effective viscosity
	ctx->eta_cr = 0.0; // creep viscosity
	ctx->DIIdif = 0.0; // diffusion creep strain rate
	ctx->DIIdis = 0.0; // dislocation creep strain rate
	ctx->DIIprl = 0.0; // Peierls creep strain rate
	ctx->DIIpl  = 0.0; // plastic strain rate
	ctx->yield  = 0.0; // yield stress

	// zero out stabilization viscosity
	svDev->eta_st = 0.0;

	// viscous initial guess
	if(ctrl->initGuess)
	{
		ctx->eta    = ctrl->eta_ref;
		ctx->eta_cr = ctrl->eta_ref;
		ctx->DIIdif = 1.0;

		PetscFunctionReturn(0);
	}

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// setup phase parameters
			ierr = setUpPhaseFD(ctx, i, aop, IOparam,  ii,  jj,  k,  ik,  jk,  kk); CHKERRQ(ierr);

			// compute phase viscosities and strain rate partitioning
			ierr = getPhaseVisc(ctx, i); CHKERRQ(ierr);

			// update stabilization viscosity
			svDev->eta_st += phRat[i]*phases->eta_st;
		}
	}

	// normalize strain rates
	if(ctx->DII)
	{
		ctx->DIIdif /= ctx->DII;
		ctx->DIIdis /= ctx->DII;
		ctx->DIIprl /= ctx->DII;
		ctx->DIIpl  /= ctx->DII;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "edgeConstEqFD"
PetscErrorCode edgeConstEqFD(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarEdge  *svEdge, // solution variables
		PetscScalar  d,      // effective shear strain rate component
		PetscScalar &s,
		AdjGrad *aop,
		ModParam *IOparam,
		PetscInt ii, 
		PetscInt jj, 
		PetscInt k, 
		PetscInt ik, 
		PetscInt jk, 
		PetscInt kk)      // Cauchy stress component
{
	// evaluate constitutive equations on the edge

	SolVarDev   *svDev;
	PetscScalar  t, eta_st;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	svDev = &svEdge->svDev;

	// evaluate deviatoric constitutive equation
	ierr = devConstEqFD(ctx, aop, IOparam,  ii,  jj,  k,  ik,  jk,  kk); CHKERRQ(ierr);

	// get stabilization viscosity
	if(ctx->ctrl->initGuess) eta_st = 0.0;
	else                     eta_st = svDev->eta_st;

	// compute stabilization stress
	s = 2.0*eta_st*svEdge->d;

	// compute history shear stress
	svEdge->s = 2.0*ctx->eta*d;

	// compute plastic strain-rate component
	t = ctx->DIIpl*d;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = t*t;

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	t = svEdge->d - svDev->I2Gdt*(svEdge->s - svEdge->h);

	// compute shear heating term contribution
	svDev->Hr = 2.0*t*svEdge->s + 2.0*svEdge->d*s;

	// compute total viscosity
	svDev->eta = ctx->eta + eta_st;

	// compute total stress
	s += svEdge->s;

	PetscFunctionReturn(0);
}
