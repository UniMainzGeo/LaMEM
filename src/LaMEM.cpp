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
 **    filename:   LaMEM.cpp
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
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

// FRAMEWORK CODE FOR LaMEM TO USE ADJOINT GRADIENT (INVERSION)
// *  developed by Georg Reuber (JGU Mainz)
// *  publication: Georg S. Reuber, Anton A. Popov, Boris J.P. Kaus, Deriving scaling laws in geodynamics using adjoint gradients, Tectonophysics, 2017
//---------------------------------------------------------------------------
// COMPUTATION OF ADJOINT INVERSION
//---------------------------------------------------------------------------
// RECIPE:
// Objective function    F(x,x(p)) = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)]      // p = parameter ; x = converged solution ; xini = comparison solution (same size as jr->gsol) ; P = Projection vector containing the proportions of solution influence
// Derivative I          dF/dx     = p*x-P*x_ini
// Adjoint operation     psi       = J^-T * dF/dx           // J = converged Jacobain matrix
// Derivative II         dr/dp     = [r(p+h) - r(p)]/h      // finite difference approximation of derivative of residual r vs parameter
// Gradients             dF/dp     = -psi^T * dr/dp
//
// ------------------------------------------------------------
// DETAILS:
// In your input file you need to define:
// IOparam.use	                   		= 1	                                // 0 = NO 1 = Free for other inversion methods 2 = Compute gradients 3 = full inversion 4 = save this forward simulation as comparison simulation
// IOparam.Tao                     		= 0                                 // 0 = Self written gradient descent using (BFGS) Hessian 1 = Use PETSCs TAO (LMVM/BLMVM)
// IOparam.mdN                     		= 4                                 // Number of parameters (same as lengths of phs & P )
// IOparam.mdI                     		= 1                                 // Number of indices (same as lengths of Ax & Ay & Az )
// IOparam.Ab                      		= 0                                 // 0 = No usage of bounds 1 = use bounds
// IOparam.Ap                      		= 1                                 // 1 = several indices ; 2 = the whole domain (will exclude mdI & Ax & Ay & Az) ; 3 = the point with maximum velocity
// IOparam.reg                     		= 1                                 // 1 = tikhonov regularization of the cost function (TN) 2 = total variation regularization (TV) (HANDLE WITH CARE)
// IOparam.OFdef                  		= 1                                 // Objective function defined by hand? - meaning you give the comparison value at specific point by hand
// IOparam.Ax					   		= {1.2}                             // Array (same length as Ay, Az) containing the x coordinates of the point where you want to compute the gradients
// IOparam.Ay					   		= {0.6}                             // Array (same length as Ax, Az) containing the y coordinates of the point where you want to compute the gradients
// IOparam.Az					   		= {0.4}                             // Array (same length as Ax, Ay) containing the z coordinates of the point where you want to compute the gradients
// IOparam.Ae					   		= {1}                               // Array (same length as Ax, Ay) containing the velocity value of the point where you want to compute the gradients (used then in cost function evaluation)
// IOparam.phs                     		= {1 2 1 2}                         // Array (same length as mdN) containing the phase of the parameter
// IOparam.typ                    		= {_RHO0_, _RHO0_, _ETA_,_ETA_}     // Array (same length as mdN) containing the parameter corresponding to the phase
// IOparam.Av                      		= {3}                               // Array (same length as Ax, Ay, Az, mdI) containing the related velocity direction in which to compute the gradient
// IOparam.Adv                     		= 1                                 // Should the point be advected?
// IOparam.tol 							= 1e-3;    							// tolerance for F/Fini after which code has converged
// (optional - tao) Lb                  = {1 0 0 0}                         // Array (same length as mdN) containing the lower bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
// (optional - tao) Ub                  = {5 8 8 8}                         // Array (same length as mdN) containing the upper bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
// (optional - GD)  IOparam.factor1     = 1e1 								// factor to multiply the gradients (should be set such that the highest gradient scales around 1/100 of its parameter)
// (optional - GD)  IOparam.factor2 	= 1.5  								// factor that increases the convergence velocity (this value is added to itself after every succesful gradient descent)
// (optional - GD)  IOparam.maxfactor2 	= 100 								// limit on the factor2
//
// ------------------------------------------------------------
// EXAMPLE IN INPUT FILE:
//  # General
//	Inv_use       = 2
//  Inv_Ap        = 1
//  Inv_OFdef     = 1
//  Inv_Tao       = 1
//  # Parameters
//  <InverseParStart>
//	   	Inv_ID  = 0
//		Inv_Typ = rho0
//		Inv_Par = 1
//	<InverseParEnd>
//  # Index
//	<InverseIndStart>
//		Inv_Ax = 4.95;
//		Inv_Ay = 0.05;
//		Inv_Az = 0.68;
//		Inv_Av = 3;
//		Inv_Ae = 1;
//	<InverseIndEnd>
//
// --> EXECUTE your input file (*.dat) (f.e. mpiexec -n 2 /home/user/lamem/bin/opt/LaMEM -ParamFile Input.dat)
//
// ------------------------------------------------------------
// Possible parameters (Inv_Typ):
//	_RHO0_, _RHON_, _RHOC_,                             // density
//	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
//	_ETA0_,	_E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
//	_BP_, _TAUP_, _GAMMA_, _Q_, _EP_, _VP_,             // Peierls creep
//	_SHEAR_, _BULK_, _KP_,                              // elasticity
//	_COHESION_, _FRICTION_, _CHSOFTID_, _FRSOFTID_,     // plasticity (Drucker-Prager)
//	_ALPHA_, _CP_, _K_, _A_                             // energy
//
// ------------------------------------------------------------
// Possible Velocities:  (Vx   Vy   Vz
//                       (1    2    3)
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
// 3) If not using tao, play around with the factors: factor1 should be set such that the smallest gradient multiplied by this factor is around 1/100 of its paramter.
//    Size of factor2 controls how fast the step size will increase (gives faster convergence but you might overstep couple of minimas).
// 
// 4) In case you want to use powerlaw viscosities make sure the reference strainrate is the same as the one that you use in the viscosity computation!!  
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
#include "scaling.h"
#include "objFunct.h"
#include "parsing.h"
#include "adjoint.h"
#include "phase.h"
//---------------------------------------------------------------------------
static char help[] = "Adjoint inversion computation .\n\n";
//--------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{

	PetscErrorCode 	ierr;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0, help); CHKERRQ(ierr);

	ModParam        IOparam;
	Scaling         scal;
	FB             *fb;
	PetscScalar    *gradar,*Ubar,*Lbar, *Par, *fcconvar, F, ts;
	Vec             val, Ub, Lb, grad, P;
	PetscInt        i, ti;
	char            str[_STR_LEN_];

	// set default to be a forward run
	IOparam.use        = 0;   		// 0 = forward run ; 1 = Neighbourhood algorithm (requires NAPlus) ; 2 = only compute adjoint gradients ; 3 = 'full' adjoint inversion with TAO ; 4 = assume this as a forward simulation and save the solution
	IOparam.mdI        = 0;
	IOparam.Ab         = 0;
	IOparam.Ap         = 2;
	IOparam.reg        = 0;
	IOparam.Adv        = 0;
	IOparam.OFdef      = 0;
	IOparam.Tao        = 1;
	IOparam.tol        = 0;
	IOparam.factor1    = 0;
	IOparam.factor1    = 0;
	IOparam.maxfactor2 = 0;
	F                  = 1e100;

	// load input file
	ierr = FBLoad(&fb); CHKERRQ(ierr);

	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_use"       , &IOparam.use,       1, 4        ); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_Ab"        , &IOparam.Ab,        1, 1        ); CHKERRQ(ierr);  // Apply bounds?
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_Ap"        , &IOparam.Ap,        1, 3        ); CHKERRQ(ierr);  // 1 = several indices ; 2 = the whole domain ; 3 = surface
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_reg"       , &IOparam.reg,       1, 2        ); CHKERRQ(ierr);  // 1 = tikhonov regularization of the cost function (TN) 2 = total variation regularization (TV)
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_Adv"       , &IOparam.Adv,       1, 1        ); CHKERRQ(ierr);  // 1 = advect the point
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_OFdef"     , &IOparam.OFdef,     1, 1        ); CHKERRQ(ierr);  // Objective function defined by hand?
	ierr = getIntParam   (fb, _OPTIONAL_, "Inv_Tao"       , &IOparam.Tao,       1, 1        ); CHKERRQ(ierr);  // Use TAO?
	ierr = getScalarParam(fb, _OPTIONAL_, "Inv_tol"       , &IOparam.tol,       1, 1        ); CHKERRQ(ierr);  // tolerance for F/Fini after which code has converged
	ierr = getScalarParam(fb, _OPTIONAL_, "Inv_factor1"   , &IOparam.factor1,   1, 1        ); CHKERRQ(ierr);  // factor to multiply the gradients (should be set such that the highest gradient scales around 1/100 of its parameter ; only used without tao)
	ierr = getScalarParam(fb, _OPTIONAL_, "Inv_factor2"   , &IOparam.factor1,   1, 1        ); CHKERRQ(ierr);  // factor that increases the convergence velocity (this value is added to itself after every succesful gradient descent ; only used without tao)
	ierr = getScalarParam(fb, _OPTIONAL_, "Inv_maxfactor2", &IOparam.maxfactor2,1, 1        ); CHKERRQ(ierr);  // limit on the factor (only used without tao)

	// Forward simulation
	if(IOparam.use == 0)
	{
		ierr = LaMEMLibMain(NULL); CHKERRQ(ierr);

		// cleanup PETSC
		ierr = PetscFinalize(); CHKERRQ(ierr);

		return 0;
	}

	IOparam.count = 1;  // iteration counter for the initial cost function

	// VECTORS
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Lb);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &Ub);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &val);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &P);    CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_  , PETSC_DETERMINE, &grad); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, 1501       , PETSC_DETERMINE, &IOparam.fcconv);   CHKERRQ(ierr);   // 1500 is the maximum inversion iterations that are accepted

	// TEMPORARY VARIABLES
	PetscInt		phsar[_MAX_PAR_];
	PetscInt      	typar[_MAX_PAR_];
	PetscScalar     W[_MAX_PAR_];
	PetscScalar     Ax[_MAX_IND_];
	PetscScalar     Ay[_MAX_IND_];
	PetscScalar     Az[_MAX_IND_];
	PetscScalar     Av[_MAX_IND_];
	PetscScalar     Ae[_MAX_IND_];
	
	ierr =  ScalingCreate(&scal, fb);



	// PARAMETERS
	// Get parameter / typ / etc.
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<InverseParStart>", "<InverseParEnd>"); CHKERRQ(ierr);

	// error checking
	if(fb->nblocks > _MAX_PAR_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many inverse parameters specified! Max allowed: %lld", (LLD)_MAX_PAR_);
	}
	if(!fb->nblocks)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define indices (mdi) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

	// read each individual parameter
	VecGetArray(P,&Par);
	VecGetArray(Ub,&Ubar);
	VecGetArray(Lb,&Lbar);
	VecGetArray(grad,&gradar);

	for(i = 0; i < fb->nblocks; i++)
	{
		ierr = getIntParam   (fb, _OPTIONAL_, "Inv_ID" , &ti, 1, max_num_phases); CHKERRQ(ierr);
		phsar[i]  = ti;                    // PHASE
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Par", &ts, 1, 1 ); CHKERRQ(ierr);
		Par[i]    = ts;                    // PARAMETER VALUES
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Uba", &ts, 1, 1 ); CHKERRQ(ierr);
		Ubar[i]   = ts;                    // UPPER BOUND
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Lba", &ts, 1, 1 ); CHKERRQ(ierr);
		Lbar[i]   = ts;                    // LOWER BOUND
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Wei", &ts, 1, 1 ); CHKERRQ(ierr);
		W[i]      = ts;                    // WEIGHTS
		ierr = getStringParam(fb, _OPTIONAL_, "Inv_Typ", str, NULL); CHKERRQ(ierr);
		if     (!strcmp(str, "rho0"))       { ti = _RHO0_; }
		else if(!strcmp(str, "rhon"))       { ti = _RHON_; }
		else if(!strcmp(str, "rhoc"))       { ti = _RHOC_; }
		else if(!strcmp(str, "eta"))        { ti = _ETA_;  }
		else if(!strcmp(str, "eta0"))       { ti = _ETA0_; }
		else if(!strcmp(str, "n"))          { ti = _N_;    }
		else if(!strcmp(str, "En"))         { ti = _EN_;   }
		else{ SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "WARNING: inversion parameter type is not yet implemented \n"); }
		typar[i]  = ti;
		gradar[i] = 0;                     // GRADIENTS

		fb->blockID++;
	}

	IOparam.phs = phsar;
	IOparam.typ = typar;
	IOparam.grd = gradar;
	IOparam.W = W;
	VecRestoreArray(P,&Par);
	VecRestoreArray(Ub,&Ubar);
	VecRestoreArray(Lb,&Lbar);
	VecRestoreArray(grad,&gradar);
	IOparam.mdN = i;

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);



	// LOCATIONS
	// Get location / value / etc.
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<InverseIndStart>", "<InverseIndEnd>"); CHKERRQ(ierr);

	// error checking
	if(fb->nblocks > _MAX_IND_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many inverse indices specified! Max allowed: %lld", (LLD)_MAX_IND_);
	}
	if(!fb->nblocks)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define parameters (mdN) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}

	// read each individual index
	for(i = 0; i < fb->nblocks; i++)
	{
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Ax", &ts, 1, 1 ); CHKERRQ(ierr);
		Ax[i] = ts /scal.length;        // X-COORDINATE
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Ay", &ts, 1, 1 ); CHKERRQ(ierr);
		Ay[i] = ts /scal.length;        // Y-COORDINATE
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Az", &ts, 1, 1 ); CHKERRQ(ierr);
		Az[i] = ts /scal.length;        // Z-COORDINATE
		ierr = getIntParam   (fb, _OPTIONAL_, "Inv_Av", &ti, 1, 3); CHKERRQ(ierr);
		Av[i] = ti;                     // VELOCITY COMPONENT
		ierr = getScalarParam(fb, _OPTIONAL_, "Inv_Ae", &ts, 1, 1 ); CHKERRQ(ierr);
		Ae[i] = ts /scal.velocity;;     // VELOCITY VALUE

		fb->blockID++;
	}

	IOparam.Ax = Ax;
	IOparam.Ay = Ay;
	IOparam.Az = Az;
	IOparam.Av = Av;
	IOparam.Ae = Ae;
	IOparam.mdI = i;

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);


	// Error checking
	if (IOparam.use != 0 && !IOparam.Ap)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have to define indices (Ap) for the inversion. Have a look into the comments in src/LaMEM.cpp");
	}
	if (IOparam.Tao == 0 && !IOparam.factor1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Without TAO you have to specify several additional hand tuned factors for the linesearch. Have a look into the comments in src/LaMEM.cpp");
	}

	//===============
	// SOLVE ADJOINT
	//===============
	// only compute the adjoint gradients or simply forward code
	if(IOparam.use == 2)
 	{
 		VecDuplicate(P,&IOparam.P);
 		VecCopy(P,IOparam.P);

 		// call LaMEM main library function
 		ierr = LaMEMLibMain(&IOparam); CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint inversion
 	else if(IOparam.use == 3)
 	{

 		VecDuplicate(P,&IOparam.P);
 		VecCopy(P,IOparam.P);

 		// if tao is used try the LMVM/BLMVM algorithms
 		if(IOparam.Tao == 1)
 		{
 	 		Tao tao;

 	 		ierr = TaoCreate(PETSC_COMM_WORLD,&tao); CHKERRQ(ierr);

 	 	 	// 1.Check if bounds are available and sets them
 	 		if (IOparam.Ab == 1)
 	 		{
 	 	 	 	ierr = TaoSetVariableBounds(tao,Lb,Ub);	 								CHKERRQ(ierr);
 	 	 	 	ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr);
 	 		}
 	 		else
 	 		{
 	 			ierr = TaoSetType(tao,TAOLMVM);CHKERRQ(ierr);                    // TAOLMVM, TAOBLMVM, TAOBMRM (bad), TAOCG (all 4), TAOTRON (might crash but is fast)
 	 		}

 	 		// 2. Set up Tao
 	 	 	ierr = TaoSetObjectiveAndGradientRoutine(tao, AdjointOptimisationTAO, &IOparam);	 	CHKERRQ(ierr);
 	 	 	ierr = TaoSetInitialVector(tao,P);	 									CHKERRQ(ierr);
 	 	 	ierr = TaoSetTolerances(tao,1e-30,1e-30,1e-30);	CHKERRQ(ierr);
 	 	 	ierr = TaoSetFunctionLowerBound(tao,1e-5);CHKERRQ(ierr);
 	 	 	ierr = TaoSetFromOptions(tao);	 										CHKERRQ(ierr);

 	 	 	// 3. Solve Tao & view result
 	 	 	ierr = TaoSolve(tao);	 												CHKERRQ(ierr);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 	 	 	TaoView(tao,PETSC_VIEWER_STDOUT_WORLD);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");

 	 	 	// 4. Clean
 	 	 	ierr = TaoDestroy(&tao);
 		}
 		else
 		// Without TAO try line search tuned gradient descent
 		{
 			ierr = AdjointOptimisation(P, F, grad, &IOparam);
 		}

 		PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------\n");
 		PetscPrintf(PETSC_COMM_WORLD,"*         INVERSION RESULT SUMMARY       *\n");
 		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 		PetscPrintf(PETSC_COMM_WORLD,"Number of inversion iterations: %d\n",IOparam.count);
 		PetscPrintf(PETSC_COMM_WORLD,"F/Fini:\n");
 		VecGetArray(IOparam.fcconv,&fcconvar);
 		for(i=0;i<IOparam.count;i++)
 		{
 			PetscPrintf(PETSC_COMM_WORLD,"%8.15f\n",fcconvar[i]);
 		}
 		VecRestoreArray(IOparam.fcconv,&fcconvar);
 		PetscPrintf(PETSC_COMM_WORLD,"\nFinal cost function:\n");
 		PetscPrintf(PETSC_COMM_WORLD,"%8.15f\n",IOparam.mfit);
 		PetscPrintf(PETSC_COMM_WORLD,"\nFinal Parameters:\n");
		VecGetArray(IOparam.P,&Par);
		for(i=0;i<IOparam.mdN;i++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%.12f\n",Par[i]);
		}
		VecRestoreArray(IOparam.P,&Par);
 		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n\n");

 	}
 	// this is a forward simulation that we want to save as comparison solution
 	else if(IOparam.use == 4)
 	{
 		// call LaMEM main library function
 		ierr = LaMEMLibMain(&IOparam); CHKERRQ(ierr);

 		// Save output
 		PetscViewer     viewer;
 		PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
 		PetscViewerSetType(viewer,PETSCVIEWERBINARY);
 		PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
 		PetscViewerFileSetName(viewer,"Forward_Solution.bin");
 	 	VecView(IOparam.xini,viewer);
 	 	PetscViewerDestroy(&viewer);
 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n        Forward Solution succesfully saved\n------------------------------------------\n");
 	}

	ierr = PetscMemzero(&IOparam, sizeof(ModParam)); CHKERRQ(ierr);

	// cleanup PETSC
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointOptimisation"
PetscErrorCode AdjointOptimisation(Vec P, PetscScalar F, Vec grad, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	PetscInt 		i, j, a, LScount;
	PetscScalar 	*Par, *Paroldar, *gradar, *dPtemp, *fcconvar, b;
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

	// Initialize parameters
	a = -1.0;
	b = IOparam->factor2;   // initial value for factor2

	// Initialize cost functions
	F = 1e100;
	Fold = 1e100;

	while(F>IOparam->tol)
	{
		// Give the updated values to the code
		VecCopy(P,IOparam->P);

		// Reset line search counter
		LScount = 1;

		// call LaMEM main library function
		ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);

		// Save intial cost function & create initial Hessian
		if(IOparam->count==1)
		{
			IOparam->mfitini = IOparam->mfit;
		}

		// Save cost function
		F = IOparam->mfit;

		// If cost function in this timestep is larger then before perform bisection line search
		while(F>Fold)
		{
			b = 1;
			PetscPrintf(PETSC_COMM_WORLD,"\n- - - - - - - - - - - - - - - - - - - - - - - - - - - \n");
			PetscPrintf(PETSC_COMM_WORLD,"              LINE SEARCH IT %d                       \n",LScount);

			VecGetArray(P,&Par);
			VecGetArray(Pold,&Paroldar);
			VecGetArray(dP,&dPtemp);
			for(i=0;i<IOparam->mdN;i++)
			{
				for(j=0;j<IOparam->mdN;j++)
				{
					dPtemp[i] = dPtemp[i]/2;
				}
			}

			// Update parameter
			for(i=0;i<IOparam->mdN;i++)
			{
				Par[i] = Paroldar[i] + dPtemp[i];
			}
			VecRestoreArray(P,&Par);
			VecRestoreArray(Pold,&Paroldar);
			VecRestoreArray(dP,&dPtemp);

			// Give the updated values to the code
			VecCopy(P,IOparam->P);

			// call LaMEM main library function
			ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);

			F = IOparam->mfit;

			LScount+=1;
			if(LScount>10)
			{
				PetscPrintf(PETSC_COMM_WORLD,"******************************************************\n");
				PetscPrintf(PETSC_COMM_WORLD,"*             YOUR SOLUTION DIVERGED                 *\n");
				PetscPrintf(PETSC_COMM_WORLD,"******************************************************\n\n");

				// Return parameters for final output
				VecCopy(P,IOparam->P);

				PetscFunctionReturn(0);
			}
		}

		// Zero out perturbation
		VecDuplicate(IOparam->P,&dP);

		// restore parameter values
		VecDuplicate(IOparam->P,&P);
		VecCopy(IOparam->P,P);

		VecGetArray(grad,&gradar);
		for(j = 0; j < IOparam->mdN; j++)
		{
			gradar[j] = IOparam->grd[j];
		}
		VecRestoreArray(grad,&gradar);

		PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"%d. IT INVERSION RESULT: line search its = %d ; F / FINI = %.12f\n\n",IOparam->count,LScount-1,IOparam->mfit/IOparam->mfitini);
		PetscPrintf(PETSC_COMM_WORLD,"FOLD = %.20f \n   F = %.20f\n\n",Fold,F);

		// BEFORE UPDATING the par vector store the old gradient & Parameter vector
		VecCopy(P,Pold);
		VecCopy(grad,gradold);
		Fold = F;

		VecGetArray(grad,&gradar);
		VecGetArray(P,&Par);
		VecGetArray(dP,&dPtemp);
		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				dPtemp[i] = dPtemp[i] + (gradar[i]);
			}
		}

		VecRestoreArray(dP,&dPtemp);
		PetscReal min;
		ierr =  VecMin(dP,NULL,&min);
		if(min<0)
		{
			min = -min;
		}
		VecGetArray(dP,&dPtemp);

		PetscPrintf(PETSC_COMM_WORLD,"Factor2 = %2.1f\n\n",b);
		for(i=0;i<IOparam->mdN;i++)
		{
			dPtemp[i] = b*dPtemp[i]*a;
		}

		b = b+IOparam->factor2;
		if(b>IOparam->maxfactor2)
		{
			b = IOparam->maxfactor2;
		}

		// Display the current state of the parameters
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%D. Diff parameter value = %.12f\n",j+1,dPtemp[j]);
		}

		PetscPrintf(PETSC_COMM_WORLD,"\n");

		// Update parameter
		for(i=0;i<IOparam->mdN;i++)
		{
			Par[i] = Par[i] + dPtemp[i];
		}
		VecRestoreArray(grad,&gradar);
		VecRestoreArray(P,&Par);
		VecRestoreArray(dP,&dPtemp);

		// Display the current state of the parameters
		VecGetArray(P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%D. Parameter value = %.12f\n",j+1,Par[j]);
		}
		VecRestoreArray(P,&Par);

		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n\n");

		// Give the updated values to the code  (actually unfortunately necessary here and at the top of this function - need to rearrange that)
		VecCopy(P,IOparam->P);

		VecGetArray(IOparam->fcconv,&fcconvar);
		fcconvar[IOparam->count] = IOparam->mfit/IOparam->mfitini;
		VecRestoreArray(IOparam->fcconv,&fcconvar);

		// count
		IOparam->count += 1;
		if(IOparam->count>1500)
		{
			PetscPrintf(PETSC_COMM_WORLD,"\n\n\nEXCEEDED 1500 FUNCTION EVALUATIONS (consider changing inversion options)\n\n\n");
			PetscFunctionReturn(0);
		}

	}

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
	PetscScalar *Par, *gradar, *fcconvar;
	ModParam    *IOparam;
	IOparam     = (ModParam*)ctx;

	// get parameter values
	VecDuplicate(P,&IOparam->P);
	VecCopy(P,IOparam->P);

	// call LaMEM main library function
	ierr = LaMEMLibMain(IOparam); CHKERRQ(ierr);

	// restore parameter values
	VecDuplicate(IOparam->P,&P);
	VecCopy(IOparam->P,P);

	// Store the gradient & misfit
	VecGetArray(grad,&gradar);
	for(j = 0; j < IOparam->mdN; j++)
	{
		gradar[j] = IOparam->grd[j];
	}
	VecRestoreArray(grad,&gradar);

	*F = IOparam->mfit;

	// Save intial cost function
	if(IOparam->count==1){IOparam->mfitini = IOparam->mfit;}

	// Display the current state of the parameters
	VecGetArray(P,&Par);
	for(j = 0; j < IOparam->mdN; j++)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%D. Parameter value = %.12f\n",j+1,Par[j]);
	}
	VecRestoreArray(P,&Par);

	// Relative cost function
	PetscPrintf(PETSC_COMM_WORLD,"mfit / mfit0 = %.12f\n------------------------------------------\n\n",IOparam->mfit/IOparam->mfitini);

	VecGetArray(IOparam->fcconv,&fcconvar);
	fcconvar[IOparam->count] = IOparam->mfit/IOparam->mfitini;
	VecRestoreArray(IOparam->fcconv,&fcconvar);

	// count
	IOparam->count += 1;
	if(IOparam->count>1500)
	{
		PetscPrintf(PETSC_COMM_WORLD,"\n\n\nEXCEEDED 1500 FUNCTION EVALUATIONS (consider changing inversion options)\n\n\n");
		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}

















