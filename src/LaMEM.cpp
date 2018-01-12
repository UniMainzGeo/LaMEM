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

















