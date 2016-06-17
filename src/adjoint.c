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
 **    filename:   input.c
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
// USAGE:
// In your .dat file you need to define:
// ComputeAdjointGradients	       = 1	                    // 0 = NO 1 = Compute gradients 2 = full inversion 3 = save this forward simulation as comparison simulation
// Adjoint_x					   = 1.2                    // Array (same length as Adjoint_y, Adjoint_z) containing the x coordinates of the point where you want to compute the gradients
// Adjoint_y					   = 0.6                    // Array (same length as Adjoint_x, Adjoint_z) containing the y coordinates of the point where you want to compute the gradients
// Adjoint_z					   = 0.4                    // Array (same length as Adjoint_x, Adjoint_y) containing the z coordinates of the point where you want to compute the gradients
// AdjointPhases                   = 1 2 1 2                // Array (same length as AdjointParameters) containing the phase of the parameter
// AdjointParameters               = 2 1 1 1                // Array (same length as AdjointPhases) containing the parameter corresponding to the phase
// AdjointVel                      = 3                      // Array (same length as Adjoint_x, Adjoint_y, Adjoint_z) containing the related velocity direction in which to compute the gradient
// (optional) AdjointLowerBound               = 1 0 0 0		// Array (same length as AdjointParameters) containing the lower bounds for the parameters (only taken into account if '-tao_type blmvm')
// (optional) AdjointUpperBound               = 5 8 8 8		// Array (same length as AdjointParameters) containing the upper bounds for the parameters (only taken into account if '-tao_type blmvm')
//
//                              Density         Elasticity  Diff creep    Dis creep    Peierl creep
// Possible parameters: (rho rho_n rho_c beta)   (K Kp G)   (Bd Ed Vd)  (Bn n  En Vn) (taup gamma q )
//                      ( 1     2     3    4 )   (5  6 7)   ( 8  9 10)  (11 12 13 14) ( 15    16  17)
//
// Remark: If you use parameter 8 or 11 you will get to see the gradient or parameter value for eta (NOT for Bd or Bn itself) which is 1/Bd or 1/Bn.
//
// Possible Velocities:  (Vx   Vy   Vz
//                       (1    2    3)
//
// LINEAR SOLVER:
// You can control the behaviour of the KSP object for the adjoint with the prefix "as_" (probably the same as "js_")
//
// FULL INVERSION REMARKS:
// 1) In case you want to perform the full adjoint inversion (ComputeAdjointGradients = 2) make sure that you have a comparison file with a petsc vector the same size as jr->gsol
//    and called Forward_Solution.bin. Most likely you want to run a forward simulation and then change the parameters to do so just run your forward model with ComputeAdjointGradients = 3
//    which automatically saves this file. Perturb the values within the input script and solve again with ComputeAdjointGradients = 2.
//
// 2) You should make sure that you use '-tao_type lmvm', since this code is written for this method!
//    You can control the behaviour of the TAO object for the adjoint with the prefix "tao_" (example: '-tao_type lmvm' ; '-tao_fatol 1e-15' ; '-tao_converged_reason')
//
// 3) You can also perform the bounded inversion by setting '-tao_type blmvm' and defining upper (AdjointUpperBound) and lower bound (AdjointLowerBound).
//
// IMPORTANT REMARKS:
// 1) Since the Adjoint needs the Jacobian matrix for computing the gradients it's crucial to make sure that
//    you compute the Jacobian matrix in the timesteps where you want to compute the gradients
//    (f.e. a linear problem would need a low value for -snes_atol [1e-20] and a low max iteration count -snes_max_it [2] to
//    guarantee the computation of the Jacobian + the option '-snes_type ksponly')
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "paraViewOutBin.h"
#include "paraViewOutSurf.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "multigrid.h"
#include "advect.h"
#include "marker.h"
#include "paraViewOutMark.h"
#include "input.h"
#include "matProps.h"
#include "objFunct.h"
#include "AVDView.h"
#include "break.h"
#include "parsing.h"
#include "adjoint.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointDestroy"
PetscErrorCode AdjointDestroy(AdjGrad *aop)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Destroy the Adjoint gradients structures
	ierr = PetscMemzero(aop, sizeof(AdjGrad)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
 //---------------------------------------------------------------------------
 #undef __FUNCT__
 #define __FUNCT__ "AdjointOptimization"
 PetscErrorCode AdjointOptimization(JacRes *jr, AdjGrad *aop, UserCtx *user, NLSol *nl,SNES snes)
 {
	 /* TODO:
		- It's only valid vor visco-elasicity since the Jacobian is assumed to be symmetric (future: Jacobian needs to be transposed)
		- boundary of 1e-10 if Perturb parameter becomes too small in the finite differences
		- the only working comparison variable is the velocity
		- still not sure if it should be 1/(2*Bn) or 1/Bn
	  */

 	PetscErrorCode ierr;
 	PetscFunctionBegin;

 	Tao 			tao;
 	PetscReal 		*F;
	Vec         	x,P,grad;

	// Give all needed paramters to the aop handle (not efficient but Tao needs it like that ...)
 	aop->jr = jr;
 	aop->nl = nl;
 	aop->user = user;
 	aop->snes = snes;

 	// Initialize adjoint iteration counter
 	aop->count = 1;

	// Create everything
 	ierr = TaoCreate(PETSC_COMM_WORLD,&tao);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumInd, PETSC_DETERMINE, &aop->vx); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumInd, PETSC_DETERMINE, &aop->vy); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumInd, PETSC_DETERMINE, &aop->vz); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumPar, PETSC_DETERMINE, &grad); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumPar, PETSC_DETERMINE, &P); CHKERRQ(ierr);

 	// Check if bounds are available and sets them
	if (user->AdjointUseBounds == 1)
	{
		Vec 			lb,ub;
		PetscScalar 	*lbb,*ubb;
		PetscInt 		j;

		// Create the bound vectors
		ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumPar, PETSC_DETERMINE, &lb); CHKERRQ(ierr);
		ierr = VecCreateMPI(PETSC_COMM_WORLD, user->AdjointNumPar, PETSC_DETERMINE, &ub); CHKERRQ(ierr);

		VecGetArray(lb,&lbb);
		VecGetArray(ub,&ubb);

		// Set the boundaries from the user input
		for(j = 0; j < user->AdjointNumPar; j++)
		{
			lbb[j] = user->AdjointLowerBound[j];
			ubb[j] = user->AdjointUpperBound[j];
		}
		VecRestoreArray(lb,&lbb);
		VecRestoreArray(ub,&ubb);

		// Apply the bounds to the TAO object
 	 	ierr = TaoSetVariableBounds(tao,lb,ub);	 								CHKERRQ(ierr);

 	 	ierr = VecDestroy(&lb);
 	 	ierr = VecDestroy(&ub);
	 }

	//===============
	// SOLVE ADJOINT
	//===============
	// only compute the adjoint gradients
 	if(user->ComputeAdjointGradients == 1)
 	{
 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);	 	CHKERRQ(ierr);

 	 	// Get the solution
 	 	ierr = VecDuplicate(jr->gsol, &x);	 	CHKERRQ(ierr);
 	 	ierr = VecCopy(jr->gsol,x); 			CHKERRQ(ierr);

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, user, aop); 		CHKERRQ(ierr);

 		// Compute the gradients
 		ierr = AdjointObjectiveFunction(tao, P, F, aop);
 		ierr = AdjointComputeGradients(tao, P, grad, aop);
 	}
 	// compute 'full' adjoint inversion
 	else if(user->ComputeAdjointGradients == 2)
 	{
 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);	 	CHKERRQ(ierr);

 	 	// Get the solution
 	 	ierr = VecDuplicate(jr->gsol, &x);	 	CHKERRQ(ierr);
 	 	ierr = VecCopy(jr->gsol,x); 			CHKERRQ(ierr);

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, user, aop); 		CHKERRQ(ierr);

 		// 1. We need to read the comparison solution
 		if (aop->count == 1)
 		{
 		 	// Get the initial parameters & the comparison solution
 		 	ierr = AdjointExtractParameterVec(nl,P,user);
 	 		ierr = AdjointLoadCompareData(aop);
 		}

 		// 2. Set up Tao
 	 	ierr = TaoSetObjectiveRoutine(tao, AdjointObjectiveFunction, aop);	 	CHKERRQ(ierr);
 	 	ierr = TaoSetGradientRoutine(tao, AdjointComputeGradients, aop);	 	CHKERRQ(ierr);
 	 	ierr = TaoSetInitialVector(tao,P);	 									CHKERRQ(ierr);
 	 	ierr = TaoSetFromOptions(tao);	 										CHKERRQ(ierr);

 	 	// 3. Solve Tao & view result
 	 	ierr = TaoSolve(tao);	 												CHKERRQ(ierr);
 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 	 	TaoView(tao,PETSC_VIEWER_STDOUT_WORLD);
 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 	 	ierr = AdjointPerturbParameterVec(nl,P, user);
 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 	}
 	// this is a forward simulation that we want to save as comparison solution
 	else if(user->ComputeAdjointGradients == 3)
 	{
 		PetscViewer     viewer;
 		PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
 		PetscViewerSetType(viewer,PETSCVIEWERBINARY);
 		PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
 		PetscViewerFileSetName(viewer,"Forward_Solution.bin");
 	 	VecView(jr->gsol,viewer);
 	 	PetscViewerDestroy(&viewer);
 	 	PetscPrintf(PETSC_COMM_WORLD,"Forward Solution succesfully saved\n------------------------------------------\n");
 	}

 	// Clean
 	ierr = TaoDestroy(&tao);
 	ierr = VecDestroy(&x);
 	ierr = VecDestroy(&P);
 	ierr = VecDestroy(&grad);

 	PetscFunctionReturn(0);
 }
 //---------------------------------------------------------------------------
 #undef __FUNCT__
 #define __FUNCT__ "AdjointObjectiveFunction"
 PetscErrorCode AdjointObjectiveFunction(Tao tao, Vec P, PetscReal *F, void *ctx)
 {
 	PetscErrorCode ierr;
 	PetscFunctionBegin;

	AdjGrad 		*aop;
	JacRes 			*jr;
	UserCtx			*user;
	NLSol           *nl;
	PetscReal       *Fini;

 	aop = (AdjGrad*)ctx;

 	jr 		= aop->jr;
 	user 	= aop->user;
 	nl 		= aop->nl;

 	ierr = VecDuplicate(jr->gsol, &aop->dF);	 	CHKERRQ(ierr);

	//=============================
	// COMPUTE OBJECTIVE FUNCTION
	//=============================
 	// only compute the gradients (F = P*x (unused) // dF/dx = P)
 	if (user->ComputeAdjointGradients == 1)
 	{
 		ierr = VecCopy(aop->pro,aop->dF); 		CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint inversion
 	else if(user->ComputeAdjointGradients == 2)
	{
 		Vec xini;
 		PetscScalar Ad;

 		PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE COST FUNCTION\n******************************************\n");

	 	// Copy temporary comparison solution
	 	ierr = VecDuplicate(jr->gsol, &xini);	 	CHKERRQ(ierr);
	 	ierr = VecCopy(aop->xini,xini); 			CHKERRQ(ierr);

 		if (aop->count == 1){}
 		else
 		{
 			// Solve with the new parameters
 			ierr = AdjointPerturbParameterVec(nl,P,user);

 			PetscLogDouble     	cputime_start_nonlinear, cputime_end_nonlinear;
 			SNES 				snes;

 			snes 	= aop->snes;

			PetscTime(&cputime_start_nonlinear);

			PetscBool 			snes_convergence;
			snes_convergence 	=	PETSC_FALSE;

			// compute inverse elastic viscosities (dependent on dt)
			ierr = JacResGetI2Gdt(jr); CHKERRQ(ierr);

			// solve nonlinear system with SNES
			ierr = SNESSolve(snes, NULL, jr->gsol); CHKERRQ(ierr);

			// print analyze convergence/divergence reason & iteration count
			ierr = SNESPrintConvergedReason(snes, &snes_convergence); CHKERRQ(ierr);

			if (!snes_convergence)
			{
				PetscPrintf(PETSC_COMM_WORLD, " **** Nonlinear solver failed to converge *** \n");
			}

			PetscTime(&cputime_end_nonlinear);

			PetscPrintf(PETSC_COMM_WORLD, " Nonlinear solve took %g (sec)\n", cputime_end_nonlinear - cputime_start_nonlinear);

			// view nonlinear residual
			ierr = JacResViewRes(jr); CHKERRQ(ierr);
 		}

 		// Incorporate projection vector (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr = VecAYPX(xini,-1,jr->gsol);        		CHKERRQ(ierr);
 		ierr = VecPointwiseMult(xini, xini,aop->pro);	CHKERRQ(ierr);

 		// Compute objective function value (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr 	= VecDot(xini,xini,&Ad);
 		Ad 		/= 2;
 		*F 		= Ad;

 		// Save the initial cost function
 		if (aop->count == 1){aop->Fini = Ad;}

 	 	PetscPrintf(PETSC_COMM_WORLD,"Current Cost function = %.12f ; F/F0 = %.12f\n",(double)(*F), Ad/aop->Fini);

 		// Compute it's derivative (dF/dx = P*x-P*x_ini)
 		ierr = VecCopy(xini,aop->dF); 		CHKERRQ(ierr);

 		// Destroy & increase count
 		ierr = VecDestroy(&xini);
 		aop->count += 1;
	}
 	else
 	{
 	 	PetscPrintf(PETSC_COMM_WORLD,"ERROR: ComputeAdjointGradient value is not defined ; Choose between [1-3]\n");
 	 	PetscFunctionReturn(1);
 	}

 	PetscFunctionReturn(0);
 }
 /*//---------------------------------------------------------------------------
 #undef __FUNCT__
 #define __FUNCT__ "AdjointFormHessian"
 PetscErrorCode AdjointFormHessian(Tao tao, Mat H, Mat Hre, void *ctx)
 {
 	PetscErrorCode ierr;
 	PetscFunctionBegin;

	// You should really use lmvm (or the bounded version blmvm) tao_solve because this is using the bfgs approximated Hessian which is reasonable.
	// Anyway deriving the Hessian for our problem might be a mathematical overkill

 	PetscFunctionReturn(0);
 }*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointComputeGradients"
PetscErrorCode AdjointComputeGradients(Tao tao, Vec P, Vec grad, void *ctx)
{
	PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE GRADIENTS\n******************************************\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	AdjGrad 			*aop;
	JacRes 				*jr;
	UserCtx				*user;
	NLSol 				*nl;
	SNES 				snes;
	FDSTAG              *fs;
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, j, CurPhase, CurPar, lrank, grank;
	PetscScalar         grd, Perturb, coord_local[3], *gradvec, *vx, *vy, *vz;
	Vec 				rpl, sol, psi, drdp, res, Perturb_vec;
	PC                  ipc;
	Scaling             *scal;

 	aop     = (AdjGrad*)ctx;
 	jr 		= aop->jr;
 	user 	= aop->user;
 	nl 		= aop->nl;
	snes 	= aop->snes;

	fs = jr->fs;

	// get the current velocities at comparison point
	ierr = AdjointPointInPro(jr, user, aop); 		CHKERRQ(ierr);

	// Set perturbation paramter for the finite differences back to 1e-6
	aop->Perturb = 1e-6;

	// Profile time
	PetscLogDouble     cputime_start, cputime_end;
	PetscTime(&cputime_start);

	scal = &jr->scal;

	// Create all needed vectors in the same size as the solution vector
	ierr = VecDuplicate(jr->gsol, &psi);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &rpl);		 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &res);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sol); 		 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &drdp);	 	 CHKERRQ(ierr);
	ierr = VecCopy(jr->gsol,sol); 				 CHKERRQ(ierr);
	ierr = VecCopy(jr->gres,res); 				 CHKERRQ(ierr);

	//========
	// SOLVE
	//========
	// Solve the adjoint equation (psi = J^-T * dF/dx)
	ierr = SNESGetKSP(snes, &ksp);         		CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"as_"); 		CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         		CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            		CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          		CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,nl->J,nl->P);	CHKERRQ(ierr);
	ierr = KSPSolve(ksp,aop->dF,psi);			CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);	CHKERRQ(ierr);

	VecGetArray(grad,&gradvec);

	//=================
	// PARAMETER LOOP
	//=================
	for(j = 0; j < user->AdjointNumPar; j++)
	{
		// Get residual since it is overwritten in VecAYPX
		ierr = VecDuplicate(jr->gres, &res);  	CHKERRQ(ierr);
		ierr = VecCopy(jr->gres,res); 			CHKERRQ(ierr);

		// Get current phase and parameter
		CurPhase = user->AdjointPhases[j];
		CurPar   = user->AdjointParameters[j];

		// Perturb the current parameter in the current phase
		ierr = AdjointGradientPerturbParameter(nl, CurPar, CurPhase, aop);      CHKERRQ(ierr);

		// get the actual used perturbation parameter which is 1e-6*parameter
		Perturb = aop->Perturb;
		ierr = VecDuplicate(jr->gsol, &Perturb_vec); CHKERRQ(ierr);
		ierr = VecSet(Perturb_vec,Perturb);			 CHKERRQ(ierr);

		// Compute residual with the converged Jacobian (dr/dp = [r(p+h) - r(p)]/h)
		ierr = FormResidual(snes, sol, rpl, nl); 	       CHKERRQ(ierr);
		ierr = VecAYPX(res,-1,rpl);        CHKERRQ(ierr);
		ierr = VecPointwiseDivide(drdp,res,Perturb_vec);   CHKERRQ(ierr);

		// Compute the gradient (dF/dp = -psi^T * dr/dp) & Save gradient
		ierr = VecDot(drdp,psi,&grd);     CHKERRQ(ierr);
		gradvec[j] 	= -1 * grd;

		// Reset perturbed parameter
		ierr = AdjointGradientResetParameter(nl, CurPar, CurPhase, aop);       CHKERRQ(ierr);

		// Print result
		PetscPrintf(PETSC_COMM_WORLD,"%D.Gradient = %.12f ; CurPar = %d ; CurPhase = %d\n",j+1,gradvec[j]*scal->velocity,CurPar,CurPhase);

		// Destroy overwritten residual vector
		ierr = VecDestroy(&res);
	}

	VecRestoreArray(grad,&gradvec);
	VecGetArray(aop->vx,&vx);
	VecGetArray(aop->vy,&vy);
	VecGetArray(aop->vz,&vz);

	// Print the solution variable at the user defined index
	for (i=0; i<user->AdjointNumInd; i++)
	{
		coord_local[0] = user->Adjoint_x[i];
		coord_local[1] = user->Adjoint_y[i];
		coord_local[2] = user->Adjoint_z[i];

		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

		// If lrank is not 13 the point is not on this processor
		if(lrank == 13)
		{

			if (user->AdjointVel[i] == 1)
			{
				PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vx] (dimensional) = %.12f\n",vx[i]*scal->velocity);
			}
			else if (user->AdjointVel[i] == 2)
			{
				PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vy] (dimensional) = %.12f\n",vy[i]*scal->velocity);
			}
			else if (user->AdjointVel[i] == 3)
			{
				PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vz] (dimensional) = %.12f\n",vz[i]*scal->velocity);
			}
		}
	}

	VecRestoreArray(aop->vx,&vx);
	VecRestoreArray(aop->vy,&vy);
	VecRestoreArray(aop->vz,&vz);

	// Clean
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&drdp);
	ierr = VecDestroy(&rpl);
	ierr = VecDestroy(&Perturb_vec);

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful & took %g s\n",cputime_end - cputime_start);

	if (user->ComputeAdjointGradients == 1)
	{
		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPointInPro"
PetscErrorCode AdjointPointInPro(JacRes *jr, UserCtx *user, AdjGrad *aop)
{
	PetscErrorCode      ierr;
	FDSTAG              *fs;
	Vec                 lproX, lproY, lproZ, gproX, gproY, gproZ, pro;
	PetscScalar         coord_local[3], *temppro, ***llproX, ***llproY, ***llproZ, *dggproX, *dggproY, *dggproZ;
	PetscScalar         *vx, *vy, *vz;
	PetscInt            ii, sx, sy, sz, nx, ny, nz, I, J, K, II, JJ, KK, lrank, grank;
	PetscScalar         xb, yb, zb, xe, ye, ze, xc, yc, zc, *iter, *ncx, *ncy, *ncz, *ccx, *ccy, *ccz, ***lvx, ***lvy, ***lvz;

	PetscFunctionBegin;

	fs = jr->fs;

	VecGetArray(aop->vx,&vx);
	VecGetArray(aop->vy,&vy);
	VecGetArray(aop->vz,&vz);

 	ierr = VecDuplicate(jr->gsol, &pro);	 	CHKERRQ(ierr);

	// Access the local velocities
	ierr = DMDAVecGetArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_X, &gproX); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gproY); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gproZ); CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lproX); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lproY); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lproZ); CHKERRQ(ierr);

	// Zero out entries
	VecZeroEntries(gproX);
	VecZeroEntries(gproY);
	VecZeroEntries(gproZ);
	VecZeroEntries(lproX);
	VecZeroEntries(lproY);
	VecZeroEntries(lproZ);

	//=================
	// INDEXING LOOP
	//=================
	for(ii = 0; ii < user->AdjointNumInd; ii++)
	{
		// Create coordinate vector
		coord_local[0] = user->Adjoint_x[ii];
		coord_local[1] = user->Adjoint_y[ii];
		coord_local[2] = user->Adjoint_z[ii];

		// get global & local ranks of a marker
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
			I = FindPointInCell(ncx, 0, nx, coord_local[0]);
			J = FindPointInCell(ncy, 0, ny, coord_local[1]);
			K = FindPointInCell(ncz, 0, nz, coord_local[2]);

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

			// Vx (interpolate x velocity)
			if(user->AdjointVel[ii] == 1)
			{
				vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);

				// get relative coordinates
				xe = (coord_local[0] - ncx[I ])/(ncx[I +1] - ncx[I ]); xb = 1.0 - xe;

				llproX[sz+KK  ][sy+JJ  ][sx+I  ] = xb;
				llproX[sz+KK  ][sy+JJ  ][sx+I+1] = xe;
				llproX[sz+KK  ][sy+JJ+1][sx+I  ] = xb;
				llproX[sz+KK  ][sy+JJ+1][sx+I+1] = xe;
				llproX[sz+KK+1][sy+JJ  ][sx+I  ] = xb;
				llproX[sz+KK+1][sy+JJ  ][sx+I+1] = xe;
				llproX[sz+KK+1][sy+JJ+1][sx+I  ] = xb;
				llproX[sz+KK+1][sy+JJ+1][sx+I+1] = xe;
			}
			// Vy (interpolate y velocity)
			else if(user->AdjointVel[ii] == 2)
			{
				vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);

				// get relative coordinates
				ye = (coord_local[1] - ncy[J ])/(ncy[J +1] - ncy[J ]); yb = 1.0 - ye;

				llproY[sz+KK  ][sy+J  ][sx+II  ] = yb;
				llproY[sz+KK  ][sy+J  ][sx+II+1] = yb;
				llproY[sz+KK  ][sy+J+1][sx+II  ] = ye;
				llproY[sz+KK  ][sy+J+1][sx+II+1] = ye;
				llproY[sz+KK+1][sy+J  ][sx+II  ] = yb;
				llproY[sz+KK+1][sy+J  ][sx+II+1] = yb;
				llproY[sz+KK+1][sy+J+1][sx+II  ] = ye;
				llproY[sz+KK+1][sy+J+1][sx+II+1] = ye;
			}
			// Vz (interpolate z velocity)
			else if(user->AdjointVel[ii] == 3)
			{
				vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

				// get relative coordinates
				ze = (coord_local[2] - ncz[K ])/(ncz[K +1] - ncz[K ]); zb = 1.0 - ze;

				llproZ[sz+K  ][sy+JJ  ][sx+II  ] = zb;
				llproZ[sz+K  ][sy+JJ  ][sx+II+1] = zb;
				llproZ[sz+K  ][sy+JJ+1][sx+II  ] = zb;
				llproZ[sz+K  ][sy+JJ+1][sx+II+1] = zb;
				llproZ[sz+K+1][sy+JJ  ][sx+II  ] = ze;
				llproZ[sz+K+1][sy+JJ  ][sx+II+1] = ze;
				llproZ[sz+K+1][sy+JJ+1][sx+II  ] = ze;
				llproZ[sz+K+1][sy+JJ+1][sx+II+1] = ze;
			}
			else
			{
				PetscPrintf(PETSC_COMM_WORLD,"ERROR: Velocity direction is not defined ; Choose between [1-3]\n");
				PetscFunctionReturn(1);
			}
			ierr = DMDAVecRestoreArray(fs->DA_X, lproX, &llproX);      CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(fs->DA_Y, lproY, &llproY);      CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(fs->DA_Z, lproZ, &llproZ);      CHKERRQ(ierr);
		}
	}

	LOCAL_TO_GLOBAL(fs->DA_X, lproX, gproX);
	LOCAL_TO_GLOBAL(fs->DA_Y, lproY, gproY);
	LOCAL_TO_GLOBAL(fs->DA_Z, lproZ, gproZ);

	ierr = VecGetArray(gproX, &dggproX);      CHKERRQ(ierr);
	ierr = VecGetArray(gproY, &dggproY);      CHKERRQ(ierr);
	ierr = VecGetArray(gproZ, &dggproZ);      CHKERRQ(ierr);

	// Put the proportion into the Projection vector where the user defined the computation coordinates (P)
	ierr = VecGetArray(pro, &temppro);      CHKERRQ(ierr);
	iter = temppro;

	ierr  = PetscMemcpy(iter, dggproX, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, dggproY, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, dggproZ, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	// restore & destroy
	ierr = VecRestoreArray(pro, &temppro);         CHKERRQ(ierr);

	ierr = VecRestoreArray(gproX, &dggproX);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gproY, &dggproY);       CHKERRQ(ierr);
	ierr = VecRestoreArray(gproZ, &dggproZ);       CHKERRQ(ierr);

	ierr = VecRestoreArray(aop->vx,&vx); 	      CHKERRQ(ierr);
	ierr = VecRestoreArray(aop->vy,&vy); 	      CHKERRQ(ierr);
	ierr = VecRestoreArray(aop->vz,&vz); 	      CHKERRQ(ierr);

	ierr = VecCopy(pro,aop->pro);

	ierr = VecDestroy(&lproX);
	ierr = VecDestroy(&lproY);
	ierr = VecDestroy(&lproZ);
	ierr = VecDestroy(&gproX);
	ierr = VecDestroy(&gproY);
	ierr = VecDestroy(&gproZ);
	ierr = VecDestroy(&pro);

	ierr = DMDAVecRestoreArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointGradientPerturbParameter"
PetscErrorCode AdjointGradientPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop)
{
	PetscScalar ini, perturb;

	PetscFunctionBegin;

	// Get the perturbation value
	perturb = aop->Perturb;

	// Perturb the current parameter in the current phase (more to be included)
	if(CurPar==1)			// rho
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho;
		nl->pc->pm->jr->phases[CurPhase].rho +=  perturb;
	}
	else if (CurPar==2)	    // rho_n
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho_n;
		nl->pc->pm->jr->phases[CurPhase].rho_n +=  perturb;
	}
	else if (CurPar==3)	    // rho_c
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho_c;
		nl->pc->pm->jr->phases[CurPhase].rho_c +=  perturb;
	}
	else if (CurPar==4)	    // beta
	{
		ini = nl->pc->pm->jr->phases[CurPhase].beta;
		nl->pc->pm->jr->phases[CurPhase].beta +=  perturb;
	}
	else if (CurPar==5)	    // K
	{
		ini = nl->pc->pm->jr->phases[CurPhase].K;
		nl->pc->pm->jr->phases[CurPhase].K +=  perturb;
	}
	else if (CurPar==6)	    // Kp
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Kp;
		nl->pc->pm->jr->phases[CurPhase].Kp +=  perturb;
	}
	else if (CurPar==7)	    // G
	{
		ini = nl->pc->pm->jr->phases[CurPhase].G;
		nl->pc->pm->jr->phases[CurPhase].G +=  perturb;
	}
	else if (CurPar==8)	    // Bd
	{
		// This kind of perturbs the whole NEWTONIAN viscosity, consider perturbing the parameters directly
		ini = nl->pc->pm->jr->phases[CurPhase].Bd;
		PetscScalar BdTemp;
		BdTemp = (1.0/(2*ini)) + perturb;
		nl->pc->pm->jr->phases[CurPhase].Bd =  (1.0/(2*BdTemp));
	}
	else if (CurPar==9)	    // Ed
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Ed;
		nl->pc->pm->jr->phases[CurPhase].Ed +=  perturb;
	}
	else if (CurPar==10)	// Vd
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Vd;
		nl->pc->pm->jr->phases[CurPhase].Vd +=  perturb;
	}
	else if (CurPar==11)	// Bn
	{
		// This kind of perturbs the whole DISLOCATION viscosity, consider perturbing the parameters directly
		ini = nl->pc->pm->jr->phases[CurPhase].Bn;
		PetscScalar BnTemp;
		BnTemp = (1.0/(2*ini)) + perturb;
		nl->pc->pm->jr->phases[CurPhase].Bn =  (1.0/(2*BnTemp));
	}
	else if (CurPar==12)	// n
	{
		ini = nl->pc->pm->jr->phases[CurPhase].n;
		nl->pc->pm->jr->phases[CurPhase].n +=  perturb;
	}
	else if (CurPar==13)	// En
	{
		ini = nl->pc->pm->jr->phases[CurPhase].En;
		nl->pc->pm->jr->phases[CurPhase].En +=  perturb;
	}
	else if (CurPar==14)	// Vn
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Vn;
		nl->pc->pm->jr->phases[CurPhase].Vn +=  perturb;
	}
	else if (CurPar==15)	// taup
	{
		ini = nl->pc->pm->jr->phases[CurPhase].taup;
		nl->pc->pm->jr->phases[CurPhase].taup +=  perturb;
	}
	else if (CurPar==16)	// gamma
	{
		ini = nl->pc->pm->jr->phases[CurPhase].gamma;
		nl->pc->pm->jr->phases[CurPhase].gamma +=  perturb;
	}
	else if (CurPar==17)	// q
	{
		ini = nl->pc->pm->jr->phases[CurPhase].q;
		nl->pc->pm->jr->phases[CurPhase].q +=  perturb;
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"ERROR: Definition of the current parameter is not defined ; Choose between [1-17]\n");
		PetscFunctionReturn(1);
	}

	// Store initial value of current parameter
	aop->Ini = ini;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointGradientResetParameter"
PetscErrorCode AdjointGradientResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop)
{
	PetscScalar ini;

	PetscFunctionBegin;

	// Get initial value of currently perturbed parameter
	ini = aop->Ini;

	// Set the current parameter back to its original value
	if (CurPar==1)        { nl->pc->pm->jr->phases[CurPhase].rho   = ini;
	}else if (CurPar==2)  { nl->pc->pm->jr->phases[CurPhase].rho_n = ini;
	}else if (CurPar==3)  { nl->pc->pm->jr->phases[CurPhase].rho_c = ini;
	}else if (CurPar==4)  { nl->pc->pm->jr->phases[CurPhase].beta  = ini;
	}else if (CurPar==5)  { nl->pc->pm->jr->phases[CurPhase].K     = ini;
	}else if (CurPar==6)  { nl->pc->pm->jr->phases[CurPhase].Kp    = ini;
	}else if (CurPar==7)  { nl->pc->pm->jr->phases[CurPhase].G     = ini;
	}else if (CurPar==8)  { nl->pc->pm->jr->phases[CurPhase].Bd    = ini;
	}else if (CurPar==9)  { nl->pc->pm->jr->phases[CurPhase].Ed    = ini;
	}else if (CurPar==10) { nl->pc->pm->jr->phases[CurPhase].Vd    = ini;
	}else if (CurPar==11) { nl->pc->pm->jr->phases[CurPhase].Bn    = ini;
	}else if (CurPar==12) { nl->pc->pm->jr->phases[CurPhase].n     = ini;
	}else if (CurPar==13) { nl->pc->pm->jr->phases[CurPhase].En    = ini;
	}else if (CurPar==14) { nl->pc->pm->jr->phases[CurPhase].Vn    = ini;
	}else if (CurPar==15) { nl->pc->pm->jr->phases[CurPhase].taup  = ini;
	}else if (CurPar==16) { nl->pc->pm->jr->phases[CurPhase].gamma = ini;
	}else if (CurPar==17) { nl->pc->pm->jr->phases[CurPhase].q     = ini;}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointLoadCompareData"
PetscErrorCode AdjointLoadCompareData(void *ctx)
{
	PetscFunctionBegin;

	PetscErrorCode  ierrp;
	PetscViewer     viewer;
	JacRes 			*jr;
	AdjGrad 		*aop;

 	aop = (AdjGrad*)ctx;

 	jr = aop->jr;

	VecDuplicate(jr->gsol, &aop->xini);

	// Load the comparison solution vector
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Forward_Solution.bin",FILE_MODE_READ,&viewer);
	ierrp = VecLoad(aop->xini,viewer);       CHKERRQ(ierrp);

	if (ierrp)
	{
		PetscPrintf(PETSC_COMM_WORLD,"ERROR: Could not load the initial solution (xini)\n");
		PetscFunctionReturn(1);
	}

	// Destroy
	PetscViewerDestroy(&viewer);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointExtractParameterVec"
PetscErrorCode AdjointExtractParameterVec(NLSol *nl,Vec P, UserCtx *user)
{
	PetscInt       j, CurPhase, CurPar;
	PetscScalar    *PTemp;

	PetscFunctionBegin;

	VecGetArray(P,&PTemp);

	for(j = 0; j < user->AdjointNumPar; j++)
	{
		// Get current phase and parameter
		CurPhase = user->AdjointPhases[j];
		CurPar   = user->AdjointParameters[j];

		// Get the current parameter
		if (CurPar==1)        {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].rho;
		}else if (CurPar==2)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].rho_n;
		}else if (CurPar==3)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].rho_c;
		}else if (CurPar==4)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].beta;
		}else if (CurPar==5)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].K;
		}else if (CurPar==6)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].Kp;
		}else if (CurPar==7)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].G;
		}else if (CurPar==8)  {PTemp[j] = 1/(2*nl->pc->pm->jr->phases[CurPhase].Bd);
		}else if (CurPar==9)  {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].Ed;
		}else if (CurPar==10) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].Vd;
		}else if (CurPar==11) {PTemp[j] = 1/(2*nl->pc->pm->jr->phases[CurPhase].Bn);
		}else if (CurPar==12) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].n;
		}else if (CurPar==13) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].En;
		}else if (CurPar==14) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].Vn;
		}else if (CurPar==15) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].taup;
		}else if (CurPar==16) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].gamma;
		}else if (CurPar==17) {PTemp[j] = nl->pc->pm->jr->phases[CurPhase].q;}

		PetscPrintf(PETSC_COMM_WORLD,"%D. Initial Parameter = %.12f ; CurPar = %d ; CurPhase = %d\n",j,PTemp[j],CurPar,CurPhase);
	}

	VecRestoreArray(P,&PTemp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPerturbParameterVec"
PetscErrorCode AdjointPerturbParameterVec(NLSol *nl,Vec P, UserCtx *user)
{
	PetscInt       j, CurPhase, CurPar;
	PetscScalar    *PTemp;

	PetscFunctionBegin;

	VecGetArray(P,&PTemp);

	for(j = 0; j < user->AdjointNumPar; j++)
	{
		// Get current phase and parameter
		CurPhase = user->AdjointPhases[j];
		CurPar   = user->AdjointParameters[j];

		// Set the current parameter back to its original value
		if (CurPar==1)        {nl->pc->pm->jr->phases[CurPhase].rho 	= PTemp[j];
		}else if (CurPar==2)  {nl->pc->pm->jr->phases[CurPhase].rho_n 	= PTemp[j];
		}else if (CurPar==3)  {nl->pc->pm->jr->phases[CurPhase].rho_c 	= PTemp[j];
		}else if (CurPar==4)  {nl->pc->pm->jr->phases[CurPhase].beta 	= PTemp[j];
		}else if (CurPar==5)  {nl->pc->pm->jr->phases[CurPhase].K 		= PTemp[j];
		}else if (CurPar==6)  {nl->pc->pm->jr->phases[CurPhase].Kp 		= PTemp[j];
		}else if (CurPar==7)  {nl->pc->pm->jr->phases[CurPhase].G 		= PTemp[j];
		}else if (CurPar==8)  {nl->pc->pm->jr->phases[CurPhase].Bd 		= 1/(2*PTemp[j]);
		}else if (CurPar==9)  {nl->pc->pm->jr->phases[CurPhase].Ed 		= PTemp[j];
		}else if (CurPar==10) {nl->pc->pm->jr->phases[CurPhase].Vd 		= PTemp[j];
		}else if (CurPar==11) {nl->pc->pm->jr->phases[CurPhase].Bn 		= 1/(2*PTemp[j]);
		}else if (CurPar==12) {nl->pc->pm->jr->phases[CurPhase].n 		= PTemp[j];
		}else if (CurPar==13) {nl->pc->pm->jr->phases[CurPhase].En 		= PTemp[j];
		}else if (CurPar==14) {nl->pc->pm->jr->phases[CurPhase].Vn 		= PTemp[j];
		}else if (CurPar==15) {nl->pc->pm->jr->phases[CurPhase].taup 	= PTemp[j];
		}else if (CurPar==16) {nl->pc->pm->jr->phases[CurPhase].gamma 	= PTemp[j];
		}else if (CurPar==17) {nl->pc->pm->jr->phases[CurPhase].q 		= PTemp[j];}

		PetscPrintf(PETSC_COMM_WORLD,"%D. Current Parameter = %.12f ; CurPar = %d ; CurPhase = %d\n",j,PTemp[j],CurPar,CurPhase);
	}

	VecRestoreArray(P,&PTemp);

	PetscFunctionReturn(0);
}
