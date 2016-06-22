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
 #define __FUNCT__ "AdjointObjectiveAndGradientFunction"
 PetscErrorCode AdjointObjectiveAndGradientFunction(AdjGrad *aop, JacRes *jr, NLSol *nl, ModParam *IOparam, SNES snes)
 {
 	PetscErrorCode ierr;
 	PetscFunctionBegin;

	Vec         	x;

	// Create everything
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vx); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vy); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vz); CHKERRQ(ierr);

 	ierr = VecDuplicate(jr->gsol, &aop->dF);	 									CHKERRQ(ierr);

	//========================================
	// COMPUTE OBJECTIVE FUNCTION & GRADIENT
	//========================================
 	// only compute the gradients (F = P*x (unused) // dF/dx = P)
 	if (IOparam->use == 2)
 	{
 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);	 	CHKERRQ(ierr);

 	 	// Get the solution
 	 	ierr = VecDuplicate(jr->gsol, &x);	 			CHKERRQ(ierr);
 	 	ierr = VecCopy(jr->gsol,x); 					CHKERRQ(ierr);

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, aop, IOparam); 	CHKERRQ(ierr);

 		ierr = VecCopy(aop->pro,aop->dF); 				CHKERRQ(ierr);

 		// Get the gradients
 		ierr = AdjointComputeGradients(jr, aop, nl, snes, IOparam); 		CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint inversion
 	else if(IOparam->use == 3)
	{
 		Vec xini;
 		PetscScalar Ad;

 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);	 			CHKERRQ(ierr);

 	 	if(IOparam->count == 1)
 	 	{
 	 		// We need to read the comparison solution
 	 	 	PetscErrorCode  ierrp;
 	 	 	PetscViewer     viewer;

 	 	 	ierr = VecDuplicate(jr->gsol, &IOparam->xini);	 	CHKERRQ(ierr);

 	 	 	// Load the comparison solution vector
 	 	 	PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Forward_Solution.bin",FILE_MODE_READ,&viewer);
 	 	 	ierrp = VecLoad(IOparam->xini,viewer);       		CHKERRQ(ierrp);

 	 	 	if (ierrp){PetscPrintf(PETSC_COMM_WORLD,"ERROR: Could not load the initial solution (xini)\n");PetscFunctionReturn(1);}

 	 	 	// Destroy
 	 	 	PetscViewerDestroy(&viewer);
 	 	}

 	 	// Get the solution
 	 	ierr = VecDuplicate(jr->gsol, &x);	 	CHKERRQ(ierr);
 	 	ierr = VecCopy(jr->gsol,x); 			CHKERRQ(ierr);

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, aop, IOparam); 		CHKERRQ(ierr);

 		PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE COST FUNCTION\n******************************************\n");

	 	// Copy temporary comparison solution
	 	ierr = VecDuplicate(jr->gsol, &xini);	 		CHKERRQ(ierr);
	 	ierr = VecCopy(IOparam->xini,xini); 			CHKERRQ(ierr);

 		// Incorporate projection vector (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr = VecAYPX(xini,-1,jr->gsol);        		CHKERRQ(ierr);
 		ierr = VecPointwiseMult(xini, xini,aop->pro);	CHKERRQ(ierr);

 		// Compute objective function value (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr 	           = VecDot(xini,xini,&Ad);
 		Ad 		          /= 2;
 		IOparam->mfit 	   = Ad;

 	 	PetscPrintf(PETSC_COMM_WORLD,"Current Cost function = %.12f\n",IOparam->mfit);

 		// Compute it's derivative (dF/dx = P*x-P*x_ini)
 		ierr = VecCopy(xini,aop->dF); 		CHKERRQ(ierr);

 		// Get the gradients
 		ierr = AdjointComputeGradients(jr, aop, nl, snes, IOparam); 		CHKERRQ(ierr);

 		// Destroy & increase count
 		ierr = VecDestroy(&xini);
	}
 	else
 	{
 	 	PetscPrintf(PETSC_COMM_WORLD,"ERROR: ComputeAdjointGradient value is not defined ; Choose between [1-2]\n");
 	 	PetscFunctionReturn(1);
 	}

 	ierr = VecDestroy(&x);

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

 	PetscFunctionReturn(0);
 }*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointComputeGradients"
PetscErrorCode AdjointComputeGradients(JacRes *jr, AdjGrad *aop, NLSol *nl, SNES snes, ModParam *IOparam)
{
	PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE GRADIENTS\n******************************************\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG              *fs;
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, j, CurPhase, CurPar, lrank, grank;
	PetscScalar         grd, Perturb, coord_local[3], *vx, *vy, *vz;
	Vec 				rpl, sol, psi, drdp, res, Perturb_vec;
	PC                  ipc;
	Scaling             *scal;

	fs = jr->fs;

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


	//=================
	// PARAMETER LOOP
	//=================
	for(j = 0; j < IOparam->mdN; j++)
	{
		// Get residual since it is overwritten in VecAYPX
		ierr = VecDuplicate(jr->gres, &res);  	CHKERRQ(ierr);
		ierr = VecCopy(jr->gres,res); 			CHKERRQ(ierr);

		// Get current phase and parameter
		CurPhase = IOparam->phs[j];
		CurPar   = IOparam->typ[j];

		// Perturb the current parameter in the current phase
		ierr = AdjointGradientPerturbParameter(nl, CurPar, CurPhase, aop, scal);      CHKERRQ(ierr);

		// get the actual used perturbation parameter which is 1e-6*parameter
		Perturb = aop->Perturb;
		ierr = VecDuplicate(jr->gsol, &Perturb_vec); 		CHKERRQ(ierr);
		ierr = VecSet(Perturb_vec,Perturb);			 		CHKERRQ(ierr);

		// Compute residual with the converged Jacobian (dr/dp = [r(p+h) - r(p)]/h)
		ierr = FormResidual(snes, sol, rpl, nl); 	       	CHKERRQ(ierr);
		ierr = VecAYPX(res,-1,rpl);        					CHKERRQ(ierr);
		ierr = VecPointwiseDivide(drdp,res,Perturb_vec);   	CHKERRQ(ierr);

		// Compute the gradient (dF/dp = -psi^T * dr/dp) & Save gradient
		ierr = VecDot(drdp,psi,&grd);     					CHKERRQ(ierr);
		IOparam->grd[j] 	= -1 * grd;

		// Reset perturbed parameter
		ierr = AdjointGradientResetParameter(nl, CurPar, CurPhase, aop);       CHKERRQ(ierr);

		// Print result
		PetscPrintf(PETSC_COMM_WORLD,"%D.Gradient (dimensional) = %.40f ; CurPar = %d ; CurPhase = %d\n",j+1, IOparam->grd[j]*aop->CurScal, CurPar, CurPhase);

		// Destroy overwritten residual vector
		ierr = VecDestroy(&res);
	}

	VecGetArray(aop->vx,&vx);
	VecGetArray(aop->vy,&vy);
	VecGetArray(aop->vz,&vz);

	if(IOparam->mdI<10)
	{
		// get the current velocities at comparison point
		ierr = AdjointPointInPro(jr, aop, IOparam); 		CHKERRQ(ierr);

		// Print the solution variable at the user defined index (if they are suffently few)
		for (i=0; i<IOparam->mdI; i++)
		{

			coord_local[0] = IOparam->Ax[i];
			coord_local[1] = IOparam->Ay[i];
			coord_local[2] = IOparam->Az[i];

			// get global & local ranks of a marker
			ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

			// If lrank is not 13 the point is not on this processor
			if(lrank == 13)
			{

				if (IOparam->Av[i] == 1)
				{
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vx] (dimensional) = %.12f\n",vx[i]*scal->velocity);
				}
				else if (IOparam->Av[i] == 2)
				{
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vy] (dimensional) = %.12f\n",vy[i]*scal->velocity);
				}
				else if (IOparam->Av[i] == 3)
				{
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vz] (dimensional) = %.12f\n",vz[i]*scal->velocity);
				}
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
	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful & took %g s\n******************************************\n",cputime_end - cputime_start);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPointInPro"
PetscErrorCode AdjointPointInPro(JacRes *jr, AdjGrad *aop, ModParam *IOparam)
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

 	ierr = VecDuplicate(jr->gsol, &pro);	 		 CHKERRQ(ierr);

	// Access the local velocities
	ierr = DMDAVecGetArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_X, &gproX); 	 CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &gproY); 	 CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &gproZ); 	 CHKERRQ(ierr);

	ierr = DMCreateLocalVector (fs->DA_X, &lproX); 	 CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &lproY); 	 CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &lproZ); 	 CHKERRQ(ierr);

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
	for(ii = 0; ii < IOparam->mdI; ii++)
	{
		// Create coordinate vector
		coord_local[0] = IOparam->Ax[ii];
		coord_local[1] = IOparam->Ay[ii];
		coord_local[2] = IOparam->Az[ii];

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
			if(IOparam->Av[ii] == 1)
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
			else if(IOparam->Av[ii] == 2)
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
			else if(IOparam->Av[ii] == 3)
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
	ierr = VecGetArray(pro, &temppro);        CHKERRQ(ierr);
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

	ierr = VecRestoreArray(aop->vx,&vx); 	       CHKERRQ(ierr);
	ierr = VecRestoreArray(aop->vy,&vy); 	       CHKERRQ(ierr);
	ierr = VecRestoreArray(aop->vz,&vz); 	       CHKERRQ(ierr);

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
PetscErrorCode AdjointGradientPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop, Scaling *scal)
{
	PetscScalar         ini, perturb, curscal;

	PetscFunctionBegin;

	// Get the perturbation value & scaling
	perturb = aop->Perturb;
	curscal = aop->CurScal;

	// Perturb the current parameter in the current phase (more to be included)
	if(CurPar==_RHO0_)			// rho
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho;
		nl->pc->pm->jr->phases[CurPhase].rho +=  perturb;
		curscal = (scal->velocity)/(scal->density);
	}
	else if (CurPar==_RHON_)	    // rho_n
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho_n;
		nl->pc->pm->jr->phases[CurPhase].rho_n +=  perturb;
		curscal = (scal->velocity)/1;
	}
	else if (CurPar==_RHOC_)	    // rho_c
	{
		ini = nl->pc->pm->jr->phases[CurPhase].rho_c;
		nl->pc->pm->jr->phases[CurPhase].rho_c +=  perturb;
		curscal = (scal->velocity)*(scal->length_si);
	}
	else if (CurPar==_K_)	    // K
	{
		ini = nl->pc->pm->jr->phases[CurPhase].K;
		nl->pc->pm->jr->phases[CurPhase].K +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_KP_)	    // Kp
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Kp;
		nl->pc->pm->jr->phases[CurPhase].Kp +=  perturb;
		curscal = (scal->velocity)/1;
	}
	else if (CurPar==_SHEAR_)	    // G
	{
		ini = nl->pc->pm->jr->phases[CurPhase].G;
		nl->pc->pm->jr->phases[CurPhase].G +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_ETA_)	    // Bd
	{
		// This kind of perturbs the whole NEWTONIAN viscosity, consider perturbing the parameters directly
		ini = nl->pc->pm->jr->phases[CurPhase].Bd;
		PetscScalar BdTemp;
		BdTemp = (1.0/(2*ini)) + perturb;
		nl->pc->pm->jr->phases[CurPhase].Bd =  (1.0/(2*BdTemp));
		curscal = (scal->velocity)/(scal->viscosity);
	}
	else if (CurPar==_ED_)	    // Ed
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Ed;
		nl->pc->pm->jr->phases[CurPhase].Ed +=  perturb;
		curscal = (scal->velocity)/(1);   // Not sure
	}
	else if (CurPar==_VD_)	// Vd
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Vd;
		nl->pc->pm->jr->phases[CurPhase].Vd +=  perturb;
		curscal = (scal->velocity)*(scal->stress_si);
	}
	else if (CurPar==_ETA0_)	// Bn
	{
		// This kind of perturbs the whole DISLOCATION viscosity, consider perturbing the parameters directly
		ini = nl->pc->pm->jr->phases[CurPhase].Bn;
		PetscScalar BnTemp;
		BnTemp = (1.0/(2*ini)) + perturb;
		nl->pc->pm->jr->phases[CurPhase].Bn =  (1.0/(2*BnTemp));
		curscal = (scal->velocity)/(scal->viscosity);
	}
	else if (CurPar==_N_)	// n
	{
		ini = nl->pc->pm->jr->phases[CurPhase].n;
		nl->pc->pm->jr->phases[CurPhase].n +=  perturb;
		curscal = (scal->velocity)/(1);
	}
	else if (CurPar==_EN_)	// En
	{
		ini = nl->pc->pm->jr->phases[CurPhase].En;
		nl->pc->pm->jr->phases[CurPhase].En +=  perturb;
		curscal = (scal->velocity)/(1);    // Not sure
	}
	else if (CurPar==_VN_)	// Vn
	{
		ini = nl->pc->pm->jr->phases[CurPhase].Vn;
		nl->pc->pm->jr->phases[CurPhase].Vn +=  perturb;
		curscal = (scal->velocity)*(scal->stress_si);
	}
	else if (CurPar==_TAUP_)	// taup
	{
		ini = nl->pc->pm->jr->phases[CurPhase].taup;
		nl->pc->pm->jr->phases[CurPhase].taup +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_GAMMA_)	// gamma
	{
		ini = nl->pc->pm->jr->phases[CurPhase].gamma;
		nl->pc->pm->jr->phases[CurPhase].gamma +=  perturb;
		curscal = (scal->velocity)/(1);
	}
	else if (CurPar==_Q_)	// q
	{
		ini = nl->pc->pm->jr->phases[CurPhase].q;
		nl->pc->pm->jr->phases[CurPhase].q +=  perturb;
		curscal = (scal->velocity)/(1);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"ERROR: Definition of the current parameter is not defined ; Choose between [0-15]\n");
		PetscFunctionReturn(1);
	}

	// Store initial value of current parameter & scaling
	aop->Ini 		= ini;
	aop->CurScal 	= curscal;

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
	if (CurPar==_RHO0_)        	{nl->pc->pm->jr->phases[CurPhase].rho   = ini;
	}else if (CurPar==_RHON_)  	{nl->pc->pm->jr->phases[CurPhase].rho_n = ini;
	}else if (CurPar==_RHOC_)  	{nl->pc->pm->jr->phases[CurPhase].rho_c = ini;
	}else if (CurPar==_K_)  	{nl->pc->pm->jr->phases[CurPhase].K     = ini;
	}else if (CurPar==_KP_)  	{nl->pc->pm->jr->phases[CurPhase].Kp    = ini;
	}else if (CurPar==_SHEAR_)  {nl->pc->pm->jr->phases[CurPhase].G     = ini;
	}else if (CurPar==_ETA_)  	{nl->pc->pm->jr->phases[CurPhase].Bd    = ini;
	}else if (CurPar==_ED_)  	{nl->pc->pm->jr->phases[CurPhase].Ed    = ini;
	}else if (CurPar==_VD_) 	{nl->pc->pm->jr->phases[CurPhase].Vd    = ini;
	}else if (CurPar==_ETA0_) 	{nl->pc->pm->jr->phases[CurPhase].Bn    = ini;
	}else if (CurPar==_N_) 		{nl->pc->pm->jr->phases[CurPhase].n     = ini;
	}else if (CurPar==_EN_) 	{nl->pc->pm->jr->phases[CurPhase].En    = ini;
	}else if (CurPar==_VN_) 	{nl->pc->pm->jr->phases[CurPhase].Vn    = ini;
	}else if (CurPar==_TAUP_) 	{nl->pc->pm->jr->phases[CurPhase].taup  = ini;
	}else if (CurPar==_GAMMA_) 	{nl->pc->pm->jr->phases[CurPhase].gamma = ini;
	}else if (CurPar==_Q_) 		{nl->pc->pm->jr->phases[CurPhase].q     = ini;}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
