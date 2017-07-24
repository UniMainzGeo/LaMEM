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
#include "objFunct.h"
#include "parsing.h"
#include "constEq.h"
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
 PetscErrorCode AdjointObjectiveAndGradientFunction(AdjGrad *aop, JacRes *jr, NLSol *nl, ModParam *IOparam, SNES snes, FreeSurf *surf)
 {

	Scaling             *scal;

 	PetscErrorCode ierr;
 	PetscFunctionBegin;

 	scal = jr->scal;

 	ierr = VecDuplicate(jr->gsol, &aop->dF);                                     CHKERRQ(ierr);

	//========================================
	// COMPUTE OBJECTIVE FUNCTION & GRADIENT
	//========================================
 	// only compute the gradients (F = P*x (unused) // dF/dx = P)
 	if (IOparam->use == 2)
 	{
 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);	 	CHKERRQ(ierr);

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, aop, IOparam, surf); 	CHKERRQ(ierr);

 		ierr = VecCopy(aop->pro,aop->dF); 				CHKERRQ(ierr);

 		// Get the gradients
 		ierr = AdjointComputeGradients(jr, aop, nl, snes, IOparam, surf);        CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint inversion
 	else if(IOparam->use == 3)
	{
 		Vec xini;
 		PetscScalar Ad;

 	 	// Create projection vector
 	 	ierr = VecDuplicate(jr->gsol, &aop->pro);                                CHKERRQ(ierr);

 	 	if(IOparam->count == 1)
 	 	{
 	 		// We need to read the comparison solution
 	 	 	PetscErrorCode  ierrp;
 	 	 	PetscViewer     viewer;

 	 	 	ierr = VecDuplicate(jr->gsol, &IOparam->xini);                       CHKERRQ(ierr);

 	 	 	// Load the comparison solution vector
 	 	 	PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Forward_Solution.bin",FILE_MODE_READ,&viewer);
 	 	 	ierrp = VecLoad(IOparam->xini,viewer);                               CHKERRQ(ierrp);

 	 	 	if (ierrp){PetscPrintf(PETSC_COMM_WORLD,"ADJOINT ERROR: Could not load the initial solution (xini)\n");PetscFunctionReturn(1);}

 	 	 	// Destroy
 	 	 	PetscViewerDestroy(&viewer);
 	 	}

 		// Put the proportion into the Projection vector where the user defined the computation coordinates (P) & get the velocities
 		ierr = AdjointPointInPro(jr, aop, IOparam, surf);                       CHKERRQ(ierr);

 		PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE COST FUNCTION\n******************************************\n");

	 	// Copy temporary comparison solution
	 	ierr = VecDuplicate(jr->gsol, &xini);                                   CHKERRQ(ierr);
	 	ierr = VecCopy(IOparam->xini,xini);                                     CHKERRQ(ierr);

 		// Incorporate projection vector (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr = VecAYPX(xini,-1,jr->gsol);                                       CHKERRQ(ierr);
 		ierr = VecPointwiseMult(xini, xini,aop->pro);                           CHKERRQ(ierr);

 		// Compute objective function value (F = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)])
 		ierr 	           = VecDot(xini,xini,&Ad);
 		Ad 		          /= 2;
 		IOparam->mfit 	   = Ad*pow(scal->velocity,2); // Dimensional misfit function

 		// Perform Tikhonov regularization (TN)
 		if (IOparam->reg==1)
 		{
 			PetscInt       i;
 			PetscScalar    R, *p;

 			VecGetArray(IOparam->P,&p);
 			for (i==0 ; i<IOparam->mdN ; i++)
 			{
 				R += (IOparam->W[i]/2) * (p[i]*p[i]);
 			}
 			VecRestoreArray(IOparam->P,&p);
 			IOparam->mfit += R;
 		}
 		else if (IOparam->reg==2)
 		{
 			PetscInt       i;
 			PetscScalar    R, *p;

 			VecGetArray(IOparam->P,&p);
 			for (i==0 ; i<IOparam->mdN ; i++)
 			{
 				R += (IOparam->W[i]/2) * sqrt(p[i]*p[i] + 100);
 			}
 			VecRestoreArray(IOparam->P,&p);
 			IOparam->mfit += R;
 		}

 	 	PetscPrintf(PETSC_COMM_WORLD,"Current Cost function = %.20f\n",IOparam->mfit);

 		// Compute it's derivative (dF/dx = P*x-P*x_ini)
 		ierr = VecCopy(xini,aop->dF); 		CHKERRQ(ierr);

 		// Get the gradients
 		ierr = AdjointComputeGradients(jr, aop, nl, snes, IOparam, surf);        CHKERRQ(ierr);

 		// Destroy
 		ierr = VecDestroy(&xini);
	}
 	else
 	{
 	 	PetscPrintf(PETSC_COMM_WORLD,"ADJOINT ERROR: ComputeAdjointGradient value is not defined ; Choose between [1-2]\n");
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

	// You should use lmvm (or the bounded version blmvm) tao_solve because this is using the bfgs approximated Hessian which is reasonable.

 	PetscFunctionReturn(0);
 }*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointComputeGradients"
PetscErrorCode AdjointComputeGradients(JacRes *jr, AdjGrad *aop, NLSol *nl, SNES snes, ModParam *IOparam, FreeSurf *surf)
{
	PetscPrintf(PETSC_COMM_WORLD,"******************************************\n      COMPUTATION OF THE GRADIENTS\n******************************************\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG              *fs;
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, j, CurPhase, CurPar, lrank, grank, fd;
	PetscScalar         dt, grd, Perturb, coord_local[3], *vx, *vy, *vz;
	Vec 				rpl, sol, psi, drdp, res, Perturb_vec;
	PC                  ipc;
	Scaling             *scal;

	fs = jr->fs;
	dt = jr->ts->dt;
	fd = 0;          // Counts FD approximations

	// Profile time
	PetscLogDouble     cputime_start, cputime_end;
	PetscTime(&cputime_start);

	scal = jr->scal;

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
	ierr = KSPSolve(ksp,aop->dF,psi);	CHKERRQ(ierr);
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
		
		// Compute residual with the converged Jacobian analytically
		if (CurPar==_RHO0_)
		{
			ierr = AdjointFormResidual(snes, sol, drdp, nl, CurPar, CurPhase);          CHKERRQ(ierr);
			aop->CurScal = (scal->velocity)/(scal->density);
		}
		else if(CurPar==_ETA0_)  // FD and adjoint can also compute the gradient for BN (look into the perturbation/residual routine later this file)
		{
			ierr = AdjointFormResidual(snes, sol, drdp, nl, CurPar, CurPhase);          CHKERRQ(ierr);
			aop->CurScal = (scal->velocity)/(scal->viscosity);
		}
		else if(CurPar==_N_)
		{
			ierr = AdjointFormResidual(snes, sol, drdp, nl, CurPar, CurPhase);          CHKERRQ(ierr);
			aop->CurScal = (scal->velocity)/(1);
		}
		else if(CurPar==_EN_)
		{
			ierr = AdjointFormResidual(snes, sol, drdp, nl, CurPar, CurPhase);          CHKERRQ(ierr);
			aop->CurScal = (scal->velocity)/(1);
		}
		else if(CurPar==_MFR_)
		{
			ierr = AdjointFormResidual(snes, sol, drdp, nl, CurPar, CurPhase);          CHKERRQ(ierr);
			aop->CurScal = (scal->velocity)/(1);
		}
		else  // Use FD approximation (dr/dp = [r(p+h) - r(p)]/h)
		{
			PetscPrintf(PETSC_COMM_WORLD,"    Adjoint Warning: No analytical residual computation known for %d (FD used)\n",CurPar);
			
			// Set perturbation paramter for the finite differences
			aop->Perturb = 1e-6;

			// Perturb the current parameter in the current phase
			ierr = AdjointGradientPerturbParameter(nl, CurPar, CurPhase, aop, scal);   CHKERRQ(ierr);

			// get the actual used perturbation parameter which is 1e-6*parameter
			Perturb = aop->Perturb;
			ierr = VecDuplicate(jr->gsol, &Perturb_vec);          CHKERRQ(ierr);
			ierr = VecSet(Perturb_vec,Perturb);                   CHKERRQ(ierr);
			
			ierr = FormResidual(snes, sol, rpl, nl);              CHKERRQ(ierr);
			ierr = VecAYPX(res,-1,rpl);                           CHKERRQ(ierr);
			ierr = VecPointwiseDivide(drdp,res,Perturb_vec);      CHKERRQ(ierr);
			
			// Reset perturbed parameter
			ierr = AdjointGradientResetParameter(nl, CurPar, CurPhase, aop);           CHKERRQ(ierr);
			
			fd += fd;  // Needed to free vectors later
		}
		
		// Compute the gradient (dF/dp = -psi^T * dr/dp) & Save gradient
		ierr = VecDot(drdp,psi,&grd);                       CHKERRQ(ierr);
		IOparam->grd[j] 	= -grd*aop->CurScal;            CHKERRQ(ierr);

		// Print result
		PetscPrintf(PETSC_COMM_WORLD,"%D.Gradient (dimensional) = %.40f ; CurPar = %d ; CurPhase = %d\n",j+1, IOparam->grd[j], CurPar, CurPhase);

		// Destroy overwritten residual vector
		ierr = VecDestroy(&res);
	}

	if(IOparam->mdI<11 && IOparam->Ap == 1)
	{
		VecGetArray(aop->vx,&vx);
		VecGetArray(aop->vy,&vy);
		VecGetArray(aop->vz,&vz);

		// get the current velocities at comparison point
		ierr = AdjointPointInPro(jr, aop, IOparam, surf);    CHKERRQ(ierr);

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
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vx] (dimensional) = %.30f at Location x = %g , y = %g , z = %g\n",vx[i]*scal->velocity,IOparam->Ax[i]*scal->length,IOparam->Ay[i]*scal->length,IOparam->Az[i]*scal->length);
				}
				else if (IOparam->Av[i] == 2)
				{
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vy] (dimensional) = %.30f at Location x = %g , y = %g , z = %g\n",vy[i]*scal->velocity,IOparam->Ax[i]*scal->length,IOparam->Ay[i]*scal->length,IOparam->Az[i]*scal->length);
				}
				else if (IOparam->Av[i] == 3)
				{
					PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vz] (dimensional) = %.30f\n",vz[i]*scal->velocity);
					// PetscPrintf(PETSC_COMM_SELF,"Computation variable [Vz] (dimensional) = %.30f at Location x = %g , y = %g , z = %g\n",vz[i]*scal->velocity,IOparam->Ax[i]*scal->length,IOparam->Ay[i]*scal->length,IOparam->Az[i]*scal->length);
				}

				if (IOparam->Adv == 1)     // advect the point?
				{
					IOparam->Ax[i] += vx[i]*dt;
					IOparam->Ay[i] += vy[i]*dt;
					IOparam->Az[i] += vz[i]*dt;
				}
				else
				{
					// simply don't advect
				}
			}
		}
		VecRestoreArray(aop->vx,&vx);
		VecRestoreArray(aop->vy,&vy);
		VecRestoreArray(aop->vz,&vz);
	}

	// Clean
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&drdp);
	ierr = VecDestroy(&rpl);
	if(fd > 0)
	{
		ierr = VecDestroy(&Perturb_vec);
	}

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful & took %g s\n******************************************\n",cputime_end - cputime_start);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPointInPro"
PetscErrorCode AdjointPointInPro(JacRes *jr, AdjGrad *aop, ModParam *IOparam, FreeSurf *surf)
{
	PetscErrorCode      ierr;
	FDSTAG              *fs;
	Vec                 lproX, lproY, lproZ, gproX, gproY, gproZ, pro;
	PetscScalar         coord_local[3], *temppro, ***llproX, ***llproY, ***llproZ, *dggproX, *dggproY, *dggproZ;
	PetscScalar         *vx, *vy, *vz;
	PetscInt            j, i, z, w, ii, sx, sy, sz, nx, ny, nz, I, J, K, II, JJ, KK, lrank, grank, p, level;
	PetscScalar         xb, yb, zb, xe, ye, ze, xc, yc, zc, *iter, *ncx, *ncy, *ncz, *ccx, *ccy, *ccz, ***lvx, ***lvy, ***lvz, ***vgrid, ***topo, ***vsurf;
	PetscReal           val;
	Discret1D           *dsz;
	InterpFlags         iflag;

	PetscFunctionBegin;

	fs = jr->fs;

 	ierr = VecDuplicate(jr->gsol, &pro);             CHKERRQ(ierr);

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

	// Zero out entries
	VecZeroEntries(gproX);
	VecZeroEntries(gproY);
	VecZeroEntries(gproZ);
	VecZeroEntries(lproX);
	VecZeroEntries(lproY);
	VecZeroEntries(lproZ);

	// If we want only a few indices we need to interpolate
	if(IOparam->Ap == 1)
	{
		// Create everything
	 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vx); CHKERRQ(ierr);
	 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vy); CHKERRQ(ierr);
	 	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam->mdI, PETSC_DETERMINE, &aop->vz); CHKERRQ(ierr);

		VecGetArray(aop->vx,&vx);
		VecGetArray(aop->vy,&vy);
		VecGetArray(aop->vz,&vz);

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
				
				if(IOparam->Av[ii] == 1)
				{
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					xe = ((coord_local[0] - ncx[I ])/(ncx[I +1] - ncx[I ]))/2; xb = (1.0 - xe)/2;
	
					llproX[sz+KK  ][sy+JJ  ][sx+I  ] = xb;
					llproX[sz+KK  ][sy+JJ  ][sx+I+1] = xe;
					llproX[sz+KK  ][sy+JJ+1][sx+I  ] = xb;
					llproX[sz+KK  ][sy+JJ+1][sx+I+1] = xe;
					llproX[sz+KK+1][sy+JJ  ][sx+I  ] = xb;
					llproX[sz+KK+1][sy+JJ  ][sx+I+1] = xe;
					llproX[sz+KK+1][sy+JJ+1][sx+I  ] = xb;
					llproX[sz+KK+1][sy+JJ+1][sx+I+1] = xe;
				}
				else if(IOparam->Av[ii] == 2)
				{
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					ye = ((coord_local[1] - ncy[J ])/(ncy[J +1] - ncy[J ]))/2; yb = (1.0 - 2*ye)/2;
	
					llproY[sz+KK  ][sy+J  ][sx+II  ] = yb;
					llproY[sz+KK  ][sy+J  ][sx+II+1] = yb;
					llproY[sz+KK  ][sy+J+1][sx+II  ] = ye;
					llproY[sz+KK  ][sy+J+1][sx+II+1] = ye;
					llproY[sz+KK+1][sy+J  ][sx+II  ] = yb;
					llproY[sz+KK+1][sy+J  ][sx+II+1] = yb;
					llproY[sz+KK+1][sy+J+1][sx+II  ] = ye;
					llproY[sz+KK+1][sy+J+1][sx+II+1] = ye;
				}
				else if(IOparam->Av[ii] == 3)
				{
					vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);
					vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);
					vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

					// get relative coordinates
					ze = ((coord_local[2] - ncz[K ])/(ncz[K +1] - ncz[K ]))/2; zb = (1.0 - 2*ze)/2;
	
					llproZ[sz+K  ][sy+JJ  ][sx+II  ] = zb;
					llproZ[sz+K  ][sy+JJ  ][sx+II+1] = zb;
					llproZ[sz+K  ][sy+JJ+1][sx+II  ] = zb;
					llproZ[sz+K  ][sy+JJ+1][sx+II+1] = zb;
					llproZ[sz+K+1][sy+JJ  ][sx+II  ] = ze;
					llproZ[sz+K+1][sy+JJ  ][sx+II+1] = ze;
					llproZ[sz+K+1][sy+JJ+1][sx+II  ] = ze;
					llproZ[sz+K+1][sy+JJ+1][sx+II+1] = ze;
				}

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
			if(IOparam->Av[ii] == 1)
			{
				ierr = VecSet(lproX,1);
			}
			else if(IOparam->Av[ii] == 2)
			{
				ierr = VecSet(lproY,1);
			}
			else if(IOparam->Av[ii] == 3)
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
						K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, z);
			
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
						K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, z);
			
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
						K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, z);
			
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
#define __FUNCT__ "AdjointFormResidual"
PetscErrorCode AdjointFormResidual(SNES snes, Vec x, Vec f, void *ctx, PetscInt CurPar, PetscInt CurPhase )
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear unused parameters
	if(snes) snes = NULL;

	// access context
	nl = (NLSol*)ctx;
	jr = nl->pc->pm->jr;

	/*// apply pressure limit at the first visco-plastic timestep and iteration
    if(jr->ts->istep == 1 && jr->ctrl->pLimPlast == PETSC_TRUE)
    {
    	jr->matLim.presLimFlg = PETSC_TRUE;
	}*/

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x); CHKERRQ(ierr);

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute residual
	if (CurPar==_RHO0_)
	{
		ierr = AdjointJacResGetResidual_ViscPowerlaw(jr, CurPar, CurPhase); CHKERRQ(ierr);
	}
	else if(CurPar==_ETA0_)
	{
		ierr = AdjointJacResGetResidual_ViscPowerlaw(jr, CurPar, CurPhase); CHKERRQ(ierr);
	}
	else if(CurPar==_N_)
	{
		ierr = AdjointJacResGetResidual_ViscPowerlaw(jr, CurPar, CurPhase); CHKERRQ(ierr);
	}
	else if(CurPar==_EN_)
	{
		ierr = AdjointJacResGetResidual_ViscPowerlaw(jr, CurPar, CurPhase); CHKERRQ(ierr);
	}
	else if(CurPar==_MFR_)
	{
		ierr = AdjointJacResGetResidual_ViscPowerlaw(jr, CurPar, CurPhase); CHKERRQ(ierr);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"ADJOINT ERROR: This residual computation is not known\n");
	}

	// copy residuals to global vector
	ierr = JacResCopyRes(jr, f); CHKERRQ(ierr);

	// deactivate pressure limit after it has been activated
	// jr->matLim.presLimFlg = PETSC_FALSE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointGradientPerturbParameter"
PetscErrorCode AdjointGradientPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop, Scaling *scal)
{
	PetscScalar         ini, perturb, curscal;
	Material_t         *mat;

	PetscFunctionBegin;

	mat = nl->pc->pm->jr->dbm->phases;

	// Get the perturbation value & scaling
	perturb = aop->Perturb;
	curscal = aop->CurScal;

	// Perturb the current parameter in the current phase (more to be included)
	if(CurPar==_RHO0_)			// rho
	{
		ini = mat[CurPhase].rho;
		perturb = ini*perturb;
		mat[CurPhase].rho +=  perturb;
		curscal = (scal->velocity)/(scal->density);
	}
	else if (CurPar==_RHON_)	    // rho_n
	{
		ini = mat[CurPhase].rho_n;
		perturb = ini*perturb;
		mat[CurPhase].rho_n +=  perturb;
		curscal = (scal->velocity)/1;
	}
	else if (CurPar==_RHOC_)	    // rho_c
	{
		ini = mat[CurPhase].rho_c;
		perturb = ini*perturb;
		mat[CurPhase].rho_c +=  perturb;
		curscal = (scal->velocity)*(scal->length_si);
	}
	else if (CurPar==_K_)	    // K
	{
		ini = mat[CurPhase].K;
		perturb = ini*perturb;
		mat[CurPhase].K +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_KP_)	    // Kp
	{
		ini = mat[CurPhase].Kp;
		perturb = ini*perturb;
		mat[CurPhase].Kp +=  perturb;
		curscal = (scal->velocity)/1;
	}
	else if (CurPar==_SHEAR_)	    // G
	{
		ini = mat[CurPhase].G;
		perturb = ini*perturb;
		mat[CurPhase].G +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_ETA_)	    // Bd
	{
		// This kind of perturbs the whole NEWTONIAN viscosity, consider perturbing the parameters directly
		ini = mat[CurPhase].Bd;
		PetscScalar BdTemp;
		perturb = perturb*(1.0/(2*ini));
		BdTemp = (1.0/(2*ini)) + perturb;//(perturb*(1.0/(2*ini)));
		mat[CurPhase].Bd =  (1.0/(2*BdTemp));
		curscal = (scal->velocity)/(scal->viscosity);
	}
	else if (CurPar==_ED_)	    // Ed
	{
		ini = mat[CurPhase].Ed;
		perturb = ini*perturb;
		mat[CurPhase].Ed +=  perturb;
		curscal = (scal->velocity)/(1);   // Not sure
	}
	else if (CurPar==_VD_)	// Vd
	{
		ini = mat[CurPhase].Vd;
		perturb = ini*perturb;
		mat[CurPhase].Vd +=  perturb;
		curscal = (scal->velocity)*(scal->stress_si);
	}
	else if (CurPar==_ETA0_)	// Bn
	{
		// This kind of perturbs the whole DISLOCATION viscosity, consider perturbing the parameters directly
		ini = mat[CurPhase].Bn;
		PetscScalar BnTemp;
		// -- Uncomment to compute gradient for ETA0 --
		perturb = perturb* (pow(mat[CurPhase].Bn * pow(2,mat[CurPhase].n) * pow(nl->pc->pm->jr->ctrl.DII_ref, mat[CurPhase].n-1) , -1/mat[CurPhase].n));
		BnTemp = (pow(mat[CurPhase].Bn * pow(2,mat[CurPhase].n) * pow(nl->pc->pm->jr->ctrl.DII_ref, mat[CurPhase].n-1) , -1/mat[CurPhase].n))  + perturb;
		mat[CurPhase].Bn = pow (2.0*BnTemp, -mat[CurPhase].n) * pow(nl->pc->pm->jr->ctrl.DII_ref, 1 - mat[CurPhase].n);
		// -- Uncomment to compute gradient for BN --
		//perturb = ini*perturb + 1e-12;
		//mat[CurPhase].Bn +=  perturb;
		curscal = (scal->velocity)/(scal->viscosity);
	}
	else if (CurPar== _N_)	// n
	{
		ini = mat[CurPhase].n;
		aop->Ini2 = mat[CurPhase].Bn;
		PetscScalar ViscTemp = (pow((aop->Ini2 * pow(2,mat[CurPhase].n) * pow(nl->pc->pm->jr->ctrl.DII_ref, mat[CurPhase].n-1)) , -1/mat[CurPhase].n));
		perturb = ini*perturb;
		mat[CurPhase].n +=  perturb;
		// We also accordingly need to perturb the inverse viscosity in this case
		mat[CurPhase].Bn = pow (2.0*ViscTemp, -mat[CurPhase].n) * pow(nl->pc->pm->jr->ctrl.DII_ref, 1 - mat[CurPhase].n);
		curscal = (scal->velocity)/(1);
	}
	else if (CurPar==_EN_)	// En
	{
		ini = mat[CurPhase].En;
		perturb = ini*perturb;
		mat[CurPhase].En +=  perturb;
		curscal = (scal->velocity)/(1);    // Not sure
	}
	else if (CurPar==_VN_)	// Vn
	{
		ini = mat[CurPhase].Vn;
		perturb = ini*perturb;
		mat[CurPhase].Vn +=  perturb;
		curscal = (scal->velocity)*(scal->stress_si);
	}
	else if (CurPar==_TAUP_)	// taup
	{
		ini = mat[CurPhase].taup;
		perturb = ini*perturb;
		mat[CurPhase].taup +=  perturb;
		curscal = (scal->velocity)/(scal->stress_si);
	}
	else if (CurPar==_GAMMA_)	// gamma
	{
		ini = mat[CurPhase].gamma;
		perturb = ini*perturb;
		mat[CurPhase].gamma +=  perturb;
		curscal = (scal->velocity)/(1);
	}
	else if (CurPar==_Q_)	// q
	{
		ini = mat[CurPhase].q;
		perturb = ini*perturb;
		mat[CurPhase].q +=  perturb;
		curscal = (scal->velocity)/(1);
	}
	else if (CurPar==_FRICTION_)	// fr
	{
		ini = mat[CurPhase].fr;
		perturb = ini*perturb;
		mat[CurPhase].fr +=  perturb;
		curscal = (scal->velocity)/scal->angle;
	}
	else if (CurPar==_COHESION_)	// ch
	{
		ini = mat[CurPhase].ch;
		perturb = ini*perturb;
		mat[CurPhase].ch +=  perturb;
		curscal = (scal->velocity)/scal->stress_si;
	}
	else if (CurPar==_CP_)	// Cp
	{
		ini = mat[CurPhase].Cp;
		perturb = ini*perturb;
		mat[CurPhase].Cp +=  perturb;
		curscal = (scal->velocity)/scal->cpecific_heat;
	}
	else if (CurPar==_A_)	// A
	{
		ini = mat[CurPhase].A;
		perturb = ini*perturb;
		mat[CurPhase].A +=  perturb;
		curscal = (scal->velocity)/scal->heat_production;
	}
	else if (CurPar==_MFR_)	// Melt fraction
	{
		ini = mat[CurPhase].mf;
		perturb = ini*perturb;
		mat[CurPhase].mf +=  perturb;
		curscal = (scal->velocity)/1;
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"ADJOINT ERROR: Definition of the current parameter is not defined ; Choose between [0-15]\n");
		PetscFunctionReturn(1);
	}

	// Store initial value of current parameter & scaling
	aop->Ini 		= ini;
	aop->CurScal 	= curscal;
	aop->Perturb    = perturb;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
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
	}else if (CurPar==_K_)  		{phases[CurPhase].K      = ini;
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
	else if (CurPar==_MFR_)			{phases[CurPhase].mf     = ini;}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointJacResGetResidual_ViscPowerlaw"
PetscErrorCode AdjointJacResGetResidual_ViscPowerlaw(JacRes *jr, PetscInt CurPar, PetscInt CurPhase)
{
	// Compute residual of nonlinear momentum and mass conservation
	// equations, based on pre-computed components of effective
	// strain-rate tensor, current values of pressure and temperature.
	// Missing components of the second invariant of the effective strain-rate
	// tensor (squares of the corresponding strain rate components) are averaged
	// form the hosting nodes using arithmetic mean.
	// DII = (0.5*D_ij*D_ij)^0.5
	// NOTE: we interpolate and average D_ij*D_ij terms instead of D_ij

	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	Soft_t     *soft;
	Controls   *ctrl;
	PetscInt    iter, numPhases, AirPhase;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, mcx, mcy, mcz;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar gx, gy, gz, tx, ty, tz, sxx, syy, szz, sxy, sxz, syz;
	PetscScalar J2Inv, theta, rho, IKdt, Tc, pc, pShift, pn, dt, fssa, xc, yc, *grav;
	PetscScalar ph, n, eII, Bn, Vn, RT, pp, ep, ef, En, e0, eta0, etamax, mf, Me_mu;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc, ***bcp;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***T, ***p_lith, ***p_pore;
	PetscScalar eta_creep, eta_vp;
	PetscScalar depth, pc_lith, pc_pore, biot, ptotal, avg_topo;
//	PetscScalar alpha, Tn,

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;
	bc = jr->bc;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	mx  = fs->dsx.tnods - 1;
	my  = fs->dsy.tnods - 1;
	mz  = fs->dsz.tnods - 1;

	// access residual context variables
	numPhases =  jr->dbm->numPhases; // number phases
	phases    =  jr->dbm->phases;    // phase parameters
	soft      =  jr->dbm->matSoft;   // material softening laws
	ctrl      = &jr->ctrl;           // control parameters
	dt        =  jr->ts->dt;         // time step
	fssa      =  ctrl->FSSA;         // density gradient penalty parameter
	grav      =  ctrl->grav;         // gravity acceleration
	pShift    =  ctrl->pShift;       // pressure shift
	biot      =  ctrl->biot;         // Biot pressure parameter
	AirPhase  =  jr->surf->AirPhase; // sticky air phase number
	avg_topo  =  jr->surf->avg_topo; // average surface topography

	// clear local residual vectors
	ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->gc);  CHKERRQ(ierr);

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
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

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

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

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

		// compute depth below the free surface
		if(AirPhase != -1) depth = avg_topo - COORD_CELL(k, sz, fs->dsz);
		else               depth = 0.0;
		if(depth < 0.0)    depth = 0.0;

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_vp, numPhases, phases, soft, svCell->phRat, ctrl, pc_lith, pc_pore, dt, pc-pShift, Tc, jr->Pd); CHKERRQ(ierr);

		// store creep viscosity
		svCell->eta_creep = eta_creep;
		svCell->eta_vp    = eta_vp;

		// compute stress, plastic strain rate and shear heating term on cell
		ierr = GetStressCell(svCell, XX, YY, ZZ); CHKERRQ(ierr);

		// get total pressure (effective pressure, computed by LaMEM, plus pore pressure)
		ptotal = pc + biot*pc_pore;

		// compute total Cauchy stresses
		sxx = svCell->sxx - ptotal;
		syy = svCell->syy - ptotal;
		szz = svCell->szz - ptotal;

		// evaluate volumetric constitutive equations
		ierr = VolConstEq(svBulk, numPhases, phases, svCell->phRat, ctrl, depth, dt, pc-pShift , Tc, jr->Pd); CHKERRQ(ierr);

		// access
		theta = svBulk->theta; // volumetric strain rate
		rho   = svBulk->rho;   // effective density
		IKdt  = svBulk->IKdt;  // inverse bulk viscosity
//		alpha = svBulk->alpha; // effective thermal expansion
		pn    = svBulk->pn;    // pressure history
//		Tn    = svBulk->Tn;    // temperature history

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
		
		// Get all necessary parameters for residual computation
		ph 		 = svCell->phRat[CurPhase];
		n 		 = phases[CurPhase].n;
		Bn 		 = phases[CurPhase].Bn;
		eII 	 =  svDev->DII;				
		if  (eII == 0.0) eII = ctrl->DII_ref;
		Vn 		 =  phases[CurPhase].Vn;
		En 		 =  phases[CurPhase].En;
		if(T) RT =  ctrl->Rugc*Tc;
		else  RT =  -1.0;
		pp 		 =  pc-pShift;
		e0 		 =  ctrl->DII_ref;
		eta0 	 =  pow(Bn*pow(2,n)*pow(e0,n-1),-1/n);
		etamax   =  1/ctrl->inv_eta_max;
		mf       =  phases[CurPhase].mf;
		Me_mu    =  phases[CurPhase].Me_Mu0;
		
		if (CurPar==_RHO0_)  // Only influences the center node
		{
			if(svCell->phRat[CurPhase] > 0)
			{
				fx[k][j][i] -= ph*(grav[0]*(0.5-(vx[k][j][i] * fssa*dt)/bdx));   fx[k][j][i+1] += ph*(grav[0]*(-(vx[k][j][i+1] * fssa*dt)/fdx - 0.5));
				fy[k][j][i] -= ph*(grav[1]*(0.5-(vy[k][j][i] * fssa*dt)/bdy));   fy[k][j+1][i] += ph*(grav[1]*(-(vy[k][j+1][i] * fssa*dt)/fdy - 0.5));
				fz[k][j][i] -= ph*(grav[2]*(0.5-(vz[k][j][i] * fssa*dt)/bdz));   fz[k+1][j][i] += ph*(grav[2]*(-(vz[k+1][j][i] * fssa*dt)/fdz - 0.5));
			}
		}
		else if(CurPar==_ETA0_)
		{
			// BN
			/*if(svCell->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				fx[k][j][i] -= ph*( (-XX*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdx*n) );   fx[k][j][i+1] += ph*( (-XX*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdx*n) );
				fy[k][j][i] -= ph*( (-YY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdy*n) );   fy[k][j+1][i] += ph*( (-YY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdy*n) );
				fz[k][j][i] -= ph*( (-ZZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdz*n) );   fz[k+1][j][i] += ph*( (-ZZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdz*n) );
			}*/
			// ETA0
			if(svCell->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k][j][i] -= ph*( (XX*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i+1] += ph*( (XX*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fy[k][j][i] -= ph*( (YY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j+1][i] += ph*( (YY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j][i] -= ph*( (ZZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k+1][j][i] += ph*( (ZZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
		}
		else if (CurPar==_N_)
		{
			if(svCell->phRat[CurPhase] > 0)
			{
				ep = exp((-En+pp*Vn)/RT);
				ef = pow(2,-n)*ep*pow(eta0,-n)*pow(e0,1-n);
				fx[k][j][i] -= ph*( (4*pow(eII,1+1/n)*XX*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fx[k][j][i+1] += ph*( (4*pow(eII,1+1/n)*XX*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
				fy[k][j][i] -= ph*( (4*pow(eII,1+1/n)*YY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fy[k][j+1][i] += ph*( (4*pow(eII,1+1/n)*YY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
				fz[k][j][i] -= ph*( (4*pow(eII,1+1/n)*ZZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fz[k+1][j][i] += ph*( (4*pow(eII,1+1/n)*ZZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
			}
		}
		else if(CurPar==_EN_)
		{
			if(svCell->phRat[CurPhase] > 0)
			{
				ep = exp((-pp*Vn-En)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k][j][i] -= ph*( (XX*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i+1] += ph*( (XX*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fy[k][j][i] -= ph*( (YY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j+1][i] += ph*( (YY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j][i] -= ph*( (ZZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k+1][j][i] += ph*( (ZZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
		}
		else if(CurPar==_MFR_)
		{
			if(svCell->phRat[CurPhase] > 0)
			{
				ep = ((1-mf)/mf);
				ef = exp(2.5+(1-mf)*pow(ep,0.48));
				fx[k][j][i] -= ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XX - vx[k][j][i]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/bdx) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));   fx[k][j][i+1] += ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XX - vx[k][j][i+1]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/fdx) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));
				fy[k][j][i] -= ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YY - vy[k][j][i]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/bdy) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));   fy[k][j+1][i] += ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YY - vy[k][j+1][i]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/fdy) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));
				fz[k][j][i] -= ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*ZZ - vz[k][j][i]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/bdz) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));   fz[k+1][j][i] += ph*(((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*ZZ - vz[k+1][j][i]*grav[0]*dt*fssa*(svBulk->rho_pf-svBulk->rho_pd))/fdz) + ((grav[0]*(svBulk->rho_pf-svBulk->rho_pd))/2));
			}
		}
		else
		{
			// momentum
			fx[k][j][i] -= 0;   fx[k][j][i+1] += 0;
			fy[k][j][i] -= 0;   fy[k][j+1][i] += 0;
			fz[k][j][i] -= 0;   fz[k+1][j][i] += 0;
		}


//****************************************
// ADHOC (HARD-CODED PRESSURE CONSTRAINTS)
//****************************************

//		if(k == 0)   fz[k][j][i]   += -p[k-1][j][i]/bdz;
//		if(k == mcz) fz[k+1][j][i] -= -p[k+1][j][i]/fdz;

		// mass - currently T-dependency is deactivated
//		gc[k][j][i] = -IKdt*(pc - pn) - theta + alpha*(Tc - Tn)/dt;
        
        gc[k][j][i] = -IKdt*(pc - pn) - theta ;
        
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
		svDev  = &svEdge->svDev;

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

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

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

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_vp, numPhases, phases, soft, svEdge->phRat, ctrl, pc_lith, pc_pore, dt, pc-pShift, Tc, jr->Pd); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, XY); CHKERRQ(ierr);

		// access xy component of the Cauchy stress
		sxy = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		
		// Get all necessary parameters for residual computation
		ph 		 = svEdge->phRat[CurPhase];
		n 		 = phases[CurPhase].n;
		Bn 		 = phases[CurPhase].Bn;
		eII 	 =  svDev->DII;				
		if  (eII == 0.0) eII = ctrl->DII_ref;
		Vn 		 =  phases[CurPhase].Vn;
		En 		 =  phases[CurPhase].En;
		if(T) RT =  ctrl->Rugc*Tc;
		else  RT =  -1.0;
		pp 		 =  pc-pShift;
		e0 		 =  ctrl->DII_ref;
		eta0 	 =  pow(Bn*pow(2,n)*pow(e0,n-1),-1/n);
		etamax   =  1/ctrl->inv_eta_max;
		mf       =  phases[CurPhase].mf;
		Me_mu    =  phases[CurPhase].Me_Mu0;

		if (CurPar==_RHO0_)
		{
			// Density only defined at the center
		}
		else if (CurPar==_ETA0_)
		{
			// ETA0
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k][j-1][i] -= ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i] += ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fy[k][j][i-1] -= ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j][i] += ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
			// BN
			/*if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				fx[k][j-1][i] -= ph*( (-XY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdy*n) );   fx[k][j][i] += ph*( (-XY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdy*n) );
				fy[k][j][i-1] -= ph*( (-XY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdx*n) );   fy[k][j][i] += ph*( (-XY*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdx*n) );
			}*/
		}
		else if(CurPar==_N_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En+pp*Vn)/RT);
				ef = pow(2,-n)*ep*pow(eta0,-n)*pow(e0,1-n);
				fx[k][j-1][i] -= ph*( (4*pow(eII,1+1/n)*XY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fx[k][j][i] += ph*( (4*pow(eII,1+1/n)*XY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
				fy[k][j][i-1] -= ph*( (4*pow(eII,1+1/n)*XY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fy[k][j][i] += ph*( (4*pow(eII,1+1/n)*XY*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
			}
		}
		else if(CurPar==_EN_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-pp*Vn-En)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k][j-1][i] -= ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i] += ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fy[k][j][i-1] -= ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j][i] += ph*( (XY*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
		}
		else if(CurPar==_MFR_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = ((1-mf)/mf);
				ef = exp(2.5+(1-mf)*pow(ep,0.48));
				fx[k][j-1][i] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XY) /bdy);   fx[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XY) /fdy);
				fy[k][j][i-1] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XY) /bdx);   fy[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XY) /fdx);
			}
		}
		else
		{
			fx[k][j-1][i] -= 0;   fx[k][j][i] += 0;
			fy[k][j][i-1] -= 0;   fy[k][j][i] += 0;
		}

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
		svDev  = &svEdge->svDev;

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

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

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

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_vp, numPhases, phases, soft, svEdge->phRat, ctrl, pc_lith, pc_pore, dt, pc-pShift, Tc, jr->Pd); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, XZ); CHKERRQ(ierr);

		// access xz component of the Cauchy stress
		sxz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);
		
		// Get all necessary parameters for residual computation
		ph 		 = svEdge->phRat[CurPhase];
		n 		 = phases[CurPhase].n;
		Bn 		 = phases[CurPhase].Bn;
		eII 	 =  svDev->DII;				
		if  (eII == 0.0) eII = ctrl->DII_ref;
		Vn 		 =  phases[CurPhase].Vn;
		En 		 =  phases[CurPhase].En;
		if(T) RT =  ctrl->Rugc*Tc;
		else  RT =  -1.0;
		pp 		 =  pc-pShift;
		e0 		 =  ctrl->DII_ref;
		eta0 	 =  pow(Bn*pow(2,n)*pow(e0,n-1),-1/n);
		etamax   =  1/ctrl->inv_eta_max;
		mf       =  phases[CurPhase].mf;
		Me_mu    =  phases[CurPhase].Me_Mu0;

		// momentum
		if (CurPar==_RHO0_)
		{
			// Density only defined at the center
		}
		else if (CurPar==_ETA0_)
		{
			// ETA0
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k-1][j][i] -= ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i] += ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j][i-1] -= ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k][j][i] += ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdx*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
			// BN
			/*if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				fx[k-1][j][i] -= ph*( (-XZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdz*n) );   fx[k][j][i] += ph*( (-XZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdz*n) );
				fz[k][j][i-1] -= ph*( (-XZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdx*n) );   fz[k][j][i] += ph*( (-XZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdx*n) );
			}*/
		}
		else if(CurPar==_N_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En+pp*Vn)/RT);
				ef = pow(2,-n)*ep*pow(eta0,-n)*pow(e0,1-n);
				fx[k-1][j][i] -= ph*( (4*pow(eII,1+1/n)*XZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fx[k][j][i] += ph*( (4*pow(eII,1+1/n)*XZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
				fz[k][j][i-1] -= ph*( (4*pow(eII,1+1/n)*XZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fz[k][j][i] += ph*( (4*pow(eII,1+1/n)*XZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdx*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
			}
		}
		else if(CurPar==_EN_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-pp*Vn-En)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fx[k-1][j][i] -= ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fx[k][j][i] += ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j][i-1] -= ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k][j][i] += ph*( (XZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdx*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
		}
		else if(CurPar==_MFR_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = ((1-mf)/mf);
				ef = exp(2.5+(1-mf)*pow(ep,0.48));
				fx[k-1][j][i] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XZ) /bdz);   fx[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XZ) /fdz);
				fz[k][j][i-1] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XZ) /bdx);   fz[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*XZ) /fdx);
			}
		}
		else
		{
			fx[k-1][j][i] -= 0;   fx[k][j][i] += 0;
			fz[k][j][i-1] -= 0;   fz[k][j][i] += 0;
		}

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
		svDev  = &svEdge->svDev;

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

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

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

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_vp, numPhases, phases, soft, svEdge->phRat, ctrl, pc_lith, pc_pore, dt, pc-pShift, Tc, jr->Pd); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, YZ); CHKERRQ(ierr);

		// access yz component of the Cauchy stress
		syz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);
		
		// Get all necessary parameters for residual computation
		ph 		 =  svEdge->phRat[CurPhase];
		n 		 =  phases[CurPhase].n;
		Bn 		 =  phases[CurPhase].Bn;
		eII 	 =  svDev->DII;				
		if  (eII == 0.0) eII = ctrl->DII_ref;
		Vn 		 =  phases[CurPhase].Vn;
		En 		 =  phases[CurPhase].En;
		if(T) RT =  ctrl->Rugc*Tc;
		else  RT =  -1.0;
		pp 		 =  pc-pShift;
		e0 		 =  ctrl->DII_ref;
		eta0 	 =  pow(Bn*pow(2,n)*pow(e0,n-1),-1/n);
		etamax   =  1/ctrl->inv_eta_max;
		mf       =  phases[CurPhase].mf;
		Me_mu    =  phases[CurPhase].Me_Mu0;

		// update momentum residuals
		if (CurPar==_RHO0_)
		{
			// Density only defined at the center
		}
		else if (CurPar==_ETA0_)
		{
			// ETA0
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fy[k-1][j][i] -= ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j][i] += ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdz*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j-1][i] -= ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(bdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k][j][i] += ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n-1)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(fdy*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
			// BN
			/*if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En-pp*Vn)/RT);
				fy[k-1][j][i] -= ph*( (-YZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdz*n) );   fy[k][j][i] += ph*( (-YZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdz*n) );
				fz[k][j-1][i] -= ph*( (-YZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(bdy*n) );   fz[k][j][i] += ph*( (-YZ*pow(eII,1/n-1)*ep*pow(Bn*ep,-1/n-1))/(fdy*n) );
			}*/
		}
		else if(CurPar==_N_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-En+pp*Vn)/RT);
				ef = pow(2,-n)*ep*pow(eta0,-n)*pow(e0,1-n);
				fy[k-1][j][i] -= ph*( (4*pow(eII,1+1/n)*YZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fy[k][j][i] += ph*( (4*pow(eII,1+1/n)*YZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdz*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
				fz[k][j-1][i] -= ph*( (4*pow(eII,1+1/n)*YZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(bdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );   fz[k][j][i] += ph*( (4*pow(eII,1+1/n)*YZ*pow(etamax,2)*pow(ef,1/n)*(-log(eII)+n*(log(eta0)+log(2*e0))+log(ef)))/(fdy*pow(n,2)*pow(pow(eII,1/n)+2*eII*etamax*pow(ef,1/n),2)) );
			}
		}
		else if(CurPar==_EN_)
		{
			if(svEdge->phRat[CurPhase] > 0)
			{
				ep = exp((-pp*Vn-En)/RT);
				ef = pow(2,-n)*pow(eta0,-n)*pow(e0,1-n)*ep;
				fy[k-1][j][i] -= ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fy[k][j][i] += ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdz*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
				fz[k][j-1][i] -= ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*bdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );   fz[k][j][i] += ph*( (YZ*pow(2,2-n)*pow(eII,1-1/n)*pow(eta0,-n)*pow(e0,1-n)*ep*pow(ef,1/n-1))/(RT*fdy*n*pow(2*pow(eII,1-1/n)*pow(ef,1/n)+1/etamax,2)) );
			}
		}
		else if(CurPar==_MFR_)
		{
			ep = ((1-mf)/mf);
			ef = exp(2.5+(1-mf)*pow(ep,0.48));
			if(svEdge->phRat[CurPhase] > 0)
			{
				fy[k-1][j][i] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YZ) /bdz);   fy[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YZ) /fdz);
				fz[k][j-1][i] -= ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YZ) /bdy);   fz[k][j][i] += ph*((2*( Me_mu*ef*((0.48*(-(1-mf)/pow(mf,2)-(1/mf))*(1-mf))/(pow(ep,0.52)) - pow(ep,0.48)) )*YZ) /fdy);
			}
		}
		else
		{
			fy[k-1][j][i] -= 0;   fy[k][j][i] += 0;
			fz[k][j-1][i] -= 0;   fz[k][j][i] += 0;
		}
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

	PetscFunctionReturn(0);
}



