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
// COMPUTATION OF ADJOINT GRADIENTS
//---------------------------------------------------------------------------
// RECIPE:
// Objective function    F(x,x(p)) = P * x^T                // p = parameter ; x = converged solution ; Projection vector containing the proportions of solution influence
// Derivative I          dF/dx     = P
// Adjoint operation     psi       = J^-T * dF/dx           // J = converged Jacobain matrix
// Derivative II         dr/dp     = [r(p+h) - r(p)]/h      // finite difference approximation of derivative of residual r vs parameter
// Gradients             dF/dp     = -psi^T * dr/dp
//
// USAGE:
// In your .dat file you need to define:
// ComputeAdjointGradients	       = 1	                    // 1 = Compute gradients 0 = NO
// Adjoint_x					   = 1.2                    // Array (same length as Adjoint_y, Adjoint_z) containing the x coordinates of the point where you want to compute the gradients
// Adjoint_y					   = 0.6                    // Array (same length as Adjoint_x, Adjoint_z) containing the y coordinates of the point where you want to compute the gradients
// Adjoint_z					   = 0.4                    // Array (same length as Adjoint_x, Adjoint_y) containing the z coordinates of the point where you want to compute the gradients
// AdjointPhases                   = 1 2 1 2                // Array (same length as AdjointParameters) containing the phase of the parameter
// AdjointParameters               = 2 1 1 1                // Array (same length as AdjointPhases) containing the parameter corresponding to the phase
// AdjointVel                      = 3                      // Array (same length as Adjoint_x, Adjoint_y, Adjoint_z) containing the related velocity direction in which to compute the gradient
//
//                              Density         Elasticity  Diff creep    Dis creep    Peierl creep
// Possible parameters: (rho rho_n rho_c beta)   (K Kp G)   (Bd Ed Vd)  (Bn n  En Vn) (taup gamma q )
//                      ( 1     2     3    4 )   (5  6 7)   ( 8  9 10)  (11 12 13 14) ( 15    16  17)
//
// Possible Velocities:  (Vx   Vy   Vz
//                       (1    2    3)
//
// You can control the behaviour of the KSP object for the adjoint with the prefix "as_"
//
// IMPORTANT: Since the Adjoint needs the Jacobian matrix for computing the gradients it's crucial to make sure that
// you compute the Jacobian matrix in the timesteps where you want to compute the gradients
// (f.e. a linear problem would need a low value for -snes_atol [1e-20] and a low max iteration count -snes_max_it [2] to
// guarantee the computation of the Jacobian + the option '-snes_type ksponly'
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

	ierr = PetscMemzero(aop, sizeof(AdjGrad)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CreateAdjoint"
PetscErrorCode CreateAdjoint(JacRes *jr, UserCtx *user, AdjGrad *aop, NLSol *nl,SNES snes)
{

	/* TODO:
		- It's only valid vor visco-elasicity since the Jacobian is assumed to be symmetric (future: Jacobian needs to be transposed)
		- boundary of 1e-10 if Perturb parameter becomes too small
		- the only working comparison variable is the velocity
	*/

	PetscPrintf(PETSC_COMM_WORLD,"Computation of adjoint gradients\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Initialize all variables
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, j, CurPhase, CurPar;
	PetscScalar         grd, Perturb, vx[_MAX_AdjointIndices_], vy[_MAX_AdjointIndices_], vz[_MAX_AdjointIndices_];
	Vec 				rpl, sol, psi, pro, drdp, res, Perturb_vec;
	PC                  ipc;
	Scaling             *scal;

	// Set perturbation paramter for the finite differences back to 1e-6
	aop->Perturb = 1e-6;

	// Profile time
	PetscLogDouble     cputime_start, cputime_end;
	PetscTime(&cputime_start);

	scal = &jr->scal;

	// Create all needed vectors in the same size as the solution vector
	ierr = VecDuplicate(jr->gsol, &pro);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &psi);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &rpl);		 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &res);	 	 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sol); 		 CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &drdp);	 	 CHKERRQ(ierr);
	ierr = VecCopy(jr->gsol,sol); 				 CHKERRQ(ierr);
	ierr = VecCopy(jr->gres,res); 				 CHKERRQ(ierr);

	//===============
	// Indexing
	//===============
	// Put the proportion into the Projection vector where the user defined the computation coordinates (dF/dx = P)
	ierr = AdjointPointInPro(jr, user, vx, vy, vz, pro); 		CHKERRQ(ierr);

	//===============
	// SOLVE
	//===============
	// Solve the adjoint equation (psi = J^-T * dF/dx)
	ierr = SNESGetKSP(snes, &ksp);         		CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"as_"); 		CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         		CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            		CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          		CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,nl->J,nl->P);	CHKERRQ(ierr);
	ierr = KSPSolve(ksp,pro,psi);				CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);	CHKERRQ(ierr);

	//===============
	// Parameter loop
	//===============
	for(j = 0; j < user->AdjointNumPar; j++)
	{
		// Get residual since it is overwritten in VecAYPX
		ierr = VecDuplicate(jr->gres, &res);  	CHKERRQ(ierr);
		ierr = VecCopy(jr->gres,res); 			CHKERRQ(ierr);

		// Get current phase and parameter
		CurPhase = user->AdjointPhases[j];
		CurPar   = user->AdjointParameters[j];

		// Perturb the current parameter in the current phase
		ierr = AdjointPerturbParameter(nl, CurPar, CurPhase, aop);      CHKERRQ(ierr);

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
		aop->grad[j] 	= -1 * grd;

		// Reset perturbed parameter
		ierr = AdjointResetParameter(nl, CurPar, CurPhase, aop);       CHKERRQ(ierr);

		// Print result
		PetscPrintf(PETSC_COMM_WORLD,"%D.Gradient = %.12f ; CurPar = %d ; CurPhase = %d\n",j+1,aop->grad[j]*scal->velocity,CurPar,CurPhase);

		// Destroy overwritten residual vector
		ierr = VecDestroy(&res);
	}

	// Print the solution variable at the user defined index
	for (i=0; i<user->AdjointNumInd; i++)
	{
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Computation variable = %.12f\n",vz[i]*scal->velocity);
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	}

	// Clean
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&pro);
	ierr = VecDestroy(&drdp);
	ierr = VecDestroy(&rpl);
	ierr = VecDestroy(&Perturb_vec);

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful & took %g s\n------------------------------------------\n",cputime_end - cputime_start);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPointInPro"
PetscErrorCode AdjointPointInPro(JacRes *jr, UserCtx *user, PetscScalar *vx, PetscScalar *vy, PetscScalar *vz, Vec pro)
{
	PetscErrorCode      ierr;
	FDSTAG              *fs;
	Vec                 lproX, lproY, lproZ, gproX, gproY, gproZ;
	PetscScalar         coord_local[3], *temppro, ***llproX, ***llproY, ***llproZ, *dggproX, *dggproY, *dggproZ;
	PetscInt            ii, sx, sy, sz, nx, ny, nz, I, J, K, II, JJ, KK, lrank, grank;
	PetscScalar         xb, yb, zb, xe, ye, ze, xc, yc, zc, *iter, *ncx, *ncy, *ncz, *ccx, *ccy, *ccz, ***lvx, ***lvy, ***lvz;

	PetscFunctionBegin;

	fs = jr->fs;

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

	VecZeroEntries(gproX);
	VecZeroEntries(gproY);
	VecZeroEntries(gproZ);
	VecZeroEntries(lproX);
	VecZeroEntries(lproY);
	VecZeroEntries(lproZ);

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

			if(user->AdjointVel[ii] == 1)
			{   // Vx
				// interpolate x velocity
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
			else if(user->AdjointVel[ii] == 2)
			{   // Vy
				// interpolate y velocity
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
			else if(user->AdjointVel[ii] == 3)
			{   // Vz
				// interpolate z velocity
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

	// Put the proportion into the Projection vector where the user defined the computation coordinates (dF/dx = P)
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

	ierr = VecDestroy(&lproX);
	ierr = VecDestroy(&lproY);
	ierr = VecDestroy(&lproZ);
	ierr = VecDestroy(&gproX);
	ierr = VecDestroy(&gproY);
	ierr = VecDestroy(&gproZ);

	ierr = DMDAVecRestoreArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointPerturbParameter"
PetscErrorCode AdjointPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop)
{
	PetscScalar ini, perturb;

	PetscFunctionBegin;

	// Get the perturbation value
	perturb = aop->Perturb;

	// Perturb the parameter in the current phase (more to be included)
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
		BdTemp = (1.0/(ini)) + perturb;   // might actually be 1/(2*ini) or 1/(ini/2)
		nl->pc->pm->jr->phases[CurPhase].Bd =  (1.0/(BdTemp));
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
		BnTemp = (1.0/(ini)) + perturb;   // might actually be 1/2*ini or 1/(ini/2)
		nl->pc->pm->jr->phases[CurPhase].Bn =  (1.0/(BnTemp));
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
#define __FUNCT__ "AdjointResetParameter"
PetscErrorCode AdjointResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop)
{
	PetscScalar ini;

	PetscFunctionBegin;

	// Get initial value of currently perturbed parameter
	ini = aop->Ini;

	// Set the used parameter back to its original value
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
