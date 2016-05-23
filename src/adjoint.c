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
// Adjoint_x					   = 1.2                    // Array containing the x coordinates of the point where you want to compute the gradients
// Adjoint_y					   = 0.6                    // Array containing the y coordinates of the point where you want to compute the gradients
// Adjoint_z					   = 0.4                    // Array containing the z coordinates of the point where you want to compute the gradients
// AdjointPhases                   = 1 2 1 2                // Array (same length as AdjointParameters) containing the phase of the parameter
// AdjointParameters               = 2 1 1 1                // Array (same length as AdjointPhases) containing the parameter corresponding to the phase
//
//                              Density         Elasticity  Diff creep    Dis creep    Peierl creep
// Possible parameters: (rho rho_n rho_c beta)   (K Kp G)   (Bd Ed Vd)  (Bn n  En Vn) (taup gamma q )
//                      ( 1     2     3    4 )   (5  6 7)   ( 8  9 10)  (11 12 13 14) ( 15    16  17)
//
// You can control the behaviour of the KSP object for the adjoint with the prefix "as_"
// IMPORTANT: Since the Adjoint needs the Jacobian matrix for computing the gradients it's crucial to make sure that
// you compute the Jacobian matrix in the timesteps where you want to compute the gradients
// (f.e. a linear problem would need a low value for -snes_atol [1e-20] and a low max iteration count -snes_max_it [2] to guarantee the computation of the Jacobian
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
PetscErrorCode CreateAdjoint(JacRes *jr, UserCtx *user, AdjGrad *aop, NLSol *nl,SNES snes, AdvCtx *actx)
{

	/* TODO:
		- It's only valid vor visco-elasicity since the Jacobian is assumed to be symmetric
		- So far only working in sequentiel mode (maybe the summation of the vector is worng in parallel?)
		- boundary of 1e-10 if Perturb parameter becomes too small
		- thereotical capable of multiple points but so far only tested for ONE point!
	*/

	PetscPrintf(PETSC_COMM_WORLD,"Computation of adjoint gradients\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Initialize all variables
	// PCStokes            pc;
	FDSTAG              *fs;
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, ii, j, lrank, grank, sx, sy, sz, nx, ny, nz, I, J, K, II, JJ, KK;
	PetscScalar         grd, xb, yb, zb, xe, ye, ze, Perturb, coord_local[3], vx[100], vy[100], vz[100], xc, yc, zc, *iter, *ncx, *ncy, *ncz, *ccx, *ccy, *ccz, ***lvx, ***lvy, ***lvz, *temppro, *tempproX, *tempproY, *tempproZ;       // Temporary vectors for the computation
	PetscScalar         rhoIni, rho_nIni, rho_cIni, betaIni, KIni, KpIni, GIni, BdIni, EdIni, VdIni, BnIni, nIni, EnIni, VnIni, taupIni, gammaIni, qIni;   // Initialize parameters for which the gradients can be computed
	Vec 				rpl, sol, psi, pro, proX, proY, proZ, drdp, res, Perturb_vec;
	PC                  ipc;

	// access context
	fs = actx->fs;

	// Set perturbation paramter for the finite differences
	Perturb = 1e-6;

	// Create all needed vectors in the same size as the solution vector
	ierr = VecDuplicate(jr->gsol, &pro); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gfx, &proX); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gfy, &proY); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gfz, &proZ); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &psi); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &rpl); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &res); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sol); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &drdp); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &Perturb_vec); CHKERRQ(ierr);
	ierr = VecCopy(jr->gsol,sol); CHKERRQ(ierr);
	ierr = VecCopy(jr->gres,res); CHKERRQ(ierr);
	ierr = VecSet(Perturb_vec,Perturb);

	// scan markers
	for(ii = 0; ii < user->AdjointNumInd; ii++) {
		// Create coordinate vector
		coord_local[0] = user->Adjoint_x[ii];
		coord_local[1] = user->Adjoint_y[ii];
		coord_local[2] = user->Adjoint_z[ii];

		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, coord_local, &lrank, &grank); CHKERRQ(ierr);

		if(lrank == -1) {
		} else {
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

			// Access the local velocities
			ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

			// get coordinates of cell center
			xc = ccx[I];
			yc = ccy[J];
			zc = ccz[K];

			// map marker on the cells of X, Y, Z & center grids
			if(coord_local[0] > xc) { II = I; } else { II = I-1; }
			if(coord_local[1] > yc) { JJ = J; } else { JJ = J-1; }
			if(coord_local[2] > zc) { KK = K; } else { KK = K-1; }

			ierr = VecGetArray(proX, &tempproX);      CHKERRQ(ierr);
			ierr = VecGetArray(proY, &tempproY);      CHKERRQ(ierr);
			ierr = VecGetArray(proZ, &tempproZ);      CHKERRQ(ierr);

			if(user->AdjointVel[ii] == 1){                  // Vx
				// interpolate x velocity
				vx[ii] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ncx, ccy, ccz);

				// get relative coordinates
				xe = (coord_local[0] - ncx[I])/(ncx[I+1] - ncx[I]); xb = 1.0 - xe;
				ye = (coord_local[1] - ccy[J])/(ccy[J+1] - ccy[J]); yb = 1.0 - ye;
				ze = (coord_local[2] - ccz[K])/(ccz[K+1] - ccz[K]); zb = 1.0 - ze;

				tempproX[sx+I  ] = xb;
				tempproX[sx+I+1] = xe;

				tempproY[sy+J  ] = yb;
				tempproY[sy+J+1] = ye;

				tempproZ[sz+K  ] = zb;
				tempproZ[sz+K+1] = ze;

			}else if(user->AdjointVel[ii] == 2){                  // Vy
				// interpolate y velocity
				vy[ii] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ncy, ccz);

				// get relative coordinates
				xe = (coord_local[0] - ccx[I])/(ccx[I+1] - ccx[I]); xb = 1.0 - xe;
				ye = (coord_local[1] - ncy[J])/(ncy[J+1] - ncy[J]); yb = 1.0 - ye;
				ze = (coord_local[2] - ccz[K])/(ccz[K+1] - ccz[K]); zb = 1.0 - ze;

				tempproX[sx+I  ] = xb;
				tempproX[sx+I+1] = xe;

				tempproY[sy+J  ] = yb;
				tempproY[sy+J+1] = ye;

				tempproZ[sz+K  ] = zb;
				tempproZ[sz+K+1] = ze;

			}else if(user->AdjointVel[ii] == 3){                  // Vz
				// interpolate z velocity
				vz[ii] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, coord_local[0], coord_local[1], coord_local[2], ccx, ccy, ncz);

				// get relative coordinates
				xe = (coord_local[0] - ccx[I])/(ccx[I+1] - ccx[I]); xb = 1.0 - xe;
				ye = (coord_local[1] - ccy[J])/(ccy[J+1] - ccy[J]); yb = 1.0 - ye;
				ze = (coord_local[2] - ncz[K])/(ncz[K+1] - ncz[K]); zb = 1.0 - ze;

				tempproX[sx+I  ] = xb;
				tempproX[sx+I+1] = xe;

				tempproY[sy+J  ] = yb;
				tempproY[sy+J+1] = ye;

				tempproZ[sz+K  ] = zb;
				tempproZ[sz+K+1] = ze;
			}

			// Put the proportion into the Projection vector where the user defined the computation coordinates (dF/dx = P)
			ierr = VecGetArray(pro, &temppro);      CHKERRQ(ierr);
			iter = temppro;

			ierr  = PetscMemcpy(iter, tempproX, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
			iter += fs->nXFace;

			ierr  = PetscMemcpy(iter, tempproY, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
			iter += fs->nYFace;

			ierr  = PetscMemcpy(iter, tempproZ, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
			iter += fs->nZFace;

			// restore access
			ierr = VecRestoreArray(proX, &tempproX);      CHKERRQ(ierr);
			ierr = VecRestoreArray(proY, &tempproY);      CHKERRQ(ierr);
			ierr = VecRestoreArray(proZ, &tempproZ);      CHKERRQ(ierr);
			ierr = VecRestoreArray(pro, &temppro);       CHKERRQ(ierr);
		}
	}

	// Solve the adjoint equation (psi = J^-T * dF/dx)
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"as_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,nl->J,nl->P);	CHKERRQ(ierr);
	ierr = KSPSolve(ksp,pro,psi);	CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);	CHKERRQ(ierr);

	//===============
	// Parameter loop
	//===============
	for(j=0; j<user->AdjointNumPar; j++){

		// Get residual since it is overwritten in VecAYPX
		ierr = VecDuplicate(jr->gres, &res); CHKERRQ(ierr);
		ierr = VecDuplicate(jr->gsol, &drdp); CHKERRQ(ierr);
		ierr = VecDuplicate(jr->gres, &rpl); CHKERRQ(ierr);
		ierr = VecCopy(jr->gres,res); CHKERRQ(ierr);

		// Get current phase and parameter
		PetscInt CurPhase = user->AdjointPhases[j];
		PetscInt CurPar   = user->AdjointParameters[j];

		// Perturb the parameter in the current phase (more to be included)
		if(CurPar==1) {
			// PetscScalar         rhoIni;
			rhoIni = nl->pc->pm->jr->phases[CurPhase].rho;
			nl->pc->pm->jr->phases[CurPhase].rho +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].rho);
		}else if (CurPar==2) {
			// PetscScalar         rho_nIni;
			rho_nIni = nl->pc->pm->jr->phases[CurPhase].rho_n;
			nl->pc->pm->jr->phases[CurPhase].rho_n +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].rho_n);
		}else if (CurPar==3) {
			// PetscScalar         rho_cIni;
			rho_cIni = nl->pc->pm->jr->phases[CurPhase].rho_c;
			nl->pc->pm->jr->phases[CurPhase].rho_c +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].rho_c);
		}else if (CurPar==4) {
			// PetscScalar         betaIni;
			betaIni = nl->pc->pm->jr->phases[CurPhase].beta;
			nl->pc->pm->jr->phases[CurPhase].beta +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].beta);
		}else if (CurPar==5) {
			// PetscScalar         KIni;
			KIni = nl->pc->pm->jr->phases[CurPhase].K;
			nl->pc->pm->jr->phases[CurPhase].K +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].K);
		}else if (CurPar==6) {
			// PetscScalar         KpIni;
			KpIni = nl->pc->pm->jr->phases[CurPhase].Kp;
			nl->pc->pm->jr->phases[CurPhase].Kp +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Kp);
		}else if (CurPar==7) {
			// PetscScalar         GIni;
			GIni = nl->pc->pm->jr->phases[CurPhase].G;
			nl->pc->pm->jr->phases[CurPhase].G +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].G);
		}else if (CurPar==8) {
			// PetscScalar         BdIni;
			BdIni = nl->pc->pm->jr->phases[CurPhase].Bd;
			nl->pc->pm->jr->phases[CurPhase].Bd +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Bd);
		}else if (CurPar==9) {
			// PetscScalar         EdIni;
			EdIni = nl->pc->pm->jr->phases[CurPhase].Ed;
			nl->pc->pm->jr->phases[CurPhase].Ed +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Ed);
		}else if (CurPar==10) {
			// PetscScalar         VdIni;
			VdIni = nl->pc->pm->jr->phases[CurPhase].Vd;
			nl->pc->pm->jr->phases[CurPhase].Vd +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Vd);
		}else if (CurPar==11) {
			// PetscScalar         BnIni;
			BnIni = nl->pc->pm->jr->phases[CurPhase].Bn;
			nl->pc->pm->jr->phases[CurPhase].Bn +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Bn);
		}else if (CurPar==12) {
			// PetscScalar         nIni;
			nIni = nl->pc->pm->jr->phases[CurPhase].n;
			nl->pc->pm->jr->phases[CurPhase].n +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].n);
		}else if (CurPar==13) {
			// PetscScalar         EnIni;
			EnIni = nl->pc->pm->jr->phases[CurPhase].En;
			nl->pc->pm->jr->phases[CurPhase].En +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].En);
		}else if (CurPar==14) {
			// PetscScalar         VnIni;
			VnIni = nl->pc->pm->jr->phases[CurPhase].Vn;
			nl->pc->pm->jr->phases[CurPhase].Vn +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].Vn);
		}else if (CurPar==15) {
			// PetscScalar         taupIni;
			taupIni = nl->pc->pm->jr->phases[CurPhase].taup;
			nl->pc->pm->jr->phases[CurPhase].taup +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].taup);
		}else if (CurPar==16) {
			// PetscScalar         gammaIni;
			gammaIni = nl->pc->pm->jr->phases[CurPhase].gamma;
			nl->pc->pm->jr->phases[CurPhase].gamma +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].gamma);
		}else if (CurPar==17) {
			// PetscScalar         qIni;
			qIni = nl->pc->pm->jr->phases[CurPhase].q;
			nl->pc->pm->jr->phases[CurPhase].q +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].q);
		}/*else if (CurPar==18) {
			// PetscScalar         etaIni;
			etaIni = nl->pc->pm->jr->phases[CurPhase].eta;
			nl->pc->pm->jr->phases[CurPhase].eta +=  (Perturb*nl->pc->pm->jr->phases[CurPhase].eta);
		}*/

		// Compute residual with the converged Jacobian and preconditioner (dr/dp = [r(p+h) - r(p)]/h)
		ierr = FormResidual(snes, sol, rpl, nl); 	       CHKERRQ(ierr);
		ierr = VecAYPX(res,-1,rpl);                        CHKERRQ(ierr);
		ierr = VecPointwiseDivide(drdp,res,Perturb_vec);   CHKERRQ(ierr);

		// Compute the gradient (dF/dp = -psi^T * dr/dp)
		ierr = VecDot(psi, drdp, &grd);     CHKERRQ(ierr);

		// Save gradient
		aop->grad[j] 	= -1 * grd;

		// VecView(pro,	PETSC_VIEWER_STDOUT_WORLD );     CHKERRQ(ierr);

		// Print result
		PetscPrintf(PETSC_COMM_WORLD,"%D.Gradient = %g ; CurPar = %d ; CurPhase = %d\n",j+1,aop->grad[j],CurPar,CurPhase);

		// Set all the used parameter back to its original value
		if (CurPar==1) {
			nl->pc->pm->jr->phases[CurPhase].rho = rhoIni;
		}else if (CurPar==2) {
			nl->pc->pm->jr->phases[CurPhase].rho_n = rho_nIni;
		}else if (CurPar==3) {
			nl->pc->pm->jr->phases[CurPhase].rho_c = rho_cIni;
		}else if (CurPar==4) {
			nl->pc->pm->jr->phases[CurPhase].beta = betaIni;
		}else if (CurPar==5) {
			nl->pc->pm->jr->phases[CurPhase].K = KIni;
		}else if (CurPar==6) {
			nl->pc->pm->jr->phases[CurPhase].Kp = KpIni;
		}else if (CurPar==7) {
			nl->pc->pm->jr->phases[CurPhase].G = GIni;
		}else if (CurPar==8) {
			nl->pc->pm->jr->phases[CurPhase].Bd = BdIni;
		}else if (CurPar==9) {
			nl->pc->pm->jr->phases[CurPhase].Ed = EdIni;
		}else if (CurPar==10) {
			nl->pc->pm->jr->phases[CurPhase].Vd = VdIni;
		}else if (CurPar==11) {
			nl->pc->pm->jr->phases[CurPhase].Bn = BnIni;
		}else if (CurPar==12) {
			nl->pc->pm->jr->phases[CurPhase].n = nIni;
		}else if (CurPar==13) {
			nl->pc->pm->jr->phases[CurPhase].En = EnIni;
		}else if (CurPar==14) {
			nl->pc->pm->jr->phases[CurPhase].Vn = VnIni;
		}else if (CurPar==15) {
			nl->pc->pm->jr->phases[CurPhase].taup = taupIni;
		}else if (CurPar==16) {
			nl->pc->pm->jr->phases[CurPhase].gamma = gammaIni;
		}else if (CurPar==17) {
			nl->pc->pm->jr->phases[CurPhase].q = qIni;
		}/*else if (CurPar==18) {
			nl->pc->pm->jr->phases[CurPhase].eta = etaIni;
		}*/

	}

	// Print the solution variable at the user defined index
	for (i=0; i<user->AdjointNumInd; i++) {
		PetscPrintf(PETSC_COMM_WORLD,"Computation variable = %g \n",vz[i]);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful\n------------------------------------------\n");

	// Clean
	ierr = VecDestroy(&rpl);
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&pro);
	ierr = VecDestroy(&proX);
	ierr = VecDestroy(&proY);
	ierr = VecDestroy(&proZ);
	ierr = VecDestroy(&res);
	ierr = VecDestroy(&drdp);
	ierr = VecDestroy(&Perturb_vec);

	PetscFunctionReturn(0);
}
