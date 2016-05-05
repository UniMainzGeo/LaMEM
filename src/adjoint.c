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
// Objective function    F(x,x(p)) = P * x^T                // p = parameter ; x = converged solution ; Projection vector containing 1 at comparison points
// Derivative I          dF/dx     = P
// Adjoint operation     psi       = J^-T * dF/dx           // J = converged Jacobain matrix
// Derivative II         dr/dp     = [r(p+h) - r(p)]/h      // finite difference approximation of derivative of residual r vs parameter
// Gradients             dF/dp     = -psi^T * dr/dp
//
// USAGE:
// In your *.dat file you need to define:
// ComputeAdjointGradients	       = 1	                    // 1 = Compute gradients 0 = NO
// AdjointIndex					   = 99 100 101             // Array containing the indices in the solution vector at which you want to compute the gradients in space (best to find them in matlab when creating the setup)
// AdjointPhases                   = 1 2 1 2                // Array (same length as AdjointParameters) containing the phase of the parameter
// AdjointParameters               = 2 1 1 1                // Array (same length as AdjointPhases) containing the parameter corresponding to the phase
//
//                              Density         Elasticity  Diff creep    Dis creep    Peierl creep
// Possible parameters: (rho rho_n rho_c beta)   (K Kp G)   (Bd Ed Vd)  (Bn n  En Vn) (taup gamma q )
//                      ( 1     2     3    4 )   (5  6 7)   ( 8  9 10)  (11 12 13 14) ( 15    16  17)
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "matFree.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "tools.h"
#include "adjoint.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointDestroy"
PetscErrorCode AdjointDestroy(AdjGrad *aop)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Nothing to clear so far (maybe for the future)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CreateAdjoint"
PetscErrorCode CreateAdjoint(JacRes *jr, UserCtx *user, AdjGrad *aop, NLSol *nl,SNES *snes)
{

	/* TODO:
		- It's only valid vor visco-elasicity since the Jacobian is assumed to be symmetric
		- For some strange reason gives different results in parallel (maybe the summation of the vector is worng in parallel?)
	*/

	PetscPrintf(PETSC_COMM_WORLD,"Computation of adjoint gradients\n");

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Initialize all variables
	KSP                 ksp;
	KSPConvergedReason  reason;
	PetscInt            i, j, sizeSol;
	PetscScalar         p, Perturb, *tempres, *temprpl, *temppsi;       // Temporary vectors for the computation
	PetscScalar         rhoIni, rho_nIni, rho_cIni, betaIni, KIni, KpIni, GIni, BdIni, EdIni, VdIni, BnIni, nIni, EnIni, VnIni, taupIni, gammaIni, qIni;   // Initialize parameters for which the gradients can be computed
	Vec 				rpl,sol,psi,pro;

	// Set perturbation paramter for the finite differences
	Perturb = 1e-6;
	p       = 1.0;

	// Get size of the solution vector
	ierr = VecGetSize(jr->gsol, &sizeSol);

	// Initialize temporary computation arrays
	PetscScalar tempdrdp[sizeSol], tempap[sizeSol];

	// Create all needed vectors in the same size as the solution vector
	ierr = VecDuplicate(jr->gsol, &pro); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &psi); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &rpl); CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gsol, &sol); CHKERRQ(ierr);
	ierr = VecCopy(jr->gsol,sol); CHKERRQ(ierr);

	// Put 1 into the Projection vector where the user defined the computation indices (dF/dx = P)
	for (i=0; i<=user->AdjointNumInd; i++) {
		ierr = VecSetValues(pro,1,&user->AdjointIndex[i],&p,INSERT_VALUES);
	}
	ierr = VecAssemblyBegin(pro);
	ierr = VecAssemblyEnd(pro);

	// Solve the main adjoint equation with KSP (psi = J^-T * dF/dx)
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,nl->J,nl->P);
	KSPSetType(ksp,KSPCG);
	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	KSPSolve(ksp,pro,psi);
	KSPGetConvergedReason(ksp,&reason);

	// Extract global residual and psi vector for computation
	VecGetArray(jr->gres,&tempres);
	VecGetArray(psi,&temppsi);

	//===============
	// Parameter loop
	//===============
	for(j=0; j<user->AdjointNumPar; j++){

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
		}

		// Compute residual with the converged Jacobian and preconditioner (dr/dp = [r(p+h) - r(p)]/h)
		ierr = FormResidual(snes, sol, rpl, nl); CHKERRQ(ierr);

		// Exract new residual
		VecGetArray(rpl,&temprpl);

		// Compute the gradient (dF/dp = psi^T * dr/dp [minus is added in the next line])
		for (i=0; i<=sizeSol; i++) {
			tempdrdp[i] 	= (tempres[i] - temprpl[i])/Perturb;
			tempap[i]   	= temppsi[i] * tempdrdp[i];
		}

		// Sum the gradient vector up to get the gradient
		for (i=0; i<=sizeSol; i++) {
			aop->grad[j] += tempap[i];
		}

		// Save gradient
		aop->grad[j] 	= -aop->grad[j];    // Include the minus sign of the equation (dF/dp = -dF/dp)

		// Rebuild the residual vector
		VecRestoreArray(rpl,&temprpl);

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
		}
	}

	PetscPrintf(PETSC_COMM_WORLD,"Computation was succesful\n------------------------------------------\n");

	// Restore parameter independent vectors
	VecRestoreArray(jr->gres,&tempres);
	VecRestoreArray(psi,&temppsi);

	// Clean
	ierr = VecDestroy(&rpl);
	ierr = VecDestroy(&psi);
	ierr = VecDestroy(&sol);
	ierr = VecDestroy(&pro);

	PetscFunctionReturn(0);
}
