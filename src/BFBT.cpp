/*
 * BFBT.cpp
 *
 *  Created on: 04.03.2020
 *      Author: daniel
 */

//---------------------------------------------------------------------------
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matrix.h"
#include "multigrid.h"
#include "lsolve.h"
#include "BFBT.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTCreate"
PetscErrorCode PMatBFBTCreate(PMat pm)
{
	PMatBlock      *P;
	JacRes         *jr;
	FDSTAG         *fs;
	DOFIndex       *dof;
	PetscInt        lnv, stv;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);

	P   = (PMatBlock*)pm->data;
	jr  = pm->jr;
	fs  = jr->fs;
	dof = &fs->dof;
	lnv = dof->lnv;
	stv = dof->stv;

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create DMDA
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &P->DA_P); CHKERRQ(ierr);

	// create matrix
	ierr = DMCreateMatrix(P->DA_P, &P->K); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// create scaling matrix
	ierr = MatAIJCreateDiag(lnv, stv, &P->C); CHKERRQ(ierr);

	// allocate work vectors
	ierr = VecDuplicate(P->xv, &P->wv0);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv2);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv3);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv4);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv5);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv7);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xp, &P->wp1);  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xp, &P->wp6);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTAssemble"
PetscErrorCode PMatBFBTAssemble(PMat pm)
{
	PMatBlock  *P;
	JacRes     *jr;
	FDSTAG     *fs;
	BCCtx      *bc;
	Vec         lvEtaCen;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, mcx, mcy, mcz, iter;
	PetscScalar bEtaX, fEtaX, bEtaY, fEtaY, bEtaZ, fEtaZ, eta;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz, dx, dy, dz;
	PetscScalar v[7], cfp[6]/*, cfv[6]*/;
	MatStencil  row[1], col[7];
	PetscScalar ***lEta, ***bcvx,  ***bcvy,  ***bcvz, ***bcp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context variables
	P 	 = (PMatBlock*)pm->data;
	jr   = pm->jr;
	fs   = jr->fs;
	bc   = jr->bc;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	mx  = fs->dsx.tnods - 1;
	my  = fs->dsy.tnods - 1;
	mz  = fs->dsz.tnods - 1;

	//===============
	// CELL VISCOSITY
	//===============

	ierr = DMGetLocalVector(fs->DA_CEN, &lvEtaCen); CHKERRQ(ierr);

	ierr = VecZeroEntries(lvEtaCen); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		lEta[k][j][i] = jr->svCell[iter++].svDev.eta;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, lvEtaCen)

	//===============
	// SCALING MATRIX
	//===============

	ierr = DMDAVecGetArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	// set iterator
	iter = fs->dof.stv;

	// X-points
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;

		eta = (lEta[k][j][I1] + lEta[k][j][I2])/2.0;

		ierr = MatSetValue(P->C, iter, iter, sqrt(eta), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Y-points
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		eta = (lEta[k][J1][i] + lEta[k][J2][i])/2.0;

		ierr = MatSetValue(P->C, iter, iter, sqrt(eta), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Z-points
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		eta = (lEta[K1][j][i] + lEta[K2][j][i])/2.0;

		ierr = MatSetValue(P->C, iter, iter, sqrt(eta), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// assemble scaling matrix
	ierr = MatAIJAssemble(P->C, bc->vNumSPC, bc->vSPCList, 1.0); CHKERRQ(ierr);

	//=======================
	// PRECONDITIONING MATRIX
	//=======================

	// clear matrix
	ierr = MatZeroEntries(P->K); CHKERRQ(ierr);

	// access bc vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// check index bounds and TPC multipliers
		Im1 = i-1; cfp[0] = 1.0; if(Im1 < 0)   { Im1++; if(bcp[k][j][i-1] != DBL_MAX) cfp[0] = -1.0; }
		Ip1 = i+1; cfp[1] = 1.0; if(Ip1 > mcx) { Ip1--; if(bcp[k][j][i+1] != DBL_MAX) cfp[1] = -1.0; }
		Jm1 = j-1; cfp[2] = 1.0; if(Jm1 < 0)   { Jm1++; if(bcp[k][j-1][i] != DBL_MAX) cfp[2] = -1.0; }
		Jp1 = j+1; cfp[3] = 1.0; if(Jp1 > mcy) { Jp1--; if(bcp[k][j+1][i] != DBL_MAX) cfp[3] = -1.0; }
		Km1 = k-1; cfp[4] = 1.0; if(Km1 < 0)   { Km1++; if(bcp[k-1][j][i] != DBL_MAX) cfp[4] = -1.0; }
		Kp1 = k+1; cfp[5] = 1.0; if(Kp1 > mcz) { Kp1--; if(bcp[k+1][j][i] != DBL_MAX) cfp[5] = -1.0; }

		// compute viscosity scaling coefficients
		eta   = lEta[k][j][i];
		bEtaX = sqrt((eta + lEta[k][j][Im1])/2.0);   fEtaX = sqrt((eta + lEta[k][j][Ip1])/2.0);
		bEtaY = sqrt((eta + lEta[k][Jm1][i])/2.0);   fEtaY = sqrt((eta + lEta[k][Jp1][i])/2.0);
		bEtaZ = sqrt((eta + lEta[Km1][j][i])/2.0);   fEtaZ = sqrt((eta + lEta[Kp1][j][i])/2.0);

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC & SPC multipliers
		v[0] = -bEtaX/bdx/dx*cfp[0];
		v[1] = -fEtaX/fdx/dx*cfp[1];
		v[2] = -bEtaY/bdy/dy*cfp[2];
		v[3] = -fEtaY/fdy/dy*cfp[3];
		v[4] = -bEtaZ/bdz/dz*cfp[4];
		v[5] = -fEtaZ/fdz/dz*cfp[5];
		v[6] = (bEtaX/bdx + fEtaX/fdx)/dx
		+      (bEtaY/bdy + fEtaY/fdy)/dy
		+      (bEtaZ/bdz + fEtaZ/fdz)/dz;

		// set matrix coefficients
		ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// assemble preconditioning matrix
	ierr = MatAIJAssemble(P->K,  bc->pNumSPC, bc->pSPCList, 1.0); CHKERRQ(ierr);

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	// restore viscosity vector
	ierr = DMRestoreLocalVector (fs->DA_CEN, &lvEtaCen); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTDestroy"
PetscErrorCode PMatBFBTDestroy(PMat pm)
{
	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);

	P = (PMatBlock*)pm->data;

	ierr = DMDestroy (&P->DA_P); CHKERRQ(ierr);
	ierr = MatDestroy(&P->K); 	 CHKERRQ(ierr);
	ierr = MatDestroy(&P->C); 	 CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv0);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wp1);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv2);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv3);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv4);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv5);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wp6);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv7);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFBTApply"
PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y)
{
	PCStokes    pc;
//	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

//	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;


	// ACHTUNG! placeholder
	ierr = MatMult(P->iS, x, y);     CHKERRQ(ierr); // xp = (S^-1)*rp


/*
		//=======================
		// BLOCK w-BFBT
		//=======================

		PMat pm;
		JacRes *jr;
		pm = pc->pm;
		jr = pm->jr;

		//assemble C                                           (get global viscosity like residual in JacRes)
		ierr = CopyViscosityToScalingVector(jr->eta_gfx, jr->eta_gfy, jr->eta_gfz, P->C); CHKERRQ(ierr);

		// rv = f
		// wp = B*A⁻1*rv
		ierr = KSPSolve(bf->vksp, P->rv, P->wv0); 		CHKERRQ(ierr); // wv0 = (Avv⁻1)*rv     | A=Avv | B=Apv | B^T=Avp
		ierr = MatMult(P->Apv, P->wv0, P->wp); 			CHKERRQ(ierr); // wp  = Apv*wv0

		// p = S⁻1*wp        S⁻1 = (BCB^T)⁻1 * BCACB^T * (BCB^T)⁻1 = K⁻1 * BCACB^T * K⁻1
		// K = BCB^T
		ierr = KSPSolve(bf->pksp, P->wp, P->wp1); 		CHKERRQ(ierr); // wp1 = K⁻1*wp     <=> K*wp1 = wp
		ierr = MatMult(P->Avp, P->wp1, P->wv2); 		CHKERRQ(ierr); // wv2 = Avp*wp1
		ierr = VecPointwiseMult(P->wv3, P->C, P->wv2); 	CHKERRQ(ierr); // wv3 = C*wv2
		ierr = MatMult(P->Avv, P->wv3, P->wv4); 		CHKERRQ(ierr); // wv4 = Avv * wv3
		ierr = VecPointwiseMult(P->wv5, P->C, P->wv4); 	CHKERRQ(ierr); // wv5 = C*wv4
		ierr = MatMult(P->Apv, P->wv5, P->wp6); 		CHKERRQ(ierr); // wp6 = Apv*wv5
		ierr = KSPSolve(bf->pksp, P->wp6, P->xp); 		CHKERRQ(ierr); // xp  = K⁻1*wp6     <=> K*xp = wp6

		// u = A⁻1*(wv-B^T*p)
		ierr = MatMult(P->Avp, P->xp, P->wv7); 			CHKERRQ(ierr); // wv7 = B^T*xp
		ierr = VecWAXPY(P->rv, -1.0, P->wv7, P->wv); 	CHKERRQ(ierr); // rv  = wv-wv7
		ierr = KSPSolve(bf->vksp, P->rv, P->xv); 		CHKERRQ(ierr); // xv  = (A⁻1)*rv



 */


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------











/*

	PMatBlock *P;
	JacRes 	  *jr;

	FDSTAG     *fs;
	BCCtx      *bc;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bvx, fvx, bvy, fvy, bvz, fvz;
	PetscScalar sbvx, sfvx, sbvy, sfvy, sbvz, sfvz;
	PetscScalar ibvx, ifvx, ibvy, ifvy, ibvz, ifvz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;

	PetscScalar v[7], cf[6];
	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcvx, ***bcvy, ***bcvz, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);


	// access residual context variables
	P 	 = (PMatBlock*)pm->data;
	jr   = pm->jr;
	fs   = jr->fs;
	bc   = jr->bc;
	num  = bc->pNumSPC;
	list = bc->pSPCList;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

//	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_VISC_TOTAL)

	// clear matrix coefficients
	ierr = MatZeroEntries(P->K); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, bc->bcvz,  &bcvz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcvx[k][j][i-1] != DBL_MAX) cf[0] = 0.0; }
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcvx[k][j][i+1] != DBL_MAX) cf[1] = 0.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcvy[k][j-1][i] != DBL_MAX) cf[2] = 0.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcvy[k][j+1][i] != DBL_MAX) cf[3] = 0.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcvz[k-1][j][i] != DBL_MAX) cf[4] = 0.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcvz[k+1][j][i] != DBL_MAX) cf[5] = 0.0; }

		// compute viscosity on cell faces
		bvx = (lk[k][j][i] + lk[k][j][Im1])/2.0;      fvx = (lk[k][j][i] + lk[k][j][Ip1])/2.0;
		bvy = (lk[k][j][i] + lk[k][Jm1][i])/2.0;      fvy = (lk[k][j][i] + lk[k][Jp1][i])/2.0;
		bvz = (lk[k][j][i] + lk[Km1][j][i])/2.0;      fvz = (lk[k][j][i] + lk[Kp1][j][i])/2.0;

		// compute square root of the viscosity
		sbvx = sqrt(bvx);	sfvx = sqrt(fvx);
		sbvy = sqrt(bvy);	sfvy = sqrt(fvy);
		sbvz = sqrt(bvz);	sfvz = sqrt(fvz);

		// compute inverse of the square root
		ibvx = 1/sbvx;		ifvx = 1/sfvx;
		ibvy = 1/sbvy;		ifvy = 1/sfvy;
		ibvz = 1/sbvz;		ifvz = 1/sfvz;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC multipliers
		v[0] = -ibvx/bdx/dx*cf[0];
		v[1] = -ifvx/fdx/dx*cf[1];
		v[2] = -ibvy/bdy/dy*cf[2];
		v[3] = -ifvy/fdy/dy*cf[3];
		v[4] = -ibvz/bdz/dz*cf[4];
		v[5] = -ifvz/fdz/dz*cf[5];
		v[6] =  (ibvx/bdx + ifvx/fdx)/dx
		+       (ibvy/bdy + ifvy/fdy)/dy
		+       (ibvz/bdz + ifvz/fdz)/dz;

		// set matrix coefficients
		//ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

		// NOTE! since only TPC are active, no SPC modification is necessary
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz);  CHKERRQ(ierr);

	// assemble temperature matrix
	ierr = MatAIJAssemble(P->K, num, list, 1.0); CHKERRQ(ierr);

*/





