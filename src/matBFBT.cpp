/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...................   wBFBT PRECONDITIONER MATRIX   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matData.h"
#include "matrix.h"
//---------------------------------------------------------------------------
PetscErrorCode wBFBTCreate(wBFBTData *P, MatData *md)
{
	FDSTAG         *fs;
	DOFIndex       *dof;
	DMBoundaryType BC_TYPE_X;
	PetscInt       periodic;
	PetscInt       lnv, stv;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// store evaluation context
	P->md = md;

	// access context variables
	fs  = md->fs;
	dof = &fs->dof;
	lnv = dof->lnv;
	stv = dof->stv;

	// set periodic flag
	periodic = fs->periodic;

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// set boundary type in x direction
	if(periodic) { BC_TYPE_X = DM_BOUNDARY_PERIODIC; }
	else         { BC_TYPE_X = DM_BOUNDARY_NONE;     }

	// create DMDA
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		BC_TYPE_X, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &P->DA_P); CHKERRQ(ierr);

	// set proper interpolation type for multigrid
	ierr = DMDASetInterpolationType(P->DA_P, DMDA_Q0); CHKERRQ(ierr);

	// create matrix
	ierr = DMCreateMatrix(P->DA_P, &P->K); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// create scaling matrix
	ierr = MatAIJCreateDiag(lnv, stv, &P->C); CHKERRQ(ierr);

	// allocate work vector
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &P->w); CHKERRQ(ierr);
	ierr = VecSetFromOptions(P->w);                                     CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode wBFBTDestroy(wBFBTData *P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = DMDestroy (&P->DA_P); CHKERRQ(ierr);
	ierr = MatDestroy(&P->K); 	 CHKERRQ(ierr);
	ierr = MatDestroy(&P->C); 	 CHKERRQ(ierr);
	ierr = VecDestroy(&P->w);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode wBFBTAssemble(wBFBTData *P)
{
	MatData    *md;
	FDSTAG     *fs;
	MatStencil  row[1], col[7];
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    mnx, mny, mnz, mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v[7], d[6], g[6], cfe[6], cfp[6], cfv[6];
	PetscScalar ***eta, ***bcvx,  ***bcvy,  ***bcvz, ***bcp;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz, dx, dy, dz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context variables
	md = P->md;
	fs = md->fs;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	mnx  = fs->dsx.tnods - 1;
	mny  = fs->dsy.tnods - 1;
	mnz  = fs->dsz.tnods - 1;

	//===============
	// SCALING MATRIX
	//===============

	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta, &eta); CHKERRQ(ierr);

	// set iterator
	iter = fs->dof.stv;

	// X-points
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		I1 = i;   if(I1 == mnx) I1--;
		I2 = i-1; if(I2 == -1)  I2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((eta[k][j][I1] + eta[k][j][I2])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Y-points
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		J1 = j;   if(J1 == mny) J1--;
		J2 = j-1; if(J2 == -1)  J2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((eta[k][J1][i] + eta[k][J2][i])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Z-points
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		K1 = k;   if(K1 == mnz) K1--;
		K2 = k-1; if(K2 == -1)  K2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((eta[K1][j][i] + eta[K2][j][i])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// assemble scaling matrix
	ierr = MatAIJAssemble(P->C, md->vNumSPC, md->vSPCListMat, 1.0); CHKERRQ(ierr);

	//=======================
	// PRECONDITIONING MATRIX
	//=======================

	// clear matrix
	ierr = MatZeroEntries(P->K); CHKERRQ(ierr);

	// access boundary constraints
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,  &bcp);  CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// check index bounds, get pressure TPC multipliers
		Im1 = i-1; cfp[0] = 1.0; if(Im1 < 0)   { Im1++; if(bcp[k][j][i-1] != DBL_MAX) cfp[0] = -1.0; }
		Ip1 = i+1; cfp[1] = 1.0; if(Ip1 > mcx) { Ip1--; if(bcp[k][j][i+1] != DBL_MAX) cfp[1] = -1.0; }
		Jm1 = j-1; cfp[2] = 1.0; if(Jm1 < 0)   { Jm1++; if(bcp[k][j-1][i] != DBL_MAX) cfp[2] = -1.0; }
		Jp1 = j+1; cfp[3] = 1.0; if(Jp1 > mcy) { Jp1--; if(bcp[k][j+1][i] != DBL_MAX) cfp[3] = -1.0; }
		Km1 = k-1; cfp[4] = 1.0; if(Km1 < 0)   { Km1++; if(bcp[k-1][j][i] != DBL_MAX) cfp[4] = -1.0; }
		Kp1 = k+1; cfp[5] = 1.0; if(Kp1 > mcz) { Kp1--; if(bcp[k+1][j][i] != DBL_MAX) cfp[5] = -1.0; }

		// get velocity SPC multipliers
		cfv[0] = 1.0; if(bcvx[k][j][i]   != DBL_MAX) cfv[0] = 0.0;
		cfv[1] = 1.0; if(bcvx[k][j][i+1] != DBL_MAX) cfv[1] = 0.0;
		cfv[2] = 1.0; if(bcvy[k][j][i]   != DBL_MAX) cfv[2] = 0.0;
		cfv[3] = 1.0; if(bcvy[k][j+1][i] != DBL_MAX) cfv[3] = 0.0;
		cfv[4] = 1.0; if(bcvz[k][j][i]   != DBL_MAX) cfv[4] = 0.0;
		cfv[5] = 1.0; if(bcvz[k+1][j][i] != DBL_MAX) cfv[5] = 0.0;

		// get viscosity scaling coefficients
		cfe[0] = sqrt((eta[k][j][i] + eta[k][j][Im1])/2.0);
		cfe[1] = sqrt((eta[k][j][i] + eta[k][j][Ip1])/2.0);
		cfe[2] = sqrt((eta[k][j][i] + eta[k][Jm1][i])/2.0);
		cfe[3] = sqrt((eta[k][j][i] + eta[k][Jp1][i])/2.0);
		cfe[4] = sqrt((eta[k][j][i] + eta[Km1][j][i])/2.0);
		cfe[5] = sqrt((eta[k][j][i] + eta[Kp1][j][i])/2.0);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get gradient vector & divergence vector scaled by viscosity and velocity SPC multipliers
		g[0] =  1.0/bdx;   d[0] =  cfv[0]/dx/cfe[0];
		g[1] = -1.0/fdx;   d[1] = -cfv[1]/dx/cfe[1];
		g[2] =  1.0/bdy;   d[2] =  cfv[2]/dy/cfe[2];
		g[3] = -1.0/fdy;   d[3] = -cfv[3]/dy/cfe[3];
		g[4] =  1.0/bdz;   d[4] =  cfv[4]/dz/cfe[4];
		g[5] = -1.0/fdz;   d[5] = -cfv[5]/dz/cfe[5];

		// compute matrix row, apply pressure TPC multipliers
		v[0] = -g[0]*d[0]*cfp[0];
		v[1] = -g[1]*d[1]*cfp[1];
		v[2] = -g[2]*d[2]*cfp[2];
		v[3] = -g[3]*d[3]*cfp[3];
		v[4] = -g[4]*d[4]*cfp[4];
		v[5] = -g[5]*d[5]*cfp[5];
		v[6] =  g[0]*d[0] + g[1]*d[1] + g[2]*d[2] + g[3]*d[3] + g[4]*d[4] + g[5]*d[5];

		// set diagonal for fully constrained cells
		if(!v[6]) v[6] = 1.0;

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set matrix coefficients
		ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// assemble preconditioning matrix
	ierr = MatAIJAssemble(P->K, md->pNumSPC, md->pSPCListMat, 1.0); CHKERRQ(ierr);

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta, &eta); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,  &bcp);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

