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
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "multigrid.h"
#include "matData.h"
#include "matrix.h"
#include "matFree.h"
#include "tools.h"
//---------------------------------------------------------------------------
// MG Level functions
//---------------------------------------------------------------------------
PetscErrorCode MGLevelCreate(MGLevel *lvl, MGLevel *fine, MatData *md)
{
	PetscInt ln=0, lnfine=0;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(!fine)
	{
		// setup top level
		lvl->md = md;
		lvl->R  = NULL;
		lvl->P  = NULL;
	}
	else
	{
		// allocate evaluation context
		ierr = PetscMalloc(sizeof(MatData), &lvl->md); CHKERRQ(ierr);

		// clear memory
		ierr = PetscMemzero(lvl->md, sizeof(MatData)); CHKERRQ(ierr);

		// coarsen matrix context
		ierr = MatDataCoarsen(lvl->md, fine->md); CHKERRQ(ierr);

		// get matrix sizes
		if     (lvl->md->idxmod == _IDX_COUPLED_) { ln = lvl->md->fs->dof.ln;  lnfine = fine->md->fs->dof.ln;  }
		else if(lvl->md->idxmod == _IDX_BLOCK_)   { ln = lvl->md->fs->dof.lnv; lnfine = fine->md->fs->dof.lnv; }

		// WARNING! CONSTANT SIZE PREALLOCATION (ADD VARIABLE PREALLOCATION)

		// preallocate restriction & prolongation matrices
		ierr = MatAIJCreate(ln,     lnfine, 12, NULL, 4, NULL, &lvl->R); CHKERRQ(ierr);
		ierr = MatAIJCreate(lnfine, ln,     8,  NULL, 7, NULL, &lvl->P); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGLevelDestroy(MGLevel *lvl)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(lvl->md->coarse)
	{
		ierr = MatDataDestroy(lvl->md); CHKERRQ(ierr);
		ierr = PetscFree     (lvl->md); CHKERRQ(ierr);

		ierr = MatDestroy(&lvl->R); CHKERRQ(ierr);
		ierr = MatDestroy(&lvl->P); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine)
{
	Mat         R;
	MatData     *mdlvl, *mdfine;
	PetscScalar v[12], bc[12], vs[12];
	PetscInt    idx[12];
	PetscInt    row, I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	R      = lvl->R;
	mdlvl  = lvl->md;
	mdfine = fine->md;

	// clear restriction matrix coefficients
	ierr = MatZeroEntries(R); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(mdfine->fs->DA_X,   mdfine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Y,   mdfine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Z,   mdfine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_CEN, mdfine->ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(mdfine->fs->DA_X,   mdfine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Y,   mdfine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Z,   mdfine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_CEN, mdfine->bcp,  &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(mdlvl->fs->DA_X,   mdlvl->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Y,   mdlvl->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Z,   mdlvl->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_CEN, mdlvl->bcp,  &cbcp);  CHKERRQ(ierr);

	// get global index of the first row in coarse grid
	if     (mdlvl->idxmod == _IDX_COUPLED_) { row = mdlvl->fs->dof.st;  }
	else if(mdlvl->idxmod == _IDX_BLOCK_)   { row = mdlvl->fs->dof.stv; }

	// set velocity stencil weights
	vs[0]  = 1.0/16.0;
	vs[1]  = 1.0/16.0;
	vs[2]  = 1.0/16.0;
	vs[3]  = 1.0/16.0;
	vs[4]  = 1.0/8.0;
	vs[5]  = 1.0/8.0;
	vs[6]  = 1.0/8.0;
	vs[7]  = 1.0/8.0;
	vs[8]  = 1.0/16.0;
	vs[9]  = 1.0/16.0;
	vs[10] = 1.0/16.0;
	vs[11] = 1.0/16.0;

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(mdlvl->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)ivx[K  ][J  ][I-1];
		idx[1]  = (PetscInt)ivx[K  ][J+1][I-1];
		idx[2]  = (PetscInt)ivx[K+1][J  ][I-1];
		idx[3]  = (PetscInt)ivx[K+1][J+1][I-1];
		idx[4]  = (PetscInt)ivx[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivx[K  ][J+1][I  ];
		idx[6]  = (PetscInt)ivx[K+1][J  ][I  ];
		idx[7]  = (PetscInt)ivx[K+1][J+1][I  ];
		idx[8]  = (PetscInt)ivx[K  ][J  ][I+1];
		idx[9]  = (PetscInt)ivx[K  ][J+1][I+1];
		idx[10] = (PetscInt)ivx[K+1][J  ][I+1];
		idx[11] = (PetscInt)ivx[K+1][J+1][I+1];

		// get fine grid boundary conditions
		bc[0]   =         fbcvx[K  ][J  ][I-1];
		bc[1]   =         fbcvx[K  ][J+1][I-1];
		bc[2]   =         fbcvx[K+1][J  ][I-1];
		bc[3]   =         fbcvx[K+1][J+1][I-1];
		bc[4]   =         fbcvx[K  ][J  ][I  ];
		bc[5]   =         fbcvx[K  ][J+1][I  ];
		bc[6]   =         fbcvx[K+1][J  ][I  ];
		bc[7]   =         fbcvx[K+1][J+1][I  ];
		bc[8]   =         fbcvx[K  ][J  ][I+1];
		bc[9]   =         fbcvx[K  ][J+1][I+1];
		bc[10]  =         fbcvx[K+1][J  ][I+1];
		bc[11]  =         fbcvx[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvx[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(mdlvl->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)ivy[K  ][J-1][I  ];
		idx[1]  = (PetscInt)ivy[K  ][J-1][I+1];
		idx[2]  = (PetscInt)ivy[K+1][J-1][I  ];
		idx[3]  = (PetscInt)ivy[K+1][J-1][I+1];
		idx[4]  = (PetscInt)ivy[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivy[K  ][J  ][I+1];
		idx[6]  = (PetscInt)ivy[K+1][J  ][I  ];
		idx[7]  = (PetscInt)ivy[K+1][J  ][I+1];
		idx[8]  = (PetscInt)ivy[K  ][J+1][I  ];
		idx[9]  = (PetscInt)ivy[K  ][J+1][I+1];
		idx[10] = (PetscInt)ivy[K+1][J+1][I  ];
    	idx[11] = (PetscInt)ivy[K+1][J+1][I+1];

    	// get fine grid boundary conditions
		bc[0]   =         fbcvy[K  ][J-1][I  ];
		bc[1]   =         fbcvy[K  ][J-1][I+1];
		bc[2]   =         fbcvy[K+1][J-1][I  ];
		bc[3]   =         fbcvy[K+1][J-1][I+1];
		bc[4]   =         fbcvy[K  ][J  ][I  ];
		bc[5]   =         fbcvy[K  ][J  ][I+1];
		bc[6]   =         fbcvy[K+1][J  ][I  ];
		bc[7]   =         fbcvy[K+1][J  ][I+1];
		bc[8]   =         fbcvy[K  ][J+1][I  ];
		bc[9]   =         fbcvy[K  ][J+1][I+1];
		bc[10]  =         fbcvy[K+1][J+1][I  ];
    	bc[11]  =         fbcvy[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvy[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(mdlvl->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)ivz[K-1][J  ][I  ];
		idx[1]  = (PetscInt)ivz[K-1][J  ][I+1];
		idx[2]  = (PetscInt)ivz[K-1][J+1][I  ];
		idx[3]  = (PetscInt)ivz[K-1][J+1][I+1];
		idx[4]  = (PetscInt)ivz[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivz[K  ][J  ][I+1];
		idx[6]  = (PetscInt)ivz[K  ][J+1][I  ];
		idx[7]  = (PetscInt)ivz[K  ][J+1][I+1];
		idx[8]  = (PetscInt)ivz[K+1][J  ][I  ];
		idx[9]  = (PetscInt)ivz[K+1][J  ][I+1];
		idx[10] = (PetscInt)ivz[K+1][J+1][I  ];
    	idx[11] = (PetscInt)ivz[K+1][J+1][I+1];

    	// get fine grid boundary conditions
		bc[0]   =         fbcvz[K-1][J  ][I  ];
		bc[1]   =         fbcvz[K-1][J  ][I+1];
		bc[2]   =         fbcvz[K-1][J+1][I  ];
		bc[3]   =         fbcvz[K-1][J+1][I+1];
		bc[4]   =         fbcvz[K  ][J  ][I  ];
		bc[5]   =         fbcvz[K  ][J  ][I+1];
		bc[6]   =         fbcvz[K  ][J+1][I  ];
		bc[7]   =         fbcvz[K  ][J+1][I+1];
		bc[8]   =         fbcvz[K+1][J  ][I  ];
		bc[9]   =         fbcvz[K+1][J  ][I+1];
		bc[10]  =         fbcvz[K+1][J+1][I  ];
    	bc[11]  =         fbcvz[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvz[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(mdlvl->idxmod == _IDX_COUPLED_)
	{
		// set pressure weights
		vs[0] = 1.0/8.0;
		vs[1] = 1.0/8.0;
		vs[2] = 1.0/8.0;
		vs[3] = 1.0/8.0;
		vs[4] = 1.0/8.0;
		vs[5] = 1.0/8.0;
		vs[6] = 1.0/8.0;
		vs[7] = 1.0/8.0;

		//-----------------------
		// P-points (coarse grid)
		//-----------------------
		ierr = DMDAGetCorners(mdlvl->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get fine grid indices
			I = 2*i;
			J = 2*j;
			K = 2*k;

			// get fine grid stencil
			idx[0] = (PetscInt)ip[K  ][J  ][I  ];
			idx[1] = (PetscInt)ip[K  ][J  ][I+1];
			idx[2] = (PetscInt)ip[K  ][J+1][I  ];
			idx[3] = (PetscInt)ip[K  ][J+1][I+1];
			idx[4] = (PetscInt)ip[K+1][J  ][I  ];
			idx[5] = (PetscInt)ip[K+1][J  ][I+1];
			idx[6] = (PetscInt)ip[K+1][J+1][I  ];
			idx[7] = (PetscInt)ip[K+1][J+1][I+1];

			// get fine grid boundary conditions
			bc[0]  =         fbcp[K  ][J  ][I  ];
			bc[1]  =         fbcp[K  ][J  ][I+1];
			bc[2]  =         fbcp[K  ][J+1][I  ];
			bc[3]  =         fbcp[K  ][J+1][I+1];
			bc[4]  =         fbcp[K+1][J  ][I  ];
			bc[5]  =         fbcp[K+1][J  ][I+1];
			bc[6]  =         fbcp[K+1][J+1][I  ];
			bc[7]  =         fbcp[K+1][J+1][I+1];

			// setup row of restriction matrix
			getRowRestrict(cbcp[k][j][i], 8, idx, bc, v, vs);

			// store full matrix row
			ierr = MatSetValues(R, 1, &row, 8, idx, v, INSERT_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_X,   mdfine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Y,   mdfine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Z,   mdfine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_CEN, mdfine->ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(mdfine->fs->DA_X,   mdfine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Y,   mdfine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Z,   mdfine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_CEN, mdfine->bcp,  &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_X,    mdlvl->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Y,    mdlvl->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Z,    mdlvl->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_CEN,  mdlvl->bcp,  &cbcp);  CHKERRQ(ierr);

	// assemble restriction matrix
	ierr = MatAIJAssemble(R, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGLevelSetupProlong(MGLevel *lvl, MGLevel *fine)
{
	Mat         P;
	MatData     *mdlvl, *mdfine;
	PetscInt    n, idx[8];
	PetscScalar v[8], bc[8], vsf[8], vsr[4], *vs;
	PetscInt    row, I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	P      = lvl->P;
	mdlvl  = lvl->md;
	mdfine = fine->md;

	// clear prolongation matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);

	// access index vectors in coarse grid
	ierr = DMDAVecGetArray(mdlvl->fs->DA_X,   mdlvl->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Y,   mdlvl->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Z,   mdlvl->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_CEN, mdlvl->ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(mdfine->fs->DA_X,   mdfine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Y,   mdfine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_Z,   mdfine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdfine->fs->DA_CEN, mdfine->bcp,  &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(mdlvl->fs->DA_X,   mdlvl->bcvx,  &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Y,   mdlvl->bcvy,  &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_Z,   mdlvl->bcvz,  &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(mdlvl->fs->DA_CEN, mdlvl->bcp,   &cbcp);  CHKERRQ(ierr);

	// get global index of the first row in the fine grid
	if     (mdfine->idxmod == _IDX_COUPLED_) { row = mdfine->fs->dof.st;  }
	else if(mdfine->idxmod == _IDX_BLOCK_)   { row = mdfine->fs->dof.stv; }

	// set reduced velocity stencil coefficients (even)
	vsr[0] = 9.0/16.0;
	vsr[1] = 3.0/16.0;
	vsr[2] = 3.0/16.0;
	vsr[3] = 1.0/16.0;

	// set full velocity stencil coefficients (odd)
	vsf[0] = 9.0/32.0;
	vsf[1] = 3.0/32.0;
	vsf[2] = 3.0/32.0;
	vsf[3] = 1.0/32.0;
	vsf[4] = 9.0/32.0;
	vsf[5] = 3.0/32.0;
	vsf[6] = 3.0/32.0;
	vsf[7] = 1.0/32.0;

	//---------------------
	// X-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(mdfine->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(j % 2) J1 = J+1; else J1 = J-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivx[K ][J ][I];
		idx[1] = (PetscInt)ivx[K ][J1][I];
		idx[2] = (PetscInt)ivx[K1][J ][I];
		idx[3] = (PetscInt)ivx[K1][J1][I];
		bc [0] =         cbcvx[K ][J ][I];
		bc [1] =         cbcvx[K ][J1][I];
		bc [2] =         cbcvx[K1][J ][I];
		bc [3] =         cbcvx[K1][J1][I];

		if(i % 2)
		{
			// extend stencil (odd)
			n      = 8;
			vs     = vsf;
			I1     = I+1;
			idx[4] = (PetscInt)ivx[K ][J ][I1];
			idx[5] = (PetscInt)ivx[K ][J1][I1];
			idx[6] = (PetscInt)ivx[K1][J ][I1];
			idx[7] = (PetscInt)ivx[K1][J1][I1];
			bc [4] =         cbcvx[K ][J ][I1];
			bc [5] =         cbcvx[K ][J1][I1];
			bc [6] =         cbcvx[K1][J ][I1];
			bc [7] =         cbcvx[K1][J1][I1];
		}

		// setup row of prolongation matrix
		getRowProlong(row, fbcvx[k][j][i], n, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Y-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(mdfine->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivy[K ][J][I ];
		idx[1] = (PetscInt)ivy[K ][J][I1];
		idx[2] = (PetscInt)ivy[K1][J][I ];
		idx[3] = (PetscInt)ivy[K1][J][I1];
		bc [0] =         cbcvy[K ][J][I ];
		bc [1] =         cbcvy[K ][J][I1];
		bc [2] =         cbcvy[K1][J][I ];
		bc [3] =         cbcvy[K1][J][I1];

		if(j % 2)
		{
			// extend stencil (odd)
			n      = 8;
			vs     = vsf;
			J1     = J+1;
			idx[4] = (PetscInt)ivy[K ][J1][I ];
			idx[5] = (PetscInt)ivy[K ][J1][I1];
			idx[6] = (PetscInt)ivy[K1][J1][I ];
			idx[7] = (PetscInt)ivy[K1][J1][I1];
			bc [4] =         cbcvy[K ][J1][I ];
			bc [5] =         cbcvy[K ][J1][I1];
			bc [6] =         cbcvy[K1][J1][I ];
			bc [7] =         cbcvy[K1][J1][I1];
		}

		// setup row of prolongation matrix
		getRowProlong(row, fbcvy[k][j][i], n, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Z-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(mdfine->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(j % 2) J1 = J+1; else J1 = J-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivz[K][J ][I ];
		idx[1] = (PetscInt)ivz[K][J ][I1];
		idx[2] = (PetscInt)ivz[K][J1][I ];
		idx[3] = (PetscInt)ivz[K][J1][I1];
		bc [0] =         cbcvz[K][J ][I ];
		bc [1] =         cbcvz[K][J ][I1];
		bc [2] =         cbcvz[K][J1][I ];
		bc [3] =         cbcvz[K][J1][I1];

		if(k % 2)
		{
			// extend stencil (odd)
			n      = 8;
			vs     = vsf;
			K1     = K+1;
			idx[4] = (PetscInt)ivz[K1][J ][I ];
			idx[5] = (PetscInt)ivz[K1][J ][I1];
			idx[6] = (PetscInt)ivz[K1][J1][I ];
			idx[7] = (PetscInt)ivz[K1][J1][I1];
			bc [4] =         cbcvz[K1][J ][I ];
			bc [5] =         cbcvz[K1][J ][I1];
			bc [6] =         cbcvz[K1][J1][I ];
			bc [7] =         cbcvz[K1][J1][I1];
		}

		// setup row of prolongation matrix
		getRowProlong(row, fbcvz[k][j][i], n, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(mdfine->idxmod== _IDX_COUPLED_)
	{
		// set pressure interpolation stencil (direct injection)
		vsr[0] = 1.0;

		//---------------------
		// P-points (fine grid)
		//---------------------
		ierr = DMDAGetCorners(mdfine->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get coarse grid indices
			I = i/2;
			J = j/2;
			K = k/2;

			idx[0] = (PetscInt)ip[K][J][I];
			bc [0] =         cbcp[K][J][I];

			// setup row of prolongation matrix
			getRowProlong(row, fbcp[k][j][i], 1, bc, v, vsr);

			// store full matrix row
			ierr = MatSetValues(P, 1, &row, 1, idx, v, INSERT_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_X,  mdlvl->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Y,  mdlvl->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Z,  mdlvl->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_CEN,mdlvl->ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(mdfine->fs->DA_X,   mdfine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Y,   mdfine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_Z,   mdfine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdfine->fs->DA_CEN, mdfine->bcp,  &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_X,   mdlvl->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Y,   mdlvl->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_Z,   mdlvl->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(mdlvl->fs->DA_CEN, mdlvl->bcp,  &cbcp);  CHKERRQ(ierr);

	// assemble prolongation matrix
	ierr = MatAIJAssemble(P, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void getRowRestrict(
		PetscScalar parent,
		PetscInt    n,
		PetscInt    idx[],
		PetscScalar bc[],
		PetscScalar v[],
		PetscScalar vs[])
{
	PetscInt j, pdof;

	// constrained DOF
	if(parent != DBL_MAX)
	{
		// get parent DOF index
		pdof = (PetscInt)parent;

		// zero out entire row, set parent DOF to unit
		for(j = 0; j < n; j++)
		{
			if(idx[j] == pdof) v[j] = 1.0;
			else               v[j] = 0.0;
		}

	}
	// active DOF
	else
	{
		// set stencil coefficients, zero out constrained DOF
		for(j = 0; j < n; j++)
		{
			if(bc[j] != DBL_MAX) v[j] = 0.0;
			else                 v[j] = vs[j];
		}
	}
}
//---------------------------------------------------------------------------
void getRowProlong(
		PetscInt    parent,
		PetscScalar parent_bc,
		PetscInt    n,
		PetscScalar bc[],
		PetscScalar v[],
		PetscScalar vs[])
{
	PetscInt    j;
	PetscScalar pdof;

	// constrained DOF
	if(parent_bc != DBL_MAX)
	{
		// get parent DOF index
		pdof = (PetscScalar)parent;

		// zero out entire row, set parent DOF to unit
		for(j = 0; j < n; j++)
		{
			if(bc[j] == pdof) v[j] = 1.0;
			else              v[j] = 0.0;
		}
	}
	// active DOF
	else
	{
		// set stencil coefficients, zero out constrained DOF
		for(j = 0; j < n; j++)
		{
			if(bc[j] != DBL_MAX) v[j] = 0.0;
			else                 v[j] = vs[j];
		}
	}
}
//---------------------------------------------------------------------------
// MG interpolation context functions
//---------------------------------------------------------------------------
PetscErrorCode MGInterpCreate(MGInterp *mgi, MatData *coarse, MatData *fine)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set evaluation context
	mgi->coarse = coarse;
	mgi->fine   = fine;

	// create work vectors
	ierr = VecCreateMPI(PETSC_COMM_WORLD, coarse->fs->dof.ln, PETSC_DETERMINE, &mgi->wc); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, fine->fs->dof.ln,   PETSC_DETERMINE, &mgi->wf); CHKERRQ(ierr);

	ierr = VecSetFromOptions(mgi->wc); CHKERRQ(ierr);
	ierr = VecSetFromOptions(mgi->wf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGInterpDestroy(MGInterp *mgi)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = VecDestroy(&mgi->wc); CHKERRQ(ierr);
	ierr = VecDestroy(&mgi->wf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// MG -functions
//---------------------------------------------------------------------------
PetscErrorCode MGCreate(MG *mg, MatData *md)
{
	PetscInt  i, l;
	MGLevel   *fine;
	char      pc_type[_str_len_];
	PetscBool opt_set;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get preconditioner type
	ierr = PetscOptionsGetString(NULL, NULL, "-gmg_pc_type", pc_type, _str_len_, &opt_set); CHKERRQ(ierr);

	// check whether multigrid is requested
	if(opt_set != PETSC_TRUE || strcmp(pc_type, "mg"))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "-gmg_pc_type option is not defined of specified incorrectly (use -gmg_pc_type mg)");
	}

	// clear object
	ierr = PetscMemzero(mg, sizeof(MG)); CHKERRQ(ierr);

	// check multigrid mesh restrictions & get actual number of levels
	ierr = MGGetNumLevels(mg, md); CHKERRQ(ierr);

	// allocate levels
	ierr = PetscMalloc(sizeof(MGLevel)*(size_t)mg->nlvl, &mg->lvls); CHKERRQ(ierr);

	// create levels
	fine = NULL;

	for(i = 0; i < mg->nlvl; i++)
	{
		ierr = MGLevelCreate(&mg->lvls[i], fine, md); CHKERRQ(ierr);

		fine = &mg->lvls[i];
	}

	// create Galerkin multigrid preconditioner
	ierr = PCCreate(PETSC_COMM_WORLD, &mg->pc);          CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(mg->pc, "gmg_");           CHKERRQ(ierr);
	ierr = PCSetType(mg->pc, PCMG);                      CHKERRQ(ierr);
	ierr = PCMGSetLevels(mg->pc, mg->nlvl, NULL);        CHKERRQ(ierr);
	ierr = PCMGSetType(mg->pc, PC_MG_MULTIPLICATIVE);    CHKERRQ(ierr);
	ierr = PCMGSetGalerkin(mg->pc, PC_MG_GALERKIN_BOTH); CHKERRQ(ierr);
	ierr = PCSetFromOptions(mg->pc);                     CHKERRQ(ierr);

	// attach restriction/prolongation matrices to the preconditioner
	for(i = 1, l = mg->nlvl-1; i < mg->nlvl; i++, l--)
	{
		ierr = PCMGSetRestriction  (mg->pc, l, mg->lvls[i].R); CHKERRQ(ierr);
		ierr = PCMGSetInterpolation(mg->pc, l, mg->lvls[i].P); CHKERRQ(ierr);
	}

	// set coarse solver setup flag
	mg->crs_setup = PETSC_FALSE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGDestroy(MG *mg)
{
	PetscInt  i;
	PetscBool flg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// view preconditioner if required
	ierr = PetscOptionsHasName(NULL, NULL, "-gmg_pc_view", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PCView(mg->pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	for(i = 0; i < mg->nlvl; i++)
	{
		ierr = MGLevelDestroy(&mg->lvls[i]); CHKERRQ(ierr);
	}

	ierr = PetscFree(mg->lvls); CHKERRQ(ierr);

	ierr = PCDestroy(&mg->pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode MGSetupCoarse(MG *mg, Mat A)
{
	KSP        ksp;
	PC         pc;
	Mat        mat;
	MGLevel   *lvl;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// skip already configured solver
	if(mg->crs_setup == PETSC_TRUE)
	{
		PetscFunctionReturn(0);
	}

	// get coarse level
	lvl = mg->lvls + mg->nlvl - 1;

	// set dummy coarse solver
	ierr = PCMGGetCoarseSolve(mg->pc, &ksp); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPPREONLY);      CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);               CHKERRQ(ierr);
	ierr = PCSetType(pc, PCNONE);            CHKERRQ(ierr);

	// force setup operators
	ierr = PCSetOperators(mg->pc, A, A);     CHKERRQ(ierr);
	ierr = PCSetUp(mg->pc);                  CHKERRQ(ierr);

	// set near null space on coarse level matrix
	ierr = KSPGetOperators(ksp, &mat, NULL); CHKERRQ(ierr);
	ierr = MatAIJSetNullSpace(mat, lvl->md); CHKERRQ(ierr);

	// set actual coarse solver options
	ierr = KSPSetOptionsPrefix(ksp, "crs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);           CHKERRQ(ierr);

	// set setup flag
	mg->crs_setup = PETSC_TRUE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGSetup(MG *mg, Mat A)
{
	// Matrices are re-assembled here, just in case
	// they will be made matrix- distance- dependent.
	// Currently they depend only on boundary conditions,
	// so changing boundary condition would also require re-assembly.

	MGLevel *lvl, *fine;
	PetscInt i;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	for(i = 1; i < mg->nlvl; i++)
	{
		lvl =  &mg->lvls[i];
		fine = &mg->lvls[i-1];

		ierr = MatDataRestrict     (lvl->md, fine->md); CHKERRQ(ierr);
		ierr = MGLevelSetupRestrict(lvl,     fine);     CHKERRQ(ierr);
		ierr = MGLevelSetupProlong (lvl,     fine);     CHKERRQ(ierr);

		if(i == 2)
		{
			ierr = TestInterp(lvl->md, fine->md, lvl->R, lvl->P);  CHKERRQ(ierr);
		}


	}

	// setup coarse grid solver if necessary
	ierr = MGSetupCoarse(mg, A); CHKERRQ(ierr);

	// tell to recompute preconditioner
	ierr = PCSetOperators(mg->pc, A, A); CHKERRQ(ierr);

	// force setup operators
	ierr = PCSetUp(mg->pc); CHKERRQ(ierr);

	// store matrices in the file if requested
	ierr = MGDumpMat(mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);



	/*

	https://petsc.org/release/src/dm/impls/stag/tutorials/ex4.c.html
	  PetscCall(PCSetType(pc_faces, PCMG));
258:     PetscCall(PCMGSetLevels(pc_faces, ctx->n_levels, NULL));
259:     for (PetscInt level = 0; level < ctx->n_levels; ++level) {
260:       KSP ksp_level;
261:       PC  pc_level;

263:       // Smoothers
264:       PetscCall(PCMGGetSmoother(pc_faces, level, &ksp_level));
265:       PetscCall(KSPGetPC(ksp_level, &pc_level));
266:       PetscCall(KSPSetOperators(ksp_level, A_faces[level], A_faces[level]));
267:       if (level > 0) PetscCall(PCSetType(pc_level, PCJACOBI));

269:       // Transfer Operators
270:       if (level > 0) {
271:         Mat restriction, interpolation;
272:         DM  dm_level   = ctx->levels[level]->dm_faces;
273:         DM  dm_coarser = ctx->levels[level - 1]->dm_faces;

275:         PetscCall(DMCreateInterpolation(dm_coarser, dm_level, &interpolation, NULL));
276:         PetscCall(PCMGSetInterpolation(pc_faces, level, interpolation));
277:         PetscCall(MatDestroy(&interpolation));
278:         PetscCall(DMCreateRestriction(dm_coarser, dm_level, &restriction));
279:         PetscCall(PCMGSetRestriction(pc_faces, level, restriction));
280:         PetscCall(MatDestroy(&restriction));
281:       }
282:     }
283:   }
	 */
}
//---------------------------------------------------------------------------
PetscErrorCode MGApply(PC pc, Vec x, Vec y)
{
	MG *mg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PCShellGetContext(pc, (void**)&mg); CHKERRQ(ierr);

	// apply multigrid preconditioner
	ierr = PCApply(mg->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGDumpMat(MG *mg)
{
	Mat         A;
	KSP         ksp;
	PetscBool   flg;
	PetscInt    l;
	PetscViewer viewer;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// view multigrid matrices if required
	ierr = PetscOptionsHasName(NULL, NULL, "-gmg_dump", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Dumping multigrid matrices to MATLAB\n"); CHKERRQ(ierr);

		viewer = PETSC_VIEWER_BINARY_(PETSC_COMM_WORLD);

		//===================================
		// OUTPUT IN THE ORDER FINE -> COARSE
		//===================================

		for(l = mg->nlvl-1; l >= 0; l--)
		{
			// level matrix
			ierr = PCMGGetSmoother(mg->pc, l, &ksp); CHKERRQ(ierr);
			ierr = KSPGetOperators(ksp, &A, NULL);   CHKERRQ(ierr);
			ierr = MatView(A, viewer);               CHKERRQ(ierr);

			if(l != 0)
			{
				// restriction
				ierr = PCMGGetRestriction(mg->pc, l, &A); CHKERRQ(ierr);
				ierr = MatView(A, viewer);                CHKERRQ(ierr);

				// prolongation
				ierr = PCMGGetInterpolation(mg->pc, l, &A); CHKERRQ(ierr);
				ierr = MatView(A, viewer);                  CHKERRQ(ierr);
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MGGetNumLevels(MG *mg, MatData *md)
{
	// check multigrid mesh restrictions, get actual number of coarsening steps

	FDSTAG   *fs;
	PetscBool opt_set;
	PetscInt  nx, ny, nz, Nx, Ny, Nz, ncors, nlevels;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// check discretization in all directions
	ierr = Discret1DCheckMG(&fs->dsx, "x", &nx); CHKERRQ(ierr);
	ierr = Discret1DCheckMG(&fs->dsy, "y", &ny); CHKERRQ(ierr);
	ierr = Discret1DCheckMG(&fs->dsz, "z", &nz); CHKERRQ(ierr);

	ncors = nx;
	if(ny < ncors) ncors = ny;
	if(nz < ncors) ncors = nz;

	// check number of levels requested on the command line
	ierr = PetscOptionsGetInt(NULL, NULL, "-gmg_pc_mg_levels", &nlevels, &opt_set); CHKERRQ(ierr);

	if(opt_set != PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Number of multigrid levels is not specified. Use option -gmg_pc_mg_levels. Max # of levels: %lld", (LLD)(ncors+1));
	}
	else if(nlevels < 2 || nlevels > ncors+1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect # of multigrid levels specified. Requested: %lld. Max. possible: %lld", (LLD)nlevels, (LLD)(ncors+1));
	}

	// set actual number of coarsening steps
	ncors = nlevels-1;

	// print grid statistics
	nx = fs->dsx.ncels >> ncors;
	ny = fs->dsy.ncels >> ncors;
	nz = fs->dsz.ncels >> ncors;

	Nx = nx*fs->dsx.nproc;
	Ny = ny*fs->dsy.nproc;
	Nz = nz*fs->dsz.nproc;
	ierr = PetscPrintf(PETSC_COMM_WORLD, "   Global coarse grid [nx,ny,nz] : [%lld, %lld, %lld]\n", (LLD)Nx, (LLD)Ny, (LLD)Nz); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "   Local coarse grid  [nx,ny,nz] : [%lld, %lld, %lld]\n", (LLD)nx, (LLD)ny, (LLD)nz); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "   Number of multigrid levels    :  %lld\n", (LLD)nlevels);                            CHKERRQ(ierr);

	// store number of levels
	mg->nlvl = nlevels;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

// test functions

PetscErrorCode comareVecs(Vec va, Vec vb)
{
	Vec       diff;
	PetscReal nrm;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = VecDuplicate(va, &diff); CHKERRQ(ierr);

	ierr = VecWAXPY(diff, -1.0, va, vb);

	ierr = VecNorm(diff, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "   *** \n");
	PetscPrintf(PETSC_COMM_WORLD, "   *** \n");
	PetscPrintf(PETSC_COMM_WORLD, "   *** \n");

	PetscPrintf(PETSC_COMM_WORLD, "   Difference            :  %g\n", nrm);

	ierr = VecDestroy(&diff); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

PetscErrorCode genRandVec(MatData *md, Vec *v)
{
	PetscRandom  rctx;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);

	ierr = VecCreateMPI(PETSC_COMM_WORLD, md->fs->dof.ln, PETSC_DETERMINE, v); CHKERRQ(ierr);

	ierr = VecSetRandom((*v), rctx); CHKERRQ(ierr);

	ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode VecSetBC(MatData *md, Vec v)
{
	PetscInt     i, num, *list;
	PetscScalar  *va;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = VecGetArray(v, &va); CHKERRQ(ierr);

	// zero out constrained residuals (velocity)
	num   = md->vNumSPC;
	list  = md->vSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = 0.0;

	// zero out constrained residuals (pressure)
	num   = md->pNumSPC;
	list  = md->pSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = 0.0;

	ierr = VecRestoreArray(v, &va); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode TestInterp(MatData *coarse, MatData *fine, Mat R, Mat P)
{
	Vec zc, zf;
	Vec rc, rf;
	Vec wc, wf;
	Vec wcmf, wfmf;
	Mat RMF, PMF;
	PetscInt nc, nf;
	MGInterp mgi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	nc = coarse->fs->dof.ln;
	nf = fine->fs->dof.ln;

	ierr = MGInterpCreate(&mgi, coarse, fine); CHKERRQ(ierr);

	// create restriction operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, nc, nf, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &RMF); CHKERRQ(ierr);
	ierr = MatSetUp(RMF); CHKERRQ(ierr);

	ierr = MatShellSetOperation(RMF, MATOP_MULT_ADD, (void(*)(void))MatFreeApplyRestrict); CHKERRQ(ierr);
	ierr = MatShellSetContext(RMF, (void*)&mgi); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(RMF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (RMF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// create prolongation operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, nf, nc, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &PMF); CHKERRQ(ierr);
	ierr = MatSetUp(PMF); CHKERRQ(ierr);

	ierr = MatShellSetOperation(PMF, MATOP_MULT_ADD, (void(*)(void))MatFreeApplyProlong); CHKERRQ(ierr);
	ierr = MatShellSetContext(PMF, (void*)&mgi); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(PMF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (PMF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// zero
	ierr = VecDuplicate(mgi.wc, &zc); CHKERRQ(ierr);
	ierr = VecDuplicate(mgi.wf, &zf); CHKERRQ(ierr);

	ierr = VecZeroEntries(zc); CHKERRQ(ierr);
	ierr = VecZeroEntries(zf); CHKERRQ(ierr);

	// random
	ierr = VecDuplicate(mgi.wc, &rc); CHKERRQ(ierr);
	ierr = VecDuplicate(mgi.wf, &rf); CHKERRQ(ierr);

	ierr = genRandVec(coarse, &rc); CHKERRQ(ierr);
	ierr = genRandVec(fine  , &rf); CHKERRQ(ierr);

	ierr = VecSetBC(coarse, rc); CHKERRQ(ierr);
	ierr = VecSetBC(fine,   rf); CHKERRQ(ierr);

	// work
	ierr = VecDuplicate(mgi.wc, &wc); CHKERRQ(ierr);
	ierr = VecDuplicate(mgi.wf, &wf); CHKERRQ(ierr);

	ierr = VecDuplicate(mgi.wc, &wcmf); CHKERRQ(ierr);
	ierr = VecDuplicate(mgi.wf, &wfmf); CHKERRQ(ierr);

	// restriction
	ierr = MatMult   (R,   rf,     wc);   CHKERRQ(ierr);
	ierr = MatMultAdd(RMF, rf, zc, wcmf); CHKERRQ(ierr);

	// prolongation
	ierr = MatMult   (P,   rc,     wf);   CHKERRQ(ierr);
	ierr = MatMultAdd(PMF, rc, zf, wfmf); CHKERRQ(ierr);

	// comparison
	ierr = comareVecs(wc, wcmf); CHKERRQ(ierr);
	ierr = comareVecs(wf, wfmf); CHKERRQ(ierr);

	ierr = MGInterpDestroy(&mgi); CHKERRQ(ierr);

	ierr = MatDestroy(&RMF); CHKERRQ(ierr);
	ierr = MatDestroy(&PMF); CHKERRQ(ierr);

	ierr = VecDestroy(&zc);   CHKERRQ(ierr);
	ierr = VecDestroy(&zf);   CHKERRQ(ierr);
	ierr = VecDestroy(&rc);   CHKERRQ(ierr);
	ierr = VecDestroy(&rf);   CHKERRQ(ierr);
	ierr = VecDestroy(&wc);   CHKERRQ(ierr);
	ierr = VecDestroy(&wf);   CHKERRQ(ierr);
	ierr = VecDestroy(&wcmf); CHKERRQ(ierr);
	ierr = VecDestroy(&wfmf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode TestInterpBC(MatData *coarse, MatData *fine, Mat R, Mat P)
{
	MGInterp mgi;

	Vec rc, rcbc;
	Vec wf, wfbc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MGInterpCreate(&mgi, coarse, fine); CHKERRQ(ierr);

	ierr = genRandVec(coarse, &rc); CHKERRQ(ierr);

	ierr = VecDuplicate(rc, &rcbc); CHKERRQ(ierr);

	ierr = VecCopy(rc, rcbc); CHKERRQ(ierr);

	ierr = VecSetBC(coarse, rcbc); CHKERRQ(ierr);

	ierr = VecDuplicate(mgi.wf, &wf); CHKERRQ(ierr);
	ierr = VecDuplicate(mgi.wf, &wfbc); CHKERRQ(ierr);


	ierr = MatMult(P, rc,   wf);   CHKERRQ(ierr);
	ierr = MatMult(P, rcbc, wfbc); CHKERRQ(ierr);

	ierr = VecSetBC(fine, wf);   CHKERRQ(ierr);
	ierr = VecSetBC(fine, wfbc); CHKERRQ(ierr);

	ierr = comareVecs(wf, wfbc); CHKERRQ(ierr);

	ierr = VecDestroy(&rc);   CHKERRQ(ierr);
	ierr = VecDestroy(&rcbc); CHKERRQ(ierr);
	ierr = VecDestroy(&wf);    CHKERRQ(ierr);
	ierr = VecDestroy(&wfbc); CHKERRQ(ierr);

	ierr = MGInterpDestroy(&mgi); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

