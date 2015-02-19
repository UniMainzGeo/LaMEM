//---------------------------------------------------------------------------
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "matrix.h"
#include "multigrid.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * remove hierarchy of grids & bc-objects (use info from fine level)
// * preallocate all restriction & interpolation operators
// * implement fine-level preconditionier completely matrix-free
// * coordinate- viscosity- residual-dependent restriction & interpolation
//---------------------------------------------------------------------------
// MG -functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelCreate"
PetscErrorCode MGLevelCreate(MGLevel *lvl, MGLevel *fine, FDSTAG *fs, BCCtx *bc)
{
	PetscInt         i, ln, lnfine;
	PetscInt         Nx,   Ny,   Nz;
	PetscInt         Px,   Py,   Pz;
	const PetscInt  *plx, *ply, *plz;
	PetscInt        *lx,  *ly,  *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!fine)
	{
		// copy data from the staggered grid on the finest level
		lvl->DA_CEN = fs->DA_CEN;
		lvl->DA_X   = fs->DA_X;
		lvl->DA_Y   = fs->DA_Y;
		lvl->DA_Z   = fs->DA_Z;
		lvl->dof    = fs->dof;
		lvl->bcvx   = bc->bcvx;
		lvl->bcvy   = bc->bcvy;
		lvl->bcvz   = bc->bcvz;
		lvl->bcp    = bc->bcp;
		lvl->R      = NULL;
		lvl->P      = NULL;
	}
	else
	{
		// get number of cells & processors in the fine grid
		ierr = DMDAGetInfo(fine->DA_CEN, 0, &Nx, &Ny, &Nz, &Px, &Py, &Pz, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

		// get number of cells per processor in fine grid
		ierr = DMDAGetOwnershipRanges(fine->DA_CEN, &plx, &ply, &plz); CHKERRQ(ierr);

		ierr = makeIntArray(&lx, plx, Px); CHKERRQ(ierr);
		ierr = makeIntArray(&ly, ply, Py); CHKERRQ(ierr);
		ierr = makeIntArray(&lz, plz, Pz); CHKERRQ(ierr);

		// coarsen uniformly in every direction
		Nx /= 2;  for(i = 0; i < Px; i++) lx[i] /= 2;
		Ny /= 2;  for(i = 0; i < Py; i++) ly[i] /= 2;
		Nz /= 2;  for(i = 0; i < Pz; i++) lz[i] /= 2;

		// central points (DA_CEN) with boundary ghost points (1-layer stencil box)
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_CEN); CHKERRQ(ierr);

		// X face (DA_X) with boundary ghost points (1-layer stencil box)
		lx[Px-1]++;
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx+1, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_X); CHKERRQ(ierr);
		lx[Px-1]--;

		// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
		ly[Py-1]++;
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny+1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_Y); CHKERRQ(ierr);
		ly[Py-1]--;

		// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
		lz[Pz-1]++;
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny, Nz+1, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_Z); CHKERRQ(ierr);

		// clear temporary storage
		ierr = PetscFree(lx); CHKERRQ(ierr);
		ierr = PetscFree(ly); CHKERRQ(ierr);
		ierr = PetscFree(lz); CHKERRQ(ierr);

		// create index arrays
		ierr = DOFIndexCreate(&lvl->dof, lvl->DA_CEN, lvl->DA_X, lvl->DA_Y, lvl->DA_Z); CHKERRQ(ierr);

		// create restricted boundary condition vectors
		ierr = DMCreateLocalVector(fs->DA_X,   &lvl->bcvx);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_Y,   &lvl->bcvy);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_Z,   &lvl->bcvz);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_CEN, &lvl->bcp);   CHKERRQ(ierr);

		// compute index arrays
		ierr = DOFIndexCompute(&lvl->dof, fine->dof.idxmod); CHKERRQ(ierr);

		// get matrix sizes
		if     (lvl->dof.idxmod == IDXCOUPLED)   { ln = lvl->dof.ln;  lnfine = fine->dof.ln;  }
		else if(lvl->dof.idxmod == IDXUNCOUPLED) { ln = lvl->dof.lnv; lnfine = fine->dof.lnv; }

		// preallocate restriction & prolongation matrices
		// WARNING! CONSTANT SIZE PREALLOCATION
		// ADD VARIABLE PREALLOCATION TO THESE MATRICES
		ierr = MatAIJCreate(ln,     lnfine, 12, NULL, 4, NULL, &lvl->R); CHKERRQ(ierr);
		ierr = MatAIJCreate(lnfine, ln,     8,  NULL, 7, NULL, &lvl->P); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelDestroy"
PetscErrorCode MGLevelDestroy(MGLevel *lvl)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(lvl->R)
	{
		ierr = DMDestroy(&lvl->DA_CEN);    CHKERRQ(ierr);
		ierr = DMDestroy(&lvl->DA_X);      CHKERRQ(ierr);
		ierr = DMDestroy(&lvl->DA_Y);      CHKERRQ(ierr);
		ierr = DMDestroy(&lvl->DA_Z);      CHKERRQ(ierr);
		ierr = DOFIndexDestroy(&lvl->dof); CHKERRQ(ierr);
		ierr = VecDestroy(&lvl->bcvx);     CHKERRQ(ierr);
		ierr = VecDestroy(&lvl->bcvy);     CHKERRQ(ierr);
		ierr = VecDestroy(&lvl->bcvz);     CHKERRQ(ierr);
		ierr = VecDestroy(&lvl->bcp);      CHKERRQ(ierr);
		ierr = MatDestroy(&lvl->R);        CHKERRQ(ierr);
		ierr = MatDestroy(&lvl->P);        CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelRestrictBC"
PetscErrorCode MGLevelRestrictBC(MGLevel *lvl, MGLevel *fine)
{
	// restrict boundary condition vectors from fine to coarse level

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// mark all variables unconstrained
	ierr = VecSet(lvl->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcp,  DBL_MAX); CHKERRQ(ierr);



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelSetupRestrict"
PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine)
{
	Mat         R;
	PetscScalar v[12], vs[12], bc[12];
	PetscInt    idx[12];
	PetscInt    mx, my, mz;
	PetscInt    row, I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fivx,  ***fivy,  ***fivz,  ***fip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	R = lvl->R;

	// clear restriction matrix coefficients
	ierr = MatZeroEntries(R); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->dof.ivx, &fivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->dof.ivy, &fivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->dof.ivz, &fivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->dof.ip,  &fip);  CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->dof.ivx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->dof.ivy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->dof.ivz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->dof.ip,  &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,   lvl->dof.ivx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,   lvl->dof.ivy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,   lvl->dof.ivz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN, lvl->dof.ip,  &cbcp);  CHKERRQ(ierr);

	// get total number of nodes in the fine grid
	ierr = DMDAGetInfo(fine->DA_CEN, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	mx++;
	my++;
	mz++;

	// get global index of the first row in coarse grid
	if     (lvl->dof.idxmod == IDXCOUPLED)   { row = lvl->dof.st;  }
	else if(lvl->dof.idxmod == IDXUNCOUPLED) { row = lvl->dof.stv; }

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
	ierr = DMDAGetCorners(lvl->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)fivx[K  ][J  ][I-1];
		idx[1]  = (PetscInt)fivx[K  ][J+1][I-1];
		idx[2]  = (PetscInt)fivx[K+1][J  ][I-1];
		idx[3]  = (PetscInt)fivx[K+1][J+1][I-1];
		idx[4]  = (PetscInt)fivx[K  ][J  ][I  ];
		idx[5]  = (PetscInt)fivx[K  ][J+1][I  ];
		idx[6]  = (PetscInt)fivx[K+1][J  ][I  ];
		idx[7]  = (PetscInt)fivx[K+1][J+1][I  ];
		idx[8]  = (PetscInt)fivx[K  ][J  ][I+1];
		idx[9]  = (PetscInt)fivx[K  ][J+1][I+1];
		idx[10] = (PetscInt)fivx[K+1][J  ][I+1];
		idx[11] = (PetscInt)fivx[K+1][J+1][I+1];

		// get fine grid boundary conditions
		bc[0]  = (PetscInt)fbcvx[K  ][J  ][I-1];
		bc[1]  = (PetscInt)fbcvx[K  ][J+1][I-1];
		bc[2]  = (PetscInt)fbcvx[K+1][J  ][I-1];
		bc[3]  = (PetscInt)fbcvx[K+1][J+1][I-1];
		bc[4]  = (PetscInt)fbcvx[K  ][J  ][I  ];
		bc[5]  = (PetscInt)fbcvx[K  ][J+1][I  ];
		bc[6]  = (PetscInt)fbcvx[K+1][J  ][I  ];
		bc[7]  = (PetscInt)fbcvx[K+1][J+1][I  ];
		bc[8]  = (PetscInt)fbcvx[K  ][J  ][I+1];
		bc[9]  = (PetscInt)fbcvx[K  ][J+1][I+1];
		bc[10] = (PetscInt)fbcvx[K+1][J  ][I+1];
		bc[11] = (PetscInt)fbcvx[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvx[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(lvl->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)fivy[K  ][J-1][I  ];
		idx[1]  = (PetscInt)fivy[K  ][J-1][I+1];
		idx[2]  = (PetscInt)fivy[K+1][J-1][I  ];
		idx[3]  = (PetscInt)fivy[K+1][J-1][I+1];
		idx[4]  = (PetscInt)fivy[K  ][J  ][I  ];
		idx[5]  = (PetscInt)fivy[K  ][J  ][I+1];
		idx[6]  = (PetscInt)fivy[K+1][J  ][I  ];
		idx[7]  = (PetscInt)fivy[K+1][J  ][I+1];
		idx[8]  = (PetscInt)fivy[K  ][J+1][I  ];
		idx[9]  = (PetscInt)fivy[K  ][J+1][I+1];
		idx[10] = (PetscInt)fivy[K+1][J+1][I  ];
    	idx[11] = (PetscInt)fivy[K+1][J+1][I+1];

    	// get fine grid boundary conditions
		bc[0]  = (PetscInt)fbcvy[K  ][J-1][I  ];
		bc[1]  = (PetscInt)fbcvy[K  ][J-1][I+1];
		bc[2]  = (PetscInt)fbcvy[K+1][J-1][I  ];
		bc[3]  = (PetscInt)fbcvy[K+1][J-1][I+1];
		bc[4]  = (PetscInt)fbcvy[K  ][J  ][I  ];
		bc[5]  = (PetscInt)fbcvy[K  ][J  ][I+1];
		bc[6]  = (PetscInt)fbcvy[K+1][J  ][I  ];
		bc[7]  = (PetscInt)fbcvy[K+1][J  ][I+1];
		bc[8]  = (PetscInt)fbcvy[K  ][J+1][I  ];
		bc[9]  = (PetscInt)fbcvy[K  ][J+1][I+1];
		bc[10] = (PetscInt)fbcvy[K+1][J+1][I  ];
    	bc[11] = (PetscInt)fbcvy[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvy[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(lvl->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// get fine grid stencil
		idx[0]  = (PetscInt)fivz[K-1][J  ][I  ];
		idx[1]  = (PetscInt)fivz[K-1][J  ][I+1];
		idx[2]  = (PetscInt)fivz[K-1][J+1][I  ];
		idx[3]  = (PetscInt)fivz[K-1][J+1][I+1];
		idx[4]  = (PetscInt)fivz[K  ][J  ][I  ];
		idx[5]  = (PetscInt)fivz[K  ][J  ][I+1];
		idx[6]  = (PetscInt)fivz[K  ][J+1][I  ];
		idx[7]  = (PetscInt)fivz[K  ][J+1][I+1];
		idx[8]  = (PetscInt)fivz[K+1][J  ][I  ];
		idx[9]  = (PetscInt)fivz[K+1][J  ][I+1];
		idx[10] = (PetscInt)fivz[K+1][J+1][I  ];
    	idx[11] = (PetscInt)fivz[K+1][J+1][I+1];

    	// get fine grid boundary conditions
		bc[0]  = (PetscInt)fbcvz[K-1][J  ][I  ];
		bc[1]  = (PetscInt)fbcvz[K-1][J  ][I+1];
		bc[2]  = (PetscInt)fbcvz[K-1][J+1][I  ];
		bc[3]  = (PetscInt)fbcvz[K-1][J+1][I+1];
		bc[4]  = (PetscInt)fbcvz[K  ][J  ][I  ];
		bc[5]  = (PetscInt)fbcvz[K  ][J  ][I+1];
		bc[6]  = (PetscInt)fbcvz[K  ][J+1][I  ];
		bc[7]  = (PetscInt)fbcvz[K  ][J+1][I+1];
		bc[8]  = (PetscInt)fbcvz[K+1][J  ][I  ];
		bc[9]  = (PetscInt)fbcvz[K+1][J  ][I+1];
		bc[10] = (PetscInt)fbcvz[K+1][J+1][I  ];
    	bc[11] = (PetscInt)fbcvz[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(cbcvz[k][j][i], 12, idx, bc, v, vs);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(lvl->dof.idxmod == IDXCOUPLED)
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
		ierr = DMDAGetCorners(lvl->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get fine grid indices
			I = 2*i;
			J = 2*j;
			K = 2*k;

			// get fine grid stencil
			idx[0] = (PetscInt)fip[K  ][J  ][I  ];
			idx[1] = (PetscInt)fip[K  ][J  ][I+1];
			idx[2] = (PetscInt)fip[K  ][J+1][I  ];
			idx[3] = (PetscInt)fip[K  ][J+1][I+1];
			idx[4] = (PetscInt)fip[K+1][J  ][I  ];
			idx[5] = (PetscInt)fip[K+1][J  ][I+1];
			idx[6] = (PetscInt)fip[K+1][J+1][I  ];
			idx[7] = (PetscInt)fip[K+1][J+1][I+1];

			// get fine grid boundary conditions
			bc[0] = (PetscInt)fbcp[K  ][J  ][I  ];
			bc[1] = (PetscInt)fbcp[K  ][J  ][I+1];
			bc[2] = (PetscInt)fbcp[K  ][J+1][I  ];
			bc[3] = (PetscInt)fbcp[K  ][J+1][I+1];
			bc[4] = (PetscInt)fbcp[K+1][J  ][I  ];
			bc[5] = (PetscInt)fbcp[K+1][J  ][I+1];
			bc[6] = (PetscInt)fbcp[K+1][J+1][I  ];
			bc[7] = (PetscInt)fbcp[K+1][J+1][I+1];

			// setup row of restriction matrix
			getRowRestrict(cbcp[k][j][i], 8, idx, bc, v, vs);

			// store full matrix row
			ierr = MatSetValues(R, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->dof.ivx, &fivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->dof.ivy, &fivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->dof.ivz, &fivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->dof.ip,  &fip);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->dof.ivx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->dof.ivy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->dof.ivz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->dof.ip,  &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->dof.ivx,  &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->dof.ivy,  &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->dof.ivz,  &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->dof.ip,   &cbcp);  CHKERRQ(ierr);

	// assemble restriction matrix
	ierr = MatAIJAssemble(R, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelSetupProlong"
PetscErrorCode MGLevelSetupProlong(MGLevel *lvl, MGLevel *fine)
{
	Mat P;
	PetscScalar v[8];
	PetscInt    idx[8];
	PetscInt    mx, my, mz;
	PetscInt    row, I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	P = lvl->P;

	// clear prolongation matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);

	// access index vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,   lvl->dof.ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,   lvl->dof.ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,   lvl->dof.ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN, lvl->dof.ip,   &ip);  CHKERRQ(ierr);

	// get total number of nodes in the coarse grid
	ierr = DMDAGetInfo(lvl->DA_CEN, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	mx++;
	my++;
	mz++;

	// get global index of the first row in the fine grid
	if     (fine->dof.idxmod == IDXCOUPLED)   { row = fine->dof.st;  }
	else if(fine->dof.idxmod == IDXUNCOUPLED) { row = fine->dof.stv; }

	// set velocity weights
	v[0] = 9.0/32.0;
	v[1] = 3.0/32.0;
	v[2] = 3.0/32.0;
	v[3] = 1.0/32.0;
	v[4] = 9.0/32.0;
	v[5] = 3.0/32.0;
	v[6] = 3.0/32.0;
	v[7] = 1.0/32.0;

	//---------------------
	// X-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(fine->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I;

		if(j % 2) J1 = J+1; else J1 = J-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(J1 == -1 || J1 == my) J1 = J;
		if(K1 == -1 || K1 == mz) K1 = K;

		idx[0] = (PetscInt)ivx[K ][J ][I ];
		idx[1] = (PetscInt)ivx[K ][J1][I ];
		idx[2] = (PetscInt)ivx[K1][J ][I ];
		idx[3] = (PetscInt)ivx[K1][J1][I ];
		idx[4] = (PetscInt)ivx[K ][J ][I1];
		idx[5] = (PetscInt)ivx[K ][J1][I1];
		idx[6] = (PetscInt)ivx[K1][J ][I1];
		idx[7] = (PetscInt)ivx[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Y-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(fine->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(j % 2) J1 = J+1; else J1 = J;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(I1 == -1 || I1 == mx) I1 = I;
		if(K1 == -1 || K1 == mz) K1 = K;

		idx[0] = (PetscInt)ivy[K ][J ][I ];
		idx[1] = (PetscInt)ivy[K ][J ][I1];
		idx[2] = (PetscInt)ivy[K1][J ][I ];
		idx[3] = (PetscInt)ivy[K1][J ][I1];
		idx[4] = (PetscInt)ivy[K ][J1][I ];
		idx[5] = (PetscInt)ivy[K ][J1][I1];
		idx[6] = (PetscInt)ivy[K1][J1][I ];
		idx[7] = (PetscInt)ivy[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Z-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(fine->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(k % 2) K1 = K+1; else K1 = K;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(j % 2) J1 = J+1; else J1 = J-1;

		if(I1 == -1 || I1 == mx) I1 = I;
		if(J1 == -1 || J1 == my) J1 = J;

		idx[0] = (PetscInt)ivz[K ][J ][I ];
		idx[1] = (PetscInt)ivz[K ][J ][I1];
		idx[2] = (PetscInt)ivz[K ][J1][I ];
		idx[3] = (PetscInt)ivz[K ][J1][I1];
		idx[4] = (PetscInt)ivz[K1][J ][I ];
		idx[5] = (PetscInt)ivz[K1][J ][I1];
		idx[6] = (PetscInt)ivz[K1][J1][I ];
		idx[7] = (PetscInt)ivz[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(fine->dof.idxmod == IDXCOUPLED)
	{
		// set pressure weights (direct injection)
		v[0] = 1.0;

		//---------------------
		// P-points (fine grid)
		//---------------------
		ierr = DMDAGetCorners(fine->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get coarse grid indices
			I = i/2;
			J = j/2;
			K = k/2;

			idx[0] = (PetscInt)ip[K][J][I];

			// store full matrix row
			ierr = MatSetValues(P, 1, &row, 1, idx, v, ADD_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP

	}
	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_X,   lvl->dof.ivx, &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,   lvl->dof.ivy, &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,   lvl->dof.ivz, &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN, lvl->dof.ip,  &ip);  CHKERRQ(ierr);

	// assemble prolongation matrix
	ierr = MatAIJAssemble(P, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// MG -functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCreate"
PetscErrorCode MGCreate(MG *mg, FDSTAG *fs, BCCtx *bc)
{
	KSP       ksp;
	PetscInt  i, l;
	MGLevel   *fine;
	char      pc_type[MAX_NAME_LEN];
	PetscBool opt_set;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get preconditioner type
	ierr = PetscOptionsGetString(NULL, "-gmg_pc_type", pc_type, MAX_NAME_LEN, &opt_set); CHKERRQ(ierr);

	// check whether multigrid is requested
	if(opt_set != PETSC_TRUE || strcmp(pc_type, "mg"))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "-gmg_pc_type option is not defined of specified incorrectly (use -gmg_pc_type mg)");
	}

	// clear object
	ierr = PetscMemzero(mg, sizeof(MG)); CHKERRQ(ierr);

	// store fine grid contexts
	mg->fs = fs;
	mg->bc = bc;

	// check multigrid mesh restrictions & get actual number of levels
	ierr = MGGetNumLevels(mg); CHKERRQ(ierr);

	// allocate levels
	ierr = PetscMalloc(sizeof(MGLevel)*(size_t)mg->nlvl, &mg->lvls); CHKERRQ(ierr);

	// create levels
	fine = NULL;

	for(i = 0; i < mg->nlvl; i++)
	{
		ierr = MGLevelCreate(&mg->lvls[i], fine, fs, bc); CHKERRQ(ierr);

		fine = &mg->lvls[i];
	}

	// create Galerkin multigrid preconditioner
	ierr = PCCreate(PETSC_COMM_WORLD, &mg->pc);       CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(mg->pc, "gmg_");        CHKERRQ(ierr);
	ierr = PCSetType(mg->pc, PCMG);                   CHKERRQ(ierr);
	ierr = PCMGSetLevels(mg->pc, mg->nlvl, NULL);     CHKERRQ(ierr);
	ierr = PCMGSetType(mg->pc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
	ierr = PCMGSetGalerkin(mg->pc, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = PCSetFromOptions(mg->pc);                  CHKERRQ(ierr);

	// setup coarse solver
	ierr = PCMGGetCoarseSolve(mg->pc, &ksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp, "crs_");  CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);            CHKERRQ(ierr);

	// attach restriction/prolongation matrices to the preconditioner
	for(i = 1, l = mg->nlvl-1; i < mg->nlvl; i++, l--)
	{
		ierr = PCMGSetRestriction  (mg->pc, l, mg->lvls[i].R); CHKERRQ(ierr);
		ierr = PCMGSetInterpolation(mg->pc, l, mg->lvls[i].P); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGDestroy"
PetscErrorCode MGDestroy(MG *mg)
{

	PetscInt  i;
	PetscBool flg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// view preconditioner if required
	ierr = PetscOptionsHasName(NULL, "-gmg_pc_view", &flg); CHKERRQ(ierr);

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
#undef __FUNCT__
#define __FUNCT__ "MGSetup"
PetscErrorCode MGSetup(MG *mg, Mat A)
{
	// Matrices are re-assembled here, just in case
	// they will be made matrix- distance- dependent.
	// Currently they depend only on boundary conditions,
	// so changing boundary condition would also require re-assembly.

	PetscInt  i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	for(i = 1; i < mg->nlvl; i++)
	{
		ierr = MGLevelRestrictBC   (&mg->lvls[i], &mg->lvls[i-1]); CHKERRQ(ierr);
		ierr = MGLevelSetupRestrict(&mg->lvls[i], &mg->lvls[i-1]); CHKERRQ(ierr);
		ierr = MGLevelSetupProlong (&mg->lvls[i], &mg->lvls[i-1]); CHKERRQ(ierr);
	}

	// tell to recompute preconditioner
	ierr = PCSetOperators(mg->pc, A, A); CHKERRQ(ierr);

	// force setup operators
	ierr = PCSetUp(mg->pc); CHKERRQ(ierr);

	// store matrices in the file if requested
	ierr = MGDumpMat(mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGApply"
PetscErrorCode MGApply(PC pc, Vec x, Vec y)
{
	MG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PCShellGetContext(pc, (void**)&mg); CHKERRQ(ierr);

	// apply multigrid preconditioner
	ierr = PCApply(mg->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGDumpMat"
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
	ierr = PetscOptionsHasName(NULL, "-gmg_dump", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Dumping multigrid matrices to MATLAB\n"); CHKERRQ(ierr);

		viewer = PETSC_VIEWER_BINARY_(PetscObjectComm((PetscObject)mg->pc));

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
#undef __FUNCT__
#define __FUNCT__ "MGGetNumLevels"
PetscErrorCode MGGetNumLevels(MG *mg)
{
	// check multigrid mesh restrictions, get actual number of coarsening steps

	FDSTAG   *fs;
	PetscBool opt_set;
	PetscInt  nx, ny, nz, Nx, Ny, Nz, ncors, nlevels;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = mg->fs;

	// check discretization in all directions
	ierr = Discret1DCheckMG(&fs->dsx, "x", &nx); CHKERRQ(ierr);                ncors = nx;
	ierr = Discret1DCheckMG(&fs->dsy, "y", &ny); CHKERRQ(ierr); if(ny < ncors) ncors = ny;
	ierr = Discret1DCheckMG(&fs->dsz, "z", &nz); CHKERRQ(ierr); if(nz < ncors) ncors = nz;

	// check number of levels requested on the command line
	ierr = PetscOptionsGetInt(NULL, "-gmg_pc_mg_levels", &nlevels, &opt_set); CHKERRQ(ierr);

	if(opt_set != PETSC_TRUE)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Number of multigrid levels is not specified. Use option -gmg_pc_mg_levels. Max # of levels: %lld", (LLD)(ncors+1));
	}
	else if(nlevels < 2 || nlevels > ncors+1)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect # of multigrid levels specified. Requested: %lld. Max. possible: %lld", (LLD)nlevels, (LLD)(ncors+1));
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

	ierr = PetscPrintf(PETSC_COMM_WORLD, " Total coarse grid size    [nx, ny, nz] : [%lld, %lld, %lld]\n", (LLD)Nx, (LLD)Ny, (LLD)Nz); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, " Coarse grid per processor [nx, ny, nz] : [%lld, %lld, %lld]\n", (LLD)nx, (LLD)ny, (LLD)nz); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of multigrid levels             : %lld\n", (LLD)nlevels);                            CHKERRQ(ierr);

	// store number of levels
	mg->nlvl = nlevels;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void getRowRestrict(PetscScalar parent, PetscInt n, PetscInt idx[], PetscScalar bc[], PetscScalar v[], PetscScalar vs[])
{
	PetscInt j, jj;

	// constrained DOF
	if(parent != DBL_MAX)
	{
		// get parent DOF index
		jj = (PetscInt)parent;

		// zero out entire row, set parent DOF to unit
		for(j = 0; j < n; j++)
		{
			if(idx[j] == jj) v[j] = 1.0;
			else             v[j] = 0.0;
		}

	}
	// active DOF
	else
	{
		// set stencil coefficient, zero out constrained DOF
		for(j = 0; j < n; j++)
		{
			if(bc[j] != DBL_MAX) v[j] = 0.0;
			else                 v[j] = vs[j];
		}
	}
}
//---------------------------------------------------------------------------

void getRowProlong(PetscScalar parent, PetscInt n, PetscInt idx[], PetscScalar bc[], PetscScalar v[], PetscScalar vs[])
{
	PetscInt j, jj;

	// constrained DOF
	if(parent != DBL_MAX)
	{
		// get parent DOF index
		jj = (PetscInt)parent;

		// zero out entire row, set parent DOF to unit
		for(j = 0; j < n; j++)
		{
			if(idx[j] == jj) v[j] = 1.0;
			else             v[j] = 0.0;
		}

	}
	// active DOF
	else
	{
		// set stencil coefficient, zero out constrained DOF
		for(j = 0; j < n; j++)
		{
			if(bc[j] != DBL_MAX) v[j] = 0.0;
			else                 v[j] = vs[j];
		}
	}
}


