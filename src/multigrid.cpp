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
 **    filename:   multigrid.c
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
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "multigrid.h"
#include "matrix.h"
#include "JacRes.h"
#include "bc.h"
#include "tools.h"
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
	PetscInt         i, ln=0, lnfine=0, refine_y;
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

		// get refinement factor in y-direction of central array (in 2D, don't refine in y-direction)
		ierr = DMDAGetRefinementFactor(fine->DA_CEN, NULL, &refine_y,NULL); CHKERRQ(ierr);

		// get number of cells per processor in fine grid
		ierr = DMDAGetOwnershipRanges(fine->DA_CEN, &plx, &ply, &plz); CHKERRQ(ierr);

		ierr = makeIntArray(&lx, plx, Px); CHKERRQ(ierr);
		ierr = makeIntArray(&ly, ply, Py); CHKERRQ(ierr);
		ierr = makeIntArray(&lz, plz, Pz); CHKERRQ(ierr);

		// coarsen uniformly in every direction
		Nx /= 2;  for(i = 0; i < Px; i++) lx[i] /= 2;
		if (refine_y==1){
			Ny /= 1;  for(i = 0; i < Py; i++) ly[i] /= 1;		// don't refine in y-direction
		}
		else
		{
			Ny /= 2;  for(i = 0; i < Py; i++) ly[i] /= 2;
		}

		Nz /= 2;  for(i = 0; i < Pz; i++) lz[i] /= 2;

		// central points (DA_CEN) with boundary ghost points (1-layer stencil box)
		ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_CEN); CHKERRQ(ierr);

		// X face (DA_X) with boundary ghost points (1-layer stencil box)
		lx[Px-1]++;
		ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx+1, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_X); CHKERRQ(ierr);
		lx[Px-1]--;

		// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
		ly[Py-1]++;
		ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny+1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_Y); CHKERRQ(ierr);
		ly[Py-1]--;

		// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
		lz[Pz-1]++;
		ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
			Nx, Ny, Nz+1, Px, Py, Pz, 1, 1, lx, ly, lz, &lvl->DA_Z); CHKERRQ(ierr);

		// clear temporary storage
		ierr = PetscFree(lx); CHKERRQ(ierr);
		ierr = PetscFree(ly); CHKERRQ(ierr);
		ierr = PetscFree(lz); CHKERRQ(ierr);

		// create index arrays
		ierr = DOFIndexCreate(&lvl->dof, lvl->DA_CEN, lvl->DA_X, lvl->DA_Y, lvl->DA_Z); CHKERRQ(ierr);

		// create restricted boundary condition vectors
		ierr = DMCreateLocalVector(lvl->DA_X,   &lvl->bcvx);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->DA_Y,   &lvl->bcvy);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->DA_Z,   &lvl->bcvz);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->DA_CEN, &lvl->bcp);   CHKERRQ(ierr);

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

	// create viscosity vectors
	ierr = DMCreateLocalVector(lvl->DA_CEN, &lvl->eta);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(lvl->DA_X,   &lvl->etax); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(lvl->DA_Y,   &lvl->etay); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(lvl->DA_Z,   &lvl->etaz); CHKERRQ(ierr);

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

	ierr = VecDestroy(&lvl->eta);          CHKERRQ(ierr);
	ierr = VecDestroy(&lvl->etax);         CHKERRQ(ierr);
	ierr = VecDestroy(&lvl->etay);         CHKERRQ(ierr);
	ierr = VecDestroy(&lvl->etaz);         CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelInitEta"
PetscErrorCode MGLevelInitEta(MGLevel *lvl, JacRes *jr)
{
	// initialize viscosity on fine grid

	PetscScalar ***eta;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize viscosity
	ierr = VecSet(lvl->eta, -1.0); CHKERRQ(ierr);

	// access viscosity vector
	ierr = DMDAVecGetArray(lvl->DA_CEN, lvl->eta, &eta); CHKERRQ(ierr);

	//----------
	// P-points
	//----------
	iter = 0;
	ierr = DMDAGetCorners (lvl->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		eta[k][j][i] = jr->svCell[iter++].svDev.eta;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_CEN, lvl->eta, &eta); CHKERRQ(ierr);

	// exchange ghost point values
	LOCAL_TO_LOCAL(lvl->DA_CEN, lvl->eta)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelAverageEta"
PetscErrorCode MGLevelAverageEta(MGLevel *lvl)
{
	// average viscosity from cell centers to velocity nodes

	PetscScalar b_eta, f_eta, n;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***eta, ***etax, ***etay, ***etaz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set viscosity
	ierr = VecSet(lvl->etax, -1.0); CHKERRQ(ierr);
	ierr = VecSet(lvl->etay, -1.0); CHKERRQ(ierr);
	ierr = VecSet(lvl->etaz, -1.0); CHKERRQ(ierr);

	// access viscosity vectors
	ierr = DMDAVecGetArray(lvl->DA_CEN, lvl->eta,  &eta);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_X,   lvl->etax, &etax); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,   lvl->etay, &etay); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,   lvl->etaz, &etaz); CHKERRQ(ierr);

	//---------
	// X-points
	//---------

	ierr = DMDAGetCorners(lvl->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		n     = 0.0;
		b_eta = eta[k][j][i-1]; if(b_eta != -1.0) { n += 1.0; } else { b_eta = 0.0; }
		f_eta = eta[k][j][i];   if(f_eta != -1.0) { n += 1.0; } else { f_eta = 0.0; }

		etax[k][j][i] = (b_eta + f_eta)/n;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	ierr = DMDAGetCorners(lvl->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		n     = 0.0;
		b_eta = eta[k][j-1][i]; if(b_eta != -1.0) { n += 1.0; } else { b_eta = 0.0; }
		f_eta = eta[k][j][i];   if(f_eta != -1.0) { n += 1.0; } else { f_eta = 0.0; }

		etay[k][j][i] = (b_eta + f_eta)/n;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	ierr = DMDAGetCorners(lvl->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		n     = 0.0;
		b_eta = eta[k-1][j][i]; if(b_eta != -1.0) { n += 1.0; } else { b_eta = 0.0; }
		f_eta = eta[k][j][i];   if(f_eta != -1.0) { n += 1.0; } else { f_eta = 0.0; }

		etaz[k][j][i] = (b_eta + f_eta)/n;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_CEN, lvl->eta,  &eta);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_X,   lvl->etax, &etax); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,   lvl->etay, &etay); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,   lvl->etaz, &etaz); CHKERRQ(ierr);

	// exchange ghost point viscosity
	LOCAL_TO_LOCAL(lvl->DA_X, lvl->etax)
	LOCAL_TO_LOCAL(lvl->DA_Y, lvl->etay)
	LOCAL_TO_LOCAL(lvl->DA_Z, lvl->etaz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelRestrictEta"
PetscErrorCode MGLevelRestrictEta(MGLevel *lvl, MGLevel *fine)
{
	// restrict inverse viscosity from fine to coarse level

	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, refine_y;
	PetscScalar sum, ***ceta, ***feta;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize viscosities
	ierr = VecSet(lvl->eta, -1.0); CHKERRQ(ierr);

	// access viscosity vector in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->eta,  &ceta); CHKERRQ(ierr);

	// access viscosity vector in fine grid
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	//-----------------------
	// P-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(lvl->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	ierr = DMDAGetRefinementFactor(fine->DA_CEN, NULL, &refine_y,NULL); CHKERRQ(ierr);	// refinement in y

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = refine_y*j;
		K = 2*k;

		// compute average viscosity
		sum = feta[K  ][J  ][I  ]
		+     feta[K  ][J  ][I+1]
		+     feta[K  ][J+1][I  ]
		+     feta[K  ][J+1][I+1]
		+     feta[K+1][J  ][I  ]
		+     feta[K+1][J  ][I+1]
		+     feta[K+1][J+1][I  ]
		+     feta[K+1][J+1][I+1];

		ceta[k][j][i] = sum/8.0;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->eta,  &ceta); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	// exchange ghost points
	LOCAL_TO_LOCAL(lvl->DA_CEN, lvl->eta)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelRestrictBC"
PetscErrorCode MGLevelRestrictBC(MGLevel *lvl, MGLevel *fine, PetscBool no_restric_bc)
{
	// restrict boundary condition vectors from fine to coarse level

	PetscInt    I, J, K, refine_y;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// mark all variables unconstrained
	ierr = VecSet(lvl->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcp,  DBL_MAX); CHKERRQ(ierr);

	// check activation
	if(no_restric_bc == PETSC_TRUE) PetscFunctionReturn(0);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	// Refinement factor in Y
	ierr = DMDAGetRefinementFactor(fine->DA_CEN, NULL, &refine_y,NULL); CHKERRQ(ierr);	// refinement in y

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(lvl->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = refine_y*j;
		K = 2*k;

		// restrict constraint
		if(fbcvx[K  ][J  ][I] != DBL_MAX
		&& fbcvx[K  ][J+1][I] != DBL_MAX
		&& fbcvx[K+1][J  ][I] != DBL_MAX
		&& fbcvx[K+1][J+1][I] != DBL_MAX)
		{
			// store parent DOF index
			cbcvx[k][j][i] = ivx[K][J][I];
		}
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
		J = refine_y*j;
		K = 2*k;

		// restrict constraint
		if(fbcvy[K  ][J][I  ] != DBL_MAX
		&& fbcvy[K  ][J][I+1] != DBL_MAX
		&& fbcvy[K+1][J][I  ] != DBL_MAX
		&& fbcvy[K+1][J][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvy[k][j][i] = ivy[K][J][I];
		}
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
		J = refine_y*j;
		K = 2*k;

		// restrict constraint
		if(fbcvz[K][J  ][I  ] != DBL_MAX
		&& fbcvz[K][J  ][I+1] != DBL_MAX
		&& fbcvz[K][J+1][I  ] != DBL_MAX
		&& fbcvz[K][J+1][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvz[k][j][i] = ivz[K][J][I];
		}
	}
	END_STD_LOOP

	if(lvl->dof.idxmod == IDXCOUPLED)
	{

		//-----------------------
		// P-points (coarse grid)
		//-----------------------
		ierr = DMDAGetCorners(lvl->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get fine grid indices
			I = 2*i;
			J = refine_y*j;
			K = 2*k;

			// restrict constraint
			if(fbcp[K  ][J  ][I  ] != DBL_MAX
			&& fbcp[K  ][J  ][I+1] != DBL_MAX
			&& fbcp[K  ][J+1][I  ] != DBL_MAX
			&& fbcp[K  ][J+1][I+1] != DBL_MAX
			&& fbcp[K+1][J  ][I  ] != DBL_MAX
			&& fbcp[K+1][J  ][I+1] != DBL_MAX
			&& fbcp[K+1][J+1][I  ] != DBL_MAX
			&& fbcp[K+1][J+1][I+1] != DBL_MAX)
			{
				// store parent DOF index
				cbcp[k][j][i] = ip[K][J][I];
			}
		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	// exchange ghost point constraints
	LOCAL_TO_LOCAL(lvl->DA_X,   lvl->bcvx)
	LOCAL_TO_LOCAL(lvl->DA_Y,   lvl->bcvy)
	LOCAL_TO_LOCAL(lvl->DA_Z,   lvl->bcvz)
	LOCAL_TO_LOCAL(lvl->DA_CEN, lvl->bcp)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelSetupRestrict"
PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine)
{
	Mat         R;
	PetscScalar v[12], bc[12], vs[12], eta[12];
	PetscInt    idx[12], refine_y;
	PetscInt    row, I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;
	PetscScalar ***fetax, ***fetay, ***fetaz;
	PetscScalar ***cetax, ***cetay, ***cetaz;

	PetscBool scale = PETSC_FALSE;
	PetscOptionsGetBool(NULL, NULL, "-rest", &scale, NULL);

	PetscErrorCode ierr;
	PetscFunctionBegin;

	R = lvl->R;

	// get refinement factor in y-direction of central array (in 2D, don't refine in y-direction)
	ierr = DMDAGetRefinementFactor(fine->DA_CEN, NULL, &refine_y,NULL); CHKERRQ(ierr);

	// clear restriction matrix coefficients
	ierr = MatZeroEntries(R); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	// access viscosity vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->etax,    &fetax); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->etay,    &fetay); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->etaz,    &fetaz); CHKERRQ(ierr);

	// access viscosity vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->etax,     &cetax); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->etay,     &cetay); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->etaz,     &cetaz); CHKERRQ(ierr);

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
		J = refine_y*j;
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

		// get fine grid viscosity
		eta[0]   =        fetax[K  ][J  ][I-1];
		eta[1]   =        fetax[K  ][J+1][I-1];
		eta[2]   =        fetax[K+1][J  ][I-1];
		eta[3]   =        fetax[K+1][J+1][I-1];
		eta[4]   =        fetax[K  ][J  ][I  ];
		eta[5]   =        fetax[K  ][J+1][I  ];
		eta[6]   =        fetax[K+1][J  ][I  ];
		eta[7]   =        fetax[K+1][J+1][I  ];
		eta[8]   =        fetax[K  ][J  ][I+1];
		eta[9]   =        fetax[K  ][J+1][I+1];
		eta[10]  =        fetax[K+1][J  ][I+1];
		eta[11]  =        fetax[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(scale, cbcvx[k][j][i], 12, idx, bc, v, vs, eta, cetax[k][j][i]);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

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
		J = refine_y*j;
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

    	// get fine grid viscosity
		eta[0]   =        fetay[K  ][J-1][I  ];
		eta[1]   =        fetay[K  ][J-1][I+1];
		eta[2]   =        fetay[K+1][J-1][I  ];
		eta[3]   =        fetay[K+1][J-1][I+1];
		eta[4]   =        fetay[K  ][J  ][I  ];
		eta[5]   =        fetay[K  ][J  ][I+1];
		eta[6]   =        fetay[K+1][J  ][I  ];
		eta[7]   =        fetay[K+1][J  ][I+1];
		eta[8]   =        fetay[K  ][J+1][I  ];
		eta[9]   =        fetay[K  ][J+1][I+1];
		eta[10]  =        fetay[K+1][J+1][I  ];
    	eta[11]  =        fetay[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(scale, cbcvy[k][j][i], 12, idx, bc, v, vs, eta, cetay[k][j][i]);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

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
		J = refine_y*j;
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

    	// get fine grid viscosity
		eta[0]   =        fetaz[K-1][J  ][I  ];
		eta[1]   =        fetaz[K-1][J  ][I+1];
		eta[2]   =        fetaz[K-1][J+1][I  ];
		eta[3]   =        fetaz[K-1][J+1][I+1];
		eta[4]   =        fetaz[K  ][J  ][I  ];
		eta[5]   =        fetaz[K  ][J  ][I+1];
		eta[6]   =        fetaz[K  ][J+1][I  ];
		eta[7]   =        fetaz[K  ][J+1][I+1];
		eta[8]   =        fetaz[K+1][J  ][I  ];
		eta[9]   =        fetaz[K+1][J  ][I+1];
		eta[10]  =        fetaz[K+1][J+1][I  ];
    	eta[11]  =        fetaz[K+1][J+1][I+1];

		// setup row of restriction matrix
		getRowRestrict(scale, cbcvz[k][j][i], 12, idx, bc, v, vs, eta, cetaz[k][j][i]);

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, INSERT_VALUES); CHKERRQ(ierr);

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
			J = refine_y*j;
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
			getRowRestrict(PETSC_FALSE, cbcp[k][j][i], 8, idx, bc, v, vs, NULL, 0.0);

			// store full matrix row
			ierr = MatSetValues(R, 1, &row, 8, idx, v, INSERT_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->etax,    &fetax); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->etay,    &fetay); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->etaz,    &fetaz); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->etax,     &cetax); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->etay,     &cetay); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->etaz,     &cetaz); CHKERRQ(ierr);

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
	PetscInt    n, idx[8], refine_y;
	PetscScalar v[8], bc[8], eta[8], vsf[8], vsr[4], *vs;
	PetscInt    row, I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;
	PetscScalar ***fetax, ***fetay, ***fetaz;
	PetscScalar ***cetax, ***cetay, ***cetaz;

	PetscBool scale = PETSC_FALSE;
	PetscOptionsGetBool(NULL, NULL, "-prol", &scale, NULL);

	PetscErrorCode ierr;
	PetscFunctionBegin;

	P = lvl->P;

	// clear prolongation matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);

	// get refinement factor in y-direction of central array (in 2D, don't refine in y-direction)
	ierr = DMDAGetRefinementFactor(fine->DA_CEN, NULL, &refine_y,NULL); CHKERRQ(ierr);

	// access index vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->dof.ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->bcvx,   &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->bcvy,   &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->bcvz,   &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->bcp,    &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->bcvx,    &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->bcvy,    &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->bcvz,    &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->bcp,     &cbcp);  CHKERRQ(ierr);

	// access viscosity vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->etax,   &fetax); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->etay,   &fetay); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->etaz,   &fetaz); CHKERRQ(ierr);

	// access viscosity vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->etax,    &cetax); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->etay,    &cetay); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->etaz,    &cetaz); CHKERRQ(ierr);

	// get global index of the first row in the fine grid
	if     (fine->dof.idxmod == IDXCOUPLED)   { row = fine->dof.st;  }
	else if(fine->dof.idxmod == IDXUNCOUPLED) { row = fine->dof.stv; }

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
	ierr = DMDAGetCorners(fine->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/refine_y;
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
		eta[0] =         cetax[K ][J ][I];
		eta[1] =         cetax[K ][J1][I];
		eta[2] =         cetax[K1][J ][I];
		eta[3] =         cetax[K1][J1][I];

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
			eta[4] =         cetax[K ][J ][I1];
			eta[5] =         cetax[K ][J1][I1];
			eta[6] =         cetax[K1][J ][I1];
			eta[7] =         cetax[K1][J1][I1];
		}

		// setup row of prolongation matrix
		getRowProlong(scale, row, fbcvx[k][j][i], n, bc, v, vs, eta, fetax[k][j][i]);


		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

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
		J = j/refine_y;
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
		eta[0] =         cetay[K ][J][I ];
		eta[1] =         cetay[K ][J][I1];
		eta[2] =         cetay[K1][J][I ];
		eta[3] =         cetay[K1][J][I1];

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
			eta[4] =         cetay[K ][J1][I ];
			eta[5] =         cetay[K ][J1][I1];
			eta[6] =         cetay[K1][J1][I ];
			eta[7] =         cetay[K1][J1][I1];

		}

		// setup row of prolongation matrix
		getRowProlong(scale, row, fbcvy[k][j][i], n, bc, v, vs, eta, fetay[k][j][i]);

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

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
		J = j/refine_y;
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
		eta[0] =         cetaz[K][J ][I ];
		eta[1] =         cetaz[K][J ][I1];
		eta[2] =         cetaz[K][J1][I ];
		eta[3] =         cetaz[K][J1][I1];

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
			eta[4] =         cetaz[K1][J ][I ];
			eta[5] =         cetaz[K1][J ][I1];
			eta[6] =         cetaz[K1][J1][I ];
			eta[7] =         cetaz[K1][J1][I1];
		}

		// setup row of prolongation matrix
		getRowProlong(scale, row, fbcvz[k][j][i], n, bc, v, vs, eta, fetaz[k][j][i]);

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, n, idx, v, INSERT_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP


	if(fine->dof.idxmod == IDXCOUPLED)
	{
		// set pressure interpolation stencil (direct injection)
		vsr[0] = 1.0;

		//---------------------
		// P-points (fine grid)
		//---------------------
		ierr = DMDAGetCorners(fine->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get coarse grid indices
			I = i/2;
			J = j/refine_y;
			K = k/2;

			idx[0] = (PetscInt)ip[K][J][I];
			bc [0] =         cbcp[K][J][I];

			// setup row of prolongation matrix
			getRowProlong(PETSC_FALSE, row, fbcp[k][j][i], 1, bc, v, vsr, NULL, 0.0);

			// store full matrix row
			ierr = MatSetValues(P, 1, &row, 1, idx, v, INSERT_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->dof.ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->bcvx,   &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->bcvy,   &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->bcvz,   &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->bcp,    &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->bcvx,    &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->bcvy,    &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->bcvz,    &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->bcp,     &cbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->etax,   &fetax); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->etay,   &fetay); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->etaz,   &fetaz); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->etax,    &cetax); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->etay,    &cetay); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->etaz,    &cetaz); CHKERRQ(ierr);

	// assemble prolongation matrix
	ierr = MatAIJAssemble(P, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void getRowRestrict(PetscBool scale,
	PetscScalar parent, PetscInt n, PetscInt idx[], PetscScalar bc[],
	PetscScalar v[], PetscScalar vs[], PetscScalar eta_fine[], PetscScalar eta_crs)
{
	PetscInt j, pdof;
	PetscScalar sum;

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
		// scale with local viscosities
		if(scale == PETSC_TRUE)
		{
			sum = 0.0;

			for(j = 0; j < n; j++)
			{
				v[j] *= eta_fine[j]/eta_crs;

				sum += v[j];

			}

			// normalize
			for(j = 0; j < n; j++)
			{
				v[j] /= sum;
			}

		}
	}
}
//---------------------------------------------------------------------------
void getRowProlong(PetscBool scale,
	PetscInt parent, PetscScalar parent_bc, PetscInt n, PetscScalar bc[],
	PetscScalar v[], PetscScalar vs[], PetscScalar eta_crs[], PetscScalar eta_fine)
{
	PetscInt    j;
	PetscScalar pdof, sum;

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

		// scale with local viscosities, compute row sum
		if(scale == PETSC_TRUE)
		{
			sum = 0.0;

			for(j = 0; j < n; j++)
			{
				v[j] *= eta_crs[j]/eta_fine;

				sum += v[j];
			}

			// normalize
			for(j = 0; j < n; j++)
			{
				v[j] /= sum;
			}

		}
	}
}
//---------------------------------------------------------------------------
/*
 		// check neighbor points
		nd = 1;
		no = 0;
		// x-velocity
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i-1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j+1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j-1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k+1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k-1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k][j][i-1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j+1][i-1]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, start, ln, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j][i-1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j][i-1]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
 		// pressure
		ind = (PetscInt) ip[k][j][i-1];    CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, start, ln, nd, no);
		// update counters
		d_nnz[iter] = nd;
		o_nnz[iter] = no;
		iter++;
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelAllocRestrict"
PetscErrorCode MGLevelAllocRestrict(MGLevel *lvl, MGLevel *fine)
{
	Mat         R;
	PetscInt    ind;
	PetscInt    row, I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	R = lvl->R;

	// clear restriction matrix coefficients
	ierr = MatZeroEntries(R); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);

	// get global index of the first row in coarse grid
	if     (lvl->dof.idxmod == IDXCOUPLED)   { row = lvl->dof.st;  }
	else if(lvl->dof.idxmod == IDXUNCOUPLED) { row = lvl->dof.stv; }

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

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(lvl->dof.idxmod == IDXCOUPLED)
	{
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
			idx[0] = (PetscInt)ip[K  ][J  ][I  ];
			idx[1] = (PetscInt)ip[K  ][J  ][I+1];
			idx[2] = (PetscInt)ip[K  ][J+1][I  ];
			idx[3] = (PetscInt)ip[K  ][J+1][I+1];
			idx[4] = (PetscInt)ip[K+1][J  ][I  ];
			idx[5] = (PetscInt)ip[K+1][J  ][I+1];
			idx[6] = (PetscInt)ip[K+1][J+1][I  ];
			idx[7] = (PetscInt)ip[K+1][J+1][I+1];

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->DA_X,   fine->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   fine->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   fine->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, fine->dof.ip,  &ip);    CHKERRQ(ierr);


	// assemble restriction matrix
	ierr = MatAIJAssemble(R, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGLevelAllocProlong"
PetscErrorCode MGLevelAllocProlong(MGLevel *lvl, MGLevel *fine)
{
	Mat P;
	PetscInt    ind;
	PetscInt    row, I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	P = lvl->P;

	// clear prolongation matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);

	// access index vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->DA_X,    lvl->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Y,    lvl->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_Z,    lvl->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->DA_CEN,  lvl->dof.ip,  &ip);    CHKERRQ(ierr);


	// get global index of the first row in the fine grid
	if     (fine->dof.idxmod == IDXCOUPLED)   { row = fine->dof.st;  }
	else if(fine->dof.idxmod == IDXUNCOUPLED) { row = fine->dof.stv; }


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

		if(j % 2) J1 = J+1; else J1 = J-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivx[K ][J ][I];
		idx[1] = (PetscInt)ivx[K ][J1][I];
		idx[2] = (PetscInt)ivx[K1][J ][I];
		idx[3] = (PetscInt)ivx[K1][J1][I];

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
		}

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

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivy[K ][J][I ];
		idx[1] = (PetscInt)ivy[K ][J][I1];
		idx[2] = (PetscInt)ivy[K1][J][I ];
		idx[3] = (PetscInt)ivy[K1][J][I1];

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

		}

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

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(j % 2) J1 = J+1; else J1 = J-1;

		// setup reduced stencil (even)
		n      = 4;
		vs     = vsr;
		idx[0] = (PetscInt)ivz[K][J ][I ];
		idx[1] = (PetscInt)ivz[K][J ][I1];
		idx[2] = (PetscInt)ivz[K][J1][I ];
		idx[3] = (PetscInt)ivz[K][J1][I1];

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
		}

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(fine->dof.idxmod == IDXCOUPLED)
	{

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

			// increment row number
			row++;

		}
		END_STD_LOOP

	}

	// restore access
	ierr = DMDAVecRestoreArray(lvl->DA_X,    lvl->dof.ivx,  &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Y,    lvl->dof.ivy,  &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_Z,    lvl->dof.ivz,  &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->DA_CEN,  lvl->dof.ip,   &ip);    CHKERRQ(ierr);

	// assemble prolongation matrix
	ierr = MatAIJAssemble(P, 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
// MG -functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCreate"
PetscErrorCode MGCreate(MG *mg, JacRes *jr)
{
	PetscInt  i, l;
	MGLevel   *fine;
	char      pc_type[_str_len_];
	PetscBool opt_set;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get preconditioner type
	ierr = PetscOptionsGetString(NULL, NULL, "-gmg_pc_type", pc_type, _str_len_, &opt_set); CHKERRQ(ierr);

	// check whether multigrid is requested
	if(opt_set != PETSC_TRUE || strcmp(pc_type, "mg"))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "-gmg_pc_type option is not defined of specified incorrectly (use -gmg_pc_type mg)");
	}

	// clear object
	ierr = PetscMemzero(mg, sizeof(MG)); CHKERRQ(ierr);

	// store finest grid context
	mg->jr = jr;

	// set boundary constraint restriction flag
	ierr = PetscOptionsHasName(NULL, NULL, "-gmg_no_restric_bc", &mg->no_restric_bc); CHKERRQ(ierr);

	// check multigrid mesh restrictions & get actual number of levels
	ierr = MGGetNumLevels(mg); CHKERRQ(ierr);

	// allocate levels
	ierr = PetscMalloc(sizeof(MGLevel)*(size_t)mg->nlvl, &mg->lvls); CHKERRQ(ierr);

	// create levels
	fine = NULL;

	for(i = 0; i < mg->nlvl; i++)
	{
		ierr = MGLevelCreate(&mg->lvls[i], fine, jr->fs, jr->bc); CHKERRQ(ierr);

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
#undef __FUNCT__
#define __FUNCT__ "MGDestroy"
PetscErrorCode MGDestroy(MG *mg)
{
	PetscInt  i;
	PetscBool flg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "MGSetupCoarse"
PetscErrorCode MGSetupCoarse(MG *mg, Mat A)
{
	KSP        ksp;
	PC         pc;
	Mat        mat;
	MGLevel   *lvl;
	DOFIndex  *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// skip already configured solver
	if(mg->crs_setup == PETSC_TRUE)
	{
		PetscFunctionReturn(0);
	}

	// get coarse level index object
	lvl = mg->lvls + mg->nlvl - 1;
	dof = &lvl->dof;

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
	ierr = MatAIJSetNullSpace(mat, dof);     CHKERRQ(ierr);

	// set actual coarse solver options
	ierr = KSPSetOptionsPrefix(ksp, "crs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);           CHKERRQ(ierr);

	// set setup flag
	mg->crs_setup = PETSC_TRUE;

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

	ierr = MGLevelInitEta(mg->lvls, mg->jr); CHKERRQ(ierr);
	ierr = MGLevelAverageEta(mg->lvls);      CHKERRQ(ierr);

	for(i = 1; i < mg->nlvl; i++)
	{
		ierr = MGLevelRestrictBC   (&mg->lvls[i], &mg->lvls[i-1], mg->no_restric_bc); CHKERRQ(ierr);
		ierr = MGLevelRestrictEta  (&mg->lvls[i], &mg->lvls[i-1]);                    CHKERRQ(ierr);
		ierr = MGLevelAverageEta   (&mg->lvls[i]);                                    CHKERRQ(ierr);
		ierr = MGLevelSetupRestrict(&mg->lvls[i], &mg->lvls[i-1]);                    CHKERRQ(ierr);
		ierr = MGLevelSetupProlong (&mg->lvls[i], &mg->lvls[i-1]);                    CHKERRQ(ierr);
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
#undef __FUNCT__
#define __FUNCT__ "MGGetNumLevels"
PetscErrorCode MGGetNumLevels(MG *mg)
{
	// check multigrid mesh restrictions, get actual number of coarsening steps

	FDSTAG   *fs;
	PetscBool opt_set;
	PetscInt  nx, ny, nz, Nx, Ny, Nz, ncors, nlevels, refine_y;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = mg->jr->fs;

	// get refinement factor in y-direction of central array (in 2D, don't refine in y-direction)
	refine_y = 2;
	ierr = PetscOptionsGetInt(NULL, NULL, "-da_refine_y", &refine_y, NULL); CHKERRQ(ierr);		// for cases where refinement is set to 1 (2D MG)

	// check discretization in all directions
	ierr = Discret1DCheckMG(&fs->dsx, "x", &nx); CHKERRQ(ierr);                ncors = nx;

	if(refine_y > 1)
	{
		ierr = Discret1DCheckMG(&fs->dsy, "y", &ny); CHKERRQ(ierr);
		if(ny < ncors) ncors = ny;
	}

	ierr = Discret1DCheckMG(&fs->dsz, "z", &nz); CHKERRQ(ierr); if(nz < ncors) ncors = nz;

	// check number of levels requested on the command line
	ierr = PetscOptionsGetInt(NULL, NULL, "-gmg_pc_mg_levels", &nlevels, &opt_set); CHKERRQ(ierr);

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
	if  (refine_y>1){
		ny = fs->dsy.ncels >> ncors;
	}
	else
	{
		ny = fs->dsy.ncels;
	}
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
