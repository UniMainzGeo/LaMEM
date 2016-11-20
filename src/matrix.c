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
 **    filename:   matrix.c
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
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "matrix.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * pressure Schur complement preconditioners
// * matrix-free preconditioner action
// * linear system scaling (fdstag multigrid paper)
// * temperature scaling
// * preallocation for temperature & pressure Schur PC
// * figure out why block factorization with penalty
//   doesn't work with non-homogeneous Dirichlet BC
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatAIJCreate"
PetscErrorCode MatAIJCreate(
	PetscInt m, PetscInt n,
	PetscInt d_nz, const PetscInt d_nnz[],
	PetscInt o_nz, const PetscInt o_nnz[],
	Mat *P)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create matrix
	ierr = MatCreate(PETSC_COMM_WORLD, P); CHKERRQ(ierr);
	ierr = MatSetType((*P), MATAIJ); CHKERRQ(ierr);
	ierr = MatSetSizes((*P), m, n, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);

	// preallocate matrix
	ierr = MatSeqAIJSetPreallocation((*P), d_nz, d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation((*P), d_nz, d_nnz, o_nz, o_nnz); CHKERRQ(ierr);

	// read custom options (required to resolve SuperLU_DIST issue)
	ierr = MatSetFromOptions((*P)); CHKERRQ(ierr);

	// throw an error if preallocation fails
	ierr = MatSetOption((*P), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetUp((*P)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatAIJCreateDiag"
PetscErrorCode MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P)
{
	PetscInt i, ii;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// preallocate
	ierr = MatAIJCreate(m, m, 1, NULL, 0, NULL, P); CHKERRQ(ierr);

	// put explicit zeroes on the diagonal
	for(i = 0; i < m; i++)
	{
		// get global row index
		ii = istart + i;

		ierr = MatSetValue((*P), ii, ii, 0.0, INSERT_VALUES); CHKERRQ(ierr);
	}

	// assemble
	ierr = MatAIJAssemble((*P), 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatAIJAssemble"
PetscErrorCode MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag)
{
//	PetscInt    m, n;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// freeze nonzero structure, constrain rows only locally
	ierr = MatSetOption(P, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);     CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);    CHKERRQ(ierr);

	// zero out constrained rows, form unit diagonal for the constrained block
	ierr = MatZeroRows(P, numRows, rows, diag, NULL, NULL); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatCreate"
PetscErrorCode PMatCreate(PMat *p_pm, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	//========================================================================
	// create preconditioner matrix context
	//========================================================================

	PMat pm;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PMat), &pm); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pm, sizeof(p_PMat)); CHKERRQ(ierr);

	// set type
	ierr = PMatSetFromOptions(pm); CHKERRQ(ierr);

	// set context
	pm->jr = jr;

	if(pm->type == _MONOLITHIC_)
	{
		// monolithic format
		pm->Create   = PMatMonoCreate;
		pm->Assemble = PMatMonoAssemble;
		pm->Destroy  = PMatMonoDestroy;
		pm->Picard   = PMatMonoPicard;
	}
	else if(pm->type == _BLOCK_)
	{
		// block format
		pm->Create   = PMatBlockCreate;
		pm->Assemble = PMatBlockAssemble;
		pm->Destroy  = PMatBlockDestroy;
		if(pm->pgamma != 1.0) pm->Picard = PMatBlockPicardSchur;
		else                  pm->Picard = PMatBlockPicardClean;
	}
	// create type-specific context
	ierr = pm->Create(pm); CHKERRQ(ierr);

	// return pointer
	(*p_pm) = pm;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatSetFromOptions"
PetscErrorCode PMatSetFromOptions(PMat pm)
{
	PetscBool   flg;
	PetscScalar pgamma;
	char        pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set matrix type
	ierr = PetscOptionsGetString(NULL, NULL,"-pcmat_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(!strcmp(pname, "mono"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner matrix type     : monolithic\n");
			pm->type = _MONOLITHIC_;
		}
		else if(!strcmp(pname, "block"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner matrix type     : block\n");
			pm->type = _BLOCK_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Incorrect matrix storage format: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Preconditioner matrix type     : monolithic\n");
		pm->type = _MONOLITHIC_;
	}

	// set penalty parameter
	pm->pgamma = 1.0;

	ierr = PetscOptionsGetScalar(NULL, NULL, "-pcmat_pgamma", &pgamma, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(pgamma < 1.0)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Penalty parameter [-pcmat_pgamma] is less than unit");
		}

		pm->pgamma = pgamma;
	}

	if(pm->pgamma > 1.0)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Penalty parameter (pgamma)     : %e\n", pm->pgamma);
	}

	// set cell stiffness function
	ierr = PetscOptionsHasName(NULL, NULL, "-pcmat_no_dev_proj", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Excluding deviatoric projection from preconditioner\n");
		pm->getStiffMat = getStiffMatClean;
	}
	else
	{
		pm->getStiffMat = getStiffMatDevProj;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatAssemble"
PetscErrorCode PMatAssemble(PMat pm)
{
	BCCtx  *bc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

//	PetscPrintf(PETSC_COMM_WORLD, " Starting preconditioner assembly\n");

	bc = pm->jr->bc;

	// shift constrained node indices to global index space
	ierr = BCShiftIndices(bc, _LOCAL_TO_GLOBAL_); CHKERRQ(ierr);

	ierr = pm->Assemble(pm); CHKERRQ(ierr);

	// shift constrained node indices back to local index space
	ierr = BCShiftIndices(bc,  _GLOBAL_TO_LOCAL_); CHKERRQ(ierr);

//	PetscPrintf(PETSC_COMM_WORLD, " Finished preconditioner assembly\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatDestroy"
PetscErrorCode PMatDestroy(PMat pm)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = pm->Destroy(pm); CHKERRQ(ierr);
	ierr = PetscFree(pm);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//.........................   MONOLITHIC MATRIX   ...........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatMonoCreate"
PetscErrorCode PMatMonoCreate(PMat pm)
{
	//=========================================================================
	// Count nonzero in diagonal and off-diagonal blocks
	//
	// Stencil of every velocity point contains 17 points:
	//    * 7 points from the same velocity component (star-type, self included)
	//    * 4 velocity points from each of the other two orthogonal directions
	//    * 2 pressure points
	//
	// Stencil of every pressure point contains 7 points:
	//    * 6 velocity points (2 points from each velocity direction)
	//    * 1 pressure point (self)
	//
	// We just need to check which of the neighboring points are local,
	// and which are not, and set the counters accordingly.
	// Some of the points do not exist near the boundaries.
	//
	//=========================================================================

	FDSTAG      *fs;
	DOFIndex    *dof;
	PMatMono    *P;
	PetscInt    ln, start, ind, nd, no, *d_nnz, *o_nnz;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access contexts
	fs  = pm->jr->fs;
	dof = &fs->dof;

	// allocate space
	ierr = PetscMalloc(sizeof(PMatMono), (void**)&P); CHKERRQ(ierr);

	// store context
	pm->data = (void*)P;

	// compute global indexing
	ierr = DOFIndexCompute(dof, IDXCOUPLED); CHKERRQ(ierr);

	// get number of local rows & global index of the first row
	ln    = dof->ln;
	start = dof->st;

	// allocate nonzero counter arrays
	ierr = makeIntArray(&d_nnz, NULL, ln); CHKERRQ(ierr);
	ierr = makeIntArray(&o_nnz, NULL, ln); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// clear iterator
	iter = 0;

	//---------
	// X-points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
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
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// check neighbor points
		nd = 1;
		no = 0;
		// y-velocity
		ind = (PetscInt) ivy[k][j][i+1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j][i-1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j-1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k+1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k-1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j-1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j-1][i]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
 		// x-velocity
		ind = (PetscInt) ivx[k][j-1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j-1][i+1]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, start, ln, nd, no);
		// pressure
		ind = (PetscInt) ip[k][j-1][i];    CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, start, ln, nd, no);
		// update counters
		d_nnz[iter] = nd;
		o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// check neighbor points
		nd = 1;
		no = 0;
		// z-velocity
		ind = (PetscInt) ivz[k][j][i+1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k][j][i-1];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k][j+1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k][j-1][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k-1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
 		// x-velocity
		ind = (PetscInt) ivx[k-1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k-1][j][i+1]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, start, ln, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k-1][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k-1][j+1][i]; CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j][i];     CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, start, ln, nd, no);
		// pressure
		ind = (PetscInt) ip[k-1][j][i];    CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, start, ln, nd, no);
		// update counters
		d_nnz[iter] = nd;
		o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// check neighbor points
		nd = 1;
		no = 0;
		// x-velocity
		ind = (PetscInt) ivx[k][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivx[k][j][i+1]; CHECK_DOF(ind, start, ln, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivy[k][j+1][i]; CHECK_DOF(ind, start, ln, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j][i];   CHECK_DOF(ind, start, ln, nd, no);
		ind = (PetscInt) ivz[k+1][j][i]; CHECK_DOF(ind, start, ln, nd, no);
		// update counters
		d_nnz[iter] = nd;
		o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// create matrices & vectors
	ierr = MatAIJCreate(ln, ln, 0, d_nnz, 0, o_nnz, &P->A); CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(ln, start, &P->M);              CHKERRQ(ierr);
	ierr = VecDuplicate(pm->jr->gsol, &P->w);               CHKERRQ(ierr);

	// clear work arrays
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatMonoAssemble"
PetscErrorCode PMatMonoAssemble(PMat pm)
{
	//======================================================================
	// Assemble effective viscosity preconditioning matrix
	//
	// Global ordering of the variables is interlaced:
	// all X-Y-Z-P DOF of the first processor
	// are followed by X-Y-Z-P DOF of the second one, and so on.
	//
	// Picard Jacobian (J) can be computed as follows when necessary:
	// J = P - M
	// P - preconditioner matrix (computed in this function)
	// M - inverse viscosity matrix (computed in this function)
	//======================================================================

	JacRes      *jr;
	FDSTAG      *fs;
	BCCtx       *bc;
	DOFIndex    *dof;
	PMatMono    *P;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, rho, IKdt, diag, pgamma, pt, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access contexts
	jr     = pm->jr;
	fs     = jr->fs;
	bc     = jr->bc;
	dof    = &fs->dof;
	P      = (PMatMono*)pm->data;

	// get density gradient stabilization parameters
	dt   = jr->ts->dt; // time step
	fssa = jr->FSSA;   // density gradient penalty parameter
    grav = jr->grav;   // gravity acceleration

	// get penalty parameter
	pgamma = pm->pgamma;

	// clear matrix coefficients
	ierr = MatZeroEntries(P->A); CHKERRQ(ierr);
	ierr = MatZeroEntries(P->M); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,   &bcp);   CHKERRQ(ierr);

	//---------------
	// central points
	//---------------

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		eta  = jr->svCell[iter].svDev.eta;
		IKdt = jr->svCell[iter].svBulk.IKdt;
		rho  = jr->svCell[iter].svBulk.rho;

		iter++;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute penalty term
		pt = -1.0/(pgamma*eta);

		// get pressure diagonal element (with penalty)
		diag = -IKdt + pt;

		// compute local matrix
		pm->getStiffMat(eta, diag, v, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

		// compute density gradient stabilization terms
		addDensGradStabil(fssa, v, rho, dt, grav, fdx, fdy, fdz, bdx, bdy, bdz);

		// get global indices of the points:
		// vx_(i), vx_(i+1), vy_(j), vy_(j+1), vz_(k), vz_(k+1), p
		idx[0] = (PetscInt) ivx[k][j][i];
		idx[1] = (PetscInt) ivx[k][j][i+1];
		idx[2] = (PetscInt) ivy[k][j][i];
		idx[3] = (PetscInt) ivy[k][j+1][i];
		idx[4] = (PetscInt) ivz[k][j][i];
		idx[5] = (PetscInt) ivz[k+1][j][i];
		idx[6] = (PetscInt) ip[k][j][i];

		// get boundary constraints
		pdofidx[0] = -1;   cf[0] = bcvx[k][j][i];
		pdofidx[1] = -1;   cf[1] = bcvx[k][j][i+1];
		pdofidx[2] = -1;   cf[2] = bcvy[k][j][i];
		pdofidx[3] = -1;   cf[3] = bcvy[k][j+1][i];
		pdofidx[4] = -1;   cf[4] = bcvz[k][j][i];
		pdofidx[5] = -1;   cf[5] = bcvz[k+1][j][i];
		pdofidx[6] = -1;   cf[6] = bcp[k][j][i];

		// constrain local matrix
		constrLocalMat(7, pdofidx, cf, v);

		// update global & penalty compensation matrices
		ierr = MatSetValues(P->A, 7, idx, 7, idx, v,  ADD_VALUES);    CHKERRQ(ierr);
		ierr = MatSetValue (P->M, idx[6], idx[6], pt, INSERT_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svXYEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// compute local matrix
		//       vx_(j-1)             vx_(j)               vy_(i-1)             vy_(i)
		v[0]  =  eta/dy/bdy; v[1]  = -eta/dy/bdy; v[2]  =  eta/dx/bdy; v[3]  = -eta/dx/bdy; // fx_(j-1) [sxy]
		v[4]  = -eta/dy/fdy; v[5]  =  eta/dy/fdy; v[6]  = -eta/dx/fdy; v[7]  =  eta/dx/fdy; // fx_(j)   [sxy]
		v[8]  =  eta/dy/bdx; v[9]  = -eta/dy/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fy_(i-1) [sxy]
		v[12] = -eta/dy/fdx; v[13] =  eta/dy/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fy_(i)   [sxy]

		// get global indices of the points: vx_(j-1), vx_(j), vy_(i-1), vy_(i)
		idx[0] = (PetscInt) ivx[k][j-1][i];
		idx[1] = (PetscInt) ivx[k][j][i];
		idx[2] = (PetscInt) ivy[k][j][i-1];
		idx[3] = (PetscInt) ivy[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k][j-1][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvy[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvy[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svXZEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// compute local matrix
		//       vx_(k-1)             vx_(k)               vz_(i-1)             vz_(i)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dx/bdz; v[3]  = -eta/dx/bdz; // fx_(k-1) [sxz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dx/fdz; v[7]  =  eta/dx/fdz; // fx_(k)   [sxz]
		v[8]  =  eta/dz/bdx; v[9]  = -eta/dz/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fz_(i-1) [sxz]
		v[12] = -eta/dz/fdx; v[13] =  eta/dz/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fz_(i)   [sxz]

		// get global indices of the points: vx_(k-1), vx_(k), vz_(i-1), vz_(i)
		idx[0] = (PetscInt) ivx[k-1][j][i];
		idx[1] = (PetscInt) ivx[k][j][i];
		idx[2] = (PetscInt) ivz[k][j][i-1];
		idx[3] = (PetscInt) ivz[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svYZEdge[iter++].svDev.eta;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// compute local matrix
		//       vy_(k-1)             vy_(k)               vz_(j-1)             vz_(j)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dy/bdz; v[3]  = -eta/dy/bdz; // fy_(k-1) [syz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dy/fdz; v[7]  =  eta/dy/fdz; // fy_(k)   [syz]
		v[8]  =  eta/dz/bdy; v[9]  = -eta/dz/bdy; v[10] =  eta/dy/bdy; v[11] = -eta/dy/bdy; // fz_(j-1) [syz]
		v[12] = -eta/dz/fdy; v[13] =  eta/dz/fdy; v[14] = -eta/dy/fdy; v[15] =  eta/dy/fdy; // fz_(j)   [syz]

		// get global indices of the points: vy_(k-1), vy_(k), vz_(j-1), vz_(j)
		idx[0] = (PetscInt) ivy[k-1][j][i];
		idx[1] = (PetscInt) ivy[k][j][i];
		idx[2] = (PetscInt) ivz[k][j-1][i];
		idx[3] = (PetscInt) ivz[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvy[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvy[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j-1][i];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp);  CHKERRQ(ierr);

	// assemble velocity-pressure matrix, remove constrained rows
	ierr = MatAIJAssemble(P->A, bc->numSPC, bc->SPCList, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->M, bc->numSPC, bc->SPCList, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatMonoPicard"
PetscErrorCode PMatMonoPicard(Mat J, Vec x, Vec r)
{
	// actual operation is: r = J*x = A*x - M*x

	PMatMono *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(J, (void**)&P); CHKERRQ(ierr);

	// compute action of preconditioner matrix
	ierr = MatMult(P->A, x, r); CHKERRQ(ierr);

	// compute compensation
	ierr = MatMult(P->M, x, P->w); CHKERRQ(ierr);

	// update result
	ierr = VecAXPY(r, -1.0, P->w);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatMonoDestroy"
PetscErrorCode PMatMonoDestroy(PMat pm)
{
	PMatMono *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	P = (PMatMono*)pm->data;

	ierr = MatDestroy (&P->A); CHKERRQ(ierr);
	ierr = MatDestroy (&P->M); CHKERRQ(ierr);
	ierr = VecDestroy (&P->w); CHKERRQ(ierr);
	ierr = PetscFree(P);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//...........................   BLOCK MATRIX   ..............................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBlockCreate"
PetscErrorCode PMatBlockCreate(PMat pm)
{
	JacRes      *jr;
	FDSTAG      *fs;
	DOFIndex    *dof;
	PMatBlock   *P;
	PetscInt    lnp, lnv, startv, startp, nd, no, ind;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    *Avv_d_nnz, *Avv_o_nnz;
	PetscInt    *Avp_d_nnz, *Avp_o_nnz;
	PetscInt    *Apv_d_nnz, *Apv_o_nnz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	jr  = pm->jr;
	fs  = jr->fs;
	dof = &fs->dof;

	// allocate space
	ierr = PetscMalloc(sizeof(PMatBlock), (void**)&P); CHKERRQ(ierr);

	// store context
	pm->data = (void*)P;

	// compute global indexing
	ierr = DOFIndexCompute(dof, IDXUNCOUPLED); CHKERRQ(ierr);

	// get number of local rows & global index of the first row
	lnv    = dof->lnv;
	startv = dof->stv;
	lnp    = dof->lnp;
	startp = dof->stp;

	// allocate nonzero counter arrays
	ierr = makeIntArray(&Avv_d_nnz, NULL, lnv); CHKERRQ(ierr);
	ierr = makeIntArray(&Avv_o_nnz, NULL, lnv); CHKERRQ(ierr);

	ierr = makeIntArray(&Avp_d_nnz, NULL, lnv); CHKERRQ(ierr);
	ierr = makeIntArray(&Avp_o_nnz, NULL, lnv); CHKERRQ(ierr);

	ierr = makeIntArray(&Apv_d_nnz, NULL, lnp); CHKERRQ(ierr);
	ierr = makeIntArray(&Apv_o_nnz, NULL, lnp); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// clear iterator (velocity matrices)
	iter = 0;

	//---------
	// X-points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// Avv
		nd = 1;
		no = 0;
		// x-velocity
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i-1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j+1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j-1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k+1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k-1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k][j][i-1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j+1][i-1]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j][i-1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j][i-1]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		Avv_d_nnz[iter] = nd;
		Avv_o_nnz[iter] = no;
		// Avp
		nd = 0;
		no = 0;
		// pressure
		ind = (PetscInt) ip[k][j][i-1];    CHECK_DOF(ind, startp, lnp, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, startp, lnp, nd, no);
		Avp_d_nnz[iter] = nd;
		Avp_o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// Avv
		nd = 1;
		no = 0;
		// y-velocity
		ind = (PetscInt) ivy[k][j][i+1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j][i-1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j-1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k+1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k-1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j-1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j-1][i]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
 		// x-velocity
		ind = (PetscInt) ivx[k][j-1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j-1][i+1]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, startv, lnv, nd, no);
		Avv_d_nnz[iter] = nd;
		Avv_o_nnz[iter] = no;
		// Avp
		nd = 0;
		no = 0;
		// pressure
		ind = (PetscInt) ip[k][j-1][i];    CHECK_DOF(ind, startp, lnp, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, startp, lnp, nd, no);
		Avp_d_nnz[iter] = nd;
		Avp_o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// Avv
		nd = 1;
		no = 0;
		// z-velocity
		ind = (PetscInt) ivz[k][j][i+1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k][j][i-1];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k][j+1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k][j-1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k-1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
 		// x-velocity
		ind = (PetscInt) ivx[k-1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k-1][j][i+1]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i+1];   CHECK_DOF(ind, startv, lnv, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k-1][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k-1][j+1][i]; CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j][i];     CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j+1][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		Avv_d_nnz[iter] = nd;
		Avv_o_nnz[iter] = no;
		// Avp
		nd = 0;
		no = 0;
		// pressure
		ind = (PetscInt) ip[k-1][j][i];    CHECK_DOF(ind, startp, lnp, nd, no);
		ind = (PetscInt) ip[k][j][i];      CHECK_DOF(ind, startp, lnp, nd, no);
		// update counters
		Avp_d_nnz[iter] = nd;
		Avp_o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	// clear iterator (pressure matrices)
	iter = 0;

	//---------
	// P-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// Apv
		nd = 0;
		no = 0;
		// x-velocity
		ind = (PetscInt) ivx[k][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivx[k][j][i+1]; CHECK_DOF(ind, startv, lnv, nd, no);
		// y-velocity
		ind = (PetscInt) ivy[k][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivy[k][j+1][i]; CHECK_DOF(ind, startv, lnv, nd, no);
		// z-velocity
		ind = (PetscInt) ivz[k][j][i];   CHECK_DOF(ind, startv, lnv, nd, no);
		ind = (PetscInt) ivz[k+1][j][i]; CHECK_DOF(ind, startv, lnv, nd, no);
		// update counters
		Apv_d_nnz[iter] = nd;
		Apv_o_nnz[iter] = no;
		iter++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// create matrices & vectors
	ierr = MatAIJCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &P->Avv);  CHKERRQ(ierr);
	ierr = MatAIJCreate(lnv, lnp, 0, Avp_d_nnz, 0, Avp_o_nnz, &P->Avp);  CHKERRQ(ierr);
	ierr = MatAIJCreate(lnp, lnv, 0, Apv_d_nnz, 0, Apv_o_nnz, &P->Apv);  CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(lnp, startp, &P->App);                       CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(lnp, startp, &P->iS);                        CHKERRQ(ierr);

	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &P->xv); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnp, PETSC_DETERMINE, &P->xp); CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->rv);                                  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xv, &P->wv);                                  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xp, &P->rp);                                  CHKERRQ(ierr);
	ierr = VecDuplicate(P->xp, &P->wp);                                  CHKERRQ(ierr);

	// free counter arrays
	ierr = PetscFree(Avv_d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(Avv_o_nnz); CHKERRQ(ierr);
	ierr = PetscFree(Avp_d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(Avp_o_nnz); CHKERRQ(ierr);
	ierr = PetscFree(Apv_d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(Apv_o_nnz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBlockAssemble"
PetscErrorCode PMatBlockAssemble(PMat pm)
{
	//======================================================================
	// (pgamma >= 1) - is a penalty parameter
	//
	// if pgamma > 1, Avv will contain velocity Schur complement:
	//
	//    Avv - Avp*(S^-1)*Apv
	//
	// where S is the pressure Schur complement preconditioner:
	//
	//    S = - 1/(K*dt) - 1/(pgamma*eta)
	//======================================================================

	JacRes      *jr;
	FDSTAG      *fs;
	BCCtx       *bc;
	DOFIndex    *dof;
	PMatBlock   *P;
	PetscInt    idx[7];
	PetscScalar v[49], a[36], d[6], g[6];
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, rho, IKdt, diag, pgamma, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	jr  = pm->jr;
	fs  = jr->fs;
	bc  = jr->bc;
	dof = &fs->dof;
	P   = (PMatBlock*)pm->data;

	// get density gradient stabilization parameters
	dt   = jr->ts->dt; // time step
	fssa = jr->FSSA;   // density gradient penalty parameter
    grav = jr->grav;   // gravity acceleration

	// get penalty parameter
	pgamma = pm->pgamma;

	// clear matrix coefficients
	ierr = MatZeroEntries(P->Avv); CHKERRQ(ierr);
	ierr = MatZeroEntries(P->Avp); CHKERRQ(ierr);
	ierr = MatZeroEntries(P->Apv); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);  CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,   &bcp);  CHKERRQ(ierr);

	//---------------
	// central points
	//---------------

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		eta  = jr->svCell[iter].svDev.eta;
		IKdt = jr->svCell[iter].svBulk.IKdt;
		rho  = jr->svCell[iter].svBulk.rho;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get pressure diagonal element (with penalty)
		diag = -IKdt -1.0/(pgamma*eta);

		// compute local matrix
		pm->getStiffMat(eta, diag, v, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

		// compute density gradient stabilization terms
		addDensGradStabil(fssa, v, rho, dt, grav, fdx, fdy, fdz, bdx, bdy, bdz);

		// compute velocity Schur complement
		if(pm->pgamma != 1.0) getVelSchur(v, d, g);

		// get global indices of the points:
		// vx_(i), vx_(i+1), vy_(j), vy_(j+1), vz_(k), vz_(k+1), p
		idx[0] = (PetscInt) ivx[k][j][i];
		idx[1] = (PetscInt) ivx[k][j][i+1];
		idx[2] = (PetscInt) ivy[k][j][i];
		idx[3] = (PetscInt) ivy[k][j+1][i];
		idx[4] = (PetscInt) ivz[k][j][i];
		idx[5] = (PetscInt) ivz[k+1][j][i];
		idx[6] = (PetscInt) ip[k][j][i];

		// get boundary constraints
		pdofidx[0] = -1;   cf[0] = bcvx[k][j][i];
		pdofidx[1] = -1;   cf[1] = bcvx[k][j][i+1];
		pdofidx[2] = -1;   cf[2] = bcvy[k][j][i];
		pdofidx[3] = -1;   cf[3] = bcvy[k][j+1][i];
		pdofidx[4] = -1;   cf[4] = bcvz[k][j][i];
		pdofidx[5] = -1;   cf[5] = bcvz[k+1][j][i];
		pdofidx[6] = -1;   cf[6] = bcp[k][j][i];

		// constrain local matrix
		constrLocalMat(7, pdofidx, cf, v);

		// extract operators
		getSubMat(v, a, d, g);

		// update global matrices
		ierr = MatSetValues(P->Avv, 6, idx,   6, idx,   a,    ADD_VALUES);    CHKERRQ(ierr);
		ierr = MatSetValues(P->Avp, 6, idx,   1, idx+6, g,    ADD_VALUES);    CHKERRQ(ierr);
		ierr = MatSetValues(P->Apv, 1, idx+6, 6, idx,   d,    ADD_VALUES);    CHKERRQ(ierr);
		ierr = MatSetValue (P->App, idx[6], idx[6], -IKdt,    INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue (P->iS,  idx[6], idx[6], 1.0/diag, INSERT_VALUES); CHKERRQ(ierr);

		// increment iterator
		iter++;
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svXYEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// compute local matrix
		//       vx_(j-1)             vx_(j)               vy_(i-1)             vy_(i)
		v[0]  =  eta/dy/bdy; v[1]  = -eta/dy/bdy; v[2]  =  eta/dx/bdy; v[3]  = -eta/dx/bdy; // fx_(j-1) [sxy]
		v[4]  = -eta/dy/fdy; v[5]  =  eta/dy/fdy; v[6]  = -eta/dx/fdy; v[7]  =  eta/dx/fdy; // fx_(j)   [sxy]
		v[8]  =  eta/dy/bdx; v[9]  = -eta/dy/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fy_(i-1) [sxy]
		v[12] = -eta/dy/fdx; v[13] =  eta/dy/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fy_(i)   [sxy]

		// get global indices of the points: vx_(j-1), vx_(j), vy_(i-1), vy_(i)
		idx[0] = (PetscInt) ivx[k][j-1][i];
		idx[1] = (PetscInt) ivx[k][j][i];
		idx[2] = (PetscInt) ivy[k][j][i-1];
		idx[3] = (PetscInt) ivy[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k][j-1][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvy[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvy[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svXZEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// compute local matrix
		//       vx_(k-1)             vx_(k)               vz_(i-1)             vz_(i)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dx/bdz; v[3]  = -eta/dx/bdz; // fx_(k-1) [sxz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dx/fdz; v[7]  =  eta/dx/fdz; // fx_(k)   [sxz]
		v[8]  =  eta/dz/bdx; v[9]  = -eta/dz/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fz_(i-1) [sxz]
		v[12] = -eta/dz/fdx; v[13] =  eta/dz/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fz_(i)   [sxz]

		// get global indices of the points: vx_(k-1), vx_(k), vz_(i-1), vz_(i)
		idx[0] = (PetscInt) ivx[k-1][j][i];
		idx[1] = (PetscInt) ivx[k][j][i];
		idx[2] = (PetscInt) ivz[k][j][i-1];
		idx[3] = (PetscInt) ivz[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get viscosity
		eta = jr->svYZEdge[iter++].svDev.eta;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// compute local matrix
		//       vy_(k-1)             vy_(k)               vz_(j-1)             vz_(j)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dy/bdz; v[3]  = -eta/dy/bdz; // fy_(k-1) [syz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dy/fdz; v[7]  =  eta/dy/fdz; // fy_(k)   [syz]
		v[8]  =  eta/dz/bdy; v[9]  = -eta/dz/bdy; v[10] =  eta/dy/bdy; v[11] = -eta/dy/bdy; // fz_(j-1) [syz]
		v[12] = -eta/dz/fdy; v[13] =  eta/dz/fdy; v[14] = -eta/dy/fdy; v[15] =  eta/dy/fdy; // fz_(j)   [syz]

		// get global indices of the points: vy_(k-1), vy_(k), vz_(j-1), vz_(j)
		idx[0] = (PetscInt) ivy[k-1][j][i];
		idx[1] = (PetscInt) ivy[k][j][i];
		idx[2] = (PetscInt) ivz[k][j-1][i];
		idx[3] = (PetscInt) ivz[k][j][i];

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvy[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvy[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j-1][i];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp);  CHKERRQ(ierr);

	// assemble velocity-pressure matrix blocks, remove constrained rows
	ierr = MatAIJAssemble(P->Avv, bc->vNumSPC, bc->vSPCList, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->Avp, bc->vNumSPC, bc->vSPCList, 0.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->Apv, bc->pNumSPC, bc->pSPCList, 0.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->App, bc->pNumSPC, bc->pSPCList, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->iS,  bc->pNumSPC, bc->pSPCList, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBlockPicardClean"
PetscErrorCode PMatBlockPicardClean(Mat J, Vec x, Vec r)
{
	//=======================================================================
	// Get action of the Picard Jacobian
	//
	// rv = Avv*xv + Avp*xp
	// rp = Apv*xv + App*xp
	//=======================================================================

	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	ierr = MatShellGetContext(J, (void**)&P); CHKERRQ(ierr);

	// extract solution blocks
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = MatMult(P->Apv, P->xv, P->rp); CHKERRQ(ierr); // rp = Apv*xv

	ierr = MatMult(P->App, P->xp, P->wp); CHKERRQ(ierr); // wp = App*xp

	ierr = VecAXPY(P->rp, 1.0, P->wp);    CHKERRQ(ierr); // rp = rp + wp

	ierr = MatMult(P->Avp, P->xp, P->rv); CHKERRQ(ierr); // rv = Avp*xp

	ierr = MatMult(P->Avv, P->xv, P->wv); CHKERRQ(ierr); // wv = Avv*xv

	ierr = VecAXPY(P->rv, 1.0, P->wv);    CHKERRQ(ierr); // rv = rv + wv

	// compose coupled residual
	ierr = VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBlockPicardSchur"
PetscErrorCode PMatBlockPicardSchur(Mat J, Vec x, Vec r)
{
	//=======================================================================
	// Get action of the Picard Jacobian using velocity Schur complement
	//
	// rv = Avv*xv + Avp*(xp + (S^-1)*Apv*xv)
	// rp = Apv*xv + App*xp
	//=======================================================================

	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	ierr = MatShellGetContext(J, (void**)&P); CHKERRQ(ierr);

	// extract solution blocks
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = MatMult(P->Apv, P->xv, P->rp); CHKERRQ(ierr); // rp = Apv*xv

	ierr = MatMult(P->iS, P->rp, P->wp);  CHKERRQ(ierr); // wp = (S^-1)*rp

	ierr = VecAXPY(P->wp, 1.0, P->xp);    CHKERRQ(ierr); // wp = wp + xp

	ierr = MatMult(P->Avp, P->wp, P->rv); CHKERRQ(ierr); // rv = Avp*wp

	ierr = MatMult(P->App, P->xp, P->wp); CHKERRQ(ierr); // wp = App*xp

	ierr = VecAXPY(P->rp, 1.0, P->wp);    CHKERRQ(ierr); // rp = rp + wp

	ierr = MatMult(P->Avv, P->xv, P->wv); CHKERRQ(ierr); // wv = Avv*xv

	ierr = VecAXPY(P->rv, 1.0, P->wv);    CHKERRQ(ierr); // rv = rv + wv

	// compose coupled residual
	ierr = VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBlockDestroy"
PetscErrorCode PMatBlockDestroy(PMat pm)
{
	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	P = (PMatBlock*)pm->data;

	ierr = MatDestroy (&P->Avv); CHKERRQ(ierr);
	ierr = MatDestroy (&P->Avp); CHKERRQ(ierr);
	ierr = MatDestroy (&P->Apv); CHKERRQ(ierr);
	ierr = MatDestroy (&P->App); CHKERRQ(ierr);
	ierr = MatDestroy (&P->iS);   CHKERRQ(ierr);
	ierr = VecDestroy (&P->rv);  CHKERRQ(ierr);
	ierr = VecDestroy (&P->rp);  CHKERRQ(ierr);
	ierr = VecDestroy (&P->xv);  CHKERRQ(ierr);
	ierr = VecDestroy (&P->xp);  CHKERRQ(ierr);
	ierr = VecDestroy (&P->wv);  CHKERRQ(ierr);
	ierr = VecDestroy (&P->wp);  CHKERRQ(ierr);
	ierr = PetscFree(P);         CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// SERVICE FUNCTIONS
//---------------------------------------------------------------------------
void getStiffMatDevProj(
	PetscScalar eta, PetscScalar diag,PetscScalar *v,
	PetscScalar dx,  PetscScalar dy,  PetscScalar dz,
	PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
	PetscScalar bdx, PetscScalar bdy, PetscScalar bdz)
{
	// compute cell stiffness matrix with deviatoric projection

	PetscScalar E43 = 4.0*eta/3.0;
	PetscScalar E23 = 2.0*eta/3.0;

	//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
	v[0]  =  E43/dx/bdx; v[1]  = -E43/dx/bdx; v[2]  = -E23/dy/bdx; v[3]  =  E23/dy/bdx; v[4]  = -E23/dz/bdx; v[5]  =  E23/dz/bdx; v[6]  =  1.0/bdx; // fx_(i)   [sxx]
	v[7]  = -E43/dx/fdx; v[8]  =  E43/dx/fdx; v[9]  =  E23/dy/fdx; v[10] = -E23/dy/fdx; v[11] =  E23/dz/fdx; v[12] = -E23/dz/fdx; v[13] = -1.0/fdx; // fx_(i+1) [sxx]
	v[14] = -E23/dx/bdy; v[15] =  E23/dx/bdy; v[16] =  E43/dy/bdy; v[17] = -E43/dy/bdy; v[18] = -E23/dz/bdy; v[19] =  E23/dz/bdy; v[20] =  1.0/bdy; // fy_(j)   [syy]
	v[21] =  E23/dx/fdy; v[22] = -E23/dx/fdy; v[23] = -E43/dy/fdy; v[24] =  E43/dy/fdy; v[25] =  E23/dz/fdy; v[26] = -E23/dz/fdy; v[27] = -1.0/fdy; // fy_(j+1) [syy]
	v[28] = -E23/dx/bdz; v[29] =  E23/dx/bdz; v[30] = -E23/dy/bdz; v[31] =  E23/dy/bdz; v[32] =  E43/dz/bdz; v[33] = -E43/dz/bdz; v[34] =  1.0/bdz; // fz_(k)   [szz]
	v[35] =  E23/dx/fdz; v[36] = -E23/dx/fdz; v[37] =  E23/dy/fdz; v[38] = -E23/dy/fdz; v[39] = -E43/dz/fdz; v[40] =  E43/dz/fdz; v[41] = -1.0/fdz; // fz_(k+1) [szz]
	v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  diag;    // g
}
//---------------------------------------------------------------------------
void getStiffMatClean(
	PetscScalar eta, PetscScalar diag,PetscScalar *v,
	PetscScalar dx,  PetscScalar dy,  PetscScalar dz,
	PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
	PetscScalar bdx, PetscScalar bdy, PetscScalar bdz)
{
	// compute cell stiffness matrix without deviatoric projection

	PetscScalar E2 = 2.0*eta;

	//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
	v[0]  =  E2/dx/bdx;  v[1]  = -E2/dx/bdx;  v[2]  =  0.0;        v[3]  =  0.0;        v[4]  =  0.0;        v[5]  =  0.0;        v[6]  =  1.0/bdx; // fx_(i)   [sxx]
	v[7]  = -E2/dx/fdx;  v[8]  =  E2/dx/fdx;  v[9]  =  0.0;        v[10] =  0.0;        v[11] =  0.0;        v[12] =  0.0;        v[13] = -1.0/fdx; // fx_(i+1) [sxx]
	v[14] =  0.0;        v[15] =  0.0;        v[16] =  E2/dy/bdy;  v[17] = -E2/dy/bdy;  v[18] =  0.0;        v[19] =  0.0;        v[20] =  1.0/bdy; // fy_(j)   [syy]
	v[21] =  0.0;        v[22] =  0.0;        v[23] = -E2/dy/fdy;  v[24] =  E2/dy/fdy;  v[25] =  0.0;        v[26] =  0.0;        v[27] = -1.0/fdy; // fy_(j+1) [syy]
	v[28] =  0.0;        v[29] =  0.0;        v[30] =  0.0;        v[31] =  0.0;        v[32] =  E2/dz/bdz;  v[33] = -E2/dz/bdz;  v[34] =  1.0/bdz; // fz_(k)   [szz]
	v[35] =  0.0;        v[36] =  0.0;        v[37] =  0.0;        v[38] =  0.0;        v[39] = -E2/dz/fdz;  v[40] =  E2/dz/fdz;  v[41] = -1.0/fdz; // fz_(k+1) [szz]
	v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  diag;    // g
}
//---------------------------------------------------------------------------
void addDensGradStabil(
	PetscScalar fssa, PetscScalar *v,
	PetscScalar rho,  PetscScalar dt,   PetscScalar *grav,
	PetscScalar fdx,  PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx,  PetscScalar bdy,  PetscScalar bdz)
{
	PetscScalar cf = -fssa*dt;

	// add stabilization terms
	v[0 ] -= cf*(rho*grav[0])/bdx;
	v[8 ] += cf*(rho*grav[0])/fdx;
	v[16] -= cf*(rho*grav[1])/bdy;
	v[24] += cf*(rho*grav[1])/fdy;
	v[32] -= cf*(rho*grav[2])/bdz;
	v[40] += cf*(rho*grav[2])/fdz;
/*
	fx[k][j][i]   -= vx[k][j][i]  *tx/bdx
	fx[k][j][i+1] += vx[k][j][i+1]*tx/fdx
	fy[k][j][i]   -= vy[k][j][i]  *ty/bdy
	fy[k][j+1][i] += vy[k][j+1][i]*ty/fdy
	fz[k][j][i]   -= vz[k][j][i]  *tz/bdz
	fz[k+1][j][i] += vz[k+1][j][i]*tz/fdz
*/
}
//---------------------------------------------------------------------------
void getVelSchur(PetscScalar v[], PetscScalar d[], PetscScalar g[])
{
	PetscScalar k;

	// extract divergence operator
	d[0] = v[42]; d[1] = v[43]; d[2] = v[44]; d[3] = v[45]; d[4] = v[46]; d[5] = v[47];

	// extract gradient operator
	g[0] = v[6];  g[1] = v[13]; g[2] = v[20]; g[3] = v[27]; g[4] = v[34]; g[5] = v[41];

	// get penalty term
	k = -1.0/v[48];

	// compute velocity Schur complement
	v[0]  += k*g[0]*d[0]; v[1]  += k*g[0]*d[1]; v[2]  += k*g[0]*d[2]; v[3]  += k*g[0]*d[3]; v[4]  += k*g[0]*d[4]; v[5]  += k*g[0]*d[5];
	v[7]  += k*g[1]*d[0]; v[8]  += k*g[1]*d[1]; v[9]  += k*g[1]*d[2]; v[10] += k*g[1]*d[3]; v[11] += k*g[1]*d[4]; v[12] += k*g[1]*d[5];
	v[14] += k*g[2]*d[0]; v[15] += k*g[2]*d[1]; v[16] += k*g[2]*d[2]; v[17] += k*g[2]*d[3]; v[18] += k*g[2]*d[4]; v[19] += k*g[2]*d[5];
	v[21] += k*g[3]*d[0]; v[22] += k*g[3]*d[1]; v[23] += k*g[3]*d[2]; v[24] += k*g[3]*d[3]; v[25] += k*g[3]*d[4]; v[26] += k*g[3]*d[5];
	v[28] += k*g[4]*d[0]; v[29] += k*g[4]*d[1]; v[30] += k*g[4]*d[2]; v[31] += k*g[4]*d[3]; v[32] += k*g[4]*d[4]; v[33] += k*g[4]*d[5];
	v[35] += k*g[5]*d[0]; v[36] += k*g[5]*d[1]; v[37] += k*g[5]*d[2]; v[38] += k*g[5]*d[3]; v[39] += k*g[5]*d[4]; v[40] += k*g[5]*d[5];
}
//---------------------------------------------------------------------------
void getSubMat(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[])
{
	// extract divergence operator
	d[0] = v[42]; d[1] = v[43]; d[2] = v[44]; d[3] = v[45]; d[4] = v[46]; d[5] = v[47];

	// extract gradient operator
	g[0] = v[6];  g[1] = v[13]; g[2] = v[20]; g[3] = v[27]; g[4] = v[34]; g[5] = v[41];

	// extract velocity operator
	a[0]  = v[0];  a[1]  = v[1];  a[2]  = v[2];  a[3]  = v[3];  a[4]  = v[4];  a[5]  = v[5];
	a[6]  = v[7];  a[7]  = v[8];  a[8]  = v[9];  a[9]  = v[10]; a[10] = v[11]; a[11] = v[12];
	a[12] = v[14]; a[13] = v[15]; a[14] = v[16]; a[15] = v[17]; a[16] = v[18]; a[17] = v[19];
	a[18] = v[21]; a[19] = v[22]; a[20] = v[23]; a[21] = v[24]; a[22] = v[25]; a[23] = v[26];
	a[24] = v[28]; a[25] = v[29]; a[26] = v[30]; a[27] = v[31]; a[28] = v[32]; a[29] = v[33];
	a[30] = v[35]; a[31] = v[36]; a[32] = v[37]; a[33] = v[38]; a[34] = v[39]; a[35] = v[40];
}
//---------------------------------------------------------------------------
/*
	// Pressure constraints implementation

	PetscScalar  cbot, ctop;
	PetscInt     mcz;
	cbot = 1.0; if(k == 0) 	 cbot = 2.0;
	ctop = 1.0; if(k == mcz) ctop = 2.0;
	//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
	v[0]  =  E43/dx/bdx; v[1]  = -E43/dx/bdx; v[2]  = -E23/dy/bdx; v[3]  =  E23/dy/bdx; v[4]  = -E23/dz/bdx; v[5]  =  E23/dz/bdx; v[6]  =  1.0/bdx;  // fx_(i)   [sxx]
	v[7]  = -E43/dx/fdx; v[8]  =  E43/dx/fdx; v[9]  =  E23/dy/fdx; v[10] = -E23/dy/fdx; v[11] =  E23/dz/fdx; v[12] = -E23/dz/fdx; v[13] = -1.0/fdx;  // fx_(i+1) [sxx]
	v[14] = -E23/dx/bdy; v[15] =  E23/dx/bdy; v[16] =  E43/dy/bdy; v[17] = -E43/dy/bdy; v[18] = -E23/dz/bdy; v[19] =  E23/dz/bdy; v[20] =  1.0/bdy;  // fy_(j)   [syy]
	v[21] =  E23/dx/fdy; v[22] = -E23/dx/fdy; v[23] = -E43/dy/fdy; v[24] =  E43/dy/fdy; v[25] =  E23/dz/fdy; v[26] = -E23/dz/fdy; v[27] = -1.0/fdy;  // fy_(j+1) [syy]
	v[28] = -E23/dx/bdz; v[29] =  E23/dx/bdz; v[30] = -E23/dy/bdz; v[31] =  E23/dy/bdz; v[32] =  E43/dz/bdz; v[33] = -E43/dz/bdz; v[34] =  cbot/bdz; // fz_(k)   [szz]
	v[35] =  E23/dx/fdz; v[36] = -E23/dx/fdz; v[37] =  E23/dy/fdz; v[38] = -E23/dy/fdz; v[39] = -E43/dz/fdz; v[40] =  E43/dz/fdz; v[41] = -ctop/fdz; // fz_(k+1) [szz]
	v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  0.0;      // g
*/
//---------------------------------------------------------------------------
void getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[])
{
	// apply two-point constraints on the ghost nodes
	PetscInt j;

	for(j = 0; j < n; j++)
	{
		// detect ghost point (the one that has negative global index)
		if(idx[j] == -1)
		{
			// ghost point detected, check whether potential primary DOF is constrained
			if(cf[pdofidx[j]] != DBL_MAX)
			{
				// primary DOF is constrained, put single-point constraint on the ghost node
				cf[j]      =  0.0;
				pdofidx[j] = -1;
			}
			else
			{
				// primary DOF is unconstrained (active), check constraint type
				if(cf[j] != DBL_MAX) cf[j] = -1.0; // no-slip (arbitrary)
				else                 cf[j] =  1.0; // free-slip
			}
		}
		else
		{	// internal point detected, cannot have two-point constraints
			pdofidx[j] = -1;
		}
	}
}
//---------------------------------------------------------------------------
void constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[])
{
	//=========================================================================
	// Apply linear constraints to a local matrix
	//
	// Function supports:
	//  - single-point constraints (x = c)
	//  - two-point constraints    (x = cf*y + c)
	// where:
	//    x   - constrained DOF
	//    y   - primary DOF
	//    cf  - linear combination parameter
	//    c   - constant prescribed value (may be zero for homogeneous case)
	//
	// For the two-point constraints it must be guaranteed that the primary
	// DOF is an active (unconstrained) DOF. Otherwise the two-point constraint
	// simply degenerates to a single-point constraint.
	//
	// Active (unconstrained) DOF are marked by DBL_MAX constant:
	//
	// cf[j] != DBL_MAX:
	//    j is a constrained DOF
	//    cf[j]      - linear combination parameter (not used by single-point constraint)
	//    pdofidx[j] - relative index of the primary DOF (-1 for single-point constraint)
	//
	// otherwise:
	//    j is an unconstrained (active) DOF
	//=========================================================================

	PetscInt i, j, jj, jst;

	for(i = 0, jj = 0; i < n; i++)
	{
		// skip constrained rows (handled by MatZeroRows after assembly)
		if(cf[i] != DBL_MAX) { jj += n; continue; }

		// store pointer to row beginning
		jst = jj;

		// proceed with unconstrained row, check columns
		for(j = 0; j < n; j++, jj++)
		{
			// detect constrained columns
			if(cf[j] != DBL_MAX)
			{
				// compute linear combination of primary & constrained columns
				// NOTE: two-point constraints only!
				if(pdofidx[j] != -1)
				{
					v[jst + pdofidx[j]] += cf[j]*v[jj];
				}

				// zero out constrained column
				v[jj] = 0.0;
			}
		}
	}
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecScatterBlockToMonolithic"
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode)
{
	// scatter block vectors to monolithic format forward & reverse

	PetscInt     fs,  gs,  bs;
	PetscScalar *fp, *gp, *bp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get sizes of the blocks
	ierr = VecGetLocalSize(f, &fs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(g, &gs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(b, &bs); CHKERRQ(ierr);

	if(bs != fs+gs)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Block sizes don't match monolithic format");
	}

	// access vectors
	ierr = VecGetArray(f, &fp); CHKERRQ(ierr);
	ierr = VecGetArray(g, &gp); CHKERRQ(ierr);
	ierr = VecGetArray(b, &bp); CHKERRQ(ierr);

	if(mode == SCATTER_FORWARD)
	{
		// block-to-monolithic
		ierr = PetscMemcpy(bp,    fp, (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(bp+fs, gp, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);
	}
	if(mode == SCATTER_REVERSE)
	{
		// monolithic-to-block
		ierr = PetscMemcpy(fp, bp,    (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(gp, bp+fs, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);
	}

	// restore access
	ierr = VecRestoreArray(f, &fp); CHKERRQ(ierr);
	ierr = VecRestoreArray(g, &gp); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &bp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
