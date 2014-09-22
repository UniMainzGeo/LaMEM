//---------------------------------------------------------------------------
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "matrix.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatCreate"
PetscErrorCode PMatCreate(
	PetscInt m, PetscInt n,
	PetscInt d_nz, const PetscInt d_nnz[],
	PetscInt o_nz, const PetscInt o_nnz[],
	Mat *P)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create matrix
	ierr = MatCreate(PETSC_COMM_WORLD, P); CHKERRQ(ierr);
	ierr = MatSetType((*P), MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatSetSizes((*P), m, n, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);

	// preallocate matrix
	ierr = MatSeqAIJSetPreallocation((*P), d_nz, d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation((*P), d_nz, d_nnz, o_nz, o_nnz); CHKERRQ(ierr);

	// throw an error if preallocation fails
	ierr = MatSetOption((*P), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetUp((*P)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatAssemble"
PetscErrorCode PMatAssemble(Mat P, PetscInt numRows, const PetscInt rows[])
{
	PetscInt    m, n;
	PetscScalar diag;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// freeze nonzero structure, constrain rows only locally
	ierr = MatSetOption(P, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);     CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);    CHKERRQ(ierr);

	// zero out constrained rows, form unit diagonal for the constrained block
	if(numRows)
	{	ierr = MatGetSize(P, &m, &n); CHKERRQ(ierr);
		if(m == n) diag = 1.0;
		else       diag = 0.0;
		ierr = MatZeroRows(P, numRows, rows, diag, NULL, NULL); CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BMatCreate"
PetscErrorCode BMatCreate(BMat *bmat,
	PetscInt  lnv,       PetscInt  lnp,
	PetscInt *Avv_d_nnz, PetscInt *Avv_o_nnz,
	PetscInt *Avp_d_nnz, PetscInt *Avp_o_nnz,
	PetscInt *Apv_d_nnz, PetscInt *Apv_o_nnz)
{
	Mat M[4];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create velocity-pressure matrix blocks
	ierr = PMatCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &bmat->Avv); CHKERRQ(ierr);
	ierr = PMatCreate(lnv, lnp, 0, Avp_d_nnz, 0, Avp_o_nnz, &bmat->Avp); CHKERRQ(ierr);
	ierr = PMatCreate(lnp, lnv, 0, Apv_d_nnz, 0, Apv_o_nnz, &bmat->Apv); CHKERRQ(ierr);
	ierr = PMatCreate(lnp, lnp, 1, NULL,      0, NULL,      &bmat->App); CHKERRQ(ierr);

	// setup matrix array
	M[0] = bmat->Avv;
	M[1] = bmat->Avp;
	M[2] = bmat->Apv;
	M[3] = bmat->App;

	// create nested matrix
	ierr = MatCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, 2, PETSC_NULL, (const Mat*)M, &bmat->A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BMatDestroy"
PetscErrorCode BMatDestroy(BMat *bmat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatDestroy(&bmat->Avv); CHKERRQ(ierr);
	ierr = MatDestroy(&bmat->Avp); CHKERRQ(ierr);
	ierr = MatDestroy(&bmat->Apv); CHKERRQ(ierr);
	ierr = MatDestroy(&bmat->App); CHKERRQ(ierr);
	ierr = MatDestroy(&bmat->A);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BMatClearSubMat"
PetscErrorCode BMatClearSubMat(BMat *bmat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// zero all matrices
	ierr = MatZeroEntries(bmat->Avv); CHKERRQ(ierr);
	ierr = MatZeroEntries(bmat->Avp); CHKERRQ(ierr);
	ierr = MatZeroEntries(bmat->Apv); CHKERRQ(ierr);
	ierr = MatZeroEntries(bmat->App); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BMatAssemble"
PetscErrorCode BMatAssemble(BMat *bmat, BCCtx *bc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// assemble velocity-pressure matrix blocks, remove constrained rows
	ierr = PMatAssemble(bmat->Avv, bc->numSPC,     bc->SPCList);     CHKERRQ(ierr);
	ierr = PMatAssemble(bmat->Avp, bc->numSPC,     bc->SPCList);     CHKERRQ(ierr);
	ierr = PMatAssemble(bmat->Apv, bc->numSPCPres, bc->SPCListPres); CHKERRQ(ierr);
	ierr = PMatAssemble(bmat->App, bc->numSPCPres, bc->SPCListPres); CHKERRQ(ierr);

	// assemble block matrix
	ierr = MatAssemblyBegin(bmat->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (bmat->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatCreateMonolithic"
PetscErrorCode PMatCreateMonolithic(
	FDSTAG *fs,
	Mat    *P,
	Mat    *M)
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

	PetscInt    ln, start, ind, nd, no, *d_nnz, *o_nnz;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	DOFIndex    *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	dof = &fs->cdof;

	// get number of local rows & global index of the first row
	start = dof->istart;
	ln    = dof->numdof;

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

	// create velocity-pressure matrix
	ierr = PMatCreate(ln, ln, 0, d_nnz, 0, o_nnz, P); CHKERRQ(ierr);
	ierr = PMatCreate(ln, ln, 1, NULL,  0, NULL,  M); CHKERRQ(ierr);

	// clear work arrays
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatAssembleMonolithic"
PetscErrorCode PMatAssembleMonolithic(
	JacRes  *jr,
	Mat      P,
	Mat      M)
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

	FDSTAG     *fs;
	BCCtx      *bc;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, IKdt, E43, E23;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];
	PetscScalar diag;
	DOFIndex    *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  = jr->fs;
	bc  = jr->cbc;   // coupled
	dof = &fs->cdof; // coupled

	// clear matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);
	ierr = MatZeroEntries(M); CHKERRQ(ierr);

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
		// get shear & inverse bulk viscosities
		eta  = jr->svCell[iter].svDev.eta;
		IKdt = jr->svCell[iter].svBulk.IKdt;
		iter++;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute local matrix
		E43 = 4.0*eta/3.0;
		E23 = 2.0*eta/3.0;

		// get diagonal element
		diag = -IKdt -1.0/eta;

		//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
		v[0]  =  E43/dx/bdx; v[1]  = -E43/dx/bdx; v[2]  = -E23/dy/bdx; v[3]  =  E23/dy/bdx; v[4]  = -E23/dz/bdx; v[5]  =  E23/dz/bdx; v[6]  =  1.0/bdx; // fx_(i)   [sxx]
		v[7]  = -E43/dx/fdx; v[8]  =  E43/dx/fdx; v[9]  =  E23/dy/fdx; v[10] = -E23/dy/fdx; v[11] =  E23/dz/fdx; v[12] = -E23/dz/fdx; v[13] = -1.0/fdx; // fx_(i+1) [sxx]
		v[14] = -E23/dx/bdy; v[15] =  E23/dx/bdy; v[16] =  E43/dy/bdy; v[17] = -E43/dy/bdy; v[18] = -E23/dz/bdy; v[19] =  E23/dz/bdy; v[20] =  1.0/bdy; // fy_(j)   [syy]
		v[21] =  E23/dx/fdy; v[22] = -E23/dx/fdy; v[23] = -E43/dy/fdy; v[24] =  E43/dy/fdy; v[25] =  E23/dz/fdy; v[26] = -E23/dz/fdy; v[27] = -1.0/fdy; // fy_(j+1) [syy]
		v[28] = -E23/dx/bdz; v[29] =  E23/dx/bdz; v[30] = -E23/dy/bdz; v[31] =  E23/dy/bdz; v[32] =  E43/dz/bdz; v[33] = -E43/dz/bdz; v[34] =  1.0/bdz; // fz_(k)   [szz]
		v[35] =  E23/dx/fdz; v[36] = -E23/dx/fdz; v[37] =  E23/dy/fdz; v[38] = -E23/dy/fdz; v[39] = -E43/dz/fdz; v[40] =  E43/dz/fdz; v[41] = -1.0/fdz; // fz_(k+1) [szz]
		v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  diag;    // g

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
		ierr = constrLocalMat(7, pdofidx, cf, v); CHKERRQ(ierr);

		// add to global matrix
		ierr = MatSetValues(P, 7, idx, 7, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// store inverse viscosity for compensation
		ierr = MatSetValue(M, idx[6], idx[6], -1.0/eta, INSERT_VALUES); CHKERRQ(ierr);
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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

		// add to global matrix
		ierr = MatSetValues(P, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

		// add to global matrix
		ierr = MatSetValues(P, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

		// add to global matrix
		ierr = MatSetValues(P, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx);       CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy);       CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz);       CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);        CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp);   CHKERRQ(ierr);

	// assemble velocity-pressure matrix, remove constrained rows
	ierr = PMatAssemble(P, bc->numSPC, bc->SPCList); CHKERRQ(ierr);
	ierr = PMatAssemble(M, bc->numSPC, bc->SPCList); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatCreateBlock"
PetscErrorCode PMatCreateBlock(
	FDSTAG *fs,
	BMat   *P,
	Mat    *M)
{
	PetscInt    lnp, lnv, startv, startp, nd, no, ind;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    *Avv_d_nnz, *Avv_o_nnz;
	PetscInt    *Avp_d_nnz, *Avp_o_nnz;
	PetscInt    *Apv_d_nnz, *Apv_o_nnz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	DOFIndex    *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	dof = &fs->udof; // uncoupled

	// get number of local rows & global index of the first row
	lnv    = dof->numdof;
	startv = dof->istart;
	lnp    = dof->numdofp;
	startp = dof->istartp;

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

	// create block matrix
	ierr = BMatCreate(
		P, lnv, lnp,
		Avv_d_nnz, Avv_o_nnz,
		Avp_d_nnz, Avp_o_nnz,
		Apv_d_nnz, Apv_o_nnz); CHKERRQ(ierr);

	ierr = PMatCreate(lnp, lnp, 1, NULL,  0, NULL, M); CHKERRQ(ierr);

	// clear counter arrays
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
#define __FUNCT__ "PMatAssembleBlock"
PetscErrorCode PMatAssembleBlock(
	JacRes      *jr,
	BMat        *P,
	Mat          M,
	PetscScalar  pgamma)
{
	//======================================================================
	// pgamma - is a penalty parameter
	//
	// if pgamma is nonzero:
	// M   - will contain -kappa (inverse penalty matrix)
	// Avv - will contain Avv + kappa*Avp*Apv (velocity Schur complement)
	// kappa = 1/(1/(K*dt) + 1/(pgamma*eta)
	//
	// otherwise:
	// M will contain -kappa (pressure Schur complement preconditioner)
	// Avv - will contain unmodified velocity operator
	// kappa = 1/(K*dt) + 1/eta
	//======================================================================
	FDSTAG      *fs;
	BCCtx       *bc;
	PetscInt    idx[7];
	PetscScalar v[49], a[36], d[6], g[6];
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, IKdt, E43, E23;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar kappa;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];
	DOFIndex    *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  = jr->fs;
	bc  = jr->ubc;   // uncoupled
	dof = &fs->udof; // uncoupled

	// clear matrix coefficients
	ierr = BMatClearSubMat(P); CHKERRQ(ierr);
	ierr = MatZeroEntries (M); CHKERRQ(ierr);

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
		// get shear & inverse bulk viscosities
		eta  = jr->svCell[iter].svDev.eta;
		IKdt = jr->svCell[iter].svBulk.IKdt;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute local matrix
		E43 = 4.0*eta/3.0;
		E23 = 2.0*eta/3.0;

		//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
		v[0]  =  E43/dx/bdx; v[1]  = -E43/dx/bdx; v[2]  = -E23/dy/bdx; v[3]  =  E23/dy/bdx; v[4]  = -E23/dz/bdx; v[5]  =  E23/dz/bdx; v[6]  =  1.0/bdx; // fx_(i)   [sxx]
		v[7]  = -E43/dx/fdx; v[8]  =  E43/dx/fdx; v[9]  =  E23/dy/fdx; v[10] = -E23/dy/fdx; v[11] =  E23/dz/fdx; v[12] = -E23/dz/fdx; v[13] = -1.0/fdx; // fx_(i+1) [sxx]
		v[14] = -E23/dx/bdy; v[15] =  E23/dx/bdy; v[16] =  E43/dy/bdy; v[17] = -E43/dy/bdy; v[18] = -E23/dz/bdy; v[19] =  E23/dz/bdy; v[20] =  1.0/bdy; // fy_(j)   [syy]
		v[21] =  E23/dx/fdy; v[22] = -E23/dx/fdy; v[23] = -E43/dy/fdy; v[24] =  E43/dy/fdy; v[25] =  E23/dz/fdy; v[26] = -E23/dz/fdy; v[27] = -1.0/fdy; // fy_(j+1) [syy]
		v[28] = -E23/dx/bdz; v[29] =  E23/dx/bdz; v[30] = -E23/dy/bdz; v[31] =  E23/dy/bdz; v[32] =  E43/dz/bdz; v[33] = -E43/dz/bdz; v[34] =  1.0/bdz; // fz_(k)   [szz]
		v[35] =  E23/dx/fdz; v[36] = -E23/dx/fdz; v[37] =  E23/dy/fdz; v[38] = -E23/dy/fdz; v[39] = -E43/dz/fdz; v[40] =  E43/dz/fdz; v[41] = -1.0/fdz; // fz_(k+1) [szz]
		v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  0.0;     // g

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
		ierr = constrLocalMat(7, pdofidx, cf, v); CHKERRQ(ierr);

		// extract operators, compute penalty terms and preconditioners
		if(pgamma)
		{
			// compute penalty parameter
			kappa = 1.0/(IKdt + 1.0/(pgamma*eta));

			// get velocity Schur complement
			ierr = getVelSchurComp(v, a, d, g, kappa); CHKERRQ(ierr);
		}
		else
		{	// compute pressure Schur complement preconditioner
			kappa = IKdt + 1.0/eta;

			// get sub-matrices
			ierr = getSubMats(v, a, d, g); CHKERRQ(ierr);
		}

		// update global matrices
		ierr = MatSetValues(P->Avv, 6, idx,   6, idx,   a, ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(P->Avp, 6, idx,   1, idx+6, g, ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(P->Apv, 1, idx+6, 6, idx,   d, ADD_VALUES); CHKERRQ(ierr);

		// pressure diagonal blocks
		ierr = MatSetValue (P->App, idx[6], idx[6], -IKdt,  INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue (M,      idx[6], idx[6], -kappa, INSERT_VALUES); CHKERRQ(ierr);

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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

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
		ierr = getTwoPointConstr(4, idx, pdofidx, cf); CHKERRQ(ierr);

		// constrain local matrix
		ierr = constrLocalMat(4, pdofidx, cf, v); CHKERRQ(ierr);

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

	// assemble block matrix
	ierr = BMatAssemble(P, bc);                              CHKERRQ(ierr);
	ierr = PMatAssemble(M, bc->numSPCPres, bc->SPCListPres); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getTwoPointConstr"
PetscErrorCode getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[])
{
	// apply two-point constraints on the ghost nodes
	PetscInt j;

	PetscFunctionBegin;

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
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "constrLocalMat"
PetscErrorCode constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[])
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

	PetscFunctionBegin;

	for(i = 0, jj = 0; i < n; i++)
	{
		// detect constrained rows
		// NOTE: constrained rows are handled by MatZeroRows after assembly
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
				if(pdofidx[j] != -1) v[jst + pdofidx[j]] += cf[j]*v[jj];

				// zero out constrained column
				v[jj] = 0.0;
			}
		}
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getVelSchurComp"
PetscErrorCode getVelSchurComp(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[], PetscScalar k)
{
	PetscFunctionBegin;

	// extract divergence operator
	d[0] = v[42]; d[1] = v[43]; d[2] = v[44]; d[3] = v[45]; d[4] = v[46]; d[5] = v[47];

	// extract gradient operator
	g[0] = v[6];  g[1] = v[13]; g[2] = v[20]; g[3] = v[27]; g[4] = v[34]; g[5] = v[41];

	// compute & extract velocity Schur complement
	a[0]  = v[0]  + k*g[0]*d[0]; a[1]  = v[1]  + k*g[0]*d[1]; a[2]  = v[2]  + k*g[0]*d[2]; a[3]  = v[3]  + k*g[0]*d[3]; a[4]  = v[4]  + k*g[0]*d[4]; a[5]  = v[5]  + k*g[0]*d[5];
	a[6]  = v[7]  + k*g[1]*d[0]; a[7]  = v[8]  + k*g[1]*d[1]; a[8]  = v[9]  + k*g[1]*d[2]; a[9]  = v[10] + k*g[1]*d[3]; a[10] = v[11] + k*g[1]*d[4]; a[11] = v[12] + k*g[1]*d[5];
	a[12] = v[14] + k*g[2]*d[0]; a[13] = v[15] + k*g[2]*d[1]; a[14] = v[16] + k*g[2]*d[2]; a[15] = v[17] + k*g[2]*d[3]; a[16] = v[18] + k*g[2]*d[4]; a[17] = v[19] + k*g[2]*d[5];
	a[18] = v[21] + k*g[3]*d[0]; a[19] = v[22] + k*g[3]*d[1]; a[20] = v[23] + k*g[3]*d[2]; a[21] = v[24] + k*g[3]*d[3]; a[22] = v[25] + k*g[3]*d[4]; a[23] = v[26] + k*g[3]*d[5];
	a[24] = v[28] + k*g[4]*d[0]; a[25] = v[29] + k*g[4]*d[1]; a[26] = v[30] + k*g[4]*d[2]; a[27] = v[31] + k*g[4]*d[3]; a[28] = v[32] + k*g[4]*d[4]; a[29] = v[33] + k*g[4]*d[5];
	a[30] = v[35] + k*g[5]*d[0]; a[31] = v[36] + k*g[5]*d[1]; a[32] = v[37] + k*g[5]*d[2]; a[33] = v[38] + k*g[5]*d[3]; a[34] = v[39] + k*g[5]*d[4]; a[35] = v[40] + k*g[5]*d[5];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getSubMats"
PetscErrorCode getSubMats(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[])
{
	PetscFunctionBegin;

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

	PetscFunctionReturn(0);
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
