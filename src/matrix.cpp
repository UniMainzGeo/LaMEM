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
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "matData.h"
#include "matrix.h"
#include "tssolve.h"
#include "fdstag.h"
#include "tools.h"

// WARNING!!! rewrite everything here using submatrix functions

//---------------------------------------------------------------------------
PetscErrorCode PMatCreate(MatData *md, Mat *A)
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
	PetscInt    ln, start, ind, nd, no, *d_nnz, *o_nnz;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs  = md->fs;
	dof = &fs->dof;

	// get number of local rows & global index of the first row
	ln    = dof->ln;
	start = dof->st;

	// allocate nonzero counter arrays
	ierr = makeIntArray(&d_nnz, NULL, ln); CHKERRQ(ierr);
	ierr = makeIntArray(&o_nnz, NULL, ln); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// clear iterator
	iter = 0;

	//---------
	// X-points
	//---------
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// create matrix
	ierr = MatAIJCreate(ln, ln, 0, d_nnz, 0, o_nnz, A); CHKERRQ(ierr);

	// attach near null space
	ierr = MatAIJSetNullSpace((*A), md); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatAssemble(MatData *md, PetscScalar pgamma, Mat A)
{
	//======================================================================
	// Assemble effective viscosity preconditioning matrix
	//
	// Global ordering of the variables is interlaced:
	// all X-Y-Z-P DOF of the first processor
	// are followed by X-Y-Z-P DOF of the second one, and so on.
	//======================================================================

	FDSTAG      *fs;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscScalar dr;
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, rescal;
	PetscScalar eta, rho, Kb, IKdt, diag, pt, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar ***vKb,  ***vrho,  ***veta, ***vetaxy, ***vetaxz, ***vetayz;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// get density gradient stabilization parameters
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration
	rescal = md->rescal; // stencil rescaling flag

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// clear matrix coefficients
	ierr = MatZeroEntries(A); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp);   CHKERRQ(ierr);

	// access parameter vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		Kb   = vKb [k][j][i];
		rho  = vrho[k][j][i];
		eta  = veta[k][j][i];
		IKdt = 1.0/Kb/dt;

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

		// set pressure two-point constraints
		SET_PRES_TPC(bcp, i-1, j,   k,   i, 0,   cf[0])
		SET_PRES_TPC(bcp, i+1, j,   k,   i, mcx, cf[1])
		SET_PRES_TPC(bcp, i,   j-1, k,   j, 0,   cf[2])
		SET_PRES_TPC(bcp, i,   j+1, k,   j, mcy, cf[3])
		SET_PRES_TPC(bcp, i,   j,   k-1, k, 0,   cf[4])
		SET_PRES_TPC(bcp, i,   j,   k+1, k, mcz, cf[5])

		// compute local matrix
		getStiffMat(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

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
		ierr = MatSetValues(A, 7, idx, 7, idx, v,  ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxy[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k][j-1][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvy[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvy[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxz[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetayz[k][j][i];

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvy[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvy[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j-1][i];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	// assemble velocity-pressure matrix, remove constrained rows
	ierr = MatAIJAssemble(A, md->numSPC, md->SPCListMat, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatComputeDiag(MatData *md, PetscScalar pgamma, Mat D)
{
	//======================================================================
	// get preconditioning matrix diagonal
	//======================================================================

	// WARNING! avoid using matrix! use vectors instead (matrix-free)

	FDSTAG      *fs;
	PetscInt    idx[7];
	PetscScalar v[16];
	PetscScalar dr;
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, rescal;
	PetscScalar eta, rho, Kb, IKdt, diag, pt, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar ***vKb,  ***vrho,  ***veta, ***vetaxy, ***vetaxz, ***vetayz;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// get density gradient stabilization parameters
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration
	rescal = md->rescal; // stencil rescaling flag

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// clear matrix coefficients
	ierr = MatZeroEntries(D); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp);   CHKERRQ(ierr);

	// access parameter vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		Kb   = vKb [k][j][i];
		rho  = vrho[k][j][i];
		eta  = veta[k][j][i];
		IKdt = 1.0/Kb/dt;

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

		// set pressure two-point constraints
		SET_PRES_TPC(bcp, i-1, j,   k,   i, 0,   cf[0])
		SET_PRES_TPC(bcp, i+1, j,   k,   i, mcx, cf[1])
		SET_PRES_TPC(bcp, i,   j-1, k,   j, 0,   cf[2])
		SET_PRES_TPC(bcp, i,   j+1, k,   j, mcy, cf[3])
		SET_PRES_TPC(bcp, i,   j,   k-1, k, 0,   cf[4])
		SET_PRES_TPC(bcp, i,   j,   k+1, k, mcz, cf[5])

		// compute local matrix
		getStiffMatDiag(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

		// compute density gradient stabilization terms
		addDensGradStabilDiag(fssa, v, rho, dt, grav, fdx, fdy, fdz, bdx, bdy, bdz);

		// get global indices of the points:
		// vx_(i), vx_(i+1), vy_(j), vy_(j+1), vz_(k), vz_(k+1), p
		idx[0] = (PetscInt) ivx[k][j][i];
		idx[1] = (PetscInt) ivx[k][j][i+1];
		idx[2] = (PetscInt) ivy[k][j][i];
		idx[3] = (PetscInt) ivy[k][j+1][i];
		idx[4] = (PetscInt) ivz[k][j][i];
		idx[5] = (PetscInt) ivz[k+1][j][i];
		idx[6] = (PetscInt) ip[k][j][i];

		// update diagonal
		ierr = MatSetValue(D, idx[0], idx[0], v[0], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[1], idx[1], v[1], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[2], idx[2], v[2], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[3], idx[3], v[3], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[4], idx[4], v[4], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[5], idx[5], v[5], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[6], idx[6], v[6], ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxy[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k][j-1][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvy[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvy[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		ierr = MatSetValue(D, idx[0], idx[0], v[0],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[1], idx[1], v[5],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[2], idx[2], v[10], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[3], idx[3], v[15], ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxz[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		ierr = MatSetValue(D, idx[0], idx[0], v[0],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[1], idx[1], v[5],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[2], idx[2], v[10], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[3], idx[3], v[15], ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetayz[k][j][i];

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvy[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvy[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j-1][i];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		ierr = MatSetValue(D, idx[0], idx[0], v[0],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[1], idx[1], v[5],  ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[2], idx[2], v[10], ADD_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(D, idx[3], idx[3], v[15], ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	// assemble diagonal preconditioner matrix, remove constrained rows
	ierr = MatAIJAssemble(D, md->numSPC, md->SPCListMat, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//.........................   MONOLITHIC MATRIX   ...........................
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoCreate(
		PMatMono    *P,
		MatData     *md,
		PetscScalar  pgamma)
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
	PetscInt    ln, start, ind, nd, no, *d_nnz, *o_nnz;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// store evaluation context
	P->md     = md;
	P->pgamma = pgamma;

	// access context
	fs  = md->fs;
	dof = &fs->dof;

	// get number of local rows & global index of the first row
	ln    = dof->ln;
	start = dof->st;

	// allocate nonzero counter arrays
	ierr = makeIntArray(&d_nnz, NULL, ln); CHKERRQ(ierr);
	ierr = makeIntArray(&o_nnz, NULL, ln); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// clear iterator
	iter = 0;

	//---------
	// X-points
	//---------
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// create matrices & vectors
	ierr = MatAIJCreate(ln, ln, 0, d_nnz, 0, o_nnz, &P->A); CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(ln, start, &P->M);              CHKERRQ(ierr);

	// create work vector
	ierr = VecCreateMPI(PETSC_COMM_WORLD, dof->ln, PETSC_DETERMINE, &P->w); CHKERRQ(ierr);
	ierr = VecSetFromOptions(P->w);

	// clear work arrays
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	// attach near null space
	ierr = MatAIJSetNullSpace(P->A, md); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoDestroy(PMatMono *P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatDestroy (&P->A); CHKERRQ(ierr);
	ierr = MatDestroy (&P->M); CHKERRQ(ierr);
	ierr = VecDestroy (&P->w); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoAssemble(PMatMono *P)
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

	MatData     *md;
	FDSTAG      *fs;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscScalar dr;
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, rescal;
	PetscScalar eta, rho, Kb, IKdt, diag, pgamma, pt, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar ***vKb,  ***vrho,  ***veta, ***vetaxy, ***vetaxz, ***vetayz;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	md = P->md;
	fs = md->fs;

	// get density gradient stabilization parameters
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration
	rescal = md->rescal; // stencil rescaling flag

	// get penalty parameter
	pgamma = P->pgamma;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// clear matrix coefficients
	ierr = MatZeroEntries(P->A); CHKERRQ(ierr);
	ierr = MatZeroEntries(P->M); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp);   CHKERRQ(ierr);

	// access parameter vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		Kb   = vKb [k][j][i];
		rho  = vrho[k][j][i];
		eta  = veta[k][j][i];
		IKdt = 1.0/Kb/dt;

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

		// set pressure two-point constraints
		SET_PRES_TPC(bcp, i-1, j,   k,   i, 0,   cf[0])
		SET_PRES_TPC(bcp, i+1, j,   k,   i, mcx, cf[1])
		SET_PRES_TPC(bcp, i,   j-1, k,   j, 0,   cf[2])
		SET_PRES_TPC(bcp, i,   j+1, k,   j, mcy, cf[3])
		SET_PRES_TPC(bcp, i,   j,   k-1, k, 0,   cf[4])
		SET_PRES_TPC(bcp, i,   j,   k+1, k, mcz, cf[5])

		// compute local matrix
		getStiffMat(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

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
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxy[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k][j-1][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvy[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvy[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[1], cf[0], dr);

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
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxz[k][j][i];

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvx[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvx[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j][i-1];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dx, fdx, bdx, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetayz[k][j][i];

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// get boundary constraints
		pdofidx[0] = 1;   cf[0] = bcvy[k-1][j][i];
		pdofidx[1] = 0;   cf[1] = bcvy[k][j][i];
		pdofidx[2] = 3;   cf[2] = bcvz[k][j-1][i];
		pdofidx[3] = 2;   cf[3] = bcvz[k][j][i];

		// stencil rescaling
		RESCALE_STENCIL(rescal, dy, fdy, bdy, cf[3], cf[2], dr);
		RESCALE_STENCIL(rescal, dz, fdz, bdz, cf[1], cf[0], dr);

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

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// add to global matrix
		ierr = MatSetValues(P->A, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	// assemble velocity-pressure matrix, remove constrained rows
	ierr = MatAIJAssemble(P->A, md->numSPC, md->SPCListMat, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->M, md->numSPC, md->SPCListMat, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode PMatMonoPicard(Mat J, Vec x, Vec r)
{
	//=======================================================================
	// matrix-vector product with Picard matrix
	//
	// r = J*x = A*x - M*x
	//=======================================================================

	PMatMono *P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
//...........................   BLOCK MATRIX   ..............................
//---------------------------------------------------------------------------
PetscErrorCode PMatBlockCreate(
		PMatBlock   *P,
		MatData     *md,
		PetscScalar  pgamma,
		PetscInt     buildwBFBT,
		PetscInt     buildBvv)
{
	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscInt    lnp, lnv, startv, startp, nd, no, ind;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    *Avv_d_nnz, *Avv_o_nnz;
	PetscInt    *Avp_d_nnz, *Avp_o_nnz;
	PetscInt    *Apv_d_nnz, *Apv_o_nnz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// store evaluation context
	P->md     = md;
	P->pgamma = pgamma;

	// access context
	fs  = md->fs;
	dof = &fs->dof;

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
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// clear iterator (velocity matrices)
	iter = 0;

	//---------
	// X-points
	//---------
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);   CHKERRQ(ierr);

	// create matrices & vectors
	ierr = MatAIJCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &P->Avv);  CHKERRQ(ierr);
	ierr = MatAIJCreate(lnv, lnp, 0, Avp_d_nnz, 0, Avp_o_nnz, &P->Avp);  CHKERRQ(ierr);
	ierr = MatAIJCreate(lnp, lnv, 0, Apv_d_nnz, 0, Apv_o_nnz, &P->Apv);  CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(lnp, startp, &P->App);                       CHKERRQ(ierr);
	ierr = MatAIJCreateDiag(lnp, startp, &P->iS);                        CHKERRQ(ierr);

	if(buildBvv)
	{
		ierr = MatAIJCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &P->Bvv); CHKERRQ(ierr);
	}

	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &P->xv); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnp, PETSC_DETERMINE, &P->xp); CHKERRQ(ierr);
	ierr = VecSetFromOptions(P->xp);                                     CHKERRQ(ierr);
	ierr = VecSetFromOptions(P->xv);                                     CHKERRQ(ierr);
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

	// attach near null space
	ierr = MatAIJSetNullSpace(P->Avv, md); CHKERRQ(ierr);

	// create BFBT preconditioner
	if(buildwBFBT)
	{
		// allocate space
		ierr = PetscMalloc(sizeof(wBFBTData), (void**)&P->wbfbt); CHKERRQ(ierr);

		// clear object
		ierr = PetscMemzero(P->wbfbt, sizeof(wBFBTData)); CHKERRQ(ierr);

		// setup data
		ierr = wBFBTCreate(P->wbfbt, md); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatBlockDestroy(PMatBlock *P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatDestroy(&P->Avv); CHKERRQ(ierr);
	ierr = MatDestroy(&P->Avp); CHKERRQ(ierr);
	ierr = MatDestroy(&P->Apv); CHKERRQ(ierr);
	ierr = MatDestroy(&P->App); CHKERRQ(ierr);
	ierr = VecDestroy(&P->rv);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->rp);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->xv);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->xp);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv);  CHKERRQ(ierr);
	ierr = VecDestroy(&P->wp);  CHKERRQ(ierr);
	ierr = MatDestroy(&P->iS);  CHKERRQ(ierr);

	if(P->Bvv)
	{
		ierr = MatDestroy(&P->Bvv); CHKERRQ(ierr);
	}
	if(P->wbfbt)
	{
		ierr = wBFBTDestroy(P->wbfbt); CHKERRQ(ierr);
		ierr = PetscFree   (P->wbfbt); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode PMatBlockAssemble(PMatBlock *P)
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
	//
	// Bvv matrix will contain clean velocity block if necessary
	//======================================================================

	MatData     *md;
	FDSTAG      *fs;
	PetscInt    idx[7];
	PetscScalar v[49], vc[49], a[36], ac[36], d[6], g[6];
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, rho, Kb, IKdt, diag, pgamma, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar ***vKb,  ***vrho,  ***veta, ***vetaxy, ***vetaxz, ***vetayz;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// accees context
	md = P->md;
	fs = md->fs;

	// get density gradient stabilization parameters
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration

	// get penalty parameter
	pgamma = P->pgamma;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// clear matrix coefficients
	ierr = MatZeroEntries(P->Avv);  CHKERRQ(ierr);
	ierr = MatZeroEntries(P->Avp);  CHKERRQ(ierr);
	ierr = MatZeroEntries(P->Apv);  CHKERRQ(ierr);

	if(P->Bvv)
	{
		ierr = MatZeroEntries(P->Bvv); CHKERRQ(ierr);
	}

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip);  CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp);  CHKERRQ(ierr);

	// access parameter vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get density, shear & inverse bulk viscosities
		Kb   = vKb [k][j][i];
		rho  = vrho[k][j][i];
		eta  = veta[k][j][i];
		IKdt = 1.0/Kb/dt;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get pressure diagonal element (with penalty)
		diag = -IKdt - 1.0/(pgamma*eta);

		// set pressure two-point constraints
		SET_PRES_TPC(bcp, i-1, j,   k,   i, 0,   cf[0])
		SET_PRES_TPC(bcp, i+1, j,   k,   i, mcx, cf[1])
		SET_PRES_TPC(bcp, i,   j-1, k,   j, 0,   cf[2])
		SET_PRES_TPC(bcp, i,   j+1, k,   j, mcy, cf[3])
		SET_PRES_TPC(bcp, i,   j,   k-1, k, 0,   cf[4])
		SET_PRES_TPC(bcp, i,   j,   k+1, k, mcz, cf[5])

		// compute local matrix
		getStiffMat(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

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

		if(P->Bvv)
		{
			// copy clean stiffness matrix
			ierr = PetscMemcpy(vc, v, sizeof(v)); CHKERRQ(ierr);

			// constrain local matrix
			constrLocalMat(7, pdofidx, cf, vc);

			// extract operators
			getSubMat(vc, ac, d, g);

			// update global matrix
			ierr = MatSetValues(P->Bvv, 6, idx, 6, idx, ac, ADD_VALUES); CHKERRQ(ierr);
		}

		// compute velocity Schur complement
		if(P->pgamma != 1.0) getVelSchur(v, d, g);

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
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxy[k][j][i];

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

		if(P->Bvv)
		{
			ierr = MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
		}
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetaxz[k][j][i];

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

		if(P->Bvv)
		{
			ierr = MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
		}
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get effective viscosity
		eta = vetayz[k][j][i];

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

		if(P->Bvv)
		{
			ierr = MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES); CHKERRQ(ierr);
		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp);  CHKERRQ(ierr);

	// assemble velocity-pressure matrix blocks, remove constrained rows
	ierr = MatAIJAssemble(P->Avv, md->vNumSPC, md->vSPCListMat, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->Avp, md->vNumSPC, md->vSPCListMat, 0.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->Apv, md->pNumSPC, md->pSPCListMat, 0.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->App, md->pNumSPC, md->pSPCListMat, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->iS,  md->pNumSPC, md->pSPCListMat, 1.0); CHKERRQ(ierr);

	if(P->Bvv)
	{
		ierr = MatAIJAssemble(P->Bvv, md->vNumSPC, md->vSPCListMat, 1.0); CHKERRQ(ierr);
	}

	// assemble BFBT preconditioner
	if(P->wbfbt)
	{
		ierr = wBFBTAssemble(P->wbfbt); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatBlockPicard(Mat J, Vec x, Vec r)
{
	//=======================================================================
	// matrix-vector product with Picard matrix
	//
	// rv = Bvv*xv + Avp*xp
	// rp = Apv*xv + App*xp
	//=======================================================================
	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(J, (void**)&P); CHKERRQ(ierr);

	// extract solution blocks
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = MatMult(P->Apv,  P->xv, P->rp); CHKERRQ(ierr); // rp = Apv*xv

	ierr = MatMult(P->App,  P->xp, P->wp); CHKERRQ(ierr); // wp = App*xp

	ierr = VecAXPY(P->rp,   1.0,   P->wp); CHKERRQ(ierr); // rp = rp + wp

	ierr = MatMult(P->Avp,  P->xp, P->rv); CHKERRQ(ierr); // rv = Avp*xp

	if(P->Bvv)
	{
		ierr = MatMult(P->Bvv,  P->xv, P->wv); CHKERRQ(ierr); // wv = Bvv*xv
	}
	else
	{
		ierr = MatMult(P->Avv,  P->xv, P->wv); CHKERRQ(ierr); // wv = Avv*xv
	}

	ierr = VecAXPY(P->rv,   1.0,   P->wv); CHKERRQ(ierr); // rv = rv + wv

	// compose coupled residual
	ierr = VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
