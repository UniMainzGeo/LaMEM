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

//---------------------------------------------------------------------------
PetscErrorCode PMatCreate(MatData *md, Mat *A, PetscInt set_null_space)
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

	
	PetscFunctionBeginUser;

	// access context
	fs  = md->fs;
	dof = &fs->dof;

	// get number of local rows & global index of the first row
	ln    = dof->ln;
	start = dof->st;

	// allocate nonzero counter arrays
	PetscCall(makeIntArray(&d_nnz, NULL, ln));
	PetscCall(makeIntArray(&o_nnz, NULL, ln));

	// access index vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip));

	// clear iterator
	iter = 0;

	//---------
	// X-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip));

	// create matrix
	PetscCall(MatAIJCreate(ln, ln, 0, d_nnz, 0, o_nnz, A));

	// free counter arrays
	PetscCall(PetscFree(d_nnz));
	PetscCall(PetscFree(o_nnz));

	// attach near null space
	if(set_null_space)
	{
		PetscCall(MatAIJSetNullSpace((*A), md));
	}

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
	PetscCall(MatZeroEntries(A));

	// access index vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip));

	// access boundary constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp));

	// access parameter vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz));

	//---------------
	// central points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(A, 7, idx, 7, idx, v,  ADD_VALUES));
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES));
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES));
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(A, 4, idx, 4, idx, v, ADD_VALUES));
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip));

	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp));

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &vKb));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &vrho));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &veta));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &vetaxy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &vetaxz));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &vetayz));

	// assemble velocity-pressure matrix, remove constrained rows
	PetscCall(MatAIJAssemble(A, md->numSPC, md->SPCListMat, 1.0));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatDiagComp(MatData *md, PetscScalar pgamma, Mat M)
{
	// assemble diagonal compensation matrix for Picard matrix-vector product

	FDSTAG      *fs;
	PetscInt    idx;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, pt;
	PetscScalar ***ip;
	PetscScalar ***veta;

	
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// clear matrix coefficients
	PetscCall(MatZeroEntries(M));

	// access pressure index and shear viscosity vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,  &ip));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->eta, &veta));

	//---------------
	// central points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get shear viscosity
		eta  = veta[k][j][i];

		// compute penalty term
		pt = -1.0/(pgamma*eta);

		// get global indices of pressure point
		idx = (PetscInt) ip[k][j][i];

		// update penalty compensation matrices
		PetscCall(MatSetValue(M, idx, idx, pt, INSERT_VALUES));
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,  &ip));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->eta, &veta));

	// assemble diagonal compensation matrix, remove constrained rows
	PetscCall(MatAIJAssemble(M, md->numSPC, md->SPCListMat, 0.0));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//.........................   MONOLITHIC MATRIX   ...........................
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoCreate(
		PMatMono    *P,
		MatData     *md,
		PetscScalar  pgamma,
		PetscInt     set_null_space)
{
	DOFIndex *dof;

	
	PetscFunctionBeginUser;

	// store evaluation context
	P->md     = md;
	P->pgamma = pgamma;

	// access context
	dof = &md->fs->dof;

	// create matrix
	PetscCall(PMatCreate(md, &P->A, set_null_space));

	// create diagonal compensation matrix
	PetscCall(MatAIJCreateDiag(dof->ln, dof->st, &P->M));

	// create work vector
	PetscCall(VecCreateMPI(PETSC_COMM_WORLD, dof->ln, PETSC_DETERMINE, &P->w));
	PetscCall(VecSetFromOptions(P->w));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoDestroy(PMatMono *P)
{
	
	PetscFunctionBeginUser;

	PetscCall(MatDestroy(&P->A));
	PetscCall(MatDestroy(&P->M));
	PetscCall(VecDestroy(&P->w));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatMonoAssemble(PMatMono *P)
{
	
	PetscFunctionBeginUser;

	PetscCall(PMatAssemble(P->md, P->pgamma, P->A));
	PetscCall(PMatDiagComp(P->md, P->pgamma, P->M));

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

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(J, (void**)&P));

	// compute action of preconditioner matrix
	PetscCall(MatMult(P->A, x, r));

	// compute compensation
	PetscCall(MatMult(P->M, x, P->w));

	// update result
	PetscCall(VecAXPY(r, -1.0, P->w));

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
		PetscInt     buildBvv,
		PetscInt     set_null_space)
{
	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscInt    lnp, lnv, startv, startp, nd, no, ind;
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    *Avv_d_nnz, *Avv_o_nnz;
	PetscInt    *Avp_d_nnz, *Avp_o_nnz;
	PetscInt    *Apv_d_nnz, *Apv_o_nnz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	
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
	PetscCall(makeIntArray(&Avv_d_nnz, NULL, lnv));
	PetscCall(makeIntArray(&Avv_o_nnz, NULL, lnv));

	PetscCall(makeIntArray(&Avp_d_nnz, NULL, lnv));
	PetscCall(makeIntArray(&Avp_o_nnz, NULL, lnv));

	PetscCall(makeIntArray(&Apv_d_nnz, NULL, lnp));
	PetscCall(makeIntArray(&Apv_o_nnz, NULL, lnp));

	// access index vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip));

	// clear iterator (velocity matrices)
	iter = 0;

	//---------
	// X-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip));

	// create matrices & vectors
	PetscCall(MatAIJCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &P->Avv));
	PetscCall(MatAIJCreate(lnv, lnp, 0, Avp_d_nnz, 0, Avp_o_nnz, &P->Avp));
	PetscCall(MatAIJCreate(lnp, lnv, 0, Apv_d_nnz, 0, Apv_o_nnz, &P->Apv));
	PetscCall(MatAIJCreateDiag(lnp, startp, &P->App));

	if(buildBvv)
	{
		PetscCall(MatAIJCreate(lnv, lnv, 0, Avv_d_nnz, 0, Avv_o_nnz, &P->Bvv));
	}

	PetscCall(VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &P->xv));
	PetscCall(VecCreateMPI(PETSC_COMM_WORLD, lnp, PETSC_DETERMINE, &P->xp));
	PetscCall(VecSetFromOptions(P->xp));
	PetscCall(VecSetFromOptions(P->xv));
	PetscCall(VecDuplicate(P->xv, &P->rv));
	PetscCall(VecDuplicate(P->xv, &P->wv));
	PetscCall(VecDuplicate(P->xp, &P->rp));
	PetscCall(VecDuplicate(P->xp, &P->wp));

	// free counter arrays
	PetscCall(PetscFree(Avv_d_nnz));
	PetscCall(PetscFree(Avv_o_nnz));
	PetscCall(PetscFree(Avp_d_nnz));
	PetscCall(PetscFree(Avp_o_nnz));
	PetscCall(PetscFree(Apv_d_nnz));
	PetscCall(PetscFree(Apv_o_nnz));

	// attach near null space
	if(set_null_space)
	{
		PetscCall(MatAIJSetNullSpace(P->Avv, md));
	}

	// create BFBT preconditioner
	if(buildwBFBT)
	{
		// allocate space
		PetscCall(PetscMalloc(sizeof(wBFBTData), (void**)&P->wbfbt));

		// clear object
		PetscCall(PetscMemzero(P->wbfbt, sizeof(wBFBTData)));

		// setup data
		PetscCall(wBFBTCreate(P->wbfbt, md));
	}
	else
	{
		PetscCall(MatAIJCreateDiag(lnp, startp, &P->iS));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PMatBlockDestroy(PMatBlock *P)
{
	
	PetscFunctionBeginUser;

	PetscCall(MatDestroy(&P->Avv));
	PetscCall(MatDestroy(&P->Avp));
	PetscCall(MatDestroy(&P->Apv));
	PetscCall(MatDestroy(&P->App));
	PetscCall(VecDestroy(&P->rv));
	PetscCall(VecDestroy(&P->rp));
	PetscCall(VecDestroy(&P->xv));
	PetscCall(VecDestroy(&P->xp));
	PetscCall(VecDestroy(&P->wv));
	PetscCall(VecDestroy(&P->wp));

	if(P->Bvv)
	{
		PetscCall(MatDestroy(&P->Bvv));
	}

	if(P->wbfbt)
	{
		PetscCall(wBFBTDestroy(P->wbfbt));
		PetscCall(PetscFree   (P->wbfbt));
	}

	if(P->iS)
	{
		PetscCall(MatDestroy(&P->iS));
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
	PetscCall(MatZeroEntries(P->Avv));
	PetscCall(MatZeroEntries(P->Avp));
	PetscCall(MatZeroEntries(P->Apv));

	if(P->Bvv)
	{
		PetscCall(MatZeroEntries(P->Bvv));
	}

	// access index vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,   &ip));

	// access boundary constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp));

	// access parameter vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz));

	//---------------
	// central points
	//---------------

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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
			PetscCall(PetscMemcpy(vc, v, sizeof(v)));

			// constrain local matrix
			constrLocalMat(7, pdofidx, cf, vc);

			// extract operators
			getSubMat(vc, ac, d, g);

			// update global matrix
			PetscCall(MatSetValues(P->Bvv, 6, idx, 6, idx, ac, ADD_VALUES));
		}

		// compute velocity Schur complement
		if(P->pgamma != 1.0) getVelSchur(v, d, g);

		// constrain local matrix
		constrLocalMat(7, pdofidx, cf, v);

		// extract operators
		getSubMat(v, a, d, g);

		// update global matrices
		PetscCall(MatSetValues(P->Avv, 6, idx,   6, idx,   a,    ADD_VALUES));
		PetscCall(MatSetValues(P->Avp, 6, idx,   1, idx+6, g,    ADD_VALUES));
		PetscCall(MatSetValues(P->Apv, 1, idx+6, 6, idx,   d,    ADD_VALUES));
		PetscCall(MatSetValue (P->App, idx[6], idx[6], -IKdt,    INSERT_VALUES));

		if(P->iS)
		{
			PetscCall(MatSetValue (P->iS,  idx[6], idx[6], 1.0/diag, INSERT_VALUES));
		}
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES));

		if(P->Bvv)
		{
			PetscCall(MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES));
		}
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES));

		if(P->Bvv)
		{
			PetscCall(MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES));
		}
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
		PetscCall(MatSetValues(P->Avv, 4, idx, 4, idx, v, ADD_VALUES));

		if(P->Bvv)
		{
			PetscCall(MatSetValues(P->Bvv, 4, idx, 4, idx, v, ADD_VALUES));
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->ivx,  &ivx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->ivy,  &ivy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->ivz,  &ivz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,   &ip));

	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp));

	// assemble velocity-pressure matrix blocks, remove constrained rows
	PetscCall(MatAIJAssemble(P->Avv, md->vNumSPC, md->vSPCListMat, 1.0));
	PetscCall(MatAIJAssemble(P->Avp, md->vNumSPC, md->vSPCListMat, 0.0));
	PetscCall(MatAIJAssemble(P->Apv, md->pNumSPC, md->pSPCListMat, 0.0));
	PetscCall(MatAIJAssemble(P->App, md->pNumSPC, md->pSPCListMat, 1.0));

	if(P->Bvv)
	{
		PetscCall(MatAIJAssemble(P->Bvv, md->vNumSPC, md->vSPCListMat, 1.0));
	}

	// assemble BFBT preconditioner
	if(P->wbfbt)
	{
		PetscCall(wBFBTAssemble(P->wbfbt));
	}

	// assemble inverse viscosity preconditioner
	if(P->iS)
	{
		PetscCall(MatAIJAssemble(P->iS,  md->pNumSPC, md->pSPCListMat, 1.0));
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

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(J, (void**)&P));

	// extract solution blocks
	PetscCall(VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_REVERSE));

	PetscCall(MatMult(P->Apv,  P->xv, P->rp)); // rp = Apv*xv

	PetscCall(MatMult(P->App,  P->xp, P->wp)); // wp = App*xp

	PetscCall(VecAXPY(P->rp,   1.0,   P->wp)); // rp = rp + wp

	PetscCall(MatMult(P->Avp,  P->xp, P->rv)); // rv = Avp*xp

	if(P->Bvv)
	{
		PetscCall(MatMult(P->Bvv,  P->xv, P->wv)); // wv = Bvv*xv
	}
	else
	{
		PetscCall(MatMult(P->Avv,  P->xv, P->wv)); // wv = Avv*xv
	}

	PetscCall(VecAXPY(P->rv,   1.0,   P->wv)); // rv = rv + wv

	// compose coupled residual
	PetscCall(VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_FORWARD));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
