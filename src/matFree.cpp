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
//...............   MATRIX-FREE JACOBIAN AND PRECONDITIONER  ................
//---------------------------------------------------------------------------

#include "LaMEM.h"
#include "matFree.h"
#include "fdstag.h"
#include "JacRes.h"
#include "matrix.h"
#include "bc.h"
#include "tssolve.h"

//---------------------------------------------------------------------------
PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec y)
{
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(A, (void**)&jr); CHKERRQ(ierr);

	// copy solution from global to local vectors, do not enforce boundary constraints
	ierr = JacResCopyVelNoBC(jr, x); CHKERRQ(ierr);

	// compute matrix-free matrix-vector product
	ierr = JacResPicardMatFree(jr); CHKERRQ(ierr);

	// assemble residuals to global vector
	ierr = JacResCopyRes(jr, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode JacResPicardMatFree(JacRes *jr)
{
	FDSTAG     *fs;
	BCCtx      *bc;
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar dx, dy, dz, tx, ty, tz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar sxx, syy, szz, sxy, sxz, syz;
	PetscScalar dxx, dyy, dzz, dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar eta, theta, tr, rho, IKdt, pc, dt, fssa, *grav;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc, ***p;

	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar cf[6];
//	PetscInt    rescal;
//	PetscScalar dr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs     = jr->fs;
	bc     = jr->bc;
	dt     = jr->ts->dt;      // time step
	fssa   = jr->ctrl.FSSA;   // density gradient penalty parameter
	grav   = jr->ctrl.grav;   // gravity acceleration
//    rescal = jr->ctrl.rescal; // stencil rescaling flag

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

    // clear local residual vectors
	ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->gc);  CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,   &bcp);   CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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

		// compute velocity gradients
		dxx = (vx[k][j][i+1] - vx[k][j][i])/dx;
		dyy = (vy[k][j+1][i] - vy[k][j][i])/dy;
		dzz = (vz[k+1][j][i] - vz[k][j][i])/dz;

		// compute & store volumetric strain rate
		theta = dxx + dyy + dzz;

		// compute & store total deviatoric strain rates
		tr   = theta/3.0;
		dxx -= tr;
		dyy -= tr;
		dzz -= tr;

		// access current pressure
		pc = p[k][j][i];

		// set pressure two-point constraints
		SET_PRES_TPC(bcp, i-1, j,   k,   i, 0,   cf[0])
		SET_PRES_TPC(bcp, i+1, j,   k,   i, mcx, cf[1])
		SET_PRES_TPC(bcp, i,   j-1, k,   j, 0,   cf[2])
		SET_PRES_TPC(bcp, i,   j+1, k,   j, mcy, cf[3])
		SET_PRES_TPC(bcp, i,   j,   k-1, k, 0,   cf[4])
		SET_PRES_TPC(bcp, i,   j,   k+1, k, mcz, cf[5])

		// compute deviatoric stresses
		sxx = 2.0*eta*dxx;
		syy = 2.0*eta*dyy;
		szz = 2.0*eta*dzz;

		// compute stabilization terms (lumped approximation)
		tx = -fssa*dt*(rho*grav[0]);
		ty = -fssa*dt*(rho*grav[1]);
		tz = -fssa*dt*(rho*grav[2]);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
		fx[k][j][i] -= ((sxx - cf[0]*pc) + vx[k][j][i]*tx)/bdx;   fx[k][j][i+1] += ((sxx - cf[1]*pc) + vx[k][j][i+1]*tx)/fdx;
		fy[k][j][i] -= ((syy - cf[2]*pc) + vy[k][j][i]*ty)/bdy;   fy[k][j+1][i] += ((syy - cf[3]*pc) + vy[k][j+1][i]*ty)/fdy;
		fz[k][j][i] -= ((szz - cf[4]*pc) + vz[k][j][i]*tz)/bdz;   fz[k+1][j][i] += ((szz - cf[5]*pc) + vz[k+1][j][i]*tz)/fdz;

		// mass
		gc[k][j][i] = -IKdt*pc - theta;
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// initialize scaling factors
		cf[0] = 1.0;
		cf[1] = 1.0;
		cf[2] = 1.0;
		cf[3] = 1.0;

		// set velocity two-point constraints
		SET_VEL_TPC(bcvx, i,   j-1, k, j, 0,   cf[1], cf[0])
		SET_VEL_TPC(bcvx, i,   j,   k, j, mny, cf[0], cf[1])
		SET_VEL_TPC(bcvy, i-1, j,   k, i, 0,   cf[3], cf[2])
		SET_VEL_TPC(bcvy, i,   j,   k, i, mnx, cf[2], cf[3])

		// get effective viscosity
		eta = jr->svXYEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// compute velocity gradients
		dvxdy = (cf[1]*vx[k][j][i] - cf[0]*vx[k][j-1][i])/dy;
		dvydx = (cf[3]*vy[k][j][i] - cf[2]*vy[k][j][i-1])/dx;

		// compute stress
		sxy = eta*(dvxdy + dvydx);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// momentum
		fx[k][j-1][i] -= sxy/bdy;   fx[k][j][i] += sxy/fdy;
		fy[k][j][i-1] -= sxy/bdx;   fy[k][j][i] += sxy/fdx;

	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// initialize scaling factors
		cf[0] = 1.0;
		cf[1] = 1.0;
		cf[2] = 1.0;
		cf[3] = 1.0;

		// set velocity two-point constraints
		SET_VEL_TPC(bcvx, i,   j,   k-1, k, 0,   cf[1], cf[0])
		SET_VEL_TPC(bcvx, i,   j,   k,   k, mnz, cf[0], cf[1])
		SET_VEL_TPC(bcvz, i-1, j,   k,   i, 0,   cf[3], cf[2])
		SET_VEL_TPC(bcvz, i,   j,   k,   i, mnx, cf[2], cf[3])

		// get effective viscosity
		eta = jr->svXZEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvxdz = (cf[1]*vx[k][j][i] - cf[0]*vx[k-1][j][i])/dz;
		dvzdx = (cf[3]*vz[k][j][i] - cf[2]*vz[k][j][i-1])/dx;

		// compute stress
		sxz = eta*(dvxdz + dvzdx);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// momentum
		fx[k-1][j][i] -= sxz/bdz;   fx[k][j][i] += sxz/fdz;
		fz[k][j][i-1] -= sxz/bdx;   fz[k][j][i] += sxz/fdx;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// initialize scaling factors
		cf[0] = 1.0;
		cf[1] = 1.0;
		cf[2] = 1.0;
		cf[3] = 1.0;

		// set velocity two-point constraints
		SET_VEL_TPC(bcvy, i,   j,   k-1, k, 0,   cf[1], cf[0])
		SET_VEL_TPC(bcvy, i,   j,   k,   k, mnz, cf[0], cf[1])
		SET_VEL_TPC(bcvz, i,   j-1, k,   j, 0,   cf[3], cf[2])
		SET_VEL_TPC(bcvz, i,   j,   k,   j, mny, cf[2], cf[3])

		// get effective viscosity
		eta = jr->svYZEdge[iter++].svDev.eta;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvydz = (cf[1]*vy[k][j][i] - cf[0]*vy[k-1][j][i])/dz;
		dvzdy = (cf[3]*vz[k][j][i] - cf[2]*vz[k][j-1][i])/dy;

		// compute stress
		syz = eta*(dvydz + dvzdy);

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// update momentum residuals
		fy[k-1][j][i] -= syz/bdz;   fy[k][j][i] += syz/fdz;
		fz[k][j-1][i] -= syz/bdy;   fz[k][j][i] += syz/fdy;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


