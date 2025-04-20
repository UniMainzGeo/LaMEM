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
#include "matData.h"
#include "matrix.h"

//---------------------------------------------------------------------------
PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec f)
{
	MatData *md;
	FDSTAG  *fs;
	Vec      lvx, lvy, lvz, gp;
	Vec      lfx, lfy, lfz, gc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(A, (void**)&md); CHKERRQ(ierr);

	fs = md->fs;

	// get temporary local vectors
	ierr = DMGetLocalVector(fs->DA_X, &lvx); CHKERRQ(ierr);
	ierr = DMGetLocalVector(fs->DA_Y, &lvy); CHKERRQ(ierr);
	ierr = DMGetLocalVector(fs->DA_Z, &lvz); CHKERRQ(ierr);
	ierr = DMGetLocalVector(fs->DA_X, &lfx); CHKERRQ(ierr);
	ierr = DMGetLocalVector(fs->DA_Y, &lfy); CHKERRQ(ierr);
	ierr = DMGetLocalVector(fs->DA_Z, &lfz); CHKERRQ(ierr);

	// get temporary global vectors
	ierr = DMGetGlobalVector(fs->DA_CEN, &gp);  CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_CEN, &gc);  CHKERRQ(ierr);

	// access solution vector
	ierr = MatFreeGetSol(md, x, lvx, lvy, lvz, gp); CHKERRQ(ierr);

	// compute matrix-vector product
	ierr = MatFreeGetPicard(md, lvx, lvy, lvz, gp, lfx, lfy, lfz, gc); CHKERRQ(ierr);

	// assemble residual
	ierr = MatFreeAssembleRes(md, f, lfx, lfy, lfz, gc); CHKERRQ(ierr);

	// restore temporary local vectors
	ierr = DMRestoreLocalVector(fs->DA_X, &lvx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_Y, &lvy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_Z, &lvz); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_X, &lfx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_Y, &lfy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_Z, &lfz); CHKERRQ(ierr);

	// restore temporary global vectors
	ierr = DMRestoreGlobalVector(fs->DA_CEN, &gp);  CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_CEN, &gc);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode MatFreeGetPicard(MatData *md,
		Vec lvx, Vec lvy, Vec lvz, Vec gp,
		Vec lfx, Vec lfy, Vec lfz, Vec gc)
{
	FDSTAG     *fs;
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar dx, dy, dz, tx, ty, tz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar sxx, syy, szz, sxy, sxz, syz;
	PetscScalar dxx, dyy, dzz, dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;

	PetscScalar eta, theta, tr, rho, Kb, IKdt, pc, dt, fssa, *grav;

	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***c, ***p;

	PetscScalar ***vKb,  ***vrho,  ***veta, ***vetaxy, ***vetaxz, ***vetayz;


	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscScalar cf[6];
//	PetscInt    rescal;
//	PetscScalar dr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs     = md->fs;   // grid context
	dt     = md->dt;   // time step
	fssa   = md->fssa; // density gradient penalty parameter
	grav   = md->grav; // gravity acceleration
//    rescal = jr->ctrl.rescal; // stencil rescaling flag

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

    // clear residual components
	ierr = VecZeroEntries(lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(lfz); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_X,   lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, gc,   &c);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, gp,   &p);   CHKERRQ(ierr);

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


	//-------------------------------
	// central points
	//-------------------------------
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
		c[k][j][i] = -IKdt*pc - theta;
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
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
		eta = vetaxy[k][j][i];

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
		eta = vetaxz[k][j][i];

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
		eta = vetayz[k][j][i];

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
	ierr = DMDAVecGetArray(fs->DA_X,   lvx,       &vx);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   lvy,       &vy);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   lvz,       &vz);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   lfx,       &fx);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   lfy,       &fy);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   lfz,       &fz);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, gc,        &c);      CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, gp,        &p);      CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   md->bcvx,  &bcvx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->bcvy,  &bcvy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->bcvz,  &bcvz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->bcp,   &bcp);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeGetSol(MatData *md, Vec x, Vec lvx, Vec lvy, Vec lvz, Vec gp)
{
	// split coupled solution vector into components, do not enforce boundary constraints

	FDSTAG            *fs;
	PetscScalar       *vx, *vy, *vz, *p;
	Vec                gvx, gvy, gvz;
	const PetscScalar *sol, *iter;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary vectors
	ierr = DMGetGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	// access temporary vectors
	ierr = VecGetArray    (gvx, &vx);  CHKERRQ(ierr);
	ierr = VecGetArray    (gvy, &vy);  CHKERRQ(ierr);
	ierr = VecGetArray    (gvz, &vz);  CHKERRQ(ierr);
	ierr = VecGetArray    (gp,  &p);   CHKERRQ(ierr);
	ierr = VecGetArrayRead(x,   &sol); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = sol;

	ierr  = PetscMemcpy(vx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(vy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(vz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr = PetscMemcpy(p,   iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray    (gvx, &vx);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gvy, &vy);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gvz, &vz);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gp,  &p);   CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(x,   &sol); CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_X,   gvx, lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y,   gvy, lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z,   gvz, lvz)

	// restore temporary global vectors
	ierr = DMRestoreGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeAssembleRes(MatData *md, Vec f, Vec lfx, Vec lfy, Vec lfz, Vec gc)
{
	// assemble residual components into coupled vector, enforce boundary constraints

	FDSTAG      *fs;
	PetscInt     i, num, *list;
	PetscScalar *fx, *fy, *fz, *c, *res, *iter;
	Vec          gfx, gfy, gfz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary global vectors
	ierr = DMGetGlobalVector(fs->DA_X, &gfx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Y, &gfy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Z, &gfz); CHKERRQ(ierr);

	// assemble residual components
	LOCAL_TO_GLOBAL(fs->DA_X, lfx, gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, lfy, gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, lfz, gfz)

	// access vectors
	ierr = VecGetArray(gfx, &fx);  CHKERRQ(ierr);
	ierr = VecGetArray(gfy, &fy);  CHKERRQ(ierr);
	ierr = VecGetArray(gfz, &fz);  CHKERRQ(ierr);
	ierr = VecGetArray(gc,  &c);   CHKERRQ(ierr);
	ierr = VecGetArray(f,   &res); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = res;

	ierr  = PetscMemcpy(iter, fx, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, fy, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, fz, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr  = PetscMemcpy(iter, c,  (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// zero out constrained residuals (velocity)
	num   = md->vNumSPC;
	list  = md->vSPCList;

	for(i = 0; i < num; i++) res[list[i]] = 0.0;

	// zero out constrained residuals (pressure)
	num   = md->pNumSPC;
	list  = md->pSPCList;

	for(i = 0; i < num; i++) res[list[i]] = 0.0;

	// restore access
	ierr = VecRestoreArray(gfx, &fx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(gfy, &fy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(gfz, &fz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(gc,  &c);   CHKERRQ(ierr);
	ierr = VecRestoreArray(f,   &res); CHKERRQ(ierr);

	// restore temporary global vectors
	ierr = DMRestoreGlobalVector(fs->DA_X, &gfx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Y, &gfy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Z, &gfz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
