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
#include "multigrid.h"

//---------------------------------------------------------------------------
// interface functions
//---------------------------------------------------------------------------

PetscErrorCode MatFreeApplyPicard(Mat A, Vec x, Vec f)
{
	// this function corresponds to MATOP_MULT operation (f = A*x)

	MatData     *md;
	PetscScalar  cfInvEta;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(A, (void**)&md); CHKERRQ(ierr);

	// do not add inverse viscosity term to pressure diagonal matrix (Picard operator)
	cfInvEta = 0.0;

	ierr = MatFreeComputeLinearOperator(md, x, f, cfInvEta); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyPreconditioner(Mat A, Vec x, Vec f)
{
	// this function corresponds to MATOP_MULT operation (f = A*x)

	MatData     *md;
	PetscScalar  cfInvEta;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(A, (void**)&md); CHKERRQ(ierr);

	// add inverse viscosity term to pressure diagonal matrix (preconditioner operator)
	cfInvEta = 1.0;

	ierr = MatFreeComputeLinearOperator(md, x, f, cfInvEta); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyRestrict(Mat R, Vec vf, Vec vc)
{
	// this function corresponds to MATOP_MULT operation (vc = R*vf)

	MGInterp *mgi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(R, (void**)&mgi); CHKERRQ(ierr);

	ierr = MatFreeComputeRestrict(mgi, vf, vc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeUpdateRestrict(Mat R, Vec vf, Vec vcb, Vec vc)
{
	// this function corresponds to MATOP_MULT_ADD operation (vc = vcb + R*vf)

	MGInterp *mgi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(R, (void**)&mgi); CHKERRQ(ierr);

	ierr = MatFreeComputeRestrict(mgi, vf, mgi->wc); CHKERRQ(ierr);

	// update coarse grid base vector
	if(vcb == vc)
	{
		ierr = VecAXPY(vc, 1.0, mgi->wc); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecWAXPY(vc, 1.0, vcb, mgi->wc); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyProlong(Mat P, Vec vc, Vec vf)
{
	// this function corresponds to MATOP_MULT operation (vf = P*vc)

	MGInterp *mgi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(P, (void**)&mgi); CHKERRQ(ierr);

	ierr = MatFreeComputeProlong(mgi, vc, vf);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeUpdateProlong(Mat P, Vec vc, Vec vfb, Vec vf)
{
	// this function corresponds to MATOP_MULT_ADD operation (vf = vfb + P*vc)

	MGInterp *mgi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(P, (void**)&mgi); CHKERRQ(ierr);

	ierr = MatFreeComputeProlong(mgi, vc, mgi->wf);

	// update fine grid base vector
	if(vfb == vf)
	{
		ierr = VecAXPY(vf, 1.0, mgi->wf); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecWAXPY(vf, 1.0, vfb, mgi->wf); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
// main computation functions
//---------------------------------------------------------------------------

PetscErrorCode MatFreeComputeLinearOperator(MatData *md, Vec x, Vec f, PetscScalar cfInvEta)
{
	FDSTAG  *fs;
	Vec      vx, vy, vz, p;
	Vec      fx, fy, fz, c;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// get temporary solution vectors
	ierr = DMGetLocalVector (fs->DA_X,   &vx); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fs->DA_Y,   &vy); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fs->DA_Z,   &vz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_CEN, &p);  CHKERRQ(ierr);

	// get temporary residual vectors
	ierr = DMGetLocalVector (fs->DA_X,   &fx); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fs->DA_Y,   &fy); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fs->DA_Z,   &fz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_CEN, &c);  CHKERRQ(ierr);

	// split solution vector
	ierr = MatFreeSplitVec(md, x, vx, vy, vz, p); CHKERRQ(ierr);

	// compute matrix-vector product
	ierr = MatFreeGetLinearOperator(md, vx, vy, vz, p, fx, fy, fz, c, cfInvEta); CHKERRQ(ierr);

	// assemble residual
	ierr = MatFreeAssembleVec(md, f, fx, fy, fz, c); CHKERRQ(ierr);

	// restore temporary vectors
	ierr = DMRestoreLocalVector (fs->DA_X,   &vx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fs->DA_Y,   &vy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fs->DA_Z,   &vz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_CEN, &p);  CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fs->DA_X,   &fx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fs->DA_Y,   &fy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fs->DA_Z,   &fz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_CEN, &c);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeComputeRestrict(MGInterp *mgi, Vec vf, Vec vc)
{
	MatData  *coarse, *fine;
	Vec       fx, fy, fz, fp;
	Vec       cx, cy, cz, cp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	coarse = mgi->coarse;
	fine   = mgi->fine;

	// get temporary fine grid vectors
	ierr = DMGetLocalVector (fine->fs->DA_X,     &fx); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fine->fs->DA_Y,     &fy); CHKERRQ(ierr);
	ierr = DMGetLocalVector (fine->fs->DA_Z,     &fz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fine->fs->DA_CEN,   &fp); CHKERRQ(ierr);

	// get temporary coarse grid vectors
	ierr = DMGetGlobalVector(coarse->fs->DA_X,   &cx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(coarse->fs->DA_Y,   &cy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(coarse->fs->DA_Z,   &cz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(coarse->fs->DA_CEN, &cp); CHKERRQ(ierr);

	// split fine grid vector
	ierr = MatFreeSplitVec(fine, vf, fx, fy, fz, fp); CHKERRQ(ierr);

	// compute restriction to coarse grid
	ierr =  MatFreeGetRestrict(coarse, fine, fx, fy, fz, fp, cx, cy, cz, cp); CHKERRQ(ierr);

	// combine coarse grid vector
	ierr = MatFreeCombineVec(coarse, vc, cx, cy, cz, cp); CHKERRQ(ierr);

	// restore temporary vectors
	ierr = DMRestoreLocalVector (fine->fs->DA_X,     &fx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fine->fs->DA_Y,     &fy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (fine->fs->DA_Z,     &fz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fine->fs->DA_CEN,   &fp); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(coarse->fs->DA_X,   &cx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(coarse->fs->DA_Y,   &cy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(coarse->fs->DA_Z,   &cz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(coarse->fs->DA_CEN, &cp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeComputeProlong(MGInterp *mgi, Vec vc, Vec vf)
{
	MatData  *coarse, *fine;
	Vec       fx, fy, fz, fp;
	Vec       cx, cy, cz, cp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	coarse = mgi->coarse;
	fine   = mgi->fine;

	// get temporary fine grid vectors
	ierr = DMGetGlobalVector(fine->fs->DA_X,     &fx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fine->fs->DA_Y,     &fy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fine->fs->DA_Z,     &fz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fine->fs->DA_CEN,   &fp); CHKERRQ(ierr);

	// get temporary coarse grid vectors
	ierr = DMGetLocalVector (coarse->fs->DA_X,   &cx); CHKERRQ(ierr);
	ierr = DMGetLocalVector (coarse->fs->DA_Y,   &cy); CHKERRQ(ierr);
	ierr = DMGetLocalVector (coarse->fs->DA_Z,   &cz); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(coarse->fs->DA_CEN, &cp); CHKERRQ(ierr);

	// split coarse grid vector
	ierr = MatFreeSplitVec(coarse, vc, cx, cy, cz, cp); CHKERRQ(ierr);

	// compute prolongation on fine grid
	ierr = MatFreeGetProlong(coarse, fine, fx, fy, fz, fp, cx, cy, cz, cp); CHKERRQ(ierr);

	// combine fine grid vector
	ierr = MatFreeCombineVec(fine, vf, fx, fy, fz, fp); CHKERRQ(ierr);

	// restore temporary vectors
	ierr = DMRestoreGlobalVector(fine->fs->DA_X,     &fx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fine->fs->DA_Y,     &fy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fine->fs->DA_Z,     &fz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fine->fs->DA_CEN,   &fp); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (coarse->fs->DA_X,   &cx); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (coarse->fs->DA_Y,   &cy); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector (coarse->fs->DA_Z,   &cz); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(coarse->fs->DA_CEN, &cp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
// helper functions
//---------------------------------------------------------------------------

PetscErrorCode MatFreeSplitVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp)
{
	// split vector into blocks, assign ghost points

	FDSTAG            *fs;
	PetscScalar       *vx, *vy, *vz, *vp;
	Vec                gvx, gvy, gvz;
	const PetscScalar *va, *iter;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary vectors
	ierr = DMGetGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	// access temporary vectors
	ierr = VecGetArray    (gvx, &vx); CHKERRQ(ierr);
	ierr = VecGetArray    (gvy, &vy); CHKERRQ(ierr);
	ierr = VecGetArray    (gvz, &vz); CHKERRQ(ierr);
	ierr = VecGetArray    (gvp, &vp); CHKERRQ(ierr);
	ierr = VecGetArrayRead(v,   &va); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = va;

	ierr  = PetscMemcpy(vx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(vy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(vz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr = PetscMemcpy(vp,  iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray    (gvx, &vx);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gvy, &vy);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gvz, &vz);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (gvp, &vp);  CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(v,   &va);  CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_X, gvx, lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y, gvy, lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z, gvz, lvz)

	// restore temporary global vectors
	ierr = DMRestoreGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeAssembleVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp)
{
	// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions

	FDSTAG *fs;
	Vec     gvx, gvy, gvz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary global vectors
	ierr = DMGetGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	// assemble ghost point contributions
	LOCAL_TO_GLOBAL(fs->DA_X, lvx, gvx)
	LOCAL_TO_GLOBAL(fs->DA_Y, lvy, gvy)
	LOCAL_TO_GLOBAL(fs->DA_Z, lvz, gvz)

	// combine vector blocks
	ierr = MatFreeCombineVec(md, v, gvx, gvy, gvz, gvp); CHKERRQ(ierr);

	// restore temporary global vectors
	ierr = DMRestoreGlobalVector(fs->DA_X, &gvx); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Y, &gvy); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(fs->DA_Z, &gvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode MatFreeCombineVec(MatData *md, Vec v, Vec gvx, Vec gvy, Vec gvz, Vec gvp)
{
	// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions

	FDSTAG      *fs;
	PetscInt     i, num, *list;
	PetscScalar *vx, *vy, *vz, *vp, *va, *iter;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = md->fs;

	// access vectors
	ierr = VecGetArray(gvx, &vx); CHKERRQ(ierr);
	ierr = VecGetArray(gvy, &vy); CHKERRQ(ierr);
	ierr = VecGetArray(gvz, &vz); CHKERRQ(ierr);
	ierr = VecGetArray(gvp, &vp); CHKERRQ(ierr);
	ierr = VecGetArray(v,   &va); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = va;

	ierr  = PetscMemcpy(iter, vx, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(iter, vy, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(iter, vz, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr  = PetscMemcpy(iter, vp,  (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// zero out constrained residuals (velocity)
	num   = md->vNumSPC;
	list  = md->vSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = 0.0;

	// zero out constrained residuals (pressure)
	num   = md->pNumSPC;
	list  = md->pSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = 0.0;

	// restore access
	ierr = VecRestoreArray(gvx, &vx); CHKERRQ(ierr);
	ierr = VecRestoreArray(gvy, &vy); CHKERRQ(ierr);
	ierr = VecRestoreArray(gvz, &vz); CHKERRQ(ierr);
	ierr = VecRestoreArray(gvp, &vp); CHKERRQ(ierr);
	ierr = VecRestoreArray(v,   &va); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
// low-level functions
//---------------------------------------------------------------------------

PetscErrorCode MatFreeGetLinearOperator(MatData *md,
		Vec lvx, Vec lvy, Vec lvz, Vec gp,
		Vec lfx, Vec lfy, Vec lfz, Vec gc,
		PetscScalar cfInvEta)
{
	// cfInvEta - inverse viscosity term prefactor
	// 0.0      - Picard operator
	// 1.0      - preconditioner operator

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs     = md->fs;     // grid context
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration
//	rescal = md->rescal; // stencil rescaling flag

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
		c[k][j][i] = -(IKdt + cfInvEta/eta)*pc - theta;
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
	ierr = DMDAVecRestoreArray(fs->DA_X,   lvx,       &vx);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   lvy,       &vy);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   lvz,       &vz);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   lfx,       &fx);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   lfy,       &fy);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   lfz,       &fz);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, gc,        &c);      CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, gp,        &p);      CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->bcvx,  &bcvx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->bcvy,  &bcvy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->bcvz,  &bcvz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->bcp,   &bcp);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &vKb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &vrho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &veta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &vetaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &vetaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &vetayz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeGetRestrict(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp)
{
	PetscScalar vs[12], sum;
	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fxa, ***fya, ***fza, ***fpa;
	PetscScalar ***cxa, ***cya, ***cza, ***cpa;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fx, &fxa);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fy, &fya);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fz, &fza);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fp, &fpa);   CHKERRQ(ierr);

	// access vectors in coarse grid
	ierr = DMDAVecGetArray(coarse->fs->DA_X,   cx, &cxa); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Y,   cy, &cya); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Z,   cz, &cza); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, cp, &cpa); CHKERRQ(ierr);

	// set velocity weights
	vs[0 ] = 1.0/16.0;
	vs[1 ] = 1.0/16.0;
	vs[2 ] = 1.0/16.0;
	vs[3 ] = 1.0/16.0;
	vs[4 ] = 1.0/8.0;
	vs[5 ] = 1.0/8.0;
	vs[6 ] = 1.0/8.0;
	vs[7 ] = 1.0/8.0;
	vs[8 ] = 1.0/16.0;
	vs[9 ] = 1.0/16.0;
	vs[10] = 1.0/16.0;
	vs[11] = 1.0/16.0;

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict from fine grid stencil
		sum = fxa[K  ][J  ][I-1]*vs[0 ]
		+     fxa[K  ][J+1][I-1]*vs[1 ]
		+     fxa[K+1][J  ][I-1]*vs[2 ]
		+     fxa[K+1][J+1][I-1]*vs[3 ]
		+     fxa[K  ][J  ][I  ]*vs[4 ]
		+     fxa[K  ][J+1][I  ]*vs[5 ]
		+     fxa[K+1][J  ][I  ]*vs[6 ]
		+     fxa[K+1][J+1][I  ]*vs[7 ]
		+     fxa[K  ][J  ][I+1]*vs[8 ]
		+     fxa[K  ][J+1][I+1]*vs[9 ]
		+     fxa[K+1][J  ][I+1]*vs[10]
		+     fxa[K+1][J+1][I+1]*vs[11];

		// store coarse grid value
		cxa[k][j][i] = sum;
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict from fine grid stencil
		sum = fya[K  ][J-1][I  ]*vs[0 ]
		+     fya[K  ][J-1][I+1]*vs[1 ]
		+     fya[K+1][J-1][I  ]*vs[2 ]
		+     fya[K+1][J-1][I+1]*vs[3 ]
		+     fya[K  ][J  ][I  ]*vs[4 ]
		+     fya[K  ][J  ][I+1]*vs[5 ]
		+     fya[K+1][J  ][I  ]*vs[6 ]
		+     fya[K+1][J  ][I+1]*vs[7 ]
		+     fya[K  ][J+1][I  ]*vs[8 ]
		+     fya[K  ][J+1][I+1]*vs[9 ]
		+     fya[K+1][J+1][I  ]*vs[10]
    	+     fya[K+1][J+1][I+1]*vs[11];

		// store coarse grid value
		cya[k][j][i] = sum;
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict from fine grid stencil
		sum = fza[K-1][J  ][I  ]*vs[0 ]
		+     fza[K-1][J  ][I+1]*vs[1 ]
		+     fza[K-1][J+1][I  ]*vs[2 ]
		+     fza[K-1][J+1][I+1]*vs[3 ]
		+     fza[K  ][J  ][I  ]*vs[4 ]
		+     fza[K  ][J  ][I+1]*vs[5 ]
		+     fza[K  ][J+1][I  ]*vs[6 ]
		+     fza[K  ][J+1][I+1]*vs[7 ]
		+     fza[K+1][J  ][I  ]*vs[8 ]
		+     fza[K+1][J  ][I+1]*vs[9 ]
		+     fza[K+1][J+1][I  ]*vs[10]
    	+     fza[K+1][J+1][I+1]*vs[11];

		// store coarse grid value
		cza[k][j][i] = sum;
	}
	END_STD_LOOP

	if(coarse->idxmod == _IDX_COUPLED_)
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
		ierr = DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get fine grid indices
			I = 2*i;
			J = 2*j;
			K = 2*k;

			// restrict from fine grid stencil
			sum = fpa[K  ][J  ][I  ]*vs[0]
			+     fpa[K  ][J  ][I+1]*vs[1]
			+     fpa[K  ][J+1][I  ]*vs[2]
			+     fpa[K  ][J+1][I+1]*vs[3]
			+     fpa[K+1][J  ][I  ]*vs[4]
			+     fpa[K+1][J  ][I+1]*vs[5]
			+     fpa[K+1][J+1][I  ]*vs[6]
			+     fpa[K+1][J+1][I+1]*vs[7];

			// store coarse grid value
			cpa[k][j][i] = sum;
		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fx, &fxa);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fy, &fya);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fz, &fza);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fp, &fpa);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(coarse->fs->DA_X,   cx, &cxa); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Y,   cy, &cya); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Z,   cz, &cza); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, cp, &cpa); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeGetProlong(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp)
{
	PetscScalar vsf[8], vsr[4], sum;
	PetscInt    I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fxa, ***fya, ***fza, ***fpa;
	PetscScalar ***cxa, ***cya, ***cza, ***cpa;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fx, &fxa);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fy, &fya);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fz, &fza);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fp, &fpa);   CHKERRQ(ierr);

	// access vectors in coarse grid
	ierr = DMDAVecGetArray(coarse->fs->DA_X,   cx, &cxa); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Y,   cy, &cya); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Z,   cz, &cza); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, cp, &cpa); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fine->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(j % 2) J1 = J+1; else J1 = J-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(i % 2)
		{
			// extend stencil (odd)
			I1  = I+1;
			sum = cxa[K ][J ][I ]*vsf[0]
			+     cxa[K ][J1][I ]*vsf[1]
			+     cxa[K1][J ][I ]*vsf[2]
			+     cxa[K1][J1][I ]*vsf[3]
			+     cxa[K ][J ][I1]*vsf[4]
			+     cxa[K ][J1][I1]*vsf[5]
			+     cxa[K1][J ][I1]*vsf[6]
			+     cxa[K1][J1][I1]*vsf[7];
		}
		else
		{
			// reduced stencil (even)
			sum = cxa[K ][J ][I]*vsr[0]
			+     cxa[K ][J1][I]*vsr[1]
			+     cxa[K1][J ][I]*vsr[2]
			+     cxa[K1][J1][I]*vsr[3];
		}

		// store fine grid value
		fxa[k][j][i] = sum;
	}
	END_STD_LOOP

	//---------------------
	// Y-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(fine->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(j % 2)
		{
			// extend stencil (odd)
			J1  = J+1;
			sum = cya[K ][J ][I ]*vsf[0]
			+     cya[K ][J ][I1]*vsf[1]
			+     cya[K1][J ][I ]*vsf[2]
			+     cya[K1][J ][I1]*vsf[3]
			+     cya[K ][J1][I ]*vsf[4]
			+     cya[K ][J1][I1]*vsf[5]
			+     cya[K1][J1][I ]*vsf[6]
			+     cya[K1][J1][I1]*vsf[7];
		}
		else
		{
			// reduced stencil (even)
			sum = cya[K ][J][I ]*vsr[0]
			+     cya[K ][J][I1]*vsr[1]
			+     cya[K1][J][I ]*vsr[2]
			+     cya[K1][J][I1]*vsr[3];
		}

		// store fine grid value
		fya[k][j][i] = sum;
	}
	END_STD_LOOP

	//---------------------
	// Z-points (fine grid)
	//---------------------
	ierr = DMDAGetCorners(fine->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(j % 2) J1 = J+1; else J1 = J-1;

		if(k % 2)
		{
			// extend stencil (odd)
			K1  = K+1;
			sum = cza[K ][J ][I ]*vsf[0]
			+     cza[K ][J ][I1]*vsf[1]
			+     cza[K ][J1][I ]*vsf[2]
			+     cza[K ][J1][I1]*vsf[3]
			+     cza[K1][J ][I ]*vsf[4]
			+     cza[K1][J ][I1]*vsf[5]
			+     cza[K1][J1][I ]*vsf[6]
			+     cza[K1][J1][I1]*vsf[7];

		}
		else
		{
			// reduced stencil (even)
			sum = cza[K][J ][I ]*vsr[0]
			+     cza[K][J ][I1]*vsr[1]
			+     cza[K][J1][I ]*vsr[2]
			+     cza[K][J1][I1]*vsr[3];
		}

		// store fine grid value
		fza[k][j][i] = sum;
	}
	END_STD_LOOP

	if(fine->idxmod== _IDX_COUPLED_)
	{
		//---------------------
		// P-points (fine grid)
		//---------------------
		ierr = DMDAGetCorners(fine->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get coarse grid indices
			I = i/2;
			J = j/2;
			K = k/2;

			// store fine grid value (direct injection)
			fpa[k][j][i] = cpa[K][J][I];
		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fx, &fxa);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fy, &fya);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fz, &fza);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fp, &fpa);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(coarse->fs->DA_X,   cx, &cxa); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Y,   cy, &cya); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Z,   cz, &cza); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, cp, &cpa); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
