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

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(A, (void**)&md));

	// do not add inverse viscosity term to pressure diagonal matrix (Picard operator)
	cfInvEta = 0.0;

	PetscCall(MatFreeComputeLinearOperator(md, x, f, cfInvEta));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyPreconditioner(Mat A, Vec x, Vec f)
{
	// this function corresponds to MATOP_MULT operation (f = A*x)

	MatDataPC   *mdpc;
	PetscScalar  cfInvEta;

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(A, (void**)&mdpc));

	// add inverse viscosity term to pressure diagonal matrix (preconditioner operator)
	cfInvEta = 1.0;

	PetscCall(MatFreeComputeLinearOperator(mdpc->md, x, f, cfInvEta));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeGetDiagonal(Mat A, Vec v)
{
	// this function corresponds to MATOP_GET_DIAGONAL operation (v = diag(A))

	MatDataPC *mdpc;

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(A, (void**)&mdpc));

	PetscCall(VecCopy(mdpc->D, v));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyRestrict(Mat R, Vec vf, Vec vc)
{
	// this function corresponds to MATOP_MULT operation (vc = R*vf)

	MGInterp *mgi;

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(R, (void**)&mgi));

	PetscCall(MatFreeComputeRestrict(mgi, vf, vc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeUpdateRestrict(Mat R, Vec vf, Vec vcb, Vec vc)
{
	// this function corresponds to MATOP_MULT_ADD operation (vc = vcb + R*vf)

	MGInterp *mgi;

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(R, (void**)&mgi));

	PetscCall(MatFreeComputeRestrict(mgi, vf, mgi->wc));

	// update coarse grid base vector
	if(vcb == vc)
	{
		PetscCall(VecAXPY(vc, 1.0, mgi->wc));
	}
	else
	{
		PetscCall(VecWAXPY(vc, 1.0, vcb, mgi->wc));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeApplyProlong(Mat P, Vec vc, Vec vf)
{
	// this function corresponds to MATOP_MULT operation (vf = P*vc)

	MGInterp *mgi;

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(P, (void**)&mgi));

	PetscCall(MatFreeComputeProlong(mgi, vc, vf));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeUpdateProlong(Mat P, Vec vc, Vec vfb, Vec vf)
{
	// this function corresponds to MATOP_MULT_ADD operation (vf = vfb + P*vc)

	MGInterp *mgi;

	
	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(P, (void**)&mgi));

	PetscCall(MatFreeComputeProlong(mgi, vc, mgi->wf));

	// update fine grid base vector
	if(vfb == vf)
	{
		PetscCall(VecAXPY(vf, 1.0, mgi->wf));
	}
	else
	{
		PetscCall(VecWAXPY(vf, 1.0, vfb, mgi->wf));
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

	
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// get temporary solution vectors
	PetscCall(DMGetLocalVector (fs->DA_X,   &vx));
	PetscCall(DMGetLocalVector (fs->DA_Y,   &vy));
	PetscCall(DMGetLocalVector (fs->DA_Z,   &vz));
	PetscCall(DMGetGlobalVector(fs->DA_CEN, &p));

	// get temporary residual vectors
	PetscCall(DMGetLocalVector (fs->DA_X,   &fx));
	PetscCall(DMGetLocalVector (fs->DA_Y,   &fy));
	PetscCall(DMGetLocalVector (fs->DA_Z,   &fz));
	PetscCall(DMGetGlobalVector(fs->DA_CEN, &c));

	// split solution vector
	PetscCall(MatFreeSplitVec(md, x, vx, vy, vz, p));

	// compute matrix-vector product
	PetscCall(MatFreeEvaluateLinearOperator(md, vx, vy, vz, p, fx, fy, fz, c, cfInvEta));

	// assemble residual
	PetscCall(MatFreeAssembleVec(md, f, fx, fy, fz, c));

	// restore temporary vectors
	PetscCall(DMRestoreLocalVector (fs->DA_X,   &vx));
	PetscCall(DMRestoreLocalVector (fs->DA_Y,   &vy));
	PetscCall(DMRestoreLocalVector (fs->DA_Z,   &vz));
	PetscCall(DMRestoreGlobalVector(fs->DA_CEN, &p));
	PetscCall(DMRestoreLocalVector (fs->DA_X,   &fx));
	PetscCall(DMRestoreLocalVector (fs->DA_Y,   &fy));
	PetscCall(DMRestoreLocalVector (fs->DA_Z,   &fz));
	PetscCall(DMRestoreGlobalVector(fs->DA_CEN, &c));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeComputeRestrict(MGInterp *mgi, Vec vf, Vec vc)
{
	MatData  *coarse, *fine;
	Vec       fx, fy, fz, fp;
	Vec       cx, cy, cz, cp;

	
	PetscFunctionBeginUser;

	// access context
	coarse = mgi->coarse;
	fine   = mgi->fine;

	// get temporary fine grid vectors
	PetscCall(DMGetLocalVector (fine->fs->DA_X,     &fx));
	PetscCall(DMGetLocalVector (fine->fs->DA_Y,     &fy));
	PetscCall(DMGetLocalVector (fine->fs->DA_Z,     &fz));
	PetscCall(DMGetGlobalVector(fine->fs->DA_CEN,   &fp));

	// get temporary coarse grid vectors
	PetscCall(DMGetGlobalVector(coarse->fs->DA_X,   &cx));
	PetscCall(DMGetGlobalVector(coarse->fs->DA_Y,   &cy));
	PetscCall(DMGetGlobalVector(coarse->fs->DA_Z,   &cz));
	PetscCall(DMGetGlobalVector(coarse->fs->DA_CEN, &cp));

	// split fine grid vector
	PetscCall(MatFreeSplitVec(fine, vf, fx, fy, fz, fp));

	// compute restriction to coarse grid
	PetscCall(MatFreeEvaluateRestrict(coarse, fine, fx, fy, fz, fp, cx, cy, cz, cp));

	// combine coarse grid vector
	PetscCall(MatFreeCombineVec(coarse, vc, cx, cy, cz, cp));

	// restore temporary vectors
	PetscCall(DMRestoreLocalVector (fine->fs->DA_X,     &fx));
	PetscCall(DMRestoreLocalVector (fine->fs->DA_Y,     &fy));
	PetscCall(DMRestoreLocalVector (fine->fs->DA_Z,     &fz));
	PetscCall(DMRestoreGlobalVector(fine->fs->DA_CEN,   &fp));
	PetscCall(DMRestoreGlobalVector(coarse->fs->DA_X,   &cx));
	PetscCall(DMRestoreGlobalVector(coarse->fs->DA_Y,   &cy));
	PetscCall(DMRestoreGlobalVector(coarse->fs->DA_Z,   &cz));
	PetscCall(DMRestoreGlobalVector(coarse->fs->DA_CEN, &cp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeComputeProlong(MGInterp *mgi, Vec vc, Vec vf)
{
	MatData  *coarse, *fine;
	Vec       fx, fy, fz, fp;
	Vec       cx, cy, cz, cp;

	
	PetscFunctionBeginUser;

	// access context
	coarse = mgi->coarse;
	fine   = mgi->fine;

	// get temporary fine grid vectors
	PetscCall(DMGetGlobalVector(fine->fs->DA_X,     &fx));
	PetscCall(DMGetGlobalVector(fine->fs->DA_Y,     &fy));
	PetscCall(DMGetGlobalVector(fine->fs->DA_Z,     &fz));
	PetscCall(DMGetGlobalVector(fine->fs->DA_CEN,   &fp));

	// get temporary coarse grid vectors
	PetscCall(DMGetLocalVector (coarse->fs->DA_X,   &cx));
	PetscCall(DMGetLocalVector (coarse->fs->DA_Y,   &cy));
	PetscCall(DMGetLocalVector (coarse->fs->DA_Z,   &cz));
	PetscCall(DMGetGlobalVector(coarse->fs->DA_CEN, &cp));

	// split coarse grid vector
	PetscCall(MatFreeSplitVec(coarse, vc, cx, cy, cz, cp));

	// compute prolongation on fine grid
	PetscCall(MatFreeEvaluateProlong(coarse, fine, fx, fy, fz, fp, cx, cy, cz, cp));

	// combine fine grid vector
	PetscCall(MatFreeCombineVec(fine, vf, fx, fy, fz, fp));

	// restore temporary vectors
	PetscCall(DMRestoreGlobalVector(fine->fs->DA_X,     &fx));
	PetscCall(DMRestoreGlobalVector(fine->fs->DA_Y,     &fy));
	PetscCall(DMRestoreGlobalVector(fine->fs->DA_Z,     &fz));
	PetscCall(DMRestoreGlobalVector(fine->fs->DA_CEN,   &fp));
	PetscCall(DMRestoreLocalVector (coarse->fs->DA_X,   &cx));
	PetscCall(DMRestoreLocalVector (coarse->fs->DA_Y,   &cy));
	PetscCall(DMRestoreLocalVector (coarse->fs->DA_Z,   &cz));
	PetscCall(DMRestoreGlobalVector(coarse->fs->DA_CEN, &cp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------
PetscErrorCode MatFreeComputeDiagonal(MatData *md, Vec d)
{
	FDSTAG  *fs;
	Vec      dx, dy, dz, dp;

	
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;

	// get temporary diagoanal vectors
	PetscCall(DMGetLocalVector (fs->DA_X,   &dx));
	PetscCall(DMGetLocalVector (fs->DA_Y,   &dy));
	PetscCall(DMGetLocalVector (fs->DA_Z,   &dz));
	PetscCall(DMGetGlobalVector(fs->DA_CEN, &dp));

	// compute diagonal entries
	PetscCall(MatFreeEvaluateDiagonal(md, dx, dy, dz, dp));

	// assemble diagonal
	PetscCall(MatFreeAssembleVec(md, d, dx, dy, dz, dp, 1.0));

	// restore temporary vectors
	PetscCall(DMRestoreLocalVector (fs->DA_X,   &dx));
	PetscCall(DMRestoreLocalVector (fs->DA_Y,   &dy));
	PetscCall(DMRestoreLocalVector (fs->DA_Z,   &dz));
	PetscCall(DMRestoreGlobalVector(fs->DA_CEN, &dp));

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

	
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary vectors
	PetscCall(DMGetGlobalVector(fs->DA_X, &gvx));
	PetscCall(DMGetGlobalVector(fs->DA_Y, &gvy));
	PetscCall(DMGetGlobalVector(fs->DA_Z, &gvz));

	// access temporary vectors
	PetscCall(VecGetArray    (gvx, &vx));
	PetscCall(VecGetArray    (gvy, &vy));
	PetscCall(VecGetArray    (gvz, &vz));
	PetscCall(VecGetArray    (gvp, &vp));
	PetscCall(VecGetArrayRead(v,   &va));

	// copy vectors component-wise
	iter = va;

	PetscCall(PetscMemcpy(vx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)));
	iter += fs->nXFace;

	PetscCall(PetscMemcpy(vy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)));
	iter += fs->nYFace;

	PetscCall(PetscMemcpy(vz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)));
	iter += fs->nZFace;

	PetscCall(PetscMemcpy(vp,  iter, (size_t)fs->nCells*sizeof(PetscScalar)));

	// restore access
	PetscCall(VecRestoreArray    (gvx, &vx));
	PetscCall(VecRestoreArray    (gvy, &vy));
	PetscCall(VecRestoreArray    (gvz, &vz));
	PetscCall(VecRestoreArray    (gvp, &vp));
	PetscCall(VecRestoreArrayRead(v,   &va));

	PetscCall(VecZeroEntries(lvx));
	PetscCall(VecZeroEntries(lvy));
	PetscCall(VecZeroEntries(lvz));

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_X, gvx, lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y, gvy, lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z, gvz, lvz)

	// restore temporary global vectors
	PetscCall(DMRestoreGlobalVector(fs->DA_X, &gvx));
	PetscCall(DMRestoreGlobalVector(fs->DA_Y, &gvy));
	PetscCall(DMRestoreGlobalVector(fs->DA_Z, &gvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeAssembleVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp, PetscScalar setVal)
{
	// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions

	FDSTAG *fs;
	Vec     gvx, gvy, gvz;

	
	PetscFunctionBeginUser;

	fs = md->fs;

	// get temporary global vectors
	PetscCall(DMGetGlobalVector(fs->DA_X, &gvx));
	PetscCall(DMGetGlobalVector(fs->DA_Y, &gvy));
	PetscCall(DMGetGlobalVector(fs->DA_Z, &gvz));

	// assemble ghost point contributions
	LOCAL_TO_GLOBAL(fs->DA_X, lvx, gvx)
	LOCAL_TO_GLOBAL(fs->DA_Y, lvy, gvy)
	LOCAL_TO_GLOBAL(fs->DA_Z, lvz, gvz)

	// combine vector blocks
	PetscCall(MatFreeCombineVec(md, v, gvx, gvy, gvz, gvp, setVal));

	// restore temporary global vectors
	PetscCall(DMRestoreGlobalVector(fs->DA_X, &gvx));
	PetscCall(DMRestoreGlobalVector(fs->DA_Y, &gvy));
	PetscCall(DMRestoreGlobalVector(fs->DA_Z, &gvz));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode MatFreeCombineVec(MatData *md, Vec v, Vec gvx, Vec gvy, Vec gvz, Vec gvp, PetscScalar setVal)
{
	// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions

	FDSTAG      *fs;
	PetscInt     i, num, *list;
	PetscScalar *vx, *vy, *vz, *vp, *va, *iter;

	
	PetscFunctionBeginUser;

	fs = md->fs;

	// access vectors
	PetscCall(VecGetArray(gvx, &vx));
	PetscCall(VecGetArray(gvy, &vy));
	PetscCall(VecGetArray(gvz, &vz));
	PetscCall(VecGetArray(gvp, &vp));
	PetscCall(VecGetArray(v,   &va));

	// copy vectors component-wise
	iter = va;

	PetscCall(PetscMemcpy(iter, vx, (size_t)fs->nXFace*sizeof(PetscScalar)));
	iter += fs->nXFace;

	PetscCall(PetscMemcpy(iter, vy, (size_t)fs->nYFace*sizeof(PetscScalar)));
	iter += fs->nYFace;

	PetscCall(PetscMemcpy(iter, vz, (size_t)fs->nZFace*sizeof(PetscScalar)));
	iter += fs->nZFace;

	PetscCall(PetscMemcpy(iter, vp,  (size_t)fs->nCells*sizeof(PetscScalar)));

	// zero out constrained residuals (velocity)
	num   = md->vNumSPC;
	list  = md->vSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = setVal;

	// zero out constrained residuals (pressure)
	num   = md->pNumSPC;
	list  = md->pSPCListVec;

	for(i = 0; i < num; i++) va[list[i]] = setVal;

	// restore access
	PetscCall(VecRestoreArray(gvx, &vx));
	PetscCall(VecRestoreArray(gvy, &vy));
	PetscCall(VecRestoreArray(gvz, &vz));
	PetscCall(VecRestoreArray(gvp, &vp));
	PetscCall(VecRestoreArray(v,   &va));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// low-level functions
//---------------------------------------------------------------------------
PetscErrorCode MatFreeEvaluateLinearOperator(MatData *md,
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

	
	PetscFunctionBeginUser;

	fs     = md->fs;     // grid context
	dt     = md->dt;     // time step
	fssa   = md->fssa;   // density gradient penalty parameter
	grav   = md->grav;   // gravity acceleration

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

    // clear residual components
	PetscCall(VecZeroEntries(lfx));
	PetscCall(VecZeroEntries(lfy));
	PetscCall(VecZeroEntries(lfz));

	// access work vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   lvx,  &vx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   lvy,  &vy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   lvz,  &vz));
	PetscCall(DMDAVecGetArray(fs->DA_X,   lfx,  &fx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   lfy,  &fy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   lfz,  &fz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, gc,   &c));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, gp,   &p));

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

	//-------------------------------
	// central points
	//-------------------------------
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
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   lvx,       &vx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   lvy,       &vy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   lvz,       &vz));
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   lfx,       &fx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   lfy,       &fy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   lfz,       &fz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, gc,        &c));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, gp,        &p));
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeEvaluateRestrict(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp)
{
	PetscScalar vs[12], sum;
	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fxa, ***fya, ***fza, ***fpa;
	PetscScalar ***cxa, ***cya, ***cza, ***cpa;

	
	PetscFunctionBeginUser;

	// access vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fx, &fxa));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fy, &fya));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fz, &fza));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fp, &fpa));

	// access vectors in coarse grid
	PetscCall(DMDAVecGetArray(coarse->fs->DA_X,   cx, &cxa));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Y,   cy, &cya));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Z,   cz, &cza));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, cp, &cpa));

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
	PetscCall(DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

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

	//-----------------------
	// P-points (coarse grid)
	//-----------------------

	// set pressure weights
	vs[0] = 1.0/8.0;
	vs[1] = 1.0/8.0;
	vs[2] = 1.0/8.0;
	vs[3] = 1.0/8.0;
	vs[4] = 1.0/8.0;
	vs[5] = 1.0/8.0;
	vs[6] = 1.0/8.0;
	vs[7] = 1.0/8.0;

	PetscCall(DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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

	// restore access
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fx, &fxa));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fy, &fya));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fz, &fza));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fp, &fpa));

	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_X,   cx, &cxa));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Y,   cy, &cya));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Z,   cz, &cza));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, cp, &cpa));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeEvaluateProlong(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp)
{
	PetscScalar vsf[8], vsr[4], sum;
	PetscInt    I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fxa, ***fya, ***fza, ***fpa;
	PetscScalar ***cxa, ***cya, ***cza, ***cpa;

	
	PetscFunctionBeginUser;

	// access vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fx, &fxa));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fy, &fya));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fz, &fza));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fp, &fpa));

	// access vectors in coarse grid
	PetscCall(DMDAVecGetArray(coarse->fs->DA_X,   cx, &cxa));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Y,   cy, &cya));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Z,   cz, &cza));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, cp, &cpa));

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
	PetscCall(DMDAGetCorners(fine->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fine->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

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
	PetscCall(DMDAGetCorners(fine->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

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

	//---------------------
	// P-points (fine grid)
	//---------------------
	PetscCall(DMDAGetCorners(fine->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

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


	// restore access
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fx, &fxa));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fy, &fya));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fz, &fza));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fp, &fpa));

	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_X,   cx, &cxa));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Y,   cy, &cya));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Z,   cz, &cza));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, cp, &cpa));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatFreeEvaluateDiagonal(MatData *md,
		Vec ldx, Vec ldy, Vec ldz, Vec gdp)
{
	// get diagonal of the preconditioner matrix

	FDSTAG      *fs;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, rho, Kb, IKdt, diag, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***vdx,  ***vdy,  ***vdz,  ***vdp;
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

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

    // clear diagonal components
	PetscCall(VecZeroEntries(ldx));
	PetscCall(VecZeroEntries(ldy));
	PetscCall(VecZeroEntries(ldz));

	// access work vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   ldx,  &vdx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   ldy,  &vdy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   ldz,  &vdz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, gdp,  &vdp));

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

		// get pressure diagonal element with inverse viscosity term
		diag = -IKdt - 1.0/eta;

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

		// update diagonal entries
		vdx[k  ][j  ][i  ] += v[0 ];
		vdx[k  ][j  ][i+1] += v[8 ];
		vdy[k  ][j  ][i  ] += v[16];
		vdy[k  ][j+1][i  ] += v[24];
		vdz[k  ][j  ][i  ] += v[32];
		vdz[k+1][j  ][i  ] += v[40];
		vdp[k  ][j  ][i  ]  = v[48];
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

		// compute local matrix
		//       vx_(j-1)             vx_(j)               vy_(i-1)             vy_(i)
		v[0]  =  eta/dy/bdy; v[1]  = -eta/dy/bdy; v[2]  =  eta/dx/bdy; v[3]  = -eta/dx/bdy; // fx_(j-1) [sxy]
		v[4]  = -eta/dy/fdy; v[5]  =  eta/dy/fdy; v[6]  = -eta/dx/fdy; v[7]  =  eta/dx/fdy; // fx_(j)   [sxy]
		v[8]  =  eta/dy/bdx; v[9]  = -eta/dy/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fy_(i-1) [sxy]
		v[12] = -eta/dy/fdx; v[13] =  eta/dy/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fy_(i)   [sxy]

		// set ghost point flags
		idx[0] = 0; if(j == 0)   { idx[0] = -1; }
		idx[1] = 0; if(j == mny) { idx[1] = -1; }
		idx[2] = 0; if(i == 0)   { idx[2] = -1; }
		idx[3] = 0; if(i == mnx) { idx[3] = -1; }

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		vdx[k][j-1][i  ] += v[0];
		vdx[k][j  ][i  ] += v[5];
		vdy[k][j  ][i-1] += v[10];
		vdy[k][j  ][i  ] += v[15];
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

		// compute local matrix
		//       vx_(k-1)             vx_(k)               vz_(i-1)             vz_(i)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dx/bdz; v[3]  = -eta/dx/bdz; // fx_(k-1) [sxz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dx/fdz; v[7]  =  eta/dx/fdz; // fx_(k)   [sxz]
		v[8]  =  eta/dz/bdx; v[9]  = -eta/dz/bdx; v[10] =  eta/dx/bdx; v[11] = -eta/dx/bdx; // fz_(i-1) [sxz]
		v[12] = -eta/dz/fdx; v[13] =  eta/dz/fdx; v[14] = -eta/dx/fdx; v[15] =  eta/dx/fdx; // fz_(i)   [sxz]

		// set ghost point flags
		idx[0] = 0; if(k == 0)   { idx[0] = -1; }
		idx[1] = 0; if(k == mnz) { idx[1] = -1; }
		idx[2] = 0; if(i == 0)   { idx[2] = -1; }
		idx[3] = 0; if(i == mnx) { idx[3] = -1; }

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		vdx[k-1][j][i  ] += v[0];
		vdx[k  ][j][i  ] += v[5];
		vdz[k  ][j][i-1] += v[10];
		vdz[k  ][j][i  ] += v[15];
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

		// compute local matrix
		//       vy_(k-1)             vy_(k)               vz_(j-1)             vz_(j)
		v[0]  =  eta/dz/bdz; v[1]  = -eta/dz/bdz; v[2]  =  eta/dy/bdz; v[3]  = -eta/dy/bdz; // fy_(k-1) [syz]
		v[4]  = -eta/dz/fdz; v[5]  =  eta/dz/fdz; v[6]  = -eta/dy/fdz; v[7]  =  eta/dy/fdz; // fy_(k)   [syz]
		v[8]  =  eta/dz/bdy; v[9]  = -eta/dz/bdy; v[10] =  eta/dy/bdy; v[11] = -eta/dy/bdy; // fz_(j-1) [syz]
		v[12] = -eta/dz/fdy; v[13] =  eta/dz/fdy; v[14] = -eta/dy/fdy; v[15] =  eta/dy/fdy; // fz_(j)   [syz]

		// set ghost point flags
		idx[0] = 0; if(k == 0)   { idx[0] = -1; }
		idx[1] = 0; if(k == mnz) { idx[1] = -1; }
		idx[2] = 0; if(j == 0)   { idx[2] = -1; }
		idx[3] = 0; if(j == mny) { idx[3] = -1; }

		// apply two-point constraints on the ghost nodes
		getTwoPointConstr(4, idx, pdofidx, cf);

		// constrain local matrix
		constrLocalMat(4, pdofidx, cf, v);

		// update diagonal
		vdy[k-1][j  ][i] += v[0];
		vdy[k  ][j  ][i] += v[5];
		vdz[k  ][j-1][i] += v[10];
		vdz[k  ][j  ][i] += v[15];
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   ldx,  &vdx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   ldy,  &vdy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   ldz,  &vdz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, gdp,  &vdp));

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
