//---------------------------------------------------------------------------
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "constEq.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResClear"
PetscErrorCode JacResClear(JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(jr, sizeof(JacRes)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCreate"
PetscErrorCode JacResCreate(
	JacRes   *jr,
	FDSTAG   *fs,
	BCCtx    *bc)
{
	DOFIndex    *dof;
	PetscScalar *svBuff;
	PetscInt     i, n, svBuffSz, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set external handles
	jr->fs = fs;
	jr->bc = bc;

	// set indexing object
	dof = &fs->dof;

	// phases must be initialized before calling this function
	numPhases = jr->numPhases;

	//========================
	// create solution vectors
	//========================

	// coupled solution vectors
	ierr = VecCreateMPI(PETSC_COMM_WORLD, dof->ln, PETSC_DETERMINE, &jr->gsol); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, dof->ln, PETSC_DETERMINE, &jr->gres); CHKERRQ(ierr);

	// global velocity components
	ierr = DMCreateGlobalVector(fs->DA_X, &jr->gvx); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &jr->gvy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &jr->gvz); CHKERRQ(ierr);

	// local velocity components
	ierr = DMCreateLocalVector(fs->DA_X, &jr->lvx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y, &jr->lvy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z, &jr->lvz);  CHKERRQ(ierr);

	// global momentum residual components
	ierr = DMCreateGlobalVector(fs->DA_X, &jr->gfx); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &jr->gfy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &jr->gfz); CHKERRQ(ierr);

	// local momentum residual components
	ierr = DMCreateLocalVector(fs->DA_X, &jr->lfx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y, &jr->lfy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z, &jr->lfz);  CHKERRQ(ierr);

	// local strain rate components
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->ldxx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->ldyy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->ldzz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_XY,  &jr->ldxy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_XZ,  &jr->ldxz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_YZ,  &jr->ldyz); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(fs->DA_XY,  &jr->gdxy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XZ,  &jr->gdxz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_YZ,  &jr->gdyz); CHKERRQ(ierr);

	// global pressure & temperature
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gp); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gT); CHKERRQ(ierr);

	// local pressure & temperature
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->lp); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->lT); CHKERRQ(ierr);

	// global continuity & energy residuals
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gc);  CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->ge);  CHKERRQ(ierr);

	//======================================
	// allocate space for solution variables
	//======================================

	ierr = PetscMalloc(sizeof(SolVarCell)*(size_t)fs->nCells, &jr->svCell);   CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(SolVarEdge)*(size_t)fs->nXYEdg, &jr->svXYEdge); CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(SolVarEdge)*(size_t)fs->nXZEdg, &jr->svXZEdge); CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(SolVarEdge)*(size_t)fs->nYZEdg, &jr->svYZEdge); CHKERRQ(ierr);

	ierr = PetscMemzero(jr->svCell,   sizeof(SolVarCell)*(size_t)fs->nCells); CHKERRQ(ierr);
	ierr = PetscMemzero(jr->svXYEdge, sizeof(SolVarEdge)*(size_t)fs->nXYEdg); CHKERRQ(ierr);
	ierr = PetscMemzero(jr->svXZEdge, sizeof(SolVarEdge)*(size_t)fs->nXZEdg); CHKERRQ(ierr);
	ierr = PetscMemzero(jr->svYZEdge, sizeof(SolVarEdge)*(size_t)fs->nYZEdg); CHKERRQ(ierr);

	// compute total size per processor of the solution variables storage buffer
	svBuffSz = numPhases*(fs->nCells + fs->nXYEdg + fs->nXZEdg + fs->nYZEdg);

	// allocate buffer for solution variables (phRat)
	ierr = makeScalArray(&jr->svBuff, NULL, svBuffSz);

	// setup pointers
	svBuff = jr->svBuff;

	n = fs->nCells;
	for(i = 0; i < n; i++) { jr->svCell[i].phRat   = svBuff; svBuff += numPhases; }

	n = fs->nXYEdg;
	for(i = 0; i < n; i++) { jr->svXYEdge[i].phRat = svBuff; svBuff += numPhases; }

	n = fs->nXZEdg;
	for(i = 0; i < n; i++) { jr->svXZEdge[i].phRat = svBuff; svBuff += numPhases; }

	n = fs->nYZEdg;
	for(i = 0; i < n; i++) { jr->svYZEdge[i].phRat = svBuff; svBuff += numPhases; }


	// create scatter context
//	ierr = FDSTAGCreateScatter(fs, jrctx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResDestroy"
PetscErrorCode JacResDestroy(JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// solution vectors
	ierr = VecDestroy(&jr->gsol);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gres);    CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gvx);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gvy);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gvz);     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->lvx);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lvy);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lvz);     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gfx);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gfy);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gfz);     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->lfx);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lfy);     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lfz);     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->ldxx);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldyy);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldzz);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldxy);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldxz);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldyz);    CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gdxy);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdxz);    CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdyz);    CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gp);      CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gT);      CHKERRQ(ierr);

	ierr = VecDestroy(&jr->lp);      CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lT);      CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gc);      CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ge);      CHKERRQ(ierr);

	// global to local scatter context
//	ierr = VecScatterDestroy(&jr->g2lctx); CHKERRQ(ierr);

	// solution variables
	ierr = PetscFree(jr->svCell);   CHKERRQ(ierr);
	ierr = PetscFree(jr->svXYEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svXZEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svYZEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svBuff);   CHKERRQ(ierr);

	// phase parameters
	ierr = PetscFree(jr->phases);  CHKERRQ(ierr);
	ierr = PetscFree(jr->matSoft); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetI2Gdt"
PetscErrorCode JacResGetI2Gdt(JacRes *jr)
{
	// compute average inverse elastic viscosity in the integration points
	// WARNING! this should be replaced by the effective elastic strain rates

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	PetscInt    i, n;
	PetscScalar dt;

	PetscFunctionBegin;

	fs = jr->fs;
	dt = jr->ts.dt;

	//=============
	// cell centers
	//=============
	n = fs->nCells;
	for(i = 0; i < n; i++)
	{	// access solution variables
		svCell = &jr->svCell[i];
		// compute & store inverse viscosity
		svCell->svDev.I2Gdt = GetI2Gdt(jr->numPhases, jr->phases, svCell->phRat, dt);
	}
	//===========
	// xy - edges
	//===========
	n = fs->nXYEdg;
	for(i = 0; i < n; i++)
	{	// access solution variables
		svEdge = &jr->svXYEdge[i];
		// compute & store inverse viscosity
		svEdge->svDev.I2Gdt = GetI2Gdt(jr->numPhases, jr->phases, svEdge->phRat, dt);
	}
	//===========
	// xz - edges
	//===========
	n = fs->nXZEdg;
	for(i = 0; i < n; i++)
	{	// access solution variables
		svEdge = &jr->svXZEdge[i];
		// compute & store inverse viscosity
		svEdge->svDev.I2Gdt = GetI2Gdt(jr->numPhases, jr->phases, svEdge->phRat, dt);
	}
	//===========
	// yz - edges
	//===========
	n = fs->nYZEdg;
	for(i = 0; i < n; i++)
	{	// access solution variables
		svEdge = &jr->svYZEdge[i];
		// compute & store inverse viscosity
		svEdge->svDev.I2Gdt = GetI2Gdt(jr->numPhases, jr->phases, svEdge->phRat, dt);
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetEffStrainRate"
PetscErrorCode JacResGetEffStrainRate(JacRes *jr)
{

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;

	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar dx, dy, dz, xx, yy, zz, xy, xz, yz, theta, tr;
	PetscScalar ***vx,  ***vy,  ***vz;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access local (ghosted) velocity components
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	// access global strain-rate components
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);

	//-------------------------------
	// central points (dxx, dyy, dzz)
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// compute velocity gradients
		xx = (vx[k][j][i+1] - vx[k][j][i])/dx;
		yy = (vy[k][j+1][i] - vy[k][j][i])/dy;
		zz = (vz[k+1][j][i] - vz[k][j][i])/dz;

		// compute & store volumetric strain rate
		theta = xx + yy + zz;
		svBulk->theta = theta;

		// compute & store total deviatoric strain rates

		tr  = theta/3.0;
		xx -= tr;
		yy -= tr;
		zz -= tr;

		svCell->dxx = xx;
		svCell->dyy = yy;
		svCell->dzz = zz;

		// compute & store effective deviatoric strain rates
		dxx[k][j][i] = xx + svCell->hxx*svDev->I2Gdt;
		dyy[k][j][i] = yy + svCell->hyy*svDev->I2Gdt;
		dzz[k][j][i] = zz + svCell->hzz*svDev->I2Gdt;

	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points (dxy)
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// compute velocity gradients
		dvxdy = (vx[k][j][i] - vx[k][j-1][i])/dy;
		dvydx = (vy[k][j][i] - vy[k][j][i-1])/dx;

		// compute & store total strain rate
		xy = 0.5*(dvxdy + dvydx);
		svEdge->d = xy;

		// compute & store effective deviatoric strain rate
		dxy[k][j][i] = xy + svEdge->h*svDev->I2Gdt;

	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (dxz)
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvxdz = (vx[k][j][i] - vx[k-1][j][i])/dz;
		dvzdx = (vz[k][j][i] - vz[k][j][i-1])/dx;

		// compute & store total strain rate
        xz = 0.5*(dvxdz + dvzdx);
        svEdge->d = xz;

		// compute & store effective deviatoric strain rate
		dxz[k][j][i] = xz + svEdge->h*svDev->I2Gdt;

	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (dyz)
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvydz = (vy[k][j][i] - vy[k-1][j][i])/dz;
		dvzdy = (vz[k][j][i] - vz[k][j-1][i])/dy;

		// compute & store total strain rate
		yz = 0.5*(dvydz + dvzdy);
		svEdge->d = yz;

		// compute & store effective deviatoric strain rate
		dyz[k][j][i] = yz + svEdge->h*svDev->I2Gdt;

	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);

	// communicate boundary strain-rate values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldzz);
	LOCAL_TO_LOCAL(fs->DA_XY,  jr->ldxy);
	LOCAL_TO_LOCAL(fs->DA_XZ,  jr->ldxz);
	LOCAL_TO_LOCAL(fs->DA_YZ,  jr->ldyz);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetVorticity"
PetscErrorCode JacResGetVorticity(JacRes *jr)
{
	// Compute components of the vorticity pseudo-vector
	// (instantaneous rotation rates around three coordinate axis).
	// Take care of rotation direction and sign convention.
	// Throughout LaMEM, right-handed coordinate system is assumed!

	FDSTAG     *fs;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar ***gwx, ***gwy, ***gwz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,  jr->lvx,  &lvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,  jr->lvy,  &lvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,  jr->lvz,  &lvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &gwz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &gwy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &gwx);  CHKERRQ(ierr);

	//-------------------------------
	// xy edge points (wz)
	//-------------------------------

	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvxdy = (lvx[k][j][i] - lvx[k][j-1][i])/SIZE_NODE(j, sy, fs->dsy);
		dvydx = (lvy[k][j][i] - lvy[k][j][i-1])/SIZE_NODE(i, sx, fs->dsx);

		// positive (counter-clockwise) rotation around Z axis X -> Y
		gwz[k][j][i] = dvydx - dvxdy;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (wy)
	//-------------------------------

	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvxdz = (lvx[k][j][i] - lvx[k-1][j][i])/SIZE_NODE(k, sz, fs->dsz);
		dvzdx = (lvz[k][j][i] - lvz[k][j][i-1])/SIZE_NODE(i, sx, fs->dsx);

		// positive (counter-clockwise) rotation around Y axis Z -> X
		gwy[k][j][i] = dvxdz - dvzdx;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (wx)
	//-------------------------------

	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvydz = (lvy[k][j][i] - lvy[k-1][j][i])/SIZE_NODE(k, sz, fs->dsz);
		dvzdy = (lvz[k][j][i] - lvz[k][j-1][i])/SIZE_NODE(j, sy, fs->dsy);

		// positive (counter-clockwise) rotation around X axis Y -> Z
		gwx[k][j][i] = dvzdy - dvydz;
	}
	END_STD_LOOP

	// restore velocity & strain rate component vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,  jr->lvx,  &lvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,  jr->lvy,  &lvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,  jr->lvz,  &lvz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &gwz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &gwy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &gwx);  CHKERRQ(ierr);

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);
	LOCAL_TO_LOCAL(fs->DA_XZ, jr->ldxz);
	LOCAL_TO_LOCAL(fs->DA_YZ, jr->ldyz);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetResidual"
PetscErrorCode JacResGetResidual(JacRes *jr)
{
	// Compute residual of nonlinear momentum and mass conservation
	// equations, based on pre-computed components of effective
	// strain-rate tensor, current values of pressure and temperature.
	// Missing components of the second invariant of the effective strain-rate
	// tensor (squares of the corresponding strain rate components) are averaged
	// form the hosting nodes using arithmetic mean.
	// DII = (0.5*D_ij*D_ij)^0.5
	// NOTE: we interpolate and average D_ij*D_ij terms instead of D_ij
	// WARNING: ADD DENSITY STABILIZATON TERMS

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	MatParLim  *matLim;
	PetscInt    iter, numPhases;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar gx, gy, gz, sxx, syy, szz, sxy, sxz, syz;
	PetscScalar J2Inv, theta, rho, IKdt, alpha, Tc, pc, Tn, pn, dt;
	PetscScalar ***fx,  ***fy,  ***fz,  ***gc;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***T;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

//	PetscInt mcz = fs->dsz.tcels - 1;

	// initialize maximum node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// access residual context variables
	numPhases =  jr->numPhases; // number phases
	phases    =  jr->phases;    // phase parameters
	matLim    = &jr->matLim;    // phase parameters limiters
	dt        =  jr->ts.dt;     // time step

	// clear local residual vectors
	ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		//=================
		// SECOND INVARIANT
		//=================

		// access strain rates
		XX = dxx[k][j][i];
		YY = dyy[k][j][i];
		ZZ = dzz[k][j][i];

		// x-y plane, i-j indices
		XY1 = dxy[k][j][i];
		XY2 = dxy[k][j+1][i];
		XY3 = dxy[k][j][i+1];
		XY4 = dxy[k][j+1][i+1];

		// x-z plane, i-k indices
		XZ1 = dxz[k][j][i];
		XZ2 = dxz[k+1][j][i];
		XZ3 = dxz[k][j][i+1];
		XZ4 = dxz[k+1][j][i+1];

		// y-z plane, j-k indices
		YZ1 = dyz[k][j][i];
		YZ2 = dyz[k+1][j][i];
		YZ3 = dyz[k][j+1][i];
		YZ4 = dyz[k+1][j+1][i];

		// compute second invariant
		J2Inv = 0.5*(XX*XX + YY*YY + ZZ*ZZ) +
		0.25*(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25*(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
		0.25*(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure
		pc = p[k][j][i];

		// current temperature
		Tc = T[k][j][i];

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, numPhases, phases, svCell->phRat, matLim, dt, pc, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on cell
		ierr = GetStressCell(svCell, matLim, XX, YY, ZZ); CHKERRQ(ierr);

		// compute total Cauchy stresses
		sxx = svCell->sxx - pc;
		syy = svCell->syy - pc;
		szz = svCell->szz - pc;

		// evaluate volumetric constitutive equations
		ierr = VolConstEq(svBulk, numPhases, phases, svCell->phRat, matLim, dt, pc, Tc); CHKERRQ(ierr);

		// access
		theta = svBulk->theta; // volumetric strain rate
		rho   = svBulk->rho;   // effective density
		IKdt  = svBulk->IKdt;  // inverse bulk viscosity
		alpha = svBulk->alpha; // effective thermal expansion
		pn    = svBulk->pn;    // pressure history
		Tn    = svBulk->Tn;    // temperature history

		// compute gravity terms

		// WARNING! correct signs of body forces in the entire code!

//		gx = -rho*jr->grav[0]/2.0;
//		gy = -rho*jr->grav[1]/2.0;
//		gz = -rho*jr->grav[2]/2.0;
		gx =  rho*jr->grav[0]/2.0;
		gy =  rho*jr->grav[1]/2.0;
		gz =  rho*jr->grav[2]/2.0;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
//		fx[k][j][i] += sxx/bdx + gx;   fx[k][j][i+1] -= sxx/fdx - gx;
//		fy[k][j][i] += syy/bdy + gy;   fy[k][j+1][i] -= syy/fdy - gy;
//		fz[k][j][i] += szz/bdz + gz;   fz[k+1][j][i] -= szz/fdz - gz;
		fx[k][j][i] -= sxx/bdx + gx;   fx[k][j][i+1] += sxx/fdx - gx;
		fy[k][j][i] -= syy/bdy + gy;   fy[k][j+1][i] += syy/fdy - gy;
		fz[k][j][i] -= szz/bdz + gz;   fz[k+1][j][i] += szz/fdz - gz;

//****************************************
// ADHOC (HARD-CODED PRESSURE CONSTRAINTS)
//****************************************

//		if(k == 0)   fz[k][j][i]   += -p[k-1][j][i]/bdz;
//		if(k == mcz) fz[k+1][j][i] -= -p[k+1][j][i]/fdz;

		// mass
		gc[k][j][i] = -IKdt*(pc - pn) - theta + alpha*(Tc - Tn)/dt;

	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// access strain rates
		XY = dxy[k][j][i];

		// x-y plane, i-j indices (i & j - bounded)
		XX1 = dxx[k][J1][I1];
		XX2 = dxx[k][J1][I2];
		XX3 = dxx[k][J2][I1];
		XX4 = dxx[k][J2][I2];

		// x-y plane, i-j indices (i & j - bounded)
		YY1 = dyy[k][J1][I1];
		YY2 = dyy[k][J1][I2];
		YY3 = dyy[k][J2][I1];
		YY4 = dyy[k][J2][I2];

		// x-y plane, i-j indices (i & j - bounded)
		ZZ1 = dzz[k][J1][I1];
		ZZ2 = dzz[k][J1][I2];
		ZZ3 = dzz[k][J2][I1];
		ZZ4 = dzz[k][J2][I2];

		// y-z plane j-k indices (j - bounded)
		XZ1 = dxz[k][J1][i];
		XZ2 = dxz[k+1][J1][i];
		XZ3 = dxz[k][J2][i];
		XZ4 = dxz[k+1][J2][i];

		// x-z plane i-k indices (i - bounded)
		YZ1 = dyz[k][j][I1];
		YZ2 = dyz[k+1][j][I1];
		YZ3 = dyz[k][j][I2];
		YZ4 = dyz[k+1][j][I2];

		// compute second invariant
		J2Inv = XY*XY +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
		0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (x-y plane, i-j indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k][j-1][i] + p[k][j-1][i-1]);

		// current temperature (x-y plane, i-j indices)
		Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k][j-1][i] + T[k][j-1][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, numPhases, phases, svEdge->phRat, matLim, dt, pc, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, XY); CHKERRQ(ierr);

		// access xy component of the Cauchy stress
		sxy = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// momentum
//		fx[k][j-1][i] += sxy/bdy;   fx[k][j][i] -= sxy/fdy;
//		fy[k][j][i-1] += sxy/bdx;   fy[k][j][i] -= sxy/fdx;
		fx[k][j-1][i] -= sxy/bdy;   fx[k][j][i] += sxy/fdy;
		fy[k][j][i-1] -= sxy/bdx;   fy[k][j][i] += sxy/fdx;

	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// access strain rates
		XZ = dxz[k][j][i];

		// x-z plane, i-k indices (i & k - bounded)
		XX1 = dxx[K1][j][I1];
		XX2 = dxx[K1][j][I2];
		XX3 = dxx[K2][j][I1];
		XX4 = dxx[K2][j][I2];

		// x-z plane, i-k indices (i & k - bounded)
		YY1 = dyy[K1][j][I1];
		YY2 = dyy[K1][j][I2];
		YY3 = dyy[K2][j][I1];
		YY4 = dyy[K2][j][I2];

		// x-z plane, i-k indices (i & k - bounded)
		ZZ1 = dzz[K1][j][I1];
		ZZ2 = dzz[K1][j][I2];
		ZZ3 = dzz[K2][j][I1];
		ZZ4 = dzz[K2][j][I2];

		// y-z plane, j-k indices (k - bounded)
		XY1 = dxy[K1][j][i];
		XY2 = dxy[K1][j+1][i];
		XY3 = dxy[K2][j][i];
		XY4 = dxy[K2][j+1][i];

		// xy plane, i-j indices (i - bounded)
		YZ1 = dyz[k][j][I1];
		YZ2 = dyz[k][j+1][I1];
		YZ3 = dyz[k][j][I2];
		YZ4 = dyz[k][j+1][I2];

		// compute second invariant
		J2Inv = XZ*XZ +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (x-z plane, i-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k-1][j][i] + p[k-1][j][i-1]);

		// current temperature (x-z plane, i-k indices)
		Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k-1][j][i] + T[k-1][j][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, numPhases, phases, svEdge->phRat, matLim, dt, pc, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, XZ); CHKERRQ(ierr);

		// access xz component of the Cauchy stress
		sxz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// momentum
//		fx[k-1][j][i] += sxz/bdz;   fx[k][j][i] -= sxz/fdz;
//		fz[k][j][i-1] += sxz/bdx;   fz[k][j][i] -= sxz/fdx;
		fx[k-1][j][i] -= sxz/bdz;   fx[k][j][i] += sxz/fdz;
		fz[k][j][i-1] -= sxz/bdx;   fz[k][j][i] += sxz/fdx;

	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// access strain rates
		YZ = dyz[k][j][i];

		// y-z plane, j-k indices (j & k - bounded)
		XX1 = dxx[K1][J1][i];
		XX2 = dxx[K1][J2][i];
		XX3 = dxx[K2][J1][i];
		XX4 = dxx[K2][J2][i];

		// y-z plane, j-k indices (j & k - bounded)
		YY1 = dyy[K1][J1][i];
		YY2 = dyy[K1][J2][i];
		YY3 = dyy[K2][J1][i];
		YY4 = dyy[K2][J2][i];

		// y-z plane, j-k indices (j & k - bounded)
		ZZ1 = dzz[K1][J1][i];
		ZZ2 = dzz[K1][J2][i];
		ZZ3 = dzz[K2][J1][i];
		ZZ4 = dzz[K2][J2][i];

		// x-z plane, i-k indices (k -bounded)
		XY1 = dxy[K1][j][i];
		XY2 = dxy[K1][j][i+1];
		XY3 = dxy[K2][j][i];
		XY4 = dxy[K2][j][i+1];

		// x-y plane, i-j indices (j - bounded)
		XZ1 = dxz[k][J1][i];
		XZ2 = dxz[k][J1][i+1];
		XZ3 = dxz[k][J2][i];
		XZ4 = dxz[k][J2][i+1];

		// compute second invariant
		J2Inv = YZ*YZ +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (y-z plane, j-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j-1][i] + p[k-1][j][i] + p[k-1][j-1][i]);

		// current temperature (y-z plane, j-k indices)
		Tc = 0.25*(T[k][j][i] + T[k][j-1][i] + T[k-1][j][i] + T[k-1][j-1][i]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, numPhases, phases, svEdge->phRat, matLim, dt, pc, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, YZ); CHKERRQ(ierr);

		// access yz component of the Cauchy stress
		syz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// update momentum residuals
//		fy[k-1][j][i] += syz/bdz;   fy[k][j][i] -= syz/fdz;
//		fz[k][j-1][i] += syz/bdy;   fz[k][j][i] -= syz/fdy;
		fy[k-1][j][i] -= syz/bdz;   fy[k][j][i] += syz/fdz;
		fz[k][j-1][i] -= syz/bdy;   fz[k][j][i] += syz/fdy;

	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopySol"
PetscErrorCode JacResCopySol(JacRes *jr, Vec x)
{
	FDSTAG      *fs;
	BCCtx       *bc;
	DOFIndex    *dof;
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, I, J, K, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcp;
	PetscScalar ***ivx,   ***ivy,   ***ivz,  ***ip;
	PetscScalar ***lvx, ***lvy, ***lvz, ***lp;
	PetscScalar *vx, *vy, *vz, *p, *sol, *iter, pmdof, bcval;
	PetscScalar *vals;
	PetscInt    num, *list;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;
	dof = &fs->dof;

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vectors
	ierr = VecGetArray(jr->gvx, &vx);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gvy, &vy);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gvz, &vz);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gp,  &p);   CHKERRQ(ierr);
	ierr = VecGetArray(x,       &sol); CHKERRQ(ierr);

	// enforce single point constraints (velocity)
	num   = bc->vNumSPC;
	list  = bc->vSPCList;
	vals  = bc->vSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	// enforce single point constraints (pressure)
	num   = bc->pNumSPC;
	list  = bc->pSPCList;
	vals  = bc->pSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	// copy vectors component-wise
	iter = sol;

	ierr  = PetscMemcpy(vx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(vy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(vz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr  = PetscMemcpy(p,  iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jr->gvx, &vx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gvy, &vy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gvz, &vz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gp,  &p);   CHKERRQ(ierr);
	ierr = VecRestoreArray(x, &sol);       CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_X,   jr->gvx, jr->lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y,   jr->gvy, jr->lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z,   jr->gvz, jr->lvz)
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gp,  jr->lp)

	// access local solution vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraints vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,   &bcp); CHKERRQ(ierr);

	// enforce two-point constraints

	//---------
	// X points
	//---------
	GET_NODE_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost points only
		if(ivx[k][j][i] == -1)
		{
			// get primary dof
			I = i; if(I < 0) I = 0; if(I > mnx) I = mnx;
			J = j; if(J < 0) J = 0; if(J > mcy) J = mcy;
			K = k; if(K < 0) K = 0; if(K > mcz) K = mcz;

			pmdof = lvx[K][J][I];
			bcval = bcvx[k][j][i];

			// no-gradient or prescribed value
			if(bcval == DBL_MAX) lvx[k][j][i] = pmdof;
			else                 lvx[k][j][i] = 2.0*bcval - pmdof;
		}
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_NODE_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost points only
		if(ivy[k][j][i] == -1)
		{
			// get primary dof
			I = i; if(I < 0) I = 0; if(I > mcx) I = mcx;
			J = j; if(J < 0) J = 0; if(J > mny) J = mny;
			K = k; if(K < 0) K = 0; if(K > mcz) K = mcz;

			pmdof = lvy[K][J][I];
			bcval = bcvy[k][j][i];

			// no-gradient or prescribed value
			if(bcval == DBL_MAX) lvy[k][j][i] = pmdof;
			else                 lvy[k][j][i] = 2.0*bcval - pmdof;
		}
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_NODE_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost points only
		if(ivz[k][j][i] == -1)
		{
			// get primary dof
			I = i; if(I < 0) I = 0; if(I > mcx) I = mcx;
			J = j; if(J < 0) J = 0; if(J > mcy) J = mcy;
			K = k; if(K < 0) K = 0; if(K > mnz) K = mnz;

			pmdof = lvz[K][J][I];
			bcval = bcvz[k][j][i];

			// no-gradient or prescribed value
			if(bcval == DBL_MAX) lvz[k][j][i] = pmdof;
			else                 lvz[k][j][i] = 2.0*bcval - pmdof;
		}
	}
	END_STD_LOOP

	//----------------
	// central points
	//---------------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost points only
		if(ip[k][j][i] == -1)
		{
			// get primary dof
			I = i; if(I < 0) I = 0; if(I > mcx) I = mcx;
			J = j; if(J < 0) J = 0; if(J > mcy) J = mcy;
			K = k; if(K < 0) K = 0; if(K > mcz) K = mcz;

			pmdof = lp[K][J][I];
			bcval = bcp[k][j][i];

			// no-gradient or prescribed value
			if(bcval == DBL_MAX) lp[k][j][i] = pmdof;
			else                 lp[k][j][i] = 2.0*bcval - pmdof;

		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyRes"
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f)
{
	FDSTAG      *fs;
	BCCtx       *bc;
	PetscInt    i, num, *list;
	PetscScalar *fx, *fy, *fz, *c, *res, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  = jr->fs;
	bc  = jr->bc;

	// access vectors
	ierr = VecGetArray(jr->gfx, &fx); CHKERRQ(ierr);
	ierr = VecGetArray(jr->gfy, &fy); CHKERRQ(ierr);
	ierr = VecGetArray(jr->gfz, &fz); CHKERRQ(ierr);
	ierr = VecGetArray(jr->gc,  &c);  CHKERRQ(ierr);
	ierr = VecGetArray(f, &res);      CHKERRQ(ierr);

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
	num   = bc->vNumSPC;
	list  = bc->vSPCList;

	for(i = 0; i < num; i++) res[list[i]] = 0.0;

	// zero out constrained residuals (pressure)
	num   = bc->pNumSPC;
	list  = bc->pSPCList;

	for(i = 0; i < num; i++) res[list[i]] = 0.0;

	// restore access
	ierr = VecRestoreArray(jr->gfx,  &fx); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gfy,  &fy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gfz,  &fz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gc,   &c);  CHKERRQ(ierr);
	ierr = VecRestoreArray(f, &res);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResViewRes"
PetscErrorCode JacResViewRes(JacRes *jr)
{
	// #define MAX3(a,b,c) (a > b ? (a > c ? a : c) : (b > c ? b : c))

	// WARNING! show assembled residual (with bc)

	PetscBool   flg;
	PetscScalar dmin, dmax, d2, fx, fy, fz, f2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// view residuals if required
	ierr = PetscOptionsHasName(NULL, "-res_log", &flg); CHKERRQ(ierr);

	if(flg != PETSC_TRUE) PetscFunctionReturn(0);

	// compute norms
	ierr = VecMin (jr->gc, NULL,   &dmin); CHKERRQ(ierr);
	ierr = VecMax (jr->gc, NULL,   &dmax); CHKERRQ(ierr);
	ierr = VecNorm(jr->gc, NORM_2, &d2);   CHKERRQ(ierr);

	ierr = VecNorm(jr->gfx, NORM_2, &fx);  CHKERRQ(ierr);
	ierr = VecNorm(jr->gfy, NORM_2, &fy);  CHKERRQ(ierr);
	ierr = VecNorm(jr->gfz, NORM_2, &fz);  CHKERRQ(ierr);

	f2 = sqrt(fx*fx + fy*fy + fz*fz);

	// print
	PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "Residual summary: \n");
	PetscPrintf(PETSC_COMM_WORLD, "  Continuity: \n");
	PetscPrintf(PETSC_COMM_WORLD, "    Div_min  = %12.12e \n", dmin);
	PetscPrintf(PETSC_COMM_WORLD, "    Div_max  = %12.12e \n", dmax);
	PetscPrintf(PETSC_COMM_WORLD, "    |Div|_2  = %12.12e \n", d2);
	PetscPrintf(PETSC_COMM_WORLD, "  Momentum: \n" );
	PetscPrintf(PETSC_COMM_WORLD, "    |mRes|_2 = %12.12e \n", f2);
	PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar JacResGetTime(JacRes *jr)
{
	return	jr->ts.time*jr->scal.time;
}
//---------------------------------------------------------------------------
PetscInt JacResGetStep(JacRes *jr)
{
	return	jr->ts.istep;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetMatParLim"
PetscErrorCode SetMatParLim(MatParLim *matLim, UserContext *usr)
{
	// initialize material parameter limits

	PetscFunctionBegin;

	matLim->eta_min      = usr->LowerViscosityCutoff;
	matLim->eta_max      = usr->UpperViscosityCutoff;
	matLim->TRef         = 0.0;
	matLim->Rugc         = usr->GasConstant;
	matLim->eta_atol     = 0.0;
	matLim->eta_rtol     = 1e-8;
	matLim->DII_atol     = 0.0;
	matLim->DII_rtol     = 1e-8;
	matLim->DII_ref      = 1.0;
	matLim->minCh        = 0.0;
	matLim->minFr        = 0.0;
	matLim->tauUlt       = DBL_MAX;
	matLim->shearHeatEff = 1.0;
	matLim->quasiHarmAvg = PETSC_FALSE;
	matLim->initGuessFlg = PETSC_TRUE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetCourantStep"
PetscErrorCode JacResGetCourantStep(JacRes *jr)
{
	//-------------------------------------
	// compute length of the next time step
	//-------------------------------------

	FDSTAG      *fs;
	TSSol       *ts;
	PetscScalar dt, lidtmax, gidtmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs =  jr->fs;
	ts = &jr->ts;

	lidtmax = 0.0;

	// determine maximum local inverse time step
	ierr = getMaxInvStep1DLocal(&fs->dsx, fs->DA_X, jr->gvx, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsy, fs->DA_Y, jr->gvy, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsz, fs->DA_Z, jr->gvz, 0, &lidtmax); CHKERRQ(ierr);

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lidtmax, &gidtmax, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gidtmax = lidtmax;
	}

	// compute time step
	gidtmax /= ts->Cmax;

	if(gidtmax < 1.0/ts->dtmax) dt = ts->dtmax;
	else                        dt = 1.0/gidtmax;

	// store new time step
	ts->pdt = ts->dt;
	ts->dt  = dt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getMaxInvStep1DLocal"
PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax)
{
	PetscScalar v, h, vmax, idt, idtmax;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, idx, ijk[3], jj, ln;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	idtmax = (*_idtmax);

	if(ds->h_uni < 0.0)
	{
		// compute time step on variable spacing grid
		PetscScalar ***va;

		ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, gv, &va);                     CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get velocity
			v = va[k][j][i];

			// prepare node index buffer
			ijk[0] = i-sx;
			ijk[1] = j-sy;
			ijk[2] = k-sz;

			// anisotropic direction-dependent criterion
			if(v >= 0.0)  idx = ijk[dir];
			else          idx = ijk[dir]-1;

			// get mesh step
			h = ds->ncoor[idx+1] - ds->ncoor[idx];

			// get inverse time step (safe to compute)
			idt = v/h;

			// update maximum inverse time step
			if(idt > idtmax) idtmax = idt;
		}
		END_STD_LOOP

		ierr = DMDAVecRestoreArray(da, gv, &va); CHKERRQ(ierr);
	}
	else
	{
		// compute time step on uniform spacing grid
		PetscScalar *va;

		// get maximum local velocity
		ierr = VecGetLocalSize(gv, &ln); CHKERRQ(ierr);
		ierr = VecGetArray(gv, &va);     CHKERRQ(ierr);

		vmax = 0.0;
		for(jj = 0; jj < ln; jj++) { v = PetscAbsScalar(va[jj]); if(v > vmax) vmax = v;	}

		ierr = VecRestoreArray(gv, &va); CHKERRQ(ierr);

		// get inverse time step
		idt = vmax/ds->h_uni;

		// update maximum inverse time step
		if(idt > idtmax) idtmax = idt;
	}

	// return result
	(*_idtmax) = idtmax;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitTemp"
PetscErrorCode JacResInitTemp(JacRes *jr)
{
	PetscScalar *T, ***lT;
	FDSTAG      *fs;
	PetscInt     jj, n, bc;
	PetscInt     mcx, mcy, mcz;
	PetscInt     i, j, k, I, J, K, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// copy temperatures from context storage to global vector
	ierr = VecGetArray(jr->gT, &T);  CHKERRQ(ierr);

	for(jj = 0, n = fs->nCells; jj < n; jj++) T[jj] = jr->svCell[jj].svBulk.Tn;

	ierr = VecRestoreArray(jr->gT, &T);  CHKERRQ(ierr);

	// scatter temperature to inter-processor ghost points
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gT, jr->lT)

	// set temperature in boundary ghost points
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		bc = 0;

		I = i; if(I < 0) { I = 0; bc = 1; } if(I > mcx) { I = mcx; bc = 1; }
		J = j; if(J < 0) { J = 0; bc = 1; } if(J > mcy) { J = mcy; bc = 1; }
		K = k; if(K < 0) { K = 0; bc = 1; } if(K > mcz) { K = mcz; bc = 1; }

		if(bc) lT[k][j][i] = lT[K][J][I];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGScatterSol"
PetscErrorCode FDSTAGScatterSol(FDSTAG *fs, JacResCtx *jrctx)
{
	// scatter solution from coupled vector to component vectors, enforce constraints

	PetscScalar *TPCVals, *TPCLinComPar;
	PetscScalar *vx, *vy, *vz, *p, *sol, *iter;
	PetscInt     i, numTPC, *TPCList, *TPCPrimeDOF;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// scatter solution vector
	ierr = VecScatterBegin(jrctx->g2lctx, jrctx->gsol, jrctx->lsol, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd  (jrctx->g2lctx, jrctx->gsol, jrctx->lsol, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

	// access 1D images of local vectors
	ierr = VecGetArray(jrctx->lvx,  &vx);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lvy,  &vy);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lvz,  &vz);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gp,   &p);   CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lsol, &sol); CHKERRQ(ierr);

	// enforce two-point constraints
	numTPC       = jrctx->numTPC;       // number of two-point constraints (TPC)
	TPCList      = jrctx->TPCList;      // local indices of TPC (ghosted layout)
	TPCPrimeDOF  = jrctx->TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
	TPCVals      = jrctx->TPCVals;      // values of TPC
	TPCLinComPar = jrctx->TPCLinComPar; // linear combination parameters

	for(i = 0; i < numTPC; i++) sol[TPCList[i]] = TPCLinComPar[i]*sol[TPCPrimeDOF[i]] + TPCVals[i];

	// copy solution components from coupled storage into local vectors
	iter = sol;

	ierr  = PetscMemcpy(vx, iter, (size_t)fs->nXFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFaceGh;

	ierr  = PetscMemcpy(vy, iter, (size_t)fs->nYFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFaceGh;

	ierr  = PetscMemcpy(vz, iter, (size_t)fs->nZFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFaceGh;

	ierr  = PetscMemcpy(p,  iter, (size_t)fs->nCellsGh*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jrctx->lvx,  &vx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lvy,  &vy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lvz,  &vz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gp,   &p);   CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lsol, &sol); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGScatterRes"
PetscErrorCode FDSTAGScatterRes(FDSTAG *fs, JacResCtx *jrctx)
{
	// copy residuals from component vectors to coupled vector, enforce constraints
	PetscInt     i, n, shift, *setList;
	PetscScalar *fx, *fy, *fz, *c, *res, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access 1D images of local vectors
	ierr = VecGetArray(jrctx->lfx,  &fx);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lfy,  &fy);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lfz,  &fz);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gc,   &c);   CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->lres, &res); CHKERRQ(ierr);

	// copy local vectors into coupled storage component-wise
	iter = res;

	ierr  = PetscMemcpy(iter, fx, (size_t)fs->nXFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFaceGh;

	ierr  = PetscMemcpy(iter, fy, (size_t)fs->nYFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFaceGh;

	ierr  = PetscMemcpy(iter, fz, (size_t)fs->nZFaceGh*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFaceGh;

	ierr  = PetscMemcpy(iter, c,  (size_t)fs->nCellsGh*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jrctx->lfx,  &fx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lfy,  &fy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lfz,  &fz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gc,   &c);   CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->lres, &res); CHKERRQ(ierr);

	// assemble global residual
	ierr = VecZeroEntries(jrctx->gres); CHKERRQ(ierr);
	ierr = VecScatterBegin(jrctx->g2lctx, jrctx->lres, jrctx->gres, ADD_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecScatterEnd  (jrctx->g2lctx, jrctx->lres, jrctx->gres, ADD_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

	// zero out constrained residuals
	n       = jrctx->numSPC;  // number of single point constraints (SPC)
	setList = jrctx->SPCList; // list of constrained DOF global IDs
	shift   = fs->istart;     // global index of the first DOF (global to local conversion)

	ierr = VecGetArray(jrctx->gres, &res); CHKERRQ(ierr);

	for(i = 0; i < n; i++) res[setList[i] - shift] = 0.0;

	ierr = VecRestoreArray(jrctx->gres, &res); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCreateScatter"
PetscErrorCode FDSTAGCreateScatter(FDSTAG *fs, JacResCtx *jrctx)
{

	// For every point in the ghosted (local) layout, create a scatter context
	// to retrieve its value from the non-ghosted (global) layout.
	// Same scatter context can be used to perform residual assembly operation.
	// Boundary ghost points are not scattered, as their values are (should be) globally accessible.

	// NOTE! Residual assembly operation is much sparser
	// It's therefore better to create another scatter context for residual assembly

	IS          gIS;
	IS          lIS;
	PetscInt    *gidx;
	PetscInt    *lidx;
	PetscInt    ln, sum, ind, start, num;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate index arrays
	ierr = makeIntArray(&gidx, NULL, fs->numdofGh); CHKERRQ(ierr);
	ierr = makeIntArray(&lidx, NULL, fs->numdofGh); CHKERRQ(ierr);

	// compute starting index in the ghosted vector storage
	ln    = fs->numdofGh;
	ierr  = MPI_Scan(&ln, &sum, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	start = sum - ln;

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   fs->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   fs->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   fs->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, fs->ip,   &ip);   CHKERRQ(ierr);

	// Collect global indices of all the local & ghost points (no boundary)
	// Also store their global indices in the ghosted vector storage.

	num = 0;

	//---------
	// X-points
	//---------
	GET_NODE_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ind = (PetscInt)ivx[k][j][i]; CHECK_DOF_INTERNAL(ind, start, num, gidx, lidx);
		start++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_NODE_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ind = (PetscInt)ivy[k][j][i]; CHECK_DOF_INTERNAL(ind, start, num, gidx, lidx);
		start++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_NODE_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ind = (PetscInt)ivz[k][j][i]; CHECK_DOF_INTERNAL(ind, start, num, gidx, lidx);
		start++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ind = (PetscInt)ip[k][j][i]; CHECK_DOF_INTERNAL(ind, start, num, gidx, lidx);
		start++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   fs->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   fs->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   fs->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, fs->ip,   &ip);   CHKERRQ(ierr);

	// create index sets to scatter from global to local vector
	ierr = ISCreateGeneral(PETSC_COMM_WORLD, num, gidx, PETSC_USE_POINTER, &gIS); CHKERRQ(ierr);
	ierr = ISCreateGeneral(PETSC_COMM_WORLD, num, lidx, PETSC_USE_POINTER, &lIS); CHKERRQ(ierr);

	// create scatter object
	ierr = VecScatterCreate(jrctx->gsol, gIS, jrctx->lsol, lIS, &jrctx->g2lctx); CHKERRQ(ierr);

	// destroy index sets & arrays
	ierr = ISDestroy(&gIS); CHKERRQ(ierr);
	ierr = ISDestroy(&lIS); CHKERRQ(ierr);
	ierr = PetscFree(gidx); CHKERRQ(ierr);
	ierr = PetscFree(lidx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCopySol"
PetscErrorCode FDSTAGCopySol(FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx, Vec x)
{
	// NOTE! Properly setting the non-local boundary ghost points is ONLY
	// necessary for output. Shear stress\tangential velocity boundary conditions
	// should be implemented at the level of the edge shear stress points.
	// As an ad-hoc solution all boundary velocities are set here.


	PetscInt    mcx, mcy, mcz;
	PetscScalar dx, dz, dy;
	PetscInt    ix[2], iy[2], iz[2];
	PetscScalar bx[2], by[2], bz[2];
	PetscScalar vx[2], vy[2], vz[2];
	PetscInt    i, j, k, I, J, K, nx, ny, nz, sx, sy, sz, isConsNode;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcp;
	PetscScalar ***ivx,   ***ivy,   ***ivz,  ***ip;
	PetscScalar ***lvx, ***lvy, ***lvz, ***lp;
	PetscScalar *gvx, *gvy, *gvz, *gp, *sol, *iter, pmdof, bcval;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vectors
	ierr = VecGetArray(jrctx->gvx,  &gvx);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gvy,  &gvy);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gvz,  &gvz);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gp,   &gp);   CHKERRQ(ierr);
	ierr = VecGetArray(x,           &sol);  CHKERRQ(ierr);

	// copy vectors component-wise
	iter = sol;

	ierr  = PetscMemcpy(gvx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(gvy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(gvz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr  = PetscMemcpy(gp,  iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jrctx->gvx, &gvx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gvy, &gvy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gvz, &gvz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gp,  &gp);   CHKERRQ(ierr);
	ierr = VecRestoreArray(x, &sol);           CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	ierr = VecZeroEntries(jrctx->lvx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->lvy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->lvz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->lp);  CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(fs->DA_X,   jrctx->gvx, jrctx->lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y,   jrctx->gvy, jrctx->lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z,   jrctx->gvz, jrctx->lvz)
	GLOBAL_TO_LOCAL(fs->DA_CEN, jrctx->gp,  jrctx->lp)

	// access local solution vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jrctx->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jrctx->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jrctx->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jrctx->lp,  &lp);  CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   fs->dofcoupl.ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   fs->dofcoupl.ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   fs->dofcoupl.ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, fs->dofcoupl.ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraints vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	//=================
	// set ghost points
	//=================

	//-------------------------------
	// xy edge points (dxy)
	//-------------------------------

	GET_NODE_RANGE          (nx, sx, fs->dsx)
	GET_NODE_RANGE          (ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost point flags
		ix[0] = (PetscInt)ivx[k][j-1][i];
		ix[1] = (PetscInt)ivx[k][j][i];
	    iy[0] = (PetscInt)ivy[k][j][i-1];
		iy[1] = (PetscInt)ivy[k][j][i];

		// boundary velocities
		bx[0] = bcvx[k][j-1][i];
		bx[1] = bcvx[k][j][i];
		by[0] = bcvy[k][j][i-1];
		by[1] = bcvy[k][j][i];

		// internal velocities
		vx[0] = lvx[k][j-1][i];
		vx[1] = lvx[k][j][i];
		vy[0] = lvy[k][j][i-1];
		vy[1] = lvy[k][j][i];

		// mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// constrain node
		isConsNode = constrEdgeNode(ix, iy, bx, by, vx, vy, dx, dy);

		// copy constrained node stencil
		if(isConsNode)
		{
			lvx[k][j-1][i] = vx[0];
			lvx[k][j][i]   = vx[1];
			lvy[k][j][i-1] = vy[0];
			lvy[k][j][i]   = vy[1];
		}
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (dxz)
	//-------------------------------
	GET_NODE_RANGE          (nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_NODE_RANGE          (nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost point flags
		ix[0] = (PetscInt)ivx[k-1][j][i];
		ix[1] = (PetscInt)ivx[k][j][i];
		iz[0] = (PetscInt)ivz[k][j][i-1];
		iz[1] = (PetscInt)ivz[k][j][i];

		// boundary velocities
		bx[0] = bcvx[k-1][j][i];
		bx[1] = bcvx[k][j][i];
		bz[0] = bcvz[k][j][i-1];
		bz[1] = bcvz[k][j][i];

		// internal velocities
		vx[0] = lvx[k-1][j][i];
		vx[1] = lvx[k][j][i];
		vz[0] = lvz[k][j][i-1];
		vz[1] = lvz[k][j][i];

		// mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// constrain node
		isConsNode = constrEdgeNode(ix, iz, bx, bz, vx, vz, dx, dz);

		// copy constrained node stencil
		if(isConsNode)
		{
			lvx[k-1][j][i] = vx[0];
			lvx[k][j][i]   = vx[1];
			lvz[k][j][i-1] = vz[0];
			lvz[k][j][i]   = vz[1];
		}
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (dyz)
	//-------------------------------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_NODE_RANGE          (ny, sy, fs->dsy)
	GET_NODE_RANGE          (nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost point flags
		iy[0] = (PetscInt)ivy[k-1][j][i];
		iy[1] = (PetscInt)ivy[k][j][i];
		iz[0] = (PetscInt)ivz[k][j-1][i];
		iz[1] = (PetscInt)ivz[k][j][i];

		// boundary velocities
		by[0] = bcvy[k-1][j][i];
		by[1] = bcvy[k][j][i];
		bz[0] = bcvz[k][j-1][i];
		bz[1] = bcvz[k][j][i];

		// internal velocities
		vy[0] = lvy[k-1][j][i];
		vy[1] = lvy[k][j][i];
		vz[0] = lvz[k][j-1][i];
		vz[1] = lvz[k][j][i];

		// mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// constrain node
		isConsNode = constrEdgeNode(iy, iz, by, bz, vy, vz, dy, dz);

		// copy constrained node stencil
		if(isConsNode)
		{
			lvy[k-1][j][i] = vy[0];
			lvy[k][j][i]   = vy[1];
			lvz[k][j-1][i] = vz[0];
			lvz[k][j][i]   = vz[1];
		}
	}
	END_STD_LOOP

	//----------------
	// central points
	//---------------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// ghost points only
		if(ip[k][j][i] == -1)
		{
			// get primary dof
			I = i; if(I < 0) I = 0; if(I > mcx) I = mcx;
			J = j; if(J < 0) J = 0; if(J > mcy) J = mcy;
			K = k; if(K < 0) K = 0; if(K > mcz) K = mcz;

			pmdof = lp[K][J][I];
			bcval = bcp[k][j][i];

			// no-gradient or prescribed value
			if(bcval == DBL_MAX) lp[k][j][i] = pmdof;
			else                 lp[k][j][i] = 2.0*bcval - pmdof;

		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jrctx->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jrctx->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jrctx->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jrctx->lp,  &lp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   fs->dofcoupl.ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   fs->dofcoupl.ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   fs->dofcoupl.ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, fs->dofcoupl.ip,  &ip);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

 static inline PetscInt constrEdgeNode(
	PetscInt    ix[],
	PetscInt    iy[],
	PetscScalar bx[],
	PetscScalar by[],
	PetscScalar vx[],
	PetscScalar vy[],
	PetscScalar dx,
	PetscScalar dy)
{
	PetscInt    cx = -1, cy = -1, px = -1, py = -1;
	PetscScalar epsb = 0.0, eps = 0.0, cfx = 0.0, cfy = 0.0;

	// epsb - is a boundary value for the strain rate component.
	// For the free surface it should be zero.
	// Alternatively, if a nonzero boundary stress is required,
	// this strain rate should be computed from the constitutive model
	// (probably nonlinear).
	// Currently only the free surface is implemented, so epsb is set to zero.

	// eps - is a cumulative strain rate from internal and Dirichlet ghost nodes.
	// It is used for enforcing the Neumann constraints.

	// cfx, cfy - are the strain-rate terms pre-multipliers used in the expansion
	// of Neumann constraints. They assume values of -1 or 1, depending on which
	// node along each direction is constrained, first or second (correspondingly).

	// determine primary internal nodes & constrained ghost nodes
	if(ix[0] == -1) { cx = 0; px = 1; cfx = -1.0; }
	if(ix[1] == -1) { cx = 1; px = 0; cfx =  1.0; }
	if(iy[0] == -1) { cy = 0; py = 1; cfy = -1.0; }
	if(iy[1] == -1) { cy = 1; py = 0; cfy =  1.0; }

	// sort out the internal nodes
	if(cx == -1 && cy == -1) return 0;

	// compute strain rates from internal nodes
	if(cx == -1) eps += (vx[1] - vx[0])/dy/2.0;
	if(cy == -1) eps += (vy[1] - vy[0])/dx/2.0;

	// expand Dirichlet constraints, update strain rates from Dirichlet nodes
	if(cx != -1 && bx[cx] != DBL_MAX) { vx[cx] = 2.0*bx[cx] - vx[px]; eps += (vx[1] - vx[0])/dy/2.0; }
	if(cy != -1 && by[cy] != DBL_MAX) { vy[cy] = 2.0*by[cy] - vy[py]; eps += (vy[1] - vy[0])/dx/2.0; }

	// expand Neumann constraints
	if(cx != -1 && bx[cx] == DBL_MAX) { vx[cx] = vx[px] + 2.0*dy*cfx*(epsb - eps); }
	if(cy != -1 && by[cy] == DBL_MAX) { vy[cy] = vy[py] + 2.0*dx*cfy*(epsb - eps); }

	return 1;
}
//---------------------------------------------------------------------------
*/
