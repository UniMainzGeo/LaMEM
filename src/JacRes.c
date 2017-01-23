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
 **    filename:   JacRes.c
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
#include "tools.h"
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
#define __FUNCT__ "JacResSetFromOptions"
PetscErrorCode JacResSetFromOptions(JacRes *jr)
{
	PetscBool   flg;
	PetscScalar gtol;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set pressure shift flag
	ierr = PetscOptionsHasName(NULL, NULL, "-skip_press_shift", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) jr->pShiftAct = PETSC_FALSE;

	ierr = PetscOptionsHasName(NULL, NULL, "-act_temp_diff", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) jr->actTemp = PETSC_TRUE;

	// set geometry tolerance
	ierr = PetscOptionsGetScalar(NULL, NULL, "-geom_tol", &gtol, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) jr->gtol = gtol;

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
	DOFIndex       *dof;
	PetscScalar    *svBuff;
	PetscInt        i, n, svBuffSz, numPhases;
	const PetscInt *lx, *ly;

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

	// velocity components
	ierr = DMCreateGlobalVector(fs->DA_X, &jr->gvx); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &jr->gvy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &jr->gvz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_X, &jr->lvx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &jr->lvy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &jr->lvz); CHKERRQ(ierr);

	// momentum residual components
	ierr = DMCreateGlobalVector(fs->DA_X, &jr->gfx); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y, &jr->gfy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z, &jr->gfz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_X, &jr->lfx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Y, &jr->lfy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_Z, &jr->lfz); CHKERRQ(ierr);

	// strain-rate components (also used as buffer vectors)
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->ldxx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->ldyy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->ldzz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_XY,  &jr->ldxy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_XZ,  &jr->ldxz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_YZ,  &jr->ldyz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XY,  &jr->gdxy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XZ,  &jr->gdxz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_YZ,  &jr->gdyz); CHKERRQ(ierr);

	// pressure
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gp); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->lp); CHKERRQ(ierr);

	// continuity residual
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gc);  CHKERRQ(ierr);

	// lithostatic pressure
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->lp_lithos); CHKERRQ(ierr);

	// pore pressure
	ierr = DMCreateLocalVector(fs->DA_CEN, &jr->lp_pore); CHKERRQ(ierr);

	// corner buffer
	ierr = DMCreateLocalVector(fs->DA_COR,  &jr->lbcor); CHKERRQ(ierr);

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

	// default geometry tolerance
	jr->gtol = 1e-15;

	// activate pressure shift
	jr->pShift    = 0.0;
	jr->pShiftAct = PETSC_TRUE;

	// switch-off temperature diffusion
	jr->actTemp = PETSC_FALSE;

	// switch off free surface tracking
	jr->AirPhase = -1;

	// change default settings
	ierr = JacResSetFromOptions(jr); CHKERRQ(ierr);

	// setup temperature parameters
	ierr = JacResCreateTempParam(jr); CHKERRQ(ierr);

	//==========================
	// 2D integration primitives
	//==========================

	// get grid partitioning in X & Y directions
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, NULL); CHKERRQ(ierr);

	// create 2D cell center grid
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_BOX,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.nproc,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, NULL, &jr->DA_CELL_2D); CHKERRQ(ierr);

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
	ierr = VecDestroy(&jr->lp);      CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gc);      CHKERRQ(ierr);

	ierr = VecDestroy(&jr->lp_lithos); CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lp_pore); CHKERRQ(ierr);

	ierr = VecDestroy(&jr->lbcor);   CHKERRQ(ierr);

	// solution variables
	ierr = PetscFree(jr->svCell);   CHKERRQ(ierr);
	ierr = PetscFree(jr->svXYEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svXZEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svYZEdge); CHKERRQ(ierr);
	ierr = PetscFree(jr->svBuff);   CHKERRQ(ierr);

	// destroy temperature parameters
	ierr = JacResDestroyTempParam(jr); CHKERRQ(ierr);

	//==========================
	// 2D integration primitives
	//==========================
	ierr = DMDestroy(&jr->DA_CELL_2D); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitScale"
PetscErrorCode JacResInitScale(JacRes *jr, UserCtx *usr)
{
	// initialize and setup scaling object, perform scaling

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize scaling object
	ierr = ScalingCreate(&jr->scal);  CHKERRQ(ierr);

	// scale input parameters
	ScalingInput(&jr->scal, usr);

	// initialize gravity acceleration
	jr->grav[0] = 0.0;
	jr->grav[1] = 0.0;
	jr->grav[2] = usr->Gravity;

	// initialize stabilization parameter
	jr->FSSA = usr->FSSA;

	// initialize time stepping parameters
	ierr = TSSolSetUp(&jr->ts, &jr->scal, usr); CHKERRQ(ierr);

	// initialize material parameter limits
	ierr = SetMatParLim(&jr->matLim, usr); CHKERRQ(ierr);

	// scale material parameter limits
	ScalingMatParLim(&jr->scal, &jr->matLim);

	// scale material parameters
	ScalingMatProp(&jr->scal, jr->phases, jr->numPhases);

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
#define __FUNCT__ "JacResGetPressShift"
PetscErrorCode JacResGetPressShift(JacRes *jr)
{
	// get average pressure near the top surface

	FDSTAG      *fs;
	PetscScalar ***p;
	PetscScalar lpShift, gpShift;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mcz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check if requested
	if(jr->pShiftAct != PETSC_TRUE) PetscFunctionReturn(0);

	fs      = jr->fs;
	mcz     = fs->dsz.tcels - 1;
	lpShift = 0.0;

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp, &p);  CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		if(k == mcz) lpShift += p[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gp, &p);  CHKERRQ(ierr);

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lpShift, &gpShift, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gpShift = lpShift;
	}

	// store pressure shift
	jr->pShift = gpShift/(PetscScalar)(fs->dsx.tcels*fs->dsy.tcels);

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
	PetscScalar gx, gy, gz, tx, ty, tz, sxx, syy, szz, sxy, sxz, syz;
	PetscScalar J2Inv, theta, rho, IKdt, Tc, pc, pShift, pn, dt, fssa, *grav;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***T, ***p_lithos, ***p_pore;
	PetscScalar eta_creep, eta_viscoplastic;
	PetscScalar depth, pc_lithos, pc_pore;
//	PetscScalar rho_lithos;
//	PetscScalar alpha, Tn,

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

//	PetscInt mcz = fs->dsz.tcels - 1;

	// initialize maximum node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// access residual context variables
	numPhases =  jr->numPhases; 	// number phases
	phases    =  jr->phases;    	// phase parameters
	matLim    = &jr->matLim;    	// phase parameters limiters
	dt        =  jr->ts.dt;     	// time step
	fssa      =  jr->FSSA;      	// density gradient penalty parameter
	grav      =  jr->grav;      	// gravity acceleration
//	rho_lithos=  matLim->rho_lithos;// density to compute lithostatic pressure in viscosity formulation
	pShift    =  jr->pShift;    	// pressure shift

	// clear local residual vectors
	ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->gc);  CHKERRQ(ierr);

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
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lithos,&p_lithos); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore,  &p_pore); CHKERRQ(ierr);

	// compute lithostatic pressure
	ierr = JacResGetLithoStaticPressure(jr); CHKERRQ(ierr);

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

		// access current lithostatic pressure
		pc_lithos = p_lithos[k][j][i];

		// access current pore pressure
		pc_pore = p_pore[k][j][i];

		// compute depth below the free surface
		depth = jr->avg_topo - COORD_CELL(k, sz, fs->dsz);

		if(depth < 0.0) depth = 0.0;

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_viscoplastic, numPhases, phases, svCell->phRat, matLim, pc_lithos, pc_pore, dt, pc-pShift, Tc); CHKERRQ(ierr);

		// store creep viscosity
		svCell->eta_creep 			= eta_creep;
		svCell->eta_viscoplastic 	= eta_viscoplastic;

		// compute stress, plastic strain rate and shear heating term on cell
		ierr = GetStressCell(svCell, matLim, XX, YY, ZZ); CHKERRQ(ierr);

		// compute total Cauchy stresses
		sxx = svCell->sxx - pc;
		syy = svCell->syy - pc;
		szz = svCell->szz - pc;



		// evaluate volumetric constitutive equations
		ierr = VolConstEq(svBulk, numPhases, phases, svCell->phRat, matLim, depth, dt, pc-pShift , Tc); CHKERRQ(ierr);

		// access
		theta = svBulk->theta; // volumetric strain rate
		rho   = svBulk->rho;   // effective density
		IKdt  = svBulk->IKdt;  // inverse bulk viscosity
//		alpha = svBulk->alpha; // effective thermal expansion
		pn    = svBulk->pn;    // pressure history
//		Tn    = svBulk->Tn;    // temperature history

		// compute gravity terms
		gx = rho*grav[0];
		gy = rho*grav[1];
		gz = rho*grav[2];

		// compute stabilization terms (lumped approximation)
		tx = -fssa*dt*gx;
		ty = -fssa*dt*gy;
		tz = -fssa*dt*gz;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
		fx[k][j][i] -= (sxx + vx[k][j][i]*tx)/bdx + gx/2.0;   fx[k][j][i+1] += (sxx + vx[k][j][i+1]*tx)/fdx - gx/2.0;
		fy[k][j][i] -= (syy + vy[k][j][i]*ty)/bdy + gy/2.0;   fy[k][j+1][i] += (syy + vy[k][j+1][i]*ty)/fdy - gy/2.0;
		fz[k][j][i] -= (szz + vz[k][j][i]*tz)/bdz + gz/2.0;   fz[k+1][j][i] += (szz + vz[k+1][j][i]*tz)/fdz - gz/2.0;

//****************************************
// ADHOC (HARD-CODED PRESSURE CONSTRAINTS)
//****************************************

//		if(k == 0)   fz[k][j][i]   += -p[k-1][j][i]/bdz;
//		if(k == mcz) fz[k+1][j][i] -= -p[k+1][j][i]/fdz;

		// mass - currently T-dependency is deactivated
//		gc[k][j][i] = -IKdt*(pc - pn) - theta + alpha*(Tc - Tn)/dt;
        
        gc[k][j][i] = -IKdt*(pc - pn) - theta ;
        
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

		// access current lithostatic pressure (x-y plane, i-j indices)
		pc_lithos = 0.25*(p_lithos[k][j][i] + p_lithos[k][j][i-1] + p_lithos[k][j-1][i] + p_lithos[k][j-1][i-1]);

		// access current pore pressure (x-y plane, i-j indices)
		pc_pore   = 0.25*(p_pore[k][j][i] + p_pore[k][j][i-1] + p_pore[k][j-1][i] + p_pore[k][j-1][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_viscoplastic, numPhases, phases, svEdge->phRat, matLim, pc_lithos, pc_pore, dt, pc-pShift, Tc); CHKERRQ(ierr);

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

		// access current lithostatic pressure (x-z plane, i-k indices)
		pc_lithos = 0.25*(p_lithos[k][j][i] + p_lithos[k][j][i-1] + p_lithos[k-1][j][i] + p_lithos[k-1][j][i-1]);

		// access current pore pressure (x-z plane, i-k indices)
		pc_pore = 0.25*(p_pore[k][j][i] + p_pore[k][j][i-1] + p_pore[k-1][j][i] + p_pore[k-1][j][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_viscoplastic, numPhases, phases, svEdge->phRat, matLim, pc_lithos, pc_pore, dt, pc-pShift, Tc); CHKERRQ(ierr);

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

		// access current lithostatic pressure (y-z plane, j-k indices)
		pc_lithos = 0.25*(p_lithos[k][j][i] + p_lithos[k][j-1][i] + p_lithos[k-1][j][i] + p_lithos[k-1][j-1][i]);

		// access current pore pressure (y-z plane, j-k indices)
		pc_pore = 0.25*(p_pore[k][j][i] + p_pore[k][j-1][i] + p_pore[k-1][j][i] + p_pore[k-1][j-1][i]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, &eta_viscoplastic, numPhases, phases, svEdge->phRat, matLim, pc_lithos, pc_pore, dt, pc-pShift, Tc); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lithos, &p_lithos); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore,   &p_pore); CHKERRQ(ierr);

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
	// copy solution from global to local vectors, enforce boundary constraints

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = JacResCopyVel (jr, x); CHKERRQ(ierr);

	ierr = JacResCopyPres(jr, x); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyVel"
PetscErrorCode JacResCopyVel(JacRes *jr, Vec x)
{
	// copy velocity from global to local vectors, enforce boundary constraints

	FDSTAG           *fs;
	BCCtx            *bc;
	PetscInt          mcx, mcy, mcz;
	PetscInt          I, J, K, fi, fj, fk;
	PetscInt          i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar       ***bcvx,  ***bcvy,  ***bcvz;
	PetscScalar       ***lvx, ***lvy, ***lvz;
	PetscScalar       *vx, *vy, *vz, pmdof;
	const PetscScalar *sol, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vectors
	ierr = VecGetArray    (jr->gvx, &vx);  CHKERRQ(ierr);
	ierr = VecGetArray    (jr->gvy, &vy);  CHKERRQ(ierr);
	ierr = VecGetArray    (jr->gvz, &vz);  CHKERRQ(ierr);
	ierr = VecGetArrayRead(x,       &sol); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = sol;

	ierr  = PetscMemcpy(vx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(vy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(vz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray    (jr->gvx, &vx);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (jr->gvy, &vy);  CHKERRQ(ierr);
	ierr = VecRestoreArray    (jr->gvz, &vz);  CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(x,       &sol); CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_X,   jr->gvx, jr->lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y,   jr->gvy, jr->lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z,   jr->gvz, jr->lvz)

	// access local solution vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	// access boundary constraints vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);

	//==============================
	// enforce two-point constraints
	//==============================

	//---------
	// X points
	//---------

	GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lvx[k][j][i];

		J = j; fj = 0;
		K = k; fk = 0;

		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcvx, lvx, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcvx, lvx, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcvx, lvx, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcvx, lvx, K, j, i, pmdof) }

		if(fj*fk) SET_EDGE_CORNER(n, lvx, K, J, i, k, j, i, pmdof)
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_NODE_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lvy[k][j][i];

		I = i; fi = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcvy, lvy, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcvy, lvy, k, j, I, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcvy, lvy, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcvy, lvy, K, j, i, pmdof) }

		if(fi*fk) SET_EDGE_CORNER(n, lvy, K, j, I, k, j, i, pmdof)
	}
	END_STD_LOOP


	//---------
	// Z points
	//---------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_NODE_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lvz[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcvz, lvz, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcvz, lvz, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcvz, lvz, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcvz, lvz, k, J, i, pmdof) }

		if(fi*fj) SET_EDGE_CORNER(n, lvz, k, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyPres"
PetscErrorCode JacResCopyPres(JacRes *jr, Vec x)
{
	// copy pressure from global to local vectors, enforce boundary constraints

	FDSTAG            *fs;
	BCCtx             *bc;
	PetscInt          mcx, mcy, mcz;
	PetscInt          I, J, K, fi, fj, fk;
	PetscInt          i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar       ***bcp;
	PetscScalar       ***lp;
	PetscScalar       *p, pmdof;
	const PetscScalar *sol, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vectors
	ierr = VecGetArray    (jr->gp, &p);   CHKERRQ(ierr);
	ierr = VecGetArrayRead(x,      &sol); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = sol + fs->nXFace + fs->nYFace + fs->nZFace;

	ierr = PetscMemcpy(p, iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray    (jr->gp, &p);   CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(x,      &sol); CHKERRQ(ierr);

	// fill local (ghosted) version of solution vectors
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gp, jr->lp)

	// access local solution vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp, &lp);  CHKERRQ(ierr);

	// access boundary constraints vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp, &bcp); CHKERRQ(ierr);

	//==============================
	// enforce two-point constraints
	//==============================

	//--------------------------
	// central points (pressure)
	//--------------------------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lp[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcp, lp, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcp, lp, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcp, lp, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcp, lp, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcp, lp, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcp, lp, K, j, i, pmdof) }

		if(fi*fj)    SET_EDGE_CORNER(n, lp, k, J, I, k, j, i, pmdof)
		if(fi*fk)    SET_EDGE_CORNER(n, lp, K, j, I, k, j, i, pmdof)
		if(fj*fk)    SET_EDGE_CORNER(n, lp, K, J, i, k, j, i, pmdof)
		if(fi*fj*fk) SET_EDGE_CORNER(n, lp, K, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp, &bcp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyRes"
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f)
{
	// copy residuals from local to global vectors, enforce boundary constraints

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
#define __FUNCT__ "JacResCopyMomentumRes"
PetscErrorCode JacResCopyMomentumRes(JacRes *jr, Vec f)
{
	// copy momentum residuals from global to local vectors for output

	FDSTAG      *fs;
	PetscScalar *fx, *fy, *fz, *res, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  = jr->fs;

	// access vectors
	ierr = VecGetArray(jr->gfx, &fx); CHKERRQ(ierr);
	ierr = VecGetArray(jr->gfy, &fy); CHKERRQ(ierr);
	ierr = VecGetArray(jr->gfz, &fz); CHKERRQ(ierr);
	ierr = VecGetArray(f, &res);      CHKERRQ(ierr);

	// copy vectors component-wise
	iter = res;

	ierr  = PetscMemcpy(fx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(fy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(fz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	// restore access
	ierr = VecRestoreArray(jr->gfx,  &fx); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gfy,  &fy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gfz,  &fz); CHKERRQ(ierr);
	ierr = VecRestoreArray(f, &res);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyContinuityRes"
PetscErrorCode JacResCopyContinuityRes(JacRes *jr, Vec f)
{
	// copy continuity residuals from global to local vectors for output

	FDSTAG      *fs;
	PetscScalar *c, *res, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  = jr->fs;

	// access vectors
	ierr = VecGetArray(jr->gc,  &c);  CHKERRQ(ierr);
	ierr = VecGetArray(f, &res);      CHKERRQ(ierr);

	// copy vectors component-wise
	iter = res + fs->dof.lnv;

	ierr = PetscMemcpy(c,  iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jr->gc,   &c);  CHKERRQ(ierr);
	ierr = VecRestoreArray(f, &res);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResViewRes"
PetscErrorCode JacResViewRes(JacRes *jr)
{
	// show assembled residual with boundary constraints
	// WARNING! rewrite this function using coupled residual vector directly

	PetscScalar dmin, dmax, d2, e2, fx, fy, fz, f2, div_tol;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get constrained residual vectors
	ierr = JacResCopyMomentumRes  (jr, jr->gres); CHKERRQ(ierr);
	ierr = JacResCopyContinuityRes(jr, jr->gres); CHKERRQ(ierr);

	// compute norms
	ierr = VecMin (jr->gc,  NULL,   &dmin); CHKERRQ(ierr);
	ierr = VecMax (jr->gc,  NULL,   &dmax); CHKERRQ(ierr);
	ierr = VecNorm(jr->gc,  NORM_2, &d2);   CHKERRQ(ierr);

	ierr = VecNorm(jr->gfx, NORM_2, &fx);   CHKERRQ(ierr);
	ierr = VecNorm(jr->gfy, NORM_2, &fy);   CHKERRQ(ierr);
	ierr = VecNorm(jr->gfz, NORM_2, &fz);   CHKERRQ(ierr);

	f2 = sqrt(fx*fx + fy*fy + fz*fz);

	if(jr->actTemp == PETSC_TRUE)
	{
    	ierr = JacResGetTempRes(jr);         CHKERRQ(ierr);
		ierr = VecNorm(jr->ge, NORM_2, &e2); CHKERRQ(ierr);
	}

	// print
	PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "Residual summary: \n");
	PetscPrintf(PETSC_COMM_WORLD, "  Continuity: \n");
	PetscPrintf(PETSC_COMM_WORLD, "    Div_min  = %12.12e \n", dmin);
	PetscPrintf(PETSC_COMM_WORLD, "    Div_max  = %12.12e \n", dmax);
	PetscPrintf(PETSC_COMM_WORLD, "    |Div|_2  = %12.12e \n", d2);
	PetscPrintf(PETSC_COMM_WORLD, "  Momentum: \n" );
	PetscPrintf(PETSC_COMM_WORLD, "    |mRes|_2 = %12.12e \n", f2);

	if(jr->actTemp == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "  Energy: \n" );
		PetscPrintf(PETSC_COMM_WORLD, "    |eRes|_2 = %12.12e \n", e2);
	}

	PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------\n");

	// stop if divergence more than tolerance
	div_tol = 0.0;
	ierr = PetscOptionsGetScalar(NULL, NULL, "-div_tol",  &div_tol,  NULL); CHKERRQ(ierr);

	if ((div_tol) && (( dmax > div_tol ) || (f2 > div_tol)))
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, " *** Emergency stop! Maximum divergence or momentum residual is too large; solver did not converge! *** \n");
	}

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
PetscErrorCode SetMatParLim(MatParLim *matLim, UserCtx *usr)
{
	// initialize material parameter limits
	PetscBool flg;
	PetscInt  cnt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	matLim->eta_min      = usr->LowerViscosityCutoff;
	matLim->eta_max      = usr->UpperViscosityCutoff;
	matLim->eta_ref      = usr->InitViscosity;

	matLim->TRef         = 0.0;
	matLim->Rugc         = 8.3144621;
	matLim->eta_atol     = 0.0;
	matLim->eta_rtol     = 1e-8;
	matLim->DII_atol     = 0.0;
	matLim->DII_rtol     = 1e-8;
	matLim->minCh        = 0.0;
	matLim->minFr        = 0.0;
	matLim->tauUlt       = DBL_MAX;
	matLim->shearHeatEff = 1.0;

	matLim->quasiHarmAvg = PETSC_FALSE;
	matLim->cf_eta_min   = 0.0;
	matLim->n_pw         = 0.0;
	matLim->MaxSNESIterBeforeApplyPlimit = 25;		

	matLim->initGuessFlg = PETSC_TRUE;
	matLim->rho_fluid    = 1040.0;
	matLim->actPorePres  = PETSC_FALSE;
	matLim->rho_lithos 	 = 0.0;	 // lithostatic density
	matLim->theta_north  = 90.0; // by default y-axis
	matLim->warn         = PETSC_TRUE;
	matLim->jac_mat_free = PETSC_FALSE;

	if(usr->DII_ref) matLim->DII_ref = usr->DII_ref;
	else
	{
		matLim->DII_ref = 1.0;
		PetscPrintf(PETSC_COMM_WORLD," WARNING: Reference strain rate DII_ref is not defined. Use a non-dimensional reference value of DII_ref =%f \n", matLim->DII_ref);
	}

	cnt = 0;

	// switch off initial guess
	ierr = PetscOptionsHasName(NULL, NULL, "-no_init_guess", &flg); CHKERRQ(ierr);
	if(flg == PETSC_TRUE) { matLim->initGuessFlg = PETSC_FALSE; }

	// plasticity stabilization parameters
	ierr = PetscOptionsHasName(NULL, NULL, "-quasi_harmonic", &flg); CHKERRQ(ierr);
	if(flg == PETSC_TRUE) { matLim->quasiHarmAvg = PETSC_TRUE; cnt++; }
	ierr = PetscOptionsGetScalar(NULL, NULL, "-cf_eta_min",  &matLim->cf_eta_min, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) { cnt++; }

	ierr = PetscOptionsGetScalar(NULL, NULL, "-n_pw",  &matLim->n_pw, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) { cnt++; }

	if(cnt > 1)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot combine plasticity stabilization methods (-quasi_harmonic -cf_eta_min -n_pw) \n");
	}

	ierr = PetscOptionsGetInt(NULL, NULL, "-MaxSNESIterBeforeApplyPlimit",  &matLim->MaxSNESIterBeforeApplyPlimit, &flg); CHKERRQ(ierr);


	// set Jacobian flag
	ierr = PetscOptionsHasName(NULL, NULL, "-jac_mat_free", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) matLim->jac_mat_free = PETSC_TRUE;

	if(cnt &&  matLim->jac_mat_free == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Analytical Jacobian is not available for plasticity stabilizations (-jac_mat_free -quasi_harmonic -cf_eta_min -n_pw) \n");
	}

	ierr = PetscOptionsGetScalar(NULL, NULL, "-rho_fluid",  &matLim->rho_fluid, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-rho_lithos", &matLim->rho_lithos, NULL); CHKERRQ(ierr);		// specify lithostatic density on commandline (if not set, we don't use this)

	// 
	ierr = PetscOptionsHasName(NULL, NULL, "-actPorePres", &flg); CHKERRQ(ierr);
	if(flg == PETSC_TRUE)
	{ 
		matLim->actPorePres=PETSC_TRUE;
	}

	ierr = PetscOptionsGetScalar(NULL, NULL, "-theta_north", &matLim->theta_north, NULL); CHKERRQ(ierr);

	ierr = PetscOptionsHasName(NULL, NULL, "-stop_warnings", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) matLim->warn = PETSC_FALSE;

	ierr = PetscOptionsGetScalar(NULL, NULL, "-shearHeatEff", &matLim->shearHeatEff, NULL); CHKERRQ(ierr);

	if(matLim->shearHeatEff > 1.0) matLim->shearHeatEff = 1.0;
	if(matLim->shearHeatEff < 0.0) matLim->shearHeatEff = 0.0;

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
	ierr = getMaxInvStep1DLocal(&fs->dsy, fs->DA_Y, jr->gvy, 1, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsz, fs->DA_Z, jr->gvz, 2, &lidtmax); CHKERRQ(ierr);

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
    
    dt = (ts->dt)*1.1;                          // slightly increase timestep
    if (dt > 1.0/gidtmax)   dt = 1.0/gidtmax;   // if dt larger than dt_courant, use courant
    if (dt > ts->dtmax)     dt = ts->dtmax;     // if dt larger than maximum dt use maximum dt
    
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
// Infinite Strain Axis (ISA) computation functions
//---------------------------------------------------------------------------
#define gradComp(v, dx, bdx1, fdx1, bdx2, fdx2, dvdx, dvdx1, dvdx2, vc) \
	dvdx  = ( v[9] - v[4])/dx; \
	dvdx1 = ((v[4] - v[0] + v[9] - v[5])/bdx1 + (v[1] - v[4] + v[6] - v[9])/fdx1)/4.0; \
	dvdx2 = ((v[4] - v[2] + v[9] - v[7])/bdx2 + (v[3] - v[4] + v[8] - v[9])/fdx2)/4.0; \
	vc    = ( v[9] + v[4])/2.0;
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getGradientVel"
PetscErrorCode getGradientVel(
	FDSTAG *fs, PetscScalar ***lvx, PetscScalar ***lvy, PetscScalar ***lvz,
	PetscInt i, PetscInt j, PetscInt k, PetscInt sx, PetscInt sy, PetscInt sz,
	Tensor2RN *L, PetscScalar *vel, PetscScalar *pvnrm)
{
	// compute velocity gradient and normalized velocities at cell center

	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar vnrm, vx[10], vy[10], vz[10], vxc, vyc, vzc;

	// get cell sizes
	dx = SIZE_CELL(i, sx, fs->dsx);   bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
	dy = SIZE_CELL(j, sy, fs->dsy);   bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
	dz = SIZE_CELL(k, sz, fs->dsz);   bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

	// vx - stencil
	vx[0] = lvx[k]  [j-1][i];
	vx[1] = lvx[k]  [j+1][i];
	vx[2] = lvx[k-1][j]  [i];
	vx[3] = lvx[k+1][j]  [i];
	vx[4] = lvx[k]  [j]  [i];
	vx[5] = lvx[k]  [j-1][i+1];
	vx[6] = lvx[k]  [j+1][i+1];
	vx[7] = lvx[k-1][j]  [i+1];
	vx[8] = lvx[k+1][j]  [i+1];
	vx[9] = lvx[k]  [j]  [i+1];

	// vy - stencil
	vy[0] = lvy[k]  [j]  [i-1];
	vy[1] = lvy[k]  [j]  [i+1];
	vy[2] = lvy[k-1][j]  [i];
	vy[3] = lvy[k+1][j]  [i];
	vy[4] = lvy[k]  [j]  [i];
	vy[5] = lvy[k]  [j+1][i-1];
	vy[6] = lvy[k]  [j+1][i+1];
	vy[7] = lvy[k-1][j+1][i];
	vy[8] = lvy[k+1][j+1][i];
	vy[9] = lvy[k]  [j+1][i];

	// vz - stencil
	vz[0] = lvz[k]  [j]  [i-1];
	vz[1] = lvz[k]  [j]  [i+1];
	vz[2] = lvz[k]  [j-1][i];
	vz[3] = lvz[k]  [j+1][i];
	vz[4] = lvz[k]  [j]  [i];
	vz[5] = lvz[k+1][j]  [i-1];
	vz[6] = lvz[k+1][j]  [i+1];
	vz[7] = lvz[k+1][j-1][i];
	vz[8] = lvz[k+1][j+1][i];
	vz[9] = lvz[k+1][j]  [i];

	gradComp(vx, dx, bdy, fdy, bdz, fdz, L->xx, L->xy, L->xz, vxc)
	gradComp(vy, dy, bdx, fdx, bdz, fdz, L->yy, L->yx, L->yz, vyc)
	gradComp(vz, dz, bdx, fdx, bdy, fdy, L->zz, L->zx, L->zy, vzc)

	// get normalized velocities
	vnrm = vxc*vxc + vyc*vyc + vzc*vzc;

	if(vnrm)
	{
		vnrm = sqrt(vnrm);

		vel[0] = vxc/vnrm;
		vel[1] = vyc/vnrm;
		vel[2] = vzc/vnrm;
	}

	if(pvnrm) (*pvnrm) = vnrm;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetISA"
PetscErrorCode JacResGetISA(JacRes *jr)
{
	// compute Infinite Strain Axis (ISA) orientation to the North

	FDSTAG      *fs;
	Tensor2RN   L;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, code, nSDFail, gnSDFail;
	PetscScalar vel[3], ISA[3];
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar ***isx, ***isy;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// set number of spectral decomposition failures
	nSDFail = 0;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &lvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &isx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &isy); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get gradient and velocities at cell center
		ierr = getGradientVel(fs, lvx, lvy, lvz, i, j, k, sx, sy, sz, &L, vel, NULL);

		// get normalized direction of Infinite Strain Axis (ISA)
		code = getISA(&L, ISA, NULL);

		if(code == 1)
		{
			// simple shear case (ISA & velocity are codirected)
			ISA[0] = vel[0];
			ISA[1] = vel[1];
		}

		// choose common sense for all ISA orientations (important for averaging)
		if(ISA[0] < 0.0 || (ISA[0] == 0.0 && ISA[1] < 0.0))
		{
			ISA[0] = -ISA[0];
			ISA[1] = -ISA[1];
		}

		// store ISA vector for output
		isx[k][j][i] = ISA[0];
		isy[k][j][i] = ISA[1];

		// report spectral decomposition failure (to adjust tolerances)
		if(code == -2)
		{
			nSDFail++;
		}

	}
	END_STD_LOOP

	// print warning
	if(jr->matLim.warn == PETSC_TRUE)
	{
		if(ISParallel(PETSC_COMM_WORLD))
		{
			MPI_Reduce(&nSDFail, &gnSDFail, 1, MPIU_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
		}
		else
		{
			gnSDFail = nSDFail;
		}

		if(gnSDFail)
		{
			PetscPrintf(PETSC_COMM_WORLD,"*****************************************************************************\n",(LLD)gnSDFail);
			PetscPrintf(PETSC_COMM_WORLD,"Warning! ISA spectral decomposition failed in %lld points. Adjust tolerances!\n",(LLD)gnSDFail);
			PetscPrintf(PETSC_COMM_WORLD,"*****************************************************************************\n",(LLD)gnSDFail);
		}
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &isx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &isy); CHKERRQ(ierr);

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetGOL"
PetscErrorCode JacResGetGOL(JacRes *jr)
{
	// compute Grain Orientation Lag (GOL) parameter

	FDSTAG      *fs;
	Tensor2RN   L;
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, code;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar lnrm, vnrm, theta, gol, max_gol, dtheta[6], vel[3], ISA[3];
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar ***ptheta, ***pgol;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// get cell index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// set maximum value of GOL parameter
	max_gol = 1.0;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &lvx);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &lvy);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &lvz);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &pgol);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &ptheta); CHKERRQ(ierr);

	//==========================================================
	// compute strain rate norm and angle between ISA & velocity
	//==========================================================

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get gradient and velocities at cell center
		ierr = getGradientVel(fs, lvx, lvy, lvz, i, j, k, sx, sy, sz, &L, vel, &vnrm);

		// get normalized direction of Infinite Strain Axis (ISA)
		code = getISA(&L, ISA, &lnrm);

		if(code == 1)
		{
			// simple shear case (ISA & velocity are codirected)
			ISA[0] = vel[0];
			ISA[1] = vel[1];
			ISA[2] = vel[2];
		}

		// get angle between ISA and velocity
		if(code < 0 || !vnrm)
		{
			// ~~~ ISA and (or) velocity is undefined ~~~

			theta = DBL_MAX;
		}
		else
		{
			// ~~~ both ISA and velocity are defined ~~~

			// get angle between ISA and velocity
			theta = ARCCOS(vel[0]*ISA[0] + vel[1]*ISA[1] + vel[2]*ISA[2]);

		}

		// store strain rate norm & angle
		pgol  [k][j][i] = lnrm;
		ptheta[k][j][i] = theta;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &ptheta); CHKERRQ(ierr);

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	//=====================================================
	// compute GOL parameter using Eulerian advection terms
	//=====================================================

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &ptheta); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get angle in current cell
		theta = ptheta[k][j][i];

		if(theta != DBL_MAX)
		{
			// get angles in neighbor cells
			if(i == 0)   { dtheta[0] = theta; } else { dtheta[0] = ptheta[k][j][i-1]; }
			if(i == mcx) { dtheta[1] = theta; } else { dtheta[1] = ptheta[k][j][i+1]; }
			if(j == 0)   { dtheta[2] = theta; } else { dtheta[2] = ptheta[k][j-1][i]; }
			if(j == mcy) { dtheta[3] = theta; } else { dtheta[3] = ptheta[k][j+1][i]; }
			if(k == 0)   { dtheta[4] = theta; } else { dtheta[4] = ptheta[k-1][j][i]; }
			if(k == mcz) { dtheta[5] = theta; } else { dtheta[5] = ptheta[k+1][j][i]; }

			// filter undefined angles
			if(dtheta[0] == DBL_MAX) dtheta[0] = M_PI_2;
			if(dtheta[1] == DBL_MAX) dtheta[1] = M_PI_2;
			if(dtheta[2] == DBL_MAX) dtheta[2] = M_PI_2;
			if(dtheta[3] == DBL_MAX) dtheta[3] = M_PI_2;
			if(dtheta[4] == DBL_MAX) dtheta[4] = M_PI_2;
			if(dtheta[5] == DBL_MAX) dtheta[5] = M_PI_2;

			// get mesh steps for the backward and forward derivatives
			bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
			bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
			bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

			// get Eulerian material time derivative
			gol = (lvx[k][j][i]*(theta - dtheta[0])/bdx + lvx[k][j][i+1]*(dtheta[1] - theta)/fdx
			+      lvy[k][j][i]*(theta - dtheta[2])/bdy + lvy[k][j+1][i]*(dtheta[3] - theta)/fdy
			+      lvz[k][j][i]*(theta - dtheta[4])/bdz + lvz[k+1][j][i]*(dtheta[5] - theta)/fdz)/2.0;

			// get strain rate norm
			lnrm = pgol[k][j][i];

			// get GOL parameter
			gol = fabs(gol)/lnrm;

			// impose limit
			if(gol > max_gol)
			{
				gol = max_gol;

			}
		}
		else
		{
			// ISA and (or) velocity is undefined, set GOL parameter to maximum
			gol = max_gol;
		}

		// store GOL parameter
		pgol[k][j][i] = gol;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &pgol);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &ptheta); CHKERRQ(ierr);

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetSHmax"
PetscErrorCode JacResGetSHmax(JacRes *jr)
{
	// compute maximum horizontal compressive stress (SHmax) orientation

	FDSTAG      *fs;
	SolVarCell  *svCell;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v1[3], v2[3], sxx, syy, sxy, s1, s2;
	PetscScalar ***dx, ***dy, ***lsxy;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// setup shear stress vector
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lsxy); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		lsxy[k][j][i] = jr->svXYEdge[iter++].s;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lsxy); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);

	// get SHmax
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &lsxy); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		svCell = &jr->svCell[iter++];
		sxx    = svCell->sxx;
		syy    = svCell->syy;
		sxy    = (lsxy[k][j][i] + lsxy[k][j][i+1] + lsxy[k][j+1][i] + lsxy[k][j+1][i+1])/4.0;

		// maximum compressive stress orientation is the eigenvector of the SMALLEST eigenvalue
		// (stress is negative in compression)
		ierr = Tensor2RS2DSpectral(sxx, syy, sxy, &s1, &s2, v1, v2, 1e-12); CHKERRQ(ierr);

		// get common sense
		if(v2[0] < 0.0 || (v2[0] == 0.0 && v2[1] < 0.0))
		{
			v2[0] = -v2[0];
			v2[1] = -v2[1];
		}

		// store direction vector for output
		dx[k][j][i] = v2[0];
		dy[k][j][i] = v2[1];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &lsxy); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetEHmax"
PetscErrorCode JacResGetEHmax(JacRes *jr)
{
	// compute maximum horizontal extension rate (EHmax) orientation

	FDSTAG      *fs;
	SolVarCell  *svCell;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v1[3], v2[3], dxx, dyy, dxy, d1, d2;
	PetscScalar ***dx, ***dy, ***ldxy;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// setup shear strain rate vector
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &ldxy); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		ldxy[k][j][i] = jr->svXYEdge[iter++].d;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &ldxy); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);

	// get EHmax
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &ldxy); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		svCell = &jr->svCell[iter++];
		dxx    = svCell->dxx;
		dyy    = svCell->dyy;
		dxy    = (ldxy[k][j][i] + ldxy[k][j][i+1] + ldxy[k][j+1][i] + ldxy[k][j+1][i+1])/4.0;

		// maximum extension rate orientation is the eigenvector of the LARGEST eigenvalue
		// (strain rate is positive in extension)
		ierr = Tensor2RS2DSpectral(dxx, dyy, dxy, &d1, &d2, v1, v2, 1e-12); CHKERRQ(ierr);

		// get common sense
		if(v1[0] < 0.0 || (v1[0] == 0.0 && v1[1] < 0.0))
		{
			v1[0] = -v1[0];
			v1[1] = -v1[1];
		}

		// store direction vector for output
		dx[k][j][i] = v1[0];
		dy[k][j][i] = v1[1];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &ldxy); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResSetVelRotation"
PetscErrorCode JacResSetVelRotation(JacRes *jr)
{
	// set velocity for the rotation benchmark

	FDSTAG      *fs;
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz;
	PetscInt    I, J, K, mcx, mcy, mcz, sx1, sz1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// set time step
	//jr->ts.dt = 0.5; // WARNING! this should ONLY be used for benchmarking with Matlab

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels-1;
	mcy = fs->dsy.tcels-1;
	mcz = fs->dsz.tcels-1;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &lvx);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &lvy);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &lvz);    CHKERRQ(ierr);

	// Vx
	GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx);
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy);
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz);

	sz1 = fs->dsz.pstart;

	START_STD_LOOP
	{
		lvx[k][j][i] = -(fs->dsz.ccoor[k-sz1]);

		if (j == 0  ) {J = j-1; lvx[k][J][i] = -(fs->dsz.ccoor[k-sz1]); }
		if (j == mcy) {J = j+1; lvx[k][J][i] = -(fs->dsz.ccoor[k-sz1]); }
		if (k == 0  ) {K = k-1; lvx[K][j][i] = -(fs->dsz.ccoor[K-sz1]); }
		if (k == mcz) {K = k+1; lvx[K][j][i] = -(fs->dsz.ccoor[K-sz1]); }

		if ((j == 0  ) && (k == 0  )) {J = j-1; K = k-1; lvx[K][J][i] = -(fs->dsz.ccoor[K-sz1]); }
		if ((j == 0  ) && (k == mcz)) {J = j-1; K = k+1; lvx[K][J][i] = -(fs->dsz.ccoor[K-sz1]); }
		if ((j == mcy) && (k == 0  )) {J = j+1; K = k-1; lvx[K][J][i] = -(fs->dsz.ccoor[K-sz1]); }
		if ((j == mcy) && (k == mcz)) {J = j+1; K = k+1; lvx[K][J][i] = -(fs->dsz.ccoor[K-sz1]); }
	}
	END_STD_LOOP

	// Vy
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx);
	GET_NODE_RANGE_GHOST_INT(ny, sy, fs->dsy);
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz);

	START_STD_LOOP
	{
		lvy[k][j][i] =  0.0;

		if (i == 0  ) {I = i-1; lvy[k][j][I] = 0.0; }
		if (i == mcx) {I = i+1; lvy[k][j][I] = 0.0; }
		if (k == 0  ) {K = k-1; lvy[K][j][i] = 0.0; }
		if (k == mcz) {K = k+1; lvy[K][j][i] = 0.0; }

		if ((i == 0  ) && (k == 0  )) {I = i-1; K = k-1; lvy[K][j][I] = 0.0; }
		if ((i == 0  ) && (k == mcz)) {I = i-1; K = k+1; lvy[K][j][I] = 0.0; }
		if ((i == mcx) && (k == 0  )) {I = i+1; K = k-1; lvy[K][j][I] = 0.0; }
		if ((i == mcx) && (k == mcz)) {I = i+1; K = k+1; lvy[K][j][I] = 0.0; }
	}
	END_STD_LOOP

	// Vz
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx);
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy);
	GET_NODE_RANGE_GHOST_INT(nz, sz, fs->dsz);

	sx1 = fs->dsx.pstart;

	START_STD_LOOP
	{
		lvz[k][j][i] = fs->dsx.ccoor[i-sx1];

		if (i == 0  ) {I = i-1; lvz[k][j][I] = fs->dsx.ccoor[I-sx1]; }
		if (i == mcx) {I = i+1; lvz[k][j][I] = fs->dsx.ccoor[I-sx1]; }
		if (j == 0  ) {J = j-1; lvz[k][J][i] = fs->dsx.ccoor[i-sx1]; }
		if (j == mcy) {J = j+1; lvz[k][J][i] = fs->dsx.ccoor[i-sx1]; }

		if ((i == 0  ) && (j == 0  )) {I = i-1; J = j-1; lvz[k][J][I] = fs->dsx.ccoor[I-sx1]; }
		if ((i == 0  ) && (j == mcy)) {I = i-1; J = j+1; lvz[k][J][I] = fs->dsx.ccoor[I-sx1]; }
		if ((i == mcx) && (j == 0  )) {I = i+1; J = j-1; lvz[k][J][I] = fs->dsx.ccoor[I-sx1]; }
		if ((i == mcx) && (j == mcy)) {I = i+1; J = j+1; lvz[k][J][I] = fs->dsx.ccoor[I-sx1]; }
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetOverPressure"
PetscErrorCode JacResGetOverPressure(JacRes *jr, Vec lop)
{
	// compute overpressure

	FDSTAG      *fs;
	PetscScalar ***op, ***p, ***p_lithos;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs  =  jr->fs;

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// access pressure vectors
	ierr = VecZeroEntries(lop); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, lop,           &op);       CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,        &p);        CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lithos, &p_lithos); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// overpressure = dynamic - lithostatic
		op[k][j][i] = p[k][j][i] - p_lithos[k][j][i];
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lop,           &op);       CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,        &p);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lithos, &p_lithos); CHKERRQ(ierr);

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, lop)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetLithoStaticPressure"
PetscErrorCode JacResGetLithoStaticPressure(JacRes *jr)
{
	// compute lithostatic pressure

	Vec         vbuff;
	FDSTAG      *fs;
	Discret1D   *dsz;
	MPI_Request srequest, rrequest;
	PetscScalar ***lp, ***ibuff, *lbuff, dz, dp, g, rho;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, iter, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs  =  jr->fs;
	dsz = &fs->dsz;
	L   =  (PetscInt)dsz->rank;
	g   =   PetscAbsScalar(jr->grav[2]);

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// get integration/communication buffer
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vbuff); CHKERRQ(ierr);

	ierr = VecZeroEntries(vbuff); CHKERRQ(ierr);

	// open index buffer for computation
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, vbuff, &ibuff); CHKERRQ(ierr);

	// open linear buffer for send/receive
	ierr = VecGetArray(vbuff, &lbuff); CHKERRQ(ierr);

	// access lithostatic pressure
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lithos, &lp); CHKERRQ(ierr);

	// start receiving integral from top domain (next)
	if(dsz->nproc != 1 && dsz->grnext != -1)
	{
		ierr = MPI_Irecv(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
	}

	// copy density
	iter = 0;

	START_STD_LOOP
	{
		lp[k][j][i] = jr->svCell[iter++].svBulk.rho;
	}
	END_STD_LOOP

	// finish receiving
	if(dsz->nproc != 1 && dsz->grnext != -1)
	{
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	}

	// compute local integral from top to bottom
	for(k = sz + nz - 1; k >= sz; k--)
	{
		START_PLANE_LOOP
		{
			// get density, cell size, pressure increment
			rho = lp[k][j][i];
			dz  = SIZE_CELL(k, sz, (*dsz));
			dp  = rho*g*dz;

			// store  lithostatic pressure
			lp[k][j][i] = ibuff[L][j][i] + dp/2.0;

			// update lithostatic pressure integral
			ibuff[L][j][i] += dp;
		}
		END_PLANE_LOOP
	}

	// send integral to bottom domain (previous)
	if(dsz->nproc != 1 && dsz->grprev != -1)
	{
		ierr = MPI_Isend(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);

		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	}

	// restore buffer and pressure vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lithos, &lp); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff, &ibuff); CHKERRQ(ierr);

	ierr = VecRestoreArray(vbuff, &lbuff); CHKERRQ(ierr);

	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff); CHKERRQ(ierr);

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_lithos)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "JacResGetPorePressure"
PetscErrorCode JacResGetPorePressure(JacRes *jr)
{
	// compute pore pressure
	FDSTAG      *fs;
	Material_t  *phases, *mat;
	PetscInt    numPhases;
	PetscScalar *phRat, rp_cv, rp;
	PetscScalar ***lp_pore, ***lp_lith, p_hydro;
	PetscInt    i,j,k,iter,iphase, sx, sy, sz, nx, ny, nz;
	PetscScalar z_top, g, gwLevel, rhow, depth;
	PetscBool   flg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Return if not activated
	if(jr->matLim.actPorePres==PETSC_FALSE)
	{
		PetscFunctionReturn(0);		
	}

	// access context
	fs        = jr->fs;
	phases    = jr->phases;
	numPhases = jr->numPhases;
	rhow      = jr->matLim.rho_fluid;
	g         = PetscAbsScalar(jr->grav[2]);

	// groundwater level
	ierr = FDSTAGGetGlobalBox(fs,NULL,NULL,NULL,NULL,NULL,&z_top); CHKERRQ(ierr);

	// set groundwater level (for example for onshore setups) 
	ierr = PetscOptionsGetScalar(NULL, NULL, "-gwLevel",  &gwLevel, &flg); CHKERRQ(ierr);
	if(flg == PETSC_TRUE) z_top = gwLevel/jr->scal.length;

	// set groundwater level equal to free surface 
	ierr = PetscOptionsHasName(NULL, NULL,"-gwLevel_eq_fs", &flg); CHKERRQ(ierr);
    if(flg == PETSC_TRUE) z_top = jr->avg_topo;

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lp_pore); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lithos, &lp_lith); CHKERRQ(ierr);

	iter = 0;
	START_STD_LOOP
	{	
		// access phase ratio array
		phRat = jr->svCell[iter++].phRat;
		
		// compute depth of the current control volume
		depth = z_top - COORD_CELL(k, sz, fs->dsz);
		if(depth < 0.0) depth = 0.0;				// we don't want these calculations in the 'air'

		// Evaluate pore pressure ratio in control volume
		rp_cv = 0.0;
		// scan all phases
		for(iphase = 0; iphase < numPhases; iphase++)
		{
			// update present phases only
			if(phRat[iphase])
			{
				// get reference to material parameters table
				mat = &phases[iphase];

				// get and check pore pressure ratio of each phase
				if(mat->rp<0.0)      mat->rp = 0.0;
				else if(mat->rp>1.0) mat->rp = 1.0;
				rp = mat->rp;

				// compute average pore pressure ratio
				rp_cv +=  phRat[iphase] * rp;
			}
		}

		// hydrostatic pressure (based on the water column)
		p_hydro = rhow * g * PetscAbsScalar(depth);

		// compute the pore pressure as product of lithostatic pressure and porepressure ratio of the control volume					
		lp_pore[k][j][i] =  p_hydro + rp_cv * (lp_lith[k][j][i]-p_hydro);		
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lp_pore); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lithos, &lp_lith); CHKERRQ(ierr);

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_pore)

	PetscFunctionReturn(0);
}
