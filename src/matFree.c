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
 **    filename:   matFree.c
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
//...............   MATRIX-FREE JACOBIAN AND PRECONDITIONER  ................
//---------------------------------------------------------------------------

#include "LaMEM.h"
#include "matFree.h"
#include "fdstag.h"
#include "JacRes.h"
#include "tssolve.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyPicard"
PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec y)
{
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, (void**)&jr); CHKERRQ(ierr);

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x); CHKERRQ(ierr);

	ierr = JacResPicardMatFree(jr); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = JacResCopyRes(jr, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResPicardMatFree"
PetscErrorCode JacResPicardMatFree(JacRes *jr)
{
	FDSTAG     *fs;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar dx, dy, dz, tx, ty, tz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar sxx, syy, szz, sxy, sxz, syz;
	PetscScalar dxx, dyy, dzz, dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar eta, theta, tr, rho, IKdt, pc, dt, fssa, *grav;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc, ***p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs   = jr->fs;
	dt   = jr->ts->dt;      // time step
	fssa = jr->ctrl.FSSA;  // density gradient penalty parameter
    grav = jr->ctrl.grav;  // gravity acceleration

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
		tr  = theta/3.0;
		dxx -= tr;
		dyy -= tr;
		dzz -= tr;

		// access current pressure
		pc = p[k][j][i];

		// compute total Cauchy stresses
		sxx = 2.0*eta*dxx - pc;
		syy = 2.0*eta*dyy - pc;
		szz = 2.0*eta*dzz - pc;

		// compute stabilization terms (lumped approximation)
		tx = -fssa*dt*(rho*grav[0]);
		ty = -fssa*dt*(rho*grav[1]);
		tz = -fssa*dt*(rho*grav[2]);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
		fx[k][j][i] -= (sxx + vx[k][j][i]*tx)/bdx;   fx[k][j][i+1] += (sxx + vx[k][j][i+1]*tx)/fdx;
		fy[k][j][i] -= (syy + vy[k][j][i]*ty)/bdy;   fy[k][j+1][i] += (syy + vy[k][j+1][i]*ty)/fdy;
		fz[k][j][i] -= (szz + vz[k][j][i]*tz)/bdz;   fz[k+1][j][i] += (szz + vz[k+1][j][i]*tz)/fdz;

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
		// get effective viscosity
		eta = jr->svXYEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// compute velocity gradients
		dvxdy = (vx[k][j][i] - vx[k][j-1][i])/dy;
		dvydx = (vy[k][j][i] - vy[k][j][i-1])/dx;

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
		// get effective viscosity
		eta = jr->svXZEdge[iter++].svDev.eta;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvxdz = (vx[k][j][i] - vx[k-1][j][i])/dz;
		dvzdx = (vz[k][j][i] - vz[k][j][i-1])/dx;

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
		// get effective viscosity
		eta = jr->svYZEdge[iter++].svDev.eta;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvydz = (vy[k][j][i] - vy[k-1][j][i])/dz;
		dvzdy = (vz[k][j][i] - vz[k][j-1][i])/dy;

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

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyJacobian"
PetscErrorCode JacApplyJacobian(Mat A, Vec x, Vec y)
{
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, (void**)&jr); CHKERRQ(ierr);

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x); CHKERRQ(ierr);

	ierr = JacResGetJ2Derivatives(jr); CHKERRQ(ierr);

	ierr = JacResJacobianMatFree(jr); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = JacResCopyRes(jr, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetJ2Derivatives"
PetscErrorCode JacResGetJ2Derivatives(JacRes *jr)
{
	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz;
	PetscScalar ***vx,  ***vy,  ***vz;
	PetscScalar ***centerSum, ***xyEdgeSum, ***xzEdgeSum, ***yzEdgeSum;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

    // clear local residual vectors
	ierr = VecZeroEntries(jr->ldxx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldyz); CHKERRQ(ierr);

	// access local (ghosted) velocity components
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	// access vectors that store sum of the derivative-vector products
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &centerSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &xyEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &xzEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &yzEdgeSum); CHKERRQ(ierr);

	//-------------------------------
	// central points (dxx, dyy, dzz)
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get effective deviatoric strain rates
		dxx = svCell->dxx + svCell->hxx*svDev->I2Gdt;
		dyy = svCell->dyy + svCell->hyy*svDev->I2Gdt;
		dzz = svCell->dzz + svCell->hzz*svDev->I2Gdt;

		// compute & store sum of the derivative-vector products
		centerSum[k][j][i] = dxx*(vx[k][j][i+1] - vx[k][j][i])/dx
		+                    dyy*(vy[k][j+1][i] - vy[k][j][i])/dy
		+                    dzz*(vz[k+1][j][i] - vz[k][j][i])/dz;
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points (dxy)
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// get effective deviatoric strain rate
		dxy = svEdge->d + svEdge->h*svDev->I2Gdt;

		// compute & store sum of the derivative-vector products
		xyEdgeSum[k][j][i] = dxy*(vx[k][j][i] - vx[k][j-1][i])/dy
		+                    dxy*(vy[k][j][i] - vy[k][j][i-1])/dx;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (dxz)
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get effective deviatoric strain rate
		dxz = svEdge->d + svEdge->h*svDev->I2Gdt;

		// compute & store sum of the derivative-vector products
		xzEdgeSum[k][j][i] = dxz*(vx[k][j][i] - vx[k-1][j][i])/dz
		+                    dxz*(vz[k][j][i] - vz[k][j][i-1])/dx;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (dyz)
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// get effective deviatoric strain rate
		dyz = svEdge->d + svEdge->h*svDev->I2Gdt;

		// compute & store sum of the derivative-vector products
		yzEdgeSum[k][j][i] = dyz*(vy[k][j][i] - vy[k-1][j][i])/dz
		+                    dyz*(vz[k][j][i] - vz[k][j-1][i])/dy;
	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &centerSum);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &xyEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &xzEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &yzEdgeSum); CHKERRQ(ierr);

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_XY,  jr->ldxy);
	LOCAL_TO_LOCAL(fs->DA_XZ,  jr->ldxz);
	LOCAL_TO_LOCAL(fs->DA_YZ,  jr->ldyz);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResJacobianMatFree"
PetscErrorCode JacResJacobianMatFree(JacRes *jr)
{
	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	PetscInt    I1, I2, J1, J2, K1, K2, mx, my, mz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar dx, dy, dz, tx, ty, tz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar sxx, syy, szz, sxy, sxz, syz;
	PetscScalar dxx, dyy, dzz, dxy, dxz, dyz;
	PetscScalar dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar eta, I2Gdt, dEta, DII, fr, rho, IKdt;
	PetscScalar theta, tr, pc, dt, fssa, *grav, dJ2v;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc, ***p;
	PetscScalar ***centerSum, ***xyEdgeSum, ***xzEdgeSum, ***yzEdgeSum;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs   = jr->fs;
	dt   = jr->ts->dt;     // time step
	fssa = jr->ctrl.FSSA; // density gradient penalty parameter
    grav = jr->ctrl.grav; // gravity acceleration

    // initialize maximum node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

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

	// access vectors that store sum of the derivative-vector products
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &centerSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &xyEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &xzEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &yzEdgeSum); CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// get parameters
		eta   = svDev->eta;
		I2Gdt = svDev->I2Gdt;
		dEta  = svDev->dEta;
		DII   = svDev->DII;
		fr    = svDev->fr;
		rho   = svBulk->rho;
		IKdt  = svBulk->IKdt;

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
		tr  = theta/3.0;
		dxx -= tr;
		dyy -= tr;
		dzz -= tr;

		// access current pressure
		pc = p[k][j][i];

		// compute total Cauchy stresses
		sxx = 2.0*eta*dxx - pc;
		syy = 2.0*eta*dyy - pc;
		szz = 2.0*eta*dzz - pc;

		// accumulate derivative-vector products (normalized)
		dJ2v = (centerSum[k][j][i]
		// x-y plane, i-j indices
		+      (xyEdgeSum[k][j][i]
		+       xyEdgeSum[k][j+1][i]
		+       xyEdgeSum[k][j][i+1]
		+       xyEdgeSum[k][j+1][i+1]
		// x-z plane, i-k indices
		+       xzEdgeSum[k][j][i]
		+       xzEdgeSum[k+1][j][i]
		+       xzEdgeSum[k][j][i+1]
		+       xzEdgeSum[k+1][j][i+1]
		// y-z plane, j-k indices
		+       yzEdgeSum[k][j][i]
		+       yzEdgeSum[k+1][j][i]
		+       yzEdgeSum[k][j+1][i]
		+       yzEdgeSum[k+1][j+1][i])/4.0)/DII;

		// get effective deviatoric strain rates (normalized)
		dxx = (svCell->dxx + svCell->hxx*I2Gdt)/DII;
		dyy = (svCell->dyy + svCell->hyy*I2Gdt)/DII;
		dzz = (svCell->dzz + svCell->hzz*I2Gdt)/DII;

		// get total sums
		sxx += (dEta*dJ2v + fr*pc)*dxx;
		syy += (dEta*dJ2v + fr*pc)*dyy;
		szz += (dEta*dJ2v + fr*pc)*dzz;

		// compute stabilization terms (lumped approximation)
		tx = -fssa*dt*(rho*grav[0]);
		ty = -fssa*dt*(rho*grav[1]);
		tz = -fssa*dt*(rho*grav[2]);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
		fx[k][j][i] -= (sxx + vx[k][j][i]*tx)/bdx;   fx[k][j][i+1] += (sxx + vx[k][j][i+1]*tx)/fdx;
		fy[k][j][i] -= (syy + vy[k][j][i]*ty)/bdy;   fy[k][j+1][i] += (syy + vy[k][j+1][i]*ty)/fdy;
		fz[k][j][i] -= (szz + vz[k][j][i]*tz)/bdz;   fz[k+1][j][i] += (szz + vz[k+1][j][i]*tz)/fdz;

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
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svDev  = &svEdge->svDev;

		// get parameters
		eta   = svDev->eta;
		I2Gdt = svDev->I2Gdt;
		dEta  = svDev->dEta;
		DII   = svDev->DII;
		fr    = svDev->fr;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// compute velocity gradients
		dvxdy = (vx[k][j][i] - vx[k][j-1][i])/dy;
		dvydx = (vy[k][j][i] - vy[k][j][i-1])/dx;

		// compute stress
		sxy = eta*(dvxdy + dvydx);

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// accumulate derivative-vector products (normalized)
		dJ2v = (xyEdgeSum[k][j][i]
		// x-y plane, i-j indices (i & j - bounded)
		+      (centerSum[k][J1][I1]
		+       centerSum[k][J1][I2]
		+       centerSum[k][J2][I1]
		+       centerSum[k][J2][I2]
		// y-z plane j-k indices (j - bounded)
		+       xzEdgeSum[k][J1][i]
		+       xzEdgeSum[k+1][J1][i]
		+       xzEdgeSum[k][J2][i]
		+       xzEdgeSum[k+1][J2][i]
		// x-z plane i-k indices (i - bounded)
		+       yzEdgeSum[k][j][I1]
		+       yzEdgeSum[k+1][j][I1]
		+       yzEdgeSum[k][j][I2]
		+       yzEdgeSum[k+1][j][I2])/4.0)/DII;

		// get effective deviatoric strain rate (normalized)
		dxy = (svEdge->d + svEdge->h*I2Gdt)/DII;

		// access current pressure (x-y plane, i-j indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k][j-1][i] + p[k][j-1][i-1]);

		// get total sums
		sxy += (dEta*dJ2v + fr*pc)*dxy;

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
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get parameters
		eta   = svDev->eta;
		I2Gdt = svDev->I2Gdt;
		dEta  = svDev->dEta;
		DII   = svDev->DII;
		fr    = svDev->fr;

		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvxdz = (vx[k][j][i] - vx[k-1][j][i])/dz;
		dvzdx = (vz[k][j][i] - vz[k][j][i-1])/dx;

		// compute stress
		sxz = eta*(dvxdz + dvzdx);

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// accumulate derivative-vector product (normalized)
		dJ2v = (xzEdgeSum[k][j][i]
		// x-z plane, i-k indices (i & k - bounded)
		+      (centerSum[K1][j][I1]
		+       centerSum[K1][j][I2]
		+       centerSum[K2][j][I1]
		+       centerSum[K2][j][I2]
		// y-z plane, j-k indices (k - bounded)
		+       xyEdgeSum[K1][j][i]
		+       xyEdgeSum[K1][j+1][i]
		+       xyEdgeSum[K2][j][i]
		+       xyEdgeSum[K2][j+1][i]
		// x-y plane, i-j indices (i - bounded)
		+       yzEdgeSum[k][j][I1]
		+       yzEdgeSum[k][j+1][I1]
		+       yzEdgeSum[k][j][I2]
		+       yzEdgeSum[k][j+1][I2])/4.0)/DII;

		// get effective deviatoric strain rate (normalized)
		dxz = (svEdge->d + svEdge->h*I2Gdt)/DII;

		// access current pressure (x-z plane, i-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k-1][j][i] + p[k-1][j][i-1]);

		// get total sums
		sxz += (dEta*dJ2v + fr*pc)*dxz;

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
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svDev  = &svEdge->svDev;

		// get parameters
		eta   = svDev->eta;
		I2Gdt = svDev->I2Gdt;
		dEta  = svDev->dEta;
		DII   = svDev->DII;
		fr    = svDev->fr;

		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dvydz = (vy[k][j][i] - vy[k-1][j][i])/dz;
		dvzdy = (vz[k][j][i] - vz[k][j-1][i])/dy;

		// compute stress
		syz = eta*(dvydz + dvzdy);

		// check index bounds
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// accumulate derivative-vector products (normalized)
		dJ2v = (yzEdgeSum[k][j][i]
		// y-z plane, j-k indices (j & k - bounded)
		+      (centerSum[K1][J1][i]
		+       centerSum[K1][J2][i]
		+       centerSum[K2][J1][i]
		+       centerSum[K2][J2][i]
		// x-z plane, i-k indices (k -bounded)
		+       xyEdgeSum[K1][j][i]
		+       xyEdgeSum[K1][j][i+1]
		+       xyEdgeSum[K2][j][i]
		+       xyEdgeSum[K2][j][i+1]
		// x-y plane, i-j indices (j - bounded)
		+       xzEdgeSum[k][J1][i]
		+       xzEdgeSum[k][J1][i+1]
		+       xzEdgeSum[k][J2][i]
		+       xzEdgeSum[k][J2][i+1])/4.0)/DII;

		// get effective deviatoric strain rate (normalized)
		dyz = (svEdge->d + svEdge->h*I2Gdt)/DII;

		// access current pressure (y-z plane, j-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j-1][i] + p[k-1][j][i] + p[k-1][j-1][i]);

		// get total sums
		syz += (dEta*dJ2v + fr*pc)*dyz;

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// update momentum residuals
		fy[k-1][j][i] -= syz/bdz;   fy[k][j][i] += syz/fdz;
		fz[k][j-1][i] -= syz/bdy;   fz[k][j][i] += syz/fdy;
	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,   &gc);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);         CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);        CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &centerSum); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &xyEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &xzEdgeSum); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &yzEdgeSum); CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
/*
		cellSum[k][j][i] =
		vx[k][j][i+1]*(2.0*dxx/3.0 -     dyy/3.0 -     dzz/3.0)/dx -
		vx[k][j][i]  *(2.0*dxx/3.0 -     dyy/3.0 -     dzz/3.0)/dx +
		vy[k][j+1][i]*(   -dxx/3.0 + 2.0*dyy/3.0 -     dzz/3.0)/dy -
		vy[k][j][i]  *(   -dxx/3.0 + 2.0*dyy/3.0 -     dzz/3.0)/dy +
		vz[k+1][j][i]*(   -dxx/3.0 -     dyy/3.0 + 2.0*dzz/3.0)/dz -
		vz[k][j][i]  *(   -dxx/3.0 -     dyy/3.0 + 2.0*dzz/3.0)/dz;
*/

