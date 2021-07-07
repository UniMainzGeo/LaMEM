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
 **    filename:   JacResAux.c
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

#include "LaMEM.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "fdstag.h"
#include "surf.h"
#include "phase.h"
#include "tools.h"
#include "Tensor.h"
#include "parsing.h"

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
#define __FUNCT__ "JacResGetOverPressure"
PetscErrorCode JacResGetOverPressure(JacRes *jr, Vec lop)
{
	// compute overpressure

	FDSTAG      *fs;
	PetscScalar ***op, ***p, ***p_lith;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs  =  jr->fs;

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// access pressure vectors
	ierr = VecZeroEntries(lop); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, lop,         &op);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// overpressure = dynamic - lithostatic
		op[k][j][i] = p[k][j][i] - p_lith[k][j][i];
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lop,         &op);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);

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
	PetscScalar ***lp, ***ibuff, *lbuff, dz, dp, g, rho,mean_p;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, iter, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs  =  jr->fs;
	dsz = &fs->dsz;
	L   =  (PetscInt)dsz->rank;
	g   =   PetscAbsScalar(jr->ctrl.grav[2]);

	// initialize
	ierr = VecZeroEntries(jr->lp_lith); CHKERRQ(ierr);

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
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &lp); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &lp); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff, &ibuff); CHKERRQ(ierr);

	ierr = VecRestoreArray(vbuff, &lbuff); CHKERRQ(ierr);

	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff); CHKERRQ(ierr);

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_lith)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetPorePressure"
PetscErrorCode JacResGetPorePressure(JacRes *jr)
{
	// compute pore pressure
	FDSTAG      *fs;
	Controls    *ctrl;
	Material_t  *phases, *mat;
	PetscScalar ***lp_pore, ***lp_lith, *phRat;
	PetscScalar ztop, g, gwLevel=0.0, rho_fluid, depth, p_hydro, rp_cv, rp;
	PetscInt    numPhases, i, j, k, iter, iphase, sx, sy, sz, nx, ny, nz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	ierr = VecZeroEntries(jr->lp_pore); CHKERRQ(ierr);

	// return if not activated
	if(jr->ctrl.gwType == _GW_NONE_) PetscFunctionReturn(0);

	// access context
	fs        =  jr->fs;
	phases    =  jr->dbm->phases;
	numPhases =  jr->dbm->numPhases;
	ctrl      = &jr->ctrl;
	rho_fluid =  ctrl->rho_fluid;
	g         =  PetscAbsScalar(ctrl->grav[2]);

	// get top boundary coordinate
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, NULL, NULL, NULL, &ztop); CHKERRQ(ierr);

	// set ground water level
	if     (ctrl->gwType == _GW_TOP_)   gwLevel = ztop;
	else if(ctrl->gwType == _GW_SURF_)  gwLevel = jr->surf->avg_topo;
	else if(ctrl->gwType == _GW_LEVEL_) gwLevel = ctrl->gwLevel;

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lp_pore); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &lp_lith); CHKERRQ(ierr);

	iter = 0;
	START_STD_LOOP
	{
		// access phase ratio array
		phRat = jr->svCell[iter++].phRat;

		// compute depth of the current control volume
		depth = gwLevel - COORD_CELL(k, sz, fs->dsz);
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
		p_hydro = rho_fluid * g * PetscAbsScalar(depth);

		// compute the pore pressure as product of lithostatic pressure and porepressure ratio of the control volume
		lp_pore[k][j][i] =  p_hydro + rp_cv * (lp_lith[k][j][i]-p_hydro);
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lp_pore); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &lp_lith); CHKERRQ(ierr);

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_pore)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetPermea"
PetscErrorCode JacResGetPermea(JacRes *jr, PetscInt bgPhase, PetscInt step, char *outfile)
{
	FILE        *db;
	FDSTAG      *fs;
	BCCtx       *bc;
	Material_t  *phases;
	Scaling     *scal;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***vz, nZFace, lvel, gvel, dp, eta, ks, bz, ez;
	char        path[_str_len_];


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	if(!jr->ctrl.getPermea || !step) PetscFunctionReturn(0);

	// access context variables
	fs     = jr->fs;
	bc     = jr->bc;
	phases = jr->dbm->phases;
	scal   = jr->scal;

	// get grid bounds in z-direction
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);

	// get total number of z-faces
	nZFace = (PetscScalar)(fs->dsx.tcels*fs->dsy.tcels*fs->dsz.tnods);

	// get fluid viscosity (background phase is fluid)
	eta = 1.0/(2.0*phases[bgPhase].Bd);

	// get pressure gradient
	dp = (bc->pbot - bc->ptop)/(ez-bz);

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	// access z-velocity vector
	ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &vz);  CHKERRQ(ierr);

	// compute volume average absolute velocity
	lvel = 0.0;

	//---------
	// Z-points
	//---------

	START_STD_LOOP
	{
		lvel += PetscAbsScalar(vz[k][j][i]);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &vz);  CHKERRQ(ierr);

	// compute global sum
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lvel, &gvel, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gvel = lvel;
	}

	// normalize
	gvel /= nZFace;

	// compute permeability
	ks = gvel*eta/dp;

	ks = PetscAbsScalar(ks);

	// output to the file
	if(ISRankZero(PETSC_COMM_WORLD))
	{

		memset(path, 0, _str_len_);
		strcpy(path, outfile);
		strcat(path, ".darcy.dat");

		db = fopen(path, "w");

		fprintf(db,"# ==============================================\n");
		fprintf(db,"# EFFECTIVE PERMEABILITY CONSTANT: %E %s \n ", ks*scal->area_si, scal->lbl_area_si);
		fprintf(db,"# ==============================================\n");

		fclose(db);

		// Also output this to the screen (for regression testing, for example)
		PetscPrintf(PETSC_COMM_WORLD,"\n");
		PetscPrintf(PETSC_COMM_WORLD,"==========================================================================\n");
		PetscPrintf(PETSC_COMM_WORLD,"EFFECTIVE PERMEABILITY CONSTANT: %E %s\n", ks*scal->area_si, scal->lbl_area_si);
		PetscPrintf(PETSC_COMM_WORLD,"==========================================================================\n");
		
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
