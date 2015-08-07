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
 **    filename:   JacResTemp.c
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
//......................   TEMPERATURE FUNCTIONS   ..........................
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

#define SCATTER_FIELD(da, vec, FIELD) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec)

#define GET_KC \
	GetTempParam(numPhases, phases, jr->svCell[iter++].phRat, &kc, NULL, NULL); \
	buff[k][j][i] = kc;

#define GET_HRXY buff[k][j][i] = jr->svXYEdge[iter++].svDev.Hr;
#define GET_HRXZ buff[k][j][i] = jr->svXZEdge[iter++].svDev.Hr;
#define GET_HRYZ buff[k][j][i] = jr->svYZEdge[iter++].svDev.Hr;

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCreateTempParam"
PetscErrorCode JacResCreateTempParam(JacRes *jr)
{
	// setup temperature parameters

	FDSTAG *fs;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create temperature DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &jr->DA_T); CHKERRQ(ierr);

	// create temperature preconditioner matrix
	ierr = DMCreateMatrix(jr->DA_T, &jr->Att); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(jr->Att, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(jr->Att, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(jr->Att, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(jr->Att, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// temperature vectors
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gT); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->dT); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->lT); CHKERRQ(ierr);

	// energy residual
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->ge); CHKERRQ(ierr);

	// create temperature diffusion solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &jr->tksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(jr->tksp,"ts_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(jr->tksp);            CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResDestroyTempParam"
PetscErrorCode JacResDestroyTempParam(JacRes *jr)
{
	// destroy temperature parameters

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// temperature parameters
	ierr = DMDestroy (&jr->DA_T); CHKERRQ(ierr);
	ierr = MatDestroy(&jr->Att);  CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gT);   CHKERRQ(ierr);
	ierr = VecDestroy(&jr->dT);   CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lT);   CHKERRQ(ierr);

	ierr = VecDestroy(&jr->ge);   CHKERRQ(ierr);

	ierr = KSPDestroy(&jr->tksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitTemp"
PetscErrorCode JacResInitTemp(JacRes *jr)
{
	// initialize temperature from markers

	PetscScalar *T;
	FDSTAG      *fs;
	PetscInt     jj, n;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// copy temperatures from context storage to global vector
	ierr = VecGetArray(jr->gT, &T);  CHKERRQ(ierr);

	for(jj = 0, n = fs->nCells; jj < n; jj++) T[jj] = jr->svCell[jj].svBulk.Tn;

	ierr = VecRestoreArray(jr->gT, &T);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyTemp"
PetscErrorCode JacResCopyTemp(JacRes *jr)
{
	// copy temperature from global to local vectors, enforce boundary constraints

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscInt    mcx, mcy, mcz;
	PetscInt    I, J, K, fi, fj, fk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***lT, ***bcT;
	PetscScalar *T, pmdof;
	PetscScalar *vals;
	PetscInt    num, *list;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vector
	ierr = VecGetArray(jr->gT, &T); CHKERRQ(ierr);

	//=================================
	// enforce single point constraints
	//=================================

	// temperature
	num   = bc->tNumSPC;
	list  = bc->tSPCList;
	vals  = bc->tSPCVals;

	for(i = 0; i < num; i++) T[list[i]] = vals[i];

	ierr = VecRestoreArray(jr->gT, &T); CHKERRQ(ierr);

	// fill local (ghosted) version of solution vector
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gT, jr->lT)

	// access local solution vector
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	// access boundary constraints vector
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcT, &bcT); CHKERRQ(ierr);

	//==============================
	// enforce two-point constraints
	//==============================

	//-----------------------------
	// central points (temperature)
	//-----------------------------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lT[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcT, lT, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcT, lT, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcT, lT, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcT, lT, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcT, lT, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcT, lT, K, j, i, pmdof) }

		if(fi*fj)    SET_EDGE_CORNER(n, lT, k, J, I, k, j, i, pmdof)
		if(fi*fk)    SET_EDGE_CORNER(n, lT, K, j, I, k, j, i, pmdof)
		if(fj*fk)    SET_EDGE_CORNER(n, lT, K, J, i, k, j, i, pmdof)
		if(fi*fj*fk) SET_EDGE_CORNER(n, lT, K, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetTempRes"
PetscErrorCode JacResGetTempRes(JacRes *jr)
{
	// compute temperature residual vector

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	PetscInt    iter, numPhases;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
 	PetscScalar dx, dy, dz;
	PetscScalar kc, rho, Cp, A, Tc, Tn, dt, Hr;
	PetscScalar ***ge, ***T, ***lk, ***hxy, ***hxz, ***hyz, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs        = jr->fs;
	numPhases = jr->numPhases; // number phases
	phases    = jr->phases;    // phase parameters
	dt        = jr->ts.dt;     // time step

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_KC)
	SCATTER_FIELD(fs->DA_XY,  jr->ldxy, GET_HRXY)
	SCATTER_FIELD(fs->DA_XZ,  jr->ldxz, GET_HRXZ)
	SCATTER_FIELD(fs->DA_YZ,  jr->ldyz, GET_HRYZ)

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ge,   &ge);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// access
		Tc  = T[k][j][i];  // current temperature
		Tn  = svBulk->Tn;  // temperature history
		rho = svBulk->rho; // effective density

		// conductivity, heat capacity, radiogenic heat production
		GetTempParam(numPhases, phases, svCell->phRat, &kc, &Cp, &A);

		// shear heating term
		Hr = svDev->Hr +
		(hxy[k][j][i] + hxy[k][j+1][i] + hxy[k][j][i+1] + hxy[k][j+1][i+1] +
		 hxz[k][j][i] + hxz[k+1][j][i] + hxz[k][j][i+1] + hxz[k+1][j][i+1] +
		 hyz[k][j][i] + hyz[k+1][j][i] + hyz[k][j+1][i] + hyz[k+1][j+1][i])/4.0;

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average conductivities
		bkx = (kc + lk[k][j][Im1])/2.0;      fkx = (kc + lk[k][j][Ip1])/2.0;
		bky = (kc + lk[k][Jm1][i])/2.0;      fky = (kc + lk[k][Jp1][i])/2.0;
		bkz = (kc + lk[Km1][j][i])/2.0;      fkz = (kc + lk[Kp1][j][i])/2.0;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute heat fluxes
		bqx = bkx*(Tc - T[k][j][i-1])/bdx;   fqx = fkx*(T[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - T[k][j-1][i])/bdy;   fqy = fky*(T[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - T[k-1][j][i])/bdz;   fqz = fkz*(T[k+1][j][i] - Tc)/fdz;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// original balance equation:

		// rho*Cp*(Tc - Tn)/dt = (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz + Hr + rho*A

		// to get positive diagonal in the preconditioner matrix
		// put right hand side to the left, which gives the following:

		ge[k][j][i] = rho*Cp*(Tc - Tn)/dt - (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz - Hr - rho*A;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ge,   &ge);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetTempMat"
PetscErrorCode JacResGetTempMat(JacRes *jr)
{
	// assemble temperature preconditioner matrix

	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	Material_t *phases;
	PetscInt    iter, numPhases;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;
	PetscScalar v[7], cf[6], kc, rho, Cp, dt;
	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcT, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs        = jr->fs;
	bc        = jr->bc;
	numPhases = jr->numPhases; // number phases
	phases    = jr->phases;    // phase parameters
	dt        = jr->ts.dt;     // time step

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_KC)

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcT,  &bcT); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svBulk = &svCell->svBulk;

		// access
		rho = svBulk->rho; // effective density

		// conductivity, heat capacity
		GetTempParam(numPhases, phases, svCell->phRat, &kc, &Cp, NULL);

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcT[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcT[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcT[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcT[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcT[k-1][j][i] != DBL_MAX) cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcT[k+1][j][i] != DBL_MAX) cf[5] = -1.0; }

		// compute average conductivities
		bkx = (kc + lk[k][j][Im1])/2.0;      fkx = (kc + lk[k][j][Ip1])/2.0;
		bky = (kc + lk[k][Jm1][i])/2.0;      fky = (kc + lk[k][Jp1][i])/2.0;
		bkz = (kc + lk[Km1][j][i])/2.0;      fkz = (kc + lk[Kp1][j][i])/2.0;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC multipliers
		v[0] = -bkx/bdx/dx*cf[0];
		v[1] = -fkx/fdx/dx*cf[1];
		v[2] = -bky/bdy/dy*cf[2];
		v[3] = -fky/fdy/dy*cf[3];
		v[4] = -bkz/bdz/dz*cf[4];
		v[5] = -fkz/fdz/dz*cf[5];
		v[6] =  rho*Cp/dt
		+       (bkx/bdx + fkx/fdx)/dx
		+       (bky/bdy + fky/fdy)/dy
		+       (bkz/bdz + fkz/fdz)/dz;

		// set matrix coefficients
		ierr = MatSetValuesStencil(jr->Att, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

		// NOTE! since only TPC are active, no SPC modification is necessary
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT);  CHKERRQ(ierr);

	// assemble temperature matrix
	ierr = MatAssemblyBegin(jr->Att, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (jr->Att, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
Diffusion term expansion

		bqx = bkx*(Tc - T[k][j][i-1])/bdx;   fqx = fkx*(T[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - T[k][j-1][i])/bdy;   fqy = fky*(T[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - T[k-1][j][i])/bdz;   fqz = fkz*(T[k+1][j][i] - Tc)/fdz;

[1]
		-(fqx - bqx)/dx
		-(fqy - bqy)/dy
		-(fqz - bqz)/dz
[2]
		(bqx - fqx)/dx
		(bqy - fqy)/dy
		(bqz - fqz)/dz
[3]
		(bkx*(Tc - T[k][j][i-1])/bdx - fkx*(T[k][j][i+1] - Tc)/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy - fky*(T[k][j+1][i] - Tc)/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz - fkz*(T[k+1][j][i] - Tc)/fdz)/dz
[4]
		(bkx*(Tc - T[k][j][i-1])/bdx + fkx*(Tc - T[k][j][i+1])/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy + fky*(Tc - T[k][j+1][i])/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz + fkz*(Tc - T[k+1][j][i])/fdz)/dz
[5]
		bkx/bdx/dx*(Tc - T[k][j][i-1]) + fkx/fdx/dx*(Tc - T[k][j][i+1])
		bky/bdy/dy*(Tc - T[k][j-1][i]) + fky/fdy/dy*(Tc - T[k][j+1][i])
		bkz/bdz/dz*(Tc - T[k-1][j][i]) + fkz/fdz/dz*(Tc - T[k+1][j][i])
[6]
		(bkx/bdx/dx + fkx/fdx/dx)*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy/dy + fky/fdy/dy)*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz/dz + fkz/fdz/dz)*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
[7]

		(bkx/bdx + fkx/fdx)/dx*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy + fky/fdy)/dy*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz + fkz/fdz)/dz*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
*/
//---------------------------------------------------------------------------
