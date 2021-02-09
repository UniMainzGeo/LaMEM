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
 **    filename:   JacResFlow.c
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
//......................   FLUID FLOW FUNCTIONS   ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matrix.h"
#include "surf.h"
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

#define GET_KI \
	ierr = JacResGetFlowParam(jr, jr->svCell[iter++].phRat, &ki, NULL); CHKERRQ(ierr); \
	buff[k][j][i] = ki;

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetFlowParam"
PetscErrorCode JacResGetFlowParam(
		JacRes      *jr,
		PetscScalar *phRat,
		PetscScalar *ki_,   // permeability
		PetscScalar *Ss_)   // specific storage
{
	// compute effective fluid flow parameters in the cell

	PetscInt    i, numPhases;
    Material_t  *phases, *M;

	PetscScalar cf, ki, Ss;

	PetscFunctionBegin;

	// initialize
	ki        = 0.0;
	Ss        = 0.0;
	numPhases = jr->dbm->numPhases;
	phases    = jr->dbm->phases;

	// average all phases
	for(i = 0; i < numPhases; i++)
	{
		M   = &phases[i];
		cf  =  phRat[i];
		ki +=  cf*M->ki;
		Ss +=  cf*M->Ss;
	}

	// store
	if(ki_) (*ki_) = ki;
    if(Ss_) (*Ss_) = Ss;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCheckFlowParam"
PetscErrorCode JacResCheckFlowParam(JacRes *jr)
{
	// check whether fluid flow parameters are properly defined

    Material_t  *phases, *M;
	PetscInt    i, numPhases;

	PetscFunctionBegin;

	// fluid flow cases only
	if(!jr->ctrl.actFluid) PetscFunctionReturn(0);

	// initialize
	numPhases = jr->dbm->numPhases;
	phases    = jr->dbm->phases;

	// check all phases
	for(i = 0; i < numPhases; i++)
	{
		M = &phases[i];

		if(M->ki == 0.0) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define permeability of phase %lld\n", (LLD)i);
		if(M->Ss == 0.0) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define specific storage of phase %lld\n", (LLD)i);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCreateFlowParam"
PetscErrorCode JacResCreateFlowParam(JacRes *jr)
{
	// setup fluid flow parameters

	FDSTAG *fs;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// fluid flow cases only
	if(!jr->ctrl.actFluid) PetscFunctionReturn(0);

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create fluid pressure DMDA
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &jr->DA_P); CHKERRQ(ierr);

	// set proper interpolation type for multigrid
	ierr = DMDASetInterpolationType(jr->DA_P, DMDA_Q0); CHKERRQ(ierr);

	// create fluid pressure preconditioner matrix
	ierr = DMCreateMatrix(jr->DA_P, &jr->App); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// fluid pressure increment vector
	ierr = DMCreateGlobalVector(jr->DA_P, &jr->dP); CHKERRQ(ierr);

	// fluid flow residual
	ierr = DMCreateGlobalVector(jr->DA_P, &jr->gf); CHKERRQ(ierr);

	// pressure solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &jr->pksp); CHKERRQ(ierr);

	// enable geometric multigrid
	ierr = KSPSetDM(jr->pksp, jr->DA_P);           CHKERRQ(ierr);
	ierr = KSPSetDMActive(jr->pksp, PETSC_FALSE);  CHKERRQ(ierr);

	// set options
	ierr = KSPSetOptionsPrefix(jr->pksp,"ps_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(jr->pksp);            CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResDestroyFlowParam"
PetscErrorCode JacResDestroyFlowParam(JacRes *jr)
{
	// destroy fluid flow parameters

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// fluid flow cases only
	if(!jr->ctrl.actFluid) PetscFunctionReturn(0);

	// temperature parameters
	ierr = DMDestroy (&jr->DA_P); CHKERRQ(ierr);
	ierr = MatDestroy(&jr->App);  CHKERRQ(ierr);
	ierr = VecDestroy(&jr->dP);   CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gf);   CHKERRQ(ierr);
	ierr = KSPDestroy(&jr->pksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResSaveFlow"
PetscErrorCode JacResSaveFlow(JacRes *jr)
{
	// save fluid pressure to history database

	FDSTAG      *fs;
	PetscScalar ***lP;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// fluid flow cases only
	if(!jr->ctrl.actFluid) PetscFunctionReturn(0);

	// access context
	fs = jr->fs;

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		jr->svCell[iter++].svBulk.fn = lP[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResUpdateFlow"
PetscErrorCode JacResUpdateFlow(JacRes *jr)
{
	// correct fluid pressure (Newton update)

	FDSTAG      *fs;
	PetscScalar ***lP, ***dP;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lP); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_P,   jr->dP,      &dP); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		lP[k][j][i] -= dP[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lP); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_P,   jr->dP,      &dP); CHKERRQ(ierr);

	// apply boundary constraints
	ierr = JacResApplyFlowBC(jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResApplyFlowBC"
PetscErrorCode JacResApplyFlowBC(JacRes *jr)
{
	// apply fluid pressure constraints

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscScalar pmdof;
	PetscScalar ***lP, ***bcf;
	PetscInt    mcx, mcy, mcz;
	PetscInt    I, J, K, fi, fj, fk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// set fluid pressure SPC
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcf,     &bcf); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		if(bcf[k][j][i] != DBL_MAX) lP[k][j][i] = bcf[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore,  &lP);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcf,      &bcf); CHKERRQ(ierr);

	// exchange internal ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_pore)

	// set fluid pressure TPC
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore,  &lP);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcf,      &bcf); CHKERRQ(ierr);

	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lP[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcf, lP, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcf, lP, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcf, lP, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcf, lP, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcf, lP, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcf, lP, K, j, i, pmdof) }

		if(fi && fj)       SET_EDGE_CORNER(n, lP, k, J, I, k, j, i, pmdof)
		if(fi && fk)       SET_EDGE_CORNER(n, lP, K, j, I, k, j, i, pmdof)
		if(fj && fk)       SET_EDGE_CORNER(n, lP, K, J, i, k, j, i, pmdof)
		if(fi && fj && fk) SET_EDGE_CORNER(n, lP, K, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore,  &lP);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcf,      &bcf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitFlow"
PetscErrorCode JacResInitFlow(JacRes *jr)
{
	// initialize pore pressure with hydrostatic pressure
	PetscLogDouble t;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// fluid flow cases only
	if(!jr->ctrl.actFluid) PetscFunctionReturn(0);

	PrintStart(&t,"Computing steady-state fluid pressure distribution", NULL);

	// initialize
	ierr = VecZeroEntries(jr->lp_pore); CHKERRQ(ierr);

	// set boundary constraints
	ierr = BCApplyFlowBC(jr->bc); CHKERRQ(ierr);
	ierr = JacResApplyFlowBC(jr); CHKERRQ(ierr);

	// assemble residual and jacobian (steady-state)
	ierr = JacResGetFlowRes(jr, 0.0); CHKERRQ(ierr);
	ierr = JacResGetFlowMat(jr, 0.0); CHKERRQ(ierr);

	// solve for steady-state pressure distriution
	ierr = KSPSetOperators(jr->pksp, jr->App, jr->App); CHKERRQ(ierr);
	ierr = KSPSetUp(jr->pksp);                          CHKERRQ(ierr);
	ierr = KSPSolve(jr->pksp, jr->gf, jr->dP);          CHKERRQ(ierr);

	// store initial guess
	ierr = JacResUpdateFlow(jr); CHKERRQ(ierr);

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetFlowRes"
PetscErrorCode JacResGetFlowRes(JacRes *jr, PetscScalar dt)
{
	// compute fluid flow residual vector
	// STEADY STATE solution is activated by setting time step to zero

	FDSTAG     *fs;
	BCCtx      *bc;
	Controls   *ctrl;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	Vec         vki;
	PetscInt    iter, cellID, I, J, K, M, N;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, jj;
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
 	PetscScalar dx, dy, dz, hx, hy, hz;
	PetscScalar invdt, ki, Ss, pc, fn, rho, eta, gz;
	PetscScalar ***gf, ***lP, ***lk, ***buff, ***bcf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs       = jr->fs;
	bc       = jr->bc;
	mx       = fs->dsx.tcels - 1;
	my       = fs->dsy.tcels - 1;
	mz       = fs->dsz.tcels - 1;
	M        = fs->dsx.ncels;
	N        = fs->dsy.ncels;
	ctrl     = &jr->ctrl;
	rho      = ctrl->rho_fluid;
	eta      = ctrl->eta_fluid;
	gz       = PetscAbsScalar(ctrl->grav[2]);

	// compute inverse time step
	if(dt) invdt = 1.0/dt;
	else   invdt = 0.0;

	ierr = DMGetLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	SCATTER_FIELD(fs->DA_CEN, vki, GET_KI)

	// access work vectors
	ierr = DMDAVecGetArray(jr->DA_P,   jr->gf,       &gf);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore,  &lP);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, vki,          &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcf,      &bcf); CHKERRQ(ierr);

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

		// constrain residual
		if(bcf[k][j][i] != DBL_MAX)
		{
			gf[k][j][i] = 0.0;
			continue;
		}

		// access current & history pressure
		pc = lP[k][j][i];
		fn = svBulk->fn;

		// permeability, specific storage
		ierr = JacResGetFlowParam(jr, svCell->phRat, &ki, &Ss); CHKERRQ(ierr);

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average permeabilities normalized by viscosity
		bkx = (ki + lk[k][j][Im1])/2.0/eta;      fkx = (ki + lk[k][j][Ip1])/2.0/eta;
		bky = (ki + lk[k][Jm1][i])/2.0/eta;      fky = (ki + lk[k][Jp1][i])/2.0/eta;
		bkz = (ki + lk[Km1][j][i])/2.0/eta;      fkz = (ki + lk[Kp1][j][i])/2.0/eta;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute fluid fluxes
		bqx = bkx* (pc - lP[k][j][i-1])/bdx;               fqx = fkx*(lP[k][j][i+1]  - pc)/fdx;
		bqy = bky* (pc - lP[k][j-1][i])/bdy;               fqy = fky*(lP[k][j+1][i]  - pc)/fdy;
		bqz = bkz*((pc - lP[k-1][j][i])/bdz + rho*gz);     fqz = fkz*((lP[k+1][j][i] - pc)/fdz + rho*gz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// compute residual
		gf[k][j][i] = Ss*(invdt*(pc - fn)) - (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz;
	}
	END_STD_LOOP


	// apply point fluid sources
	if(!ctrl->initGuess)
	{
		for(jj = 0; jj < bc->nsource; jj++)
		{
			cellID = bc->isource[jj];

			if(cellID != -1)
			{
				// get I, J, K cell indices
				GET_CELL_IJK(cellID, I, J, K, M, N);

				// get cell sizes
				hx = SIZE_CELL(I, 0, fs->dsx);
				hy = SIZE_CELL(J, 0, fs->dsy);
				hz = SIZE_CELL(K, 0, fs->dsz);

				// add source scaled by cell volume
				gf[sz+K][sy+J][sx+I] -= bc->vsource[jj]/hx*hy*hz;
			}
		}
	}

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_P,   jr->gf,      &gf);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, vki,         &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcf,     &bcf); CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetFlowMat"
PetscErrorCode JacResGetFlowMat(JacRes *jr, PetscScalar dt)
{
	// assemble fluid pressure preconditioner matrix
	// STEADY STATE solution is activated by setting time step to zero

	FDSTAG     *fs;
	BCCtx      *bc;
	Controls   *ctrl;
	SolVarCell *svCell;
	Vec         vki;
	PetscInt    iter;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;
	PetscScalar v[7], cf[6], ki, Ss, invdt, eta;
	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcf, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs   = jr->fs;
	bc   = jr->bc;
	mx   = fs->dsx.tcels - 1;
	my   = fs->dsy.tcels - 1;
	mz   = fs->dsz.tcels - 1;
	ctrl = &jr->ctrl;
	eta  = ctrl->eta_fluid;

	// compute inverse time step
	if(dt) invdt = 1.0/dt;
	else   invdt = 0.0;

	ierr = DMGetLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, vki, GET_KI)

	// clear matrix coefficients
	ierr = MatZeroEntries(jr->App); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, vki,     &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcf, &bcf); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		// constrain matrix row
		if(bcf[k][j][i] != DBL_MAX)
		{
			row[0].k = k;
			row[0].j = j;
			row[0].i = i;
			row[0].c = 0;
			v  [0]   = 1.0;

			ierr = MatSetValuesStencil(jr->App, 1, row, 1, row, v, ADD_VALUES); CHKERRQ(ierr);

			continue;
		}

		// permeability, specific storage
		ierr = JacResGetFlowParam(jr, svCell->phRat, &ki, &Ss); CHKERRQ(ierr);

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcf[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcf[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcf[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcf[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcf[k-1][j][i] != DBL_MAX) cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcf[k+1][j][i] != DBL_MAX) cf[5] = -1.0; }

		// compute average permeabilities normalized by viscosity
		bkx = (ki + lk[k][j][Im1])/2.0/eta;      fkx = (ki + lk[k][j][Ip1])/2.0/eta;
		bky = (ki + lk[k][Jm1][i])/2.0/eta;      fky = (ki + lk[k][Jp1][i])/2.0/eta;
		bkz = (ki + lk[Km1][j][i])/2.0/eta;      fkz = (ki + lk[Kp1][j][i])/2.0/eta;

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
		v[6] =  invdt*Ss
		+       (bkx/bdx + fkx/fdx)/dx
		+       (bky/bdy + fky/fdy)/dy
		+       (bkz/bdz + fkz/fdz)/dz;

		// set matrix coefficients
		ierr = MatSetValuesStencil(jr->App, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, vki,     &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcf, &bcf); CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	// assemble fluid pressure matrix
	ierr = MatAIJAssemble(jr->App, 0, NULL, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetFlowSource"
PetscErrorCode JacResGetFlowSource(JacRes *jr)
{
	// compute fluid flow Stokes source

	FDSTAG     *fs;
	Controls   *ctrl;
	SolVarCell *svCell;
	Vec         vki;
	PetscInt    iter, fluidPhase;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
 	PetscScalar dx, dy, dz, V, lV, tV, fRat, lflux, flux, source, lbuf[2], gbuf[2];
	PetscScalar ki, Ss, pc, rho, eta, gz;
	PetscScalar ***lP, ***lk, ***buff, ***vol;
	Vec         lvol;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// relevant cases only
	if(jr->ctrl.fluidPhase == -1 || jr->ctrl.initGuess) PetscFunctionReturn(0);

	// access context
	fs         = jr->fs;
	mx         = fs->dsx.tcels - 1;
	my         = fs->dsy.tcels - 1;
	mz         = fs->dsz.tcels - 1;
	ctrl       = &jr->ctrl;
	rho        = ctrl->rho_fluid;
	eta        = ctrl->eta_fluid;
	fluidPhase = ctrl->fluidPhase;
	gz         =  PetscAbsScalar(ctrl->grav[2]);

	ierr = DMGetLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	//===============
	// compute volume
	//===============

	ierr = DMGetLocalVector(fs->DA_CEN, &lvol);     CHKERRQ(ierr);
	ierr = VecZeroEntries(lvol);                    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, lvol, &vol); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;
	lV   = 0.0;

	START_STD_LOOP
	{
		fRat = jr->svCell[iter++].phRat[fluidPhase];

		if(fRat > 0.0)
		{
			// compute cell volume, update total volume
			dx  = SIZE_CELL(i, sx, fs->dsx);
			dy  = SIZE_CELL(j, sy, fs->dsy);
			dz  = SIZE_CELL(k, sz, fs->dsz);
			V   = dx*dy*dz*fRat;
			lV += V;

			// store volume
			vol[k][j][i] = V;

		}
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvol, &vol); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, lvol)

	//=============
	// compute flux
	//=============

	SCATTER_FIELD(fs->DA_CEN, vki, GET_KI)

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, vki,         &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, lvol,        &vol); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter  = 0;
	lflux = 0.0;

	START_STD_LOOP
	{
		if(vol[k][j][i] == 0.0) continue;

		// access solution variables
		svCell = &jr->svCell[iter++];

		// access current pressure
		pc  = lP[k][j][i];

		// permeability, specific storage
		ierr = JacResGetFlowParam(jr, svCell->phRat, &ki, &Ss); CHKERRQ(ierr);

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average permeabilities normalized by viscosity
		bkx = (ki + lk[k][j][Im1])/2.0/eta;      fkx = (ki + lk[k][j][Ip1])/2.0/eta;
		bky = (ki + lk[k][Jm1][i])/2.0/eta;      fky = (ki + lk[k][Jp1][i])/2.0/eta;
		bkz = (ki + lk[Km1][j][i])/2.0/eta;      fkz = (ki + lk[Kp1][j][i])/2.0/eta;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute fluid fluxes
		bqx = bkx* (pc - lP[k][j][i-1])/bdx;            fqx = fkx* (lP[k][j][i+1] - pc)/fdx;
		bqy = bky* (pc - lP[k][j-1][i])/bdy;            fqy = fky* (lP[k][j+1][i] - pc)/fdy;
		bqz = bkz*((pc - lP[k-1][j][i])/bdz + rho*gz);  fqz = fkz*((lP[k+1][j][i] - pc)/fdz + rho*gz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// compute integral
		if(i > 0)  { if(vol[k]  [j]  [i-1] == 0.0) { lflux -= (bqx*dy*dz); } }
		if(i < mx) { if(vol[k]  [j]  [i+1] == 0.0) { lflux += (fqx*dy*dz); } }
		if(j > 0)  { if(vol[k]  [j-1][i]   == 0.0) { lflux -= (bqy*dx*dz); } }
		if(j < my) { if(vol[k]  [j+1][i]   == 0.0) { lflux += (fqy*dx*dz); } }
		if(k > 0)  { if(vol[k-1][j]  [i]   == 0.0) { lflux -= (bqz*dx*dy); } }
		if(k < mz) { if(vol[k+1][j]  [i]   == 0.0) { lflux += (fqz*dx*dy); } }

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, vki,         &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvol,        &vol); CHKERRQ(ierr);

	//===============
	// compute source
	//===============

	// compute global sum
	if(ISParallel(PETSC_COMM_WORLD))
	{
		lbuf[0] = lV;
		lbuf[1] = lflux;

		ierr = MPI_Allreduce(lbuf, gbuf, 2, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

		tV   = gbuf[0];
		flux = gbuf[1];

	}
	else
	{
		tV   = lV;
		flux = lflux;
	}

	if(tV == 0.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Fluid Stokes domain cannot have zero volume (fluid_phase)\n");
	}

	// compute source term
	source = flux/tV;

	// store source term
	ierr = DMDAVecGetArray(fs->DA_CEN, lvol, &vol); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;

	// store source term
	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		if(vol[k][j][i] != 0.0)
		{
			svCell->source = source;
		}
		else
		{
			svCell->source = 0.0;
		}
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvol, &vol); CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(fs->DA_CEN, &lvol); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(fs->DA_CEN, &vki);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetFlowFlux"
PetscErrorCode JacResGetFlowFlux(JacRes *jr, Vec lvx, Vec lvy, Vec lvz)
{
	// compute fluid flow Stokes source

	FDSTAG     *fs;
	Controls   *ctrl;
	SolVarCell *svCell;
	Vec         vki;
	PetscInt    iter;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
	PetscScalar ki, Ss, pc, rho, eta, gz;
	PetscScalar ***lP, ***lk, ***buff, ***vx, ***vy, ***vz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs       = jr->fs;
	mx       = fs->dsx.tcels - 1;
	my       = fs->dsy.tcels - 1;
	mz       = fs->dsz.tcels - 1;
	ctrl     = &jr->ctrl;
	rho      = ctrl->rho_fluid;
	eta      = ctrl->eta_fluid;
	gz       = PetscAbsScalar(ctrl->grav[2]);

	ierr = DMGetLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

	//===============
	// compute fluxes
	//===============

	SCATTER_FIELD(fs->DA_CEN, vki, GET_KI)

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, vki     ,    &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, lvx,         &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, lvy,         &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, lvz,         &vz);  CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter  = 0;

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		// access current pressure
		pc = lP[k][j][i];

		if(pc < 0.0)
		{
			vx[k][j][i] = 0.0;
			vy[k][j][i] = 0.0;
			vz[k][j][i] = 0.0;

			continue;
		}

		// permeability, specific storage
		ierr = JacResGetFlowParam(jr, svCell->phRat, &ki, &Ss); CHKERRQ(ierr);

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average permeabilities normalized by viscosity

		bkx = (ki + lk[k][j][Im1])/2.0/eta;      fkx = (ki + lk[k][j][Ip1])/2.0/eta;
		bky = (ki + lk[k][Jm1][i])/2.0/eta;      fky = (ki + lk[k][Jp1][i])/2.0/eta;
		bkz = (ki + lk[Km1][j][i])/2.0/eta;      fkz = (ki + lk[Kp1][j][i])/2.0/eta;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute fluid fluxes
		bqx = bkx* (pc - lP[k][j][i-1])/bdx;            fqx = fkx* (lP[k][j][i+1] - pc)/fdx;
		bqy = bky* (pc - lP[k][j-1][i])/bdy;            fqy = fky* (lP[k][j+1][i] - pc)/fdy;
		bqz = bkz*((pc - lP[k-1][j][i])/bdz + rho*gz);  fqz = fkz*((lP[k+1][j][i] - pc)/fdz + rho*gz);

		// store velocities
		vx[k][j][i] = -(bqx + fqx)/2.0;
		vy[k][j][i] = -(bqy + fqy)/2.0;
		vz[k][j][i] = -(bqz + fqz)/2.0;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lP);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, vki,         &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvx,         &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvy,         &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvz,         &vz);  CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(fs->DA_CEN, &vki); CHKERRQ(ierr);

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
