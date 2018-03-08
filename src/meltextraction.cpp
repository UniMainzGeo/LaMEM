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
 **    filename:   adjoint.c
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

/*
 * TODO:
 * - Pay attention local and global might still be inconsistent in Exchange Volume
 * - So far only tested for exactly 1 extraction step
 * - What happens if the melt quantity actually decreases instead of increase? This might be not treated here yet
 * - What pressure and T should the new marker have?
 * - Two things are hard-coded (ctrl+F "hard-coded" to find them)
 * - How to handle the actual volume change in injection?:
 * 		1) Delete all markers from cell and inject equally distributed (and volume equal) markers of injection phase?
 * 		2) Reduce volume of existing markers in the cell such that volume is again correct?
 * 		3) Is it phyiscally actually correct if it is just corrected for the inflated volumeas is done right now?
 * 		4) Artificially massively increase the injected volume until it is correct?
 *
 * 	POTENTIAL HELP:
 * 	- Implement a timestep criterium for the melt extraction. Only incrementally move melt up such that velocity doesn't run away.
 */

#include "LaMEM.h"
#include "cvi.h"
#include "phase.h"
#include "fdstag.h"
#include "constEq.h"
#include "advect.h"
#include "tools.h"
#include "JacRes.h"
#include "meltextraction.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionCreate"
PetscErrorCode MeltExtractionCreate(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMV);      CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMVmerge); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (jr->fs->DA_CEN, &jr->ldMV);      CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdc);       CHKERRQ(ierr);
	ierr = DMCreateLocalVector (jr->fs->DA_CEN, &jr->ldc);       CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionDestroy"
PetscErrorCode MeltExtractionDestroy(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecDestroy(&jr->gdMV);           CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdMVmerge);      CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldMV);           CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdc);            CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldc);            CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionSave"
PetscErrorCode MeltExtractionSave(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscScalar  ***ldMV;
	PetscInt      i, j, k, nx, ny, nz, sx, sy, sz, iter;
	SolVarBulk   *svBulk;
	SolVarCell   *svCell;
	FDSTAG       *fs;

	fs = jr->fs;

	ierr = DMDAVecGetArray (jr->fs->DA_CEN, jr->ldMV, &ldMV);      CHKERRQ(ierr);

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx);
	GET_CELL_RANGE(ny, sy, fs->dsy);
	GET_CELL_RANGE(nz, sz, fs->dsz);

	START_STD_LOOP
	{
		svCell = &jr->svCell[iter++];
		svBulk = &svCell->svBulk;

		ldMV[k][j][i] = svBulk->dMF;

	}END_STD_LOOP

	ierr = DMDAVecRestoreArray (jr->fs->DA_CEN, jr->ldMV, &ldMV);      CHKERRQ(ierr);

	LOCAL_TO_GLOBAL(fs->DA_CEN, jr->ldMV, jr->gdMV);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarker"
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx)
{
	//=======================================================================
	// interpolate increments of the history field of melt extraction to markers
	//=======================================================================

		FDSTAG      *fs;
		JacRes      *jr;
		Marker      *P;
		AdvVelCtx   vi;
		PetscScalar  UP;
		PetscInt     nx, ny, sx, sy, sz;
		PetscInt     jj, ID, I, J, K;
		PetscScalar *vgdc, ***vldc;

		PetscErrorCode ierr;
		PetscFunctionBegin;

		ierr = ADVelCreate(actx, &vi);  CHKERRQ(ierr);

		fs = actx->fs;
		jr = actx->jr;

		// starting indices & number of cells
		sx = fs->dsx.pstart; nx = fs->dsx.ncels;
		sy = fs->dsy.pstart; ny = fs->dsy.ncels;
		sz = fs->dsz.pstart;


		// access 1D layouts of global vectors
		ierr = VecGetArray(jr->gdc, &vgdc);  CHKERRQ(ierr);

		for(jj = 0; jj < fs->nCells; jj++) vgdc[jj] = jr->svCell[jj].svBulk.dMF;

		// restore access
		ierr = VecRestoreArray(jr->gdc, &vgdc);  CHKERRQ(ierr);

		// communicate boundary values
		GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gdc, jr->ldc);


		// access 3D layouts of local vectors
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

		// scan ALL markers
		for(jj = 0; jj < actx->nummark; jj++)
		{
			// access next marker
			P = &actx->markers[jj];

			// get consecutive index of the host cell
			ID = actx->cellnum[jj];

			// expand I, J, K cell indices
			GET_CELL_IJK(ID, I, J, K, nx, ny)

			// access buffer
			UP = vldc[sz+K][sy+J][sx+I];

			if(UP > 0)
			{
				ierr = MeltExtractionInject(actx,&vi, ID, I, J, K, UP);  CHKERRQ(ierr);
				vldc[sz+K][sy+J][sx+I] = 0; // It is all distributed in markers
			}
			else
			{
				P->Mtot += -UP;  // has to increase by the amount of melt change
				P->Mvol +=  UP;  // has to decrease
			}
		}

		// restore access
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

		ierr = ADVelDestroy(&vi);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarkerBackToGrid"
PetscErrorCode MeltExtractionInterpMarkerBackToGrid(AdvCtx *actx)
{
	// Project the following history fields from markers to grid:

	// - Mtot

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     ID, I, J, K, II, JJ, KK;
	PetscInt     ii, jj, numPhases;
	PetscInt     nx, ny, sx, sy, sz, nCells;
	PetscScalar  xp, yp, zp, xc, yc, zc, wxn, wyn, wzn, wxc, wyc, wzc, w = 0.0;
	PetscScalar  UPXY, UPXZ, UPYZ;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs        = actx->fs;
	jr        = actx->jr;
	numPhases = actx->dbm->numPhases;

	// check marker phases
	ierr = ADVCheckMarkPhases(actx); CHKERRQ(ierr);

	//======
	// CELLS
	//======

	// marker-to-grid projection (cell nodes)

	// number of cells
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nCells = fs->nCells;

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
		svCell->svBulk.mfextot = 0.0;
		svCell->svBulk.mfVol  = 0.0;
	}

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables
		svCell->svBulk.mfextot += w*P->Mtot;
		svCell->svBulk.mfVol   += w*P->Mvol;

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		// normalize history variables
		svCell->svBulk.mfextot /= w;
		svCell->svBulk.mfVol   /= w;
	}

	//======
	// EDGES
	//======

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// clear local vectors
	ierr = VecZeroEntries(jr->ldxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldyz); CHKERRQ(ierr);

	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// set interpolated fields to defaults
	UPXY = 1.0; UPXZ = 1.0; UPYZ = 1.0;

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get coordinates of cell center
		xc = fs->dsx.ccoor[I];
		yc = fs->dsy.ccoor[J];
		zc = fs->dsz.ccoor[K];

		// map marker on the control volumes of edge nodes
		if(xp > xc) { II = I+1; } else { II = I; }
		if(yp > yc) { JJ = J+1; } else { JJ = J; }
		if(zp > zc) { KK = K+1; } else { KK = K; }

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get interpolation weights in node control volumes
		wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);
		wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);
		wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);

		UPXY = P->Mtot;  UPXZ = P->Mtot;  UPYZ = P->Mtot;

		// update required fields from marker to edge nodes
		lxy[sz+K ][sy+JJ][sx+II] += wxn*wyn*wzc*UPXY;
		lxz[sz+KK][sy+J ][sx+II] += wxn*wyc*wzn*UPXZ;
		lyz[sz+KK][sy+JJ][sx+I ] += wxc*wyn*wzn*UPYZ;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// assemble global vectors
	LOCAL_TO_GLOBAL(fs->DA_XY, jr->ldxy, jr->gdxy)
	LOCAL_TO_GLOBAL(fs->DA_XZ, jr->ldxz, jr->gdxz)
	LOCAL_TO_GLOBAL(fs->DA_YZ, jr->ldyz, jr->gdyz)

	// access 1D layouts of global vectors
	ierr = VecGetArray(jr->gdxy, &gxy);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdxz, &gxz);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdyz, &gyz);  CHKERRQ(ierr);

	// copy (normalized) data to the residual context
	for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].svDev.mfextot = gxy[jj]/jr->svXYEdge[jj].ws;
	for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].svDev.mfextot = gxz[jj]/jr->svXZEdge[jj].ws;
	for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].svDev.mfextot = gyz[jj]/jr->svYZEdge[jj].ws;

	// restore access
	ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionExchangeVolume"
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr)
{
	FDSTAG      *fs;
	Discret1D   *dsz;
	SolVarCell  *svCell;
	Controls    *ctrl;
	PetscInt     i, j, k, K, sx, sy, sz, nx, ny, nz, iter;
	Vec          dgMVVec, dgMVVecmerge, dlMVVec, dlMVVecmerge;
	PetscScalar  bz, ez;
	PetscScalar  level;
	PetscScalar *vdlMVVecmerge, *vdgMVVec, *vdgMVVecmerge, **vdgMVVecmerge2, **vdgMVVec2, **vdlMVVecmerge2, **vdlMVVec2, ***vdlMV, ***vdgMV;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs    = jr->fs;
	dsz   = &fs->dsz;
	ctrl  = &jr->ctrl;

	// get local coordinate bounds
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);

	// create column communicator
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &dgMVVec);      CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &dgMVVecmerge); CHKERRQ(ierr);
	ierr = DMGetLocalVector (jr->DA_CELL_2D, &dlMVVec);      CHKERRQ(ierr);
	ierr = DMGetLocalVector (jr->DA_CELL_2D, &dlMVVecmerge); CHKERRQ(ierr);
	ierr = VecZeroEntries   (dgMVVec);                       CHKERRQ(ierr);
	ierr = VecZeroEntries   (dgMVVecmerge);                  CHKERRQ(ierr);
	ierr = VecZeroEntries   (dlMVVec);                       CHKERRQ(ierr);
	ierr = VecZeroEntries   (dlMVVecmerge);                  CHKERRQ(ierr);

	// scan all local cells
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dlMVVec, &vdlMVVec2); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldMV ,    &vdlMV);CHKERRQ(ierr);
	START_STD_LOOP
	{
		vdlMVVec2[j][i] += vdlMV[k][j][i];
	}
	END_STD_LOOP
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldMV  ,   &vdlMV);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dlMVVec, &vdlMVVec2); CHKERRQ(ierr);

	LOCAL_TO_GLOBAL(jr->DA_CELL_2D,dlMVVec,dgMVVec);

	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(dgMVVec, &vdgMVVec); CHKERRQ(ierr);
		ierr = VecGetArray(dgMVVecmerge, &vdgMVVecmerge); CHKERRQ(ierr);

		ierr = MPI_Allreduce(vdgMVVec, vdgMVVecmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);

		ierr = VecRestoreArray(dgMVVec, &vdgMVVec); CHKERRQ(ierr);
		ierr = VecRestoreArray(dgMVVecmerge, &vdgMVVecmerge); CHKERRQ(ierr);

		GLOBAL_TO_LOCAL(jr->DA_CELL_2D, dgMVVecmerge, dlMVVecmerge);
	}
	else
	{
		ierr = VecCopy(dgMVVec,dgMVVecmerge);  CHKERRQ(ierr);
		GLOBAL_TO_LOCAL(jr->DA_CELL_2D, dgMVVecmerge, dlMVVecmerge);
	}

	level = ctrl->DExt;

	// scan all local cells
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dlMVVecmerge, &vdlMVVecmerge2); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldMV ,    &vdlMV);CHKERRQ(ierr);
	START_PLANE_LOOP
	{
		if(vdlMVVecmerge2[j][i] < 0)
		{
			// check whether point belongs to domain
			if(level >= bz && level < ez)
			{
				// find containing cell
				K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, level);

				// interpolate velocity
				vdlMV[sz+K][j][i] = -vdlMVVecmerge2[j][i];
			}
		}
	}
	END_PLANE_LOOP
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dlMVVecmerge, &vdlMVVecmerge2); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldMV ,    &vdlMV);CHKERRQ(ierr);
	LOCAL_TO_GLOBAL(jr->DA_CELL_2D, dlMVVecmerge, dgMVVecmerge);


	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldMV ,    &vdlMV);CHKERRQ(ierr);
	iter = 0;
	// Finally bring it back to the center nodes
	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		// scan all phases
		//for(i = 0; i < numPhases; i++)
		//{

			// update present phases only
			//if(phRat[i])
			//{
				// get reference to material parameters table
				//mat = &phases[i];

				// Get PD data
				//if(mat->Pd_rho == 1)
				//{
					// Only interpolate positive anomalies (negative was already extracted)
					if(vdlMV[k][j][i] > 0)
					{
						svCell->svBulk.dMF  = vdlMV[k][j][i];
					}
				//}

			//}

		//}

	}END_STD_LOOP
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldMV ,    &vdlMV);CHKERRQ(ierr);
	LOCAL_TO_GLOBAL(fs->DA_CEN, jr->ldMV, jr->gdMV);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInject"
PetscErrorCode MeltExtractionInject(AdvCtx *actx, AdvVelCtx *vi, PetscInt ID, PetscInt I, PetscInt J, PetscInt K, PetscScalar UP)
{

	PetscInt    jj, ipn, n, ninj, pind, found, PhInject;
	PetscScalar xs[3], xe[3], xp[3];
	PetscScalar cf_rand;
	PetscRandom    rctx;
	Marker      *P;
	FDSTAG      *fs;
	Controls    *ctrl;
	JacRes      *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	found = 0;
	jr    = actx->jr;
	fs    = actx->fs;
	ctrl  = &jr->ctrl;

	PhInject = ctrl->PhExt;

	// get markers in cell
	n = vi->markstart[ID+1] - vi->markstart[ID];

	// scan cell markers
	for(jj = 0; jj < n; jj++)
	{
		// get marker index
		pind = vi->markind[vi->markstart[ID] + jj];

		P = &actx->markers[pind];

		if(P->phase == PhInject && UP > 0)
		{
			UP = UP - (1-P->Mvol);
			if (UP < 0)
			{
				P->Mvol += UP + (1-P->Mvol);   // has to increase
				UP = 0;
			}
			else
			{
				P->Mvol += (1-P->Mvol);   // has to increase
			}
			found = 1;
		}
	}

	// We have not found a marker of the correct phase or there is still melt to be injected
	if(found == 0 || UP > 0)
	{
		ninj = (PetscInt)ceil(UP);  // Amount of markers we have to inject

		// allocate memory for new markers
		actx->nrecv = ninj;
		ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);
		ierr = PetscMemzero(actx->recvbuf, (size_t)actx->nrecv*sizeof(Marker)); CHKERRQ(ierr);

		// initialize the random number generator
		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

		// get cell coordinates
		xs[0] = fs->dsx.ncoor[I]; xe[0] = fs->dsx.ncoor[I+1];
		xs[1] = fs->dsy.ncoor[J]; xe[1] = fs->dsy.ncoor[J+1];
		xs[2] = fs->dsz.ncoor[K]; xe[2] = fs->dsz.ncoor[K+1];

		for(ipn = 0; ipn<ninj; ipn++)
		{
			// create random coordinate within this cell
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[0] = (xe[0] - xs[0]) * cf_rand + xs[0];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[1] = (xe[1] - xs[1]) * cf_rand + xs[1];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[2] = (xe[2] - xs[2]) * cf_rand + xs[2];

			/*// calculate the closest (parent marker)
			for (ii = 0; ii < n; ii++)
			{
				X  = markers[ii].X;
				dx[0] = X[0] - xp[0];
				dx[1] = X[1] - xp[1];
				dx[2] = X[2] - xp[2];
				x  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

				if (ii == 0)    { sumind = x; sind = ii; }
				if (x < sumind) { sumind = x; sind = ii; }
			}*/

			// create new marker
			// actx->recvbuf[ipn]      = markers[sind];

			// hard-coded new marker properties for debugging
			actx->recvbuf[ipn].phase = PhInject;
			actx->recvbuf[ipn].p = 10000;
			actx->recvbuf[ipn].T = 10000;
			actx->recvbuf[ipn].APS = 5;
			actx->recvbuf[ipn].Mtot = 0;

			actx->recvbuf[ipn].S.xx = 0;
			actx->recvbuf[ipn].S.xy = 0;
			actx->recvbuf[ipn].S.xz = 0;
			actx->recvbuf[ipn].S.yy = 0;
			actx->recvbuf[ipn].S.yz = 0;
			actx->recvbuf[ipn].S.zz = 0;

			if (UP > 1)
			{
				actx->recvbuf[ipn].Mvol = 2;
				UP -= 1;
			}
			else
			{
				actx->recvbuf[ipn].Mvol = 1+UP;
				UP = 0;
			}

			actx->recvbuf[ipn].U[0] = 0;
			actx->recvbuf[ipn].U[1] = 0;
			actx->recvbuf[ipn].U[2] = 0;
			actx->recvbuf[ipn].X[0] = xp[0];
			actx->recvbuf[ipn].X[1] = xp[1];
			actx->recvbuf[ipn].X[2] = xp[2];
		}

		// destroy random context
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

		// store new markers
		actx->ndel = 0;
		ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

		// compute host cells for all the markers
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

		// print info
		// PetscPrintf(PETSC_COMM_WORLD,"Melt extraction injected %i markers.\n", ninj);

		// clear
		ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}















