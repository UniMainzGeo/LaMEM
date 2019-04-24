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
 **    filename:   subgrid.c
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
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "subgrid.h"
#include "advect.h"
#include "phase.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "bc.h"
#include "JacRes.h"
#include "surf.h"
#include "marker.h"
#include "AVD.h"
#include "cvi.h"
#include "Tensor.h"
#include "tools.h"

/*
#START_DOC#
\lamemfunction{\verb- ADVCreate -}
Create advection context

\lamemfunction{\verb- ADVDestroy -}
Destroy advection context

\lamemfunction{\verb- ADVAdvect -}
Main advection routine

#END_DOC#
*/

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSubGrid"
PetscErrorCode ADVMarkSubGrid(AdvCtx *actx)
{
	// check marker distribution and merge or clone markers based on subgrid

	FDSTAG           *fs;
	PetscInt          icell, isubcell, imark, I, J, K, i, j, ib, ie;
	PetscInt          ncx, ncy, npx, npy, npz, nmark, ncell, nclone, nmerge;
	PetscInt         *markind;
	PetscScalar       s[3], h[3], *x;
	PetscLogDouble    t0, t1;
	ipair             t;
	spair             d;
	vector <Marker>   iclone;
	vector <PetscInt> imerge;
	vector <ipair>    cell;
	vector <spair>    dist;
	vector <Marker>   mark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscTime(&t0); CHKERRQ(ierr);

	// access context
	fs    = actx->fs;
	ncx   = fs->dsx.ncels;
	ncy   = fs->dsy.ncels;
	npx   = actx->NumPartX;
	npy   = actx->NumPartY;
	npz   = actx->NumPartZ;
	ncell = npx*npy*npz;

	// reserve space for index & distance storage
	cell.reserve(_mark_buff_sz_);
	dist.reserve(_mark_buff_sz_);
	mark.reserve(_mark_buff_sz_);

	iclone.reserve(actx->nummark*_mark_buff_ratio_/100);
	imerge.reserve(actx->nummark*_mark_buff_ratio_/100);

	iclone.clear();
	imerge.clear();

	nclone = 0;
	nmerge = 0;

	// process local cells
	for(icell = 0; icell < fs->nCells; icell++)
	{
		// get number of markers per cell
		nmark = actx->markstart[icell+1] - actx->markstart[icell];

		if(!nmark)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Empty control volume");
		}

		// expand I, J, K cell indices
		GET_CELL_IJK(icell, I, J, K, ncx, ncy)

		// get cell starting coordinates
		s[0] = COORD_NODE(I, 0, fs->dsx);
		s[1] = COORD_NODE(J, 0, fs->dsy);
		s[2] = COORD_NODE(K, 0, fs->dsz);

		// get subcell sizes
		h[0] = SIZE_CELL(I, 0, fs->dsx)/(PetscScalar)npx;
		h[1] = SIZE_CELL(J, 0, fs->dsy)/(PetscScalar)npy;
		h[2] = SIZE_CELL(K, 0, fs->dsz)/(PetscScalar)npz;

		// clear index storage
		cell.clear();

		// map markers on subcells
		markind = actx->markind + actx->markstart[icell];

		for(j = 0; j < nmark; j++)
		{
			// get marker index & coordinates
			imark = markind[j];
			x      = actx->markers[imark].X;

			// compute containing subcell index
			MAP_SUBCELL(I, x[0], s[0], h[0], npx);
			MAP_SUBCELL(J, x[1], s[1], h[1], npy);
			MAP_SUBCELL(K, x[2], s[2], h[2], npz);

			GET_CELL_ID(isubcell, I, J, K, npx, npy);

			// store marker and subcell indices
			t.first  = isubcell;
			t.second = imark;
			cell.push_back(t);
		}

		// sort markers by subcells
		sort(cell.begin(), cell.end());

		// push end-of-array stamp
		t.first  = -1;
		t.second =  0;
		cell.push_back(t);

		// process subcells
		for(i = 0, ib = 0; i < ncell; i++)
		{
			// get index of next populated subcell
			isubcell = cell[ib].first;

			// check whether current subcell is empty
			if(isubcell != i)
			{
				// clone markers
				ierr = ADVMarkClone(actx, icell, i, s, h, dist, iclone); CHKERRQ(ierr);

				// update counter
				nclone++;
			}
			else
			{
				// find next populated subcell
				ie = ib; while(cell[ie].first == isubcell) ie++;

				// merge markers if required
				if(ie - ib > actx->npmax)
				{
					ierr = ADVMarkCheckMerge(actx, ib, ie, nmerge, mark, cell, iclone, imerge); CHKERRQ(ierr);
				}

				// switch to next populated subcell
				ib = ie;
			}
		}
	}

	// rearrange storage after marker resampling
	ierr = ADVCollectGarbageVec(actx, iclone, imerge); CHKERRQ(ierr);

	// compute host cells for all the markers
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,
		"Marker control [%lld]: (subgrid) cloned %lld markers and merged %lld markers in %1.4e s\n",
		(LLD)actx->iproc, (LLD)nclone, (LLD)nmerge, t1-t0);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkClone"
PetscErrorCode ADVMarkClone(
	AdvCtx          *actx,
	PetscInt         icell,
	PetscInt         isubcell,
	PetscScalar      s[3],
	PetscScalar      h[3],
	vector <spair>  &dist,
	vector <Marker> &iclone)
{
	// clone closest marker & put it in the center of an empty subcell
	// current marker storage is not modified, the following is done instead:
	//  - all newly created markers are stored for insertion in iclone

	BCCtx            *bc;
	spair             d;
	Marker            P;
	PetscScalar       xc[3], *x;
	PetscInt          I, J, K, j, npx, npy, imark, nmark, *markind;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bc      = actx->jr->bc;
	npx     = actx->NumPartX;
	npy     = actx->NumPartY;
	nmark   = actx->markstart[icell+1] - actx->markstart[icell];
	markind = actx->markind + actx->markstart[icell];

	// expand I, J, K subcell indices
	GET_CELL_IJK(isubcell, I, J, K, npx, npy)

	// get coordinates of subcell center
	COORD_SUBCELL(xc[0], I, s[0], h[0]);
	COORD_SUBCELL(xc[1], J, s[1], h[1]);
	COORD_SUBCELL(xc[2], K, s[2], h[2]);

	// clear distance storage
	dist.clear();

	// find closest marker (cell-wise approximation)
	for(j = 0; j < nmark; j++)
	{
		// get marker index & coordinates
		imark = markind[j];
		x     = actx->markers[imark].X;

		// store marker and distance
		d.first  = EDIST(x, xc);
		d.second = imark;
		dist.push_back(d);
	}

	// sort markers by distance
	sort(dist.begin(), dist.end());

	// clone closest marker
	P = actx->markers[dist.begin()->second];

	// place clone in cell center
	P.X[0] = xc[0];
	P.X[1] = xc[1];
	P.X[2] = xc[2];

	// override marker phase (if necessary)
	ierr = BCOverridePhase(bc, icell, &P); CHKERRQ(ierr);

	// store cloned marker
	iclone.push_back(P);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkCheckMerge"
PetscErrorCode ADVMarkCheckMerge(
	AdvCtx            *actx,
	PetscInt           ib,
	PetscInt           ie,
	PetscInt          &nmerge,
	vector <Marker>   &mark,
	vector <ipair>    &cell,
	vector <Marker>   &iclone,
	vector <PetscInt> &imerge)
{
	// merge markers in a densely populated subcell
	// never merge markers of different phases
	// current marker storage is not modified, the following is done instead:
	//  - indices of merged markers are flagged for removal in imerge
	//  - all newly created markers are stored for insertion in iclone
	//  - difference between original and final number of markers is added to counter

	PetscInt j, jb, je, k, sz, phase, nmark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// copy marker phase IDs
	for(j = ib; j < ie; j++)
	{
		cell[j].first = actx->markers[cell[j].second].phase;
	}

	// sort markers by phase
	sort(cell.begin() + ib, cell.begin() + ie);

	// merge markers of the same phase
	jb = ib;

	do
	{
		// get current phase ID
		phase = cell[jb].first;

		// find all markers of current phase
		je = jb; while(cell[je].first == phase && je < ie) je++;

		// get number of markers
		nmark = je - jb;

		// merge if number of markers exceeds threshold
		if(nmark > actx->npmax)
		{
			// copy markers to buffer
			mark.clear();

			for(j = jb; j < je; j++)
			{
				mark.push_back(actx->markers[cell[j].second]);
			}

			// merge markers
			ierr = ADVMarkMerge(mark, nmark, actx->npmax, sz); CHKERRQ(ierr);

			// update counter
			nmerge += nmark - actx->npmax;

			// flag merged markers
			for(j = jb, k = 0; j < je; j++, k++)
			{
				if(mark[k].phase == -1)
				{
					imerge.push_back(cell[j].second);
				}
			}

			// store new markers
			for(k = nmark; k < sz; k++)
			{
				if(mark[k].phase != -1)
				{
					iclone.push_back(mark[k]);
				}
			}
		}

		// switch to next phase
		jb = je;

	} while(je < ie);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkMerge"
PetscErrorCode ADVMarkMerge(
	vector <Marker> &mark,
	PetscInt         nmark,
	PetscInt         npmax,
	PetscInt        &sz)
{
	// recursively find and merge closest markers until required number is reached
	// put new markers in the end of the storage, mark merged markers with phase -1

	Marker       P;
	PetscInt     j, k, jmin, kmin;
	PetscScalar  d, dmin;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize storage size
	sz = nmark;

	while(nmark > npmax)
	{
		// find closest markers
		dmin  = DBL_MAX;
		jmin  = 0;
		kmin  = 0;

		for(j = 0; j < sz; j++)
		{
			if(mark[j].phase == -1) continue;

			for(k = j+1; k < sz; k++)
			{
				if(mark[k].phase == -1) continue;

				d = EDIST(mark[j].X, mark[k].X);

				if(d < dmin)
				{
					dmin = d;
					jmin = j;
					kmin = k;
				}
			}
		}

		// merge closest markers
		ierr = MarkerMerge(mark[jmin], mark[kmin], P); CHKERRQ(ierr);

		// store new marker
		mark.push_back(P);

		// mark original markers as merged
		mark[jmin].phase = -1;
		mark[kmin].phase = -1;

		// update counters
		nmark--;
		sz++;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCollectGarbageVec"
PetscErrorCode ADVCollectGarbageVec(AdvCtx *actx, vector <Marker> &recvbuf, vector <PetscInt> &idel)
{
	// rearrange storage after marker resampling

	Marker   *markers;
	PetscInt  nummark, nrecv, ndel;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access storage
	nummark = actx->nummark;
	markers = actx->markers;
	nrecv   = (PetscInt)recvbuf.size();
	ndel    = (PetscInt)idel.size();

	// close holes in marker storage
	while(nrecv && ndel)
	{
		markers[idel[ndel-1]] = recvbuf[nrecv-1];
		nrecv--;
		ndel--;
	}

	if(nrecv)
	{
		// make sure space is enough
		ierr = ADVReAllocStorage(actx, nummark + nrecv); CHKERRQ(ierr);

		// make sure we have a correct storage pointer
		markers = actx->markers;

		// put the rest in the end of marker storage
		while(nrecv)
		{
			markers[nummark++] = recvbuf[nrecv-1];
			nrecv--;
		}
	}

	if(ndel)
	{
		// collect garbage
		while(ndel)
		{
			if(idel[ndel-1] != nummark-1)
			{
				markers[idel[ndel-1]] = markers[nummark-1];
			}
			nummark--;
			ndel--;
		}
	}

	// store new number of markers
	actx->nummark = nummark;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkCrossFreeSurf"
PetscErrorCode ADVMarkCrossFreeSurf(AdvCtx *actx)
{
	// change marker phase when crossing free surface

	Marker          *P, *IP;
	FDSTAG          *fs;
	FreeSurf        *surf;
	Vec             vphase;
	PetscInt        sx, sy, sz, nx, ny;
	PetscInt        ii, jj, ID, I, J, K, L, AirPhase, phaseID, nmark, *markind, markid;
	PetscScalar     ***ltopo, ***phase, *ncx, *ncy, topo, xp, yp, zp, *X, *IX;
    spair           d;
	vector <spair>  dist;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free-surface cases only
	if(!actx->surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	surf      = actx->surf;
	fs        = actx->fs;
	L         = fs->dsz.rank;
	AirPhase  = surf->AirPhase;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// grid coordinates
	ncx = fs->dsx.ncoor;
	ncy = fs->dsy.ncoor;

	// reserve marker distance buffer
	dist.reserve(_mark_buff_sz_);

	// request local vector for reference sedimentation phases
	ierr = DMGetLocalVector(fs->DA_CEN, &vphase);

	// compute reference sedimentation phases
	ierr = ADVGetSedPhase(actx, vphase); CHKERRQ(ierr);

	// access topography & phases
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &ltopo);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN,    vphase,      &phase);  CHKERRQ(ierr);

	// scan all markers
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

		// compute surface topography at marker position
		topo = InterpLin2D(ltopo, I, J, L, sx, sy, xp, yp, ncx, ncy);

		// check whether rock marker is above the free surface
		if(P->phase != AirPhase && zp > topo)
		{
			// erosion (physical or numerical) -> rock turns into air
			P->phase = AirPhase;

			//=======================================================================
			// WARNING! At best clone history from nearest air marker
			//=======================================================================
		}

		// check whether air marker is below the free surface
		if(P->phase == AirPhase && zp < topo)
		{
			if(surf->SedimentModel > 0)
			{
				// sedimentation (physical) -> air turns into a prescribed rock
				P->phase = surf->phase;
			}
			else
			{
				// sedimentation (numerical) -> air turns into closest (reference) rock
				X = P->X;

				// get marker list in containing cell
				nmark   = actx->markstart[ID+1] - actx->markstart[ID];
				markind = actx->markind + actx->markstart[ID];

				// clear distance storage
				dist.clear();

				for(ii = 0; ii < nmark; ii++)
				{
					// get current marker
					markid = markind[ii];
					IP     = &actx->markers[markid];

					// sort out air markers
					if(IP->phase == AirPhase) continue;

					// get marker coordinates
					IX = IP->X;

					// store marker index and distance
					d.first  = EDIST(X, IX);
					d.second = markid;

					dist.push_back(d);
				}

				// find closest rock marker (if any)
				if(dist.size())
				{
					// sort rock markers by distance
					sort(dist.begin(), dist.end());

					// copy phase from closest marker
					IP = &actx->markers[dist.begin()->second];

					P->phase = IP->phase;
				}
				else
				{
					// no local rock marker found, set phase to reference
					phaseID = (PetscInt)phase[sz+K][sy+J][sx+I];

					if(phaseID < 0)
					{
						SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect sedimentation phase");
					}

					P->phase = phaseID;
				}

				//=======================================================================
				// WARNING! At best clone history from nearest rock marker
				//=======================================================================
			}
		}
	}

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &ltopo);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN,    vphase,      &phase);  CHKERRQ(ierr);

	// restore phase vector
	ierr = DMRestoreLocalVector(fs->DA_CEN, &vphase); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVGetSedPhase"
PetscErrorCode ADVGetSedPhase(AdvCtx *actx, Vec vphase)
{
	// compute reference sedimentation phases

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscInt     nCells, nMarks, numPhases, sedPhase, AirPhase, ii, jj, ID;
	PetscScalar  maxMark, ***phase;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs        = actx->fs;
	jr        = actx->jr;
	numPhases = actx->dbm->numPhases;
	AirPhase  = jr->surf->AirPhase;

	// number of cells & markers
	nCells = fs->nCells;
	nMarks = actx->nummark;

	// clear marker counters
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		for(ii = 0; ii < numPhases; ii++)
		{
			svCell->phRat[ii] = 0.0;
		}
	}

	// update marker counters
	for(jj = 0; jj < nMarks; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		svCell->phRat[P->phase] += 1.0;
	}

	// initialize phase vector
	ierr = VecSet(vphase, -1.0); CHKERRQ(ierr);

	// compute dominant sedimentation phase
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, vphase, &phase); CHKERRQ(ierr);

	iter = 0;

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		maxMark  =  0.0;
		sedPhase = -1;

		for(ii = 0; ii < numPhases; ii++)
		{
			if(ii == AirPhase) continue;

			if(svCell->phRat[ii] > maxMark)
			{
				maxMark  = svCell->phRat[ii];
				sedPhase = ii;
			}
		}

		phase[k][j][i] = (PetscScalar)sedPhase;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, vphase, &phase); CHKERRQ(ierr);

	// exchange ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, vphase)

	// propagate sedimentation phases into one layer of air cells
	ierr = DMDAVecGetArray(fs->DA_CEN, vphase, &phase); CHKERRQ(ierr);

	START_STD_LOOP
	{
		if(phase[k][j][i] == -1.0 && phase[k-1][j][i] >= 0.0)
		{
			phase[k][j][i]   =  phase[k-1][j][i];
			phase[k+1][j][i] = -2.0;
		}
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, vphase, &phase); CHKERRQ(ierr);

	// exchange ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, vphase)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
	Marker   P;
	PetscInt nmark = 5, npmax = 2, sz;
	vector   <Marker> mark;
	mark.reserve(_mark_buff_sz_);
	mark.clear();
	PetscMemzero(&P, sizeof(Marker));
	P.phase = 1; P.X[0] = 1; P.X[1] = 1; P.X[2] = 0; mark.push_back(P);
	P.phase = 1; P.X[0] = 1; P.X[1] = 5; P.X[2] = 0; mark.push_back(P);
	P.phase = 1; P.X[0] = 3; P.X[1] = 4; P.X[2] = 0; mark.push_back(P);
	P.phase = 1; P.X[0] = 4; P.X[1] = 3; P.X[2] = 0; mark.push_back(P);
	P.phase = 1; P.X[0] = 5; P.X[1] = 5; P.X[2] = 0; mark.push_back(P);
	ierr = ADVMarkMerge(mark, nmark, npmax, sz); CHKERRQ(ierr);

	AdvCtx actx;
	Marker  P;
	PetscMemzero(&actx, sizeof(AdvCtx));
	PetscMemzero(&P,    sizeof(Marker));
	ierr = ADVReAllocStorage(&actx, 100); CHKERRQ(ierr);
	P.phase = 2; P.X[0] = 5; P.X[1] = 5; P.X[2] = 0; actx.markers[0] = P;
	P.phase = 1; P.X[0] = 1; P.X[1] = 1; P.X[2] = 0; actx.markers[1] = P;
	P.phase = 2; P.X[0] = 4; P.X[1] = 3; P.X[2] = 0; actx.markers[2] = P;
	P.phase = 1; P.X[0] = 1; P.X[1] = 5; P.X[2] = 0; actx.markers[3] = P;
	P.phase = 2; P.X[0] = 3; P.X[1] = 4; P.X[2] = 0; actx.markers[4] = P;
	actx.npmax   = 2;
	actx.nummark = 5;
	vector <Marker>   iclone;
	vector <PetscInt> imerge;
	vector <ipair>    cell;
	vector <Marker>   mark;
	cell.reserve(_mark_buff_sz_);
	mark.reserve(_mark_buff_sz_);
	iclone.reserve(actx.nummark*_mark_buff_ratio_/100);
	imerge.reserve(actx.nummark*_mark_buff_ratio_/100);
	cell.clear();
	mark.clear();
	iclone.clear();
	imerge.clear();
	PetscInt nmerge = 0, ib = 0, ie = 5;
	ipair    t;
	t.first = 0; t.second = 0; cell.push_back(t);
	t.first = 0; t.second = 1; cell.push_back(t);
	t.first = 0; t.second = 2; cell.push_back(t);
	t.first = 0; t.second = 3; cell.push_back(t);
	t.first = 0; t.second = 4; cell.push_back(t);
	ierr = ADVMarkCheckMerge(&actx, ib, ie, nmerge, mark, cell, iclone, imerge); CHKERRQ(ierr);
*/
