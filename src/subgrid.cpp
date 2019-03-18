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
#include <vector>
#include <algorithm>
#include <utility>
using namespace std;
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
	// check marker distribution and delete or inject markers if necessary

	FDSTAG         *fs;
	Marker         *P;
	PetscInt       ID, I, J, K, i, ii, jj;
	PetscInt       ncx, ncy,  npx, npy, npz, nmark, ncell;
	PetscInt       ninj, nmrg, nmax;
	PetscInt       *pid;
	PetscScalar    s[3], h[3], *x;
	PetscInt       jb, je, cellid;
	PetscLogDouble t0, t1;

    vector < pair < PetscScalar, PetscInt > > dist;
    vector < pair < PetscInt,    PetscInt > > cell;
    pair          < PetscInt,    PetscInt >   t;


//	PetscInt      randNoise;           // random noise flag for marker distribution

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
	ninj  = 0;
	nmrg  = 0;

	// get max number of markers estimate per cell
	nmax = _max_nmark_*actx->npmax;
	nmax = nmax*nmax*nmax;

	// reserve space for index storage
	cell.reserve(nmax);

	// process local cells
	for(i = 0; i < fs->nCells; i++)
	{
		// get number of markers per cell
		nmark = actx->markstart[i+1] - actx->markstart[i];

		if(!nmark)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, " Empty control volume");
		}

		// expand I, J, K cell indices
		GET_CELL_IJK(i, I, J, K, ncx, ncy)

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
		pid = actx->markind + actx->markstart[i];

		for(ii = 0; ii < nmark; ii++)
		{
			// get marker & coordinates
			jj = pid[ii];
			P  = &actx->markers[jj];
			x  = P->X;

			// compute containing subcell index
			MAP_SUBCELL(I, x[0], s[0], h[0], npx);
			MAP_SUBCELL(J, x[1], s[1], h[1], npy);
			MAP_SUBCELL(K, x[2], s[2], h[2], npz);
			GET_CELL_ID(ID, I, J, K, npx, npy);

			// store marker and cell numbers
			t.first  = ID;
			t.second = jj;
			cell.push_back(t);

		}

		// sort markers by subcells
		sort(cell.begin(), cell.end());

		// push and-of-array stamp
		t.first  = -1;
		t.second =  0;
		cell.push_back(t);

		// process cells
		for(ii = 0, jb = 0; ii < ncell; ii++)
		{
			// get index of next populated cell
			cellid = cell[jb].first;

			// check whether current cell is empty
			if(cellid > ii)
			{
				// process empty subcell (insert markers)

				// ...............

				ninj++;
			}
			else
			{
				// find next populated cell or and-of-array stamp
				je = jb; while(cell[je].first == cellid) je++;

				// process non-empty subcell (merge markers of same phase)
				if(jb-je > actx->npmax)
				{
					// ...............

					nmrg++;
				}

				// switch to next populated cell
				jb = je;
			}
		}
	}

	// store new markers
//	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// compute host cells for all the markers
//	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// update arrays for marker-cell interaction
//	ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);


	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Marker control [%lld]: (subgrid) injected %lld markers and merged %lld markers in %1.4e s\n",(LLD)actx->iproc, (LLD)ninj, (LLD)nmrg, t1-t0);

	PetscFunctionReturn(0);
}


/*

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVExchange"
PetscErrorCode ADVExchange(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION EXCHANGE ROUTINE
	//
	// Exchange markers between the processors resulting from the position change
	//=======================================================================

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->advect == ADV_NONE) PetscFunctionReturn(0);


	// count number of markers to be sent to each neighbor domain
	ierr = ADVMapMarkToDomains(actx); CHKERRQ(ierr);

	// communicate number of markers with neighbor processes
	ierr = ADVExchangeNumMark(actx); CHKERRQ(ierr);

	// create send and receive buffers for asynchronous MPI communication
	ierr = ADVCreateMPIBuff(actx); CHKERRQ(ierr);

	// communicate markers with neighbor processes
	ierr = ADVExchangeMark(actx); CHKERRQ(ierr);

	// store received markers, collect garbage
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// free communication buffer
	ierr = ADVDestroyMPIBuff(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMapMarkToDomains"
PetscErrorCode ADVMapMarkToDomains(AdvCtx *actx)
{
	// count number of markers to be sent to each neighbor domain

	PetscInt     i, lrank, cnt;
	PetscMPIInt  grank;
	FDSTAG      *fs;

	PetscErrorCode  ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// clear send counters
	ierr = PetscMemzero(actx->nsendm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	// scan markers
	for(i = 0, cnt = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		if(grank == -1)
		{
			// currently all the markers must remain in the box
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "ERROR! Marker outflow is currently not implemented!");

			// otherwise, number of deleted markers should be updated here, i.e.:
			// cnt++;

			// we should also decide what to do with the outflow markers
			// periodic boundary conditions?
		}
		else if(grank != actx->iproc)
		{
			// count markers that should be sent to each neighbor
			actx->nsendm[lrank]++;
			cnt++;
		}
	}

	// store number of deleted markers
	actx->ndel = cnt;

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCreateMPIBuff"
PetscErrorCode ADVCreateMPIBuff(AdvCtx *actx)
{
	// create send and receive buffers for asynchronous MPI communication

	// NOTE! Currently the memory allocation model is fully dynamic.
	// Maybe it makes sense to introduce static model with reallocation.
	FDSTAG     *fs;
	PetscInt    i, cnt, lrank;
	PetscMPIInt grank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// compute buffer pointers
	actx->nsend = getPtrCnt(_num_neighb_, actx->nsendm, actx->ptsend);
	actx->nrecv = getPtrCnt(_num_neighb_, actx->nrecvm, actx->ptrecv);

	actx->sendbuf = NULL;
	actx->recvbuf = NULL;
	actx->idel    = NULL;

	// allocate exchange buffers & array of deleted (sent) marker indices
	if(actx->nsend) { ierr = PetscMalloc((size_t)actx->nsend*sizeof(Marker),   &actx->sendbuf); CHKERRQ(ierr); }
	if(actx->nrecv) { ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker),   &actx->recvbuf); CHKERRQ(ierr); }
	if(actx->ndel)  { ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel);    CHKERRQ(ierr); }

	// copy markers to send buffer, store their indices
	for(i = 0, cnt = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		if(grank == -1)
		{
			// currently this situation is prohibited in ADVMapMarkersDomains
			// but we already handle it here anyway

			// delete marker from the storage
			actx->idel[cnt++] = i;

		}
		else if(grank != actx->iproc)
		{
			// store marker in the send buffer
			actx->sendbuf[actx->ptsend[lrank]++] = actx->markers[i];

			// delete marker from the storage
			actx->idel[cnt++] = i;
		}
	}

	// rewind send buffer pointers
	rewindPtr(_num_neighb_, actx->ptsend);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVDestroyMPIBuff"
PetscErrorCode ADVDestroyMPIBuff(AdvCtx *actx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// destroy buffers
	ierr = PetscFree(actx->sendbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->idel);  	 CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCollectGarbage"
PetscErrorCode ADVCollectGarbage(AdvCtx *actx)
{
	// store received markers, collect garbage

	Marker   *markers, *recvbuf;
	PetscInt *idel, nummark, nrecv, ndel;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access storage
	nummark = actx->nummark;
	markers = actx->markers;

	nrecv   = actx->nrecv;
	recvbuf = actx->recvbuf;

	ndel    = actx->ndel;
	idel    = actx->idel;

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
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVUpdateMarkCell"
PetscErrorCode ADVUpdateMarkCell(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction

	FDSTAG      *fs;
	PetscInt    *numMarkCell, *m, i, p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// allocate marker counter array
	ierr = makeIntArray(&numMarkCell, NULL, fs->nCells); CHKERRQ(ierr);

	// count number of markers in the cells
	for(i = 0; i < actx->nummark; i++) numMarkCell[actx->cellnum[i]]++;

	// store starting indices of markers belonging to a cell
	actx->markstart[0] = 0;
	for(i = 1; i < fs->nCells+1; i++) actx->markstart[i] = actx->markstart[i-1]+numMarkCell[i-1];

	// allocate memory for id offset
	ierr = makeIntArray(&m, NULL, fs->nCells); CHKERRQ(ierr);

	// store marker indices belonging to a cell
	for(i = 0; i < actx->nummark; i++)
	{
		p = actx->markstart[actx->cellnum[i]];
		actx->markind[p + m[actx->cellnum[i]]] = i;
		m[actx->cellnum[i]]++;
	}

	// free memory
	ierr = PetscFree(m);           CHKERRQ(ierr);
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCheckCorners"
PetscErrorCode ADVCheckCorners(AdvCtx *actx)
{
	// check corner marker distribution
	// if empty insert one marker in center of corner
	FDSTAG         *fs;
	BCCtx          *bc;
	Marker         *P, *markers;
	NumCorner      *numcorner;
	PetscInt       i, j, ii, jj, kk, ind, nx, ny, nz, lx[3];
	PetscInt       n, p, I, J, K;
	PetscInt       ID, ineigh, Ii, Ji, Ki, indcell[27], nummark[27];
	PetscInt       ninj = 0, nind = 0, sind = 0;
	PetscScalar    dx[3], xp[3], xc[3], xs[3], xe[3], *X;
	PetscScalar    x, sumind = 0.0;
	PetscScalar    cf_rand;
	PetscRandom    rctx;
	PetscLogDouble t0,t1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bc = actx->jr->bc;

	ierr = PetscTime(&t0); CHKERRQ(ierr);

	fs     = actx->fs;
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nz     = fs->dsz.ncels;

	// allocate memory for corners
	ierr = PetscMalloc((size_t)fs->nCells*sizeof(NumCorner), &numcorner); CHKERRQ(ierr);
	ierr = PetscMemzero(numcorner, (size_t)fs->nCells*sizeof(NumCorner)); CHKERRQ(ierr);

	// count number of markers in corners
	for(i = 0; i < fs->nCells; i++)
	{
		// get I, J, K cell indices
		GET_CELL_IJK(i, I, J, K, nx, ny);

		// get coordinates of cell center
		xc[0] = fs->dsx.ccoor[I];
		xc[1] = fs->dsy.ccoor[J];
		xc[2] = fs->dsz.ccoor[K];

		// load markers in cell
		n = actx->markstart[i+1] - actx->markstart[i];
		p = actx->markstart[i];

		for(ii = 0; ii < n; ii++)
		{
			P = &actx->markers[actx->markind[p+ii]];

			// get marker coordinates
			xp[0] = P->X[0];
			xp[1] = P->X[1];
			xp[2] = P->X[2];

			// check position in cell
			if ((xp[0] < xc[0]) && (xp[1] < xc[1]) && (xp[2] < xc[2])) numcorner[i].s[0]++; // LSW
			if ((xp[0] > xc[0]) && (xp[1] < xc[1]) && (xp[2] < xc[2])) numcorner[i].s[1]++; // LSE
			if ((xp[0] < xc[0]) && (xp[1] > xc[1]) && (xp[2] < xc[2])) numcorner[i].s[2]++; // LNW
			if ((xp[0] > xc[0]) && (xp[1] > xc[1]) && (xp[2] < xc[2])) numcorner[i].s[3]++; // LNE

			if ((xp[0] < xc[0]) && (xp[1] < xc[1]) && (xp[2] > xc[2])) numcorner[i].s[4]++; // USW
			if ((xp[0] > xc[0]) && (xp[1] < xc[1]) && (xp[2] > xc[2])) numcorner[i].s[5]++; // USE
			if ((xp[0] < xc[0]) && (xp[1] > xc[1]) && (xp[2] > xc[2])) numcorner[i].s[6]++; // UNW
			if ((xp[0] > xc[0]) && (xp[1] > xc[1]) && (xp[2] > xc[2])) numcorner[i].s[7]++; // UNE
		}
	}

	// count how much to allocate memory for new markers
	for(i = 0; i < fs->nCells; i++)
	{
		for(j = 0; j < 8; j++)
		{
			if (numcorner[i].s[j] == 0) ninj++;
		}
	}

	// if no need for new markers
	if (!ninj)
	{
		// clear memory
		ierr = PetscFree(numcorner); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// initialize
	lx[0] = -1;
	lx[1] = 0;
	lx[2] = 1;

	// allocate memory for new markers
	actx->nrecv = ninj;
	ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->recvbuf, (size_t)actx->nrecv*sizeof(Marker)); CHKERRQ(ierr);

	// initialize the random number generator
	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

	// inject markers in corners
	for(i = 0; i < fs->nCells; i++)
	{
		for(j = 0; j < 8; j++)
		{
			if (numcorner[i].s[j] == 0)
			{
				// expand I, J, K cell indices
				GET_CELL_IJK(i, I, J, K, nx, ny);

				// get coordinates of cell center
				xc[0] = fs->dsx.ccoor[I];
				xc[1] = fs->dsy.ccoor[J];
				xc[2] = fs->dsz.ccoor[K];

				// get coordinates of cell node start
				xs[0] = fs->dsx.ncoor[I];
				xs[1] = fs->dsy.ncoor[J];
				xs[2] = fs->dsz.ncoor[K];

				// get coordinates of cell node end
				xe[0] = fs->dsx.ncoor[I+1];
				xe[1] = fs->dsy.ncoor[J+1];
				xe[2] = fs->dsz.ncoor[K+1];

				// GET MARKER NUMBER FROM LOCAL NEIGHBOURS
				ineigh = 0;
				n      = 0;
				for (kk = 0; kk < 3; kk++)
				{
					// get Ki index
					Ki = K + lx[kk];

					for (jj = 0; jj < 3; jj++)
					{
						// get Ji index
						Ji = J + lx[jj];

						for (ii = 0; ii < 3; ii++)
						{
							// get Ii index
							Ii = I + lx[ii];

							// check if within local domain bounds
							if ((Ki > -1) && (Ki < nz) && (Ji > -1) && (Ji < ny) && (Ii > -1) && (Ii < nx))
							{
								// get single index
								GET_CELL_ID(ID, Ii, Ji, Ki, nx, ny);

								// get markers number
								n += actx->markstart[ID+1] - actx->markstart[ID];

								// save markers number and single indices of cells
								indcell[ineigh] = ID;
								nummark[ineigh] = actx->markstart[ID+1] - actx->markstart[ID];
							}
							else
							{
								// save markers number and single indices of cells
								indcell[ineigh] = -1;
								nummark[ineigh] =  0;
							}

							// increase counter
							ineigh++;
						}
					}
				}

				// allocate memory for markers in cell
				ierr = PetscMalloc((size_t)n*sizeof(Marker),&markers); CHKERRQ(ierr);
				ierr = PetscMemzero(markers,(size_t)n*sizeof(Marker)); CHKERRQ(ierr);

				// load markers from all local neighbours
				jj = 0;
				for (ineigh = 0; ineigh < 27; ineigh++)
				{
					if ((indcell[ineigh]>-1) && (nummark[ineigh]>0))
					{
						for (ii = 0; ii < nummark[ineigh]; ii++)
						{
							// get index
							ind = actx->markind[actx->markstart[indcell[ineigh]] + ii];

							// save marker
							markers[jj] = actx->markers[ind];
							jj++;
						}
					}
				}

				// create coordinate of corner
				if (j == 0) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 1) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 2) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 3) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}

				if (j == 4) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 5) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 6) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 7) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}

				// add some random noise
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[0] += (cf_rand-0.5)*((xc[0]-xs[0])*0.5)*0.5;
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[1] += (cf_rand-0.5)*((xc[1]-xs[1])*0.5)*0.5;
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[2] += (cf_rand-0.5)*((xc[2]-xs[2])*0.5)*0.5;

				// calculate the closest (parent marker)
				for (ii = 0; ii < n; ii++)
				{
					X  = markers[ii].X;
					dx[0] = X[0] - xp[0];
					dx[1] = X[1] - xp[1];
					dx[2] = X[2] - xp[2];
					x  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

					if (ii == 0)    { sumind = x; sind = ii; }
					if (x < sumind) { sumind = x; sind = ii; }
				}

				// create new marker
				actx->recvbuf[nind]      = markers[sind];
				actx->recvbuf[nind].X[0] = xp[0];
				actx->recvbuf[nind].X[1] = xp[1];
				actx->recvbuf[nind].X[2] = xp[2];

				// override marker phase (if necessary)
				ierr = BCOverridePhase(bc, i, actx->recvbuf + nind); CHKERRQ(ierr);

				// increase counter
				nind++;

				// free memory
				ierr = PetscFree(markers); CHKERRQ(ierr);
			}
		}
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
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Marker control [%lld]: (Corners ) injected %lld markers in %1.4e s \n",(LLD)actx->iproc, (LLD)ninj, t1-t0);

	// clear
	ierr = PetscFree(numcorner);     CHKERRQ(ierr);
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

*/
//---------------------------------------------------------------------------
