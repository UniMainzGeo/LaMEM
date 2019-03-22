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
	// check marker distribution and merge or inject markers based on subgrid

	FDSTAG           *fs;
	PetscInt          cellid, markid, icell, I, J, K, i, j, ib, ie;
	PetscInt          ncx, ncy, npx, npy, npz, nmark, ncell;
	PetscInt         *markind;
	PetscScalar       s[3], h[3], xc[3], *x;
	PetscLogDouble    t0, t1;
    ipair             t;
    spair             d;
    Marker            P;
	vector <Marker>   inject;
	vector <PetscInt> imerge;
	vector <ipair>    cell;
    vector <spair>    dist;

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
	cell.reserve(1000);
	dist.reserve(1000);

	inject.reserve(actx->nummark/10);
	imerge.reserve(actx->nummark/10);

	inject.clear();
	imerge.clear();

	// process local cells
	for(icell = 0; icell < fs->nCells; icell++)
	{
		// get number of markers per cell
		nmark = actx->markstart[icell+1] - actx->markstart[icell];

		if(!nmark)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, " Empty control volume");
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
			markid = markind[j];
			x      = actx->markers[markid].X;

			// compute containing subcell index
			MAP_SUBCELL(I, x[0], s[0], h[0], npx);
			MAP_SUBCELL(J, x[1], s[1], h[1], npy);
			MAP_SUBCELL(K, x[2], s[2], h[2], npz);

			GET_CELL_ID(cellid, I, J, K, npx, npy);

			// store marker and cell numbers
			t.first  = cellid;
			t.second = markid;
			cell.push_back(t);

		}

		// sort markers by subcells
		sort(cell.begin(), cell.end());

		// push end-of-array stamp
		t.first  = -1;
		t.second =  0;
		cell.push_back(t);

		// process cells
		for(i = 0, ib = 0; i < ncell; i++)
		{
			// get index of next populated cell
			cellid = cell[ib].first;

			// check whether current cell is empty
			if(cellid != i)
			{
				// expand I, J, K indices
				GET_CELL_IJK(i, I, J, K, npx, npy)

				// get coordinates of cell center
				COORD_SUBCELL(xc[0], I, s[0], h[0]);
				COORD_SUBCELL(xc[1], J, s[1], h[1]);
				COORD_SUBCELL(xc[2], K, s[2], h[2]);

				// clear distance storage
				dist.clear();

				// find closest marker (cell-wise approximation)
				for(j = 0; j < nmark; j++)
				{
					// get marker index & coordinates
					markid = markind[j];
					x      = actx->markers[markid].X;

					// store marker and distance
					d.first  = EDIST(x, xc);
					d.second = markid;
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

				// store injected marker
				inject.push_back(P);
			}
			else
			{
				// find next populated cell or end-of-array stamp
				ie = ib; while(cell[ie].first == cellid) ie++;




				// merge markers





				// switch to next populated cell
				ib = ie;
			}
		}
	}

	// rearrange storage after marker resampling
	ierr = ADVCollectGarbageVec(actx, inject, imerge); CHKERRQ(ierr);

	// compute host cells for all the markers
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// update arrays for marker-cell interaction
	ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Marker control [%lld]: (subgrid) injected %lld markers and merged %lld markers in %1.4e s\n",(LLD)actx->iproc, (LLD)inject.size(), (LLD)imerge.size(), t1-t0);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

void checkMergeMarkers(
	Marker            *markers,
	vector <PetscInt> &imerge,
	vector <ipair>    &cell,
	PetscInt           npmax,
	PetscInt           ib,
	PetscInt           ie)
{
	PetscInt j, jb, je, phaseid;


//	1. more than npmax markers of each phase
//	2. find and merge two closest markers of each phase (recursively)

	// replace marker cell index with phase ID
	for(j = ib; j < ie; j++)
	{
		cell[j].first = markers[cell[j].second].phase;
	}

	// sort markers by phase
	sort(cell.begin() + ib, cell.begin() + ie);


	jb = ib;

	do
	{
		// get next phase ID
		phaseid = cell[jb].first;


		je = jb; while(cell[je].first == phaseid && je < ie) je++;

		if((je-jb) > npmax)
		{


	imerge.push_back(99); // **************** ALITA BATTLE ANGEL !!! *********************


	// EDIST



		}



	} while(je < ie);





/*

					// collect marker phases

					for(jj = je, jb = 0; ii < ncell; ii++)

					// process non-empty subcell (merge markers of same phase)
					if(je-jb > actx->npmax)
					{
						// ...............

						nmrg++;

					}

*/


}

//---------------------------------------------------------------------------

void MergeMarkers(
	Marker            *markers,
	vector <PetscInt> &imerge,
	vector <ipair>    &cell,
	PetscInt           jb,
	PetscInt           je)
{
	// merge two closest markers


	Marker      *pj, *pk;
	PetscInt     j, k, jmin, kmin;
	PetscScalar *xj, *xk, *uj, *uk, d, dmin;




	dmin = DBL_MAX;
	jmin = 0;
	kmin = 0;

	for(j = jb; j < je; j++)
	{
		xj = markers[cell[j].second].X;

		for(k = j+1; k < je; k++)
		{
			xk = markers[cell[k].second].X;

			d = EDIST(xj, xk);

			if(d < dmin)
			{
				dmin = d;
				jmin = j;
				kmin = k;
			}
		}
	}

	// merge markers
	pj = &markers[cell[jmin].second];
	pk = &markers[cell[kmin].second];
	xj = pj->X;
	xk = pk->X;
	uj = pj->U;
	uk = pk->U;


	xj[0]   = (xj[0]   + xk[0])  /2.0;
	xj[1]   = (xj[1]   + xk[1])  /2.0;
	xj[2]   = (xj[2]   + xk[2])  /2.0;
	pj->p   = (pj->p   + pk->p)  /2.0;
	pj->T   = (pj->T   + pk->T)  /2.0;
	pj->APS = (pj->APS + pk->APS)/2.0;
	pj->ATS = (pj->ATS + pk->ATS)/2.0;
	uj[0]   = (uj[0]   + uk[0])  /2.0;
	uj[1]   = (uj[1]   + uk[1])  /2.0;
	uj[2]   = (uj[2]   + uk[2])  /2.0;

	Tensor2RSSum2(&pj->S, 0.5, &pk->S, 0.5, &pj->S);

	// store merged marker
	imerge.push_back(cell[kmin].second);
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


/*
//---------------------------------------------------------------------------

		cell.clear();

		for(ii = 0; ii < 3; ii++)
		{
			t.first  = 5;
			t.second = 99;
			cell.push_back(t);
		}

		for(ii = 0; ii < 5; ii++)
		{
			t.first  = 1;
			t.second = 99;
			cell.push_back(t);
		}

		for(ii = 0; ii < 7; ii++)
		{
			t.first  = 3;
			t.second = 99;
			cell.push_back(t);
		}

		sort(cell.begin(), cell.end());

		t.first  = -1;
		t.second =  0;
		cell.push_back(t);

//---------------------------------------------------------------------------
*/

