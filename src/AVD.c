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
 **    filename:   AVD.c
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
//........   Routines based on the Approximate Voronoi Diagram (AVD)  .......
//---------------------------------------------------------------------------

// The algorithm computes an Approximate Voronoi Diagram (AVD) in 3D using a given set of point coordinates.
// The AVD algorithm, is described in:
//    M. Velic, D.A. May & L. Moresi (2008), "A Fast Robust Algorithm for Computing Discrete Voronoi Diagrams",

#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "AVD.h"
#include "tools.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDCreate"
PetscErrorCode AVDCreate(AVD *A)
{
	PetscInt    p, npoints;
	PetscInt    ind;
	PetscInt    i, j, k;
	PetscInt    mx,my,mz;
	PetscScalar x[3], dx[3];
	PetscScalar s[3];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize variables
	mx  = A->nx+2;
	my  = A->ny+2;
	mz  = A->nz+2;

	dx[0] = A->dx;
	dx[1] = A->dy;
	dx[2] = A->dz;

	s[0] = A->xs[0]-dx[0]*0.5;
	s[1] = A->xs[1]-dx[1]*0.5;
	s[2] = A->xs[2]-dx[2]*0.5;

	// --------------
	//   AVD CELLS
	// --------------
	// allocate memory for cells plus one layer of boundary cells
	ierr = PetscMalloc((size_t)(mx*my*mz)*sizeof(AVDCell), &A->cell); CHKERRQ(ierr);
	ierr = PetscMemzero(A->cell, (size_t)(mx*my*mz)*sizeof(AVDCell)); CHKERRQ(ierr);

	for (k=0; k<mz; k++)
	{
		// compute z - center coordinate
		x[2] = s[2] + (PetscScalar)k*dx[2];

		for (j=0; j<my; j++)
		{
			// compute y - center coordinate
			x[1] = s[1] + (PetscScalar)j*dx[1];

			for (i=0; i<mx; i++)
			{
				// compute x - center coordinate
				x[0] = s[0] + (PetscScalar)i*dx[0];

				ind = i + j * mx + k * mx*my;
				A->cell[ind].ind  = ind;
				A->cell[ind].i    = i;
				A->cell[ind].j    = j;
				A->cell[ind].k    = k;
				A->cell[ind].x[0] = x[0];
				A->cell[ind].x[1] = x[1];
				A->cell[ind].x[2] = x[2];
				A->cell[ind].done = PETSC_FALSE;
				A->cell[ind].p    = AVD_CELL_UNCLAIMED;
				A->cell[ind].col  = 0;

				// create boundary
				if ( (i==0) || (i==mx-1) ) { A->cell[ind].p = AVD_CELL_MASK; }
				if ( (j==0) || (j==my-1) ) { A->cell[ind].p = AVD_CELL_MASK; }
				if ( (k==0) || (k==mz-1) ) { A->cell[ind].p = AVD_CELL_MASK; }
			}
		}
	}

	// --------------
	//   AVD CHAIN
	// --------------
	A->buffer = 1;
	npoints   = A->npoints;

	// allocate memory for chains
	ierr = PetscMalloc((size_t)(npoints)*sizeof(AVDChain), &A->chain); CHKERRQ(ierr);
	ierr = PetscMemzero(A->chain, (size_t)(npoints)*sizeof(AVDChain)); CHKERRQ(ierr);

	for (p=0; p < npoints; p++)
	{
		// initialize dominant axis for half-centroid
		A->chain[p].xh[0] = 0.0;
		A->chain[p].xh[1] = 0.0;
		A->chain[p].xh[2] = 0.0;

		// allocate memory for chains
		A->chain[p].iclaim  = A->buffer;
		A->chain[p].ibound = A->buffer;

		ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(A->chain[p].iclaim + A->buffer), &A->chain[p].claim); CHKERRQ(ierr);
		ierr = PetscMemzero(A->chain[p].claim, sizeof(PetscInt)*(size_t)(A->chain[p].iclaim + A->buffer)); CHKERRQ(ierr);

		ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(A->chain[p].ibound + A->buffer), &A->chain[p].bound); CHKERRQ(ierr);
		ierr = PetscMemzero(A->chain[p].bound, sizeof(PetscInt)*(size_t)(A->chain[p].ibound + A->buffer)); CHKERRQ(ierr);
	}

	// --------------
	//   AVD POINTS
	// --------------
	// allocate memory for points
	ierr = PetscMalloc((size_t)(npoints)*sizeof(Marker), &A->points); CHKERRQ(ierr);
	ierr = PetscMemzero(A->points, (size_t)(npoints)*sizeof(Marker)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDDestroy"
PetscErrorCode AVDDestroy(AVD *A)
{
	PetscInt p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// --------------
	//   AVD CELLS
	// --------------
	ierr = PetscFree(A->cell); CHKERRQ(ierr);

	// --------------
	//   AVD CHAIN
	// --------------
	for (p = 0; p < A->npoints; p++)
	{
		if (A->chain[p].claim ) { ierr = PetscFree(A->chain[p].claim ); CHKERRQ(ierr); }
		if (A->chain[p].bound ) { ierr = PetscFree(A->chain[p].bound ); CHKERRQ(ierr); }
	}
	ierr = PetscFree(A->chain ); CHKERRQ(ierr);

	// --------------
	//   AVD POINTS
	// --------------
	ierr = PetscFree(A->points); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDCellInit"
PetscErrorCode AVDCellInit(AVD *A)
{
	Marker     *points;
	PetscInt    npoints;
	PetscInt    p,i,j,k;
	PetscInt    mx,my,mz,ind;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize variables
	points  = A->points;
	npoints = A->npoints;

	mx = A->nx+2;
	my = A->ny+2;
	mz = A->nz+2;

	// find positions of points inside Voronoi cells
	for (p = 0; p < npoints; p++){

		// compute cell index of the particles
		i = (PetscInt)((points[p].X[0] - (A->xs[0] - A->dx))/A->dx);
		j = (PetscInt)((points[p].X[1] - (A->xs[1] - A->dy))/A->dy);
		k = (PetscInt)((points[p].X[2] - (A->xs[2] - A->dz))/A->dz);

		// if a particle is exactly on the border then make sure it is in a valid cell inside the element
		if (i == mx-1) { i--; }
		if (j == my-1) { j--; }
		if (k == mz-1) { k--; }

		ind = i+j*mx+k*mx*my;

		if (A->cell[ind].p == AVD_CELL_MASK) {
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Inserting cells into boundary cells is not permitted \n");
		}

		A->cell[ind].p                   = p;         // particle index
		A->chain[p].nclaimed             = 1;         // number of claimed cells, currently just the one the point initially resides within
		A->chain[p].length               = 0;
		A->chain[p].done                 = PETSC_FALSE;
		A->chain[p].ind                  = ind;       // ith particle is in cell i +j*mx + k*mx*my
		A->chain[p].claim[0]             = ind;       // ith particle claimed cell it resides within, i.e. cell_index i+j*mx+k*mx*my
		A->chain[p].claim[1]             = -1;        // mark end of claimed_cells list with -1

		// update initial chain
		ierr = AVDUpdateChain(A,p); CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDClaimCells"
PetscErrorCode AVDClaimCells(AVD *A, const PetscInt ip)
{
	PetscInt    i,count;
	PetscScalar x0[3], x1[3], x2[3], dist;
	AVDChain    *bchain;
	AVDCell     *cells;
	Marker      *points;
	PetscInt    cell_num0;
	PetscInt    buffer;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	buffer = A->buffer;
	bchain = &A->chain[ip];
	cells  = A->cell;
	points = A->points;

	count  = 0;
	bchain->nclaimed = 0;

	for (i=0; i<bchain->length; i++) {
		cell_num0 = bchain->bound[i]; // cell number we are trying to claim

		// if cell unclaimed, then claim it
		if (cells[cell_num0].p == AVD_CELL_UNCLAIMED)
		{
			// re-alloc, note that we need one space more than the number of points to terminate the list
			if( count == bchain->iclaim-1 ) { ierr = AVDReAlloc(bchain, buffer); CHKERRQ(ierr); }

			// claim cell
			bchain->claim[count] = cell_num0;

			// update counters
			bchain->nclaimed++;
			count++;

			// mark cell as owned by particle ip
			cells[cell_num0].p = ip;
		}

		else if (cells[cell_num0].p != ip)
		{
			// perform distance test between points to determine ownership
			x2[0] = points[ip].X[0];
			x2[1] = points[ip].X[1];
			x2[2] = points[ip].X[2];

			x1[0] = points[cells[cell_num0].p].X[0];
			x1[1] = points[cells[cell_num0].p].X[1];
			x1[2] = points[cells[cell_num0].p].X[2];

			// cell centroid
			x0[0] = cells[cell_num0].x[0];
			x0[1] = cells[cell_num0].x[1];
			x0[2] = cells[cell_num0].x[2];

			dist = AVDDistanceTest(x0,x1,x2);
			if (dist > 0.0)
			{
				bchain->claim[count] = cell_num0;

				// update counters
				bchain->nclaimed++;
				count++;

				// mark cell as owned by particle ip
				cells[cell_num0].p = ip;
			}
		}

		// mark end of list
		bchain->claim[count] = -1;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDUpdateChain"
PetscErrorCode AVDUpdateChain(AVD *A, const PetscInt ip)
{
	PetscInt i,k;
	PetscInt count;
	PetscInt cell_num0,cell_num1,cell_num[6];
	AVDChain *bchain;
	AVDCell  *cells,*cell0;
	PetscInt mx,my,buffer;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	buffer = A->buffer;
	mx     = A->nx+2;
	my     = A->ny+2;
	bchain = &A->chain[ip];
	cells  = A->cell;

	count = 0;
	bchain->length = 0;
	for( i=0; i<bchain->nclaimed; i++) {
		cell_num0 = bchain->claim[i];
		cell0 = &cells[cell_num0];

		if (cell0->p == AVD_CELL_MASK) { continue; }

		cell_num[0] = (cell0->i  ) + (cell0->j-1)*mx + (cell0->k  )*mx*my; // S
		cell_num[1] = (cell0->i  ) + (cell0->j+1)*mx + (cell0->k  )*mx*my; // N
		cell_num[2] = (cell0->i+1) + (cell0->j  )*mx + (cell0->k  )*mx*my; // E
		cell_num[3] = (cell0->i-1) + (cell0->j  )*mx + (cell0->k  )*mx*my; // W
		cell_num[4] = (cell0->i  ) + (cell0->j  )*mx + (cell0->k+1)*mx*my; // Front
		cell_num[5] = (cell0->i  ) + (cell0->j  )*mx + (cell0->k-1)*mx*my; // Back

		// boundary protection
		for (k=0; k<6; k++) {
			if (cells[cell_num[k]].p == AVD_CELL_MASK) {
				cell_num[k] = AVD_CELL_MASK;
			}
		}

		for (k=0; k<6; k++)
		{
			cell_num1 = cell_num[k];

			// if cell does not belong to the particle add it to new boundary array and mark it as done
			if (cell_num1 != AVD_CELL_MASK)
			{
				if ( (cells[cell_num1].p != ip) && (!cells[cell_num1].done) )
				{
					// re-alloc - we need one space more than the number of points to terminate the list
					if (count == bchain->ibound-1 ) { ierr = AVDReAlloc(bchain, buffer); CHKERRQ(ierr); }

					// add new cell to boundary
					bchain->bound[count] = cell_num1;

					// increase counters
					bchain->length++;
					count++;

					// mark cell as done
					cells[cell_num1].done = PETSC_TRUE;
				}
			}
		}
	}

	// reset the processed flags
	for (i=0; i<count; i++)
	{
		cells[ bchain->bound[i] ].done = PETSC_FALSE;
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDReAlloc"
PetscErrorCode AVDReAlloc(AVDChain *chain, PetscInt buffer)
{
	PetscInt *temp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// 1. allocate memory for claimed cells
	ierr = makeIntArray(&temp, NULL, chain->iclaim + buffer); CHKERRQ(ierr);

	// copy current data
	ierr = PetscMemcpy(temp, chain->claim, (size_t)(chain->nclaimed + buffer)*sizeof(PetscInt)); CHKERRQ(ierr);

	// delete previous storage
	ierr = PetscFree(chain->claim); CHKERRQ(ierr);

	// save new capacity & storage
	chain->claim      = temp;
	chain->iclaim    += buffer;

	// 2. allocate memory for boundary cells
	ierr = makeIntArray(&temp, NULL, chain->ibound + buffer); CHKERRQ(ierr);

	// copy current data
	ierr = PetscMemcpy(temp, chain->bound, (size_t)(chain->length + buffer)*sizeof(PetscInt)); CHKERRQ(ierr);

	// delete previous storage
	ierr = PetscFree(chain->bound); CHKERRQ(ierr);

	// save new capacity & storage
	chain->bound      = temp;
	chain->ibound    += buffer;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDLoadPoints"
PetscErrorCode AVDLoadPoints(AdvCtx *actx, AVD *A, PetscInt ind)
{
	PetscInt    i, ii;
	PetscFunctionBegin;

	// load particles only within the Voronoi cell
	for (i = 0; i < A->npoints; i++)
	{
		// get index
		ii = actx->markind[actx->markstart[ind] + i];

		// save marker
		A->points[i] = actx->markers[ii];

		// save index
		A->chain [i].gind  = ii;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDInjectDeletePoints"
PetscErrorCode AVDInjectDeletePoints(AdvCtx *actx, AVD *A, PetscInt cellID)
{
	BCCtx      *bc;
	PetscInt    i, ii, n, ind;
	PetscInt    num_chain, hclaim;
	PetscInt    npoints, new_nmark = 0;
	PetscScalar xmin, xmax, ymin, ymax, zmin, zmax;
	PetscScalar xaxis, yaxis, zaxis;
	PetscScalar xp[3], xc[3], xh[3];
	PetscInt    *area, *sind, axis;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bc = actx->jr->bc;

	npoints = A->npoints;
	n  = (A->nx+2)*(A->ny+2)*(A->nz+2);
	

	// allocate memory to injected/deleted markers
	if      (npoints < A->mmin) new_nmark = A->mmin - npoints;
	else if (npoints > A->mmax) new_nmark = npoints - A->mmax;

	// allocate memory for sorting
	ierr = makeIntArray(&area, NULL, npoints); CHKERRQ(ierr);
	ierr = makeIntArray(&sind, NULL, npoints); CHKERRQ(ierr);
	

	// compute dominant axis
	for (i = 0; i < npoints; i++)
	{
		// initialize min and max with the origin of the particle
		xmin = A->points[i].X[0];
		xmax = A->points[i].X[0];
		ymin = A->points[i].X[1];
		ymax = A->points[i].X[1];
		zmin = A->points[i].X[2];
		zmax = A->points[i].X[2];

		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				if (A->cell[ii].x[0] < xmin) xmin = A->cell[ii].x[0];
				if (A->cell[ii].x[0] > xmax) xmax = A->cell[ii].x[0];
				if (A->cell[ii].x[1] < ymin) ymin = A->cell[ii].x[1];
				if (A->cell[ii].x[1] > ymax) ymax = A->cell[ii].x[1];
				if (A->cell[ii].x[2] < zmin) zmin = A->cell[ii].x[2];
				if (A->cell[ii].x[2] > zmax) zmax = A->cell[ii].x[2];
			}
		}

		// initialize dominant axis
		A->chain[i].axis = -1;

		// dominant axis
		xaxis = xmax-xmin;
		yaxis = ymax-ymin;
		zaxis = zmax-zmin;

		if ((xaxis > yaxis) && (xaxis > zaxis)) { A->chain[i].xh[0] = (xmax+xmin)*0.5; A->chain[i].axis = 0; }
		if ((yaxis > xaxis) && (yaxis > zaxis)) { A->chain[i].xh[1] = (ymax+ymin)*0.5; A->chain[i].axis = 1; }
		if ((zaxis > xaxis) && (zaxis > yaxis)) { A->chain[i].xh[2] = (zmax+zmin)*0.5; A->chain[i].axis = 2; }
	}

	// create colour - which cells to consider for the half-centroid
	for (i = 0; i < npoints; i++)
	{
		// half axis
		xh[0] = A->chain[i].xh[0];
		xh[1] = A->chain[i].xh[1];
		xh[2] = A->chain[i].xh[2];

		// point coordinate
		xp[0] = A->points[i].X[0];
		xp[1] = A->points[i].X[1];
		xp[2] = A->points[i].X[2];

		axis = A->chain[i].axis;

		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				// cell coordinate
				xc[0] = A->cell[ii].x[0];
				xc[1] = A->cell[ii].x[1];
				xc[2] = A->cell[ii].x[2];

				// mark coloring
				if (axis==-1) A->cell[ii].col = 1;
				else
				{
					if      ((xh[axis] < xp[axis]) && (xc[axis] < xh[axis])) A->cell[ii].col = 1;
					else if ((xh[axis] > xp[axis]) && (xc[axis] > xh[axis])) A->cell[ii].col = 1;
				}
			}
		}
	}

	// calculate half-centroid
	for (i = 0; i < npoints; i++)
	{
		hclaim = 0;
		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				// total claimed
				A->chain[i].tclaimed++;

				if (A->cell[ii].col == 1)
				{
					hclaim++;
					A->chain[i].xc[0] += A->cell[ii].x[0];
					A->chain[i].xc[1] += A->cell[ii].x[1];
					A->chain[i].xc[2] += A->cell[ii].x[2];
				}
			}
		}

		// centroid coordinates
		A->chain[i].xc[0] = A->chain[i].xc[0]/(PetscScalar)hclaim;
		A->chain[i].xc[1] = A->chain[i].xc[1]/(PetscScalar)hclaim;
		A->chain[i].xc[2] = A->chain[i].xc[2]/(PetscScalar)hclaim;

		// initialize variables for sorting
		sind[i] = i;
		area[i] = A->chain[i].tclaimed;
	}

	// sort in ascending order
	ierr = PetscSortIntWithArray(npoints,area,sind); CHKERRQ(ierr);
	

	// inject markers
	if      (npoints < A->mmin)
	{
		// do not insert more markers than available voronoi domains
		if (npoints < new_nmark) new_nmark = npoints;

		ind = npoints - 1;
		for (i = 0; i < new_nmark; i++)
		{
			num_chain = sind[ind];

			// inject same properties as parent marker except for position
			actx->recvbuf[actx->cinj+i]      = A->points[num_chain];
			actx->recvbuf[actx->cinj+i].X[0] = A->chain [num_chain].xc[0];
			actx->recvbuf[actx->cinj+i].X[1] = A->chain [num_chain].xc[1];
			actx->recvbuf[actx->cinj+i].X[2] = A->chain [num_chain].xc[2];

			//PetscPrintf(PETSC_COMM_SELF,"# Marker Control [%lld]: injected [%g,%g,%g]\n",(LLD)actx->iproc, A->chain [num_chain].xc[0], A->chain [num_chain].xc[1], A->chain [num_chain].xc[2]);

			// override marker phase (if necessary)
			ierr = BCOverridePhase(bc, cellID, actx->recvbuf + actx->cinj + i); CHKERRQ(ierr);

			ind--;
		}
		// update total counter
		actx->cinj +=new_nmark;
	}
	// delete markers
	else if (npoints > A->mmax)
	{
		ind = 0;
		for (i = 0; i < new_nmark; i++)
		{
			num_chain = sind[ind];
			actx->idel[actx->cdel+i] = A->chain[num_chain].gind;
			ind++;
		}
		// update total counter
		actx->cdel +=new_nmark;
	}
	

	// free memory
	ierr = PetscFree(area); CHKERRQ(ierr);
	ierr = PetscFree(sind); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDExecuteMarkerInjection"
PetscErrorCode AVDExecuteMarkerInjection(AdvCtx *actx, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind)
{

	AVD          A;
	PetscInt       i,claimed;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize some parameters
	A.nx = actx->avdx;
	A.ny = actx->avdy;
	A.nz = actx->avdz;

	A.mmin = actx->nmin;
	A.mmax = actx->nmax;

	A.npoints = npoints;

	A.xs[0] = xs[0];
	A.xs[1] = xs[1];
	A.xs[2] = xs[2];

	A.xe[0] = xe[0];
	A.xe[1] = xe[1];
	A.xe[2] = xe[2];

	A.dx = (xe[0]-xs[0])/(PetscScalar)A.nx;
	A.dy = (xe[1]-xs[1])/(PetscScalar)A.ny;
	A.dz = (xe[2]-xs[2])/(PetscScalar)A.nz;

	// AVD structures
	ierr = AVDCreate(&A); CHKERRQ(ierr);

	// load particles
	ierr = AVDLoadPoints(actx,&A,ind); CHKERRQ(ierr);

	// initialize AVD cells
	ierr = AVDCellInit(&A); CHKERRQ(ierr);

	// AVD algorithm
	claimed = 1;
	while (claimed!= 0)
	{
		claimed = 0;
		for (i = 0; i < npoints; i++)
		{
			ierr = AVDClaimCells(&A,i); CHKERRQ(ierr);
			claimed += A.chain[i].nclaimed;
			ierr = AVDUpdateChain(&A,i); CHKERRQ(ierr);
		}
	}

	// inject/delete markers
	ierr = AVDInjectDeletePoints(actx, &A, ind); CHKERRQ(ierr);

	// destroy AVD structure
	ierr = AVDDestroy(&A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// NEW MARKER CONTROL
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDMarkerControl"
PetscErrorCode AVDMarkerControl(AdvCtx *actx)
{
	// check marker distribution and delete or inject markers if necessary

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check if activated
	if(!actx->markContr) PetscFunctionReturn(0);

	PetscPrintf(PETSC_COMM_WORLD,"# NEW Marker Control Routine \n");

	// AVD routine for every control volume
	ierr = AVDMarkerControlMV(actx, _CELL_); CHKERRQ(ierr); // CELLS

	ierr = AVDMarkerControlMV(actx, _XYED_); CHKERRQ(ierr); // XY Edge

	ierr = AVDMarkerControlMV(actx, _XZED_); CHKERRQ(ierr); // XZ Edge

	ierr = AVDMarkerControlMV(actx, _YZED_); CHKERRQ(ierr); // YZ Edge

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDMarkerControlMV"
PetscErrorCode AVDMarkerControlMV(AdvCtx *actx, VolumeCase vtype)
{
	MarkerVolume  mv;
	PetscInt      dir = -1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if      (vtype == _CELL_) dir = -1;
	else if (vtype == _XYED_) dir =  2;
	else if (vtype == _XZED_) dir =  1;
	else if (vtype == _YZED_) dir =  0;

	// create MarkerVolume
	ierr = AVDCreateMV(actx, &mv, dir); CHKERRQ(ierr);

	// map markers
	ierr = AVDMapMarkersMV(actx, &mv, dir); CHKERRQ(ierr);

	// main marker control routine
	ierr = AVDCheckCellsMV(actx, &mv, dir); CHKERRQ(ierr);

	// free MarkerVolume
	ierr = AVDDestroyMV(&mv); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDCheckCellsMV"
PetscErrorCode AVDCheckCellsMV(AdvCtx *actx, MarkerVolume *mv, PetscInt dir)
{
	// check marker distribution and delete or inject markers if necessary
	PetscScalar    xs[3], xe[3];
	PetscInt       ind, i, j, k, M, N;
	PetscInt       n, ninj, ndel, nmin;
	PetscLogDouble t0,t1;
	char           lbl[_lbl_sz_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// record time
	ierr = PetscTime(&t0); CHKERRQ(ierr);

	// get number of cells
	M = mv->M;
	N = mv->N;

	// calculate storage
	ninj = 0;
	ndel = 0;
	for(ind = 0; ind < mv->ncells; ind++)
	{
		// no of markers in cell
		n = mv->markstart[ind+1] - mv->markstart[ind];

		if (n < actx->nmin)
		{
			// expand i, j, k cell indices
			GET_CELL_IJK(ind, i, j, k, M, N);

			// half-volumes
			nmin = actx->nmin;
			if ((dir == 0) && ((i == 0) | (i+1 == mv->M))) { nmin = (PetscInt) (actx->nmin/2+1); }
			if ((dir == 1) && ((j == 0) | (j+1 == mv->N))) { nmin = (PetscInt) (actx->nmin/2+1); }
			if ((dir == 2) && ((k == 0) | (k+1 == mv->P))) { nmin = (PetscInt) (actx->nmin/2+1); }

			if (n < nmin)
			{
				if ((nmin - n) > n) ninj += n;
				else                ninj += nmin - n;
			}
		}
		if (n > actx->nmax) ndel += n - actx->nmax;
	}

	// if no need for injection/deletion
	if ((!ninj) && (!ndel)) PetscFunctionReturn(0);

	actx->nrecv = ninj;
	actx->ndel  = ndel;

	// allocate memory
	if(ninj) { ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker),   &actx->recvbuf); CHKERRQ(ierr); }
	if(ndel) { ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel   ); CHKERRQ(ierr); }

	actx->cinj = 0;
	actx->cdel = 0;

	// inject/delete
	for(ind = 0; ind < mv->ncells; ind++)
	{
		// no of markers in cell
		n = mv->markstart[ind+1] - mv->markstart[ind];

		if ((n < actx->nmin) || (n > actx->nmax))
		{
			// expand i, j, k cell indices
			GET_CELL_IJK(ind, i, j, k, M, N);

			// get cell coordinates
			xs[0] = mv->xcoord[i]; xe[0] = mv->xcoord[i+1];
			xs[1] = mv->ycoord[j]; xe[1] = mv->ycoord[j+1];
			xs[2] = mv->zcoord[k]; xe[2] = mv->zcoord[k+1];

			// here calculate half volumes minimum
			nmin = actx->nmin;
			if ((dir == 0) && ((i == 0) | (i+1 == mv->M))) { nmin = (PetscInt) (actx->nmin/2+1); }
			if ((dir == 1) && ((j == 0) | (j+1 == mv->N))) { nmin = (PetscInt) (actx->nmin/2+1); }
			if ((dir == 2) && ((k == 0) | (k+1 == mv->P))) { nmin = (PetscInt) (actx->nmin/2+1); }

			// inject/delete markers
			if ((n < nmin) || (n > actx->nmax))
			{
				ierr = AVDAlgorithmMV(actx, mv, n, xs, xe, ind, nmin); CHKERRQ(ierr);
			}
		}
	}

	// store new markers
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// clear
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->idel);    CHKERRQ(ierr);

	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);

	if      (dir==-1) sprintf(lbl,"CELL");
	else if (dir== 0) sprintf(lbl,"XYED");
	else if (dir== 1) sprintf(lbl,"XZED");
	else if (dir== 2) sprintf(lbl,"YZED");

	PetscPrintf(PETSC_COMM_WORLD,"# Marker Control [%lld]: (AVD %s) injected %lld markers and deleted %lld markers in %1.4e s\n",(LLD)actx->iproc,lbl, (LLD)ninj, (LLD)ndel, t1-t0);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDMapMarkersMV"
PetscErrorCode AVDMapMarkersMV(AdvCtx *actx, MarkerVolume *mv, PetscInt dir)
{
	// creates arrays to optimize marker-cell interaction
	FDSTAG      *fs;
	PetscScalar *X;
	PetscInt     i, ID, I, J, K;
	PetscInt    *numMarkCell, *m, p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// create coordinate arrays
	if (dir == 0) // x-coord
	{
		mv->xcoord[0] = fs->dsx.ncoor[0]; mv->xcoord[mv->M] = fs->dsx.ncoor[fs->dsx.ncels];
		for(i = 1; i < mv->M; i++) mv->xcoord[i] = fs->dsx.ccoor[i-1];
	}
	else for(i = 0; i < mv->M+1; i++) mv->xcoord[i] = fs->dsx.ncoor[i];

	if (dir == 1) // y-coord
	{
		mv->ycoord[0] = fs->dsy.ncoor[0]; mv->ycoord[mv->N] = fs->dsy.ncoor[fs->dsy.ncels];
		for(i = 1; i < mv->N; i++) mv->ycoord[i] = fs->dsy.ccoor[i-1];
	}
	else for(i = 0; i < mv->N+1; i++) mv->ycoord[i] = fs->dsy.ncoor[i];

	if (dir == 2) // z-coord
	{
		mv->zcoord[0] = fs->dsz.ncoor[0]; mv->zcoord[mv->P] = fs->dsz.ncoor[fs->dsz.ncels];
		for(i = 1; i < mv->P; i++) mv->zcoord[i] = fs->dsz.ccoor[i-1];
	}
	else for(i = 0; i < mv->P+1; i++) mv->zcoord[i] = fs->dsz.ncoor[i];

	// loop over all local particles
	for(i = 0; i < actx->nummark; i++)
	{
		// get marker coordinates
		X = actx->markers[i].X;

		// find I, J, K indices by bisection algorithm
		I = FindPointInCell(mv->xcoord, 0, mv->M, X[0]);
		J = FindPointInCell(mv->ycoord, 0, mv->N, X[1]);
		K = FindPointInCell(mv->zcoord, 0, mv->P, X[2]);

		// compute and store consecutive index
		GET_CELL_ID(ID, I, J, K, mv->M, mv->N);

		mv->cellnum[i] = ID;

	}

	// allocate marker counter array
	ierr = makeIntArray(&numMarkCell, NULL, mv->ncells); CHKERRQ(ierr);

	// count number of markers in the cells
	for(i = 0; i < actx->nummark; i++) numMarkCell[mv->cellnum[i]]++;

	// store starting indices of markers belonging to a cell
	mv->markstart[0] = 0;
	for(i = 1; i < mv->ncells+1; i++) mv->markstart[i] = mv->markstart[i-1]+numMarkCell[i-1];

	// allocate memory for id offset
	ierr = makeIntArray(&m, NULL, mv->ncells); CHKERRQ(ierr);

	// store marker indices belonging to a cell
	for(i = 0; i < actx->nummark; i++)
	{
		p = mv->markstart[mv->cellnum[i]];
		mv->markind[p + m[mv->cellnum[i]]] = i;
		m[mv->cellnum[i]]++;
	}

	// free memory
	ierr = PetscFree(m);           CHKERRQ(ierr);
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDCreateMV"
PetscErrorCode AVDCreateMV(AdvCtx *actx, MarkerVolume *mv, PetscInt dir)
{
	// allocate memory and info to marker volume control structure
	FDSTAG      *fs;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	mv->ncells = 0;

	// get number of cells
	if (dir == 0) mv->M = fs->dsx.ncels+1; else mv->M = fs->dsx.ncels;
	if (dir == 1) mv->N = fs->dsy.ncels+1; else mv->N = fs->dsy.ncels;
	if (dir == 2) mv->P = fs->dsz.ncels+1; else mv->P = fs->dsz.ncels;

	// total number of cells
	mv->ncells = mv->M * mv->N * mv->P;

	// allocate memory for host cell numbers
	ierr = makeIntArray(&mv->cellnum, NULL, actx->markcap); CHKERRQ(ierr);

	// allocate memory for id marker arranging per cell
	ierr = makeIntArray(&mv->markind, NULL, actx->markcap); CHKERRQ(ierr);

	// memory for starting indices
	ierr = makeIntArray(&mv->markstart, NULL, mv->ncells+1); CHKERRQ(ierr);

	// allocate memory for local coordinates
	ierr = makeScalArray(&mv->xcoord, NULL, mv->M+1); CHKERRQ(ierr);
	ierr = makeScalArray(&mv->ycoord, NULL, mv->N+1); CHKERRQ(ierr);
	ierr = makeScalArray(&mv->zcoord, NULL, mv->P+1); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDDestroyMV"
PetscErrorCode AVDDestroyMV(MarkerVolume *mv)
{
	// free memory
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscFree(mv->cellnum);    CHKERRQ(ierr);
	ierr = PetscFree(mv->markind);    CHKERRQ(ierr);
	ierr = PetscFree(mv->markstart);  CHKERRQ(ierr);

	ierr = PetscFree(mv->xcoord);    CHKERRQ(ierr);
	ierr = PetscFree(mv->ycoord);    CHKERRQ(ierr);
	ierr = PetscFree(mv->zcoord);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDLoadPointsMV"
PetscErrorCode AVDLoadPointsMV(AdvCtx *actx, MarkerVolume *mv, AVD *A, PetscInt ind)
{
	PetscInt    i, ii;
	PetscFunctionBegin;

	// load particles only within the Voronoi cell
	for (i = 0; i < A->npoints; i++)
	{
		// get index
		ii = mv->markind[mv->markstart[ind] + i];

		// save marker
		A->points[i] = actx->markers[ii];

		// save index
		A->chain [i].gind  = ii;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDInjectPointsMV"
PetscErrorCode AVDInjectPointsMV(AdvCtx *actx, AVD *A)
{
	FDSTAG     *fs;
	BCCtx      *bc;
	PetscInt    i, ii, n, ind, I, J, K, cellID;
	PetscInt    num_chain, hclaim;
	PetscInt    npoints, new_nmark = 0;
	PetscScalar xmin, xmax, ymin, ymax, zmin, zmax;
	PetscScalar xaxis, yaxis, zaxis;
	PetscScalar xp[3], xc[3], xh[3];
	PetscInt    *area, *sind, axis;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bc = actx->jr->bc;
	fs = actx->fs;

	npoints = A->npoints;
	n  = (A->nx+2)*(A->ny+2)*(A->nz+2);

	// allocate memory for sorting
	ierr = makeIntArray(&area, NULL, npoints); CHKERRQ(ierr);
	ierr = makeIntArray(&sind, NULL, npoints); CHKERRQ(ierr);

	// compute dominant axis
	for (i = 0; i < npoints; i++)
	{
		// initialize min and max with the origin of the particle
		xmin = A->points[i].X[0];
		xmax = A->points[i].X[0];
		ymin = A->points[i].X[1];
		ymax = A->points[i].X[1];
		zmin = A->points[i].X[2];
		zmax = A->points[i].X[2];

		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				if (A->cell[ii].x[0] < xmin) xmin = A->cell[ii].x[0];
				if (A->cell[ii].x[0] > xmax) xmax = A->cell[ii].x[0];
				if (A->cell[ii].x[1] < ymin) ymin = A->cell[ii].x[1];
				if (A->cell[ii].x[1] > ymax) ymax = A->cell[ii].x[1];
				if (A->cell[ii].x[2] < zmin) zmin = A->cell[ii].x[2];
				if (A->cell[ii].x[2] > zmax) zmax = A->cell[ii].x[2];
			}
		}

		// initialize dominant axis
		A->chain[i].axis = -1;

		// dominant axis
		xaxis = xmax-xmin;
		yaxis = ymax-ymin;
		zaxis = zmax-zmin;

		if ((xaxis > yaxis) && (xaxis > zaxis)) { A->chain[i].xh[0] = (xmax+xmin)*0.5; A->chain[i].axis = 0; }
		if ((yaxis > xaxis) && (yaxis > zaxis)) { A->chain[i].xh[1] = (ymax+ymin)*0.5; A->chain[i].axis = 1; }
		if ((zaxis > xaxis) && (zaxis > yaxis)) { A->chain[i].xh[2] = (zmax+zmin)*0.5; A->chain[i].axis = 2; }
	}

	// create colour - which cells to consider for the half-centroid
	for (i = 0; i < npoints; i++)
	{
		// half axis
		xh[0] = A->chain[i].xh[0];
		xh[1] = A->chain[i].xh[1];
		xh[2] = A->chain[i].xh[2];

		// point coordinate
		xp[0] = A->points[i].X[0];
		xp[1] = A->points[i].X[1];
		xp[2] = A->points[i].X[2];

		axis = A->chain[i].axis;

		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				// cell coordinate
				xc[0] = A->cell[ii].x[0];
				xc[1] = A->cell[ii].x[1];
				xc[2] = A->cell[ii].x[2];

				// mark coloring
				if (axis==-1) A->cell[ii].col = 1;
				else
				{
					if      ((xh[axis] <= xp[axis]) && (xc[axis] <= xh[axis])) A->cell[ii].col = 1;
					else if ((xh[axis] >= xp[axis]) && (xc[axis] >= xh[axis])) A->cell[ii].col = 1;
				}
			}
		}
	}

	// calculate half-centroid
	for (i = 0; i < npoints; i++)
	{
		hclaim = 0;
		for (ii = 0; ii < n; ii++)
		{
			if (A->cell[ii].p == i)
			{
				// total claimed
				A->chain[i].tclaimed++;

				if (A->cell[ii].col == 1)
				{
					hclaim++;
					A->chain[i].xc[0] += A->cell[ii].x[0];
					A->chain[i].xc[1] += A->cell[ii].x[1];
					A->chain[i].xc[2] += A->cell[ii].x[2];
				}
			}
		}

		// centroid coordinates
		A->chain[i].xc[0] = A->chain[i].xc[0]/(PetscScalar)hclaim;
		A->chain[i].xc[1] = A->chain[i].xc[1]/(PetscScalar)hclaim;
		A->chain[i].xc[2] = A->chain[i].xc[2]/(PetscScalar)hclaim;

		// initialize variables for sorting
		sind[i] = i;
		area[i] = A->chain[i].tclaimed;
	}

	// sort in ascending order
	ierr = PetscSortIntWithArray(npoints,area,sind); CHKERRQ(ierr);

	// do not insert more markers than available voronoi domains
	new_nmark = A->mmin - npoints;
	if (npoints < new_nmark) new_nmark = npoints;

	ind = npoints - 1;
	for (i = 0; i < new_nmark; i++)
	{
		num_chain = sind[ind];

		// inject same properties as parent marker except for position
		actx->recvbuf[actx->cinj+i]      = A->points[num_chain];
		actx->recvbuf[actx->cinj+i].X[0] = A->chain [num_chain].xc[0];
		actx->recvbuf[actx->cinj+i].X[1] = A->chain [num_chain].xc[1];
		actx->recvbuf[actx->cinj+i].X[2] = A->chain [num_chain].xc[2];

		// print info
		//PetscPrintf(PETSC_COMM_SELF,"# Marker Control [%lld]: injected [%g,%g,%g]\n",(LLD)actx->iproc, A->chain [num_chain].xc[0], A->chain [num_chain].xc[1], A->chain [num_chain].xc[2]);

		// --- this is not ideal with multiple control volumes (i.e. use mv for BCOverridePhase) ---
		// find I, J, K indices by bisection algorithm
		I = FindPointInCell(fs->dsx.ncoor, 0, fs->dsx.ncels, actx->recvbuf[actx->cinj+i].X[0]);
		J = FindPointInCell(fs->dsy.ncoor, 0, fs->dsy.ncels, actx->recvbuf[actx->cinj+i].X[1]);
		K = FindPointInCell(fs->dsz.ncoor, 0, fs->dsz.ncels, actx->recvbuf[actx->cinj+i].X[2]);

		// compute and store consecutive index
		GET_CELL_ID(cellID, I, J, K, fs->dsx.ncels, fs->dsy.ncels);

		// override marker phase (if necessary) - need to calculate cellID
		ierr = BCOverridePhase(bc, cellID, actx->recvbuf + actx->cinj + i); CHKERRQ(ierr);
		// -----------------------------------------------------------------------------------------

		ind--;
	}
	// update total counter
	actx->cinj +=new_nmark;

	// free memory
	ierr = PetscFree(area); CHKERRQ(ierr);
	ierr = PetscFree(sind); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDDeletePointsMV"
PetscErrorCode AVDDeletePointsMV(AdvCtx *actx, AVD *A)
{
	PetscInt    i, ind;
	PetscInt    num_chain;
	PetscInt    npoints, new_nmark = 0;
	PetscInt    *area, *sind;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	npoints = A->npoints;
	new_nmark = npoints - A->mmax;

	// allocate memory for sorting
	ierr = makeIntArray(&area, NULL, npoints); CHKERRQ(ierr);
	ierr = makeIntArray(&sind, NULL, npoints); CHKERRQ(ierr);

	// initialize variables for sorting
	for (i = 0; i < npoints; i++)
	{
		sind[i] = i;
		area[i] = A->chain[i].tclaimed;
	}

	// sort in ascending order
	ierr = PetscSortIntWithArray(npoints,area,sind); CHKERRQ(ierr);

		ind = 0;
		for (i = 0; i < new_nmark; i++)
		{
			num_chain = sind[ind];
			actx->idel[actx->cdel+i] = A->chain[num_chain].gind;
			ind++;

			// print info
			//PetscPrintf(PETSC_COMM_SELF,"# Marker Control [%lld]: deleted [%g,%g,%g]\n",(LLD)actx->iproc, A->chain [num_chain].xc[0], A->chain [num_chain].xc[1], A->chain [num_chain].xc[2]);
		}
		// update total counter
		actx->cdel +=new_nmark;

	// free memory
	ierr = PetscFree(area); CHKERRQ(ierr);
	ierr = PetscFree(sind); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDAlgorithmMV"
PetscErrorCode AVDAlgorithmMV(AdvCtx *actx, MarkerVolume *mv, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind, PetscInt nmin)
{

	AVD          A;
	PetscInt     i,claimed;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize some parameters
	A.nx = actx->avdx;
	A.ny = actx->avdy;
	A.nz = actx->avdz;

	A.mmin = nmin;
	A.mmax = actx->nmax;

	A.npoints = npoints;

	A.xs[0] = xs[0];
	A.xs[1] = xs[1];
	A.xs[2] = xs[2];

	A.xe[0] = xe[0];
	A.xe[1] = xe[1];
	A.xe[2] = xe[2];

	A.dx = (xe[0]-xs[0])/(PetscScalar)A.nx;
	A.dy = (xe[1]-xs[1])/(PetscScalar)A.ny;
	A.dz = (xe[2]-xs[2])/(PetscScalar)A.nz;

	// create AVD structure
	ierr = AVDCreate(&A); CHKERRQ(ierr);

	// load particles
	ierr = AVDLoadPointsMV(actx,mv,&A,ind); CHKERRQ(ierr);

	// initialize AVD cells
	ierr = AVDCellInit(&A); CHKERRQ(ierr);

	// do AVD algorithm
	claimed = 1;
	while (claimed!= 0)
	{
		claimed = 0;
		for (i = 0; i < npoints; i++)
		{
			ierr = AVDClaimCells(&A,i); CHKERRQ(ierr);
			claimed += A.chain[i].nclaimed;
			ierr = AVDUpdateChain(&A,i); CHKERRQ(ierr);
		}
	}

	// inject markers
	if (A.npoints < A.mmin) { ierr = AVDInjectPointsMV(actx, &A); CHKERRQ(ierr); }

	// delete markers
	if (A.npoints > A.mmax) { ierr = AVDDeletePointsMV(actx, &A); CHKERRQ(ierr); }

	// destroy AVD structure
	ierr = AVDDestroy(&A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
