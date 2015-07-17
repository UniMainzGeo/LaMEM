//---------------------------------------------------------------------------
//........   Routines based on the Approximate Voronoi Diagram (AVD)  .......
//---------------------------------------------------------------------------

// The algorithm computes an Approximate Voronoi Diagram (AVD) in 3D using a given set of point coordinates.
// The AVD algorithm, is described in:
//    M. Velic, D.A. May & L. Moresi (2008), "A Fast Robust Algorithm for Computing Discrete Voronoi Diagrams",

#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "AVD.h"
#include "Utils.h"

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
		x[2] = s[2] + k*dx[2];

		for (j=0; j<my; j++)
		{
			// compute y - center coordinate
			x[1] = s[1] + j*dx[1];

			for (i=0; i<mx; i++)
			{
				// compute x - center coordinate
				x[0] = s[0] + i*dx[0];

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
PetscErrorCode AVDInjectDeletePoints(AdvCtx *actx, AVD *A)
{
	PetscInt    i, ii, n, ind;
	PetscInt    num_chain, hclaim;
	PetscInt    npoints, new_nmark = 0;
	PetscScalar xmin, xmax, ymin, ymax, zmin, zmax;
	PetscScalar xaxis, yaxis, zaxis;
	PetscScalar xp[3], xc[3], xh[3];
	PetscInt    *area, *sind, axis;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		A->chain[i].xc[0] = A->chain[i].xc[0]/hclaim;
		A->chain[i].xc[1] = A->chain[i].xc[1]/hclaim;
		A->chain[i].xc[2] = A->chain[i].xc[2]/hclaim;

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

	A.dx = (xe[0]-xs[0])/A.nx;
	A.dy = (xe[1]-xs[1])/A.ny;
	A.dz = (xe[2]-xs[2])/A.nz;

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
	ierr = AVDInjectDeletePoints(actx, &A); CHKERRQ(ierr);

	// destroy AVD structure
	AVDDestroy(&A);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
