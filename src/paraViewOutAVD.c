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
 **    filename:   AVDView.c
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
 *  Originally developed by Dave A. May on 6/21/11.
 *  Copyright 2011 Geophysical Fluid Dynamics. All rights reserved.
 *
 *  Adopted for use in LaMEM by Anton A. Popov
 *
 *  The algorithm computes an Approximate Voronoi Diagram (AVD) in 3D using a given set of point coordinates.
 *
 *  The AVD algorithm, is described in:
 *    M. Velic, D.A. May & L. Moresi,
 *    "A Fast Robust Algorithm for Computing Discrete Voronoi Diagrams",
 *    Journal of Mathematical Modelling and Algorithms,
 *    Volume 8, Number 3, 343-355, DOI: 10.1007/s10852-008-9097-6
 *
 *
 *  Notes:
 *    This implementation uses von-Neumann neighbourhoods for boundary chain growth.
 *    Do not be tempted to implement "diagonal" neighbourhood growth cycles - this will greatly increase the
 *    size of the boundary chain (and thus memory usage will increase and CPU time will decrease).
 */
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "phase.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "Tensor.h"
#include "advect.h"
#include "paraViewOutBin.h"
#include "paraViewOutAVD.h"
#include "tools.h"
//---------------------------------------------------------------------------
#define __AVD_DEBUG_MODE
//---------------------------------------------------------------------------
// ........................... AVDPoint3D ...................................
//---------------------------------------------------------------------------
void AVDPoint3DCreate(const PetscInt npoints, AVDPoint3D *P)
{
	AVDPoint3D points;

	points = (AVDPoint3D) malloc( sizeof(struct _p_AVDPoint3D)*(size_t)(npoints) );
	memset( points, 0, sizeof(struct _p_AVDPoint3D)*(size_t)(npoints) );

	*P = points;
}
//---------------------------------------------------------------------------
void AVDPoint3DDestroy(AVDPoint3D *P)
{
	AVDPoint3D points;

	if (!P) { return; }
	points = *P;
	free(points);
	*P = NULL;
}
//---------------------------------------------------------------------------
// ........................... AVDCell3D ....................................
//---------------------------------------------------------------------------
void AVDCell3DCreate(const PetscInt mx, const PetscInt my, const PetscInt mz, AVDCell3D *C)
{
	AVDCell3D cells;
	PetscInt i,j,k;

	cells = (AVDCell3D) malloc( sizeof(struct _p_AVDCell3D)*(size_t)(mx*my*mz) );
	memset( cells, 0, sizeof(struct _p_AVDCell3D)*(size_t)(mx*my*mz) );

	for (k=0; k<mz; k++) {
		for (j=0; j<my; j++) {
			for (i=0; i<mx; i++) {
				PetscInt ind;

				ind = i + j * mx + k * mx*my;
				cells[ind].index = ind;
				cells[ind].i = i;
				cells[ind].j = j;
				cells[ind].k = k;

				if ( (i==0) || (i==mx-1) ) { cells[ind].p = AVD_CELL_MASK; }
				if ( (j==0) || (j==my-1) ) { cells[ind].p = AVD_CELL_MASK; }
				if ( (k==0) || (k==mz-1) ) { cells[ind].p = AVD_CELL_MASK; }

			}
		}
	}

	*C = cells;
}
//---------------------------------------------------------------------------
void AVDCell3DDestroy(AVDCell3D *C)
{
	AVDCell3D cells;

	if (!C) { return; }
	cells = *C;
	free(cells);
	*C = NULL;
}

//---------------------------------------------------------------------------
// ........................... AVDChain3D ...................................
//---------------------------------------------------------------------------
void AVDChain3DCreate(const PetscInt npoints, const PetscInt buffer, AVDChain3D *CH)
{
	AVDChain3D chains;
	PetscInt p;

	chains = (AVDChain3D) malloc( sizeof(struct _p_AVDChain3D)*(size_t)(npoints) );
	memset( chains, 0, sizeof(struct _p_AVDChain3D)*(size_t)(npoints) );
	for (p=0; p<npoints; p++) {
		chains[p].new_claimed_cells_malloced = buffer;
		chains[p].new_boundary_cells_malloced = buffer;

		chains[p].new_claimed_cells = (PetscInt*) malloc (sizeof(PetscInt)*(size_t)(buffer+1) );
		chains[p].new_boundary_cells = (PetscInt*) malloc (sizeof(PetscInt)*(size_t)(buffer+1) );
	}

	*CH = chains;
}
//---------------------------------------------------------------------------
void AVDChain3DDestroy(const PetscInt npoints, AVDChain3D *CH)
{
	AVDChain3D chains;
	PetscInt p;

	if (!CH) { return; }
	chains = *CH;
	for (p=0; p<npoints; p++) {
		if (chains[p].new_claimed_cells) {
			free( chains[p].new_claimed_cells );
			chains[p].new_claimed_cells = NULL;
		}

		if (chains[p].new_boundary_cells) {
			free( chains[p].new_boundary_cells );
			chains[p].new_boundary_cells = NULL;
		}
	}
	free(chains);
	*CH = NULL;
}
//---------------------------------------------------------------------------
// ............................... AVD3D ....................................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVDViewCreate"
PetscErrorCode AVDViewCreate(AVD3D *A, AdvCtx *actx, PetscInt refine)
{
	AVD3D          avd3D;
	FDSTAG         *fs;
	PetscLogDouble tb, te;
	PetscScalar    bx, by, bz, ex, ey, ez;
	PetscInt       nx, ny, nz;
	PetscInt       i, count, claimed;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = actx->fs;

	ierr = PetscTime(&tb); CHKERRQ(ierr);

	// compute local grid resolution
	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	nx = refine*fs->dsx.ncels;
	ny = refine*fs->dsy.ncels;
	nz = refine*fs->dsz.ncels;

	// allocate space, setup structures
	AVD3DAllocate(nx, ny, nz, 1, actx->nummark, &avd3D);

	AVD3DSetDomainSize(avd3D, bx, ex, by, ey, bz, ez);

	ierr = AVD3DSetParallelExtent(avd3D, fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc); CHKERRQ(ierr);

	ierr = AVD3DLoadPoints(avd3D, actx); CHKERRQ(ierr);

	AVD3DResetCells(avd3D);

	ierr = AVD3DInit(avd3D); CHKERRQ(ierr);

	// assign all cells to points
	count   = avd3D->npoints;
	claimed = 1;

	while(claimed != 0)
	{
		claimed = 0;

		for(i = 0; i < count; i++)
		{
			AVD3DClaimCells(avd3D, i);

			claimed += avd3D->chains[i].num_claimed;

			AVD3DUpdateChain(avd3D, i);
		}
	}

	ierr = PetscTime(&te); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"AVD generator complete in %1.4e (sec)\n", te-tb);

	ierr = AVD3DReportMemory(avd3D); CHKERRQ(ierr);

	(*A) = avd3D;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void AVD3DDestroy(AVD3D *A)
{
	AVD3D aa;
	if (!A) { return; }
	aa = *A;
	if (aa->chains) {
		AVDChain3DDestroy(aa->npoints,&aa->chains);
	}
	if (aa->cells) {
		AVDCell3DDestroy(&aa->cells);
	}
	if (aa->points) {
		AVDPoint3DDestroy(&aa->points);
	}

	if (aa->ownership_ranges_i) {
		free(aa->ownership_ranges_i);
	}
	if (aa->ownership_ranges_j) {
		free(aa->ownership_ranges_j);
	}
	if (aa->ownership_ranges_k) {
		free(aa->ownership_ranges_k);
	}

	free(aa);
	*A = NULL;
}
//---------------------------------------------------------------------------
void AVD3DAllocate(
	const PetscInt mx,
	const PetscInt my,
	const PetscInt mz,
	const PetscInt buffer,
	const PetscInt npoints,
	AVD3D          *A)
{

	AVD3D avd3D;

	avd3D = (AVD3D) malloc( sizeof(struct _p_AVD3D) );
	memset( avd3D, 0, sizeof(struct _p_AVD3D) );

	avd3D->buffer = buffer;
	avd3D->mx = mx;
	avd3D->my = my;
	avd3D->mz = mz;

	avd3D->mx_mesh = mx+2;
	avd3D->my_mesh = my+2;
	avd3D->mz_mesh = mz+2;

	AVDCell3DCreate(
		(const PetscInt)avd3D->mx_mesh,
		(const PetscInt)avd3D->my_mesh,
		(const PetscInt)avd3D->mz_mesh,
		&avd3D->cells);

	avd3D->npoints = npoints;

	AVDPoint3DCreate(npoints, &avd3D->points);

	AVDChain3DCreate(npoints, (const PetscInt)avd3D->buffer, &avd3D->chains);

	(*A) = avd3D;
}
//---------------------------------------------------------------------------
PetscErrorCode AVD3DSetParallelExtent(AVD3D A, PetscInt M, PetscInt N, PetscInt P)
{
	PetscInt *tmp;
	PetscInt pid,i,j,k,sum;
	PetscErrorCode ierr;

	A->M = M;
	A->N = N;
	A->P = P;

	tmp = (PetscInt*) malloc(sizeof(PetscInt)*(size_t)(A->M*A->N*A->P+1));
	memset(tmp,0,sizeof(PetscInt)*(size_t)(A->M*A->N*A->P+1));

	A->ownership_ranges_i = (PetscInt*) malloc( sizeof(PetscInt)*(size_t)(A->M+1) );
	A->ownership_ranges_j = (PetscInt*) malloc( sizeof(PetscInt)*(size_t)(A->N+1) );
	A->ownership_ranges_k = (PetscInt*) malloc( sizeof(PetscInt)*(size_t)(A->P+1) );

	memset(tmp,0,sizeof(PetscInt)*(size_t)(A->M*A->N*A->P+1));
	ierr = MPI_Allgather(&A->mx,1,MPIU_INT,tmp,1,MPIU_INT,PETSC_COMM_WORLD); CHKERRQ(ierr);
	j = k = 0;
	sum = 0;
	for (i=0; i<A->M; i++) {
		pid = i + j*A->M + k*A->M*A->N;
		A->ownership_ranges_i[i] = sum;
		sum = sum + tmp[pid];
	} A->ownership_ranges_i[i] = sum;

	memset(tmp,0,sizeof(PetscInt)*(size_t)(A->M*A->N*A->P+1));
	ierr = MPI_Allgather(&A->my,1,MPIU_INT,tmp,1,MPIU_INT,PETSC_COMM_WORLD); CHKERRQ(ierr);
	i = k = 0;
	sum = 0;
	for (j=0; j<A->N; j++) {
		pid = i + j*A->M + k*A->M*A->N;
		A->ownership_ranges_j[j] = sum;
		sum = sum + tmp[pid];
	} A->ownership_ranges_j[j] = sum;

	memset(tmp,0,sizeof(PetscInt)*(size_t)(A->M*A->N*A->P+1));
	ierr = MPI_Allgather(&A->mz,1,MPIU_INT,tmp,1,MPIU_INT,PETSC_COMM_WORLD); CHKERRQ(ierr);
	i = j = 0;
	sum = 0;
	for (k=0; k<A->P; k++) {
		pid = i + j*A->M + k*A->M*A->N;
		A->ownership_ranges_k[k] = sum;
		sum = sum + tmp[pid];
	} A->ownership_ranges_k[k] = sum;

	A->gmx = A->ownership_ranges_i[A->M];
	A->gmy = A->ownership_ranges_j[A->N];
	A->gmz = A->ownership_ranges_k[A->P];

	free(tmp);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void AVD3DSetDomainSize(AVD3D A,
	const PetscScalar x0,
	const PetscScalar x1,
	const PetscScalar y0,
	const PetscScalar y1,
	const PetscScalar z0,
	const PetscScalar z1)
{
	A->x0 = x0;
	A->x1 = x1;
	A->y0 = y0;
	A->y1 = y1;
	A->z0 = z0;
	A->z1 = z1;

	A->dx = (x1-x0)/(PetscScalar)A->mx;
	A->dy = (y1-y0)/(PetscScalar)A->my;
	A->dz = (z1-z0)/(PetscScalar)A->mz;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVD3DLoadPoints"
PetscErrorCode AVD3DLoadPoints(AVD3D A, AdvCtx *actx)
{
	// load viewer points from markers
	// WARNING! Redundant copy of all markers coordinates & phases
	// can use markers directly

	Marker     *markers, *M;
	AVDPoint3D  points, P;
	PetscInt    i, npoints;

	PetscFunctionBegin;

	npoints = A->npoints;
	points  = A->points;
	markers = actx->markers;

	for(i = 0; i < npoints; i++)
	{
		M        = markers + i;
		P        = points  + i;
		P->x     = M->X[0];
		P->y     = M->X[1];
		P->z     = M->X[2];
		P->phase = M->phase;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void AVD3DResetCells(AVD3D A)
{
	PetscInt i,j,k;
	PetscInt mx,my,mz;

	mx = A->mx_mesh;
	my = A->my_mesh;
	mz = A->mz_mesh;

	for (k=0; k<mz; k++) {
		for (j=0; j<my; j++) {
			for (i=0; i<mx; i++) {
				PetscInt ii = i + j * mx + k * mx*my;

				A->cells[ii].p = -1;
				A->cells[ii].done = AVD_FALSE;

				if ( (i==0) || (i==mx-1) ) { A->cells[ii].p = AVD_CELL_MASK; }
				if ( (j==0) || (j==my-1) ) { A->cells[ii].p = AVD_CELL_MASK; }
				if ( (k==0) || (k==mz-1) ) { A->cells[ii].p = AVD_CELL_MASK; }

			}
		}
	}
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVD3DInit"
PetscErrorCode AVD3DInit(AVD3D A)
{
	// i = (xp - (x0-dx) )/mx_mesh

	AVDPoint3D points;
	PetscInt   p, i, j, k, npoints;
	PetscInt   mx, my, mz, ind;

	npoints = A->npoints;
	points  = A->points;

	mx = A->mx_mesh;
	my = A->my_mesh;
	mz = A->mz_mesh;

	for(p = 0; p < npoints; p++)
	{
		i = (PetscInt)((points[p].x - (A->x0 - A->dx))/A->dx);
		j = (PetscInt)((points[p].y - (A->y0 - A->dy))/A->dy);
		k = (PetscInt)((points[p].z - (A->z0 - A->dz))/A->dz);

		// if a particle is exactly on the border then make sure it is in a valid cell inside the element
		if (i == mx) { i--; }
		if (j == my) { j--; }
		if (k == mz) { k--; }

		// check bounds
		if (i==0)            SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: i==0:  %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);
		if (j==0)            SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: j==0:  %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);
		if (k==0)            SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: k==0:  %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);
		if (i==A->mx_mesh-1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: i==mx: %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);
		if (j==A->my_mesh-1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: j==my: %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);
		if (k==A->mz_mesh-1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: k==mz: %lf %lf %lf\n", points[p].x, points[p].y, points[p].z);

		ind = i+j*mx+k*mx*my;
		if (A->cells[ind].p == AVD_CELL_MASK)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "AVD3dInit: Inserting points into boundary cells - this is not permitted\n");
		}

		A->cells[i+j*mx+k*mx*my].p = p; // particle index
		A->chains[p].num_claimed = 1; // number of claimed cells, currently just the one the point initially resides within
		A->chains[p].length = 0;
		A->chains[p].total_claimed = 1; // total of claimed cells
		A->chains[p].done = AVD_FALSE;
		A->chains[p].index = i+j*mx+k*mx*my; // ith particle is in cell i +j*mx + k*mx*my
		A->chains[p].new_claimed_cells[0] = i+j*mx+k*mx*my;  // ith particle claimed cell it resides within, i.e. cell_index i+j*mx+k*mx*my
		A->chains[p].new_claimed_cells[1] = -1; // mark end of claimed_cells list with -1

		AVD3DUpdateChain(A, p);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static inline PetscScalar AVD3DDistanceTest(
	PetscScalar x0,
	PetscScalar y0,
	PetscScalar z0,
	PetscScalar x1,
	PetscScalar y1,
	PetscScalar z1,
	PetscScalar x2,
	PetscScalar y2,
	PetscScalar z2)
{
	return (x1+x2-x0-x0)*(x1-x2) + (y1+y2-y0-y0)*(y1-y2) + (z1+z2-z0-z0)*(z1-z2);
}
//---------------------------------------------------------------------------
void AVD3DClaimCells(AVD3D A, const PetscInt p_i)
{
	// Claim cells for particle p_i in the list

	PetscInt i,count;
	PetscScalar x0,y0,x1,y1,x2,y2,z0,z1,z2,dist1;
	PetscInt *temp;
	AVDChain3D bchain;
	AVDPoint3D points;
	AVDCell3D cells;
	PetscInt cell_num0;
	PetscInt buffer;
	PetscScalar dx,dy,dz;

	buffer = A->buffer;
	dx = A->dx;
	dy = A->dy;
	dz = A->dz;
	bchain = &A->chains[p_i];
	cells = A->cells;
	points = A->points;

	count = 0;
	bchain->num_claimed = 0;

	for (i=0; i<bchain->length; i++) {
		cell_num0 = bchain->new_boundary_cells[i]; // cell number we are trying to claim

#ifdef __AVD_DEBUG_MODE
		if (cell_num0<0) {
			printf("  AVD3dClaimCells(ERROR): p_i = %lld, [%lld] \n", (LLD)p_i,(LLD)cell_num0 );
			printf("  AVD3dClaimCells(ERROR):   point %f %f %f \n", A->points[p_i].x,A->points[p_i].y,A->points[p_i].z);
			exit(1);
		}

		if (cells[cell_num0].p == AVD_CELL_MASK) { printf("YOU SHOULD NEVER HAVE A MASKED CELL IN YOUR LIST\n"); exit(1); }
#endif

		if (cells[cell_num0].p == AVD_CELL_UNCLAIMED) { // if cell unclaimed, then claim it

// WARNING!!! NEVER use realloc! Either use malloc, or C++ containers

			// Realloc, note that we need one space more than the number of points to terminate the list
			if( count == bchain->new_claimed_cells_malloced-1  ){
				temp = (PetscInt*) realloc( bchain->new_claimed_cells, (size_t)(bchain->new_claimed_cells_malloced + buffer + 1)*sizeof(PetscInt) );
				bchain->new_claimed_cells = temp;
				bchain->new_claimed_cells_malloced += buffer;

				temp = (PetscInt*) realloc( bchain->new_boundary_cells, (size_t)(bchain->new_boundary_cells_malloced + buffer + 1)*sizeof(PetscInt) );
				bchain->new_boundary_cells = temp;
				bchain->new_boundary_cells_malloced += buffer;
			}
			bchain->new_claimed_cells[count] = cell_num0;
			bchain->num_claimed++;
			count++;
			cells[cell_num0].p = p_i; // mark cell as owned by particle p_i
		} else if (cells[cell_num0].p != p_i) {
			// perform distance test between points to determine ownership
			x2 = points[p_i].x;
			y2 = points[p_i].y;
			z2 = points[p_i].z;

			x1 = points[cells[cell_num0].p].x;
			y1 = points[cells[cell_num0].p].y;
			z1 = points[cells[cell_num0].p].z;

			// cell centroid
			x0 = (PetscScalar)cells[cell_num0].i*dx + (A->x0 - dx + 0.5*dx);
			y0 = (PetscScalar)cells[cell_num0].j*dy + (A->y0 - dy + 0.5*dy);
			z0 = (PetscScalar)cells[cell_num0].k*dz + (A->z0 - dz + 0.5*dz);

			dist1 = AVD3DDistanceTest(x0,y0,z0,x1,y1,z1,x2,y2,z2);
			if (dist1 > 0.0) {
				bchain->new_claimed_cells[count] = cell_num0;
				bchain->num_claimed++;
				count++;
				cells[cell_num0].p = p_i; // mark cell as owned by particle p_i
			}
		}
		bchain->new_claimed_cells[count] = -1; // mark end of list
	}
}
//---------------------------------------------------------------------------
void AVD3DUpdateChain(AVD3D A, const PetscInt p_i)
{
	PetscInt i,k;
	PetscInt count;
	PetscInt cell_num0,cell_num1,cell_num[6];
	AVDChain3D bchain;
	AVDCell3D cells,cell0;
	PetscInt mx,my,buffer;
	PetscInt *temp;

	buffer = A->buffer;
	mx = A->mx_mesh;
	my = A->my_mesh;
	bchain = &A->chains[p_i];
	cells = A->cells;

	count = 0;
	bchain->length = 0;
	for( i=0; i<bchain->num_claimed; i++) {
		cell_num0 = bchain->new_claimed_cells[i];
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
				cell_num[k] = -2;
			}
		}

		for (k=0; k<6; k++) {
			cell_num1 = cell_num[k];

			 // if cell does not already belong to the particle and hasn't been
			 // marked as being done then add it to new boundary array and mark it as done

			if (cell_num1 != -2) {
				if ( (cells[cell_num1].p != p_i) && (cells[cell_num1].done != AVD_TRUE) ) {

// WARNING!!! NEVER use realloc! Either use malloc, or C++ containers

					// Realloc, note that we need one space more than the number of points to terminate the list
					if (count == bchain->new_boundary_cells_malloced-1 ) {
						temp = (PetscInt*)realloc( bchain->new_claimed_cells, (size_t)(bchain->new_claimed_cells_malloced + buffer + 1)*sizeof(PetscInt) );
						bchain->new_claimed_cells = temp;
						bchain->new_claimed_cells_malloced += buffer;

						temp = (PetscInt*)realloc( bchain->new_boundary_cells, (size_t)(bchain->new_boundary_cells_malloced + buffer + 1)*sizeof(PetscInt) );
						bchain->new_boundary_cells = temp;
						bchain->new_boundary_cells_malloced += buffer;
					}
#ifdef __AVD_DEBUG_MODE
					if (cell_num1<0) {
						printf("  AVD3DUpdateChain(ERROR): INSERTING negative cell index \n");
						printf("  AVD3DUpdateChain(ERROR):   k=%lld :: cell0 i,j,k = %lld,%lld,%lld neighbourid [%lld]\n", (LLD)k,(LLD)(cell0->i), (LLD)(cell0->j), (LLD)(cell0->k), (LLD)cell_num1 );
						exit(0);
					}
#endif
					bchain->new_boundary_cells[count] = cell_num1;
					bchain->length++;
					count++;
					cells[cell_num1].done = AVD_TRUE;
				}
			}
		}

	}

	// reset the processed flags
	for (i=0; i<count; i++){
		cells[ bchain->new_boundary_cells[i] ].done = AVD_FALSE;
	}
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AVD3DReportMemory"
PetscErrorCode AVD3DReportMemory(AVD3D A)
{
	AVDChain3D  chains;
	PetscInt    i, npoints, ncells, nchains;
	PetscScalar lmem[3], gmem[3], sum;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	chains  = A->chains;

	// points memory [MB]
	npoints = A->npoints;
	lmem[0] = (PetscScalar)(sizeof(struct _p_AVDPoint3D)*(size_t)npoints)/1048576.0;

	// cells memory
	ncells = A->mx_mesh * A->my_mesh * A->mz_mesh;
	lmem[1] = (PetscScalar)(sizeof(struct _p_AVDCell3D)*(size_t)ncells)/1048576.0;

	// chains memory [MB]
	nchains = 0;
	for(i = 0; i < npoints; i++)
	{
		nchains += chains[i].new_boundary_cells_malloced
		+          chains[i].new_claimed_cells_malloced;
	}
	lmem[2] = (PetscScalar)(sizeof(struct _p_AVDChain3D)*(size_t)npoints
	+                       sizeof(PetscInt)            *(size_t)nchains)/1048576.0;

	// get global sum on rank zero
	gmem[0] = 0;
	gmem[1] = 0;
	gmem[2] = 0;

	ierr = MPI_Reduce(lmem, gmem, 3, MPIU_SCALAR, MPI_SUM, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

	sum  = gmem[0] + gmem[1] + gmem[2];

	// output
	PetscPrintf(PETSC_COMM_WORLD,"AVD3D memory usage: \n");
	PetscPrintf(PETSC_COMM_WORLD,"  points:  %f [MB] \n", gmem[0]);
	PetscPrintf(PETSC_COMM_WORLD,"  cells:   %f [MB] \n", gmem[1]);
	PetscPrintf(PETSC_COMM_WORLD,"  chains:  %f [MB] \n", gmem[2]);
	PetscPrintf(PETSC_COMM_WORLD,"  total:   %f [MB] \n", sum);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// .......................... AVD ParaView Output ...........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDCreate"
PetscErrorCode PVAVDCreate(PVAVD *pvavd, FB *fb)
{
	char filename[_STR_LEN_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	ierr = getIntParam(fb, _OPTIONAL_, "out_avd", &pvavd->outavd, 1, 1); CHKERRQ(ierr);

	if(!pvavd->outavd) PetscFunctionReturn(0);

	// initialize
	pvavd->outpvd = 1; // pvd file output flag
	pvavd->refine = 2; // Voronoi Diagram refinement factor

	// read
	ierr = getStringParam(fb, _OPTIONAL_, "out_file_name", filename,                  "output"); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_avd_pvd",   &pvavd->outpvd,                1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_avd_ref",   &pvavd->refine, 1, _max_avd_refine_); CHKERRQ(ierr);

	// print
	if(pvavd->outpvd)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Writing AVD phase .pvd file to disk\n");
	}

	// set file name
	sprintf(pvavd->outfile, "%s_phase", filename);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWriteTimeStep"
PetscErrorCode PVAVDWriteTimeStep(PVAVD *pvavd, const char *dirName, PetscScalar ttime, PetscInt tindx)
{
	// Create a 3D Voronoi diagram from particles with phase information
	// write the file to disk and perform scaling/unscaling of the variables

	AVD3D A;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!pvavd->outavd) PetscFunctionReturn(0);

	// create Approximate Voronoi Diagram from particles
	ierr = AVDViewCreate(&A, pvavd->actx, pvavd->refine); CHKERRQ(ierr);

	// update .pvd file if necessary
	if(pvavd->outpvd)
	{
		ierr = UpdatePVDFile(dirName, pvavd->outfile, "pvtr", &pvavd->offset, ttime, tindx); CHKERRQ(ierr);
	}

	ierr = PVAVDWritePVTR(pvavd, A, dirName); CHKERRQ(ierr);

	ierr = PVAVDWriteVTR(pvavd, A, dirName); CHKERRQ(ierr);

	// cleanup
	AVD3DDestroy(&A);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWritePVTR"
PetscErrorCode PVAVDWritePVTR(PVAVD *pvavd, AVD3D A, const char *dirName)
{
	FILE        *fp;
	char        *fname;
	PetscMPIInt inproc, irank;
	PetscInt    r2d, p, pi, pj, pk, nproc, rank;

	PetscFunctionBegin;

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	MPI_Comm_size(PETSC_COMM_WORLD, &inproc); nproc = (PetscInt)inproc;
	MPI_Comm_rank(PETSC_COMM_WORLD, &irank);  rank  = (PetscInt)irank;

	// open outfile.pvts file in the output directory (write mode)
	asprintf(&fname, "%s/%s.pvtr", dirName, pvavd->outfile);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	pk  = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj  = r2d/(A->M);
	pi  = r2d - pj*A->M;

	WriteXMLHeader(fp, "PRectilinearGrid");

	fprintf(fp, "  <PRectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\" GhostLevel=\"0\" >\n",
		0LL,(LLD)(A->gmx),
		0LL,(LLD)(A->gmy),
		0LL,(LLD)(A->gmz));

	fprintf(fp, "    <PCoordinates>\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "    </PCoordinates>\n");

	fprintf(fp, "    <PCellData>\n");

	fprintf(fp, "      <PDataArray type=\"UInt8\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" />\n");

	fprintf(fp, "    </PCellData>\n");

	fprintf(fp, "    <PPointData>\n");
	fprintf(fp, "    </PPointData>\n");

	for(p=0; p<nproc; p++)
	{
		pk = p/(A->M*A->N);
		r2d = p - pk*(A->M*A->N);
		pj = r2d/(A->M);
		pi = r2d - pj*A->M;

		fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" Source=\"%s_p%1.6lld.vtr\" />\n",
				(LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
				(LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
				(LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]),
				pvavd->outfile, (LLD)p );
	}

	fprintf(fp, "  </PRectilinearGrid>\n");

	fprintf(fp, "</VTKFile>\n");

	fclose( fp );

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWriteVTR"
PetscErrorCode PVAVDWriteVTR(PVAVD *pvavd, AVD3D A, const char *dirName)
{
	// WARNING! writing single entry at a time is too slow. Use buffers instead!

	PetscMPIInt   irank;
	FILE          *fp;
	char          *fname;
	PetscInt      i, j, k, ii;
	PetscInt      r2d, pi, pj, pk, rank;
	PetscScalar   chLen;
	float         crd;
	unsigned char phase;
	int           offset, L;

	PetscFunctionBegin;

	// access context
	chLen = pvavd->actx->jr->scal->length;

	MPI_Comm_rank(PETSC_COMM_WORLD, &irank);  rank = (PetscInt)irank;

	// open outfile_p_XXXXXX.vtr file in the output directory (write mode)
	asprintf(&fname, "%s/%s_p%1.6lld.vtr", dirName, pvavd->outfile, (LLD)rank);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	pk  = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj  = r2d/(A->M);
	pi  = r2d - pj*A->M;

	// write header
	WriteXMLHeader(fp, "RectilinearGrid");

  fprintf(fp, "  <RectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\" >\n",
		  (LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
		  (LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
		  (LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]));

	fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" >\n",
			(LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
			(LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
			(LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]));

	offset = 0;

	fprintf(fp, "    <Coordinates>\n");

	// X
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->mx+1));
	// Y
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->my+1));
	// Z
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->mz+1));

	fprintf(fp, "    </Coordinates>\n");

	fprintf(fp, "    <CellData>\n");

	// phase
	fprintf(fp, "      <DataArray type=\"UInt8\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);

	fprintf(fp, "    </CellData>\n");

	fprintf(fp, "    <PointData>\n");
	fprintf(fp, "    </PointData>\n");

	fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");

	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// X
	L = (int)sizeof(float)*(int)(A->mx+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->mx+1; i++ ) {
		crd = (float)((A->x0 + (PetscScalar)i*A->dx)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// Y
	L = (int)sizeof(float)*(int)(A->my+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->my+1; i++ ) {
		crd = (float)((A->y0 + (PetscScalar)i*A->dy)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// Z
	L = (int)sizeof(float)*(int)(A->mz+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->mz+1; i++ ) {
		crd = (float)((A->z0 + (PetscScalar)i*A->dz)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// phase
	L = (int)sizeof(unsigned char)*(int)(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(int), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++)
			{
				ii    = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				phase = (unsigned char)A->points[A->cells[ii].p].phase;
				fwrite(&phase,sizeof(unsigned char),1,fp);
			}
		}
	}
	fprintf(fp,"\n  </AppendedData>\n");

	fprintf(fp, "</VTKFile>\n");

	fclose( fp );

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
