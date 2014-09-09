/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

AVDPhaseViewer.c

Created by Dave A. May on 6/21/11.
Copyright 2011 Geophysical Fluid Dynamics. All rights reserved.

The algorithm computes an Approximate Voronoi Diagram (AVD) in 3D using a given set of point coordinates.

The AVD algorithm, is described in:
    M. Velic, D.A. May & L. Moresi,
    "A Fast Robust Algorithm for Computing Discrete Voronoi Diagrams",
    Journal of Mathematical Modelling and Algorithms,
    Volume 8, Number 3, 343-355, DOI: 10.1007/s10852-008-9097-6

Notes:
    This implementation uses von-Neumann neighbourhoods for boundary chain growth.
    Do not be tempted to implement "diagonal" neighbourhood growth cycles - this will greatly increase the
    size of the boundary chain (and thus memory usage will increase and CPU time will decrease).

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
 */

#include "LaMEM.h"
#include "AVDPhaseViewer.h"
#include "Elements.h"

#define __AVD_DEBUG_MODE

void _AVDCell3dCreate(const PetscInt mx,const PetscInt my, const PetscInt mz,AVDCell3d *C)
{
	AVDCell3d cells;
	PetscInt i,j,k;

	cells = (AVDCell3d) malloc( sizeof(struct _p_AVDCell3d)*(size_t)(mx*my*mz) );
	memset( cells, 0, sizeof(struct _p_AVDCell3d)*(size_t)(mx*my*mz) );

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

void AVDCell3dDestroy(AVDCell3d *C)
{
	AVDCell3d cells;

	if (!C) { return; }
	cells = *C;
	free(cells);
	*C = NULL;
}

void AVDCell3dReset(AVD3d A)
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

void _AVDChain3dCreate(const PetscInt npoints, const PetscInt buffer,AVDChain3d *CH)
{
	AVDChain3d chains;
	PetscInt p;

	chains = (AVDChain3d) malloc( sizeof(struct _p_AVDChain3d)*(size_t)(npoints) );
	memset( chains, 0, sizeof(struct _p_AVDChain3d)*(size_t)(npoints) );
	for (p=0; p<npoints; p++) {
		chains[p].new_claimed_cells_malloced = buffer;
		chains[p].new_boundary_cells_malloced = buffer;

		chains[p].new_claimed_cells = (PetscInt*) malloc (sizeof(PetscInt)*(size_t)buffer );
		chains[p].new_boundary_cells = (PetscInt*) malloc (sizeof(PetscInt)*(size_t)buffer );
	}

	*CH = chains;
}

void AVDChain3dDestroy(const PetscInt npoints,AVDChain3d *CH)
{
	AVDChain3d chains;
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

void AVDPoint3dCreate(const PetscInt npoints, AVDPoint3d *P)
{
	AVDPoint3d points;

	points = (AVDPoint3d) malloc( sizeof(struct _p_AVDPoint3d)*(size_t)(npoints) );
	memset( points, 0, sizeof(struct _p_AVDPoint3d)*(size_t)(npoints) );

	*P = points;
}

void AVDPoint3dDestroy(AVDPoint3d *P)
{
	AVDPoint3d points;

	if (!P) { return; }
	points = *P;
	free(points);
	*P = NULL;
}

/*
 i = (xp - (x0-dx) )/mx_mesh
*/
void AVD3dCreate(const PetscInt mx,const PetscInt my, const PetscInt mz,const PetscInt buffer,AVD3d *A)
{
	AVD3d avd3d;

	avd3d = (AVD3d) malloc( sizeof(struct _p_AVD3d) );
	memset( avd3d, 0, sizeof(struct _p_AVD3d) );

	avd3d->buffer = buffer;
	avd3d->mx = mx;
	avd3d->my = my;
	avd3d->mz = mz;

	avd3d->mx_mesh = mx+2;
	avd3d->my_mesh = my+2;
	avd3d->mz_mesh = mz+2;

	_AVDCell3dCreate((const PetscInt)avd3d->mx_mesh,(const PetscInt)avd3d->my_mesh,(const PetscInt)avd3d->mz_mesh, &avd3d->cells);

	*A = avd3d;

}

PetscErrorCode AVD3dSetParallelExtent(AVD3d A,PetscInt M,PetscInt N,PetscInt P)
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

void AVD3dDestroy(AVD3d *A)
{
	AVD3d aa;
	if (!A) { return; }
	aa = *A;
	if (aa->chains) {
		AVDChain3dDestroy(aa->npoints,&aa->chains);
	}
	if (aa->cells) {
		AVDCell3dDestroy(&aa->cells);
	}
	if (aa->points) {
		AVDPoint3dDestroy(&aa->points);
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

void AVD3dSetDomainSize(AVD3d A,const PetscScalar x0,const PetscScalar x1,const PetscScalar y_0,const PetscScalar y_1,const PetscScalar z0,const PetscScalar z1)
{
	A->x0 = x0;
	A->x1 = x1;
	A->y0 = y_0;
	A->y1 = y_1;
	A->z0 = z0;
	A->z1 = z1;

	A->dx = (x1-x0)/(PetscScalar)A->mx;
	A->dy = (y_1-y_0)/(PetscScalar)A->my;
	A->dz = (z1-z0)/(PetscScalar)A->mz;
}

void AVD3dSetPoints(AVD3d A,const PetscInt npoints,AVDPoint3d points)
{

	if (A->chains) {
		printf("Deallocating existing chains\n");
		AVDChain3dDestroy(A->npoints,&A->chains);
	}
	_AVDChain3dCreate(npoints,(const PetscInt)A->buffer,&A->chains);

	A->npoints = npoints;
	A->points = points;
}

static inline PetscScalar AVD3dDistanceTest(PetscScalar x0,PetscScalar y_0,PetscScalar z0,PetscScalar x1,PetscScalar y_1,PetscScalar z1,PetscScalar x2,PetscScalar y2,PetscScalar z2)
{
	return (x1+x2-x0-x0)*(x1-x2) + (y_1+y2-y_0-y_0)*(y_1-y2) + (z1+z2-z0-z0)*(z1-z2);
}

/* Claim cells for particle p_i in the list */
void AVD3dClaimCells(AVD3d A,const PetscInt p_i)
{
	PetscInt i,count;
	PetscScalar x0,y_0,x1,y_1,x2,y2,z0,z1,z2,dist1;
	PetscInt *temp;
	AVDChain3d bchain;
	AVDPoint3d points;
	AVDCell3d cells;
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
		cell_num0 = bchain->new_boundary_cells[i]; /* cell number we are trying to claim */

#ifdef __AVD_DEBUG_MODE
		if (cell_num0<0) {
			printf("  AVD3dClaimCells(ERROR): p_i = %lld, [%lld] \n", (LLD)p_i,(LLD)cell_num0 );
			printf("  AVD3dClaimCells(ERROR):   point %f %f %f \n", A->points[p_i].x,A->points[p_i].y,A->points[p_i].z);
			exit(0);
		}
		if (cells[cell_num0].p == AVD_CELL_MASK) { printf("YOU SHOULD NEVER HAVE A MASKED CELL IN YOUR LIST\n"); exit(1); }
#endif

		if (cells[cell_num0].p == AVD_CELL_UNCLAIMED) { /* if cell unclaimed, then claim it */
			/* Realloc, note that we need one space more than the number of points to terminate the list */
			if( count == bchain->new_claimed_cells_malloced-1  ){
				temp = (PetscInt*) realloc( bchain->new_claimed_cells, (size_t)(bchain->new_claimed_cells_malloced + buffer)*sizeof(PetscInt) );
				bchain->new_claimed_cells = temp;
				bchain->new_claimed_cells_malloced += buffer;

				temp = (PetscInt*) realloc( bchain->new_boundary_cells, (size_t)(bchain->new_boundary_cells_malloced + buffer)*sizeof(PetscInt) );
				bchain->new_boundary_cells = temp;
				bchain->new_boundary_cells_malloced += buffer;
			}
			bchain->new_claimed_cells[count] = cell_num0;
			bchain->num_claimed++;
			count++;
			cells[cell_num0].p = p_i; /* mark cell as owned by particle p_i */
		} else if (cells[cell_num0].p != p_i) {
			/* perform distance test between points to determine ownership */
			x2 = points[p_i].x;
			y2 = points[p_i].y;
			z2 = points[p_i].z;

			x1 = points[cells[cell_num0].p].x;
			y_1 = points[cells[cell_num0].p].y;
			z1 = points[cells[cell_num0].p].z;

			/* cell centroid */
			x0 = cells[cell_num0].i*dx + (A->x0 - dx + 0.5*dx);
			y_0 = cells[cell_num0].j*dy + (A->y0 - dy + 0.5*dy);
			z0 = cells[cell_num0].k*dz + (A->z0 - dz + 0.5*dz);

			dist1 = AVD3dDistanceTest(x0,y_0,z0,x1,y_1,z1,x2,y2,z2);
			if (dist1 > 0.0) {
				bchain->new_claimed_cells[count] = cell_num0;
				bchain->num_claimed++;
				count++;
				cells[cell_num0].p = p_i; /* mark cell as owned by particle p_i */
			}
		}
		bchain->new_claimed_cells[count] = -1; /* mark end of list */
	}
}

void AVD3dUpdateChain(AVD3d A,const PetscInt p_i)
{
	PetscInt i,k;
	PetscInt count;
	PetscInt cell_num0,cell_num1,cell_num[6];
	AVDChain3d bchain;
	AVDCell3d cells,cell0;
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


		/* boundary protection */
		for (k=0; k<6; k++) {
			if (cells[cell_num[k]].p == AVD_CELL_MASK) {
				cell_num[k] = -2;
			}
		}

		for (k=0; k<6; k++) {
			cell_num1 = cell_num[k];
			/*
			 if cell does not already belong to the particle and hasn't been
			 marked as being done then add it to new boundary array and mark it as done
			 */
			if (cell_num1 != -2) {
				if ( (cells[cell_num1].p != p_i) && (cells[cell_num1].done != AVD_TRUE) ) {
					/* Realloc, note that we need one space more than the number of points to terminate the list */
					if (count == bchain->new_boundary_cells_malloced-1 ) {
						temp = (PetscInt*)realloc( bchain->new_claimed_cells, (size_t)(bchain->new_claimed_cells_malloced + buffer)*sizeof(PetscInt) );
						bchain->new_claimed_cells = temp;
						bchain->new_claimed_cells_malloced += buffer;

						temp = (PetscInt*)realloc( bchain->new_boundary_cells, (size_t)(bchain->new_boundary_cells_malloced + buffer)*sizeof(PetscInt) );
						bchain->new_boundary_cells = temp;
						bchain->new_boundary_cells_malloced += buffer;
					}
#ifdef __AVD_DEBUG_MODE
					if (cell_num1<0) {
						printf("  AVD3dUpdateChain(ERROR): INSERTING negative cell index \n");
						printf("  AVD3dUpdateChain(ERROR):   k=%lld :: cell0 i,j,k = %lld,%lld,%lld neighbourid [%lld]\n", (LLD)k,(LLD)(cell0->i), (LLD)(cell0->j), (LLD)(cell0->k), (LLD)cell_num1 );
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

	/* reset the processed flags */
	for (i=0; i<count; i++){
		cells[ bchain->new_boundary_cells[i] ].done = AVD_FALSE;
	}
}

/* I would like to be able to use this function with different numbers of particles */
void AVD3dInit(AVD3d A,const PetscInt npoints,AVDPoint3d points)
{
	PetscInt p,i,j,k;
	PetscInt mx,my,mz,ind;

	//printf("AVD3dInit: \n");
	if (npoints != A->npoints) {
		printf("AVD3dInit: npoints != A->npoints\n");
		exit(0);
	}

	mx = A->mx_mesh;
	my = A->my_mesh;
	mz = A->mz_mesh;

	for (p=0; p<npoints; p++){
		/* check if point outside the domain */
		if (points[p].x < A->x0) {  printf("AVD3dInit(ERROR): xp(%1.6f) < x0(%1.6f) \n", points[p].x,A->x0); exit(1); }
		if (points[p].y < A->y0) {  printf("AVD3dInit(ERROR): yp(%1.6f) < y0(%1.6f) \n", points[p].y,A->y0); exit(1); }
		if (points[p].z < A->z0) {  printf("AVD3dInit(ERROR): zp(%1.6f) < z0(%1.6f) \n", points[p].z,A->z0); exit(1); }

		if (points[p].x > A->x1) {  printf("AVD3dInit(ERROR): xp(%1.6f) > x1(%1.6f) \n", points[p].x,A->x1); exit(1); }
		if (points[p].y > A->y1) {  printf("AVD3dInit(ERROR): yp(%1.6f) > y1(%1.6f) \n", points[p].y,A->y1); exit(1); }
		if (points[p].z > A->z1) {  printf("AVD3dInit(ERROR): zp(%1.6f) > z1(%1.6f) \n", points[p].z,A->z1); exit(1); }


		i = (PetscInt)((points[p].x - (A->x0 - A->dx))/A->dx);
		j = (PetscInt)((points[p].y - (A->y0 - A->dy))/A->dy);
		k = (PetscInt)((points[p].z - (A->z0 - A->dz))/A->dz);

		/* If a particle is exactly on the border then make sure it is in a valid cell inside the element */
		if (i == mx) { i--; }
		if (j == my) { j--; }
		if (k == mz) { k--; }

		if (i==0) { printf("AVD3dInit(ERROR): i==0: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }
		if (j==0) { printf("AVD3dInit(ERROR): j==0: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }
		if (k==0) { printf("AVD3dInit(ERROR): k==0: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }

		if (i==A->mx_mesh-1) { printf("AVD3dInit(ERROR): i==mx: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }
		if (j==A->my_mesh-1) { printf("AVD3dInit(ERROR): j==my: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }
		if (k==A->mz_mesh-1) { printf("AVD3dInit(ERROR): k==mz: %lf %lf %lf \n", points[p].x, points[p].y, points[p].z ); exit(1); }


		ind = i+j*mx+k*mx*my;
		if (A->cells[ind].p == AVD_CELL_MASK) {
			printf("AVD3dInit(ERROR): Inserting cells into boundary cells - this is not permitted\n");
			exit(1);
		}

		A->cells[i+j*mx+k*mx*my].p = p; /* particle index */
		A->chains[p].num_claimed = 1; /* number of claimed cells, currently just the one the point initially resides within */
		A->chains[p].length = 0;
		A->chains[p].total_claimed = 1; /* total of claimed cells */
		A->chains[p].done = AVD_FALSE;
		A->chains[p].index = i+j*mx+k*mx*my; /* ith particle is in cell i +j*mx + k*mx*my */
		A->chains[p].new_claimed_cells[0] = i+j*mx+k*mx*my;  /* ith particle claimed cell it resides within, i.e. cell_index i+j*mx+k*mx*my */
		A->chains[p].new_claimed_cells[1] = -1; /* mark end of claimed_cells list with -1 */

		AVD3dUpdateChain(A,p);
	}

}

void AVD3dReportMemory(AVD3d A)
{
	PetscInt c,ncells,npoints;
	PetscScalar mem_cells, mem_points, mem_chain;

	ncells = A->mx_mesh * A->my_mesh * A->mz_mesh;
	mem_cells = (PetscScalar)(sizeof(struct _p_AVDCell3d) * (size_t)ncells );

	npoints = A->npoints;
	mem_points = (PetscScalar)(sizeof(struct _p_AVDPoint3d) * (size_t)npoints );

	mem_chain = 0.0;
	for (c=0; c<npoints; c++) {
		mem_chain = mem_chain + sizeof(PetscInt) * 7.0 + sizeof(char) * 1.0 + sizeof(char*) * 2.0;
		mem_chain = mem_chain + (PetscScalar)( sizeof(PetscInt) * (size_t)A->chains[c].new_boundary_cells_malloced );
		mem_chain = mem_chain + (PetscScalar)( sizeof(PetscInt) * (size_t)A->chains[c].new_claimed_cells_malloced );
	}

	printf("AVD3d Memory usage: \n");
	printf("  points:  %1.4f (MB) \n", mem_points*1.0e-6);
	printf("  cells:   %1.4f (MB) \n", mem_cells*1.0e-6);
	printf("  chains:  %1.4f (MB) \n", mem_chain*1.0e-6);
	printf("  [total]: %1.4f (MB) \n", (mem_chain+mem_cells+mem_points)*1.0e-6);

}

PetscErrorCode AVD3dReportPMemory(AVD3d A)
{
	PetscInt c,ncells,npoints;
	PetscScalar mem_cells, mem_points, mem_chain;
	PetscScalar gmem_cells, gmem_points, gmem_chain;
	PetscScalar _gmem_cells[2], _gmem_points[2], _gmem_chain[2];
	PetscErrorCode ierr;

	ncells = A->mx_mesh * A->my_mesh * A->mz_mesh;
	mem_cells = (PetscScalar)(sizeof(struct _p_AVDCell3d) * (size_t)ncells );

	npoints = A->npoints;
	mem_points = (PetscScalar)(sizeof(struct _p_AVDPoint3d) * (size_t)npoints );

	mem_chain = 0.0;
	for (c=0; c<npoints; c++) {
		mem_chain = mem_chain + sizeof(PetscInt) * 7.0 + sizeof(char) * 1.0 + sizeof(char*) * 2.0;
		mem_chain = mem_chain + (PetscScalar)( sizeof(PetscInt) * (size_t)A->chains[c].new_boundary_cells_malloced );
		mem_chain = mem_chain + (PetscScalar)( sizeof(PetscInt) * (size_t)A->chains[c].new_claimed_cells_malloced );
	}

	ierr = MPI_Allreduce(&mem_chain,&gmem_chain,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_cells,&gmem_cells,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_points,&gmem_points,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);

	ierr = MPI_Allreduce(&mem_chain,&_gmem_chain[0],1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_cells,&_gmem_cells[0],1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_points,&_gmem_points[0],1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);

	ierr = MPI_Allreduce(&mem_chain,&_gmem_chain[1],1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_cells,&_gmem_cells[1],1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&mem_points,&_gmem_points[1],1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);


	PetscPrintf(PETSC_COMM_WORLD,"#AVD3d PMemory usage: \n");
	PetscPrintf(PETSC_COMM_WORLD,"#  points:  %1.4f [%1.4f/%1.4f] (MB) \n", gmem_points*1.0e-6,_gmem_points[0]*1.0e-6,_gmem_points[1]*1.0e-6);
	PetscPrintf(PETSC_COMM_WORLD,"#  cells:   %1.4f [%1.4f/%1.4f] (MB) \n", gmem_cells*1.0e-6,_gmem_cells[0]*1.0e-6,_gmem_cells[1]*1.0e-6);
	PetscPrintf(PETSC_COMM_WORLD,"#  chains:  %1.4f [%1.4f/%1.4f] (MB) \n", gmem_chain*1.0e-6,_gmem_chain[0]*1.0e-6,_gmem_chain[1]*1.0e-6);
	PetscPrintf(PETSC_COMM_WORLD,"#  [total]: %1.4f [%1.4f/%1.4f] (MB) \n",
							(gmem_chain+gmem_cells+gmem_points)*1.0e-6,
							(_gmem_chain[0]+_gmem_cells[0]+_gmem_points[0])*1.0e-6,
							(_gmem_chain[1]+_gmem_cells[1]+_gmem_points[1])*1.0e-6 );
	PetscFunctionReturn(0);
}

/* this is not a parallel viewer, just a subdomain viewer */
void AVDPhaseViewer_native_AppendedVTR(AVD3d A,const char name[], const char DirectoryName[], float scaling)
{
	PetscMPIInt rank;
	FILE*	fp;
	PetscInt  i,j,k;
	PetscInt offset,L;
	char *fname;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	asprintf(&fname,"./%s/%s-p%1.6d.vtr",DirectoryName, name,rank);
	if ((fp = fopen ( fname, "w")) == NULL)  {
		exit(1);
	}

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
  fprintf(fp, "  <RectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\" >\n", 0LL,(LLD)(A->mx), 0LL,(LLD)(A->my), 0LL,(LLD)(A->mz));
	fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" >\n", 0LL,(LLD)(A->mx), 0LL,(LLD)(A->my), 0LL,(LLD)(A->mz));

	offset = 0;

	fprintf(fp, "    <Coordinates>\n");
	/* X */
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->mx+1);

	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->my+1);

	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->mz+1);

	fprintf(fp, "    </Coordinates>\n");

	fprintf(fp, "    <CellData>\n");

	// Distinguish 32-bit and 64-bit integers

#if defined(PETSC_USE_64BIT_INDICES)
	// pid
	fprintf(fp, "      <DataArray type=\"Int64\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + sizeof(PetscInt) + sizeof(PetscInt)*(A->mx * A->my * A->mz);
	// phase
	fprintf(fp, "      <DataArray type=\"Int64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + sizeof(PetscInt) + sizeof(PetscInt)*(A->mx * A->my * A->mz);
#else
	// pid
	fprintf(fp, "      <DataArray type=\"Int32\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
	// phase
	fprintf(fp, "      <DataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
#endif


	fprintf(fp, "    </CellData>\n");

	fprintf(fp, "    <PointData>\n");
	fprintf(fp, "    </PointData>\n");

	fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");


	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// X
	L = (PetscInt)sizeof(float)*(A->mx+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->mx+1; i++ ) {
		float val = (float)(A->x0 + (i)*A->dx)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// Y
	L = (PetscInt)sizeof(float)*(A->my+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->my+1; i++ ) {
		float val = (float)(A->y0 + (i)*A->dy)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// Z
	L = (PetscInt)sizeof(float)*(A->mz+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->mz+1; i++ ) {
		float val = (float)(A->z0 + (i)*A->dz)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// pid
	L = (PetscInt)sizeof(PetscInt)*(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++) {
				PetscInt ii = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				PetscInt val = A->cells[ii].p;
				fwrite(&val,sizeof(PetscInt),1,fp);
			}
		}
	}

	// phase
	L = (PetscInt)sizeof(PetscInt)*(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++) {
				PetscInt ii = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				PetscInt phase;

				phase = A->points[ A->cells[ii].p ].phase;
				fwrite(&phase,sizeof(PetscInt),1,fp);
			}
		}
	}
	fprintf(fp,"\n  </AppendedData>\n");


	fprintf(fp, "</VTKFile>\n");

	fclose( fp );
}

void AVDPhaseViewer_AppendedVTR(AVD3d A,const char name[], const char DirectoryName[], float scaling)
{
	PetscMPIInt rank;
	FILE*	fp;
	char *fname;
	PetscInt  i,j,k;
	PetscInt offset,L;
	PetscInt r2d,pi,pj,pk;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	asprintf(&fname,"./%s/%s-p%1.6d.vtr",DirectoryName, name,rank);
	if ((fp = fopen ( fname, "w")) == NULL)  {
		exit(1);
	}
	free(fname);

	pk = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj = r2d/(A->M);
	pi = r2d - pj*A->M;


	fprintf(fp, "<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
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
	/* X */
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->mx+1);

	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->my+1);

	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(float)*(A->mz+1);

	fprintf(fp, "    </Coordinates>\n");

	fprintf(fp, "    <CellData>\n");

	// Distinguish 32-bit and 64-bit integers

#if defined(PETSC_USE_64BIT_INDICES)
	// pid
	fprintf(fp, "      <DataArray type=\"Int64\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
	// phase
	fprintf(fp, "      <DataArray type=\"Int64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
#else
	// pid
	fprintf(fp, "      <DataArray type=\"Int32\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
	// phase
	fprintf(fp, "      <DataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset = offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*(A->mx * A->my * A->mz);
#endif

	fprintf(fp, "    </CellData>\n");

	fprintf(fp, "    <PointData>\n");
	fprintf(fp, "    </PointData>\n");

	fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");


	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// X
	L = (PetscInt)sizeof(float)*(A->mx+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->mx+1; i++ ) {
		float val = (float)(A->x0 + (i)*A->dx)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// Y
	L = (PetscInt)sizeof(float)*(A->my+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->my+1; i++ ) {
		float val = (float)(A->y0 + (i)*A->dy)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// Z
	L = (PetscInt)sizeof(float)*(A->mz+1);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for( i=0; i<A->mz+1; i++ ) {
		float val = (float)(A->z0 + (i)*A->dz)*scaling;
		fwrite(&val,sizeof(float),1,fp);
	}

	// pid
	L = (PetscInt)sizeof(PetscInt)*(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++) {
				PetscInt ii = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				PetscInt val = A->cells[ii].p;
				fwrite(&val,sizeof(PetscInt),1,fp);
			}
		}
	}

	// phase
	L = (PetscInt)sizeof(PetscInt)*(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(PetscInt), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++) {
				PetscInt ii = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				PetscInt phase;

				phase = A->points[ A->cells[ii].p ].phase;
				fwrite(&phase,sizeof(PetscInt),1,fp);
			}
		}
	}
	fprintf(fp,"\n  </AppendedData>\n");


	fprintf(fp, "</VTKFile>\n");

	fclose( fp );
}

void AVDPhaseViewer_PVTR(AVD3d A,const char name[], const char DirectoryName[])
{
	PetscMPIInt nproc,rank;
	FILE*	fp;
	char *fname;
	PetscInt r2d,p,pi,pj,pk;

	MPI_Comm_size(PETSC_COMM_WORLD,&nproc);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	asprintf(&fname,"./%s/%s.pvtr",DirectoryName, name);
	if ((fp = fopen ( fname, "w")) == NULL)  {
		exit(1);
	}
	free(fname);

	pk = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj = r2d/(A->M);
	pi = r2d - pj*A->M;

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf(fp, "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp, "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
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

	// Distinguish 32-bit and 64-bit integers

#if defined(PETSC_USE_64BIT_INDICES)
	fprintf(fp, "      <PDataArray type=\"Int64\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Int64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" />\n");
#else
	fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"pid\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" />\n");
#endif

	fprintf(fp, "    </PCellData>\n");

	fprintf(fp, "    <PPointData>\n");
	fprintf(fp, "    </PPointData>\n");

	for (p=0; p<nproc; p++) {
		pk = p/(A->M*A->N);
		r2d = p - pk*(A->M*A->N);
		pj = r2d/(A->M);
		pi = r2d - pj*A->M;

		asprintf(&fname,"%s-p%1.6lld.vtr",name,(LLD)p);
		fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" Source=\"%s\" />\n",
				(LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
				(LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
				(LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]),
						fname );

		free(fname);
	}

	fprintf(fp, "  </PRectilinearGrid>\n");


	fprintf(fp, "</VTKFile>\n");

	fclose( fp );
}

PetscErrorCode LaMEMGetLocalBoundingBox(LaMEMVelPressureDA C,DM DA_Processors,DM da,PetscScalar lmin[],PetscScalar lmax[])
{
	DM cda;
	PetscInt i,j,k;
	PetscInt xsp,ysp,zsp,xmp,ymp,zmp;
	PetscInt iel_x,iel_y,iel_z;
	DMDACoor3d ***coords, coord_elem[MAX_nnel];
	Vec			 		gc;
	DAVPElementType element_type;
	PetscInt n,nnel;
  PetscScalar mmin[3]={PETSC_MAX_REAL,PETSC_MAX_REAL,PETSC_MAX_REAL},mmax[3]={PETSC_MIN_REAL,PETSC_MIN_REAL,PETSC_MIN_REAL};
	PetscErrorCode ierr;

	element_type = C->type;
	nnel    = C->nnel;

	ierr = DMDAGetCorners(DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);

	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){

				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables
				if( (element_type==DAVP_Q1P0) | (element_type==DAVP_Q1Q1) ) {
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}

				GetElementCoords(coord_elem, coords, i,j,k, 1);

				for (n=0; n<nnel; n++) {

					mmin[0] = PetscMin(mmin[0],PetscRealPart(coord_elem[n].x));CHKERRQ(ierr);
					mmax[0] = PetscMax(mmax[0],PetscRealPart(coord_elem[n].x));CHKERRQ(ierr);

					mmin[1] = PetscMin(mmin[1],PetscRealPart(coord_elem[n].y));CHKERRQ(ierr);
					mmax[1] = PetscMax(mmax[1],PetscRealPart(coord_elem[n].y));CHKERRQ(ierr);

					mmin[2] = PetscMin(mmin[2],PetscRealPart(coord_elem[n].z));CHKERRQ(ierr);
					mmax[2] = PetscMax(mmax[2],PetscRealPart(coord_elem[n].z));CHKERRQ(ierr);
				}

			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);

  if (lmin) {ierr = PetscMemcpy(lmin,mmin,3*sizeof(PetscScalar));CHKERRQ(ierr);}
  if (lmax) {ierr = PetscMemcpy(lmax,mmax,3*sizeof(PetscScalar));CHKERRQ(ierr);}

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AVDPhaseViewerLoadLaMEMPoints"
PetscErrorCode AVDPhaseViewerLoadLaMEMPoints(PetscScalar min[],PetscScalar max[],PetscInt nLpoints,Particles Lpoints[],PetscInt *_nump,AVDPoint3d *_points)
{
	PetscInt nump;
	AVDPoint3d points;
	PetscInt c,p;
	PetscInt gc,gnLpoints;
	PetscScalar x0,x1,y_0,y_1,z0,z1;
	PetscErrorCode ierr;

	x0 = min[0];  x1 = max[0];
	y_0 = min[1];  y_1 = max[1];
	z0 = min[2];  z1 = max[2];

	nump = nLpoints; /* might be an upper bound as some points might be outside domain */

	AVDPoint3dCreate(nump,&points);

	c = 0;
	for (p=0; p<nump; p++) {
		PetscScalar xp,yp,zp;

		xp = Lpoints[p].x;
		yp = Lpoints[p].y;
		zp = Lpoints[p].z;

		/* check if point outside the domain */
		if (xp < x0) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): xp(%1.6f) < x0(%1.6f) \n", xp,x0);*/ continue; }
		if (yp < y_0) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): yp(%1.6f) < y0(%1.6f) \n", yp,y0);*/ continue; }
		if (zp < z0) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): zp(%1.6f) < z0(%1.6f) \n", zp,z0);*/ continue; }

		if (xp > x1) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): xp(%1.6f) > x1(%1.6f) \n", xp,x1);*/ continue; }
		if (yp > y_1) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): yp(%1.6f) > y1(%1.6f) \n", yp,y1);*/ continue; }
		if (zp > z1) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): zp(%1.6f) > z1(%1.6f) \n", zp,z1);*/ continue; }

		if (Lpoints[p].phase < 0) {  /*printf("AVDPhaseViewerLoadLaMEMPoints(WARN): phase(%lld) < 0 \n",Lpoints[p].phase );*/ continue; }


		points[c].x = xp;
		points[c].y = yp;
		points[c].z = zp;

		points[c].phase = (PetscInt)Lpoints[p].phase;
		c++;
	}
	nump = c;

	ierr = MPI_Allreduce(&c,&gc,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&nLpoints,&gnLpoints,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"AVD3d: Loaded %lld of %lld LaMEM points\n",(LLD)gc,(LLD)gnLpoints );

	*_nump   = nump;
	*_points = points;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AVDPhaseViewerExecute"
PetscErrorCode AVDPhaseViewerExecute(LaMEMVelPressureDA C,DM dau,DM dap,PetscInt nLpoints,Particles Lpoints[],
		const char NAME[], const char DirectoryName[], float scaling)
{
	AVD3d A;
	AVDPoint3d points;
	PetscInt fac;
	PetscInt npoints;
	PetscInt i,count,claimed;
	PetscScalar mmin[3],mmax[3];
	PetscInt mx,my,mz,M,N,P;
	PetscLogDouble t0,t1;
	PetscErrorCode ierr;

	ierr = PetscTime(&t0);CHKERRQ(ierr);
	ierr = LaMEMGetLocalBoundingBox(C,dap,dau,mmin,mmax);CHKERRQ(ierr);

	/* get info from petsc objects */
	ierr = DMDAGetInfo(dap,0, 0,0,0, &M,&N,&P, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dap,0,0,0,&mx,&my,&mz);CHKERRQ(ierr);
	ierr = AVDPhaseViewerLoadLaMEMPoints(mmin,mmax,nLpoints,Lpoints,&npoints,&points);CHKERRQ(ierr);

	/* build up avd data structure */
	fac = 2;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-AVDPhaseViewer_factor",&fac,PETSC_NULL);CHKERRQ(ierr);
	AVD3dCreate(fac*mx,fac*my,fac*mz,1,&A);
	AVD3dSetParallelExtent(A,M,N,P);

	AVD3dSetDomainSize(A,mmin[0],mmax[0],mmin[1],mmax[1],mmin[2],mmax[2]);




	/* set points */
	AVD3dSetPoints(A,npoints,points);

	/* reset cells */
	AVDCell3dReset(A);

	AVD3dInit(A,npoints,points);
	ierr = PetscTime(&t1);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"# AVD3d: Initialisation %1.4e (sec)\n", t1-t0);


	ierr = PetscTime(&t0);CHKERRQ(ierr);
	count = npoints;
	claimed = 1;
	while (claimed != 0){
		claimed = 0 ;
		for (i=0; i<count; i++){
			AVD3dClaimCells(A,i);
			claimed += A->chains[i].num_claimed;
			AVD3dUpdateChain(A,i);
		}
	}
	ierr = PetscTime(&t1);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"# AVD3d: Voronoi generator %1.4e (sec)\n", t1-t0);


	AVD3dReportPMemory(A);
	/*AVDPhaseViewer_native_AppendedVTR(A,"avd_phase_native.vtr",DirectoryName, scaling);*/
	ierr = PetscTime(&t0);CHKERRQ(ierr);
	AVDPhaseViewer_AppendedVTR(A,NAME, DirectoryName, scaling);
	AVDPhaseViewer_PVTR(A,NAME, DirectoryName);
	ierr = PetscTime(&t1);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"# AVD3d: VTR writer %1.4e (sec)\n", t1-t0);

	AVDPoint3dDestroy(&A->points);
	AVD3dDestroy(&A);
	PetscFunctionReturn(0);
}

