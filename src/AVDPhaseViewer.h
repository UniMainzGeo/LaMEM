/*
 *  AVDPhaseViewer.c
 *
 *
 *  Created by Dave A. May on 6/21/11.
 *  Copyright 2011 Geophysical Fluid Dynamics. All rights reserved.
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
 *
 *
 */

#ifndef __AVD3d_h__
#define __AVD3d_h__

#define AVD_TRUE  'T'
#define AVD_FALSE 'F'
#define AVD_CELL_MASK -2
#define AVD_CELL_UNCLAIMED -1

typedef struct _p_AVDCell3d *AVDCell3d;
typedef struct _p_AVDChain3d *AVDChain3d;
typedef struct _p_AVDPoint3d *AVDPoint3d;
typedef struct _p_AVD3d *AVD3d;

struct _p_AVDCell3d {
	PetscInt p;     /* particle index which this cell is attributed to */
	PetscInt index; /* i + j.mx + k.mx.my */
	PetscInt i,j,k; /* cell index i,j,k*/
	char done;
};


struct _p_AVDChain3d {
	PetscInt p;
	PetscInt index;
	PetscInt length; /* current length of the boundary chain */
	PetscInt num_claimed;
	PetscInt total_claimed;
	PetscInt new_boundary_cells_malloced;
	PetscInt new_claimed_cells_malloced;
	PetscInt *new_boundary_cells;
	PetscInt *new_claimed_cells;
	char done;
};

struct _p_AVDPoint3d {
	PetscScalar x,y,z;
	PetscInt phase;
};

struct _p_AVD3d {
	PetscScalar x0,x1,y0,y1,z0,z1; /* size of domain */
	PetscScalar dx,dy,dz;
	PetscInt buffer;
	PetscInt mx,my,mz; /* user specified resolution */
	PetscInt mx_mesh,my_mesh,mz_mesh; /* computational resolution, mx_mesh = mx + 1 */
	AVDCell3d cells;
	PetscInt npoints;
	AVDChain3d chains;
	AVDPoint3d points;
	PetscInt M,N,P;
	PetscInt gmx,gmy,gmz;
	PetscInt *ownership_ranges_i; /* for pvtr output */
	PetscInt *ownership_ranges_j; /* for pvtr output */
	PetscInt *ownership_ranges_k; /* for pvtr output */
};


/* public functions */
PetscErrorCode AVDPhaseViewerExecute(LaMEMVelPressureDA C,DM dau,DM dap,PetscInt nLpoints,Particles Lpoints[],
		const char NAME[], const char DirectoryName[], float scaling);

void _AVDCell3dCreate(const PetscInt mx,const PetscInt my, const PetscInt mz,AVDCell3d *C);
void AVDCell3dDestroy(AVDCell3d *C);
void AVDCell3dReset(AVD3d A);
void _AVDChain3dCreate(const PetscInt npoints, const PetscInt buffer,AVDChain3d *CH);
void AVDChain3dDestroy(const PetscInt npoints,AVDChain3d *CH);
void AVDPoint3dCreate(const PetscInt npoints, AVDPoint3d *P);
void AVDPoint3dDestroy(AVDPoint3d *P);
void AVD3dCreate(const PetscInt mx,const PetscInt my, const PetscInt mz,const PetscInt buffer,AVD3d *A);
PetscErrorCode AVD3dSetParallelExtent(AVD3d A,PetscInt M,PetscInt N,PetscInt P);
void AVD3dDestroy(AVD3d *A);
void AVD3dSetDomainSize(AVD3d A,const PetscScalar x0,const PetscScalar x1,const PetscScalar y0,const PetscScalar y1,const PetscScalar z0,const PetscScalar z1);
void AVD3dSetPoints(AVD3d A,const PetscInt npoints,AVDPoint3d points);
void AVD3dClaimCells(AVD3d A,const PetscInt p_i);
void AVDPhaseViewer_PVTR(AVD3d A,const char name[], const char DirectoryName[]);
PetscErrorCode LaMEMGetLocalBoundingBox(LaMEMVelPressureDA C,DM DA_Processors,DM da,PetscScalar lmin[],PetscScalar lmax[]);
PetscErrorCode AVDPhaseViewerLoadLaMEMPoints(PetscScalar min[],PetscScalar max[],PetscInt nLpoints,Particles Lpoints[],PetscInt *_nump,AVDPoint3d *_points);
void AVDPhaseViewer_AppendedVTR(AVD3d A,const char name[], const char DirectoryName[], float scaling);
void AVDPhaseViewer_native_AppendedVTR(AVD3d A,const char name[], const char DirectoryName[], float scaling);
void AVD3dUpdateChain(AVD3d A,const PetscInt p_i);
void AVD3dInit(AVD3d A,const PetscInt npoints,AVDPoint3d points);
void AVD3dReportMemory(AVD3d A);
PetscErrorCode AVD3dReportPMemory(AVD3d A);


#endif
