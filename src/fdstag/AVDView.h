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

#ifndef __AVDView_h__
#define __AVDView_h__

#define AVD_TRUE  'T'
#define AVD_FALSE 'F'
#define AVD_CELL_MASK -2
#define AVD_CELL_UNCLAIMED -1

typedef struct _p_AVDCell3D *AVDCell3D;
typedef struct _p_AVDChain3D *AVDChain3D;
typedef struct _p_AVDPoint3D *AVDPoint3D;
typedef struct _p_AVD3D *AVD3D;

//---------------------------------------------------------------------------

struct _p_AVDCell3D
{
	PetscInt p;     // particle index which this cell is attributed to
	PetscInt index; // i + j.mx + k.mx.my
	PetscInt i,j,k; // cell index i,j,k
	char     done;
};

//---------------------------------------------------------------------------

void AVDCell3DCreate(const PetscInt mx, const PetscInt my, const PetscInt mz, AVDCell3D *C);

void AVDCell3DDestroy(AVDCell3D *C);

void AVDCell3DReset(AVD3D A);

//---------------------------------------------------------------------------

struct _p_AVDChain3D
{
	PetscInt  p;
	PetscInt  index;
	PetscInt  length; // current length of the boundary chain
	PetscInt  num_claimed;
	PetscInt  total_claimed;
	PetscInt  new_boundary_cells_malloced;
	PetscInt  new_claimed_cells_malloced;
	PetscInt *new_boundary_cells;
	PetscInt *new_claimed_cells;
	char      done;
};

//---------------------------------------------------------------------------

void AVDChain3DCreate(const PetscInt npoints, const PetscInt buffer, AVDChain3D *CH);

void AVDChain3DDestroy(const PetscInt npoints, AVDChain3D *CH);

//---------------------------------------------------------------------------

struct _p_AVDPoint3D
{
	PetscScalar x,y,z;
	PetscInt    phase;
};

//---------------------------------------------------------------------------

void AVDPoint3DCreate(const PetscInt npoints, AVDPoint3D *P);

void AVDPoint3DDestroy(AVDPoint3D *P);

//---------------------------------------------------------------------------

struct _p_AVD3D
{
	PetscScalar x0, x1, y0, y1, z0, z1; // size of domain
	PetscScalar dx, dy, dz;
	PetscInt    buffer;
	PetscInt    mx, my, mz; // user specified resolution
	PetscInt    mx_mesh, my_mesh, mz_mesh; // computational resolution, mx_mesh = mx + 2
	AVDCell3D   cells;
	PetscInt    npoints;
	AVDChain3D  chains;
	AVDPoint3D  points;
	PetscInt    M, N, P;
	PetscInt    gmx, gmy, gmz;
	PetscInt   *ownership_ranges_i; // for pvtr output
	PetscInt   *ownership_ranges_j; // for pvtr output
	PetscInt   *ownership_ranges_k; // for pvtr output
};

//---------------------------------------------------------------------------

void AVD3DCreate(
	const PetscInt mx,
	const PetscInt my,
	const PetscInt mz,
	const PetscInt buffer,
	AVD3D          *A);

void AVD3DDestroy(AVD3D *A);

PetscErrorCode AVD3DSetParallelExtent(AVD3D A, PetscInt M, PetscInt N, PetscInt P);

void AVD3DSetDomainSize(AVD3D A,
	const PetscScalar x0,
	const PetscScalar x1,
	const PetscScalar y0,
	const PetscScalar y1,
	const PetscScalar z0,
	const PetscScalar z1);

void AVD3DSetPoints(AVD3D A,const PetscInt npoints, AVDPoint3D points);

PetscErrorCode AVD3DInit(AVD3D A, const PetscInt npoints, AVDPoint3D points);

void AVD3DClaimCells(AVD3D A,const PetscInt p_i);

void AVD3DUpdateChain(AVD3D A,const PetscInt p_i);

PetscErrorCode AVD3DReportMemory(AVD3D A);

//---------------------------------------------------------------------------

typedef struct
{

	long int  offset;  // pvd file offset
	PetscInt  outpvd;  // pvd file output flag
	PetscBool avdAct;  // output activation flag

} AVDView;

//---------------------------------------------------------------------------

PetscErrorCode AVDViewCreate(AVDView *avdout);

PetscErrorCode AVDViewWriteStep(AVDView *avdout, AdvCtx *actx, const char DirectoryName[], PetscScalar ttime, PetscInt tindx);

PetscErrorCode AVDViewWriteVTR(AVD3D A, const char name[], const char DirectoryName[], PetscScalar scaling);

PetscErrorCode AVDViewWritePVTR(AVD3D A, const char name[], const char DirectoryName[]);

//---------------------------------------------------------------------------

PetscErrorCode AVDViewLoadPoints(AdvCtx *actx, PetscInt *_nump, AVDPoint3D *_points);

PetscErrorCode AVDViewExecute(AdvCtx *actx, const char NAME[], const char DirectoryName[]);

//---------------------------------------------------------------------------
#endif
