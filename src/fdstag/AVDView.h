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
 **    filename:   AVDView.h
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

struct _p_AVDPoint3D
{
	PetscScalar x,y,z;
	PetscInt    phase;
};

//---------------------------------------------------------------------------

void AVDPoint3DCreate(const PetscInt npoints, AVDPoint3D *P);

void AVDPoint3DDestroy(AVDPoint3D *P);

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

PetscErrorCode AVDViewCreate(AVD3D *A, AdvCtx *actx, PetscInt refine);

void AVD3DDestroy(AVD3D *A);

void AVD3DAllocate(
	const PetscInt mx,
	const PetscInt my,
	const PetscInt mz,
	const PetscInt buffer,
	const PetscInt npoints,
	AVD3D          *A);

PetscErrorCode AVD3DSetParallelExtent(AVD3D A, PetscInt M, PetscInt N, PetscInt P);

void AVD3DSetDomainSize(AVD3D A,
	const PetscScalar x0,
	const PetscScalar x1,
	const PetscScalar y0,
	const PetscScalar y1,
	const PetscScalar z0,
	const PetscScalar z1);

PetscErrorCode AVD3DLoadPoints(AVD3D A, AdvCtx *actx);

void AVD3DResetCells(AVD3D A);

PetscErrorCode AVD3DInit(AVD3D A);

void AVD3DClaimCells(AVD3D A,const PetscInt p_i);

void AVD3DUpdateChain(AVD3D A,const PetscInt p_i);

PetscErrorCode AVD3DReportMemory(AVD3D A);

//---------------------------------------------------------------------------

typedef struct
{
	AdvCtx    *actx;    // advection context
	char      *outfile; // output file name
	long int  offset;   // pvd file offset
	PetscInt  outavd;   // AVD output flag
	PetscInt  refine;   // Voronoi Diagram refinement factor
	PetscInt  outpvd;   // pvd file output flag

} PVAVD;

//---------------------------------------------------------------------------

PetscErrorCode PVAVDCreate(PVAVD *pvavd, AdvCtx *actx, const char *filename);

PetscErrorCode PVAVDDestroy(PVAVD *pvavd);

PetscErrorCode PVAVDReadFromOptions(PVAVD *pvavd);

PetscErrorCode PVAVDWriteTimeStep(PVAVD *pvavd, const char *dirName, PetscScalar ttime, PetscInt tindx);

PetscErrorCode PVAVDWritePVTR(PVAVD *pvavd, AVD3D A, const char *dirName);

PetscErrorCode PVAVDWriteVTR(PVAVD *pvavd, AVD3D A, const char *dirName);

//---------------------------------------------------------------------------
#endif
