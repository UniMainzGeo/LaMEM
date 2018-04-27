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
 **    filename:   marker.h
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
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// input polygon data
typedef struct
{
	PetscInt     dir; // normal vector of polygon plane
	PetscInt   ax[2]; // axis that span the polygon plane
	PetscInt   phase; // phase that the polygon defines
	PetscInt    type; // type can be of additive or assigning nature
	PetscInt     num; // number of polygon slices defining the volume
	PetscInt     len; // number of nodes of polygon
	PetscInt    idxs; // index of first polygon slice
	PetscInt    gidx; // global plane index (consistent with planes of markers)
	PetscInt    lidx; // local plane index (consistent with planes of markers)
	PetscInt       n; // number of polygon nodes
	PetscInt   nmark; // number of markers in current volume
	PetscScalar   *X; // coordinates of the polygon x1,y1,x2,y2,...xn,yn

} Polygon2D;

//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, UserCtx *user);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, UserCtx *user);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx, UserCtx *user);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFilePolygons (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBands        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitPipes        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitGeoth		   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFault        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitRozhko       (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitCon          (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitPrefrac      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDomes        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitRotation     (AdvCtx *actx, UserCtx *user);

//---------------------------------------------------------------------------

// service functions

PetscErrorCode ADVMarkSetTempFromFile  (AdvCtx *actx, UserCtx *user);

void ADVMarkSecIdx(AdvCtx *actx, UserCtx *user, PetscInt dir, PetscInt Nslice, PetscInt *idx);

//---------------------------------------------------------------------------

// definitions

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

//---------------------------------------------------------------------------
#endif
