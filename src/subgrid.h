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
 **    filename:   subgrid.h
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
//.................   SUBGRID MARKER RESAMPLING ROUTINES   ..................
//---------------------------------------------------------------------------
#ifndef __subgrid_h__
#define __subgrid_h__
//---------------------------------------------------------------------------

struct AdvCtx;
struct Marker;

//---------------------------------------------------------------------------

// resample markers
PetscErrorCode ADVMarkSubGrid(AdvCtx *actx);

// change marker phase when crossing free surface
PetscErrorCode ADVMarkCrossFreeSurf(AdvCtx *actx);

// compute reference sedimentation phases
PetscErrorCode ADVGetSedPhase(AdvCtx *actx, Vec vphase);

// rearrange storage after marker resampling
PetscErrorCode ADVCollectGarbageVec(AdvCtx *actx, vector <Marker> &recvbuf, vector <PetscInt> &idel);

#define MAP_SUBCELL(i, x, s, h, n) \
{ i = (PetscInt)PetscFloorReal(((x) - (s))/(h)); if(i > n - 1) { i = n - 1; } if(i < 0) { i = 0; } }

#define COORD_SUBCELL(x, i, s, h) (x) = (s) + (i)*(h) + (h)/2.0

#define EDIST(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));

//---------------------------------------------------------------------------
#endif
