/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
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

// clone closest marker & put it in the center of an empty subcell
PetscErrorCode ADVMarkClone(
	AdvCtx          *actx,
	PetscInt         icell,
	PetscInt         isubcell,
	PetscScalar     *s,
	PetscScalar     *h,
	vector <spair>  &dist,
	vector <Marker> &iclone);

// merge markers in a densely populated subcell
PetscErrorCode ADVMarkCheckMerge(
	AdvCtx            *actx,
	PetscInt           ib,
	PetscInt           ie,
	PetscInt          &nmerge,
	vector <Marker>   &mark,
	vector <ipair>    &cell,
	vector <Marker>   &iclone,
	vector <PetscInt> &imerge);

// recursively find and merge closest markers until required number is reached
PetscErrorCode ADVMarkMerge(
	vector <Marker> &mark,
	PetscInt         nmark,
	PetscInt         npmax,
	PetscInt        &sz);

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
