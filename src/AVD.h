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
 **    filename:   AVD.h
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
//..........   Routines based on Approximate Voronoi Diagram (AVD)  .........
//---------------------------------------------------------------------------
#ifndef __AVD_h__
#define __AVD_h__
//---------------------------------------------------------------------------

struct Marker;
struct AdvCtx;

//---------------------------------------------------------------------------

#define AVD_CELL_MASK      -2
#define AVD_CELL_UNCLAIMED -1

struct AVDCell
{
	PetscInt    ind;                       // single cell index
	PetscInt    i,j,k;                     // i,j,k  cell index
	PetscScalar x[3];                      // coordinates of center
	PetscInt    p;                         // marker index
	PetscBool   done;                      // flag
	PetscInt    col;                       // colour for half-centroid

} ;

struct AVDChain
{
	PetscInt    p;                         // marker index
	PetscInt    ind;                       // index
	PetscInt    length;                    // current length of chain
	PetscInt    nclaimed;                  // claimed cells in the current cycle
	PetscInt    tclaimed;                  // total no of cells claimed
	PetscInt    ibound;                    // no. of new boundary cells
	PetscInt    iclaim;                    // no. of new claimed cells
	PetscInt    *bound;                    // new boundary cells
	PetscInt    *claim;                    // new claimed cells
	PetscBool   done;                      // flag
	PetscInt    gind;                      // marker index in actx->markers
	PetscScalar xc[3];                     // centroid coordinates
	PetscScalar xh[3];                     // half-axis of the centroid
	PetscInt    axis;                      // dominant axis of centroid cell

};

struct AVD
{
	PetscInt    mmin, mmax;                // limit number of markers
	PetscScalar xs[3],xe[3];               // coordinate limits of the Voronoi diagram
	PetscScalar dx,dy,dz;                  // spacing
	PetscInt    nx,ny,nz;                  // grid cells
	PetscInt    buffer;                    // buffer
	AVDCell     *cell;                     // voronoi grid (size of nx*ny*nz)
	AVDChain    *chain;                    // voronoi chain for every point (size of npoints)
	Marker      *points;                   // points that we want to compute voronoi diagram (size of npoints)
	PetscInt    npoints;                   // no. markers

} ;

struct MarkerVolume
{
	PetscInt    *cellnum;                  // host cells local number for each marker
	PetscInt    *markind;                  // id (position) of markers clustered for every cell
	PetscInt    *markstart;                // start id in markind for every cell
	PetscInt     ncells;                   // no of total cells
	PetscScalar *xcoord, *ycoord, *zcoord; // coordinate arrays
	PetscInt     M, N, P;                  // number of cells in each direction

} ;

enum VolumeCase
{
	_CELL_, // center cell
	_XYED_, // xy-edge
	_XZED_, // xz-edge
	_YZED_  // yz-edge

};

//---------------------------------------------------------------------------

// basic AVD routines
PetscErrorCode AVDCreate     (AVD *A);
PetscErrorCode AVDDestroy    (AVD *A);
PetscErrorCode AVDCellInit   (AVD *A);
PetscErrorCode AVDClaimCells (AVD *A, const PetscInt ip);
PetscErrorCode AVDUpdateChain(AVD *A, const PetscInt ip);
PetscErrorCode AVDReAlloc    (AVDChain *chain,PetscInt buffer);

// routines for old marker control
PetscErrorCode AVDLoadPoints            (AdvCtx *actx, AVD *A, PetscInt ind);
PetscErrorCode AVDInjectDeletePoints    (AdvCtx *actx, AVD *A, PetscInt cellID);
PetscErrorCode AVDExecuteMarkerInjection(AdvCtx *actx, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind);

// new marker control (for every control volume)
PetscErrorCode AVDMarkerControl  (AdvCtx *actx);
PetscErrorCode AVDMarkerControlMV(AdvCtx *actx, VolumeCase vtype);
PetscErrorCode AVDCheckCellsMV   (AdvCtx *actx, MarkerVolume *mv, PetscInt dir);
PetscErrorCode AVDMapMarkersMV   (AdvCtx *actx, MarkerVolume *mv, PetscInt dir);
PetscErrorCode AVDCreateMV       (AdvCtx *actx, MarkerVolume *mv, PetscInt dir);
PetscErrorCode AVDAlgorithmMV    (AdvCtx *actx, MarkerVolume *mv, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind, PetscInt nmin);
PetscErrorCode AVDLoadPointsMV   (AdvCtx *actx, MarkerVolume *mv, AVD *A, PetscInt ind);
PetscErrorCode AVDInjectPointsMV (AdvCtx *actx, AVD *A);
PetscErrorCode AVDDeletePointsMV (AdvCtx *actx, AVD *A);
PetscErrorCode AVDDestroyMV      (MarkerVolume *mv);

//---------------------------------------------------------------------------
static inline PetscScalar AVDDistanceTest(PetscScalar x0[3],PetscScalar x1[3],PetscScalar x2[3])
{
	return (x1[0]+x2[0]-2*x0[0])*(x1[0]-x2[0]) + (x1[1]+x2[1]-2*x0[1])*(x1[1]-x2[1]) + (x1[2]+x2[2]-2*x0[2])*(x1[2]-x2[2]);
}
//---------------------------------------------------------------------------
#endif
