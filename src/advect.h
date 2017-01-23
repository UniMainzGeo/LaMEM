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
 **    filename:   advect.h
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
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __advect_h__
#define __advect_h__
//---------------------------------------------------------------------------

#define _cap_overhead_ 1.3
#define _max_nmark_ 5

//---------------------------------------------------------------------------

// marker initialization type enumeration
typedef enum
{
	_GEOM_,    // read geometric primitives from input file
	_FILES_,   // read coordinates, phase and temperature from files in parallel
	_POLYGONS_ // read polygons from file redundantly

} SetupType;

//---------------------------------------------------------------------------

// marker-to-edge / edge-to-marker interpolation cases
typedef enum
{
	_PHASE_,     // phase ratio
	_STRESS_,    // deviatoric stress
	_APS_,       // accumulated plastic strain
	_VORTICITY_, // vorticity pseudo-vector components
	_DISP_       // displacement

} InterpCase;

//-----------------------------------------------------------------------------

typedef enum
{
	EULER,          // euler explicit in time
	RUNGE_KUTTA_2,  // runge-kutta 2nd order in space

} AdvectionType;

//-----------------------------------------------------------------------------

typedef enum
{
	STAG,      // trilinear interp from fdstag points
	MINMOD,    // minmod interp to nodes, trilinear interp to markers + correction
	STAG_P     // empirical approach (T. Gerya)

} VelInterpType;

//---------------------------------------------------------------------------

typedef struct
{
	PetscInt s[8]; // 8 corners

} NumCorner;

//---------------------------------------------------------------------------

// Set of variables that must be tracked during the advection steps.
// During the nonlinear iteration, history accumulates on the integration points.
// During the advection a cumulative update of the history variables is added to
// the individual histories of the markers locating in the integration point
// control volume. There are four types of control volumes for different stress
// components: XY-, XZ-, YZ- edges, and cell centers. Each particle should be
// mapped on all types of control volumes to update corresponding components.

//---------------------------------------------------------------------------

// Advection context
typedef struct
{
	// staggered grid
	FDSTAG   *fs;
	JacRes   *jr;
	FreeSurf *surf;

	SetupType     msetup;                // marker initialization type
	PetscInt      NumPartX;              // markers per cell in x-direction
	PetscInt      NumPartY;              //                 ... y-direction
	PetscInt      NumPartZ;              //                 ... z-direction
	PetscInt      randNoise;             // random noise flag for marker distribution
	PetscInt      bgPhase;               // background phase ID

	PetscInt      saveMark;              // flag for saving markers
	char          saveName[MAX_NAME_LEN]; // marker output file name
	char          savePath[MAX_PATH_LEN]; // marker output directory

	AdvectionType advection;             // advection scheme
	PetscInt      newAdv;                // new advection flag (temporary)
	VelInterpType velinterp;             // velocity interpolation scheme
	PetscScalar   A;                     // FDSTAG velocity interpolation parameter

	PetscInt      markContr;             // flag to activate marker control
	PetscInt      newMarkContr;          // new marker control flag (temporary)
	PetscInt      nmin, nmax;            // minimum and maximum number of markers per cell
	PetscInt      avdx, avdy, avdz;      // AVD cells refinement factors

	//====================
	// RUN TIME PARAMETERS
	//====================
	PetscInt    AirPhase;   // air phase number
	PetscScalar Ttop;       // top surface temperature
	PetscInt    cinj, cdel; // injected & deleted marker counters

	//=============
	// COMMUNICATOR
	//=============
	MPI_Comm  icomm;   // distinct communicator for communicating markers
	PetscInt  nproc;   // total number of processors
	PetscInt  iproc;   // processor rank

	//========
	// STORAGE
	//========
	PetscInt  nummark;    // local number of markers
	PetscInt  markcap;    // capacity of marker storage
	Marker   *markers;    // storage for local markers

	//========================
	// MARKER-CELL INTERACTION
	//========================
	PetscInt *cellnum;    // host cells local number for each marker
	PetscInt *markind;    // id (position) of markers clustered for every cell
	PetscInt *markstart;  // start id in markind for every cell

	//=========
	// EXCHANGE
	//=========
	Marker   *sendbuf; // send buffer
	Marker   *recvbuf; // receive buffer

	PetscInt  nsend;                // total number of markers to be sent (local)
	PetscInt  nsendm[_num_neighb_]; // number of markers to be sent to each process
	PetscInt  ptsend[_num_neighb_]; // send buffer pointers

	PetscInt  nrecv;                // total number of markers to be received (local)
	PetscInt  nrecvm[_num_neighb_]; // number of markers to be received from each process
	PetscInt  ptrecv[_num_neighb_]; // receive buffer pointers

	PetscInt  ndel; // number of markers to be deleted from storage
	PetscInt *idel; // indices of markers to be deleted

} AdvCtx;

//---------------------------------------------------------------------------
// create advection object
PetscErrorCode ADVCreate(AdvCtx *actx, FB *fb);

// read advection object from restart database
PetscErrorCode ADVReadRestart(AdvCtx *actx, FILE *fp);

// read advection object from restart database
PetscErrorCode ADVWriteRestart(AdvCtx *actx, FILE *fp);

// create communicator and
PetscErrorCode ADVCreateData(AdvCtx *actx);

// destroy advection context
PetscErrorCode ADVDestroy(AdvCtx *actx);

// (re)allocate marker storage
PetscErrorCode ADVReAllocStorage(AdvCtx *actx, PetscInt capacity);

// perform advection step
PetscErrorCode ADVAdvect(AdvCtx *actx);

// remap markers onto the grid
PetscErrorCode ADVRemap(AdvCtx *actx);

// exchange markers between the processors resulting from the position change
PetscErrorCode ADVExchange(AdvCtx *actx);

// project history INCREMENTS from grid to markers
PetscErrorCode ADVProjHistGridToMark(AdvCtx *actx);

// interpolate field history increments to markers
PetscErrorCode ADVInterpFieldToMark(AdvCtx *actx, InterpCase icase);

// update marker positions from current velocities & time step
PetscErrorCode ADVAdvectMark(AdvCtx *actx);

// count number of markers to be sent to each neighbor domain
PetscErrorCode ADVMapMarkToDomains(AdvCtx *actx);

// communicate number of markers with neighbor processes
PetscErrorCode ADVExchangeNumMark(AdvCtx *actx);

// create send and receive buffers for asynchronous MPI communication
PetscErrorCode ADVCreateMPIBuff(AdvCtx *actx);

// communicate markers with neighbor processes
PetscErrorCode ADVExchangeMark(AdvCtx *actx);

// store received markers, collect garbage
PetscErrorCode ADVCollectGarbage(AdvCtx *actx);

// free communication buffer
PetscErrorCode ADVDestroyMPIBuff(AdvCtx *actx);

// find host cells for local markers
PetscErrorCode ADVMapMarkToCells(AdvCtx *actx);

// creates arrays to optimize marker-cell interaction
PetscErrorCode ADVUpdateMarkCell(AdvCtx *actx);

// project history fields from markers to grid
PetscErrorCode ADVProjHistMarkToGrid(AdvCtx *actx);

// marker-to-cell projection
PetscErrorCode ADVInterpMarkToCell(AdvCtx *actx);

// marker-to-edge projection
PetscErrorCode ADVInterpMarkToEdge(AdvCtx *actx, PetscInt iphase, InterpCase icase);

// inject or delete markers
PetscErrorCode ADVMarkControl(AdvCtx *actx);

PetscErrorCode ADVCheckCorners(AdvCtx *actx);

// delete marker outflow
PetscErrorCode ADVMarkDeleteOutflow(AdvCtx *actx);

// change marker phase when crossing free surface
PetscErrorCode ADVMarkCrossFreeSurf(AdvCtx *actx, PetscScalar tol);

// check marker phases
PetscErrorCode ADVCheckMarkPhases(AdvCtx *actx, PetscInt numPhases);

//-----------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------

// compute pointers from counts, return total count
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]);

// rewind pointers after using them as access iterators
void rewindPtr(PetscInt n, PetscInt ptr[]);

// compute phase ratio array
PetscErrorCode getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum);

// find ID of the cell containing point (call this function for local point only!)
static inline PetscInt FindPointInCell(
	PetscScalar *px, // node coordinates
	PetscInt     L,  // index of the leftmost node
	PetscInt     R,  // index of the rightmost node
	PetscScalar  x)  // point coordinate
{
	// get initial guess assuming uniform grid
	PetscInt M = L + (PetscInt)((x-px[L])/((px[R]-px[L])/(PetscScalar)(R-L)));

	if(M == R) return R-1;

	if(px[M]   <= x) L=M;
	if(px[M+1] >= x) R=M+1;

	while((R-L) > 1)
	{
		M = (L+R)/2;
		if(px[M] <= x) L=M;
		if(px[M] >= x) R=M;

	}
	return(L);
}
//-----------------------------------------------------------------------------
static inline PetscScalar InterpLin3D(
	PetscScalar ***lv,
	PetscInt    i,
	PetscInt    j,
	PetscInt    k,
	PetscInt    sx,
	PetscInt    sy,
	PetscInt    sz,
	PetscScalar xp,
	PetscScalar yp,
	PetscScalar zp,
	PetscScalar *cx,
	PetscScalar *cy,
	PetscScalar *cz)
{
	PetscScalar xb, yb, zb, xe, ye, ze, v;

	// get relative coordinates
	xe = (xp - cx[i])/(cx[i+1] - cx[i]); xb = 1.0 - xe;
	ye = (yp - cy[j])/(cy[j+1] - cy[j]); yb = 1.0 - ye;
	ze = (zp - cz[k])/(cz[k+1] - cz[k]); zb = 1.0 - ze;

	// interpolate & return result
	v =
	lv[sz+k  ][sy+j  ][sx+i  ]*xb*yb*zb +
	lv[sz+k  ][sy+j  ][sx+i+1]*xe*yb*zb +
	lv[sz+k  ][sy+j+1][sx+i  ]*xb*ye*zb +
	lv[sz+k  ][sy+j+1][sx+i+1]*xe*ye*zb +
	lv[sz+k+1][sy+j  ][sx+i  ]*xb*yb*ze +
	lv[sz+k+1][sy+j  ][sx+i+1]*xe*yb*ze +
	lv[sz+k+1][sy+j+1][sx+i  ]*xb*ye*ze +
	lv[sz+k+1][sy+j+1][sx+i+1]*xe*ye*ze;

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar InterpLin2D(
	PetscScalar ***lv,
	PetscInt    i,
	PetscInt    j,
	PetscInt    L,
	PetscInt    sx,
	PetscInt    sy,
	PetscScalar xp,
	PetscScalar yp,
	PetscScalar *cx,
	PetscScalar *cy)
{
	PetscScalar xb, yb, xe, ye, v;

	// get relative coordinates
	xe = (xp - cx[i])/(cx[i+1] - cx[i]); xb = 1.0 - xe;
	ye = (yp - cy[j])/(cy[j+1] - cy[j]); yb = 1.0 - ye;

	// interpolate & return result
	v =
	lv[L][sy+j  ][sx+i  ]*xb*yb +
	lv[L][sy+j  ][sx+i+1]*xe*yb +
	lv[L][sy+j+1][sx+i  ]*xb*ye +
	lv[L][sy+j+1][sx+i+1]*xe*ye;

	return v;
}
//-----------------------------------------------------------------------------
#endif
