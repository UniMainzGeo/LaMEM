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

#include "Tensor.h" // required for Marker declaration

//---------------------------------------------------------------------------

struct FB;
struct FDSTAG;
struct JacRes;
struct FreeSurf;
struct DBMat;

//---------------------------------------------------------------------------
//............   Material marker (history variables advection)   ............
//---------------------------------------------------------------------------

struct Marker
{
	PetscInt    phase; // phase identifier
	PetscScalar X[3];  // global coordinates
	PetscScalar p;     // pressure
	PetscScalar T;     // temperature
	PetscScalar APS;   // accumulated plastic strain
	PetscScalar ATS;   // accumulated total strain
	Tensor2RS   S;     // deviatoric stress
	PetscScalar U[3];  // displacement

	// WARNING! after adding new field modify marker merge routine (below)
};

// merge two markers and average history and position C = (A + B)/2
PetscErrorCode MarkerMerge(Marker &A, Marker &B, Marker &C);

//---------------------------------------------------------------------------

// marker initialization type enumeration
enum SetupType
{
	_GEOM_,    // read geometric primitives from input file
	_FILES_,   // read coordinates, phase and temperature from files in parallel
	_POLYGONS_ // read polygons from file redundantly

};

//---------------------------------------------------------------------------

// marker-to-edge / edge-to-marker interpolation cases
enum InterpCase
{
	_PHASE_,     // phase ratio
	_STRESS_,    // deviatoric stress
	_APS_,       // accumulated plastic strain
	_ATS_,       // accumulated total strain
	_VORTICITY_, // vorticity pseudo-vector components
	_DISP_       // displacement

};

//-----------------------------------------------------------------------------

enum AdvectionType
{
	ADV_NONE,       // no advection (grid-based solution)
	BASIC_EULER,    // basic Euler implementation (STAG interpolation only)
	EULER,          // Euler explicit in time
	RUNGE_KUTTA_2,  // Runge-Kutta 2nd order in space
};

//-----------------------------------------------------------------------------

enum VelInterpType
{
	STAG,      // trilinear interpolation from FDSTAG points
	MINMOD,    // MINMOD interpolation to nodes, trilinear interpolation to markers + correction
	STAG_P     // empirical approach (T. Gerya)

};

//-----------------------------------------------------------------------------

enum MarkCtrlType
{
	CTRL_NONE,  // no marker control
	CTRL_BASIC, // AVD for cells + corner insertion
	CTRL_AVD,   // pure AVD for all control volumes
	CTRL_SUB    // simple marker control method based on higher resolution grid (subgrid)

};

//---------------------------------------------------------------------------

struct NumCorner
{
	PetscInt s[8]; // 8 corners

};

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
struct AdvCtx
{
	// staggered grid
	FDSTAG   *fs;
	JacRes   *jr;
	FreeSurf *surf;
	DBMat    *dbm;

	SetupType     msetup;              // marker initialization type
	PetscInt      NumPartX;            // markers per cell in x-direction
	PetscInt      NumPartY;            //                 ... y-direction
	PetscInt      NumPartZ;            //                 ... z-direction
	PetscInt      randNoise;           // random noise flag for marker distribution
	PetscInt      randNoiseGP;         // random noise flag, subsequently applied to geometric primitives
	PetscInt      bgPhase;             // background phase ID

	PetscInt      saveMark;            // flag for saving markers
	char          saveFile[_str_len_]; // marker output file name

	AdvectionType advect;              // advection scheme
	VelInterpType interp;              // velocity interpolation scheme
	PetscScalar   A;                   // FDSTAG velocity interpolation parameter

	MarkCtrlType  mctrl;               // marker control type

	//====================
	// RUN TIME PARAMETERS
	//====================
	PetscInt    cinj, cdel;       // injected & deleted marker counters
	PetscInt    nmin, nmax;       // minimum and maximum number of markers per cell
	PetscInt    avdx, avdy, avdz; // AVD cells refinement factors
	PetscInt    npmax;            // maximum number of same phase markers per subcell

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

};

//---------------------------------------------------------------------------
// create advection object
PetscErrorCode ADVCreate(AdvCtx *actx, FB *fb);

PetscErrorCode ADVSetType(AdvCtx *actx, FB *fb);

// read advection object from restart database
PetscErrorCode ADVReadRestart(AdvCtx *actx, FILE *fp);

// read advection object from restart database
PetscErrorCode ADVWriteRestart(AdvCtx *actx, FILE *fp);

// create communicator and separator
PetscErrorCode ADVCreateData(AdvCtx *actx);

// destroy advection context
PetscErrorCode ADVDestroy(AdvCtx *actx);

// set background phase in all control volumes
PetscErrorCode ADVSetBGPhase(AdvCtx *actx);

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

// store host cell ID for every marker & list of marker IDs in every cell
PetscErrorCode ADVMapMarkToCells(AdvCtx *actx);

// project history fields from markers to grid
PetscErrorCode ADVProjHistMarkToGrid(AdvCtx *actx);

// marker-to-cell projection
PetscErrorCode ADVInterpMarkToCell(AdvCtx *actx);

// marker-to-edge projection
PetscErrorCode ADVInterpMarkToEdge(AdvCtx *actx, PetscInt iphase, InterpCase icase);

// inject or delete markers
PetscErrorCode ADVMarkControl(AdvCtx *actx);

PetscErrorCode ADVCheckCorners(AdvCtx *actx);

// check marker phases
PetscErrorCode ADVCheckMarkPhases(AdvCtx *actx);

// update history variables without advection
PetscErrorCode ADVUpdateHistADVNone(AdvCtx *actx);

// get maximum inverse time step (CFL)
PetscErrorCode ADVSelectTimeStep(AdvCtx *actx, PetscInt *restart);

//---------------------------------------------------------------------------
#endif
