/*
 * passive_tracer.h
 *
 *
 *  Created on: Jul 28, 2020
 *      Author: piccolo
 */
#ifndef passive_tracer_h_
#define passive_tracer_h_
//---------------------------------------------------------------------------

#include "Tensor.h" // required for Marker declaration
#include "LaMEM.h"
#include "AVD.h"
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "fdstag.h"
#include "bc.h"
#include "tools.h"
#include "phase_transition.h"
#include "phase.h"
#include "constEq.h"
#include "parsing.h"
#include "objFunct.h"
//---------------------------------------------------------------------------

struct FB;
struct FDSTAG;
struct JacRes;
struct FreeSurf;
struct DBMat;
/*
 * The passive tracer are placed at the cell center. They are globally identified by the ID of the cell
 * The global ID number is given by the following formulation ID = i + Npx*(ip+j*nx+jp*ny)+Npx*Npy*npx*npy*(k+kp*nz)
 * where ip,jp,kp are the current index of the processor; Npx,Npy,Npz are the total number of processor along x,y,z direction
 * nx,ny,nz are the number of nodes along x,y,z direction contained in each processor.
 * Each time the marker are advected, and interpolated revelant information: e.g. pressure and temperature
 */

struct Passive_Tracers
{
	PetscInt    nummark ;
	Vec         ID;    // global identification number
	Vec    phase; // phase identifier
	Vec    x;      // global coordinates
	Vec    y;      //
	Vec    z;      //
	Vec    p;     // pressure
	Vec    T;     // temperature
};

PetscErrorCode ADVPtrReAllocStorage(AdvCtx *actx);

PetscErrorCode ADVPassiveTracerInit(AdvCtx *actx);

PetscErrorCode ADVPtrInitCoord(AdvCtx *actx);

PetscErrorCode ADV_Assign_Phase(AdvCtx *actx);

PetscErrorCode ADVAdvectPassiveTracer(AdvCtx *actx);

PetscErrorCode ADVMarkCrossFreeSurfPassive_Tracers(AdvCtx *actx);

PetscErrorCode ADVPtrDestroy(AdvCtx *actx);

PetscErrorCode Passive_tracers_save(AdvCtx *actx);


#endif /* PASSIVE_TRACER_H */





