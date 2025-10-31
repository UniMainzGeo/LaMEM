#include "LaMEMPassiveTracer_C.h"
#include "../LaMEM.h"
#include "../AVD.h"
#include "../advect.h"
#include "../scaling.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../tools.h"
#include "../phase_transition.h"
#include "../phase.h"
#include "../constEq.h"
#include "../parsing.h"
#include "../objFunct.h"
#include "../surf.h"
#include "../subgrid.h"
#include "../tssolve.h"

extern "C" {

// Main passive tracer functions
PetscErrorCode LaMEMPassiveTracer_ADVPtrPassive_Tracer_create(void *actx, void *fb) {
    return ADVPtrPassive_Tracer_create((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMPassiveTracer_ADVPtrReCreateStorage(void *actx) {
    return ADVPtrReCreateStorage((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_ADVPassiveTracerInit(void *actx) {
    return ADVPassiveTracerInit((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_ADVPtrInitCoord(void *actx) {
    return ADVPtrInitCoord((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_ADV_Assign_Phase(void *actx) {
    return ADV_Assign_Phase((AdvCtx*)actx);
}

// Advection functions
PetscErrorCode LaMEMPassiveTracer_ADVAdvectPassiveTracer(void *actx) {
    return ADVAdvectPassiveTracer((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_ADVMarkCrossFreeSurfPassive_Tracers(void *actx) {
    return ADVMarkCrossFreeSurfPassive_Tracers((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_Check_advection_condition(void *actx, PetscInt jj, PetscInt ID, 
                                                           PetscScalar xp, PetscScalar yp, PetscScalar zp, 
                                                           PetscScalar P, PetscScalar T, PetscScalar mf) {
    return Check_advection_condition((AdvCtx*)actx, jj, ID, xp, yp, zp, P, T, mf);
}

// Cleanup and I/O functions
PetscErrorCode LaMEMPassiveTracer_ADVPtrDestroy(void *actx) {
    return ADVPtrDestroy((AdvCtx*)actx);
}

PetscErrorCode LaMEMPassiveTracer_Passive_Tracer_WriteRestart(void *actx, FILE *fp) {
    return Passive_Tracer_WriteRestart((AdvCtx*)actx, fp);
}

PetscErrorCode LaMEMPassiveTracer_ReadPassive_Tracers(void *actx, FILE *fp) {
    return ReadPassive_Tracers((AdvCtx*)actx, fp);
}

// Utility functions
PetscErrorCode LaMEMPassiveTracer_Sync_Vector(Vec x, void *actx, PetscInt nummark) {
    return Sync_Vector(x, (AdvCtx*)actx, nummark);
}

}