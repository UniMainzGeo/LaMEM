#include "LaMEMCVI_C.h"
#include "../LaMEM.h"
#include "../cvi.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../advect.h"
#include "../surf.h"
#include "../bc.h"
#include "../tssolve.h"
#include "../tools.h"

extern "C" {

// Thin wrappers around actual CVI functions
PetscErrorCode LaMEMCVI_AdvectMain(void *actx) {
    return ADVelAdvectMain((AdvCtx*)actx);
}

PetscErrorCode LaMEMCVI_InterpPT(void *actx) {
    return ADVelInterpPT((AdvCtx*)actx);
}

PetscErrorCode LaMEMCVI_AdvectScheme(void *actx, void *vi) {
    return ADVelAdvectScheme((AdvCtx*)actx, (AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_CollectIndices(void *actx, void *vi) {
    return ADVelCollectIndices((AdvCtx*)actx, (AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_Create(void *actx, void *vi) {
    return ADVelCreate((AdvCtx*)actx, (AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_Destroy(void *vi) {
    return ADVelDestroy((AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_RungeKuttaStep(void *vi, PetscScalar dt, PetscScalar a, PetscInt type) {
    return ADVelRungeKuttaStep((AdvVelCtx*)vi, dt, a, type);
}

PetscErrorCode LaMEMCVI_InterpMain(void *vi) {
    return ADVelInterpMain((AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_InterpSTAG(void *vi) {
    return ADVelInterpSTAG((AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_InterpMINMOD(void *vi) {
    return ADVelInterpMINMOD((AdvVelCtx*)vi);
}

PetscErrorCode LaMEMCVI_Exchange(void *vi) {
    return ADVelExchange((AdvVelCtx*)vi);
}

}