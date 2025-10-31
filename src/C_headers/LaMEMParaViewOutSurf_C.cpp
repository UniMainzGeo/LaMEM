#include "LaMEMParaViewOutSurf_C.h"
#include "../LaMEM.h"
#include "../paraViewOutSurf.h"
#include "../paraViewOutBin.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../fdstag.h"
#include "../surf.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMParaViewOutSurf_PVSurfCreate(void *pvsurf, void *fb) {
    return PVSurfCreate((PVSurf*)pvsurf, (FB*)fb);
}

PetscErrorCode LaMEMParaViewOutSurf_PVSurfDestroy(void *pvsurf) {
    return PVSurfDestroy((PVSurf*)pvsurf);
}

PetscErrorCode LaMEMParaViewOutSurf_PVSurfWriteTimeStep(void *pvsurf, const char *dirName, PetscScalar ttime) {
    return PVSurfWriteTimeStep((PVSurf*)pvsurf, dirName, ttime);
}

PetscErrorCode LaMEMParaViewOutSurf_PVSurfWriteVTS(void *pvsurf, const char *dirName) {
    return PVSurfWriteVTS((PVSurf*)pvsurf, dirName);
}

PetscErrorCode LaMEMParaViewOutSurf_PVSurfWritePVTS(void *pvsurf, const char *dirName) {
    return PVSurfWritePVTS((PVSurf*)pvsurf, dirName);
}

}