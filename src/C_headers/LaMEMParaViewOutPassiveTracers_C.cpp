#include "LaMEMParaViewOutPassiveTracers_C.h"
#include "../LaMEM.h"
#include "../paraViewOutBin.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../advect.h"
#include "../JacRes.h"
#include "../paraViewOutPassiveTracers.h"
#include "../tools.h"
#include "../passive_tracer.h"


extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrCreate(void *pvptr, void *fb) {
    return PVPtrCreate((PVPtr*)pvptr, (FB*)fb);
}

PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWriteTimeStep(void *pvptr, const char *dirName, PetscScalar ttime) {
    return PVPtrWriteTimeStep((PVPtr*)pvptr, dirName, ttime);
}

PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWriteVTU(void *pvptr, const char *dirName) {
    return PVPtrWriteVTU((PVPtr*)pvptr, dirName);
}

PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWritePVTU(void *pvptr, const char *dirName) {
    return PVPtrWritePVTU((PVPtr*)pvptr, dirName);
}

}