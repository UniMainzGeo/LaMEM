#include "LaMEMParaViewOutMark_C.h"
#include "../LaMEM.h"
#include "../paraViewOutMark.h"
#include "../paraViewOutBin.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../advect.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMParaViewOutMark_PVMarkCreate(void *pvmark, void *fb) {
    return PVMarkCreate((PVMark*)pvmark, (FB*)fb);
}

PetscErrorCode LaMEMParaViewOutMark_PVMarkWriteTimeStep(void *pvmark, const char *dirName, PetscScalar ttime) {
    return PVMarkWriteTimeStep((PVMark*)pvmark, dirName, ttime);
}

PetscErrorCode LaMEMParaViewOutMark_PVMarkWriteVTU(void *pvmark, const char *dirName) {
    return PVMarkWriteVTU((PVMark*)pvmark, dirName);
}

PetscErrorCode LaMEMParaViewOutMark_PVMarkWritePVTU(void *pvmark, const char *dirName) {
    return PVMarkWritePVTU((PVMark*)pvmark, dirName);
}

}