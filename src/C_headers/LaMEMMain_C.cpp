#include "LaMEMMain_C.h"
#include "../LaMEM.h"
#include "../scaling.h"
#include "../objFunct.h"
#include "../parsing.h"
#include "../adjoint.h"
#include "../phase.h"

extern "C" {

// Thin wrappers - just cast and call the actual existing functions
PetscErrorCode LaMEMMain_LibMain(void *fb, PetscLogStage *stages) {
    return LaMEMLibMain((FB*)fb, stages);
}

}