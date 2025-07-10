#ifndef LAMEMMAIN_C_H
#define LAMEMMAIN_C_H

#include <petsc.h>
#include "LaMEMMain_C.h"
#include "../LaMEM.h"
#include "../scaling.h"
#include "../objFunct.h"
#include "../parsing.h"
#include "../adjoint.h"
#include "../phase.h"


#ifdef __cplusplus
extern "C" {
#endif

// Thin wrapper for LaMEMLibMain
PetscErrorCode LaMEMMain_LibMain(void *fb, PetscLogStage *stages);

#ifdef __cplusplus
}
#endif

#endif