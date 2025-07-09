#ifndef LAMEMMAIN_C_H
#define LAMEMMAIN_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual LaMEM functions from the source
PetscErrorCode LaMEMMain_LibMain(void *fb, PetscLogStage *stages);
PetscErrorCode LaMEMMain_AdjointMain(void *IOparam);

#ifdef __cplusplus
}
#endif

#endif