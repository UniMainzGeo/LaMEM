#ifndef LAMEMCVI_C_H
#define LAMEMCVI_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual CVI functions from cvi.h
PetscErrorCode LaMEMCVI_AdvectMain(void *actx);
PetscErrorCode LaMEMCVI_InterpPT(void *actx);
PetscErrorCode LaMEMCVI_AdvectScheme(void *actx, void *vi);
PetscErrorCode LaMEMCVI_CollectIndices(void *actx, void *vi);
PetscErrorCode LaMEMCVI_Create(void *actx, void *vi);
PetscErrorCode LaMEMCVI_Destroy(void *vi);
PetscErrorCode LaMEMCVI_RungeKuttaStep(void *vi, PetscScalar dt, PetscScalar a, PetscInt type);
PetscErrorCode LaMEMCVI_InterpMain(void *vi);
PetscErrorCode LaMEMCVI_InterpSTAG(void *vi);
PetscErrorCode LaMEMCVI_InterpMINMOD(void *vi);
PetscErrorCode LaMEMCVI_Exchange(void *vi);

#ifdef __cplusplus
}
#endif

#endif