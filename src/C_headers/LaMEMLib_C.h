#ifndef LAMEMLIB_C_H
#define LAMEMLIB_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Simple wrappers for existing LaMEMLib functions
PetscErrorCode LaMEMLib_Main(void *param, PetscLogStage stages[4]);
PetscErrorCode LaMEMLib_Create(void *lm, void *param);
PetscErrorCode LaMEMLib_Solve(void *lm, void *param, PetscLogStage stages[4]);
PetscErrorCode LaMEMLib_SaveOutput(void *lm);
PetscErrorCode LaMEMLib_SaveRestart(void *lm);
PetscErrorCode LaMEMLib_LoadRestart(void *lm);
PetscErrorCode LaMEMLib_SaveGrid(void *lm);
PetscErrorCode LaMEMLib_DryRun(void *lm);
PetscErrorCode LaMEMLib_Destroy(void *lm);

#ifdef __cplusplus
}
#endif

#endif