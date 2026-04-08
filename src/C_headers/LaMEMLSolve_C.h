#ifndef LAMEMLSOLVE_C_H
#define LAMEMLSOLVE_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual lsolve functions from the source
PetscErrorCode LaMEMLSolve_PCStokesSetFromOptions(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesCreate(void **p_pc, void *pm);
PetscErrorCode LaMEMLSolve_PCStokesSetup(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesDestroy(void *pc);

// Block factorization functions
PetscErrorCode LaMEMLSolve_PCStokesBFCreate(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesBFSetFromOptions(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesBFDestroy(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesBFSetup(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesBFApply(Mat JP, Vec r, Vec x);

// Multigrid functions
PetscErrorCode LaMEMLSolve_PCStokesMGCreate(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesMGDestroy(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesMGSetup(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesMGApply(Mat JP, Vec x, Vec y);

// User-defined preconditioner functions
PetscErrorCode LaMEMLSolve_PCStokesUserCreate(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesUserAttachIS(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesUserDestroy(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesUserSetup(void *pc);
PetscErrorCode LaMEMLSolve_PCStokesUserApply(Mat JP, Vec x, Vec y);

#ifdef __cplusplus
}
#endif

#endif