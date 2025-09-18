#ifndef LAMEMNLSOLVE_C_H
#define LAMEMNLSOLVE_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual nlsolve functions from the source
PetscErrorCode LaMEMNLSolve_NLSolClear(void *nl);
PetscErrorCode LaMEMNLSolve_NLSolCreate(void *nl, void *pc, SNES *p_snes);
PetscErrorCode LaMEMNLSolve_NLSolDestroy(void *nl);
PetscErrorCode LaMEMNLSolve_DisplaySpecifiedSolverOptions(void *pc, SNES snes);

// Callback functions for SNES
PetscErrorCode LaMEMNLSolve_FormResidual(SNES snes, Vec x, Vec f, void *ctx);
PetscErrorCode LaMEMNLSolve_FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx);
PetscErrorCode LaMEMNLSolve_JacApplyMFFD(Mat A, Vec x, Vec y);
PetscErrorCode LaMEMNLSolve_SNESPrintConvergedReason(SNES snes, PetscLogDouble t_beg);
PetscErrorCode LaMEMNLSolve_SNESCoupledTest(SNES snes, PetscInt it, PetscReal xnorm, 
                                            PetscReal gnorm, PetscReal f, 
                                            SNESConvergedReason *reason, void *cctx);

#ifdef __cplusplus
}
#endif

#endif