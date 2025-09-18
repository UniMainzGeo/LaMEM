#include "LaMEMNLSolve_C.h"
#include "../LaMEM.h"
#include "../matrix.h"
#include "../fdstag.h"
#include "../tssolve.h"
#include "../multigrid.h"
#include "../lsolve.h"
#include "../nlsolve.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMNLSolve_NLSolClear(void *nl) {
    return NLSolClear((NLSol*)nl);
}

PetscErrorCode LaMEMNLSolve_NLSolCreate(void *nl, void *pc, SNES *p_snes) {
    return NLSolCreate((NLSol*)nl, (PCStokes)pc, p_snes);
}

PetscErrorCode LaMEMNLSolve_NLSolDestroy(void *nl) {
    return NLSolDestroy((NLSol*)nl);
}

PetscErrorCode LaMEMNLSolve_DisplaySpecifiedSolverOptions(void *pc, SNES snes) {
    return DisplaySpecifiedSolverOptions((PCStokes)pc, snes);
}

// Callback functions for SNES
PetscErrorCode LaMEMNLSolve_FormResidual(SNES snes, Vec x, Vec f, void *ctx) {
    return FormResidual(snes, x, f, ctx);
}

PetscErrorCode LaMEMNLSolve_FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx) {
    return FormJacobian(snes, x, Amat, Pmat, ctx);
}

PetscErrorCode LaMEMNLSolve_JacApplyMFFD(Mat A, Vec x, Vec y) {
    return JacApplyMFFD(A, x, y);
}

PetscErrorCode LaMEMNLSolve_SNESPrintConvergedReason(SNES snes, PetscLogDouble t_beg) {
    return SNESPrintConvergedReason(snes, t_beg);
}

PetscErrorCode LaMEMNLSolve_SNESCoupledTest(SNES snes, PetscInt it, PetscReal xnorm, 
                                            PetscReal gnorm, PetscReal f, 
                                            SNESConvergedReason *reason, void *cctx) {
    return SNESCoupledTest(snes, it, xnorm, gnorm, f, reason, cctx);
}

}