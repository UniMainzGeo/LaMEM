#include "LaMEMLSolve_C.h"
#include "../LaMEM.h"
#include "../fdstag.h"
#include "../matrix.h"
#include "../multigrid.h"
#include "../lsolve.h"
#include "../JacRes.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMLSolve_PCStokesSetFromOptions(void *pc) {
    return PCStokesSetFromOptions((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesCreate(void **p_pc, void *pm) {
    return PCStokesCreate((PCStokes*)p_pc, (PMat)pm);
}

PetscErrorCode LaMEMLSolve_PCStokesSetup(void *pc) {
    return PCStokesSetup((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesDestroy(void *pc) {
    return PCStokesDestroy((PCStokes)pc);
}

// Block factorization functions
PetscErrorCode LaMEMLSolve_PCStokesBFCreate(void *pc) {
    return PCStokesBFCreate((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesBFSetFromOptions(void *pc) {
    return PCStokesBFSetFromOptions((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesBFDestroy(void *pc) {
    return PCStokesBFDestroy((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesBFSetup(void *pc) {
    return PCStokesBFSetup((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesBFApply(Mat JP, Vec r, Vec x) {
    return PCStokesBFApply(JP, r, x);
}

// Multigrid functions
PetscErrorCode LaMEMLSolve_PCStokesMGCreate(void *pc) {
    return PCStokesMGCreate((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesMGDestroy(void *pc) {
    return PCStokesMGDestroy((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesMGSetup(void *pc) {
    return PCStokesMGSetup((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesMGApply(Mat JP, Vec x, Vec y) {
    return PCStokesMGApply(JP, x, y);
}

// User-defined preconditioner functions
PetscErrorCode LaMEMLSolve_PCStokesUserCreate(void *pc) {
    return PCStokesUserCreate((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesUserAttachIS(void *pc) {
    return PCStokesUserAttachIS((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesUserDestroy(void *pc) {
    return PCStokesUserDestroy((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesUserSetup(void *pc) {
    return PCStokesUserSetup((PCStokes)pc);
}

PetscErrorCode LaMEMLSolve_PCStokesUserApply(Mat JP, Vec x, Vec y) {
    return PCStokesUserApply(JP, x, y);
}

}