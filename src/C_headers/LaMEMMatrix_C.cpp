#include "LaMEMMatrix_C.h"
#include "../LaMEM.h"
#include "../matrix.h"
#include "../tssolve.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// Matrix creation and assembly functions
PetscErrorCode LaMEMMatrix_MatAIJCreate(PetscInt m, PetscInt n, PetscInt d_nz, const PetscInt d_nnz[], 
                                        PetscInt o_nz, const PetscInt o_nnz[], Mat *P) {
    return MatAIJCreate(m, n, d_nz, d_nnz, o_nz, o_nnz, P);
}

PetscErrorCode LaMEMMatrix_MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P) {
    return MatAIJCreateDiag(m, istart, P);
}

PetscErrorCode LaMEMMatrix_MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag) {
    return MatAIJAssemble(P, numRows, rows, diag);
}

PetscErrorCode LaMEMMatrix_MatAIJSetNullSpace(Mat P, void *dof) {
    return MatAIJSetNullSpace(P, (DOFIndex*)dof);
}

// PMat context management
PetscErrorCode LaMEMMatrix_PMatCreate(void **p_pm, void *jr) {
    return PMatCreate((PMat*)p_pm, (JacRes*)jr);
}

PetscErrorCode LaMEMMatrix_PMatSetFromOptions(void *pm) {
    return PMatSetFromOptions((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatAssemble(void *pm) {
    return PMatAssemble((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatDestroy(void *pm) {
    return PMatDestroy((PMat)pm);
}

// Monolithic matrix functions
PetscErrorCode LaMEMMatrix_PMatMonoCreate(void *pm) {
    return PMatMonoCreate((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatMonoAssemble(void *pm) {
    return PMatMonoAssemble((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatMonoPicard(Mat J, Vec x, Vec r) {
    return PMatMonoPicard(J, x, r);
}

PetscErrorCode LaMEMMatrix_PMatMonoDestroy(void *pm) {
    return PMatMonoDestroy((PMat)pm);
}

// Block matrix functions
PetscErrorCode LaMEMMatrix_PMatBlockCreate(void *pm) {
    return PMatBlockCreate((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatBlockAssemble(void *pm) {
    return PMatBlockAssemble((PMat)pm);
}

PetscErrorCode LaMEMMatrix_PMatBlockPicardClean(Mat J, Vec x, Vec r) {
    return PMatBlockPicardClean(J, x, r);
}

PetscErrorCode LaMEMMatrix_PMatBlockPicardSchur(Mat J, Vec x, Vec r) {
    return PMatBlockPicardSchur(J, x, r);
}

PetscErrorCode LaMEMMatrix_PMatBlockDestroy(void *pm) {
    return PMatBlockDestroy((PMat)pm);
}

// Service functions
void LaMEMMatrix_GetStiffMatDevProj(PetscScalar eta, PetscScalar diag, PetscScalar *v, PetscScalar *cf,
                                    PetscScalar dx, PetscScalar dy, PetscScalar dz,
                                    PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                    PetscScalar bdx, PetscScalar bdy, PetscScalar bdz) {
    getStiffMatDevProj(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);
}

void LaMEMMatrix_GetStiffMatClean(PetscScalar eta, PetscScalar diag, PetscScalar *v, PetscScalar *cf,
                                  PetscScalar dx, PetscScalar dy, PetscScalar dz,
                                  PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                  PetscScalar bdx, PetscScalar bdy, PetscScalar bdz) {
    getStiffMatClean(eta, diag, v, cf, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);
}

void LaMEMMatrix_AddDensGradStabil(PetscScalar fssa, PetscScalar *v, PetscScalar rho, PetscScalar dt, 
                                   PetscScalar *grav, PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                   PetscScalar bdx, PetscScalar bdy, PetscScalar bdz) {
    addDensGradStabil(fssa, v, rho, dt, grav, fdx, fdy, fdz, bdx, bdy, bdz);
}

void LaMEMMatrix_GetVelSchur(PetscScalar v[], PetscScalar d[], PetscScalar g[]) {
    getVelSchur(v, d, g);
}

void LaMEMMatrix_GetSubMat(PetscScalar v[], PetscScalar a[], PetscScalar d[], PetscScalar g[]) {
    getSubMat(v, a, d, g);
}

void LaMEMMatrix_GetTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[]) {
    getTwoPointConstr(n, idx, pdofidx, cf);
}

void LaMEMMatrix_ConstrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[]) {
    constrLocalMat(n, pdofidx, cf, v);
}

// Vector scatter functions
PetscErrorCode LaMEMMatrix_VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, PetscInt mode) {
    ScatterMode scatter_mode;
    switch(mode) {
        case 0: scatter_mode = SCATTER_FORWARD; break;
        case 1: scatter_mode = SCATTER_REVERSE; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return VecScatterBlockToMonolithic(f, g, b, scatter_mode);
}

}