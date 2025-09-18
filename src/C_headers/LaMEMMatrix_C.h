#ifndef LAMEMMATRIX_C_H
#define LAMEMMATRIX_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Matrix creation and assembly functions
PetscErrorCode LaMEMMatrix_MatAIJCreate(PetscInt m, PetscInt n, PetscInt d_nz, const PetscInt d_nnz[], 
                                        PetscInt o_nz, const PetscInt o_nnz[], Mat *P);
PetscErrorCode LaMEMMatrix_MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P);
PetscErrorCode LaMEMMatrix_MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag);
PetscErrorCode LaMEMMatrix_MatAIJSetNullSpace(Mat P, void *dof);

// PMat context management
PetscErrorCode LaMEMMatrix_PMatCreate(void **p_pm, void *jr);
PetscErrorCode LaMEMMatrix_PMatSetFromOptions(void *pm);
PetscErrorCode LaMEMMatrix_PMatAssemble(void *pm);
PetscErrorCode LaMEMMatrix_PMatDestroy(void *pm);

// Monolithic matrix functions
PetscErrorCode LaMEMMatrix_PMatMonoCreate(void *pm);
PetscErrorCode LaMEMMatrix_PMatMonoAssemble(void *pm);
PetscErrorCode LaMEMMatrix_PMatMonoPicard(Mat J, Vec x, Vec r);
PetscErrorCode LaMEMMatrix_PMatMonoDestroy(void *pm);

// Block matrix functions
PetscErrorCode LaMEMMatrix_PMatBlockCreate(void *pm);
PetscErrorCode LaMEMMatrix_PMatBlockAssemble(void *pm);
PetscErrorCode LaMEMMatrix_PMatBlockPicardClean(Mat J, Vec x, Vec r);
PetscErrorCode LaMEMMatrix_PMatBlockPicardSchur(Mat J, Vec x, Vec r);
PetscErrorCode LaMEMMatrix_PMatBlockDestroy(void *pm);

// Service functions
void LaMEMMatrix_GetStiffMatDevProj(PetscScalar eta, PetscScalar diag, PetscScalar *v, PetscScalar *cf,
                                    PetscScalar dx, PetscScalar dy, PetscScalar dz,
                                    PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                    PetscScalar bdx, PetscScalar bdy, PetscScalar bdz);
void LaMEMMatrix_GetStiffMatClean(PetscScalar eta, PetscScalar diag, PetscScalar *v, PetscScalar *cf,
                                  PetscScalar dx, PetscScalar dy, PetscScalar dz,
                                  PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                  PetscScalar bdx, PetscScalar bdy, PetscScalar bdz);
void LaMEMMatrix_AddDensGradStabil(PetscScalar fssa, PetscScalar *v, PetscScalar rho, PetscScalar dt, 
                                   PetscScalar *grav, PetscScalar fdx, PetscScalar fdy, PetscScalar fdz,
                                   PetscScalar bdx, PetscScalar bdy, PetscScalar bdz);
void LaMEMMatrix_GetVelSchur(PetscScalar v[], PetscScalar d[], PetscScalar g[]);
void LaMEMMatrix_GetSubMat(PetscScalar v[], PetscScalar a[], PetscScalar d[], PetscScalar g[]);
void LaMEMMatrix_GetTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[]);
void LaMEMMatrix_ConstrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[]);

// Vector scatter functions
PetscErrorCode LaMEMMatrix_VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, PetscInt mode);

#ifdef __cplusplus
}
#endif

#endif