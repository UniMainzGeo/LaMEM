#ifndef LAMEMTENSOR_C_H
#define LAMEMTENSOR_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Elastic stress rotation functions
void LaMEMTensor_GetRotationMatrix(void *R, PetscScalar dt, PetscScalar wx, PetscScalar wy, PetscScalar wz);
void LaMEMTensor_RotateStress(void *R, void *S, void *SR);

// Tensor copying and utility functions
void LaMEMTensor_Tensor2RSCopy(void *A, void *B);
void LaMEMTensor_Tensor2RNClear(void *A);
void LaMEMTensor_Tensor2RSClear(void *A);

// Tensor comparison and norm functions
PetscInt LaMEMTensor_Tensor2RNCheckEq(void *A, void *B, PetscScalar tol);
void LaMEMTensor_Tensor2RNNorm(void *A, PetscScalar *pk);
void LaMEMTensor_Tensor2RSNorm(void *A, PetscScalar *pk);

// Tensor arithmetic operations
void LaMEMTensor_Tensor2RNDivide(void *A, PetscScalar k);
void LaMEMTensor_Tensor2RNTrace(void *A);
void LaMEMTensor_Tensor2RNSym(void *A, void *B);
void LaMEMTensor_Tensor2RNProduct(void *A, void *B, void *C);
void LaMEMTensor_Tensor2RNTranspose(void *A, void *B);
void LaMEMTensor_Tensor2RNCopy(void *A, void *B);
void LaMEMTensor_Tensor2RNCopySym(void *A, void *B);
void LaMEMTensor_Tensor2RNUnit(void *A);
void LaMEMTensor_Tensor2RNSum3(void *A, PetscScalar ka, void *B, PetscScalar kb, 
                               void *C, PetscScalar kc, void *R);

// Tensor display functions
void LaMEMTensor_Tensor2RNView(void *A, const char *msg);
void LaMEMTensor_Tensor2RSView(void *A, const char *msg);

// Eigenvalue and spectral decomposition functions
PetscInt LaMEMTensor_Tensor2RNEigen(void *L, PetscScalar tol, PetscScalar eval[]);
PetscInt LaMEMTensor_Tensor2RSSpectral(void *A, PetscScalar eval[], PetscScalar evect[], 
                                       PetscScalar ttol, PetscScalar ltol, PetscInt itmax);
PetscErrorCode LaMEMTensor_Tensor2RS2DSpectral(PetscScalar axx, PetscScalar ayy, PetscScalar axy,
                                               PetscScalar *pa1, PetscScalar *pa2,
                                               PetscScalar v1[], PetscScalar v2[], PetscScalar tol);

// Infinite Strain Axis (ISA) calculation
PetscInt LaMEMTensor_getISA(void *pL, PetscScalar ISA[], PetscScalar *plnrm);

#ifdef __cplusplus
}
#endif

#endif