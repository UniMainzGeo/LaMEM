#include "LaMEMTensor_C.h"
#include "../LaMEM.h"
#include "../Tensor.h"
#include "../tools.h"

extern "C" {

// Elastic stress rotation functions
void LaMEMTensor_GetRotationMatrix(void *R, PetscScalar dt, PetscScalar wx, PetscScalar wy, PetscScalar wz) {
    GetRotationMatrix((Tensor2RN*)R, dt, wx, wy, wz);
}

void LaMEMTensor_RotateStress(void *R, void *S, void *SR) {
    RotateStress((Tensor2RN*)R, (Tensor2RS*)S, (Tensor2RS*)SR);
}

// Tensor copying and utility functions
void LaMEMTensor_Tensor2RSCopy(void *A, void *B) {
    Tensor2RSCopy((Tensor2RS*)A, (Tensor2RS*)B);
}

void LaMEMTensor_Tensor2RNClear(void *A) {
    Tensor2RNClear((Tensor2RN*)A);
}

void LaMEMTensor_Tensor2RSClear(void *A) {
    Tensor2RSClear((Tensor2RS*)A);
}

// Tensor comparison and norm functions
PetscInt LaMEMTensor_Tensor2RNCheckEq(void *A, void *B, PetscScalar tol) {
    return Tensor2RNCheckEq((Tensor2RN*)A, (Tensor2RN*)B, tol);
}

void LaMEMTensor_Tensor2RNNorm(void *A, PetscScalar *pk) {
    Tensor2RNNorm((Tensor2RN*)A, pk);
}

void LaMEMTensor_Tensor2RSNorm(void *A, PetscScalar *pk) {
    Tensor2RSNorm((Tensor2RS*)A, pk);
}

// Tensor arithmetic operations
void LaMEMTensor_Tensor2RNDivide(void *A, PetscScalar k) {
    Tensor2RNDivide((Tensor2RN*)A, k);
}

void LaMEMTensor_Tensor2RNTrace(void *A) {
    Tensor2RNTrace((Tensor2RN*)A);
}

void LaMEMTensor_Tensor2RNSym(void *A, void *B) {
    Tensor2RNSym((Tensor2RN*)A, (Tensor2RN*)B);
}

void LaMEMTensor_Tensor2RNProduct(void *A, void *B, void *C) {
    Tensor2RNProduct((Tensor2RN*)A, (Tensor2RN*)B, (Tensor2RN*)C);
}

void LaMEMTensor_Tensor2RNTranspose(void *A, void *B) {
    Tensor2RNTranspose((Tensor2RN*)A, (Tensor2RN*)B);
}

void LaMEMTensor_Tensor2RNCopy(void *A, void *B) {
    Tensor2RNCopy((Tensor2RN*)A, (Tensor2RN*)B);
}

void LaMEMTensor_Tensor2RNCopySym(void *A, void *B) {
    Tensor2RNCopySym((Tensor2RN*)A, (Tensor2RS*)B);
}

void LaMEMTensor_Tensor2RNUnit(void *A) {
    Tensor2RNUnit((Tensor2RN*)A);
}

void LaMEMTensor_Tensor2RNSum3(void *A, PetscScalar ka, void *B, PetscScalar kb, 
                               void *C, PetscScalar kc, void *R) {
    Tensor2RNSum3((Tensor2RN*)A, ka, (Tensor2RN*)B, kb, (Tensor2RN*)C, kc, (Tensor2RN*)R);
}

// Tensor display functions
void LaMEMTensor_Tensor2RNView(void *A, const char *msg) {
    Tensor2RNView((Tensor2RN*)A, msg);
}

void LaMEMTensor_Tensor2RSView(void *A, const char *msg) {
    Tensor2RSView((Tensor2RS*)A, msg);
}

// Eigenvalue and spectral decomposition functions
PetscInt LaMEMTensor_Tensor2RNEigen(void *L, PetscScalar tol, PetscScalar eval[]) {
    return Tensor2RNEigen((Tensor2RN*)L, tol, eval);
}

PetscInt LaMEMTensor_Tensor2RSSpectral(void *A, PetscScalar eval[], PetscScalar evect[], 
                                       PetscScalar ttol, PetscScalar ltol, PetscInt itmax) {
    return Tensor2RSSpectral((Tensor2RS*)A, eval, evect, ttol, ltol, itmax);
}

PetscErrorCode LaMEMTensor_Tensor2RS2DSpectral(PetscScalar axx, PetscScalar ayy, PetscScalar axy,
                                               PetscScalar *pa1, PetscScalar *pa2,
                                               PetscScalar v1[], PetscScalar v2[], PetscScalar tol) {
    return Tensor2RS2DSpectral(axx, ayy, axy, pa1, pa2, v1, v2, tol);
}

// Infinite Strain Axis (ISA) calculation
PetscInt LaMEMTensor_getISA(void *pL, PetscScalar ISA[], PetscScalar *plnrm) {
    return getISA((Tensor2RN*)pL, ISA, plnrm);
}

}