#ifndef LAMEMAVD_C_H
#define LAMEMAVD_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around existing AVD functions
PetscErrorCode LaMEMAVD_Create(void *A);
PetscErrorCode LaMEMAVD_Destroy(void *A);
PetscErrorCode LaMEMAVD_CellInit(void *A);
PetscErrorCode LaMEMAVD_ClaimCells(void *A, PetscInt ip);
PetscErrorCode LaMEMAVD_UpdateChain(void *A, PetscInt ip);
PetscErrorCode LaMEMAVD_ReAlloc(void *chain, PetscInt buffer);
PetscErrorCode LaMEMAVD_LoadPoints(void *actx, void *A, PetscInt ind);
PetscErrorCode LaMEMAVD_InjectDeletePoints(void *actx, void *A, PetscInt cellID);
PetscErrorCode LaMEMAVD_ExecuteMarkerInjection(void *actx, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind);
PetscErrorCode LaMEMAVD_MarkerControl(void *actx);
PetscErrorCode LaMEMAVD_MarkerControlMV(void *actx, PetscInt vtype);
PetscErrorCode LaMEMAVD_CheckCellsMV(void *actx, void *mv, PetscInt dir);
PetscErrorCode LaMEMAVD_CreateMV(void *actx, void *mv, PetscInt dir);
PetscErrorCode LaMEMAVD_DestroyMV(void *mv);
PetscErrorCode LaMEMAVD_MapMarkersMV(void *actx, void *mv, PetscInt dir);
PetscErrorCode LaMEMAVD_LoadPointsMV(void *actx, void *mv, void *A, PetscInt ind);
PetscErrorCode LaMEMAVD_InjectPointsMV(void *actx, void *A);
PetscErrorCode LaMEMAVD_DeletePointsMV(void *actx, void *A);
PetscErrorCode LaMEMAVD_AlgorithmMV(void *actx, void *mv, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind, PetscInt nmin);
PetscInt LaMEMAVD_FindPointInCell(PetscScalar *px, PetscInt L, PetscInt R, PetscScalar x);
PetscScalar LaMEMAVD_DistanceTest(PetscScalar x0[3], PetscScalar x1[3], PetscScalar x2[3]);

#ifdef __cplusplus
}
#endif

#endif