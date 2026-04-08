#ifndef LAMEMDIKE_C_H
#define LAMEMDIKE_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around existing dike functions
PetscErrorCode LaMEMDike_DBDikeCreate(void *dbdike, void *dbm, void *fb, void *jr, PetscBool PrintOutput);
PetscErrorCode LaMEMDike_DBReadDike(void *dbdike, void *dbm, void *fb, void *jr, PetscBool PrintOutput);
PetscErrorCode LaMEMDike_GetDikeContr(void *ctx, PetscScalar *phRat, PetscInt AirPhase, 
                                      PetscScalar *dikeRHS, PetscScalar *y_c, PetscInt J);
PetscErrorCode LaMEMDike_k_heatsource(void *jr, void *phases, PetscScalar *Tc, PetscScalar *phRat, 
                                      PetscScalar *k, PetscScalar *rho_A, PetscScalar *y_c, PetscInt J);
PetscErrorCode LaMEMDike_LocateDikeZones(void *actx);
PetscErrorCode LaMEMDike_ComputeSxxMagP(void *jr, PetscInt nD);
PetscErrorCode LaMEMDike_SmoothSxxEff(void *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2);
PetscErrorCode LaMEMDike_SetDikeZones(void *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2);
PetscErrorCode LaMEMDike_ReadRestart(void *dbdike, void *dbm, void *jr, void *ts, FILE *fp);
PetscErrorCode LaMEMDike_WriteRestart(void *jr, FILE *fp);
PetscErrorCode LaMEMDike_Destroy(void *jr);

#ifdef __cplusplus
}
#endif

#endif