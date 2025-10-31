#ifndef LAMEMJACRES_C_H
#define LAMEMJACRES_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around existing JacRes functions
PetscErrorCode LaMEMJacRes_Create(void *jr, void *fb);
PetscErrorCode LaMEMJacRes_CreateData(void *jr);
PetscErrorCode LaMEMJacRes_ReadRestart(void *jr, FILE *fp);
PetscErrorCode LaMEMJacRes_WriteRestart(void *jr, FILE *fp);
PetscErrorCode LaMEMJacRes_Destroy(void *jr);

// Main computation functions
PetscErrorCode LaMEMJacRes_FormResidual(void *jr, Vec x, Vec f);
PetscErrorCode LaMEMJacRes_GetI2Gdt(void *jr);
PetscErrorCode LaMEMJacRes_GetPressShift(void *jr);
PetscErrorCode LaMEMJacRes_GetEffStrainRate(void *jr);
PetscErrorCode LaMEMJacRes_GetVorticity(void *jr);
PetscErrorCode LaMEMJacRes_GetResidual(void *jr);

// Solution copying functions
PetscErrorCode LaMEMJacRes_CopySol(void *jr, Vec x);
PetscErrorCode LaMEMJacRes_CopyVel(void *jr, Vec x);
PetscErrorCode LaMEMJacRes_CopyPres(void *jr, Vec x);

// Pressure initialization
PetscErrorCode LaMEMJacRes_InitPres(void *jr, void *ts);
PetscErrorCode LaMEMJacRes_InitLithPres(void *jr, void *actx, void *ts);

// Residual copying and viewing
PetscErrorCode LaMEMJacRes_CopyRes(void *jr, Vec f);
PetscErrorCode LaMEMJacRes_CopyMomentumRes(void *jr, Vec f);
PetscErrorCode LaMEMJacRes_CopyContinuityRes(void *jr, Vec f);
PetscErrorCode LaMEMJacRes_ViewRes(void *jr);

// Lithostatic and pore pressure
PetscErrorCode LaMEMJacRes_GetLithoStaticPressure(void *jr);
PetscErrorCode LaMEMJacRes_GetPorePressure(void *jr);

// Temperature-related functions (if available)
PetscErrorCode LaMEMJacRes_CheckTempParam(void *jr);
PetscErrorCode LaMEMJacRes_CreateTempParam(void *jr);
PetscErrorCode LaMEMJacRes_DestroyTempParam(void *jr);
PetscErrorCode LaMEMJacRes_GetTempRes(void *jr, PetscScalar dt);

#ifdef __cplusplus
}
#endif

#endif