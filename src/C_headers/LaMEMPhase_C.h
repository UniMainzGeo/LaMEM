#ifndef LAMEMPHASE_C_H
#define LAMEMPHASE_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Experiment type constants for rheological profiles
#define LAMEM_UNIAXIAL     0
#define LAMEM_SIMPLESHEAR  1
#define LAMEM_NONE         2

// Main material database functions
PetscErrorCode LaMEMPhase_DBMatCreate(void *dbm, void *fb, PetscBool PrintOutput);
PetscErrorCode LaMEMPhase_DBMatReadSoft(void *dbm, void *fb, PetscBool PrintOutput);
PetscErrorCode LaMEMPhase_DBMatReadPhase(void *dbm, void *fb, PetscBool PrintOutput);
PetscErrorCode LaMEMPhase_DBMatOverwriteWithGlobalVariables(void *dbm, void *fb);

// Material property printing and utilities
void LaMEMPhase_MatPrintScalParam(PetscScalar par, const char key[], const char label[],
                                  void *scal, const char title[], PetscInt *print_title);
PetscErrorCode LaMEMPhase_PrintMatProp(void *MatProp);

// Predefined rheological profile functions
PetscErrorCode LaMEMPhase_GetProfileName(void *fb, void *scal, char name[], const char key[]);
PetscErrorCode LaMEMPhase_SetDiffProfile(void *m, char name[]);
PetscErrorCode LaMEMPhase_SetDislProfile(void *m, char name[]);
PetscErrorCode LaMEMPhase_SetPeirProfile(void *m, char name[]);

// Experimental correction functions
PetscErrorCode LaMEMPhase_CorrExpPreFactor(PetscScalar *B, PetscScalar n, PetscInt type, PetscInt MPa);
PetscErrorCode LaMEMPhase_CorrExpStressStrainRate(PetscScalar *D, PetscScalar *S, PetscInt type, PetscInt MPa);

#ifdef __cplusplus
}
#endif

#endif