#ifndef LAMEMJACRESAUX_C_H
#define LAMEMJACRESAUX_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual JacResAux functions from the source
PetscErrorCode LaMEMJacResAux_GetGradientVel(void *fs, PetscScalar ***lvx, PetscScalar ***lvy, PetscScalar ***lvz,
                                             PetscInt i, PetscInt j, PetscInt k, PetscInt sx, PetscInt sy, PetscInt sz,
                                             void *L, PetscScalar *vel, PetscScalar *pvnrm);
PetscErrorCode LaMEMJacResAux_GetSHmax(void *jr);
PetscErrorCode LaMEMJacResAux_GetEHmax(void *jr);
PetscErrorCode LaMEMJacResAux_GetOverPressure(void *jr, Vec lop);
PetscErrorCode LaMEMJacResAux_GetLithoStaticPressure(void *jr);
PetscErrorCode LaMEMJacResAux_GetPorePressure(void *jr);
PetscErrorCode LaMEMJacResAux_GetPermea(void *jr, PetscInt bgPhase, PetscInt step, char *outfile);

#ifdef __cplusplus
}
#endif

#endif