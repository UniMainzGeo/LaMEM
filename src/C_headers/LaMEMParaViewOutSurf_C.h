#ifndef LAMEMPARAVIEWOUTSURF_C_H
#define LAMEMPARAVIEWOUTSURF_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual surface ParaView output functions
PetscErrorCode LaMEMParaViewOutSurf_PVSurfCreate(void *pvsurf, void *fb);
PetscErrorCode LaMEMParaViewOutSurf_PVSurfDestroy(void *pvsurf);
PetscErrorCode LaMEMParaViewOutSurf_PVSurfWriteTimeStep(void *pvsurf, const char *dirName, PetscScalar ttime);
PetscErrorCode LaMEMParaViewOutSurf_PVSurfWriteVTS(void *pvsurf, const char *dirName);
PetscErrorCode LaMEMParaViewOutSurf_PVSurfWritePVTS(void *pvsurf, const char *dirName);

#ifdef __cplusplus
}
#endif

#endif