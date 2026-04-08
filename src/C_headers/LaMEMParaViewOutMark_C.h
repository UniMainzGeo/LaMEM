#ifndef LAMEMPARAVIEWOUTMARK_C_H
#define LAMEMPARAVIEWOUTMARK_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual marker ParaView output functions
PetscErrorCode LaMEMParaViewOutMark_PVMarkCreate(void *pvmark, void *fb);
PetscErrorCode LaMEMParaViewOutMark_PVMarkWriteTimeStep(void *pvmark, const char *dirName, PetscScalar ttime);
PetscErrorCode LaMEMParaViewOutMark_PVMarkWriteVTU(void *pvmark, const char *dirName);
PetscErrorCode LaMEMParaViewOutMark_PVMarkWritePVTU(void *pvmark, const char *dirName);

#ifdef __cplusplus
}
#endif

#endif