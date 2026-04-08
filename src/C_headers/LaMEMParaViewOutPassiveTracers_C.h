#ifndef LAMEMPARAVIEWOUTPASSIVETRACERS_C_H
#define LAMEMPARAVIEWOUTPASSIVETRACERS_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual passive tracer ParaView output functions
PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrCreate(void *pvptr, void *fb);
PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWriteTimeStep(void *pvptr, const char *dirName, PetscScalar ttime);
PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWriteVTU(void *pvptr, const char *dirName);
PetscErrorCode LaMEMParaViewOutPassiveTracers_PVPtrWritePVTU(void *pvptr, const char *dirName);

#ifdef __cplusplus
}
#endif

#endif