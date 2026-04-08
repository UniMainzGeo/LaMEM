#ifndef LAMEMPARAVIEWOUTAVD_C_H
#define LAMEMPARAVIEWOUTAVD_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// AVD constants
#define LAMEM_AVD_TRUE  1
#define LAMEM_AVD_FALSE 0
#define LAMEM_AVD_CELL_MASK -999
#define LAMEM_AVD_CELL_UNCLAIMED -1

// AVDPoint3D functions
void LaMEMParaViewOutAVD_AVDPoint3DCreate(const PetscInt npoints, void **P);
void LaMEMParaViewOutAVD_AVDPoint3DDestroy(void **P);

// AVDCell3D functions
void LaMEMParaViewOutAVD_AVDCell3DCreate(const PetscInt mx, const PetscInt my, const PetscInt mz, void **C);
void LaMEMParaViewOutAVD_AVDCell3DDestroy(void **C);

// AVDChain3D functions
void LaMEMParaViewOutAVD_AVDChain3DCreate(const PetscInt npoints, const PetscInt buffer, void **CH);
void LaMEMParaViewOutAVD_AVDChain3DDestroy(const PetscInt npoints, void **CH);

// Main AVD3D functions
PetscErrorCode LaMEMParaViewOutAVD_AVDViewCreate(void **A, void *actx, PetscInt refine);
void LaMEMParaViewOutAVD_AVD3DDestroy(void **A);
void LaMEMParaViewOutAVD_AVD3DAllocate(const PetscInt mx, const PetscInt my, const PetscInt mz,
                                       const PetscInt buffer, const PetscInt npoints, void **A);
PetscErrorCode LaMEMParaViewOutAVD_AVD3DSetParallelExtent(void *A, PetscInt M, PetscInt N, PetscInt P);
void LaMEMParaViewOutAVD_AVD3DSetDomainSize(void *A, const PetscScalar x0, const PetscScalar x1,
                                            const PetscScalar y0, const PetscScalar y1,
                                            const PetscScalar z0, const PetscScalar z1);
PetscErrorCode LaMEMParaViewOutAVD_AVD3DLoadPoints(void *A, void *actx);
void LaMEMParaViewOutAVD_AVD3DResetCells(void *A);
PetscErrorCode LaMEMParaViewOutAVD_AVD3DInit(void *A);
void LaMEMParaViewOutAVD_AVD3DClaimCells(void *A, const PetscInt p_i);
void LaMEMParaViewOutAVD_AVD3DUpdateChain(void *A, const PetscInt p_i);

// PVAVD (ParaView AVD output) functions
PetscErrorCode LaMEMParaViewOutAVD_PVAVDCreate(void *pvavd, void *fb);
PetscErrorCode LaMEMParaViewOutAVD_PVAVDWriteTimeStep(void *pvavd, const char *dirName, PetscScalar ttime);
PetscErrorCode LaMEMParaViewOutAVD_PVAVDWritePVTR(void *pvavd, void *A, const char *dirName);
PetscErrorCode LaMEMParaViewOutAVD_PVAVDWriteVTR(void *pvavd, void *A, const char *dirName);

#ifdef __cplusplus
}
#endif

#endif