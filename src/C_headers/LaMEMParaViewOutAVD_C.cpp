#include "LaMEMParaViewOutAVD_C.h"
#include "../LaMEM.h"
#include "../paraViewOutAVD.h"
#include "../paraViewOutBin.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../fdstag.h"
#include "../advect.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// AVDPoint3D functions
void LaMEMParaViewOutAVD_AVDPoint3DCreate(const PetscInt npoints, void **P) {
    AVDPoint3DCreate(npoints, (AVDPoint3D*)P);
}

void LaMEMParaViewOutAVD_AVDPoint3DDestroy(void **P) {
    AVDPoint3DDestroy((AVDPoint3D*)P);
}

// AVDCell3D functions
void LaMEMParaViewOutAVD_AVDCell3DCreate(const PetscInt mx, const PetscInt my, const PetscInt mz, void **C) {
    AVDCell3DCreate(mx, my, mz, (AVDCell3D*)C);
}

void LaMEMParaViewOutAVD_AVDCell3DDestroy(void **C) {
    AVDCell3DDestroy((AVDCell3D*)C);
}

// AVDChain3D functions
void LaMEMParaViewOutAVD_AVDChain3DCreate(const PetscInt npoints, const PetscInt buffer, void **CH) {
    AVDChain3DCreate(npoints, buffer, (AVDChain3D*)CH);
}

void LaMEMParaViewOutAVD_AVDChain3DDestroy(const PetscInt npoints, void **CH) {
    AVDChain3DDestroy(npoints, (AVDChain3D*)CH);
}

// Main AVD3D functions
PetscErrorCode LaMEMParaViewOutAVD_AVDViewCreate(void **A, void *actx, PetscInt refine) {
    return AVDViewCreate((AVD3D*)A, (AdvCtx*)actx, refine);
}

void LaMEMParaViewOutAVD_AVD3DDestroy(void **A) {
    AVD3DDestroy((AVD3D*)A);
}

void LaMEMParaViewOutAVD_AVD3DAllocate(const PetscInt mx, const PetscInt my, const PetscInt mz,
                                       const PetscInt buffer, const PetscInt npoints, void **A) {
    AVD3DAllocate(mx, my, mz, buffer, npoints, (AVD3D*)A);
}

PetscErrorCode LaMEMParaViewOutAVD_AVD3DSetParallelExtent(void *A, PetscInt M, PetscInt N, PetscInt P) {
    return AVD3DSetParallelExtent((AVD3D)A, M, N, P);
}

void LaMEMParaViewOutAVD_AVD3DSetDomainSize(void *A, const PetscScalar x0, const PetscScalar x1,
                                            const PetscScalar y0, const PetscScalar y1,
                                            const PetscScalar z0, const PetscScalar z1) {
    AVD3DSetDomainSize((AVD3D)A, x0, x1, y0, y1, z0, z1);
}

PetscErrorCode LaMEMParaViewOutAVD_AVD3DLoadPoints(void *A, void *actx) {
    return AVD3DLoadPoints((AVD3D)A, (AdvCtx*)actx);
}

void LaMEMParaViewOutAVD_AVD3DResetCells(void *A) {
    AVD3DResetCells((AVD3D)A);
}

PetscErrorCode LaMEMParaViewOutAVD_AVD3DInit(void *A) {
    return AVD3DInit((AVD3D)A);
}

void LaMEMParaViewOutAVD_AVD3DClaimCells(void *A, const PetscInt p_i) {
    AVD3DClaimCells((AVD3D)A, p_i);
}

void LaMEMParaViewOutAVD_AVD3DUpdateChain(void *A, const PetscInt p_i) {
    AVD3DUpdateChain((AVD3D)A, p_i);
}

// PVAVD (ParaView AVD output) functions
PetscErrorCode LaMEMParaViewOutAVD_PVAVDCreate(void *pvavd, void *fb) {
    return PVAVDCreate((PVAVD*)pvavd, (FB*)fb);
}

PetscErrorCode LaMEMParaViewOutAVD_PVAVDWriteTimeStep(void *pvavd, const char *dirName, PetscScalar ttime) {
    return PVAVDWriteTimeStep((PVAVD*)pvavd, dirName, ttime);
}

PetscErrorCode LaMEMParaViewOutAVD_PVAVDWritePVTR(void *pvavd, void *A, const char *dirName) {
    return PVAVDWritePVTR((PVAVD*)pvavd, (AVD3D)A, dirName);
}

PetscErrorCode LaMEMParaViewOutAVD_PVAVDWriteVTR(void *pvavd, void *A, const char *dirName) {
    return PVAVDWriteVTR((PVAVD*)pvavd, (AVD3D)A, dirName);
}

}