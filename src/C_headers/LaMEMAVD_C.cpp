#include "LaMEMAVD_C.h"
#include "../LaMEM.h"
#include "../AVD.h"
#include "../advect.h"
#include "../scaling.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../tools.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMAVD_Create(void *A) {
    return AVDCreate((AVD*)A);
}

PetscErrorCode LaMEMAVD_Destroy(void *A) {
    return AVDDestroy((AVD*)A);
}

PetscErrorCode LaMEMAVD_CellInit(void *A) {
    return AVDCellInit((AVD*)A);
}

PetscErrorCode LaMEMAVD_ClaimCells(void *A, PetscInt ip) {
    return AVDClaimCells((AVD*)A, ip);
}

PetscErrorCode LaMEMAVD_UpdateChain(void *A, PetscInt ip) {
    return AVDUpdateChain((AVD*)A, ip);
}

PetscErrorCode LaMEMAVD_ReAlloc(void *chain, PetscInt buffer) {
    return AVDReAlloc((AVDChain*)chain, buffer);
}

PetscErrorCode LaMEMAVD_LoadPoints(void *actx, void *A, PetscInt ind) {
    return AVDLoadPoints((AdvCtx*)actx, (AVD*)A, ind);
}

PetscErrorCode LaMEMAVD_InjectDeletePoints(void *actx, void *A, PetscInt cellID) {
    return AVDInjectDeletePoints((AdvCtx*)actx, (AVD*)A, cellID);
}

PetscErrorCode LaMEMAVD_ExecuteMarkerInjection(void *actx, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind) {
    return AVDExecuteMarkerInjection((AdvCtx*)actx, npoints, xs, xe, ind);
}

PetscErrorCode LaMEMAVD_MarkerControl(void *actx) {
    return AVDMarkerControl((AdvCtx*)actx);
}

PetscErrorCode LaMEMAVD_MarkerControlMV(void *actx, PetscInt vtype) {
    // Convert integer to enum if needed, otherwise pass directly
    return AVDMarkerControlMV((AdvCtx*)actx, (VolumeCase)vtype);
}

PetscErrorCode LaMEMAVD_CheckCellsMV(void *actx, void *mv, PetscInt dir) {
    return AVDCheckCellsMV((AdvCtx*)actx, (MarkerVolume*)mv, dir);
}

PetscErrorCode LaMEMAVD_CreateMV(void *actx, void *mv, PetscInt dir) {
    return AVDCreateMV((AdvCtx*)actx, (MarkerVolume*)mv, dir);
}

PetscErrorCode LaMEMAVD_DestroyMV(void *mv) {
    return AVDDestroyMV((MarkerVolume*)mv);
}

PetscErrorCode LaMEMAVD_MapMarkersMV(void *actx, void *mv, PetscInt dir) {
    return AVDMapMarkersMV((AdvCtx*)actx, (MarkerVolume*)mv, dir);
}

PetscErrorCode LaMEMAVD_LoadPointsMV(void *actx, void *mv, void *A, PetscInt ind) {
    return AVDLoadPointsMV((AdvCtx*)actx, (MarkerVolume*)mv, (AVD*)A, ind);
}

PetscErrorCode LaMEMAVD_InjectPointsMV(void *actx, void *A) {
    return AVDInjectPointsMV((AdvCtx*)actx, (AVD*)A);
}

PetscErrorCode LaMEMAVD_DeletePointsMV(void *actx, void *A) {
    return AVDDeletePointsMV((AdvCtx*)actx, (AVD*)A);
}

PetscErrorCode LaMEMAVD_AlgorithmMV(void *actx, void *mv, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind, PetscInt nmin) {
    return AVDAlgorithmMV((AdvCtx*)actx, (MarkerVolume*)mv, npoints, xs, xe, ind, nmin);
}

//PetscInt LaMEMAVD_FindPointInCell(PetscScalar *px, PetscInt L, PetscInt R, PetscScalar x) {
//    return FindPointInCell(px, L, R, x);
//}

PetscScalar LaMEMAVD_DistanceTest(PetscScalar x0[3], PetscScalar x1[3], PetscScalar x2[3]) {
    return AVDDistanceTest(x0, x1, x2);
}

}