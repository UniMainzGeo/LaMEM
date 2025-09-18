#include "LaMEMDike_C.h"
#include "../LaMEM.h"
#include "../phase.h"
#include "../parsing.h"
#include "../JacRes.h"
#include "../dike.h"
#include "../constEq.h"
#include "../bc.h"
#include "../tssolve.h"
#include "../scaling.h"
#include "../fdstag.h"
#include "../tools.h"
#include "../surf.h"
#include "../advect.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMDike_DBDikeCreate(void *dbdike, void *dbm, void *fb, void *jr, PetscBool PrintOutput) {
    return DBDikeCreate((DBPropDike*)dbdike, (DBMat*)dbm, (FB*)fb, (JacRes*)jr, PrintOutput);
}

PetscErrorCode LaMEMDike_DBReadDike(void *dbdike, void *dbm, void *fb, void *jr, PetscBool PrintOutput) {
    return DBReadDike((DBPropDike*)dbdike, (DBMat*)dbm, (FB*)fb, (JacRes*)jr, PrintOutput);
}

PetscErrorCode LaMEMDike_GetDikeContr(void *ctx, PetscScalar *phRat, PetscInt AirPhase, 
                                      PetscScalar *dikeRHS, PetscScalar *y_c, PetscInt J) {
    return GetDikeContr((ConstEqCtx*)ctx, phRat, AirPhase, *dikeRHS, *y_c, J);
}

PetscErrorCode LaMEMDike_k_heatsource(void *jr, void *phases, PetscScalar *Tc, PetscScalar *phRat, 
                                      PetscScalar *k, PetscScalar *rho_A, PetscScalar *y_c, PetscInt J) {
    return Dike_k_heatsource((JacRes*)jr, (Material_t*)phases, *Tc, phRat, *k, *rho_A, *y_c, J);
}

PetscErrorCode LaMEMDike_LocateDikeZones(void *actx) {
    return Locate_Dike_Zones((AdvCtx*)actx);
}

PetscErrorCode LaMEMDike_ComputeSxxMagP(void *jr, PetscInt nD) {
    return Compute_sxx_magP((JacRes*)jr, nD);
}

PetscErrorCode LaMEMDike_SmoothSxxEff(void *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2) {
    return Smooth_sxx_eff((JacRes*)jr, nD, nPtr, j1, j2);
}

PetscErrorCode LaMEMDike_SetDikeZones(void *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2) {
    return Set_dike_zones((JacRes*)jr, nD, nPtr, j1, j2);
}

PetscErrorCode LaMEMDike_ReadRestart(void *dbdike, void *dbm, void *jr, void *ts, FILE *fp) {
    return DynamicDike_ReadRestart((DBPropDike*)dbdike, (DBMat*)dbm, (JacRes*)jr, (TSSol*)ts, fp);
}

PetscErrorCode LaMEMDike_WriteRestart(void *jr, FILE *fp) {
    return DynamicDike_WriteRestart((JacRes*)jr, fp);
}

PetscErrorCode LaMEMDike_Destroy(void *jr) {
    return DynamicDike_Destroy((JacRes*)jr);
}

}