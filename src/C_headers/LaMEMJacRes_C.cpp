#include "LaMEMJacRes_C.h"
#include "../LaMEM.h"
#include "../JacRes.h"
#include "../parsing.h"
#include "../tssolve.h"
#include "../scaling.h"
#include "../fdstag.h"
#include "../surf.h"
#include "../bc.h"
#include "../phase.h"
#include "../constEq.h"
#include "../tools.h"
#include "../advect.h"
#include "../dike.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMJacRes_Create(void *jr, void *fb) {
    return JacResCreate((JacRes*)jr, (FB*)fb);
}

PetscErrorCode LaMEMJacRes_CreateData(void *jr) {
    return JacResCreateData((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_ReadRestart(void *jr, FILE *fp) {
    return JacResReadRestart((JacRes*)jr, fp);
}

PetscErrorCode LaMEMJacRes_WriteRestart(void *jr, FILE *fp) {
    return JacResWriteRestart((JacRes*)jr, fp);
}

PetscErrorCode LaMEMJacRes_Destroy(void *jr) {
    return JacResDestroy((JacRes*)jr);
}

// Main computation functions
PetscErrorCode LaMEMJacRes_FormResidual(void *jr, Vec x, Vec f) {
    return JacResFormResidual((JacRes*)jr, x, f);
}

PetscErrorCode LaMEMJacRes_GetI2Gdt(void *jr) {
    return JacResGetI2Gdt((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetPressShift(void *jr) {
    return JacResGetPressShift((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetEffStrainRate(void *jr) {
    return JacResGetEffStrainRate((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetVorticity(void *jr) {
    return JacResGetVorticity((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetResidual(void *jr) {
    return JacResGetResidual((JacRes*)jr);
}

// Solution copying functions
PetscErrorCode LaMEMJacRes_CopySol(void *jr, Vec x) {
    return JacResCopySol((JacRes*)jr, x);
}

PetscErrorCode LaMEMJacRes_CopyVel(void *jr, Vec x) {
    return JacResCopyVel((JacRes*)jr, x);
}

PetscErrorCode LaMEMJacRes_CopyPres(void *jr, Vec x) {
    return JacResCopyPres((JacRes*)jr, x);
}

// Pressure initialization
PetscErrorCode LaMEMJacRes_InitPres(void *jr, void *ts) {
    return JacResInitPres((JacRes*)jr, (TSSol*)ts);
}

PetscErrorCode LaMEMJacRes_InitLithPres(void *jr, void *actx, void *ts) {
    return JacResInitLithPres((JacRes*)jr, (AdvCtx*)actx, (TSSol*)ts);
}

// Residual copying and viewing
PetscErrorCode LaMEMJacRes_CopyRes(void *jr, Vec f) {
    return JacResCopyRes((JacRes*)jr, f);
}

PetscErrorCode LaMEMJacRes_CopyMomentumRes(void *jr, Vec f) {
    return JacResCopyMomentumRes((JacRes*)jr, f);
}

PetscErrorCode LaMEMJacRes_CopyContinuityRes(void *jr, Vec f) {
    return JacResCopyContinuityRes((JacRes*)jr, f);
}

PetscErrorCode LaMEMJacRes_ViewRes(void *jr) {
    return JacResViewRes((JacRes*)jr);
}

// Lithostatic and pore pressure
PetscErrorCode LaMEMJacRes_GetLithoStaticPressure(void *jr) {
    return JacResGetLithoStaticPressure((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetPorePressure(void *jr) {
    return JacResGetPorePressure((JacRes*)jr);
}

// Temperature-related functions (comment out if functions don't exist)
PetscErrorCode LaMEMJacRes_CheckTempParam(void *jr) {
    return JacResCheckTempParam((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_CreateTempParam(void *jr) {
    return JacResCreateTempParam((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_DestroyTempParam(void *jr) {
    return JacResDestroyTempParam((JacRes*)jr);
}

PetscErrorCode LaMEMJacRes_GetTempRes(void *jr, PetscScalar dt) {
    return JacResGetTempRes((JacRes*)jr, dt);
}

}