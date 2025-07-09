#include "LaMEMJacResAux_C.h"
#include "../LaMEM.h"
#include "../scaling.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../surf.h"
#include "../phase.h"
#include "../tools.h"
#include "../Tensor.h"
#include "../parsing.h"

extern "C" {

// Thin wrappers - just cast and call the actual existing functions
PetscErrorCode LaMEMJacResAux_GetGradientVel(void *fs, PetscScalar ***lvx, PetscScalar ***lvy, PetscScalar ***lvz,
                                             PetscInt i, PetscInt j, PetscInt k, PetscInt sx, PetscInt sy, PetscInt sz,
                                             void *L, PetscScalar *vel, PetscScalar *pvnrm) {
    return getGradientVel((FDSTAG*)fs, lvx, lvy, lvz, i, j, k, sx, sy, sz, (Tensor2RN*)L, vel, pvnrm);
}

PetscErrorCode LaMEMJacResAux_GetSHmax(void *jr) {
    return JacResGetSHmax((JacRes*)jr);
}

PetscErrorCode LaMEMJacResAux_GetEHmax(void *jr) {
    return JacResGetEHmax((JacRes*)jr);
}

PetscErrorCode LaMEMJacResAux_GetOverPressure(void *jr, Vec lop) {
    return JacResGetOverPressure((JacRes*)jr, lop);
}

PetscErrorCode LaMEMJacResAux_GetLithoStaticPressure(void *jr) {
    return JacResGetLithoStaticPressure((JacRes*)jr);
}

PetscErrorCode LaMEMJacResAux_GetPorePressure(void *jr) {
    return JacResGetPorePressure((JacRes*)jr);
}

PetscErrorCode LaMEMJacResAux_GetPermea(void *jr, PetscInt bgPhase, PetscInt step, char *outfile) {
    return JacResGetPermea((JacRes*)jr, bgPhase, step, outfile);
}

}