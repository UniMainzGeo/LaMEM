#include "LaMEMBC_C.h"
#include "../LaMEM.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../fdstag.h"
#include "../tools.h"
#include "../advect.h"
#include "../phase.h"
#include "../constEq.h"
#include "../surf.h"

extern "C" {

// Boundary condition creation and management
PetscErrorCode LaMEMBC_Create(void *bc, void *fb) {
    return BCCreate((BCCtx*)bc, (FB*)fb);
}

PetscErrorCode LaMEMBC_ReadRestart(void *bc, void *fp) {
    return BCReadRestart((BCCtx*)bc, (FILE*)fp);
}

PetscErrorCode LaMEMBC_WriteRestart(void *bc, void *fp) {
    return BCWriteRestart((BCCtx*)bc, (FILE*)fp);
}

PetscErrorCode LaMEMBC_CreateData(void *bc) {
    return BCCreateData((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_Destroy(void *bc) {
    return BCDestroy((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ReadFixCell(void *bc, void *fb) {
    return BCReadFixCell((BCCtx*)bc, (FB*)fb);
}

// Main boundary condition application
PetscErrorCode LaMEMBC_Apply(void *bc) {
    return BCApply((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplySPC(void *bc) {
    return BCApplySPC((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ShiftIndices(void *bc, PetscInt stype) {
    ShiftType shift_type;
    switch(stype) {
        case 0: shift_type = _LOCAL_TO_GLOBAL_; break;  // LAMEM_LOCAL_TO_GLOBAL
        case 1: shift_type = _GLOBAL_TO_LOCAL_; break;  // LAMEM_GLOBAL_TO_LOCAL
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return BCShiftIndices((BCCtx*)bc, shift_type);
}

// Specific constraint applications
PetscErrorCode LaMEMBC_ApplyPres(void *bc) {
    return BCApplyPres((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyTemp(void *bc) {
    return BCApplyTemp((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyVelDefault(void *bc) {
    return BCApplyVelDefault((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyBezier(void *bc) {
    return BCApplyBezier((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyBoundVel(void *bc) {
    return BCApplyBoundVel((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyVelBox(void *bc) {
    return BCApplyVelBox((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyVelCylinder(void *bc) {
    return BCApplyVelCylinder((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyPhase(void *bc) {
    return BCApplyPhase((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyCells(void *bc) {
    return BCApplyCells((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ListSPC(void *bc) {
    return BCListSPC((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyVelTPC(void *bc) {
    return BCApplyVelTPC((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_PlumeInflow(void *bc) {
    return BC_Plume_inflow((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyPermeablePressure(void *bc) {
    return BCApply_Permeable_Pressure((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_ApplyPlumePresure(void *bc) {
    return BCApplyPres_Plume_Pressure((BCCtx*)bc);
}

// Bezier block functions
PetscErrorCode LaMEMBC_BlockCreate(void *bcb, void *scal, void *fb) {
    return BCBlockCreate((BCBlock*)bcb, (Scaling*)scal, (FB*)fb);
}

PetscErrorCode LaMEMBC_BlockGetPosition(void *bcb, PetscScalar t, PetscInt *f, PetscScalar x[]) {
    return BCBlockGetPosition((BCBlock*)bcb, t, f, x);
}

PetscErrorCode LaMEMBC_BlockGetPolygon(void *bcb, PetscScalar Xb[], PetscScalar *cpoly) {
    return BCBlockGetPolygon((BCBlock*)bcb, Xb, cpoly);
}

// Velocity box functions
PetscErrorCode LaMEMBC_VelBoxCreate(void *velbox, void *scal, void *fb) {
    return VelBoxCreate((VelBox*)velbox, (Scaling*)scal, (FB*)fb);
}

PetscErrorCode LaMEMBC_VelBoxPrint(void *velbox, void *scal, PetscInt cnt) {
    return VelBoxPrint((VelBox*)velbox, (Scaling*)scal, cnt);
}

// Velocity cylinder functions
PetscErrorCode LaMEMBC_VelCylinderCreate(void *velcyl, void *scal, void *fb) {
    return VelCylinderCreate((VelCylinder*)velcyl, (Scaling*)scal, (FB*)fb);
}

PetscErrorCode LaMEMBC_VelCylinderPrint(void *velcyl, void *scal, PetscInt cnt) {
    return VelCylinderPrint((VelCylinder*)velcyl, (Scaling*)scal, cnt);
}

// Service functions
PetscErrorCode LaMEMBC_GetBGStrainRates(void *bc, PetscScalar *Exx, PetscScalar *Eyy, PetscScalar *Ezz, 
                                        PetscScalar *Exy, PetscScalar *Eyz, PetscScalar *Exz,
                                        PetscScalar *Rxx, PetscScalar *Ryy, PetscScalar *Rzz) {
    return BCGetBGStrainRates((BCCtx*)bc, Exx, Eyy, Ezz, Exy, Eyz, Exz, Rxx, Ryy, Rzz);
}

PetscErrorCode LaMEMBC_GetVelins(void *bc) {
    return BCGetVelins((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_GetTempBound(void *bc, PetscScalar *Tbot) {
    return BCGetTempBound((BCCtx*)bc, Tbot);
}

PetscErrorCode LaMEMBC_StretchGrid(void *bc) {
    return BCStretchGrid((BCCtx*)bc);
}

PetscErrorCode LaMEMBC_OverridePhase(void *bc, PetscInt cellID, void *P) {
    return BCOverridePhase((BCCtx*)bc, cellID, (Marker*)P);
}

PetscErrorCode LaMEMBC_GetAverageLithostatic(void *bc) {
    return GetAverageLithostatic((BCCtx*)bc);
}

PetscScalar LaMEMBC_GetDensity(void *bc, PetscInt Phase, PetscScalar T, PetscScalar p) {
    return GetDensity((BCCtx*)bc, Phase, T, p);
}

}