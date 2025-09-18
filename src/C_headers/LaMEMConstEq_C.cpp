#include "LaMEMConstEq_C.h"
#include "../LaMEM.h"
#include "../constEq.h"
#include "../tssolve.h"
#include "../surf.h"
#include "../phase.h"
#include "../JacRes.h"
#include "../meltParam.h"
#include "../tools.h"
#include "../phase_transition.h"
#include "../scaling.h"
#include "../parsing.h"
#include "../bc.h"
#include "../dike.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMConstEq_SetUpConstEq(void *ctx, void *jr) {
    return setUpConstEq((ConstEqCtx*)ctx, (JacRes*)jr);
}

PetscErrorCode LaMEMConstEq_SetUpCtrlVol(void *ctx, PetscScalar *phRat, void *svDev, void *svBulk, 
                                         PetscScalar p, PetscScalar p_lith, PetscScalar p_pore, 
                                         PetscScalar T, PetscScalar DII, PetscScalar z, PetscScalar Le) {
    return setUpCtrlVol((ConstEqCtx*)ctx, phRat, (SolVarDev*)svDev, (SolVarBulk*)svBulk, 
                        p, p_lith, p_pore, T, DII, z, Le);
}

PetscErrorCode LaMEMConstEq_SetUpPhase(void *ctx, PetscInt ID) {
    return setUpPhase((ConstEqCtx*)ctx, ID);
}

PetscErrorCode LaMEMConstEq_DevConstEq(void *ctx) {
    return devConstEq((ConstEqCtx*)ctx);
}

PetscErrorCode LaMEMConstEq_GetPhaseVisc(void *ctx, PetscInt ID) {
    return getPhaseVisc((ConstEqCtx*)ctx, ID);
}

PetscScalar LaMEMConstEq_GetConsEqRes(PetscScalar eta, void *pctx) {
    return getConsEqRes(eta, pctx);
}

PetscScalar LaMEMConstEq_ApplyStrainSoft(void *soft, PetscInt ID, PetscScalar APS, PetscScalar Le, PetscScalar par) {
    return applyStrainSoft((Soft_t*)soft, ID, APS, Le, par);
}

PetscScalar LaMEMConstEq_GetI2Gdt(PetscInt numPhases, void *phases, PetscScalar *phRat, PetscScalar dt) {
    return getI2Gdt(numPhases, (Material_t*)phases, phRat, dt);
}

PetscErrorCode LaMEMConstEq_VolConstEq(void *ctx) {
    return volConstEq((ConstEqCtx*)ctx);
}

PetscErrorCode LaMEMConstEq_CellConstEq(void *ctx, void *svCell, PetscScalar dxx, PetscScalar dyy, PetscScalar dzz,
                                        PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz, 
                                        PetscScalar *gres, PetscScalar *rho, PetscScalar *dikeRHS) {
    return cellConstEq((ConstEqCtx*)ctx, (SolVarCell*)svCell, dxx, dyy, dzz, 
                       *sxx, *syy, *szz, *gres, *rho, *dikeRHS);
}

PetscErrorCode LaMEMConstEq_EdgeConstEq(void *ctx, void *svEdge, PetscScalar d, PetscScalar *s) {
    return edgeConstEq((ConstEqCtx*)ctx, (SolVarEdge*)svEdge, d, *s);
}

PetscErrorCode LaMEMConstEq_CheckConvConstEq(void *ctx) {
    return checkConvConstEq((ConstEqCtx*)ctx);
}

PetscErrorCode LaMEMConstEq_SetDataPhaseDiagram(void *pd, PetscScalar p, PetscScalar T, char pdn[]) {
    return setDataPhaseDiagram((PData*)pd, p, T, pdn);
}

}