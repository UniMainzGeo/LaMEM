#ifndef LAMEMCONSTEQ_C_H
#define LAMEMCONSTEQ_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around existing constEq functions
PetscErrorCode LaMEMConstEq_SetUpConstEq(void *ctx, void *jr);
PetscErrorCode LaMEMConstEq_SetUpCtrlVol(void *ctx, PetscScalar *phRat, void *svDev, void *svBulk, 
                                         PetscScalar p, PetscScalar p_lith, PetscScalar p_pore, 
                                         PetscScalar T, PetscScalar DII, PetscScalar z, PetscScalar Le);
PetscErrorCode LaMEMConstEq_SetUpPhase(void *ctx, PetscInt ID);
PetscErrorCode LaMEMConstEq_DevConstEq(void *ctx);
PetscErrorCode LaMEMConstEq_GetPhaseVisc(void *ctx, PetscInt ID);
PetscScalar LaMEMConstEq_GetConsEqRes(PetscScalar eta, void *pctx);
PetscScalar LaMEMConstEq_ApplyStrainSoft(void *soft, PetscInt ID, PetscScalar APS, PetscScalar Le, PetscScalar par);
PetscScalar LaMEMConstEq_GetI2Gdt(PetscInt numPhases, void *phases, PetscScalar *phRat, PetscScalar dt);
PetscErrorCode LaMEMConstEq_VolConstEq(void *ctx);
PetscErrorCode LaMEMConstEq_CellConstEq(void *ctx, void *svCell, PetscScalar dxx, PetscScalar dyy, PetscScalar dzz,
                                        PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz, 
                                        PetscScalar *gres, PetscScalar *rho, PetscScalar *dikeRHS);
PetscErrorCode LaMEMConstEq_EdgeConstEq(void *ctx, void *svEdge, PetscScalar d, PetscScalar *s);
PetscErrorCode LaMEMConstEq_CheckConvConstEq(void *ctx);
PetscErrorCode LaMEMConstEq_SetDataPhaseDiagram(void *pd, PetscScalar p, PetscScalar T, char pdn[]);

#ifdef __cplusplus
}
#endif

#endif