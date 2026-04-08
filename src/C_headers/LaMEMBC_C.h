#ifndef LAMEMBC_C_H
#define LAMEMBC_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Shift types
#define LAMEM_LOCAL_TO_GLOBAL  0
#define LAMEM_GLOBAL_TO_LOCAL  1

// Boundary condition creation and management
PetscErrorCode LaMEMBC_Create(void *bc, void *fb);
PetscErrorCode LaMEMBC_ReadRestart(void *bc, void *fp);
PetscErrorCode LaMEMBC_WriteRestart(void *bc, void *fp);
PetscErrorCode LaMEMBC_CreateData(void *bc);
PetscErrorCode LaMEMBC_Destroy(void *bc);
PetscErrorCode LaMEMBC_ReadFixCell(void *bc, void *fb);

// Main boundary condition application
PetscErrorCode LaMEMBC_Apply(void *bc);
PetscErrorCode LaMEMBC_ApplySPC(void *bc);
PetscErrorCode LaMEMBC_ShiftIndices(void *bc, PetscInt stype);

// Specific constraint applications
PetscErrorCode LaMEMBC_ApplyPres(void *bc);
PetscErrorCode LaMEMBC_ApplyTemp(void *bc);
PetscErrorCode LaMEMBC_ApplyVelDefault(void *bc);
PetscErrorCode LaMEMBC_ApplyBezier(void *bc);
PetscErrorCode LaMEMBC_ApplyBoundVel(void *bc);
PetscErrorCode LaMEMBC_ApplyVelBox(void *bc);
PetscErrorCode LaMEMBC_ApplyVelCylinder(void *bc);
PetscErrorCode LaMEMBC_ApplyPhase(void *bc);
PetscErrorCode LaMEMBC_ApplyCells(void *bc);
PetscErrorCode LaMEMBC_ListSPC(void *bc);
PetscErrorCode LaMEMBC_ApplyVelTPC(void *bc);
PetscErrorCode LaMEMBC_PlumeInflow(void *bc);
PetscErrorCode LaMEMBC_ApplyPermeablePressure(void *bc);
PetscErrorCode LaMEMBC_ApplyPlumePresure(void *bc);

// Bezier block functions
PetscErrorCode LaMEMBC_BlockCreate(void *bcb, void *scal, void *fb);
PetscErrorCode LaMEMBC_BlockGetPosition(void *bcb, PetscScalar t, PetscInt *f, PetscScalar x[]);
PetscErrorCode LaMEMBC_BlockGetPolygon(void *bcb, PetscScalar Xb[], PetscScalar *cpoly);

// Velocity box functions
PetscErrorCode LaMEMBC_VelBoxCreate(void *velbox, void *scal, void *fb);
PetscErrorCode LaMEMBC_VelBoxPrint(void *velbox, void *scal, PetscInt cnt);

// Velocity cylinder functions
PetscErrorCode LaMEMBC_VelCylinderCreate(void *velcyl, void *scal, void *fb);
PetscErrorCode LaMEMBC_VelCylinderPrint(void *velcyl, void *scal, PetscInt cnt);

// Service functions
PetscErrorCode LaMEMBC_GetBGStrainRates(void *bc, PetscScalar *Exx, PetscScalar *Eyy, PetscScalar *Ezz, 
                                        PetscScalar *Exy, PetscScalar *Eyz, PetscScalar *Exz,
                                        PetscScalar *Rxx, PetscScalar *Ryy, PetscScalar *Rzz);
PetscErrorCode LaMEMBC_GetVelins(void *bc);
PetscErrorCode LaMEMBC_GetTempBound(void *bc, PetscScalar *Tbot);
PetscErrorCode LaMEMBC_StretchGrid(void *bc);
PetscErrorCode LaMEMBC_OverridePhase(void *bc, PetscInt cellID, void *P);
PetscErrorCode LaMEMBC_GetAverageLithostatic(void *bc);
PetscScalar LaMEMBC_GetDensity(void *bc, PetscInt Phase, PetscScalar T, PetscScalar p);

#ifdef __cplusplus
}
#endif

#endif