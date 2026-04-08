#ifndef LAMEMMULTIGRID_C_H
#define LAMEMMULTIGRID_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// MGLevel functions
PetscErrorCode LaMEMMultigrid_MGLevelCreate(void *lvl, void *fine, void *fs, void *bc);
PetscErrorCode LaMEMMultigrid_MGLevelDestroy(void *lvl);
PetscErrorCode LaMEMMultigrid_MGLevelInitEta(void *lvl, void *jr);
PetscErrorCode LaMEMMultigrid_MGLevelAverageEta(void *lvl);
PetscErrorCode LaMEMMultigrid_MGLevelRestrictEta(void *lvl, void *fine);
PetscErrorCode LaMEMMultigrid_MGLevelRestrictBC(void *lvl, void *fine, PetscBool no_restric_bc);
PetscErrorCode LaMEMMultigrid_MGLevelSetupRestrict(void *lvl, void *fine);
PetscErrorCode LaMEMMultigrid_MGLevelSetupProlong(void *lvl, void *fine);

// Main MG functions
PetscErrorCode LaMEMMultigrid_MGCreate(void *mg, void *jr);
PetscErrorCode LaMEMMultigrid_MGDestroy(void *mg);
PetscErrorCode LaMEMMultigrid_MGSetup(void *mg, Mat A);
PetscErrorCode LaMEMMultigrid_MGSetupCoarse(void *mg, Mat A);
PetscErrorCode LaMEMMultigrid_MGApply(PC pc, Vec x, Vec y);
PetscErrorCode LaMEMMultigrid_MGDumpMat(void *mg);
PetscErrorCode LaMEMMultigrid_MGGetNumLevels(void *mg);

// Helper functions for matrix row setup
void LaMEMMultigrid_GetRowRestrict(PetscBool scale, PetscScalar parent, PetscInt n, 
                                   PetscInt idx[], PetscScalar bc[], PetscScalar v[], 
                                   PetscScalar vs[], PetscScalar eta_fine[], PetscScalar eta_crs);
void LaMEMMultigrid_GetRowProlong(PetscBool scale, PetscInt parent, PetscScalar parent_bc, 
                                  PetscInt n, PetscScalar bc[], PetscScalar v[], PetscScalar vs[], 
                                  PetscScalar eta_crs[], PetscScalar eta_fine);

#ifdef __cplusplus
}
#endif

#endif