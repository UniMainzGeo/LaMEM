#include "LaMEMMultigrid_C.h"
#include "../LaMEM.h"
#include "../fdstag.h"
#include "../multigrid.h"
#include "../matrix.h"
#include "../JacRes.h"
#include "../bc.h"
#include "../tools.h"

extern "C" {

// MGLevel functions
PetscErrorCode LaMEMMultigrid_MGLevelCreate(void *lvl, void *fine, void *fs, void *bc) {
    return MGLevelCreate((MGLevel*)lvl, (MGLevel*)fine, (FDSTAG*)fs, (BCCtx*)bc);
}

PetscErrorCode LaMEMMultigrid_MGLevelDestroy(void *lvl) {
    return MGLevelDestroy((MGLevel*)lvl);
}

PetscErrorCode LaMEMMultigrid_MGLevelInitEta(void *lvl, void *jr) {
    return MGLevelInitEta((MGLevel*)lvl, (JacRes*)jr);
}

PetscErrorCode LaMEMMultigrid_MGLevelAverageEta(void *lvl) {
    return MGLevelAverageEta((MGLevel*)lvl);
}

PetscErrorCode LaMEMMultigrid_MGLevelRestrictEta(void *lvl, void *fine) {
    return MGLevelRestrictEta((MGLevel*)lvl, (MGLevel*)fine);
}

PetscErrorCode LaMEMMultigrid_MGLevelRestrictBC(void *lvl, void *fine, PetscBool no_restric_bc) {
    return MGLevelRestrictBC((MGLevel*)lvl, (MGLevel*)fine, no_restric_bc);
}

PetscErrorCode LaMEMMultigrid_MGLevelSetupRestrict(void *lvl, void *fine) {
    return MGLevelSetupRestrict((MGLevel*)lvl, (MGLevel*)fine);
}

PetscErrorCode LaMEMMultigrid_MGLevelSetupProlong(void *lvl, void *fine) {
    return MGLevelSetupProlong((MGLevel*)lvl, (MGLevel*)fine);
}

// Main MG functions
PetscErrorCode LaMEMMultigrid_MGCreate(void *mg, void *jr) {
    return MGCreate((MG*)mg, (JacRes*)jr);
}

PetscErrorCode LaMEMMultigrid_MGDestroy(void *mg) {
    return MGDestroy((MG*)mg);
}

PetscErrorCode LaMEMMultigrid_MGSetup(void *mg, Mat A) {
    return MGSetup((MG*)mg, A);
}

PetscErrorCode LaMEMMultigrid_MGSetupCoarse(void *mg, Mat A) {
    return MGSetupCoarse((MG*)mg, A);
}

PetscErrorCode LaMEMMultigrid_MGApply(PC pc, Vec x, Vec y) {
    return MGApply(pc, x, y);
}

PetscErrorCode LaMEMMultigrid_MGDumpMat(void *mg) {
    return MGDumpMat((MG*)mg);
}

PetscErrorCode LaMEMMultigrid_MGGetNumLevels(void *mg) {
    return MGGetNumLevels((MG*)mg);
}

// Helper functions for matrix row setup
void LaMEMMultigrid_GetRowRestrict(PetscBool scale, PetscScalar parent, PetscInt n, 
                                   PetscInt idx[], PetscScalar bc[], PetscScalar v[], 
                                   PetscScalar vs[], PetscScalar eta_fine[], PetscScalar eta_crs) {
    getRowRestrict(scale, parent, n, idx, bc, v, vs, eta_fine, eta_crs);
}

void LaMEMMultigrid_GetRowProlong(PetscBool scale, PetscInt parent, PetscScalar parent_bc, 
                                  PetscInt n, PetscScalar bc[], PetscScalar v[], PetscScalar vs[], 
                                  PetscScalar eta_crs[], PetscScalar eta_fine) {
    getRowProlong(scale, parent, parent_bc, n, bc, v, vs, eta_crs, eta_fine);
}

}