#include "LaMEMPhase_C.h"
#include "../LaMEM.h"
#include "../phase.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../objFunct.h"
#include "../JacRes.h"
#include "../phase_transition.h"

extern "C" {

// Main material database functions
PetscErrorCode LaMEMPhase_DBMatCreate(void *dbm, void *fb, PetscBool PrintOutput) {
    return DBMatCreate((DBMat*)dbm, (FB*)fb, PrintOutput);
}

PetscErrorCode LaMEMPhase_DBMatReadSoft(void *dbm, void *fb, PetscBool PrintOutput) {
    return DBMatReadSoft((DBMat*)dbm, (FB*)fb, PrintOutput);
}

PetscErrorCode LaMEMPhase_DBMatReadPhase(void *dbm, void *fb, PetscBool PrintOutput) {
    return DBMatReadPhase((DBMat*)dbm, (FB*)fb, PrintOutput);
}

PetscErrorCode LaMEMPhase_DBMatOverwriteWithGlobalVariables(void *dbm, void *fb) {
    return DBMatOverwriteWithGlobalVariables((DBMat*)dbm, (FB*)fb);
}

// Material property printing and utilities
void LaMEMPhase_MatPrintScalParam(PetscScalar par, const char key[], const char label[],
                                  void *scal, const char title[], PetscInt *print_title) {
    MatPrintScalParam(par, key, label, (Scaling*)scal, title, print_title);
}

PetscErrorCode LaMEMPhase_PrintMatProp(void *MatProp) {
    return PrintMatProp((Material_t*)MatProp);
}

// Predefined rheological profile functions
PetscErrorCode LaMEMPhase_GetProfileName(void *fb, void *scal, char name[], const char key[]) {
    return GetProfileName((FB*)fb, (Scaling*)scal, name, key);
}

PetscErrorCode LaMEMPhase_SetDiffProfile(void *m, char name[]) {
    return SetDiffProfile((Material_t*)m, name);
}

PetscErrorCode LaMEMPhase_SetDislProfile(void *m, char name[]) {
    return SetDislProfile((Material_t*)m, name);
}

PetscErrorCode LaMEMPhase_SetPeirProfile(void *m, char name[]) {
    return SetPeirProfile((Material_t*)m, name);
}

// Experimental correction functions
PetscErrorCode LaMEMPhase_CorrExpPreFactor(PetscScalar *B, PetscScalar n, PetscInt type, PetscInt MPa) {
    return CorrExpPreFactor(*B, n, (ExpType)type, MPa);
}

PetscErrorCode LaMEMPhase_CorrExpStressStrainRate(PetscScalar *D, PetscScalar *S, PetscInt type, PetscInt MPa) {
    return CorrExpStressStrainRate(*D, *S, (ExpType)type, MPa);
}

}