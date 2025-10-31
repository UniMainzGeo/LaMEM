#include "LaMEMParsing_C.h"
#include "../LaMEM.h"
#include "../parsing.h"
#include "../tools.h"


extern "C" {

// File buffer functions
PetscErrorCode LaMEMParsing_FBLoad(void **pfb, PetscBool DisplayOutput, char *restartFileName) {
    return FBLoad((FB**)pfb, DisplayOutput, restartFileName);
}

PetscErrorCode LaMEMParsing_FBDestroy(void **pfb) {
    return FBDestroy((FB**)pfb);
}

PetscErrorCode LaMEMParsing_FBParseBuffer(void *fb) {
    return FBParseBuffer((FB*)fb);
}

PetscErrorCode LaMEMParsing_FBFindBlocks(void *fb, PetscInt ptype, const char *keybeg, const char *keyend) {
    return FBFindBlocks((FB*)fb, (ParamType)ptype, keybeg, keyend);
}

PetscErrorCode LaMEMParsing_FBFreeBlocks(void *fb) {
    return FBFreeBlocks((FB*)fb);
}

char** LaMEMParsing_FBGetLineRanges(void *fb, PetscInt *lnbeg, PetscInt *lnend) {
    return FBGetLineRanges((FB*)fb, lnbeg, lnend);
}

// Array parameter reading functions
PetscErrorCode LaMEMParsing_FBGetIntArray(void *fb, const char *key, PetscInt *nvalues, 
                                          PetscInt *values, PetscInt num, PetscBool *found) {
    return FBGetIntArray((FB*)fb, key, nvalues, values, num, found);
}

PetscErrorCode LaMEMParsing_FBGetScalarArray(void *fb, const char *key, PetscInt *nvalues, 
                                             PetscScalar *values, PetscInt num, PetscBool *found) {
    return FBGetScalarArray((FB*)fb, key, nvalues, values, num, found);
}

PetscErrorCode LaMEMParsing_FBGetString(void *fb, const char *key, char *str, PetscBool *found) {
    return FBGetString((FB*)fb, key, str, found);
}

// Wrapper parameter reading functions
PetscErrorCode LaMEMParsing_getIntParam(void *fb, PetscInt ptype, const char *key, 
                                        PetscInt *val, PetscInt num, PetscInt maxval) {
    return getIntParam((FB*)fb, (ParamType)ptype, key, val, num, maxval);
}

PetscErrorCode LaMEMParsing_getScalarParam(void *fb, PetscInt ptype, const char *key, 
                                           PetscScalar *val, PetscInt num, PetscScalar scal) {
    return getScalarParam((FB*)fb, (ParamType)ptype, key, val, num, scal);
}

PetscErrorCode LaMEMParsing_getStringParam(void *fb, PetscInt ptype, const char *key, 
                                           char *str, const char *_default) {
    return getStringParam((FB*)fb, (ParamType)ptype, key, str, _default);
}

// PETSc options functions
PetscErrorCode LaMEMParsing_PetscOptionsReadFromFile(void *fb, PetscBool DisplayOutput) {
    return PetscOptionsReadFromFile((FB*)fb, DisplayOutput);
}

PetscErrorCode LaMEMParsing_PetscOptionsReadRestart(FILE *fp) {
    return PetscOptionsReadRestart(fp);
}

PetscErrorCode LaMEMParsing_PetscOptionsWriteRestart(FILE *fp) {
    return PetscOptionsWriteRestart(fp);
}

PetscErrorCode LaMEMParsing_PetscOptionsGetCheckString(const char key[], char str[], PetscBool *set) {
    return PetscOptionsGetCheckString(key, str, set);
}

// Default solver options
PetscErrorCode LaMEMParsing_StokesSetDefaultSolverOptions(void *fb) {
    return StokesSetDefaultSolverOptions((FB*)fb);
}

}