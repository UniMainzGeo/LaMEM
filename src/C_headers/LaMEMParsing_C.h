#ifndef LAMEMPARSING_C_H
#define LAMEMPARSING_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// File buffer functions
PetscErrorCode LaMEMParsing_FBLoad(void **pfb, PetscBool DisplayOutput, char *restartFileName);
PetscErrorCode LaMEMParsing_FBDestroy(void **pfb);
PetscErrorCode LaMEMParsing_FBParseBuffer(void *fb);
PetscErrorCode LaMEMParsing_FBFindBlocks(void *fb, PetscInt ptype, const char *keybeg, const char *keyend);
PetscErrorCode LaMEMParsing_FBFreeBlocks(void *fb);
char** LaMEMParsing_FBGetLineRanges(void *fb, PetscInt *lnbeg, PetscInt *lnend);

// Array parameter reading functions
PetscErrorCode LaMEMParsing_FBGetIntArray(void *fb, const char *key, PetscInt *nvalues, 
                                          PetscInt *values, PetscInt num, PetscBool *found);
PetscErrorCode LaMEMParsing_FBGetScalarArray(void *fb, const char *key, PetscInt *nvalues, 
                                             PetscScalar *values, PetscInt num, PetscBool *found);
PetscErrorCode LaMEMParsing_FBGetString(void *fb, const char *key, char *str, PetscBool *found);

// Wrapper parameter reading functions
PetscErrorCode LaMEMParsing_getIntParam(void *fb, PetscInt ptype, const char *key, 
                                        PetscInt *val, PetscInt num, PetscInt maxval);
PetscErrorCode LaMEMParsing_getScalarParam(void *fb, PetscInt ptype, const char *key, 
                                           PetscScalar *val, PetscInt num, PetscScalar scal);
PetscErrorCode LaMEMParsing_getStringParam(void *fb, PetscInt ptype, const char *key, 
                                           char *str, const char *_default);

// PETSc options functions
PetscErrorCode LaMEMParsing_PetscOptionsReadFromFile(void *fb, PetscBool DisplayOutput);
PetscErrorCode LaMEMParsing_PetscOptionsReadRestart(FILE *fp);
PetscErrorCode LaMEMParsing_PetscOptionsWriteRestart(FILE *fp);
PetscErrorCode LaMEMParsing_PetscOptionsGetCheckString(const char key[], char str[], PetscBool *set);

// Default solver options
PetscErrorCode LaMEMParsing_StokesSetDefaultSolverOptions(void *fb);

// Constants for parsing modes
#define LAMEM_REQUIRED 1
#define LAMEM_OPTIONAL 2

#ifdef __cplusplus
}
#endif

#endif