#ifndef LAMEMMAIN_C_H
#define LAMEMMAIN_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrapper for LaMEM's main logic
PetscErrorCode LaMEMMain(int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif