#ifndef LAMEMSCALING_C_H
#define LAMEMSCALING_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Unit system constants
#define LAMEM_UNITS_NONE  0
#define LAMEM_UNITS_SI    1  
#define LAMEM_UNITS_GEO   2

// Main scaling function - only one that exists in scaling.cpp
PetscErrorCode LaMEMScaling_ScalingCreate(void *scal, void *fb, PetscBool PrintOutput);

#ifdef __cplusplus
}
#endif

#endif