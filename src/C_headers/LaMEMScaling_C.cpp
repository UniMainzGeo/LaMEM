#include "LaMEMScaling_C.h"
#include "../LaMEM.h"
#include "../scaling.h"
#include "../parsing.h"

extern "C" {

// Main scaling function - thin wrapper around the only function in scaling.cpp
PetscErrorCode LaMEMScaling_ScalingCreate(void *scal, void *fb, PetscBool PrintOutput) {
    return ScalingCreate((Scaling*)scal, (FB*)fb, PrintOutput);
}

}