#include "LaMEMMain_C.h"
#include "../LaMEM.h"
#include "../scaling.h"
#include "../objFunct.h"
#include "../parsing.h"
#include "../adjoint.h"
#include "../phase.h"

extern "C" {

// Thin wrapper for LaMEM's main logic
PetscErrorCode LaMEMMain_main(int argc, char **argv)
{
    // Call the main function in LaMEM.cpp
    return LaMEM_main(argc, argv);
}

}