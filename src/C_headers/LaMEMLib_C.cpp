#include "LaMEMLib_C.h"
#include "../LaMEM.h"
#include "../phase.h"
#include "../dike.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../tools.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../interpolate.h"
#include "../surf.h"
#include "../paraViewOutBin.h"
#include "../paraViewOutSurf.h"
#include "../multigrid.h"
#include "../matrix.h"
#include "../lsolve.h"
#include "../nlsolve.h"
#include "../multigrid.h"
#include "../Tensor.h"
#include "../advect.h"
#include "../marker.h"
#include "../paraViewOutMark.h"
#include "../paraViewOutAVD.h"
#include "../objFunct.h"
#include "../adjoint.h"
#include "../paraViewOutPassiveTracers.h"
#include "../LaMEMLib.h"
#include "../phase_transition.h"
#include "../passive_tracer.h"

extern "C" {

PetscErrorCode LaMEMLib_Main(void *param, PetscLogStage stages[4]) {
    return LaMEMLibMain(param, stages);
}

PetscErrorCode LaMEMLib_Create(void *lm, void *param) {
    return LaMEMLibCreate((LaMEMLib*)lm, param);
}

PetscErrorCode LaMEMLib_Solve(void *lm, void *param, PetscLogStage stages[4]) {
    return LaMEMLibSolve((LaMEMLib*)lm, param, stages);
}

PetscErrorCode LaMEMLib_SaveRestart(void *lm) {
    return LaMEMLibSaveRestart((LaMEMLib*)lm);
}

PetscErrorCode LaMEMLib_LoadRestart(void *lm) {
    return LaMEMLibLoadRestart((LaMEMLib*)lm);
}

PetscErrorCode LaMEMLib_SaveGrid(void *lm) {
    return LaMEMLibSaveGrid((LaMEMLib*)lm);
}

PetscErrorCode LaMEMLib_DryRun(void *lm) {
    return LaMEMLibDryRun((LaMEMLib*)lm);
}

PetscErrorCode LaMEMLib_Destroy(void *lm) {
    return LaMEMLibDestroy((LaMEMLib*)lm);
}

}