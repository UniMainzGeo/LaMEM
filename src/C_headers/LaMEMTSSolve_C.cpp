#include "LaMEMTSSolve_C.h"
#include "../LaMEM.h"
#include "../tssolve.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tools.h"

extern "C" {

// Main time stepping control functions
PetscErrorCode LaMEMTSSolve_TSSolCreate(void *ts, void *fb) {
    return TSSolCreate((TSSol*)ts, (FB*)fb);
}

PetscInt LaMEMTSSolve_TSSolIsDone(void *ts) {
    return TSSolIsDone((TSSol*)ts);
}

PetscErrorCode LaMEMTSSolve_TSSolStepForward(void *ts) {
    return TSSolStepForward((TSSol*)ts);
}

PetscInt LaMEMTSSolve_TSSolIsRestart(void *ts) {
    return TSSolIsRestart((TSSol*)ts);
}

PetscInt LaMEMTSSolve_TSSolIsOutput(void *ts) {
    return TSSolIsOutput((TSSol*)ts);
}

// CFL and time step control
PetscErrorCode LaMEMTSSolve_TSSolGetCFLStep(void *ts, PetscScalar gidtmax, PetscInt *restart) {
    return TSSolGetCFLStep((TSSol*)ts, gidtmax, restart);
}

// Time step scheduling functions
PetscErrorCode LaMEMTSSolve_TSSolGetPeriodSteps(PetscScalar dt_start, PetscScalar dt_end, 
                                                PetscScalar span, PetscScalar *dt, PetscInt *n) {
    return TSSolGetPeriodSteps(dt_start, dt_end, span, dt, *n);
}

PetscErrorCode LaMEMTSSolve_TSSolMakeSchedule(void *ts) {
    return TSSolMakeSchedule((TSSol*)ts);
}

PetscErrorCode LaMEMTSSolve_TSSolAdjustSchedule(void *ts, PetscScalar dt_cfl, 
                                                PetscInt istep, PetscScalar *schedule) {
    return TSSolAdjustSchedule((TSSol*)ts, dt_cfl, istep, schedule);
}

}