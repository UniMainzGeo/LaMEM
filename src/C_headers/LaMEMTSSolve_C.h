#ifndef LAMEMTSSOLVE_C_H
#define LAMEMTSSOLVE_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Main time stepping control functions
PetscErrorCode LaMEMTSSolve_TSSolCreate(void *ts, void *fb);
PetscInt LaMEMTSSolve_TSSolIsDone(void *ts);
PetscErrorCode LaMEMTSSolve_TSSolStepForward(void *ts);
PetscInt LaMEMTSSolve_TSSolIsRestart(void *ts);
PetscInt LaMEMTSSolve_TSSolIsOutput(void *ts);

// CFL and time step control
PetscErrorCode LaMEMTSSolve_TSSolGetCFLStep(void *ts, PetscScalar gidtmax, PetscInt *restart);

// Time step scheduling functions
PetscErrorCode LaMEMTSSolve_TSSolGetPeriodSteps(PetscScalar dt_start, PetscScalar dt_end, 
                                                PetscScalar span, PetscScalar *dt, PetscInt *n);
PetscErrorCode LaMEMTSSolve_TSSolMakeSchedule(void *ts);
PetscErrorCode LaMEMTSSolve_TSSolAdjustSchedule(void *ts, PetscScalar dt_cfl, 
                                                PetscInt istep, PetscScalar *schedule);

#ifdef __cplusplus
}
#endif

#endif