#ifndef LAMEMPASSIVETRACER_C_H
#define LAMEMPASSIVETRACER_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Passive tracer condition constants
#define LAMEM_ALWAYS     0
#define LAMEM_MELT_FR    1
#define LAMEM_TEMP_PTR   2
#define LAMEM_TIME_PTR   3
#define LAMEM_PRES_PTR   4

// Main passive tracer functions
PetscErrorCode LaMEMPassiveTracer_ADVPtrPassive_Tracer_create(void *actx, void *fb);
PetscErrorCode LaMEMPassiveTracer_ADVPtrReCreateStorage(void *actx);
PetscErrorCode LaMEMPassiveTracer_ADVPassiveTracerInit(void *actx);
PetscErrorCode LaMEMPassiveTracer_ADVPtrInitCoord(void *actx);
PetscErrorCode LaMEMPassiveTracer_ADV_Assign_Phase(void *actx);

// Advection functions
PetscErrorCode LaMEMPassiveTracer_ADVAdvectPassiveTracer(void *actx);
PetscErrorCode LaMEMPassiveTracer_ADVMarkCrossFreeSurfPassive_Tracers(void *actx);
PetscErrorCode LaMEMPassiveTracer_Check_advection_condition(void *actx, PetscInt jj, PetscInt ID, 
                                                           PetscScalar xp, PetscScalar yp, PetscScalar zp, 
                                                           PetscScalar P, PetscScalar T, PetscScalar mf);

// Cleanup and I/O functions
PetscErrorCode LaMEMPassiveTracer_ADVPtrDestroy(void *actx);
PetscErrorCode LaMEMPassiveTracer_Passive_Tracer_WriteRestart(void *actx, FILE *fp);
PetscErrorCode LaMEMPassiveTracer_ReadPassive_Tracers(void *actx, FILE *fp);

// Utility functions
PetscErrorCode LaMEMPassiveTracer_Sync_Vector(Vec x, void *actx, PetscInt nummark);

#ifdef __cplusplus
}
#endif

#endif