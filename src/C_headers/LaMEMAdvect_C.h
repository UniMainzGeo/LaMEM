#ifndef LAMEMADVECT_C_H
#define LAMEMADVECT_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Advection types
#define LAMEM_ADV_NONE          0  // no advection (grid-based solution)
#define LAMEM_BASIC_EULER       1  // basic Euler implementation
#define LAMEM_EULER             2  // Euler explicit in time
#define LAMEM_RUNGE_KUTTA_2     3  // Runge-Kutta 2nd order in space

// Velocity interpolation types
#define LAMEM_STAG              0  // trilinear interpolation from FDSTAG points
#define LAMEM_MINMOD            1  // MINMOD interpolation
#define LAMEM_STAG_P            2  // empirical approach

// Marker control types
#define LAMEM_CTRL_NONE         0  // no marker control
#define LAMEM_CTRL_BASIC        1  // AVD for cells + corner insertion
#define LAMEM_CTRL_AVD          2  // pure AVD for all control volumes
#define LAMEM_CTRL_SUB          3  // simple marker control method

// Setup types
#define LAMEM_GEOM              0  // read geometric primitives from input file
#define LAMEM_FILES             1  // read coordinates, phase and temperature from files
#define LAMEM_POLYGONS          2  // read polygons from file

// Interpolation cases
#define LAMEM_INTERP_PHASE      0  // phase ratio
#define LAMEM_INTERP_STRESS     1  // deviatoric stress
#define LAMEM_INTERP_APS        2  // accumulated plastic strain
#define LAMEM_INTERP_ATS        3  // accumulated total strain
#define LAMEM_INTERP_VORTICITY  4  // vorticity pseudo-vector components
#define LAMEM_INTERP_DISP       5  // displacement

// Core advection functions
PetscErrorCode LaMEMAdvect_Create(void *actx, void *fb);
PetscErrorCode LaMEMAdvect_SetType(void *actx, void *fb);
PetscErrorCode LaMEMAdvect_Destroy(void *actx);

// Restart functions
PetscErrorCode LaMEMAdvect_ReadRestart(void *actx, void *fp);
PetscErrorCode LaMEMAdvect_WriteRestart(void *actx, void *fp);

// Data management
PetscErrorCode LaMEMAdvect_CreateData(void *actx);
PetscErrorCode LaMEMAdvect_SetBGPhase(void *actx);
PetscErrorCode LaMEMAdvect_ReAllocStorage(void *actx, PetscInt capacity);

// Main advection operations
PetscErrorCode LaMEMAdvect_Advect(void *actx);
PetscErrorCode LaMEMAdvect_Remap(void *actx);
PetscErrorCode LaMEMAdvect_Exchange(void *actx);

// History projection and interpolation
PetscErrorCode LaMEMAdvect_ProjHistGridToMark(void *actx);
PetscErrorCode LaMEMAdvect_InterpFieldToMark(void *actx, PetscInt icase);
PetscErrorCode LaMEMAdvect_ProjHistMarkToGrid(void *actx);

// Marker movement and positioning
PetscErrorCode LaMEMAdvect_AdvectMark(void *actx);
PetscErrorCode LaMEMAdvect_MapMarkToDomains(void *actx);
PetscErrorCode LaMEMAdvect_MapMarkToCells(void *actx);

// MPI communication
PetscErrorCode LaMEMAdvect_ExchangeNumMark(void *actx);
PetscErrorCode LaMEMAdvect_CreateMPIBuff(void *actx);
PetscErrorCode LaMEMAdvect_ExchangeMark(void *actx);
PetscErrorCode LaMEMAdvect_DestroyMPIBuff(void *actx);

// Periodic boundary conditions
PetscErrorCode LaMEMAdvect_ApplyPeriodic(void *actx);

// Garbage collection and cleanup
PetscErrorCode LaMEMAdvect_CollectGarbage(void *actx);

// Interpolation functions
PetscErrorCode LaMEMAdvect_InterpMarkToCell(void *actx);
PetscErrorCode LaMEMAdvect_InterpMarkToEdge(void *actx, PetscInt iphase, PetscInt icase);

// Marker control and quality checks
PetscErrorCode LaMEMAdvect_MarkControl(void *actx);
PetscErrorCode LaMEMAdvect_CheckCorners(void *actx);
PetscErrorCode LaMEMAdvect_CheckMarkPhases(void *actx);

// Special cases and utilities
PetscErrorCode LaMEMAdvect_UpdateHistADVNone(void *actx);
PetscErrorCode LaMEMAdvect_SelectTimeStep(void *actx, PetscInt *restart);
PetscErrorCode LaMEMAdvect_MarkerAdiabatic(void *actx);

// Marker manipulation functions
PetscErrorCode LaMEMAdvect_MarkerMerge(void *A, void *B, void *C);

// Getter functions for marker information
PetscErrorCode LaMEMAdvect_GetNumMarkers(void *actx, PetscInt *nummark);
PetscErrorCode LaMEMAdvect_GetMarkerCapacity(void *actx, PetscInt *markcap);
PetscErrorCode LaMEMAdvect_GetMarkerData(void *actx, PetscInt idx, PetscInt *phase, PetscScalar X[3], PetscScalar *p, PetscScalar *T);

// Setter functions for advection parameters
PetscErrorCode LaMEMAdvect_SetAdvectionType(void *actx, PetscInt advect_type);
PetscErrorCode LaMEMAdvect_SetVelInterpType(void *actx, PetscInt interp_type);
PetscErrorCode LaMEMAdvect_SetMarkCtrlType(void *actx, PetscInt ctrl_type);
PetscErrorCode LaMEMAdvect_SetMarkersPerCell(void *actx, PetscInt NumPartX, PetscInt NumPartY, PetscInt NumPartZ);

#ifdef __cplusplus
}
#endif

#endif