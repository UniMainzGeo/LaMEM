#include "LaMEMAdvect_C.h"
#include "../LaMEM.h"
#include "../advect.h"
#include "../phase.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../surf.h"
#include "../marker.h"
#include "../AVD.h"
#include "../cvi.h"
#include "../subgrid.h"
#include "../tools.h"
#include "../phase_transition.h"
#include "../passive_tracer.h"

extern "C" {

// Core advection functions
PetscErrorCode LaMEMAdvect_Create(void *actx, void *fb) {
    return ADVCreate((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMAdvect_SetType(void *actx, void *fb) {
    return ADVSetType((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMAdvect_Destroy(void *actx) {
    return ADVDestroy((AdvCtx*)actx);
}

// Restart functions
PetscErrorCode LaMEMAdvect_ReadRestart(void *actx, void *fp) {
    return ADVReadRestart((AdvCtx*)actx, (FILE*)fp);
}

PetscErrorCode LaMEMAdvect_WriteRestart(void *actx, void *fp) {
    return ADVWriteRestart((AdvCtx*)actx, (FILE*)fp);
}

// Data management
PetscErrorCode LaMEMAdvect_CreateData(void *actx) {
    return ADVCreateData((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_SetBGPhase(void *actx) {
    return ADVSetBGPhase((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_ReAllocStorage(void *actx, PetscInt capacity) {
    return ADVReAllocStorage((AdvCtx*)actx, capacity);
}

// Main advection operations
PetscErrorCode LaMEMAdvect_Advect(void *actx) {
    return ADVAdvect((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_Remap(void *actx) {
    return ADVRemap((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_Exchange(void *actx) {
    return ADVExchange((AdvCtx*)actx);
}

// History projection and interpolation
PetscErrorCode LaMEMAdvect_ProjHistGridToMark(void *actx) {
    return ADVProjHistGridToMark((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_InterpFieldToMark(void *actx, PetscInt icase) {
    InterpCase case_enum;
    switch(icase) {
        case LAMEM_INTERP_PHASE:     case_enum = _PHASE_; break;
        case LAMEM_INTERP_STRESS:    case_enum = _STRESS_; break;
        case LAMEM_INTERP_APS:       case_enum = _APS_; break;
        case LAMEM_INTERP_ATS:       case_enum = _ATS_; break;
        case LAMEM_INTERP_VORTICITY: case_enum = _VORTICITY_; break;
        case LAMEM_INTERP_DISP:      case_enum = _DISP_; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return ADVInterpFieldToMark((AdvCtx*)actx, case_enum);
}

PetscErrorCode LaMEMAdvect_ProjHistMarkToGrid(void *actx) {
    return ADVProjHistMarkToGrid((AdvCtx*)actx);
}

// Marker movement and positioning
PetscErrorCode LaMEMAdvect_AdvectMark(void *actx) {
    return ADVAdvectMark((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_MapMarkToDomains(void *actx) {
    return ADVMapMarkToDomains((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_MapMarkToCells(void *actx) {
    return ADVMapMarkToCells((AdvCtx*)actx);
}

// MPI communication
PetscErrorCode LaMEMAdvect_ExchangeNumMark(void *actx) {
    return ADVExchangeNumMark((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_CreateMPIBuff(void *actx) {
    return ADVCreateMPIBuff((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_ExchangeMark(void *actx) {
    return ADVExchangeMark((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_DestroyMPIBuff(void *actx) {
    return ADVDestroyMPIBuff((AdvCtx*)actx);
}

// Periodic boundary conditions
PetscErrorCode LaMEMAdvect_ApplyPeriodic(void *actx) {
    return ADVApplyPeriodic((AdvCtx*)actx);
}

// Garbage collection and cleanup
PetscErrorCode LaMEMAdvect_CollectGarbage(void *actx) {
    return ADVCollectGarbage((AdvCtx*)actx);
}

// Interpolation functions
PetscErrorCode LaMEMAdvect_InterpMarkToCell(void *actx) {
    return ADVInterpMarkToCell((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_InterpMarkToEdge(void *actx, PetscInt iphase, PetscInt icase) {
    InterpCase case_enum;
    switch(icase) {
        case LAMEM_INTERP_PHASE:     case_enum = _PHASE_; break;
        case LAMEM_INTERP_STRESS:    case_enum = _STRESS_; break;
        case LAMEM_INTERP_APS:       case_enum = _APS_; break;
        case LAMEM_INTERP_ATS:       case_enum = _ATS_; break;
        case LAMEM_INTERP_VORTICITY: case_enum = _VORTICITY_; break;
        case LAMEM_INTERP_DISP:      case_enum = _DISP_; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return ADVInterpMarkToEdge((AdvCtx*)actx, iphase, case_enum);
}

// Marker control and quality checks
PetscErrorCode LaMEMAdvect_MarkControl(void *actx) {
    return ADVMarkControl((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_CheckCorners(void *actx) {
    return ADVCheckCorners((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_CheckMarkPhases(void *actx) {
    return ADVCheckMarkPhases((AdvCtx*)actx);
}

// Special cases and utilities
PetscErrorCode LaMEMAdvect_UpdateHistADVNone(void *actx) {
    return ADVUpdateHistADVNone((AdvCtx*)actx);
}

PetscErrorCode LaMEMAdvect_SelectTimeStep(void *actx, PetscInt *restart) {
    return ADVSelectTimeStep((AdvCtx*)actx, restart);
}

PetscErrorCode LaMEMAdvect_MarkerAdiabatic(void *actx) {
    return ADVMarkerAdiabatic((AdvCtx*)actx);
}

// Marker manipulation functions
PetscErrorCode LaMEMAdvect_MarkerMerge(void *A, void *B, void *C) {
    return MarkerMerge(*(Marker*)A, *(Marker*)B, *(Marker*)C);
}

// Getter functions for marker information
PetscErrorCode LaMEMAdvect_GetNumMarkers(void *actx, PetscInt *nummark) {
    if (!actx || !nummark) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    *nummark = adv->nummark;
    return 0;
}

PetscErrorCode LaMEMAdvect_GetMarkerCapacity(void *actx, PetscInt *markcap) {
    if (!actx || !markcap) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    *markcap = adv->markcap;
    return 0;
}

PetscErrorCode LaMEMAdvect_GetMarkerData(void *actx, PetscInt idx, PetscInt *phase, PetscScalar X[3], PetscScalar *p, PetscScalar *T) {
    if (!actx) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    if (idx < 0 || idx >= adv->nummark) return PETSC_ERR_ARG_OUTOFRANGE;
    
    Marker *marker = &adv->markers[idx];
    if (phase) *phase = marker->phase;
    if (X) {
        X[0] = marker->X[0];
        X[1] = marker->X[1];
        X[2] = marker->X[2];
    }
    if (p) *p = marker->p;
    if (T) *T = marker->T;
    return 0;
}

// Setter functions for advection parameters
PetscErrorCode LaMEMAdvect_SetAdvectionType(void *actx, PetscInt advect_type) {
    if (!actx) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    
    switch(advect_type) {
        case LAMEM_ADV_NONE:        adv->advect = ADV_NONE; break;
        case LAMEM_BASIC_EULER:     adv->advect = BASIC_EULER; break;
        case LAMEM_EULER:           adv->advect = EULER; break;
        case LAMEM_RUNGE_KUTTA_2:   adv->advect = RUNGE_KUTTA_2; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return 0;
}

PetscErrorCode LaMEMAdvect_SetVelInterpType(void *actx, PetscInt interp_type) {
    if (!actx) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    
    switch(interp_type) {
        case LAMEM_STAG:    adv->interp = STAG; break;
        case LAMEM_MINMOD:  adv->interp = MINMOD; break;
        case LAMEM_STAG_P:  adv->interp = STAG_P; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return 0;
}

PetscErrorCode LaMEMAdvect_SetMarkCtrlType(void *actx, PetscInt ctrl_type) {
    if (!actx) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    
    switch(ctrl_type) {
        case LAMEM_CTRL_NONE:   adv->mctrl = CTRL_NONE; break;
        case LAMEM_CTRL_BASIC:  adv->mctrl = CTRL_BASIC; break;
        case LAMEM_CTRL_AVD:    adv->mctrl = CTRL_AVD; break;
        case LAMEM_CTRL_SUB:    adv->mctrl = CTRL_SUB; break;
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return 0;
}

PetscErrorCode LaMEMAdvect_SetMarkersPerCell(void *actx, PetscInt NumPartX, PetscInt NumPartY, PetscInt NumPartZ) {
    if (!actx) return PETSC_ERR_ARG_NULL;
    AdvCtx *adv = (AdvCtx*)actx;
    
    adv->NumPartX = NumPartX;
    adv->NumPartY = NumPartY;
    adv->NumPartZ = NumPartZ;
    return 0;
}

}