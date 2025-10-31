#include "LaMEMPhaseTransition_C.h"
#include "../LaMEM.h"
#include "../AVD.h"
#include "../advect.h"
#include "../scaling.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../tools.h"
#include "../phase_transition.h"
#include "../phase.h"
#include "../constEq.h"
#include "../parsing.h"
#include "../objFunct.h"
#include "../surf.h"
#include "../tssolve.h"
#include "../dike.h"

extern "C" {

// Main phase transition database functions
PetscErrorCode LaMEMPhaseTransition_DBMatReadPhaseTr(void *dbm, void *fb) {
    return DBMatReadPhaseTr((DBMat*)dbm, (FB*)fb);
}

PetscErrorCode LaMEMPhaseTransition_Set_Constant_Phase_Transition(void *ph, void *dbm, void *fb) {
    return Set_Constant_Phase_Transition((Ph_trans_t*)ph, (DBMat*)dbm, (FB*)fb);
}

PetscErrorCode LaMEMPhaseTransition_Set_Box_Phase_Transition(void *ph, void *dbm, void *fb) {
    return Set_Box_Phase_Transition((Ph_trans_t*)ph, (DBMat*)dbm, (FB*)fb);
}

PetscErrorCode LaMEMPhaseTransition_Set_NotInAirBox_Phase_Transition(void *ph, void *dbm, void *fb) {
    return Set_NotInAirBox_Phase_Transition((Ph_trans_t*)ph, (DBMat*)dbm, (FB*)fb);
}

PetscErrorCode LaMEMPhaseTransition_Set_Clapeyron_Phase_Transition(void *ph, void *dbm, void *fb) {
    return Set_Clapeyron_Phase_Transition((Ph_trans_t*)ph, (DBMat*)dbm, (FB*)fb);
}

PetscErrorCode LaMEMPhaseTransition_SetClapeyron_Eq(void *ph) {
    return SetClapeyron_Eq((Ph_trans_t*)ph);
}

PetscErrorCode LaMEMPhaseTransition_Overwrite_density(void *dbm) {
    return Overwrite_density((DBMat*)dbm);
}

// Dynamic phase transition initialization and management
PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_Init(void *jr) {
    return DynamicPhTr_Init((JacRes*)jr);
}

PetscErrorCode LaMEMPhaseTransition_MovingBox(void *PhaseTrans, void *ts, void *jr) {
    return MovingBox((Ph_trans_t*)PhaseTrans, (TSSol*)ts, (JacRes*)jr);
}

PetscErrorCode LaMEMPhaseTransition_LinkNotInAirBoxes(void *PhaseTrans, void *jr) {
    return LinkNotInAirBoxes((Ph_trans_t*)PhaseTrans, (JacRes*)jr);
}

// Main phase transition execution
PetscErrorCode LaMEMPhaseTransition_Phase_Transition(void *actx) {
    return Phase_Transition((AdvCtx*)actx);
}

// Individual phase transition check functions
PetscInt LaMEMPhaseTransition_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                         void *ctrl, void *scal, void *svCell, PetscInt *ph_out, 
                                         PetscScalar *T_out, PetscInt *InsideAbove, PetscScalar time, 
                                         void *jr, PetscInt cellID) {
    return Transition((Ph_trans_t*)PhaseTrans, (Marker*)P, PH1, PH2, *(Controls*)ctrl, (Scaling*)scal, 
                      (SolVarCell*)svCell, ph_out, T_out, InsideAbove, time, (JacRes*)jr, cellID);
}

PetscInt LaMEMPhaseTransition_Check_Constant_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                              void *ctrl, void *svCell, PetscInt *ph_out, 
                                                              PetscInt *InAbove, PetscScalar time) {
    return Check_Constant_Phase_Transition((Ph_trans_t*)PhaseTrans, (Marker*)P, PH1, PH2, *(Controls*)ctrl, 
                                            (SolVarCell*)svCell, ph_out, InAbove, time);
}

PetscInt LaMEMPhaseTransition_Check_Box_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                         void *scal, PetscInt *ph_out, PetscScalar *T_out, 
                                                         PetscInt *InAbove) {
    return Check_Box_Phase_Transition((Ph_trans_t*)PhaseTrans, (Marker*)P, PH1, PH2, (Scaling*)scal, 
                                       ph_out, T_out, InAbove);
}

PetscInt LaMEMPhaseTransition_Check_NotInAirBox_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                                 void *scal, PetscInt *ph_out, PetscScalar *T_out, 
                                                                 void *jr, PetscInt cellID) {
    return Check_NotInAirBox_Phase_Transition((Ph_trans_t*)PhaseTrans, (Marker*)P, PH1, PH2, (Scaling*)scal, 
                                               ph_out, T_out, (JacRes*)jr, cellID);
}

PetscInt LaMEMPhaseTransition_Check_Clapeyron_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                               void *ctrl, PetscInt *ph_out, PetscInt *InAbove) {
    return Check_Clapeyron_Phase_Transition((Ph_trans_t*)PhaseTrans, (Marker*)P, PH1, PH2, *(Controls*)ctrl, 
                                             ph_out, InAbove);
}

// Utility functions
PetscInt LaMEMPhaseTransition_Check_Phase_above_below(PetscInt *phase_array, void *P, PetscInt num_phas) {
    return Check_Phase_above_below(phase_array, (Marker*)P, num_phas);
}

// Restart I/O functions
PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_WriteRestart(void *jr, FILE *fp) {
    return DynamicPhTr_WriteRestart((JacRes*)jr, fp);
}

PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_ReadRestart(void *jr, FILE *fp) {
    return DynamicPhTr_ReadRestart((JacRes*)jr, fp);
}

PetscErrorCode LaMEMPhaseTransition_DynamicPhTrDestroy(void *dbm) {
    return DynamicPhTrDestroy((DBMat*)dbm);
}

}