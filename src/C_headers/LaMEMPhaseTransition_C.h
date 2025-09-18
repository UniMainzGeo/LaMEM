#ifndef LAMEMPHASETRANSITION_C_H
#define LAMEMPHASETRANSITION_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Phase transition type constants
#define LAMEM_CONSTANT       0
#define LAMEM_CLAPEYRON      1
#define LAMEM_BOX            2
#define LAMEM_NOTINAIRBOX    3

// Phase transition parameter constants
#define LAMEM_T_             0
#define LAMEM_PRESSURE_      1
#define LAMEM_DEPTH_         2
#define LAMEM_X_             3
#define LAMEM_Y_             4
#define LAMEM_PLASTICSTRAIN_ 5
#define LAMEM_MELTFRACTION_  6
#define LAMEM_TIME_          7

// Main phase transition database functions
PetscErrorCode LaMEMPhaseTransition_DBMatReadPhaseTr(void *dbm, void *fb);
PetscErrorCode LaMEMPhaseTransition_Set_Constant_Phase_Transition(void *ph, void *dbm, void *fb);
PetscErrorCode LaMEMPhaseTransition_Set_Box_Phase_Transition(void *ph, void *dbm, void *fb);
PetscErrorCode LaMEMPhaseTransition_Set_NotInAirBox_Phase_Transition(void *ph, void *dbm, void *fb);
PetscErrorCode LaMEMPhaseTransition_Set_Clapeyron_Phase_Transition(void *ph, void *dbm, void *fb);
PetscErrorCode LaMEMPhaseTransition_SetClapeyron_Eq(void *ph);
PetscErrorCode LaMEMPhaseTransition_Overwrite_density(void *dbm);

// Dynamic phase transition initialization and management
PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_Init(void *jr);
PetscErrorCode LaMEMPhaseTransition_MovingBox(void *PhaseTrans, void *ts, void *jr);
PetscErrorCode LaMEMPhaseTransition_LinkNotInAirBoxes(void *PhaseTrans, void *jr);

// Main phase transition execution
PetscErrorCode LaMEMPhaseTransition_Phase_Transition(void *actx);

// Individual phase transition check functions
PetscInt LaMEMPhaseTransition_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                         void *ctrl, void *scal, void *svCell, PetscInt *ph_out, 
                                         PetscScalar *T_out, PetscInt *InsideAbove, PetscScalar time, 
                                         void *jr, PetscInt cellID);
PetscInt LaMEMPhaseTransition_Check_Constant_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                              void *ctrl, void *svCell, PetscInt *ph_out, 
                                                              PetscInt *InAbove, PetscScalar time);
PetscInt LaMEMPhaseTransition_Check_Box_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                         void *scal, PetscInt *ph_out, PetscScalar *T_out, 
                                                         PetscInt *InAbove);
PetscInt LaMEMPhaseTransition_Check_NotInAirBox_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                                 void *scal, PetscInt *ph_out, PetscScalar *T_out, 
                                                                 void *jr, PetscInt cellID);
PetscInt LaMEMPhaseTransition_Check_Clapeyron_Phase_Transition(void *PhaseTrans, void *P, PetscInt PH1, PetscInt PH2, 
                                                               void *ctrl, PetscInt *ph_out, PetscInt *InAbove);

// Utility functions
PetscInt LaMEMPhaseTransition_Check_Phase_above_below(PetscInt *phase_array, void *P, PetscInt num_phas);

// Restart I/O functions
PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_WriteRestart(void *jr, FILE *fp);
PetscErrorCode LaMEMPhaseTransition_DynamicPhTr_ReadRestart(void *jr, FILE *fp);
PetscErrorCode LaMEMPhaseTransition_DynamicPhTrDestroy(void *dbm);

#ifdef __cplusplus
}
#endif

#endif