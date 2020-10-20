/*
 * phase_transition.h
 *
 *  Created on: Apr 20, 2020
 *      Author: piccolo
 */

#ifndef phase_transition_h_
#define phase_transition_h_
struct Ph_trans_t;
struct JacRes;
struct ConstEqCtx;
struct Marker;


// read phase transition law
PetscErrorCode DBMatReadPhaseTr(DBMat *dbm, FB *fb);
PetscErrorCode Set_Constant_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb,PetscInt ID);
PetscErrorCode Set_Clapeyron_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb, PetscInt ID);
PetscErrorCode SetClapeyron_Eq(Ph_trans_t *ph);
PetscErrorCode Overwrite_density(DBMat *dbm);
PetscErrorCode Phase_Transition(AdvCtx *actx);
PetscErrorCode Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1,PetscInt PH2);
PetscInt Check_Phase_above_below(PetscInt *phase_array, Marker *P,PetscInt num_phas);
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2);
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2);// softening parameter

#endif /* PHASE_TRANSITION_H_ */
