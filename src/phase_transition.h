/*
 * phase_transition.h
 *
 *  Created on: Apr 20, 2020
 *      Author: piccolo
 */

#ifndef phase_transition_h_
#define phase_transition_h_
struct Marker;
struct Ph_trans_t;
struct JacRes;

PetscErrorCode PhTr_assign_primPh(AdvCtx *actx);
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt id, PetscInt PH);
PetscErrorCode Phase_Transition(AdvCtx *actx);


#endif /* PHASE_TRANSITION_H_ */
