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

PetscErrorCode Phase_Transition(ConstEqCtx  *ctx,PetscInt ph);
PetscErrorCode Transition(Ph_trans_t *PhaseTrans, ConstEqCtx  *ctx, PetscInt id);


#endif /* PHASE_TRANSITION_H_ */
