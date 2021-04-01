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
struct Material_t;
struct Freesurf;

// read phase transition law
PetscErrorCode DBMatReadPhaseTr(DBMat *dbm, FB *fb);
PetscErrorCode Set_Constant_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb);
PetscErrorCode Set_Clapeyron_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb);
PetscErrorCode Set_Box_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb);
PetscErrorCode SetClapeyron_Eq(Ph_trans_t *ph);
PetscErrorCode Overwrite_density(DBMat *dbm);
PetscErrorCode Phase_Transition(AdvCtx *actx);
PetscErrorCode Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1,PetscInt PH2, 
			  Controls ctrl,Scaling *scal, SolVarCell *svCell, PetscInt *ph, PetscScalar *T, JacRes *jr);  //, PetscInt *InsideAbove, PetscScalar); CHANGED BY JS
PetscInt Check_Phase_above_below(PetscInt *phase_array, Marker *P,PetscInt num_phas);
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl, SolVarCell *svCell); //, PetscInt *ph, PetscInt *InsideAbove, PetscScalar time);
PetscInt Check_Box_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, PetscInt *ph, PetscScalar *T); //, PetscInt *InsideAbove);
PetscInt Check_DikeBox_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, PetscInt *ph_out, PetscScalar *T_out, JacRes *jr);
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl);//, PetscInt *ph);//, PetscInt *InsideAbove); 

// edited the Check_Constant_Phase_Transition() for *ctrl instead of ctrl, same for Check_Clapeyron_Phase_Transition()
//PetscErrorCode Check_AirPhaseRatio_Box_Transition(Material_t *mat, Controls ctrl, SolVarCell *svCell);  // NEW FOR DIKE, air check inside dike phase transition box, maybe ConstEqCtx *ctx not necessary??

#endif /* PHASE_TRANSITION_H_ */
