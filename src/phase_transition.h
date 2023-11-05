/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

/*
 *  Created on: Apr 20, 2020
 *      Author: piccolo
 */

//-----------------------------------------------------------------------------

#ifndef phase_transition_h_
#define phase_transition_h_

//-----------------------------------------------------------------------------

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
PetscErrorCode Set_NotInAirBox_Phase_Transition(Ph_trans_t *ph, DBMat *dbm, FB *fb);
PetscErrorCode SetClapeyron_Eq(Ph_trans_t *ph);
PetscErrorCode Overwrite_density(DBMat *dbm);
PetscErrorCode Phase_Transition(AdvCtx *actx);
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1,PetscInt PH2, 
			  Controls ctrl,Scaling *scal, SolVarCell *svCell, PetscInt *ph, PetscScalar *T, PetscInt *InsideAbove, PetscScalar, JacRes *jr, PetscInt cellID);
PetscInt Check_Phase_above_below(PetscInt *phase_array, Marker *P,PetscInt num_phas);
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl, SolVarCell *svCell, PetscInt *ph, PetscInt *InsideAbove, PetscScalar time);
PetscInt Check_Box_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, PetscInt *ph, PetscScalar *T, PetscInt *InsideAbove);      
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl, PetscInt *ph, PetscInt *InsideAbove);

PetscInt Check_NotInAirBox_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, PetscInt *ph_out, PetscScalar *T_out, JacRes *jr, PetscInt cellID);
PetscErrorCode MovingBox(Ph_trans_t *PhaseTrans, TSSol *ts, JacRes *jr);
PetscErrorCode LinkNotInAirBoxes(Ph_trans_t *PhaseTrans, JacRes *jr);

PetscErrorCode DynamicPhTr_Init(JacRes *jr);
PetscErrorCode DynamicPhTr_WriteRestart2(JacRes *jr, FILE *fp);
PetscErrorCode DynamicPhTr_WriteRestart(JacRes *jr, FILE *fp);
PetscErrorCode DynamicPhTr_ReadRestart(JacRes *jr, FILE *fp);
PetscErrorCode DynamicPhTrDestroy(DBMat *dbm);

//-----------------------------------------------------------------------------
#endif
