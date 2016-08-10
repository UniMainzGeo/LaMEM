/*
 * nlsolveExplicit.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Beatriz
 */

#ifndef __nlsolveExplicit_h__
#define __nlsolveExplicit_h__

PetscErrorCode GetVelocities(JacRes *jr);
//PetscErrorCode GetPressure2(JacRes *jr);
//PetscErrorCode GetStress(JacRes *jr);
//PetscErrorCode JacResGetMomentumResidual2(JacRes *jr);
//PetscErrorCode SolveEquationsWave(JacRes *jr);
PetscErrorCode FormMomentumResidualPressureAndVelocities(JacRes *jr);
//PetscErrorCode GetPressure(JacRes *jr);
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user);
PetscErrorCode CheckElasticProperties(JacRes *jr, UserCtx *user);
PetscErrorCode ShowValues(JacRes *jr, PetscInt n); //, UserCtx *user);
PetscErrorCode SaveVelocitiesForSeismicStation(JacRes *jr, UserCtx *user);
//PetscErrorCode PutSeismicSource(JacRes *jr, AdvCtx *actx, UserCtx *user);
PetscErrorCode PrintStress(JacRes *jr);
PetscErrorCode UpdateHistoryFields(JacRes *jr);
PetscErrorCode ModifyStress(JacRes *jr);

#endif
