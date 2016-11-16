/*
 * nlsolveExplicit.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Beatriz
 */

#ifndef __nlsolveExplicit_h__
#define __nlsolveExplicit_h__

PetscErrorCode GetVelocities(JacRes *jr, UserCtx *user);
PetscErrorCode FormMomentumResidualPressureAndVelocities(JacRes *jr, UserCtx *user);
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user);
PetscErrorCode ChangeTimeStep(JacRes *jr, UserCtx *user);
PetscErrorCode CheckElasticProperties(JacRes *jr, UserCtx *user);
PetscErrorCode SaveVelocitiesForSeismicStation(JacRes *jr, UserCtx *user);
//PetscErrorCode PutSeismicSource(JacRes *jr, AdvCtx *actx, UserCtx *user);
PetscErrorCode PrintStress(JacRes *jr);
//PetscErrorCode UpdateHistoryFieldsAndGetAxialStressStrain(JacRes *jr, PetscScalar *axial_stress, PetscScalar *axial_strain);
//PetscErrorCode ModifyStress(JacRes *jr);
PetscErrorCode GetCellCoordinatesSource(JacRes *jr);
PetscErrorCode GetStressFromSource(JacRes *jr, UserCtx *User, PetscInt i, PetscInt j, PetscInt k, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz);
//PetscErrorCode FormMomentumResidualPressureAndVelocities(JacRes *jr);
//PetscErrorCode GetPressure(JacRes *jr);
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user);
PetscErrorCode CheckElasticProperties(JacRes *jr, UserCtx *user);
//PetscErrorCode ShowValues(JacRes *jr, PetscInt n); //, UserCtx *user);
PetscErrorCode ShowValues(JacRes *jr, UserCtx *user, PetscInt n);
PetscErrorCode SaveVelocitiesForSeismicStation(JacRes *jr, UserCtx *user);
//PetscErrorCode PutSeismicSource(JacRes *jr, AdvCtx *actx, UserCtx *user);
PetscErrorCode PrintStress(JacRes *jr);
PetscErrorCode UpdateHistoryFieldsAndGetAxialStressStrain(JacRes *jr, PetscScalar *axial_stress); //, PetscScalar *axial_strain);
//PetscErrorCode ModifyStress(JacRes *jr);
PetscErrorCode GetCellCoordinatesSourceAndSeismicStation(JacRes *jr);
PetscErrorCode GetStressFromSource(JacRes *jr, UserCtx *User, PetscInt i, PetscInt j, PetscInt k, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz);
PetscScalar GetBoundaryDamping(UserCtx *user, PetscInt i, PetscInt j, PetscInt k);
#endif