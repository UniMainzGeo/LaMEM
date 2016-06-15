/*
 * nlsolveExplicit.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Beatriz
 */

#ifndef __nlsolveExplicit_h__
#define __nlsolveExplicit_h__

PetscErrorCode GetVelocities(JacRes *jr);
PetscErrorCode FormMomentumResidualAndTheta(Vec x, void *ctx);
PetscErrorCode GetPressure(JacRes *jr);
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user);

#endif
