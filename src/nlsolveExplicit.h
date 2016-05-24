/*
 * nlsolveExplicit.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Beatriz
 */

#ifndef __nlsolveExplicit_h__
#define __nlsolveExplicit_h__

PetscErrorCode GetVelocities(JacRes *jr);
PetscErrorCode FormMomentumResidualAndTheta(SNES snes, Vec x, Vec gK, void *ctx);
PetscErrorCode GetPressure(JacRes *jr);

#endif
