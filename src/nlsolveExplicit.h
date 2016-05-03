/*
 * nlsolveExplicit.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Beatriz
 */

#ifndef __nlsolveExplicit_h__
#define __nlsolveExplicit_h__

PetscErrorCode NLSolverExp(JacRes *jr);
PetscErrorCode FormMomentumResidual(SNES snes, Vec x, Vec f, void *ctx);
PetscErrorCode GetVelocities(JacRes *jr, Vec x, Vec f);

#endif
