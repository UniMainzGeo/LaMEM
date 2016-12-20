/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   advect.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
// COMPUTATION OF ADJOINT GRADIENTS
//---------------------------------------------------------------------------
#ifndef __adjoint_h__
#define __adjoint_h__
//---------------------------------------------------------------------------
// Structure that holds paramters for the adjoint gradient computation
typedef struct
{
	PetscScalar      Ini;                     // Initial value of perturbed parameter
	PetscScalar      Ini2;                     // If n is the parameter we need two initials
	PetscScalar      Perturb;                 // Perturbation parameter for the finite differences
	PetscScalar      CurScal;
	Vec              dF;
	Vec 			 pro;
	Vec              vx, vy, vz;
} AdjGrad;

// Compute the gradients for the adjoint inversion
PetscErrorCode AdjointObjectiveAndGradientFunction(AdjGrad *aop, JacRes *jr, NLSol *nl, ModParam *IOparam, SNES snes, FreeSurf *surf);

// Compute the gradients for the adjoint inversion
PetscErrorCode AdjointComputeGradients(JacRes *jr, AdjGrad *aop, NLSol *nl, SNES snes, ModParam *IOparam, FreeSurf *surf);

// Interpolate the adjoint points and include them into the projection vector
PetscErrorCode AdjointPointInPro(JacRes *jr, AdjGrad *aop, ModParam *IOparam, FreeSurf *surf);

// Perturb the input parameters within the gradient computation
PetscErrorCode AdjointGradientPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop, Scaling *scal);

// reset the perturbed input parameter within the gradient computation
PetscErrorCode AdjointGradientResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop);

// Manage the analytical residual functions
PetscErrorCode AdjointFormResidual(SNES snes, Vec x, Vec f, void *ctx, PetscInt CurPar, PetscInt CurPhase );

// get analytical residual for the density
PetscErrorCode AdjointJacResGetResidual_Density(JacRes *jr, PetscInt CurPar, PetscInt CurPhase);

// get analytical residual for the reference viscosity for powerlaw rheology
PetscErrorCode AdjointJacResGetResidual_ViscPowerlaw(JacRes *jr, PetscInt CurPar, PetscInt CurPhase);

// To clear the memory
PetscErrorCode AdjointDestroy(AdjGrad *aop);

#endif
