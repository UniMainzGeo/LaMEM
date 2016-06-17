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
	PetscScalar      Perturb;                 // Perturbation parameter for the finite differences
	PetscScalar      Fini;
	Vec              dF;
	Vec              xini;
	Vec 			 pro;
	Vec              vx, vy, vz;
	PetscInt         count;

	UserCtx          *user;
	JacRes           *jr;
	NLSol            *nl;
	SNES             snes;
} AdjGrad;

// Perform the adjoint inversion
PetscErrorCode AdjointOptimization(JacRes *jr, AdjGrad *aop, UserCtx *user, NLSol *nl,SNES snes);

// Compute the gradients for the adjoint inversion
PetscErrorCode AdjointComputeGradients(Tao tao, Vec P, Vec grad, void *ctx);

// Compute the objective function that will be minimized
PetscErrorCode AdjointObjectiveFunction(Tao tao, Vec x, PetscReal *F, void *ctx);

// Interpolate the adjoint points and include them into the projection vector
PetscErrorCode AdjointPointInPro(JacRes *jr, UserCtx *user, AdjGrad *aop);

// Perturb the input parameters within the gradient computation
PetscErrorCode AdjointGradientPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop);

// reset the perturbed input parameter within the gradient computation
PetscErrorCode AdjointGradientResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop);

// Perturb the parameters after every inversion iteration
PetscErrorCode AdjointPerturbParameterVec(NLSol *nl,Vec P, UserCtx *user);

// Extract the intital parameter Vector for the inversion
PetscErrorCode AdjointExtractParameterVec(NLSol *nl,Vec P, UserCtx *user);

// Read the comparison solution vector
PetscErrorCode AdjointLoadCompareData(void *ctx);

// To clear the memory
PetscErrorCode AdjointDestroy(AdjGrad *aop);

#endif
