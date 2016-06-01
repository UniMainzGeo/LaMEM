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
	PetscScalar      grad[_MAX_AdjointPars_]; // Vector containing the gradients (dF/dp = -psi * dr/dp)   // _MAX_AdjointIndices_ defined in fdstagTypes.h (needed for user input as well)
	PetscScalar      Ini;                     // Initial value of perturbed parameter
	PetscScalar      Perturb;                 // Perturbation parameter for the finite differences
} AdjGrad;

// make context for adjoint
PetscErrorCode CreateAdjoint(JacRes *jr, UserCtx *user, AdjGrad *aop, NLSol *nl,SNES snes);

// Interpolate the adjoint points and include them into the projection vector
PetscErrorCode AdjointPointInPro(JacRes *jr, UserCtx *user, PetscScalar *vx, PetscScalar *vy, PetscScalar *vz, Vec pro );

// Perturb the input parameters
PetscErrorCode AdjointPerturbParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop);

// reset the perturbed input parameter
PetscErrorCode AdjointResetParameter(NLSol *nl, PetscInt CurPar, PetscInt CurPhase, AdjGrad *aop);

// To clear the memory
PetscErrorCode AdjointDestroy(AdjGrad *aop);

#endif
