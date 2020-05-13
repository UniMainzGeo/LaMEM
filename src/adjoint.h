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
 **    filename:   adjoint.h
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

struct Scaling;
struct FDSTAG;
struct FreeSurf;
struct DBMat;
struct Tensor2RN;
struct PData;
struct Material_t;
struct JacRes;
struct Controls;
struct NLSol;
struct ModParam;
struct FB;

// Some global maxes on parameter and index numbers
#define _MAX_PAR_ 50
#define _MAX_IND_ 100

// Structure that holds paramters for the adjoint gradient computation
struct AdjGrad
{
	PetscScalar      Ini;                     // Initial value of perturbed parameter
	PetscScalar      Ini2;                    // If n is the parameter we need two initials
	PetscScalar      Perturb;                 // Perturbation parameter for the finite differences
	PetscScalar      CurScal, CurScalst;
	PetscScalar      DII_ref;                 // SUPER UNNECESSARY but DII is otherwise not accesible
	Vec              dF, dFg, dFst;
	Vec 			 pro, stpro;
	Vec              vx, vy, vz, stx, sty, stz;
	Vec              dphidu;
	Vec              gradfield;                // Used if gradient at every point is computed (same size as jr->p)
};

// Structure that holds vectors required by TAO
struct Adjoint_Vecs
{
	Vec             val, Ub, Lb, grad, P;
};

// Adjoint optimization driving routines
PetscErrorCode AdjointOptimisation(Vec P, PetscScalar F, Vec grad, void *ctx);
PetscErrorCode AdjointOptimisationTAO(Tao tao, Vec P, PetscReal *F, Vec grad, void *ctx);
PetscErrorCode LaMEMAdjointReadInputSetDefaults(ModParam **p_IOparam, FB **p_fb, Adjoint_Vecs *Adjoint_Vectors);
PetscErrorCode LaMEMAdjointMain(ModParam *IOparam, FB *fb);

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

// Gradient function for field sensitivity for rho (FD approximation)
PetscErrorCode AdjointFormResidualFieldFDRho(SNES snes, Vec x, Vec psi, NLSol *nl, AdjGrad *aop );

// To clear the memory
PetscErrorCode AdjointDestroy(AdjGrad *aop);

#endif
