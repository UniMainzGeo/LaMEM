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
 **    filename:   constEq.h
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
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#ifndef __Tensor_h__
#define __Tensor_h__

//---------------------------------------------------------------------------
//....    Non-Symmetric second rank tensor (gradient & rotation tensors) ....
//---------------------------------------------------------------------------

struct Tensor2RN
{
	PetscScalar xx, xy, xz;
	PetscScalar yx, yy, yz;
	PetscScalar zx, zy, zz;

};

//---------------------------------------------------------------------------
//.......   Symmetric second rank tensor (stress & strain tensors)   ........
//---------------------------------------------------------------------------

struct Tensor2RS
{
	PetscScalar xx, xy, xz;
	PetscScalar     yy, yz;
	PetscScalar         zz;

};

//---------------------------------------------------------------------------
// Elastic stress rotation functions
//---------------------------------------------------------------------------

// compute rotation matrix from axis & angle (Euler-Rodrigues formula)
void GetRotationMatrix(
	Tensor2RN   *R,   // rotation matrix
	PetscScalar  dt,  // time step
	PetscScalar  wx,  // vorticity vector components
	PetscScalar  wy,  // ...
	PetscScalar  wz); // ...

// rotate stress tensor
void RotateStress(Tensor2RN *R, Tensor2RS *S, Tensor2RS *SR);

// copy symmetric second order tensor B = A
void Tensor2RSCopy(Tensor2RS *A, Tensor2RS *B);

//---------------------------------------------------------------------------
// Infinite Strain Axis (ISA) calculation functions
//---------------------------------------------------------------------------
void Tensor2RNClear(Tensor2RN *A);

PetscInt Tensor2RNCheckEq(Tensor2RN *A, Tensor2RN *B, PetscScalar tol);

void Tensor2RNNorm(Tensor2RN *A, PetscScalar *pk);

void Tensor2RSNorm(Tensor2RS *A, PetscScalar *pk);

void Tensor2RNDivide(Tensor2RN *A, PetscScalar k);

void Tensor2RNTrace(Tensor2RN *A);

void Tensor2RNSym(Tensor2RN *A, Tensor2RN *B);

void Tensor2RNProduct(Tensor2RN *A, Tensor2RN *B, Tensor2RN *C);

void Tensor2RNTranspose(Tensor2RN *A, Tensor2RN *B);

void Tensor2RNCopy(Tensor2RN *A, Tensor2RN *B);

void Tensor2RNCopySym(Tensor2RN *A, Tensor2RS *B);

void Tensor2RNUnit(Tensor2RN *A);

void Tensor2RNDivide(Tensor2RN *A, PetscScalar k);

void Tensor2RNSum3(
	Tensor2RN *A, PetscScalar ka,
	Tensor2RN *B, PetscScalar kb,
	Tensor2RN *C, PetscScalar kc,
	Tensor2RN *R);

void Tensor2RNView(Tensor2RN *A, const char *msg);

void Tensor2RSView(Tensor2RS *A, const char *msg);

PetscInt Tensor2RNEigen(Tensor2RN *L, PetscScalar tol, PetscScalar eval[]);

PetscInt Tensor2RSSpectral(
	Tensor2RS   *A,      // symmetric tensor
	PetscScalar eval[],  // eigenvalues (sorted)
	PetscScalar evect[], // eigenvectors (corresponding)
	PetscScalar ttol,    // tight tolerance (convergence condition)
	PetscScalar ltol,    // loose tolerance (divergence condition)
	PetscInt    itmax);  // maximum number rotations

PetscInt getISA(Tensor2RN *pL, PetscScalar ISA[], PetscScalar *plnrm);

PetscErrorCode Tensor2RS2DSpectral(
	PetscScalar  axx,
	PetscScalar  ayy,
	PetscScalar  axy,
	PetscScalar *pa1,
	PetscScalar *pa2,
	PetscScalar  v1[],
	PetscScalar  v2[],
	PetscScalar  tol);

//---------------------------------------------------------------------------
#endif
