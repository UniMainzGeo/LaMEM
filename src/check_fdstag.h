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
 **    filename:   check_fdstag.h
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
//.........   LaMEM - FDSTAG CANONICAL INTERFACE CHECKING ROUTINES   ........
//---------------------------------------------------------------------------
#ifndef __check_fdstag_h__
#define __check_fdstag_h__

//---------------------------------------------------------------------------

PetscErrorCode DoJacTests(NLSol *nl);

PetscErrorCode DoPicardTests(NLSol *nl);

//---------------------------------------------------------------------------
/*
PetscErrorCode DoDarcyTests(NLCtx *nlctx, UserCtx *user);

PetscErrorCode DarcyPostProcess(NLCtx *nlctx, UserCtx *user);

PetscErrorCode DoMGTests(NLCtx *nlctx, PVOut *pvout);

PetscErrorCode BlockSetUpMGVelBlock(NLCtx *nlctx, MGCtx *mg);

PetscErrorCode FDSTAGOutputResTest(FDSTAG *fs, JacResCtx *jrctx, Vec f);

PetscErrorCode GetLinRes(Mat A, Vec x, Vec rhs, Vec res);

PetscErrorCode FieldSplitTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime);

PetscErrorCode MGTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime);

//---------------------------------------------------------------------------
typedef struct
{
	PetscScalar a;
	PetscScalar b;

}MaxDiff;
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar A[8];
	PetscScalar cxb, cxe;
	PetscScalar cyb, cye;
	PetscScalar czb, cze;

} BCValues;

PetscScalar InterpolateLinear3D(PetscScalar cx, PetscScalar cy, PetscScalar cz,  BCValues BC);

// Initialize velocity vectors to test strain rates.
// Selects the velocity direction x, y, z (0, 1, 2) with vectDir,
// and applies constant gradient in the gradDir direction
// between the values begVal & endVal (in the positive direction of gradDir).
// Initializes boundary ghost points in the tangential directions accordingly.

PetscErrorCode InitVelocityTest(
	JacRes      *jr,
	UserCtx     *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal);

PetscErrorCode JacResCtxClearVelocity(JacRes *jr, PetscInt vectDir);

PetscErrorCode StrainRateSingleComp(
	JacRes      *jr,
	UserCtx     *usr,
	PVOut       *pvout,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  ttime,
	PetscInt     itime);

PetscErrorCode StrainRateInterpTest(
	JacRes      *jr,
	UserCtx     *usr,
	PVOut       *pvout);

//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar A[8];
	PetscScalar cxb, cxe;
	PetscScalar cyb, cye;
	PetscScalar czb, cze;

} BCValues;

PetscScalar InterpolateLinear3D(PetscScalar cx, PetscScalar cy, PetscScalar cz,  BCValues BC);

// Initialize velocity vectors to test velocity gradients
// Selects the velocity direction x, y, z (0, 1, 2) with vectDir,
// and applies constant gradient in the gradDir direction
// between the values begVal & endVal (in the positive direction of gradDir).
// Initializes boundary ghost points in the tangential directions accordingly.

PetscErrorCode InitVelocityTest(
	JacRes      *jr,
	UserCtx     *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal,
	Tensor2RN    *L);

PetscErrorCode ClearVelocity(JacRes *jr);

PetscErrorCode VelGradSingleComp(
		JacRes      *jr,
		UserCtx     *usr,
		PetscInt     vectDir,
		PetscInt     gradDir,
		PetscScalar  begVal,
		PetscScalar  endVal);

PetscErrorCode VelGradTest(
	JacRes      *jr,
	UserCtx     *usr);
*/

//---------------------------------------------------------------------------
#endif
