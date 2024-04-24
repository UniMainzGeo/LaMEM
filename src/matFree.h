/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...............   MATRIX-FREE JACOBIAN AND PRECONDITIONER  ................
//---------------------------------------------------------------------------
#ifndef __matFree_h__
#define __matFree_h__

//---------------------------------------------------------------------------

struct JacRes;

//---------------------------------------------------------------------------

PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec y);

PetscErrorCode JacResPicardMatFree(JacRes *jr);

//---------------------------------------------------------------------------
/*
PetscErrorCode JacApplyJacobian(Mat A, Vec x, Vec y);

PetscErrorCode JacResGetJ2Derivatives(JacRes *jr);

PetscErrorCode JacResJacobianMatFree(JacRes *jr);
*/
//---------------------------------------------------------------------------
#endif
