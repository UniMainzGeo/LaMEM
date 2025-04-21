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

struct MatData;

//---------------------------------------------------------------------------

PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec f);

PetscErrorCode MatFreeGetPicard(MatData *md,
		Vec lvx, Vec lvy, Vec lvz, Vec gp,
		Vec lfx, Vec lfy, Vec lfz, Vec gc);

// access solution vector
PetscErrorCode MatFreeGetSol(MatData *md, Vec x, Vec lvx, Vec lvy, Vec lvz, Vec gp);

// assemble residual
PetscErrorCode MatFreeGetRes(MatData *md, Vec f, Vec lfx, Vec lfy, Vec lfz, Vec gc);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// set velocity two-point constraint
#define SET_VEL_TPC(bc, i, j, k, ind, lim, icf, gcf) { if(ind == lim) { gcf = 0.0; if(bc[k][j][i] != DBL_MAX) { icf = 2.0; } else { icf = 0.0; } } }

//---------------------------------------------------------------------------
#endif
