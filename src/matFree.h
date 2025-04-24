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

PetscErrorCode MatFreeApplyPicard(Mat A, Vec x, Vec f);

PetscErrorCode MatFreeApplyPreconditioner(Mat A, Vec x, Vec f);

PetscErrorCode MatFreeApplyLinearOperator(MatData *md, Vec x, Vec f, PetscScalar cfInvEta);

PetscErrorCode MatFreeApplyRestrict(Mat R, Vec vf, Vec vcb, Vec vc);

PetscErrorCode MatFreeApplyProlong(Mat P, Vec vc, Vec vfb, Vec vf);

// split vector into blocks, assign ghost points
PetscErrorCode MatFreeSplitVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp);

// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions
PetscErrorCode MatFreeAssembleVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp);

// combine blocks into a vector, enforce boundary conditions
PetscErrorCode MatFreeCombineVec(MatData *md, Vec v, Vec gvx, Vec gvy, Vec gvz, Vec gvp);

PetscErrorCode MatFreeGetLinearOperator(MatData *md,
		Vec lvx, Vec lvy, Vec lvz, Vec gp,
		Vec lfx, Vec lfy, Vec lfz, Vec gc,
		PetscScalar cfInvEta);

// cfInvEta - inverse viscosity term prefactor
// 0.0      - Picard operator
// 1.0      - preconditioner operator

PetscErrorCode MatFreeGetRestrict(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp);

PetscErrorCode MatFreeGetProlong(
		MatData *coarse, MatData *fine,
		Vec fx, Vec fy, Vec fz, Vec fp,
		Vec cx, Vec cy, Vec cz, Vec cp);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// set velocity two-point constraint
#define SET_VEL_TPC(bc, i, j, k, ind, lim, icf, gcf) { if(ind == lim) { gcf = 0.0; if(bc[k][j][i] != DBL_MAX) { icf = 2.0; } else { icf = 0.0; } } }

//---------------------------------------------------------------------------
#endif

