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

PetscErrorCode MatFreeApplyRestrict(Mat R, Vec vf, Vec vcb, Vec vc);

PetscErrorCode MatFreeApplyProlong(Mat P, Vec vc, Vec vfb, Vec vf);

// split vector into blocks, assign ghost points
PetscErrorCode MatFreeSplitVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp);

// assemble ghost point contributions, combine blocks into a vector, enforce boundary conditions
PetscErrorCode MatFreeAssembleVec(MatData *md, Vec v, Vec lvx, Vec lvy, Vec lvz, Vec gvp);

// combine blocks into a vector, enforce boundary conditions
PetscErrorCode MatFreeCombineVec(MatData *md, Vec v, Vec gvx, Vec gvy, Vec gvz, Vec gvp);

PetscErrorCode MatFreeGetPicard(MatData *md,
		Vec lvx, Vec lvy, Vec lvz, Vec gp,
		Vec lfx, Vec lfy, Vec lfz, Vec gc);

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

/*
1. Set PC_MG_GALERKIN_NONE

2. Call PCMGSetRestriction and PCMGSetInterpolation to set restriction/interpolation matrices

    a. On the top levels set the shell matrices equipped with MATOP_MULT_ADD

    b. On the bottom levels set assembled matrices

3. Call PCMGGetSmoother and KSPSetOperators to set linear operators

    a. On the top levels set the shell matrices equipped with MATOP_GET_DIAGONAL and MATOP_MULT

    b. On the bottom levels generate the operators by explicitly calling MatMatMatMult (R is not the same as P)

        Use MAT_INITIAL_MATRIX for the first time

        Use MAT_REUSE_MATRIX for the subsequent calls

PetscErrorCode MatMult(Mat mat, Vec x, Vec y); // y = A*x

PetscErrorCode MatMultAdd(Mat mat, Vec v1, Vec v2, Vec v3); // v3 = v2 + A∗v1

PetscErrorCode MatGetDiagonal(Mat mat,  Vec v); // v = diag(A)

PetscCall(MatShellSetOperation(A, MATOP_GET_DIAGONAL ,(void(*)(void))MatGetDiagonal_Laplacian2D))

PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A, Vec diag)

*/
