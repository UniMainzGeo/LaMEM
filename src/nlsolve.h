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
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __nlsolve_h__
#define __nlsolve_h__
//---------------------------------------------------------------------------
// Jacobian type
enum JacType
{
	_PICARD_, // constant effective coefficients approximation (viscosity, conductivity, stress)
	_MFFD_    // built-in finite difference approximation

};
//---------------------------------------------------------------------------

struct NLSol
{
	JacRes     *jr;      // Jacobian-residual context
	PCData      pc;      // preconditioner context
	Mat         MFFD;    // matrix-free finite difference Jacobian
	Mat         PICARD;  // Picard Jacobian
	JacType     jtype;   // actual type of Jacobian operator
	PetscInt    it;      // iteration counter
	PetscInt    it_Nwt;  // newton iteration counter
	PetscScalar refRes;  // reference residual norm
	PetscScalar rtolPic; // relative Picard residual reduction tolerance
	PetscInt    nNwtIt;  // number of Newton iterations before switch to Picard
	PetscScalar rtolNwt; // Newton divergence tolerance
};
//---------------------------------------------------------------------------

PetscErrorCode NLSolCreate(SNES *p_snes, JacRes *jr);

PetscErrorCode NLSolDestroy(SNES *p_snes);

// compute residual vector
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx);

// compute Jacobian matrix and preconditioner
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx);

PetscErrorCode JacApply(Mat A, Vec x, Vec y);

PetscErrorCode SNESCoupledTest(
	SNES                snes,
	PetscInt            it,
	PetscReal           xnorm,
	PetscReal           gnorm,
	PetscReal           f,
	SNESConvergedReason *reason,
	void                *cctx);

PetscErrorCode SNESPrintConvergedReason(SNES snes, PetscLogDouble t_beg);

//---------------------------------------------------------------------------
#endif
