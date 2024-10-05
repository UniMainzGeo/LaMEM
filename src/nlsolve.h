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
	//===================
	// assembled matrices
	//===================
	_PICARD_,   // constant effective coefficients approximation (viscosity, conductivity, stress)

	//============
	// matrix-free
	//============
	_MFFD_ // built-in finite difference approximation

};
//---------------------------------------------------------------------------
struct NLSol
{
	Mat       J;      // Jacobian matrix
	Mat       P;      // preconditioner
	Mat       MFFD;   // matrix-free finite difference Jacobian
	PCStokes  pc;     // Stokes preconditioner

	JacType     jtype;    // actual type of Jacobian operator
	PetscInt    it;       // iteration counter
	PetscInt    it_Nwt;   // newton iteration counter
	PetscScalar refRes;   // reference residual norm
	PetscInt    nPicIt;   // number of Picard iteraions before switch to Newton
	PetscScalar rtolPic;  // relative Picard residual reduction tolerance
	PetscInt    nNwtIt;   // number of Newton iterations before switch to Picard
	PetscScalar rtolNwt;  // Newton divergence tolerance
} ;

//---------------------------------------------------------------------------

PetscErrorCode NLSolClear(NLSol *nl);

PetscErrorCode NLSolCreate(NLSol *nl, PCStokes pc, SNES *p_snes);

PetscErrorCode NLSolDestroy(NLSol *nl);

//---------------------------------------------------------------------------

// compute residual vector
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx);

// compute Jacobian matrix and preconditioner
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx);

//---------------------------------------------------------------------------

PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y);

//---------------------------------------------------------------------------

PetscErrorCode SNESPrintConvergedReason(SNES snes, 	PetscLogDouble t_beg);

PetscErrorCode SNESCoupledTest(
	SNES                snes,
	PetscInt            it,
	PetscReal           xnorm,
	PetscReal           gnorm,
	PetscReal           f,
	SNESConvergedReason *reason,
	void                *cctx);

//---------------------------------------------------------------------------

PetscErrorCode DisplaySpecifiedSolverOptions(PCStokes pc, SNES snes);

//---------------------------------------------------------------------------

#endif
