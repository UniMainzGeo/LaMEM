//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __nlsolve_h__
#define __nlsolve_h__
//---------------------------------------------------------------------------

// Jacobian type
typedef enum
{
	//===================
	// assembled matrices
	//===================
	_PICARD_,   // constant effective coefficients approximation (viscosity, conductivity, stress)
//	_FDCOLOR_,  // finite difference coloring approximation with full sparsity pattern
//	_ANALYTIC_, // analytic Jacobian with full sparsity pattern
//	_APPROX_,   // analytic Jacobian truncated to Picard sparsity pattern (possibly with diagonal compensation)
//	_FDAPPROX_, // finite difference coloring approximation truncated to Picard sparsity pattern
	//============
	// matrix-free
	//============
	_MF_,  // analytic
	_MFFD_ // built-in finite difference approximation

} JacType;

//---------------------------------------------------------------------------
typedef struct
{
	Mat       J;      // Jacobian matrix
	Mat       P;      // preconditioner
	Mat       MFFD;   // matrix-free finite difference Jacobian
	PCStokes  pc;     // Stokes preconditioner
	JacType   jtype;  // actual type of Jacobian operator
	PetscInt  nPicIt; // number of Picard iteraions before switch to Newton

} NLSol;
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

PetscErrorCode SNESPrintConvergedReason(SNES snes);

//PetscErrorCode SNESBlockStopTest(SNES snes, PetscInt it, PetscReal xnorm,
//	PetscReal gnorm, PetscReal f, SNESConvergedReason *reason, void *cctx);

//---------------------------------------------------------------------------
#endif
