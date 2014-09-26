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
	PICARD,   // constant effective coefficients approximation (viscosity, conductivity, stress)
//	FDCOLOR,  // finite difference coloring approximation with full sparsity pattern
//	ANALYTIC, // analytic Jacobian with full sparsity pattern
//	APPROX,   // analytic Jacobian truncated to Picard sparsity pattern (possibly with diagonal compensation)
//	FDAPPROX, // finite difference coloring approximation truncated to Picard sparsity pattern
	//============
	// matrix-free
	//============
	MF,  // analytic
	MFFD // built-in finite difference approximation

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

//PetscErrorCode SNESBlockStopTest(SNES snes, PetscInt it, PetscReal xnorm,
//	PetscReal gnorm, PetscReal f, SNESConvergedReason *reason, void *cctx);

//PetscErrorCode SNESPrintConvergedReason(SNES snes);

//---------------------------------------------------------------------------
#endif
