//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __nlsolve_h__
#define __nlsolve_h__
//---------------------------------------------------------------------------

// Enumeration of Jacobian types
typedef enum
{
	NONE,
	//===================
	// assembled matrices
	//===================
	PICARD,   // constant effective coefficients approximation (viscosity, conductivity, stress)
//	ANALYTIC, // analytic Jacobian with full sparsity pattern
//	APPROX,   // analytic Jacobian truncated to Picard sparsity pattern (possibly with diagonal compensation)
//	FDCOLOR,  // finite difference coloring approximation with full sparsity pattern
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
	// application contexts
	Mat            Jac;        // Jacobian matrix (SHELL)
	Mat            MFFD;       // matrix-free finite difference Jacobian

	JacRes *jr;      // Jacobian & residual evaluation context

	PCStokes  pc;
	SNES           snes;   // nonlinear solver
	SNESLineSearch snesls; // line search context
	KSP            ksp;
//	PC             pc;

	JacType   jactype;    // actual type of Jacobian operator

} NLSolver;
//---------------------------------------------------------------------------
/*
PetscErrorCode NLCtxCreate(
	NLCtx       *nlctx,
	FDSTAG      *fs,
	BCCtx       *cbc,
	BCCtx       *sbc,
	JacResCtx   *jrctx);

PetscErrorCode NLCtxDestroy(NLCtx *nlctx);

//---------------------------------------------------------------------------

// compute FDSTAG block residual vector
PetscErrorCode FDSTAGFormResidual(SNES snes, Vec x, Vec f, void *ctx);

// compute FDSTAG block Jacobian matrix and preconditioner
PetscErrorCode FDSTAGFormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx);

//---------------------------------------------------------------------------

PetscErrorCode JacCreateMFFD(NLCtx *nlctx);

PetscErrorCode JacComputeMFFD(NLCtx *nlctx, SNES snes, Vec x);

PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y);

//---------------------------------------------------------------------------

PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec y);

//---------------------------------------------------------------------------

//PetscErrorCode JacCreateAnalytic(NLCtx *nlctx);

//PetscErrorCode JacCreateApprox  (NLCtx *nlctx);

//PetscErrorCode JacCreateFDColor (NLCtx *nlctx);

//PetscErrorCode JacCreateFDApprox(NLCtx *nlctx);

//PetscErrorCode JacCreateMF      (NLCtx *nlctx);

//---------------------------------------------------------------------------

//PetscErrorCode SNESBlockStopTest(SNES snes, PetscInt it, PetscReal xnorm,
//	PetscReal gnorm, PetscReal f, SNESConvergedReason *reason, void *cctx);

//PetscErrorCode SNESPrintConvergedReason(SNES snes);
*/
//---------------------------------------------------------------------------
#endif
