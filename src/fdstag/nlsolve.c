//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "interface.h"
#include "Assembly_FDSTAG.h"
#include "LaMEMLib_FDSTAG_private.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLCtxCreate"
PetscErrorCode NLCtxCreate(
	NLCtx       *nlctx,
	BlockMat    *bmat,
	FDSTAG      *fs,
	BCCtx       *cbc,
	BCCtx       *sbc,
	JacResCtx   *jrctx)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // initialize application contexts
    nlctx->bmat  = bmat,
    nlctx->fs    = fs;
    nlctx->cbc    = cbc;
    nlctx->sbc    = sbc;
    nlctx->jrctx = jrctx;

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, fs->dofcoupl.numdof, fs->dofcoupl.numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, nlctx, &nlctx->Jac); CHKERRQ(ierr);

	ierr = MatSetUp(nlctx->Jac); CHKERRQ(ierr);

	// postpone Jacobian-vector product definition
	nlctx->jactype = NONE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLCtxDestroy"
PetscErrorCode NLCtxDestroy(NLCtx *nlctx)
{
	PetscErrorCode ierr;
    PetscFunctionBegin;

	ierr = MatDestroy(&nlctx->Jac);  CHKERRQ(ierr);
	ierr = MatDestroy(&nlctx->MFFD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGFormResidual"
PetscErrorCode FDSTAGFormResidual(SNES snes, Vec x, Vec f, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// following parameters are not currently used (stop warning messages):
	if(snes) snes = NULL;

	// access context
	NLCtx       *nlctx = (NLCtx*)ctx;
	FDSTAG      *fs    = nlctx->fs;
    BCCtx       *cbc   = nlctx->cbc;
	JacResCtx   *jrctx = nlctx->jrctx;

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = FDSTAGCopySol(fs, cbc, jrctx, x); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = FDSTAGetEffStrainRate(fs, jrctx); CHKERRQ(ierr);

	// compute residual
	ierr = FDSTAGetResidual(fs, jrctx); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = FDSTAGCopyRes(fs, cbc, jrctx, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);

/*
 	NLCtx       *nlctx = (NLCtx*)ctx;
	BlockMat    *bmat  = nlctx->bmat;
	FDSTAG      *fs    = nlctx->fs;
	JacResCtx   *jrctx = nlctx->jrctx;
	UserContext *user  = nlctx->user;

	// split monolithic solution vector into velocity and pressure blocks
	ierr = BlockMatMonolithicToBlock(bmat, x); CHKERRQ(ierr);
	// copy solution into FDSTAG data structures
	ierr = FDSTAGCopySolution(fs, jrctx, user, bmat->wv, bmat->wp); CHKERRQ(ierr);
	// set boundary ghost points
//	ierr = FDSTAGSetBCGhost(fs, jrctx); CHKERRQ(ierr);
	// compute effective strain rate
	ierr = FDSTAGetEffStrainRate(fs, jrctx); CHKERRQ(ierr);
	// compute residual
	ierr = FDSTAGetResidual(fs, jrctx); CHKERRQ(ierr);
	// copy residual into velocity and pressure blocks
	ierr = FDSTAGCopyResidual(fs, jrctx, user, bmat->wv, bmat->wp);
	// assemble monolithic residual vector from velocity and pressure blocks
	ierr = BlockMatBlockToMonolithic(bmat, f); CHKERRQ(ierr);
*/
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGFormJacobian"
PetscErrorCode FDSTAGFormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	// Compute FDSTAG Jacobian matrix and preconditioner

	PetscInt it, maxit = 5;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear unused parameters
	if(Amat) Amat = NULL;
	if(Pmat) Pmat = NULL;

	// access context
	NLCtx     *nlctx = (NLCtx*)ctx;
	FDSTAG    *fs    = nlctx->fs;
    BCCtx     *cbc    = nlctx->cbc;
	JacResCtx *jrctx = nlctx->jrctx;
	BlockMat  *bmat  = nlctx->bmat;

	// assemble Picard matrix (preconditioner)
	ierr = BlockMatCompute(bmat, fs, cbc, jrctx); CHKERRQ(ierr);

	// in case no Jacobian has been set yet (start with Picard)
	if(nlctx->jactype == NONE)
	{
		// set Picard Jacobian
		ierr = MatShellSetOperation(nlctx->Jac, MATOP_MULT, (void(*)(void))JacApplyPicard); CHKERRQ(ierr);

		ierr = MatAssemblyBegin(nlctx->Jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (nlctx->Jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

		// activate Picard Jacobian type
		nlctx->jactype = PICARD;
	}
	// switch form Picard to MFFD Newton after fixed number of iterations
	else if(nlctx->jactype == PICARD)
	{
		// get iteration number
		ierr = SNESGetIterationNumber(snes, &it); CHKERRQ(ierr);

		if(it == maxit)
		{
			// create MFFD Jacobian
			ierr = JacCreateMFFD(nlctx); CHKERRQ(ierr);

			// set MFFD Jacobian
			ierr = MatShellSetOperation(nlctx->Jac, MATOP_MULT, (void(*)(void))JacApplyMFFD); CHKERRQ(ierr);

			ierr = MatAssemblyBegin(nlctx->Jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
			ierr = MatAssemblyEnd  (nlctx->Jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

			// activate MFFD Jacobian type
			nlctx->jactype = MFFD;
		}
	}

	// setup MFFD Jacobian if necessary
	if(nlctx->jactype == MFFD)
	{
		ierr = JacComputeMFFD(nlctx, snes, x); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateMFFD"
PetscErrorCode JacCreateMFFD(NLCtx *nlctx)
{
	// create matrix-free finite-difference Jacobian approximation

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	FDSTAG *fs = nlctx->fs;

	// create MFFD shell matrix
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, fs->dofcoupl.numdof, fs->dofcoupl.numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, &nlctx->MFFD); CHKERRQ(ierr);

	// Database options
	// -mat_mffd_type - wp or ds (see MATMFFD_WP or MATMFFD_DS)
	// -mat_mffd_err <error_rel> - Sets error_rel
	// -mat_mffd_unim <umin> - Sets umin (for default PETSc routine that computes h only)
	// -mat_mffd_check_positivity	-

	ierr = MatSetFromOptions(nlctx->MFFD); CHKERRQ(ierr);

	ierr = MatSetUp(nlctx->MFFD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacComputeMFFD"
PetscErrorCode JacComputeMFFD(NLCtx *nlctx, SNES snes, Vec x)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	JacResCtx *jrctx = nlctx->jrctx;

	// set residual evaluation function
	ierr = MatMFFDSetFunction(nlctx->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);

	// set base vector (current solution)
	ierr = MatMFFDSetBase(nlctx->MFFD, x, jrctx->gres); CHKERRQ(ierr);

	// assemble matrix
	ierr = MatAssemblyBegin(nlctx->MFFD, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (nlctx->MFFD, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyMFFD"
PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y)
{
	NLCtx *nlctx;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, &nlctx); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult(nlctx->MFFD, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyPicard"
PetscErrorCode JacApplyPicard(Mat A, Vec x, Vec y)
{
	NLCtx    *nlctx;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, &nlctx); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult(nlctx->bmat->A, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateAnalytic"
PetscErrorCode JacCreateAnalytic(NLCtx *nlctx)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateApprox"
PetscErrorCode JacCreateApprox(NLCtx *nlctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateFDColor"
PetscErrorCode JacCreateFDColor(NLCtx *nlctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateFDApprox"
PetscErrorCode JacCreateFDApprox(NLCtx *nlctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacCreateMF"
PetscErrorCode JacCreateMF(NLCtx *nlctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESBlockStopTest"
PetscErrorCode SNESBlockStopTest(SNES snes, PetscInt it, PetscReal xnorm,
	PetscReal gnorm, PetscReal f, SNESConvergedReason *reason, void *cctx)
{
	// monitor residual & perform stop test

	Vec         Bx, Bdx;
	PetscInt    maxit;
	PetscScalar stol;
	PetscScalar nrmSolUpVels, nrmSolVels, resVels, tolVels;
	PetscScalar nrmSolUpPres, nrmSolPres, resPres, tolPres;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// following parameters are not currently used (stop warning messages):
	if(xnorm) xnorm = 0.0;
	if(gnorm) gnorm = 0.0;
	if(f)     f     = 0.0;

	// access context
	NLCtx    *nlctx = (NLCtx*)cctx;
	BlockMat *bmat  = nlctx->bmat;

	// get relative tolerance & maximum number of iterations
	ierr = SNESGetTolerances(snes, NULL, NULL, &stol, &maxit, NULL); CHKERRQ(ierr);

	// later may assign different tolerances for velocity and pressure
	tolVels = stol;
	tolPres = stol;

	// get norms of solution update
	ierr = SNESGetSolutionUpdate(snes, &Bdx);               CHKERRQ(ierr);
	ierr = BlockMatMonolithicToBlock(bmat, Bdx);            CHKERRQ(ierr);
	ierr = VecNorm(bmat->wv, NORM_INFINITY, &nrmSolUpVels); CHKERRQ(ierr);
	ierr = VecNorm(bmat->wp, NORM_INFINITY, &nrmSolUpPres); CHKERRQ(ierr);

	// get norms of current solution
	ierr = SNESGetSolution(snes, &Bx);                      CHKERRQ(ierr);
	ierr = BlockMatMonolithicToBlock(bmat, Bx);             CHKERRQ(ierr);
	ierr = VecNorm(bmat->wv, NORM_INFINITY, &nrmSolVels);   CHKERRQ(ierr);
	ierr = VecNorm(bmat->wp, NORM_INFINITY, &nrmSolPres);   CHKERRQ(ierr);

	// fuses
	if(nrmSolUpVels == 0.0) nrmSolUpVels = 1.0;
	if(nrmSolUpPres == 0.0) nrmSolUpPres = 1.0;
	if(nrmSolVels   == 0.0) nrmSolVels   = nrmSolUpVels;
	if(nrmSolPres   == 0.0) nrmSolPres   = nrmSolUpPres;

	// compute residuals
	resVels = nrmSolUpVels/nrmSolVels;
	resPres = nrmSolUpPres/nrmSolPres;

	// print residuals and tolerances
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#========================================================================================================/\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"# SNES Iteration # /  vel.  res.:         vel.  tol.:         /  pres. res.:         pres. tol.:         /\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#  %3D             /  %14.12e  %14.12e  /  %14.12e  %14.12e  /\n", it, resVels, tolVels, resPres, tolPres);     CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#========================================================================================================/\n"); CHKERRQ(ierr);

	// check convergence
	if(resVels <= tolVels && resPres <= tolPres)
	{
//		ierr = PetscPrintf(PETSC_COMM_WORLD, "Nonlinear solution converged\n"); CHKERRQ(ierr);
		*reason = SNES_CONVERGED_SNORM_RELATIVE;
	}
	else if(it+1 == maxit)
	{
//		ierr = PetscPrintf(PETSC_COMM_WORLD, "Nonlinear solution reached maximum number iterations\n"); CHKERRQ(ierr);
		*reason = SNES_DIVERGED_MAX_IT;
	}
	else
	{
		// continue iterations
		*reason = SNES_CONVERGED_ITERATING;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESPrintConvergedReason"
PetscErrorCode SNESPrintConvergedReason(SNES snes)
{

	SNESConvergedReason reason;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = SNESGetConvergedReason(snes, &reason);

    // CONVERGENCE

    if(reason == SNES_CONVERGED_FNORM_ABS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason: ||F|| < atol \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_FNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason: ||F|| < rtol*||F_initial|| \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_SNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason: Newton computed step size small; || delta x || < stol || x ||\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_ITS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_TR_DELTA)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason: SNES_CONVERGED_TR_DELTA\n"); CHKERRQ(ierr);
	}

    // DIVERGENCE

	else if(reason == SNES_DIVERGED_FUNCTION_DOMAIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: the new x location passed the function is not in the domain of F\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FUNCTION_COUNT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: too many function evaluations\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINEAR_SOLVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: the linear solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FNORM_NAN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: residual norm is NAN\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_MAX_IT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINE_SEARCH)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: the line search failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_INNER)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: the inner solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LOCAL_MIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason: || J^T b || is small, implies converged to local minimum of F\n"); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------


/*


	//====================================================================

	// Important SNES functions:
 		SNESCreate
		SNESSetFunction
		SNESSetJacobian
		SNESGetKSP
		SNESSetFromOptions
		SNESSolve
		SNESDestroy
		SNESSetConvergenceTest
		SNESMonitorSet
		SNESSetLagJacobian
		SNESSetLagPreconditioner
		SNESVISetVariableBounds
		SNESGetSolution
		SNESGetConvergedReason
		SNESVISetVariableBounds

PetscErrorCode PreCheck(SNESLineSearch,Vec,Vec,PetscBool*,void*)
PetscErrorCode PostCheck(SNESLineSearch,Vec,Vec,Vec,PetscBool*,PetscBool*,void*)

		PetscErrorCode  SNESGetConvergedReason(SNES snes, SNESConvergedReason *reason)
	//====================================================================

		Eisenstat-Walker method to reset relative tolerances for linear solve
		according to the progress of the nonlinear solve

	//====================================================================


	//====================================================================

	// JACOBIAN FINITE DIFFERENCE COLORING APPROXIMATION (ASSEMBLED FORM FOR TESTING)

	// Sparsity structure should be assembled in the mat

	PetscErrorCode  MatGetColoring(Mat mat, const MatColoringType type, ISColoring *iscoloring)

	PetscErrorCode  MatFDColoringCreate(Mat mat, ISColoring iscoloring, MatFDColoring *color)

	PetscErrorCode  MatFDColoringSetFromOptions(MatFDColoring matfd)

	PetscErrorCode ISColoringDestroy(ISColoring *iscoloring)

	PetscErrorCode  MatFDColoringSetFunction(MatFDColoring matfd,PetscErrorCode (*f)(void),void *fctx)

	Notes: This function is usually used automatically by SNES (when one uses SNESSetJacobian() with the argument SNESDefaultComputeJacobianColor())
	and only needs to be used by someone computing a matrix via coloring directly by calling MatFDColoringApply()

	PetscErrorCode  MatFDColoringApply(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)

	//====================================================================


//---------------------------------------------------------------------------


*/
