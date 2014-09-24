//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "interface.h"
#include "Assembly_FDSTAG.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolCreate"
PetscErrorCode NLSolCreate(NLSol *nl, PCStokes pc, SNES snes)
{
	KSP             ksp;
	PC              ipc;
	SNESLineSearch  ls;
	JacRes         *jr;
	DOFIndex       *dof;
	PetscBool       flg;

    PetscErrorCode ierr;
    PetscFunctionBegin;

    // access context
	nl->pc = pc;
	jr     = pc->jr;
	dof    = &(jr->fs->cdof);

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->J); CHKERRQ(ierr);
	ierr = MatSetUp(nl->J);                              CHKERRQ(ierr);

	// create matrix-free Preconditioner operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->P); CHKERRQ(ierr);
	ierr = MatSetUp(nl->P);                              CHKERRQ(ierr);

	// create finite-difference Jacobian
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, &nl->MFFD); CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(nl->MFFD,"fd_");       CHKERRQ(ierr);
	ierr = MatSetFromOptions(nl->MFFD);               CHKERRQ(ierr);
	ierr = MatSetUp(nl->MFFD);                        CHKERRQ(ierr);

	// setup nonlinear solver
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                     CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESNEWTONLS);                         CHKERRQ(ierr);
	ierr = SNESGetLineSearch(snes, &ls);                            CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);          CHKERRQ(ierr);
	ierr = SNESSetFunction(snes, jr->gres, &FormResidual, &nl);     CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes, nl->J, nl->P, &FormJacobian, &nl); CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);                                CHKERRQ(ierr);

	// setup linear solver & preconditioner
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"js_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          CHKERRQ(ierr);

//	ierr = KSPSetConvergenceTest(ksp, &KSPBlockStopTest, &bmat, NULL);CHKERRQ(ierr);
//	ierr = SNESSetConvergenceTest(snes, SNESBlockStopTest, &nlctx, NULL); CHKERRQ(ierr);

	// set Jacobian type & initial guess
	nl->jtype = JPICARD;
	ierr = VecSet(jr->gsol, 0.0); CHKERRQ(ierr);

	// read number of Picard iterations
	ierr = PetscOptionsGetInt(PETSC_NULL, "-npicard", &nl->nPicIt, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) nl->nPicIt = 5;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolDestroy"
PetscErrorCode NLSolDestroy(NLSol *nl)
{
	PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = MatDestroy(&nl->J);    CHKERRQ(ierr);
	ierr = MatDestroy(&nl->P);    CHKERRQ(ierr);
	ierr = MatDestroy(&nl->MFFD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormResidual"
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx)
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear unused parameters
	if(snes) snes = NULL;

	// access context
	nl = (NLSol*)ctx;
	jr = nl->pc->jr;

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute residual
	ierr = JacResGetResidual(jr); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = JacResCopyRes(jr, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	// Compute FDSTAG Jacobian matrix and preconditioner

	PetscInt it;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	NLSol    *nl = (NLSol*)ctx;
	PCStokes  pc = nl->pc;
	JacRes   *jr = pc->jr;

	// setup preconditioner
	ierr = PCStokesSetup(pc);                                                CHKERRQ(ierr);
	ierr = MatShellSetOperation(Pmat, MATOP_MULT, (void(*)(void))pc->Apply); CHKERRQ(ierr);
	ierr = MatShellSetContext(Pmat, pc->data);                               CHKERRQ(ierr);

	// switch Jacobian after fixed number of iterations
	ierr = SNESGetIterationNumber(snes, &it); CHKERRQ(ierr);
	if(it == nl->nPicIt) nl->jtype = JMFFD;

	// setup Jacobian ...
	if(nl->jtype == JPICARD)
	{
		// ... Picard
		ierr = MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))pc->Picard); CHKERRQ(ierr);
		ierr = MatShellSetContext(Amat, pc->data);                                CHKERRQ(ierr);

	}
	else if(nl->jtype == JMFFD)
	{
		// ... matrix-free finite-difference (MMFD)
		ierr = MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);
		ierr = MatMFFDSetBase(nl->MFFD, x, jr->gres);                                                      CHKERRQ(ierr);
		ierr = MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))JacApplyMFFD);                        CHKERRQ(ierr);
		ierr = MatShellSetContext(Amat, (void*)&nl->MFFD);                                                 CHKERRQ(ierr);
	}

	// assemble Jacobian & preconditioner
	ierr = MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Amat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Pmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyMFFD"
PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y)
{
	Mat *MFFD;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, (void**)&MFFD); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult((*MFFD), x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
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
