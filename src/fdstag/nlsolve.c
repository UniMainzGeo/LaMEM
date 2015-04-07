//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "Utils.h"
#include "tools.h"

//---------------------------------------------------------------------------
// * add bound checking for iterative solution vector in SNES
// * automatically set -snes_type ksponly (for linear problems)
// * add line search (PETSc) and whatever load control methods (arc-length?)
// * closed-form matrix-free Jacobian
// * residual function scaling
// * adaptive setting of absolute tolerance based on previous steps residual norms
//   (also for linear solves)
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolClear"
PetscErrorCode NLSolClear(NLSol *nl)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(nl, sizeof(NLSol)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolCreate"
PetscErrorCode NLSolCreate(NLSol *nl, PCStokes pc, SNES *p_snes)
{
	SNES            snes;
	KSP             ksp;
	PC              ipc;
	SNESLineSearch  ls;
	JacRes         *jr;
	DOFIndex       *dof;
	PetscBool       flg, useCustomTest=PETSC_FALSE,custom_ksp_monitor=PETSC_FALSE;

    PetscErrorCode ierr;
    PetscFunctionBegin;

	// store context
 	nl->pc = pc;

 	// access context
	jr  = pc->pm->jr;
	dof = &(jr->fs->dof);

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->J); CHKERRQ(ierr);
	ierr = MatSetUp(nl->J);                              CHKERRQ(ierr);

	// create matrix-free Preconditioner operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->P); CHKERRQ(ierr);
	ierr = MatSetUp(nl->P);                              CHKERRQ(ierr);

	// create finite-difference Jacobian
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, &nl->MFFD); CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(nl->MFFD,"fd_");       CHKERRQ(ierr);
	ierr = MatSetFromOptions(nl->MFFD);               CHKERRQ(ierr);
	ierr = MatSetUp(nl->MFFD);                        CHKERRQ(ierr);

	// setup nonlinear solver
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                     CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESNEWTONLS);                         CHKERRQ(ierr);
	ierr = SNESGetLineSearch(snes, &ls);                            CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);          CHKERRQ(ierr);
	ierr = SNESSetFunction(snes, jr->gres, &FormResidual, nl);      CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes, nl->J, nl->P, &FormJacobian, nl);  CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);                                CHKERRQ(ierr);

	// setup linear solver & preconditioner
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"js_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          CHKERRQ(ierr);


	// activate custom test for linear iterations?
	ierr = PetscOptionsGetBool(NULL,"-js_use_custom_test",&useCustomTest,&flg); CHKERRQ(ierr);

	if ( useCustomTest )
	{
		PetscPrintf(PETSC_COMM_WORLD,"Using custom test function for residuals\n");
		nl->wsCtx.epsfrac    = 0.005; // set default to 1 percent
		nl->wsCtx.eps        = 0.0;
		nl->wsCtx.rnorm_init = 0.0;
		nl->wsCtx.winwidth   = 20;

		ierr = PetscOptionsGetInt( PETSC_NULL,"-js_ksp_difftol_winwidth",&nl->wsCtx.winwidth,PETSC_NULL ); CHKERRQ(ierr);
		ierr = PetscOptionsGetReal( PETSC_NULL,"-js_ksp_difftol_eps"    ,&nl->wsCtx.epsfrac,PETSC_NULL ); CHKERRQ(ierr);
		ierr = PetscOptionsGetBool( PETSC_NULL, "-js_custom_ksp_monitor", &custom_ksp_monitor, &flg ); CHKERRQ(ierr);
		
		if ( (nl->wsCtx.winwidth < 1) | (nl->wsCtx.winwidth > _max_win_size_))
		{
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "js_ksp_difftol_winwidth should be between to be between 1 and %lld\n",(LLD) _max_win_size_);
		}
		if ( (nl->wsCtx.epsfrac < 0.0) | (nl->wsCtx.epsfrac > 1.0) ) 
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "js_ksp_difftol_eps should be chosen to be between 0 and 1\n");
		}
		
		
		// set residual monitors
		if(custom_ksp_monitor)
		{
			ierr = KSPMonitorSet(ksp, &KSPWinStopMonitor,&nl->wsCtx, NULL); CHKERRQ(ierr);
		}
	
		
/*
		ierr = PetscPrintf( PETSC_COMM_WORLD, "Stopping conditions: \n");
		PetscPrintf( PETSC_COMM_WORLD, "rtol : %14.12e\n",(double)ctx->rtol);
		PetscPrintf( PETSC_COMM_WORLD, "atol : %14.12e\n",(double)ctx->atol);
		PetscPrintf( PETSC_COMM_WORLD, "maxit: %D\n",ctx->maxits);
		PetscPrintf( PETSC_COMM_WORLD, "difftol_eps: %14.12e\n",nl->wsCtx.epsfrac);
		PetscPrintf( PETSC_COMM_WORLD, "difftol_winwidth: %D\n",nl->wsCtx.winwidth);
*/
//		ierr = KSPStopCondConfig(ksp, &nl->wsCtx); CHKERRQ(ierr);
		ierr = KSPSetConvergenceTest(ksp, &KSPWinStopTest, &nl->wsCtx, NULL);CHKERRQ(ierr);
	}

//	ierr = SNESSetConvergenceTest(snes, SNESBlockStopTest, &nlctx, NULL); CHKERRQ(ierr);

	// initialize Jacobian controls
	nl->jtype   = _PICARD_;
	nl->nPicIt  = 5;
	nl->rtolPic = 1e-2;
	nl->nNwtIt  = 35;
	nl->rtolNwt  = 1.1;

	// override from command line
	ierr = PetscOptionsGetInt   (PETSC_NULL, "-snes_Picard_max_it",             &nl->nPicIt, &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-snes_PicardSwitchToNewton_rtol", &nl->rtolPic,&flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (PETSC_NULL, "-snes_NewtonSwitchToPicard_it",   &nl->nNwtIt, &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-snes_NewtonSwitchToPicard_rtol", &nl->rtolNwt, &flg); CHKERRQ(ierr);

	// set initial guess
	ierr = VecSet(jr->gsol, 0.0); CHKERRQ(ierr);

	// return solver
	(*p_snes) = snes;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESActEW"
PetscErrorCode SNESActEW(SNES snes)
{
	// activate Eisenstat-Walker method

    PetscErrorCode ierr;
    PetscFunctionBegin;

	ierr = SNESKSPSetUseEW(snes, PETSC_TRUE); CHKERRQ(ierr);

	ierr = SNESKSPSetParametersEW(snes,
		PETSC_DEFAULT,
		1e-3,
		1e-2,
		PETSC_DEFAULT,
		PETSC_DEFAULT,
		PETSC_DEFAULT,
		PETSC_DEFAULT); CHKERRQ(ierr);
/*
	PetscInt  version;
	PetscReal rtol_0;
	PetscReal rtol_max;
	PetscReal pgamma;
	PetscReal alpha;
	PetscReal alpha2;
	PetscReal threshold;

	ierr = SNESKSPGetParametersEW(snes,
		&version,
		&rtol_0,
		&rtol_max,
		&pgamma,
		&alpha,
		&alpha2,
		&threshold); CHKERRQ(ierr);
*/

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
	jr = nl->pc->pm->jr;

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

	Vec         r;
	NLSol       *nl;
	PCStokes    pc;
	PMat        pm;
	JacRes      *jr;
	PetscInt    it, it_newton;
	PetscScalar nrm;

	// clear unused parameters
	if(Amat) Amat = NULL;
	if(Pmat) Pmat = NULL;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	nl = (NLSol*)ctx;
	pc = nl->pc;
	pm = pc->pm;
	jr = pm->jr;
    it_newton = 0;
    
	//========================
	// Jacobian type selection
	//========================

	// get iteration counter and residual norm
	ierr = SNESGetIterationNumber(snes, &it);     CHKERRQ(ierr);
	ierr = SNESGetFunction(snes, &r, NULL, NULL); CHKERRQ(ierr);
	ierr = VecNorm(r, NORM_2, &nrm);              CHKERRQ(ierr);

    
    // initialize
	if(!it)
	{
		nl->it     = 0;
		nl->refRes = nrm;
        nl->jtype = _PICARD_;
  	}
	else if(nl->jtype == _PICARD_)
	{
		// Picard case, check to switch to Newton
		//if(nrm < nl->refRes*nl->tolPic || nl->it > nl->nPicIt)
        if(nrm < (nl->refRes*nl->rtolPic) )
        {
			PetscPrintf(PETSC_COMM_WORLD,"        ***        \n");
			PetscPrintf(PETSC_COMM_WORLD,"SWITCH TO USING MMFD JACOBIAN; ||F||/||F0||=%e,  PicIt=%i\n",nrm/nl->refRes,nl->nPicIt);
			PetscPrintf(PETSC_COMM_WORLD,"        ***        \n");

			nl->jtype = _MFFD_;
		}
	}
	else if(nl->jtype == _MFFD_)
	{
		// Newton case, check to switch to Picard
		if(nrm > nl->refRes*nl->rtolNwt || it_newton > nl->nNwtIt)
		{
			PetscPrintf(PETSC_COMM_WORLD,"        ***          \n");
			PetscPrintf(PETSC_COMM_WORLD,"SWITCH TO USING PICARD JACOBIAN  ||F||/||F0||=%e,  PicIt=%i\n \n",nrm/nl->refRes,nl->nNwtIt);
			PetscPrintf(PETSC_COMM_WORLD,"        ***          \n");

			nl->jtype = _PICARD_;
		}
	}
   
    
    if ((JacResGetStep(jr)<2) & (nl->it == 0)){
        // During the first and second timestep of a simulation, always start with picard iterations
        // that is important as plasticity is only activated during the second timestep, whereas the code might have
        // switched to MFFD already during the first timestep (and that solution is quite far off the plastic solution).
        nl->jtype = _PICARD_;
    }

    /* print info*/
    if (nl->jtype == _PICARD_)
    {
        PetscPrintf(PETSC_COMM_WORLD,"USING PICARD JACOBIAN for iteration %i,  ||F||/||F0||=%e \n",nl->it, nrm/nl->refRes);
    }
    else if(nl->jtype == _MFFD_)
    {
       PetscPrintf(PETSC_COMM_WORLD,"USING MMFD JACOBIAN for iteration %i,    ||F||/||F0||=%e \n",nl->it, nrm/nl->refRes);
        it_newton = it_newton+1;
    }
    
    
    
	// count iterations
	nl->it++;

	//========================

	// setup preconditioner
	ierr = PMatAssemble(pm);                                                  CHKERRQ(ierr);
	ierr = PCStokesSetup(pc);                                                 CHKERRQ(ierr);
	ierr = MatShellSetOperation(nl->P, MATOP_MULT, (void(*)(void))pc->Apply); CHKERRQ(ierr);
	ierr = MatShellSetContext(nl->P, pc);                                     CHKERRQ(ierr);

	// setup Jacobian ...
	if(nl->jtype == _PICARD_)
	{
		// ... Picard
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))pm->Picard); CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, pm->data);                                CHKERRQ(ierr);
	}
	else if(nl->jtype == _MFFD_)
	{
		// ... matrix-free finite-difference (MMFD)
		ierr = MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);
		ierr = MatMFFDSetBase(nl->MFFD, x, jr->gres);                                                      CHKERRQ(ierr);
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))JacApplyMFFD);                       CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, (void*)&nl->MFFD);                                                CHKERRQ(ierr);
	}

	// assemble Jacobian & preconditioner
	ierr = MatAssemblyBegin(nl->P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (nl->P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(nl->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (nl->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyMFFD"
PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y)
{
	Mat *FD;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, (void**)&FD); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult((*FD), x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESPrintConvergedReason"
PetscErrorCode SNESPrintConvergedReason(SNES snes)
{
	SNESConvergedReason reason;
	PetscInt            its;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = SNESGetIterationNumber(snes, &its);    CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(snes, &reason);  CHKERRQ(ierr);

    // CONVERGENCE

    if(reason == SNES_CONVERGED_FNORM_ABS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: ||F|| < atol \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_FNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: ||F|| < rtol*||F_initial|| \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_SNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: Newton computed step size small; || delta x || < stol || x ||\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_ITS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_TR_DELTA)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: SNES_CONVERGED_TR_DELTA\n"); CHKERRQ(ierr);
	}

    // DIVERGENCE

	else if(reason == SNES_DIVERGED_FUNCTION_DOMAIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the new x location passed the function is not in the domain of F\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FUNCTION_COUNT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: too many function evaluations\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINEAR_SOLVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the linear solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FNORM_NAN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: residual norm is NAN\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_MAX_IT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINE_SEARCH)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the line search failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_INNER)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the inner solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LOCAL_MIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: || J^T b || is small, implies converged to local minimum of F\n"); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of iterations : %lld\n", (LLD)its);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "CheckVelocityError"
PetscErrorCode CheckVelocityError(UserCtx *user)
{
	PetscScalar MaxVel, MinVel;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Error detection
	ierr = VecMax(user->sol_advect, PETSC_NULL, &MaxVel);	CHKERRQ(ierr);
	ierr = VecMin(user->sol_advect, PETSC_NULL, &MinVel); CHKERRQ(ierr);
	MaxVel = PetscMax(MaxVel, PetscAbsScalar(MinVel));

	if(isnan(MaxVel))
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "  *** Emergency stop! Maximum velocity is NaN ***  \n");
	}

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
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "KSPWinStopTest"
PetscErrorCode KSPWinStopTest(KSP ksp, PetscInt thisit, PetscScalar thisnorm, KSPConvergedReason *reason, void *mctx)
{
	PetscErrorCode    ierr;
	PetscInt          i=0,inow=0,maxits;
	PetscScalar       rtol,atol,dtol,ttol;
	WinStopCtx       *winstop = (WinStopCtx*) mctx;
	PetscBool         winnorm =  PETSC_FALSE;


	// initialise converged reason
	*reason = KSP_CONVERGED_ITERATING;

	// get tolerances
	ierr = KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits); CHKERRQ(ierr);

	// store current and last norm
	i = thisit % 2;
	winstop->rnorms[i] = thisnorm;

	// only evaluate the norm in the second iteration
	if ( thisit > 0 )
	{
		// set initial norm (mean of frist and second residual) and epsilon
		if (thisit == 1)
		{
			winstop->rnorm_init = 0.5 * (winstop->rnorms[0] + winstop->rnorms[1]);
			winstop->eps        = winstop->epsfrac * winstop->rnorm_init;
//          Maybe, epsilon should be independent of the initial residual and a rather small number of the order 1e-4
//			winstop->eps        = winstop->epsfrac;
		}

		// compute the total tolerance
		ttol = PetscMax(rtol * winstop->rnorm_init, atol);

		// compute residual difference and store it in the running window
		inow = thisit % (winstop->winwidth-1);
 		if (inow == 0) inow = winstop->winwidth-1;
 		else           inow = inow-1;
		winstop->rnormdiffs[inow] = fabs(winstop->rnorms[1]- winstop->rnorms[0]);

		// compute the criterion as soon as we have enough iterations
		if (thisit >= winstop->winwidth)
		{
//			winstop->diffnorm = getStdv(winstop->rnormdiffs, winstop->winwidth-1);
			winstop->diffnorm = getVar(winstop->rnormdiffs, winstop->winwidth-1);
			winnorm  = PETSC_TRUE;
		}
		else
		{
			winstop->diffnorm = 0.0;
			winnorm  = PETSC_FALSE;
		}



		// --- Check norms ---

		// problems?
		if (PetscIsInfOrNanScalar(thisnorm))
		{
			PetscInfo(ksp,"Linear solver has created a not a number (NaN) as the residual norm, declaring divergence \n");
			*reason = KSP_DIVERGED_NANORINF;
		}

		// ttol
		else if (thisnorm <= ttol)
		{
			// atol
			if (thisnorm < atol)
			{
			  PetscInfo3(ksp,"Linear solver has converged. Residual norm %14.12e is less than absolute tolerance %14.12e at iteration %D\n",(double)thisnorm,(double)atol,thisit);
			  *reason = KSP_CONVERGED_ATOL;
			}

			// rtol
			else
			{
			  if (winstop->rnorm_init)
			  {
				PetscInfo4(ksp,"Linear solver has converged. Residual norm %14.12e is less than relative tolerance %14.12e times initial residual norm %14.12e at iteration %D\n",(double)thisnorm,(double)rtol,(double)winstop->rnorm_init,thisit);
			  }
			  else
			  {
				PetscInfo4(ksp,"Linear solver has converged. Residual norm %14.12e is less than relative tolerance %14.12e times initial right hand side norm %14.12e at iteration %D\n",(double)thisnorm,(double)rtol,(double)winstop->rnorm_init,thisit);
			  }
			  *reason = KSP_CONVERGED_RTOL;
			}
		}

		// difftol
		else if (winnorm & (winstop->diffnorm < winstop->eps) )
		{
			PetscInfo4(ksp,"Linear solver has converged. The standard deviation of the residual differences within the running window %14.12e is less than %g % of the initial right hand side norm %14.12e at iteration %D\n",(double)winstop->diffnorm,(double)winstop->eps,(double)winstop->rnorm_init,thisit);
			*reason = KSP_CONVERGED_ITS;
		}

		// maxits
		else if ( thisit == maxits+1 )
		{
			PetscInfo2(ksp,"Linear solver has converged. The maximum number of iterations %D has been reached at iteration %D\n",maxits,thisit);
			*reason = KSP_CONVERGED_ITS;
		}



		// divergence
		else if ( thisnorm >= (dtol * winstop->rnorm_init))
		{
			PetscInfo3(ksp,"Linear solver is diverging. Initial right hand size norm %14.12e, current residual norm %14.12e at iteration %D\n",(double)winstop->rnorm_init,(double)thisnorm,thisit);
			*reason = KSP_DIVERGED_DTOL;
		}

		// otherwise, continue iterations
		else
		{
			*reason = KSP_CONVERGED_ITERATING;
		}
	
	
	
	}
	else
	{
		*reason = KSP_CONVERGED_ITERATING;	
	}


	return(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "KSPWinStopMonitor"
PetscErrorCode KSPWinStopMonitor(KSP ksp,PetscInt thisit,PetscScalar thisnorm, void *mctx)
{
	PetscErrorCode    ierr;
	WinStopCtx       *winstop = (WinStopCtx*) mctx;
	
	// Output to screen
	if (thisit == 0 && ((PetscObject)ksp)->prefix)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norms for %sksp solve.\n",((PetscObject)ksp)->prefix); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  KSP_it       KSP_res      Var(dKSP_res(1..%0.2d))        Eps\n",winstop->winwidth); CHKERRQ(ierr);
	}
	else if (thisit < winstop->winwidth)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "   %0.3d   %14.12e         ---                 ---\n",thisit, thisnorm); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "   %0.3d   %14.12e  %14.12e  %14.12e  \n",thisit, thisnorm, winstop->diffnorm, winstop->eps); CHKERRQ(ierr);
	}


	return(0);
}
