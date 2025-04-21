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
#include "LaMEM.h"
#include "matData.h"
#include "matrix.h"
#include "fdstag.h"
#include "tssolve.h"
#include "multigrid.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "JacRes.h"
#include "matFree.h"
#include "tools.h"
//---------------------------------------------------------------------------
PetscErrorCode NLSolCreate(SNES *p_snes, JacRes *jr)
{
	SNES            snes;
	KSP             ksp;
	Mat             J, P;
	PC              ipc;
	SNESLineSearch  ls;
	DOFIndex       *dof;
	SNESType        type;
	NLSol          *nl;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// allocate nonlinear solver context
	ierr = PetscMalloc(sizeof(NLSol), &nl); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(nl, sizeof(NLSol)); CHKERRQ(ierr);

	// access context
	dof = &(jr->fs->dof);

	// store context
	nl->jr = jr;

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &J); CHKERRQ(ierr);
	ierr = MatSetUp(J);                              CHKERRQ(ierr);

	// create matrix-free preconditioner operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &P); CHKERRQ(ierr);
	ierr = MatSetUp(P);                              CHKERRQ(ierr);

	// create finite-difference Jacobian
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, &nl->MFFD); CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(nl->MFFD,"fd_");       CHKERRQ(ierr);
	ierr = MatSetFromOptions(nl->MFFD);               CHKERRQ(ierr);
	ierr = MatSetUp(nl->MFFD);                        CHKERRQ(ierr);

	// create preconditioner and matrix
	ierr = PMatCreate    (&nl->pm, jr);     CHKERRQ(ierr);
	ierr = PCStokesCreate(&nl->pc, nl->pm); CHKERRQ(ierr);

	// setup nonlinear solver
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                        CHKERRQ(ierr);
	ierr = SNESSetApplicationContext(snes, (void*)nl);                 CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESNEWTONLS);                            CHKERRQ(ierr);
	ierr = SNESGetLineSearch(snes, &ls);                               CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);             CHKERRQ(ierr);
	ierr = SNESSetFunction(snes, jr->gres, &FormResidual, NULL);       CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes, J, P, &FormJacobian, NULL);           CHKERRQ(ierr);
	ierr = SNESSetConvergenceTest(snes, &SNESCoupledTest, NULL, NULL); CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);                                   CHKERRQ(ierr);

	// setup linear solver & preconditioner
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"js_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          CHKERRQ(ierr);

	// initialize Jacobian controls
	nl->jtype   = _PICARD_;
	nl->rtolPic = 1e-2;
	nl->nNwtIt  = 35;
	nl->rtolNwt = 1.1;

	// override from command line
	ierr = PetscOptionsGetScalar(NULL, NULL, "-snes_PicardSwitchToNewton_rtol", &nl->rtolPic, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (NULL, NULL, "-snes_NewtonSwitchToPicard_it",   &nl->nNwtIt,  NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-snes_NewtonSwitchToPicard_rtol", &nl->rtolNwt, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsHasName  (NULL, NULL, "-snes_picard_mat_free",           &nl->matFreePic);    CHKERRQ(ierr);

	// return solver
	(*p_snes) = snes;

	// check solver type compatibility with temperature diffsion activation
	ierr = SNESGetType(snes, &type); CHKERRQ(ierr);

	if(jr->ctrl.actTemp && !strcmp(type, SNESKSPONLY))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "act_temp_diff = 1 and -snes_type ksponly are incompatible, use -snes_max_it 1 instead\n");
	}

	// force one nonlinear iteration regardless of the initial residual
	ierr = SNESSetForceIteration(snes, PETSC_TRUE); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode NLSolDestroy(SNES *p_snes)
{
	NLSol *nl;
	Mat    J, P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = SNESGetApplicationContext((*p_snes), &nl); CHKERRQ(ierr);
	ierr = SNESGetJacobian((*p_snes), &J, &P, NULL, NULL); CHKERRQ(ierr);

	ierr = MatDestroy(&J);          CHKERRQ(ierr);
	ierr = MatDestroy(&P);          CHKERRQ(ierr);
	ierr = MatDestroy(&nl->MFFD);   CHKERRQ(ierr);
	ierr = PMatDestroy(nl->pm);     CHKERRQ(ierr);
	ierr = PCStokesDestroy(nl->pc); CHKERRQ(ierr);
	ierr = PetscFree(nl);           CHKERRQ(ierr);
	ierr = SNESDestroy(p_snes);     CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx)
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(ctx);

	ierr = SNESGetApplicationContext(snes, &nl); CHKERRQ(ierr);

	// access context
	jr = nl->jr;

	ierr = JacResFormResidual(jr, x, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	Vec         r;
	NLSol       *nl;
	PCStokes    pc;
	PMat        pm;
	JacRes      *jr;
	PetscInt    it;
	Controls   *ctrl;
	PetscScalar nrm;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(ctx);

	ierr = SNESGetApplicationContext(snes, &nl); CHKERRQ(ierr);

	// access context
	pc   =  nl->pc;
	pm   =  pc->pm;
	jr   =  nl->jr;
	ctrl = &jr->ctrl;

	//========================
	// Jacobian type selection
	//========================

	// get iteration counter and residual norm
	ierr = SNESGetIterationNumber(snes, &it);     CHKERRQ(ierr);
	ierr = SNESGetFunction(snes, &r, NULL, NULL); CHKERRQ(ierr);
	ierr = VecNorm(r, NORM_2, &nrm);              CHKERRQ(ierr);

	if(!nrm) nrm = 1.0;

	// initialize (always start with Picard)
	if(!it)
	{
		nl->it     = 0;
		nl->refRes = nrm;
		nl->jtype  = _PICARD_;
		nl->it_Nwt = 0;
	}
	else if(nl->jtype == _PICARD_)
	{
		// Picard case, check to switch to Newton
		if(nrm < nl->refRes*nl->rtolPic)
		{
			nl->jtype  = _MFFD_;
			nl->it_Nwt = 0;
		}
	}
	else if(nl->jtype == _MFFD_)
	{
		// Newton case, check to switch to Picard
		if(nrm > nl->refRes*nl->rtolNwt || nl->it_Nwt > (nl->nNwtIt-1))
		{
			nl->jtype = _PICARD_;
		}
	}

	// print info
	if(nl->jtype == _PICARD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3lld PICARD ||F||/||F0||=%e \n", (LLD)nl->it, nrm/nl->refRes);
	}
	else if(nl->jtype == _MFFD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3lld MMFD   ||F||/||F0||=%e \n", (LLD)nl->it, nrm/nl->refRes);
		nl->it_Nwt++;
	}

	// switch off pressure limit for plasticity after first iteration
	if(!ctrl->initGuess && it > 1)
	{
		ctrl->pLimPlast = 0;
	}

	// count iterations
	nl->it++;

	//=====================
	// setup preconditioner
	//=====================
	ierr = PMatAssemble(pm, jr);                                             CHKERRQ(ierr);
	ierr = PCStokesSetup(pc);                                                CHKERRQ(ierr);
	ierr = MatShellSetOperation(Pmat, MATOP_MULT, (void(*)(void))pc->Apply); CHKERRQ(ierr);
	ierr = MatShellSetContext(Pmat, pc);                                     CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY);                       CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Pmat, MAT_FINAL_ASSEMBLY);                       CHKERRQ(ierr);

	//===============
	// setup Jacobian
	//===============
	if(nl->jtype == _PICARD_)
	{
		if(nl->matFreePic)
		{
			// ... matrix-free Picard operator
			ierr = MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))JacApplyPicard); CHKERRQ(ierr);
			ierr = MatShellSetContext(Amat, (void*)jr);                                   CHKERRQ(ierr);
		}
		else
		{
			// ... assembled Picard operator
			ierr = MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))pm->Picard);     CHKERRQ(ierr);
			ierr = MatShellSetContext(Amat, pm->data);                                    CHKERRQ(ierr);
		}
	}
	else if(nl->jtype == _MFFD_)
	{
		// prepare MMFD operator
		ierr = MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);
		ierr = MatMFFDSetBase(nl->MFFD, x, jr->gres);                                                      CHKERRQ(ierr);
		ierr = MatMFFDSetType(nl->MFFD, MATMFFD_WP);                                                       CHKERRQ(ierr);

		// ... matrix-free finite-difference (MMFD)
		ierr = MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))JacApplyMFFD);                       CHKERRQ(ierr);
		ierr = MatShellSetContext(Amat, (void*)&nl->MFFD);                                                CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Amat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y)
{
	Mat *FD;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ierr = MatShellGetContext(A, (void**)&FD); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult((*FD), x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode SNESCoupledTest(
	SNES                snes,
	PetscInt            it,
	PetscReal           xnorm,
	PetscReal           gnorm,
	PetscReal           f,
	SNESConvergedReason *reason,
	void                *cctx)
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(cctx);

	ierr = SNESGetApplicationContext(snes, &nl); CHKERRQ(ierr);

	// access context
	jr = nl->jr;

	// call default convergence test
	ierr = SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, NULL); CHKERRQ(ierr);

	//=============================
	// Temperature diffusion solver
	//=============================

	if(!it) PetscFunctionReturn(0);

	if(jr->ctrl.actTemp)
	{
		ierr = JacResGetTempRes(jr, jr->ts->dt);            CHKERRQ(ierr);
		ierr = JacResGetTempMat(jr, jr->ts->dt);            CHKERRQ(ierr);
		ierr = KSPSetOperators(jr->tksp, jr->Att, jr->Att); CHKERRQ(ierr);
		ierr = KSPSetUp(jr->tksp);                          CHKERRQ(ierr);
		ierr = KSPSolve(jr->tksp, jr->ge, jr->dT);          CHKERRQ(ierr);
		ierr = JacResUpdateTemp(jr);                        CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode SNESPrintConvergedReason(SNES snes, PetscLogDouble t_beg)
{
	PetscLogDouble      t_end;
	SNESConvergedReason reason;
	PetscInt            its;
	KSP                 ksp;
	KSPConvergedReason  ksp_reason;
	PetscInt            div_severe;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set flag for severe divergence
	div_severe = 0;

	PetscTime(&t_end);

	ierr = SNESGetIterationNumber(snes, &its);    CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(snes, &reason);  CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	if(reason < 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "**************   NONLINEAR SOLVER FAILED TO CONVERGE!   ****************** \n");
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	if(reason == SNES_CONVERGED_FNORM_ABS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : ||F|| < atol \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_FNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : ||F|| < rtol*||F_initial|| \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_SNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : Newton computed step size small; || delta x || < stol || x ||\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_ITS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_ITERATING)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : SNES_CONVERGED_ITERATING\n"); CHKERRQ(ierr);
	}

	// DIVERGENCE

	else if(reason == SNES_DIVERGED_FUNCTION_DOMAIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the new x location passed the function is not in the domain of F\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FUNCTION_COUNT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : too many function evaluations\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINEAR_SOLVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the linear solve failed\n"); CHKERRQ(ierr);

		// detect severe divergence reason
		ierr = SNESGetKSP(snes, &ksp);                  CHKERRQ(ierr);
		ierr = KSPGetConvergedReason(ksp, &ksp_reason); CHKERRQ(ierr);

		if(ksp_reason == KSP_DIVERGED_BREAKDOWN
		|| ksp_reason == KSP_DIVERGED_INDEFINITE_PC
		|| ksp_reason == KSP_DIVERGED_NANORINF
		|| ksp_reason == KSP_DIVERGED_INDEFINITE_MAT)
		{
			div_severe = 1;
		}
	}
	else if(reason == SNES_DIVERGED_FNORM_NAN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : residual norm is NAN\n"); CHKERRQ(ierr);

		div_severe = 1;
	}
	else if(reason == SNES_DIVERGED_MAX_IT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINE_SEARCH)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the line search failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_INNER)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the inner solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LOCAL_MIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : || J^T b || is small, implies converged to local minimum of F\n"); CHKERRQ(ierr);
	}

	// stop if severe divergence reason detected
	if(div_severe)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Severe divergence reason detected (see above)");
	}

	PetscPrintf(PETSC_COMM_WORLD, "Number of iterations    : %lld\n", (LLD)its);

	PetscPrintf(PETSC_COMM_WORLD, "SNES solution time      : %g (sec)\n", t_end - t_beg);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

