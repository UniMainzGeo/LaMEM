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
	PC              pc;
	SNESLineSearch  ls;
	DOFIndex       *dof;
	SNESType        type;
	NLSol          *nl;
	PetscBool       flag;

	
	PetscFunctionBeginUser;

	// create nonlinear solver context
	PetscCall(PetscMalloc(sizeof(NLSol), &nl));
	PetscCall(PetscMemzero(nl, sizeof(NLSol)));

	// store residual context context
	nl->jr = jr;

	// access context
	dof = &(jr->fs->dof);

	// create matrix-free Jacobian operator
	PetscCall(MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &J));
	PetscCall(MatSetUp(J));

	// set Jacobian application operation
	PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))JacApply));

	// create matrix-free preconditioner operator
	PetscCall(MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &P));
	PetscCall(MatSetUp(P));

	// create Picard Jacobian
	PetscCall(MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->PICARD));
	PetscCall(MatSetUp(nl->PICARD));

	// create finite-difference Jacobian
	PetscCall(MatCreateMFFD(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, &nl->MFFD));
	PetscCall(MatSetOptionsPrefix(nl->MFFD,"fd_"));
	PetscCall(MatSetFromOptions(nl->MFFD));
	PetscCall(MatSetUp(nl->MFFD));

	// setup nonlinear solver
	PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
	PetscCall(SNESSetApplicationContext(snes, (void*)nl));
	PetscCall(SNESSetType(snes, SNESNEWTONLS));
	PetscCall(SNESGetLineSearch(snes, &ls));
	PetscCall(SNESLineSearchSetType(ls, SNESLINESEARCHBASIC));
	PetscCall(SNESSetFunction(snes, jr->gres, &FormResidual, NULL));
	PetscCall(SNESSetJacobian(snes, J, P, &FormJacobian, NULL));
	PetscCall(SNESSetConvergenceTest(snes, &SNESCoupledTest, NULL, NULL));
	PetscCall(SNESSetFromOptions(snes));

	// setup linear solver & preconditioner
	PetscCall(SNESGetKSP(snes, &ksp));
	PetscCall(KSPSetOptionsPrefix(ksp,"js_"));
	PetscCall(KSPSetFromOptions(ksp));
	PetscCall(KSPGetPC(ksp, &pc));
	PetscCall(PCSetType(pc, PCMAT));

	// create preconditioner context
	PetscCall(PCDataCreate(&nl->pc, jr, nl->PICARD, P));

	// initialize Jacobian controls
	nl->jtype    = _PICARD_;
	nl->rtolPic  = 1e-2;
	nl->minItPic = 5;
	nl->rtolNwt  = 1.2;
	nl->maxItNwt = 20;

	// override from command line
	PetscCall(PetscOptionsGetScalar(NULL, NULL, "-snes_picard_rtol",  &nl->rtolPic,  NULL));
	PetscCall(PetscOptionsGetInt   (NULL, NULL, "-snes_picard_minit", &nl->minItPic, NULL));
	PetscCall(PetscOptionsGetScalar(NULL, NULL, "-snes_newton_rtol",  &nl->rtolNwt,  NULL));
	PetscCall(PetscOptionsGetInt   (NULL, NULL, "-snes_newton_maxit", &nl->maxItNwt, NULL));

	// return solver
	(*p_snes) = snes;

	// check solver type compatibility with temperature diffsion activation
	PetscCall(SNESGetType(snes, &type));

	if(jr->ctrl.actTemp && !strcmp(type, SNESKSPONLY))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "act_temp_diff = 1 and -snes_type ksponly are incompatible, use -snes_max_it 1 instead\n");
	}

	// get automatic absolute tolerance initialization flags
	PetscCall(PetscOptionsHasName(NULL, NULL, "-snes_atol_auto",   &flag)); if(flag) { nl->snes_atol_auto   = 1; }
	PetscCall(PetscOptionsHasName(NULL, NULL, "-js_ksp_atol_auto", &flag)); if(flag) { nl->js_ksp_atol_auto = 1; }
	PetscCall(PetscOptionsHasName(NULL, NULL, "-ts_ksp_atol_auto", &flag)); if(flag) { nl->ts_ksp_atol_auto = 1; }

	// force one nonlinear iteration regardless of the initial residual
	PetscCall(SNESSetForceIteration(snes, PETSC_TRUE));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode NLSolDestroy(SNES *p_snes)
{
	NLSol *nl;
	Mat    J, P;

	
	PetscFunctionBeginUser;

	PetscCall(SNESGetApplicationContext((*p_snes), &nl));
	PetscCall(SNESGetJacobian((*p_snes), &J, &P, NULL, NULL));

	PetscCall(MatDestroy(&J));
	PetscCall(MatDestroy(&P));
	PetscCall(MatDestroy(&nl->MFFD));
	PetscCall(MatDestroy(&nl->PICARD));
	PetscCall(PCDataDestroy(&nl->pc));
	PetscCall(PetscFree(nl));
	PetscCall(SNESDestroy(p_snes));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx)
{
	NLSol  *nl;
	JacRes *jr;

	
	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(ctx);

	PetscCall(SNESGetApplicationContext(snes, &nl));

	// access context
	jr = nl->jr;

	PetscCall(JacResFormResidual(jr, x, f));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	Vec         r;
	NLSol       *nl;
	JacRes      *jr;
	PetscInt    it;
	Controls    *ctrl;
	PetscScalar nrm;

	
	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(ctx);

	PetscCall(SNESGetApplicationContext(snes, &nl));

	// access context
	jr   =  nl->jr;
	ctrl = &jr->ctrl;

	//========================
	// Jacobian type selection
	//========================

	// get iteration counter and residual norm
	PetscCall(SNESGetIterationNumber(snes, &it));
	PetscCall(SNESGetFunction(snes, &r, NULL, NULL));
	PetscCall(VecNorm(r, NORM_2, &nrm));

	if(!nrm) nrm = 1.0;

	// initialize (always start with Picard)
	if(!it)
	{
		nl->it     = 0;
		nl->refRes = nrm;
		nl->jtype  = _PICARD_;
		nl->itNwt  = 0;
	}
	else if(nl->jtype == _PICARD_)
	{
		// Picard case, check to switch to Newton (convergence)
		if(nrm < nl->refRes*nl->rtolPic)
		{
			nl->jtype = _MFFD_;
			nl->itNwt = 0;
		}
	}
	else if(nl->jtype == _MFFD_)
	{
		// Newton case, check to switch to Picard (divergence)
		if(nrm > nl->refRes*nl->rtolNwt || nl->itNwt > (nl->maxItNwt-1))
		{
			nl->jtype = _PICARD_;
		}
	}

	// force minimum number of Picard iterations
	if(nl->it < nl->minItPic)
	{
		nl->jtype = _PICARD_;
	}

	// print info
	if(nl->jtype == _PICARD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3" PetscInt_FMT " PICARD ||F||/||F0||=%e \n", nl->it, nrm/nl->refRes);
	}
	else if(nl->jtype == _MFFD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3" PetscInt_FMT " MMFD   ||F||/||F0||=%e \n", nl->it, nrm/nl->refRes);
		nl->itNwt++;
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
	PetscCall(PCDataSetup(&nl->pc, jr));

	// assembly preconditioner
	PetscCall(MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY));
	PetscCall(MatAssemblyEnd  (Pmat, MAT_FINAL_ASSEMBLY));

	//===============
	// setup Jacobian
	//===============
	if(nl->jtype == _PICARD_)
	{
		// assembly Picard operator
		PetscCall(MatAssemblyBegin(nl->PICARD, MAT_FINAL_ASSEMBLY));
		PetscCall(MatAssemblyEnd  (nl->PICARD, MAT_FINAL_ASSEMBLY));

		// set Picard operator
		PetscCall(MatShellSetContext(Amat, (void*)&nl->PICARD));
	}
	else if(nl->jtype == _MFFD_)
	{
		// prepare MMFD operator
		PetscCall(MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes));
		PetscCall(MatMFFDSetBase(nl->MFFD, x, nl->jr->gres));
		PetscCall(MatMFFDSetType(nl->MFFD, MATMFFD_WP));

		// assembly MMFD operator
		PetscCall(MatAssemblyBegin(nl->MFFD, MAT_FINAL_ASSEMBLY));
		PetscCall(MatAssemblyEnd  (nl->MFFD, MAT_FINAL_ASSEMBLY));

		// set MMFD operator
		PetscCall(MatShellSetContext(Amat, (void*)&nl->MFFD));
	}

	// assembly Jacobian
	PetscCall(MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY));
	PetscCall(MatAssemblyEnd  (Amat, MAT_FINAL_ASSEMBLY));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacApply(Mat A, Vec x, Vec y)
{
	Mat *J;

	PetscFunctionBeginUser;

	// access context
	PetscCall(MatShellGetContext(A, (void**)&J));

	// compute Jacobian times vector product
	PetscCall(MatMult((*J), x, y));

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
	NLSol    *nl;
	JacRes   *jr;
	KSP      js_ksp;

	PetscScalar norm;

	PetscFunctionBeginUser;

	// clear unused parameters
	UNUSED(cctx);

	PetscCall(SNESGetApplicationContext(snes, &nl));

	// access context
	jr = nl->jr;

	PetscCall(SNESGetKSP(snes, &js_ksp));

	// update absolute tolerances
	PetscCall(SNESUpdateAbsTol(snes,    nl->snes_atol_auto,   nl->snes_ref_norm,   f, it));
	PetscCall(KSPUpdateAbsTol (js_ksp,  nl->js_ksp_atol_auto, nl->js_ksp_ref_norm, f, it));

	// call default convergence test
	PetscCall(SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, NULL));

	//=============================
	// Temperature diffusion solver
	//=============================

	if(jr->ctrl.actTemp)
	{
		// get residual and tangent matrix
		PetscCall(JacResGetTempRes(jr, jr->ts->dt));
		PetscCall(JacResGetTempMat(jr, jr->ts->dt));

		// update absolute tolerance
		PetscCall(VecNorm(jr->ge, NORM_2, &norm));
		PetscCall(NLSolvePushNorm(nl->ts_ksp_ref_norm, jr->ts_ksp_ref_norm, norm));
		PetscCall(KSPUpdateAbsTol(jr->tksp, nl->ts_ksp_atol_auto, nl->ts_ksp_ref_norm, norm, it));

		// compute and apply temperature correction (not on first iteration)
		if(it)
		{
			PetscCall(KSPSetOperators(jr->tksp, jr->Att, jr->Att));
			PetscCall(KSPSetUp(jr->tksp));
			PetscCall(KSPSolve(jr->tksp, jr->ge, jr->dT));
			PetscCall(JacResUpdateTemp(jr));
		}
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

	
	PetscFunctionBeginUser;

	// set flag for severe divergence
	div_severe = 0;

	PetscTime(&t_end);

	PetscCall(SNESGetIterationNumber(snes, &its));
	PetscCall(SNESGetConvergedReason(snes, &reason));

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	if(reason < 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "**************   NONLINEAR SOLVER FAILED TO CONVERGE!   ****************** \n");
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	if(reason == SNES_CONVERGED_FNORM_ABS)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : ||F|| < atol \n");
	}
	else if(reason == SNES_CONVERGED_FNORM_RELATIVE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : ||F|| < rtol*||F_initial|| \n");
	}
	else if(reason == SNES_CONVERGED_SNORM_RELATIVE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : Newton computed step size small; || delta x || < stol || x ||\n");
	}
	else if(reason == SNES_CONVERGED_ITS)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : maximum iterations reached\n");
	}
	else if(reason == SNES_CONVERGED_ITERATING)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : SNES_CONVERGED_ITERATING\n");
	}

	// DIVERGENCE

	else if(reason == SNES_DIVERGED_FUNCTION_DOMAIN)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the new x location passed the function is not in the domain of F\n");
	}
	else if(reason == SNES_DIVERGED_FUNCTION_COUNT)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : too many function evaluations\n");
	}
	else if(reason == SNES_DIVERGED_LINEAR_SOLVE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the linear solve failed\n");

		// detect severe divergence reason
		PetscCall(SNESGetKSP(snes, &ksp));
		PetscCall(KSPGetConvergedReason(ksp, &ksp_reason));

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
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : residual norm is NAN\n");

		div_severe = 1;
	}
	else if(reason == SNES_DIVERGED_MAX_IT)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : maximum iterations reached\n");
	}
	else if(reason == SNES_DIVERGED_LINE_SEARCH)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the line search failed\n");
	}
	else if(reason == SNES_DIVERGED_INNER)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : the inner solve failed\n");
	}
	else if(reason == SNES_DIVERGED_LOCAL_MIN)
	{
		PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : || J^T b || is small, implies converged to local minimum of F\n");
	}

	// stop if severe divergence reason detected
	if(div_severe)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Severe divergence reason detected (see above)");
	}

	PetscPrintf(PETSC_COMM_WORLD, "Number of iterations    : %" PetscInt_FMT "\n", its);

	PetscPrintf(PETSC_COMM_WORLD, "SNES solution time      : %g (sec)\n", t_end - t_beg);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode SNESUpdateAbsTol(SNES snes, PetscInt set, PetscScalar &refNorm, PetscScalar norm, PetscInt it)
{
	PetscScalar rtol, abstol;

	PetscFunctionBeginUser;

	// only update tolerance if all criteria are satisfied:
	// 1 - tolerance setting is requested
	// 2 - first iteration of the SNES solve
	// 3 - current norm is larger than the reference

	if(!set || it || norm < refNorm) PetscFunctionReturn(0);

	// update reference norm
	refNorm = norm;

	// get relative tolerance
	PetscCall(SNESGetTolerances(snes, NULL, &rtol, NULL, NULL, NULL));

	// compute new absolute tolerance
	abstol = refNorm*rtol;

	// update absolute tolerance
	PetscCall(SNESSetTolerances(snes, abstol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode KSPUpdateAbsTol(KSP ksp, PetscInt set, PetscScalar &refNorm, PetscScalar norm, PetscInt it)
{
	PetscScalar rtol, abstol;

	PetscFunctionBeginUser;

	// only update tolerance if all criteria are satisfied:
	// 1 - tolerance setting is requested
	// 2 - first iteration of the SNES solve
	// 3 - current norm is larger than the reference

	if(!set || it || norm < refNorm) PetscFunctionReturn(0);

	// update reference norm
	refNorm = norm;

	// get relative tolerance
	PetscCall(KSPGetTolerances(ksp, &rtol, NULL, NULL, NULL));

	// compute new absolute tolerance
	abstol = refNorm*rtol;

	// update absolute tolerance
	PetscCall(KSPSetTolerances(ksp, PETSC_CURRENT, abstol, PETSC_CURRENT, PETSC_CURRENT));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode NLSolvePushNorm(PetscScalar ref_norm, PetscScalar ref_norm_init, PetscScalar &norm)
{
	// override norm with initial value if reference is not initialized

	PetscFunctionBeginUser;

	if(!ref_norm)
	{
		if(norm < ref_norm_init)
		{
			norm = ref_norm_init;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


