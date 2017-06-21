/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   nlsolve.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "matrix.h"
#include "fdstag.h"
#include "tssolve.h"
#include "multigrid.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "JacRes.h"
#include "matFree.h"
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
	PetscBool       flg;

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

	ierr = SNESSetConvergenceTest(snes, &SNESCoupledTest, nl, NULL); CHKERRQ(ierr);

	// initialize Jacobian controls
	nl->jtype   = _PICARD_;
	nl->nPicIt  = 5;
	nl->rtolPic = 1e-2;
	nl->nNwtIt  = 35;
	nl->rtolNwt = 1.1;

	// override from command line
	ierr = PetscOptionsGetInt   (NULL, NULL, "-snes_Picard_max_it",             &nl->nPicIt, &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-snes_PicardSwitchToNewton_rtol", &nl->rtolPic,&flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (NULL, NULL, "-snes_NewtonSwitchToPicard_it",   &nl->nNwtIt, &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-snes_NewtonSwitchToPicard_rtol", &nl->rtolNwt, &flg); CHKERRQ(ierr);

	// set initial guess
	ierr = VecSet(jr->gsol, 0.0); CHKERRQ(ierr);

	// return solver
	(*p_snes) = snes;

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

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

	ierr = JacResFormResidual(jr, x, f); CHKERRQ(ierr);

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
	Controls   *ctrl;
	PetscScalar nrm;

	// clear unused parameters
	if(Amat) Amat = NULL;
	if(Pmat) Pmat = NULL;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	nl   = (NLSol*)ctx;
	pc   =  nl->pc;
	pm   =  pc->pm;
	jr   =  pm->jr;
	ctrl = &jr->ctrl;

    it_newton = 0;

    //========================
	// Jacobian type selection
	//========================

	// get iteration counter and residual norm
	ierr = SNESGetIterationNumber(snes, &it);     CHKERRQ(ierr);
	ierr = SNESGetFunction(snes, &r, NULL, NULL); CHKERRQ(ierr);
	ierr = VecNorm(r, NORM_2, &nrm);              CHKERRQ(ierr);
    
	if(!nrm) nrm = 1.0;

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
		if(nrm < nl->refRes*nl->rtolPic)
		{
			if(ctrl->jac_mat_free)
			{
				nl->jtype = _MF_;
			}
			else
			{
				nl->jtype = _MFFD_;
			}
		}
	}
	else if(nl->jtype == _MFFD_)
	{
		// Newton case, check to switch to Picard
		if(nrm > nl->refRes*nl->rtolNwt || it_newton > nl->nNwtIt)
		{
			nl->jtype = _PICARD_;
		}
	}

	if(jr->ts->istep < 2 && nl->it == 0)
	{
		// During the first and second timestep of a simulation, always start with picard iterations
		// that is important as plasticity is only activated during the second timestep, whereas the code might have
		// switched to MFFD already during the first timestep (and that solution is quite far off the plastic solution).
		nl->jtype = _PICARD_;
	}

	// print info
	if(nl->jtype == _PICARD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3lld PICARD ||F||/||F0||=%e \n", (LLD)nl->it, nrm/nl->refRes);
	}
	else if(nl->jtype == _MFFD_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3lld MMFD   ||F||/||F0||=%e \n", (LLD)nl->it, nrm/nl->refRes);
		it_newton++;
	}
	else if(nl->jtype == _MF_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%3lld MF     ||F||/||F0||=%e \n", (LLD)nl->it, nrm/nl->refRes);
		it_newton++;
	}

	// switch off pressure limit for plasticity after first iteration
	if(!ctrl->initGuess && it > 1)
	{
		ctrl->pLimPlast = 0;
	}

	// count iterations
	nl->it++;

	// setup preconditioner
	ierr = PMatAssemble(pm);                                                  CHKERRQ(ierr);
	ierr = PCStokesSetup(pc);                                                 CHKERRQ(ierr);
	ierr = MatShellSetOperation(nl->P, MATOP_MULT, (void(*)(void))pc->Apply); CHKERRQ(ierr);
	ierr = MatShellSetContext(nl->P, pc);                                     CHKERRQ(ierr);

	// setup Jacobian
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
	else if(nl->jtype == _MF_)
	{
		// ... matrix-free closed-form
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))JacApplyJacobian); CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, (void*)jr);                                     CHKERRQ(ierr);
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
PetscErrorCode SNESPrintConvergedReason(SNES snes, 	PetscLogDouble t_beg)
{
	PetscLogDouble      t_end;
	SNESConvergedReason reason;
	PetscInt            its;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
	else if(reason == SNES_CONVERGED_TR_DELTA)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Convergence Reason : SNES_CONVERGED_TR_DELTA\n"); CHKERRQ(ierr);
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
	}
	else if(reason == SNES_DIVERGED_FNORM_NAN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "SNES Divergence Reason  : residual norm is NAN\n"); CHKERRQ(ierr);
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

	PetscPrintf(PETSC_COMM_WORLD, "Number of iterations    : %lld\n", (LLD)its);

	PetscPrintf(PETSC_COMM_WORLD, "SNES solution time      : %g (sec)\n", t_end - t_beg);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESCoupledTest"
PetscErrorCode SNESCoupledTest(
	SNES                snes,
	PetscInt            it,
	PetscReal           xnorm,
	PetscReal           gnorm,
	PetscReal           f,
	SNESConvergedReason *reason,
	void                *cctx)
{

//	PetscErrorCode  SNESGetFunction(SNES snes,Vec *r,PetscErrorCode (**f)(SNES,Vec,Vec,void*),void **ctx)
//	PetscErrorCode  SNESGetSolution(SNES snes,Vec *x)

	// currently just calls temperature diffusion solver
	// together with standard convergence test
	// later should include temperature convergence test as well

	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
   	nl = (NLSol*)cctx;
   	jr = nl->pc->pm->jr;

	// call default convergence test
	ierr = SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, NULL); CHKERRQ(ierr);

	//=============================
	// Temperature diffusion solver
	//=============================

	if(!it) PetscFunctionReturn(0);

    if(jr->ctrl.actTemp)
    {
    	ierr = JacResGetTempRes(jr);                        CHKERRQ(ierr);
    	ierr = JacResGetTempMat(jr);                        CHKERRQ(ierr);
    	ierr = KSPSetOperators(jr->tksp, jr->Att, jr->Att); CHKERRQ(ierr);
    	ierr = KSPSetUp(jr->tksp);                          CHKERRQ(ierr);
    	ierr = KSPSolve(jr->tksp, jr->ge, jr->dT);          CHKERRQ(ierr);
    	ierr = JacResUpdateTemp(jr);                        CHKERRQ(ierr);
     }

	PetscFunctionReturn(0);
}
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
*/
//---------------------------------------------------------------------------
