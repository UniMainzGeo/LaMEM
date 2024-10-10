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
	PC              ipc;
	SNESLineSearch  ls;
	DOFIndex       *dof;
	PetscBool       flg;
	SNESType        type;
	NLSol          *nl;



	PetscBool    ksp_mat_free;






	PCStokes pc;

    PetscErrorCode ierr;
    PetscFunctionBeginUser;


	ierr = PetscMalloc(sizeof(NLSol), &nl); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(nl, sizeof(NLSol)); CHKERRQ(ierr);


	// create Stokes preconditioner, matrix and nonlinear solver
	ierr = PCStokesCreate(&pc, jr);     CHKERRQ(ierr);



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


	// PetscErrorCode SNESSetApplicationContext(SNES snes, void *usrP)



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


//

	ierr = PetscOptionsHasName(NULL, NULL, "-js_mat_free ", &ksp_mat_free); CHKERRQ(ierr);



	// return solver
	(*p_snes) = snes;

	// display specified solver options
	ierr = DisplaySpecifiedSolverOptions(pc, snes); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
 	PetscFunctionBeginUser;

	ierr = SNESGetApplicationContext((*p_snes), &nl); CHKERRQ(ierr);

 	ierr = MatDestroy(&nl->J);      CHKERRQ(ierr);
	ierr = MatDestroy(&nl->P);      CHKERRQ(ierr);
	ierr = MatDestroy(&nl->MFFD);   CHKERRQ(ierr);
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
	if(ctx) ctx = NULL;

	ierr = SNESGetApplicationContext(snes, &nl); CHKERRQ(ierr);

	// clear unused parameters
	if(snes) snes = NULL;

	// access context
	nl = (NLSol*)ctx;
	jr = nl->pc->pm->jr;

	ierr = JacResFormResidual(jr, x, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	// Compute FDSTAG Jacobian matrix and preconditioner

	Vec         r;
	NLSol       *nl;
	PCStokes    pc;
	PMat        pm;
	JacRes      *jr;
	PetscInt    it;
	Controls   *ctrl;
	PetscScalar nrm;

	// clear unused parameters
	if(Amat) Amat = NULL;
	if(Pmat) Pmat = NULL;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	nl   = (NLSol*)ctx;
	pc   =  nl->pc;
	pm   =  pc->pm;
	jr   =  pm->jr;
	ctrl = &jr->ctrl;

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
		nl->it_Nwt++;
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
		// ... matrix-free Picard
//		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))JacApplyPicard); CHKERRQ(ierr);
//		ierr = MatShellSetContext(nl->J, (void*)jr);                                   CHKERRQ(ierr);

		// ... Picard
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))pm->Picard);     CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, pm->data);                                    CHKERRQ(ierr);

	}
	else if(nl->jtype == _MFFD_)
	{
		// ... matrix-free finite-difference (MMFD)
		ierr = MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);
		ierr = MatMFFDSetBase(nl->MFFD, x, jr->gres);                                                      CHKERRQ(ierr);
		ierr = MatMFFDSetType(nl->MFFD, MATMFFD_WP);                                                       CHKERRQ(ierr);
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
PetscErrorCode SNESPrintConvergedReason(SNES snes, 	PetscLogDouble t_beg)
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
	PetscFunctionBeginUser;

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
PetscErrorCode DisplaySpecifiedSolverOptions(PCStokes pc, SNES snes)
{
	PetscErrorCode 	ierr;
	KSP 			ksp_coarse, ksp_levels, ksp;
	PC 				pc_coarse, pc_levels;
	char      		pname[_str_len_];
	PCStokesMG 		*mg;
	FDSTAG          *fs;
	PCType 			pc_type;
	KSPType 	    ksp_type;
	PCStokesUser 	*user;
	PetscScalar 	scalar;
	PetscInt 		integer, refine_y;
	PetscBool		found;
	MatSolverType   solver_type;

 	PetscFunctionBeginUser;

	/* 	This routine prints the solver options that are specified on the command-line or in the PetscOptions of the LaMEM input script
	 	More complete options can be displayed with the command-line options
		-js_ksp_view
		and -in case of multigrid-
		-gmg_pc_view
		This however clutters the output a bit
	*/

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");


	// Report (linear and  nonlinear) solver options currently used
	// get solver context
	user = (PCStokesUser*)pc->data;

	PetscPrintf(PETSC_COMM_WORLD, "Solver parameters specified: \n");
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	KSPGetType(ksp, &ksp_type);
	PetscPrintf(PETSC_COMM_WORLD, "   Outermost Krylov solver       : %s \n", ksp_type);
	if (pc->type == _STOKES_MG_){

		mg 		= 	(PCStokesMG*)pc->data; // retrieve MG object
		ierr 	= 	PCMGGetSmoother(mg->mg.pc, 1, &ksp_levels); 	CHKERRQ(ierr);
		ierr 	= 	KSPGetPC(ksp_levels, &pc_levels);        	CHKERRQ(ierr);

		// Multigrid solver
		PetscPrintf(PETSC_COMM_WORLD, "   Solver type                   : multigrid \n");

		/* Do we have a 2D setup & only refine in x/z direction? */
		fs  = pc->pm->jr->fs;
		ierr = DMDAGetRefinementFactor(fs->DA_CEN, NULL, &refine_y, NULL); CHKERRQ(ierr);
		if (refine_y==1){
			PetscPrintf(PETSC_COMM_WORLD, "   Multigrid refinement in x/z \n");
		}

		/* Multigrid parameters for the smootheners at the various levels */
		ierr = PetscOptionsGetString(NULL, NULL,"-gmg_mg_levels_ksp_type", pname, _str_len_, &found); CHKERRQ(ierr);
		if (found){
			PetscPrintf(PETSC_COMM_WORLD, "   Multigrid smoother levels KSP : %s \n", pname);
		}
		else{ // default
			ierr = KSPGetType(ksp_levels, &ksp_type); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "   Multigrid smoother levels KSP : %s \n", ksp_type);
		}

		// depending on the smoother, there may be more options
		if (!strcmp(pname, "richardson")){
			// Options that go with the Richardson solver
			ierr = PetscOptionsGetScalar(NULL, NULL,"-gmg_mg_levels_ksp_richardson_scale", &scalar, &found); CHKERRQ(ierr);
			if (found){PetscPrintf(PETSC_COMM_WORLD, "   Multigrid dampening parameter : %f \n", scalar); }
		}

		// preconditioner
		ierr = PetscOptionsGetString(NULL, NULL,"-gmg_mg_levels_pc_type", pname, _str_len_, &found); CHKERRQ(ierr);
		if (found){
			PetscPrintf(PETSC_COMM_WORLD, "   Multigrid smoother levels PC  : %s \n", pname);
		}
		else{
			ierr = PCGetType(pc_levels, &pc_type); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "   Multigrid smoother levels PC  : %s \n", pc_type);
		}

		ierr = PetscOptionsGetInt(NULL, NULL,"-gmg_mg_levels_ksp_max_it", &integer, &found); CHKERRQ(ierr);
		if (found){PetscPrintf(PETSC_COMM_WORLD, "   Number of smoothening steps   : %lld \n", (LLD) integer); }
		/* ----- */

		/* Coarse grid parameters */

		// Extract PETSc default parameters
		ierr = PCMGGetCoarseSolve(mg->mg.pc, &ksp_coarse); CHKERRQ(ierr);
		ierr = KSPGetPC(ksp_coarse, &pc_coarse);               CHKERRQ(ierr);

		ierr = PetscOptionsGetString(NULL, NULL,"-crs_ksp_type", pname, _str_len_, &found); CHKERRQ(ierr);
		if (found){
			PetscPrintf(PETSC_COMM_WORLD, "   Coarse level KSP              : %s \n", pname);
		}
		else{ // default
			ierr = KSPGetType(ksp_coarse, &ksp_type); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "   Coarse level KSP              : %s \n", ksp_type);
		}

		ierr = PetscOptionsGetString(NULL, NULL,"-crs_pc_type", pname, _str_len_, &found); CHKERRQ(ierr);
		if (found){
			PetscPrintf(PETSC_COMM_WORLD, "   Coarse level PC               : %s \n", pname);
		}
		else{ // default
			ierr = PCGetType(pc_coarse, &pc_type); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "   Coarse level PC               : %s \n", pc_type);
		}

		if (!strcmp(pname, PCLU)){
			// direct solver @ coarse level
			ierr = PetscOptionsGetString(NULL, NULL,"-crs_pc_factor_mat_solver_type", pname, _str_len_, &found); CHKERRQ(ierr);
			if (found){
				PetscPrintf(PETSC_COMM_WORLD, "   Coarse level solver package   : %s \n", pname);
			}
			else{ // default
				ierr = PCFactorGetMatSolverType(pc_coarse, &solver_type); CHKERRQ(ierr);
				PetscPrintf(PETSC_COMM_WORLD, "   Coarse level solver package   : %s \n", solver_type);
			}
		}
		else if (!strcmp(pname, PCREDUNDANT))
		{
			//redundant solver @ coarse level
			ierr = PetscOptionsGetInt(NULL, NULL,"-crs_pc_redundant_number", &integer, &found); CHKERRQ(ierr);
			if (found){PetscPrintf(PETSC_COMM_WORLD, "   Number of redundant solvers   : %lld \n", (LLD) integer); }
			ierr = PetscOptionsGetString(NULL, NULL,"-crs_redundant_pc_factor_mat_solver_type", pname, _str_len_, &found); CHKERRQ(ierr);
			if (found){
				PetscPrintf(PETSC_COMM_WORLD, "   Redundant solver package      : %s \n", pname);
			}

		}
		// we can add more options here if interested [e.g. for telescope]
		/* ----- */

	}
	else if (pc->type == _STOKES_USER_){

		// Direct solver
		ierr = PetscOptionsGetString(NULL, NULL,"-jp_pc_type", pname, _str_len_, &found); CHKERRQ(ierr);
		if (found){
			if(!strcmp(pname, "lu")){
				if (ISParallel(PETSC_COMM_WORLD)){
					PetscPrintf(PETSC_COMM_WORLD, "   Solver type                   : parallel direct/lu \n");
				}
				else{
					PetscPrintf(PETSC_COMM_WORLD, "   Solver type                   : serial direct/lu \n");
				}
			}
		}

		ierr = PCGetType(user->pc, &pc_type);  CHKERRQ(ierr);
		ierr = PCFactorGetMatSolverType(user->pc, &solver_type);  CHKERRQ(ierr);

		if(!solver_type)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Solver package                : petsc default\n");
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Solver package                : %s \n", solver_type);
		}

	}


	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

