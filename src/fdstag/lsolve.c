//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
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
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesSetFromOptions"
PetscErrorCode PCStokesSetFromOptions(PCStokes pc)
{
	PetscBool found;
	char      pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetString(PETSC_NULL,"-jp_type", pname, MAX_NAME_LEN, &found); CHKERRQ(ierr);

	if(found == PETSC_TRUE)
	{
		if     (!strcmp(pname, "bf"))   pc->type = _STOKES_BF_;
		else if(!strcmp(pname, "mg"))   pc->type = _STOKES_MG_;
		else if(!strcmp(pname, "user")) pc->type = _STOKES_USER_;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#Incorrect Jacobian preconditioner type: %s \n", pname);
	}
	else pc->type =_STOKES_USER_;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesCreate"
PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm)
{
	//========================================================================
	// create Stokes preconditioner context
	//========================================================================

	PCStokes pc;
	PMatType pm_type;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PCStokes), &pc); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pc, sizeof(p_PCStokes)); CHKERRQ(ierr);

	// set type
	ierr = PCStokesSetFromOptions(pc); CHKERRQ(ierr);

	if(pc->type == _STOKES_BF_)
	{
		// Block Factorization
		pc->Create  = PCStokesBFCreate;
		pc->Setup   = PCStokesBFSetup;
		pc->Destroy = PCStokesBFDestroy;
		pc->Apply   = PCStokesBFApply;
		pm_type     = _BLOCK_;
	}
	else if(pc->type == _STOKES_MG_)
	{
		// Galerkin multigrid
		pc->Create  = PCStokesMGCreate;
		pc->Setup   = PCStokesMGSetup;
		pc->Destroy = PCStokesMGDestroy;
		pc->Apply   = PCStokesMGApply;
		pm_type     = _MONOLITHIC_;
	}
	else if(pc->type == _STOKES_USER_)
	{
		// user-defined
		pc->Create  = PCStokesUserCreate;
		pc->Setup   = PCStokesUserSetup;
		pc->Destroy = PCStokesUserDestroy;
		pc->Apply   = PCStokesUserApply;
		pm_type     = _MONOLITHIC_;
	}

	// check matrix type
	if(pm->type != pm_type) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes preconditioner matrix type used\n");

	// set matrix
	pc->pm = pm;

	// create preconditioner
	ierr = pc->Create(pc); CHKERRQ(ierr);

	// return preconditioner
	(*p_pc) = pc;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesSetup"
PetscErrorCode PCStokesSetup(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = pc->Setup(pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesDestroy"
PetscErrorCode PCStokesDestroy(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = pc->Destroy(pc); CHKERRQ(ierr);
	ierr = PetscFree(pc);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFCreate"
PetscErrorCode PCStokesBFCreate(PCStokes pc)
{
	PC          vpc;
	PCStokesBF *bf;
	JacRes     *jr;
	PetscBool   flg;
	char        pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesBF), (void**)&bf); CHKERRQ(ierr);

	// store context
	pc->data = (void*)bf;

	// access context
	jr = pc->pm->jr;

	// set velocity solver type
	ierr = PetscOptionsGetString(PETSC_NULL,"-bf_vs_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if     (!strcmp(pname, "mg"))   bf->vtype = _VEL_MG_;
		else if(!strcmp(pname, "user")) bf->vtype = _VEL_USER_;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#Incorrect velocity solver type: %s \n", pname);
	}
	else bf->vtype = _VEL_USER_;

	// create velocity solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &bf->vksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(bf->vksp,"vs_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(bf->vksp);            CHKERRQ(ierr);

	// create & set velocity multigrid preconditioner
	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGCreate(&bf->vmg, jr->fs, jr->ubc, IDXUNCOUPLED); CHKERRQ(ierr);

		ierr = KSPGetPC(bf->vksp, &vpc);         CHKERRQ(ierr);
		ierr = PCSetType(vpc, PCSHELL);          CHKERRQ(ierr);
		ierr = PCShellSetContext(vpc, &bf->vmg); CHKERRQ(ierr);
		ierr = PCShellSetApply(vpc, MGApply);    CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFDestroy"
PetscErrorCode PCStokesBFDestroy(PCStokes pc)
{
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;

	ierr = KSPDestroy(&bf->vksp);  CHKERRQ(ierr);

	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGDestroy(&bf->vmg); CHKERRQ(ierr);
	}

	ierr = PetscFree(bf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFSetup"
PetscErrorCode PCStokesBFSetup(PCStokes pc)
{
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;

	ierr = KSPSetOperators(bf->vksp, P->Avv, P->Avv); CHKERRQ(ierr);

	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGSetup(&bf->vmg, P->Avv); CHKERRQ(ierr);
	}

	ierr = KSPSetUp(bf->vksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFApply"
PetscErrorCode PCStokesBFApply(Mat JP, Vec r, Vec x)
{
	//======================================================================
	// r - residual vector      (input)
	// x - approximate solution (output)
	//======================================================================

	PCStokes    pc;
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;

	// extract residual blocks
	ierr = VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_REVERSE); CHKERRQ(ierr);

	ierr = VecPointwiseDivide(P->xp, P->rp, P->S); CHKERRQ(ierr); // xp = (S^-1)*rp

//	this fucking sign has tremendous influence on convergence rate! it's better this way!
//	ierr = VecScale(P->xp, -1.0);                  CHKERRQ(ierr); // xp = -xp

	ierr = MatMult(P->Avp, P->xp, P->wv);          CHKERRQ(ierr); // wv = Avp*xp

	ierr = VecAXPY(P->rv, -1.0, P->wv);            CHKERRQ(ierr); // rv = rv - wv

	ierr = KSPSolve(bf->vksp, P->rv, P->xv);       CHKERRQ(ierr); // xv = (Avv^-1)*rv

		// compose approximate solution
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGCreate"
PetscErrorCode PCStokesMGCreate(PCStokes pc)
{
	PCStokesMG *mg;
	JacRes     *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesMG), (void**)&mg); CHKERRQ(ierr);

	// store context
	pc->data = (void*)mg;

	// access context
	jr = pc->pm->jr;

	// create context
	ierr = MGCreate(&mg->mg, jr->fs, jr->cbc, IDXCOUPLED); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGDestroy"
PetscErrorCode PCStokesMGDestroy(PCStokes pc)
{
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	mg = (PCStokesMG*)pc->data;

	ierr = MGDestroy(&mg->mg); CHKERRQ(ierr);
	ierr = PetscFree(mg);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGSetup"
PetscErrorCode PCStokesMGSetup(PCStokes pc)
{
	PCStokesMG *mg;
	PMatMono   *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// acces context
	mg = (PCStokesMG*)pc->data;
	P  = (PMatMono*)  pc->pm->data;

	ierr = MGSetup(&mg->mg, P->A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGApply"
PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y)
{
	PCStokes    pc;
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	mg = (PCStokesMG*)pc->data;

	// apply multigrid preconditioner
	ierr = PCApply(mg->mg.pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//............................. USER-DEFINED ................................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserCreate"
PetscErrorCode PCStokesUserCreate(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesUser), (void**)&user); CHKERRQ(ierr);

	// store context
	pc->data = (void*)user;

	// create user-defined preconditioner
	ierr = PCCreate(PETSC_COMM_WORLD, &user->pc); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(user->pc, "jp_");   CHKERRQ(ierr);
	ierr = PCSetFromOptions(user->pc);            CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserDestroy"
PetscErrorCode PCStokesUserDestroy(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	user = (PCStokesUser*)pc->data;

	ierr = PCDestroy(&user->pc); CHKERRQ(ierr);
	ierr = PetscFree(user);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserSetup"
PetscErrorCode PCStokesUserSetup(PCStokes pc)
{
	PCStokesUser *user;
	PMatMono     *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	user = (PCStokesUser*)pc->data;
	P    = (PMatMono*)    pc->pm->data;

	// compute preconditioner
	ierr = PCSetOperators(user->pc, P->A, P->A);  CHKERRQ(ierr);
	ierr = PCSetUp(user->pc);                     CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserApply"
PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y)
{
	PCStokes      pc;
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	user = (PCStokesUser*)pc->data;

	// apply user-defined preconditioner
	ierr = PCApply(user->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
//---------------------------------------------------------------------------
// read custom stop tolerances
{
	// clear matrix norms
	bmat->nrmVV = 0.0; bmat->nrmVP = 0.0;
	bmat->nrmPV = 0.0; bmat->nrmPP = 0.0;

 	// set default tolerances
	bmat->rtolV=1e-8;
	bmat->atolV=1e-16;
	bmat->rtolP=1e-8;
	bmat->atolP=1e-16;

	// read stop tolerances
	ierr = PetscOptionsGetReal(PETSC_NULL, "-v_rtol", &bmat->rtolV, PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL, "-v_atol", &bmat->atolV, PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL, "-p_rtol", &bmat->rtolP, PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL, "-p_atol", &bmat->atolP, PETSC_NULL); CHKERRQ(ierr);

	// print norm type
	PetscPrintf(PETSC_COMM_WORLD, "StokesResidual: Using L_inf \n");

	// print stop tolerances
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Stopping conditions: \n");
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Velocity: RTOL: %e	ATOL: %e\n", bmat->rtolV, bmat->atolV); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Pressure: RTOL: %e	ATOL: %e\n", bmat->rtolP, bmat->atolP); CHKERRQ(ierr);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "KSPBlockStopTest"
PetscErrorCode KSPBlockStopTest(KSP ksp, PetscInt n, PetscScalar rnorm, KSPConvergedReason *reason, void *mctx)
{
	// monitor residual & perform stop test

	Vec         Bs, Br;
	PetscInt    maxits;
	PetscScalar nrmSolVels, nrmSolPres, nrmResVels, nrmResPres, tolVels, tolPres;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// stop warning messages for unused parameters
	if(rnorm) rnorm = 0.0;

	// access context
	BlockMat *bmat  = (BlockMat*)mctx;

	// get maximum iterations
	ierr = KSPGetTolerances(ksp, NULL, NULL, NULL, &maxits); CHKERRQ(ierr);

	// get norms of velocity & pressure solutions
	ierr = KSPBuildSolution(ksp, bmat->x, &Bs);           CHKERRQ(ierr);
	ierr = BlockMatMonolithicToBlock(bmat, Bs);           CHKERRQ(ierr);
	ierr = VecNorm(bmat->wv, NORM_INFINITY, &nrmSolVels); CHKERRQ(ierr);
	ierr = VecNorm(bmat->wp, NORM_INFINITY, &nrmSolPres); CHKERRQ(ierr);

	// get norms of velocity & pressure residulas
	ierr = KSPBuildResidual(ksp, bmat->b, bmat->x, &Br);  CHKERRQ(ierr);
	ierr = BlockMatMonolithicToBlock(bmat, Br);           CHKERRQ(ierr);
	ierr = VecNorm(bmat->wv, NORM_INFINITY, &nrmResVels); CHKERRQ(ierr);
	ierr = VecNorm(bmat->wp, NORM_INFINITY, &nrmResPres); CHKERRQ(ierr);

	// compute velocity & pressure tolerances (solution-dependent)
	tolVels = bmat->rtolV*(bmat->nrmVV*nrmSolVels + bmat->nrmVP*nrmSolPres + bmat->nrmRHSV) + bmat->atolV;
	tolPres = bmat->rtolP*(bmat->nrmPV*nrmSolVels + bmat->nrmPP*nrmSolPres + bmat->nrmRHSP) + bmat->atolP;

	// print header at first iteration
	if(n == 0 && ((PetscObject)ksp)->prefix)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norms for %s solve.\n",((PetscObject)ksp)->prefix); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  # KSP Iterations      /  Residual:           Tolerance:          /  ... for each %3D blocks  [%s] \n",
				2, ((PetscObject)ksp)->prefix); CHKERRQ(ierr);
	}

	// print residuals and tolerances
	ierr = PetscPrintf(PETSC_COMM_WORLD,"%3D KSP Residual norms  ",n); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "/  %14.12e  %14.12e  ",   nrmResVels, tolVels); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "/  %14.12e  %14.12e  \n", nrmResPres, tolPres); CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD, "/  [%s] \n", ((PetscObject)ksp)->prefix); CHKERRQ(ierr);

	// perform stop test, set convergence reason
	if(nrmSolVels && nrmResVels <= tolVels && nrmSolPres && nrmResPres <= tolPres)
	{
//		PetscPrintf(PETSC_COMM_WORLD, "Linear solution converged\n");
		*reason = KSP_CONVERGED_RTOL;
	}
	else if(n+1 == maxits)
	{
//		PetscPrintf(PETSC_COMM_WORLD, "Linear solution reached maximum number iterations\n");
		*reason = KSP_DIVERGED_ITS;
	}
	else
	{
		// continue iterations
		*reason = KSP_CONVERGED_ITERATING;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PowellHestenes"
PetscErrorCode PowellHestenes(BlockMat *bmat, Vec r, Vec x)
{
	KSP 			ksp;
	Vec 			f, u, fh, p, dp, d;
	PetscScalar 	MaximumDivergence, div_max;
	PetscInt 		PH_MaxInnerIterations, iter;
	PetscBool 		flg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read options from command line
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-ph_tol", &MaximumDivergence,     &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) MaximumDivergence = 1e-8;

	// create work vectors
	ierr = VecDuplicate(bmat->wv, &f);   CHKERRQ(ierr);
	ierr = VecDuplicate(bmat->wv, &u);   CHKERRQ(ierr);
	ierr = VecDuplicate(bmat->wv, &fh);  CHKERRQ(ierr);
	ierr = VecDuplicate(bmat->wp, &p);   CHKERRQ(ierr);
	ierr = VecDuplicate(bmat->wp, &dp);  CHKERRQ(ierr);
	ierr = VecDuplicate(bmat->wp, &d); CHKERRQ(ierr);

	// Create KSP solver environment
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);          CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, bmat->Avv, bmat->Avv); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);                     CHKERRQ(ierr);

	// copy rhs vector
	ierr = BlockMatMonolithicToBlock(bmat, r); CHKERRQ(ierr);
	ierr = VecCopy(bmat->wv, f);               CHKERRQ(ierr);

// debug ========================================
//	ierr = VecCopy(bmat->wp, d); CHKERRQ(ierr);
//	ierr = VecNorm(d, NORM_INFINITY, &div_max); CHKERRQ(ierr);
//	PetscPrintf(PETSC_COMM_WORLD, "Initial divergence norm = %1.10e \n", div_max);
//===============================================

	// set initial guess
	ierr = VecZeroEntries(p); CHKERRQ(ierr);

	// perform Powell-Hesteness iterations
	iter = 1;

	do
	{	// compute modified force vector: fh = f - Avp*p
		ierr = VecCopy(f, fh);           CHKERRQ(ierr);
		ierr = MatMult(bmat->Avp, p, u); CHKERRQ(ierr);
		ierr = VecAXPY(fh, -1.0, u);     CHKERRQ(ierr);

		// solve for velocity: u = Avv\fh
		ierr = KSPSolve(ksp, fh, u);  CHKERRQ(ierr);

		// compute divergence: d = Apv*u
		ierr = MatMult(bmat->Apv, u, d); CHKERRQ(ierr);

		// compute pressure increment: dp = kIM*div
		ierr = VecPointwiseMult(dp, bmat->kIM, d); CHKERRQ(ierr);

		// update pressure: p += dP
		ierr = VecAXPY(p, 1.0, dp); CHKERRQ(ierr);

		// compute divergence norm
		ierr = VecNorm(d, NORM_INFINITY, &div_max); CHKERRQ(ierr);

		// increment iteration count
		iter++;

		// print progress
		PetscPrintf(PETSC_COMM_WORLD, "Powell-Hesteness iteration %4D, max(Div) = %1.10e \n", iter, div_max);

		// check iteration count
		if(iter == PH_MaxInnerIterations)
		{
//			break;
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Error! Powell-Hesteness iterations did not converge; Check your setup and gamma");
		}

	} while(div_max > MaximumDivergence);

	// copy velocity & pressure
	ierr = VecCopy(u, bmat->wv);               CHKERRQ(ierr);
	ierr = VecCopy(p, bmat->wp);               CHKERRQ(ierr);
	ierr = BlockMatBlockToMonolithic(bmat, x); CHKERRQ(ierr);

	// cleaning up
	ierr = VecDestroy(&f);   CHKERRQ(ierr);
	ierr = VecDestroy(&u);   CHKERRQ(ierr);
	ierr = VecDestroy(&fh);  CHKERRQ(ierr);
	ierr = VecDestroy(&p);   CHKERRQ(ierr);
	ierr = VecDestroy(&dp);  CHKERRQ(ierr);
	ierr = VecDestroy(&d);   CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
