//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
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
//---------------------------------------------------------------------------
// * implement preconditioners in PETSc
// * add default solver options
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
		if(!strcmp(pname, "bf"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : block factorization\n");
            
			pc->type = _STOKES_BF_;
		}
		else if(!strcmp(pname, "mg"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : coupled Galerkin geometric multigrid\n");
			pc->type = _STOKES_MG_;
		}
		else if(!strcmp(pname, "user"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : user-defined\n");
			pc->type = _STOKES_USER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect Jacobian preconditioner type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : user-defined\n");
		pc->type = _STOKES_USER_;
	}

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

	// read options
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
	if(pm->type != pm_type) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect Stokes preconditioner matrix type used");

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

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesBF), (void**)&bf); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(bf, sizeof(PCStokesBF)); CHKERRQ(ierr);

	// store context
	pc->data = (void*)bf;

	// read options
	ierr = PCStokesBFSetFromOptions(pc); CHKERRQ(ierr);

	// access context
	jr = pc->pm->jr;

	// create velocity solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &bf->vksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(bf->vksp,"vs_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(bf->vksp);            CHKERRQ(ierr);

	// create & set velocity multigrid preconditioner
	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGCreate(&bf->vmg, jr);           CHKERRQ(ierr);
		ierr = KSPGetPC(bf->vksp, &vpc);         CHKERRQ(ierr);
		ierr = PCSetType(vpc, PCSHELL);          CHKERRQ(ierr);
		ierr = PCShellSetContext(vpc, &bf->vmg); CHKERRQ(ierr);
		ierr = PCShellSetApply(vpc, MGApply);    CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFSetFromOptions"
PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc)
{
	PCStokesBF *bf;

	PetscBool   flg;
	char        pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;

	// set factorization type
	ierr = PetscOptionsGetString(PETSC_NULL,"-bf_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(!strcmp(pname, "upper"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : upper \n");

			bf->type = _UPPER_;
		}
		else if(!strcmp(pname, "lower"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : lower \n");

			bf->type = _LOWER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect block factorization type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : upper \n");

		bf->type = _UPPER_;
	}

	// set velocity solver type
	ierr = PetscOptionsGetString(PETSC_NULL,"-bf_vs_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(!strcmp(pname, "mg"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : Galerkin geometric multigrid\n");

			bf->vtype = _VEL_MG_;
		}
		else if(!strcmp(pname, "user"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : user-defined\n");

			bf->vtype = _VEL_USER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect velocity solver type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : user-defined\n");

		bf->vtype = _VEL_USER_;
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

	if(bf->type == _UPPER_)
	{
		//=======================
		// BLOCK UPPER TRIANGULAR
		//=======================

		// Schur complement already contains negative sign (no negative sign here)
		ierr = MatMult(P->iS, P->rp, P->xp);     CHKERRQ(ierr); // xp = (S^-1)*rp

		ierr = MatMult(P->Avp, P->xp, P->wv);    CHKERRQ(ierr); // wv = Avp*xp

		ierr = VecAXPY(P->rv, -1.0, P->wv);      CHKERRQ(ierr); // rv = rv - wv

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv
	}
	else if(bf->type == _LOWER_)
	{
		//=======================
		// BLOCK LOWER TRIANGULAR
		//=======================

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv

		ierr = MatMult(P->Apv, P->xv, P->wp);    CHKERRQ(ierr); // wp = Apv*xv

		ierr = VecAXPY(P->rp, -1.0, P->wp);      CHKERRQ(ierr); // rp = rp - wp

		// Schur complement already contains negative sign (no negative sign here)
		ierr = MatMult(P->iS, P->rp, P->xp);     CHKERRQ(ierr); // xp = (S^-1)*rp
	}

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
	ierr = MGCreate(&mg->mg, jr); CHKERRQ(ierr);

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
*/
