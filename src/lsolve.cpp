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
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "matrix.h"
#include "multigrid.h"
#include "lsolve.h"
#include "JacRes.h"
//---------------------------------------------------------------------------
// * implement preconditioners in PETSc
// * add default solver options
//---------------------------------------------------------------------------
PetscErrorCode PCStokesSetFromOptions(PCStokes pc)
{
	PetscBool found;
	char      pname[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL, NULL,"-jp_type", pname, _str_len_, &found); CHKERRQ(ierr);

	if(found == PETSC_TRUE)
	{
		if(!strcmp(pname, "bf"))
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : block factorization\n");
			pc->type = _STOKES_BF_;
		}
		else if(!strcmp(pname, "mg"))
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : coupled Galerkin geometric multigrid\n");
			pc->type = _STOKES_MG_;
		}
		else if(!strcmp(pname, "user"))
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : user-defined\n");
			pc->type = _STOKES_USER_;
		}
		else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect Jacobian preconditioner type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : user-defined\n");
		pc->type = _STOKES_USER_;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm)
{
	//========================================================================
	// create Stokes preconditioner context
	//========================================================================

	PCStokes pc;
	PMatType pm_type;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PCStokes), &pc); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pc, sizeof(p_PCStokes)); CHKERRQ(ierr);

	// read options
	ierr = PCStokesSetFromOptions(pc); CHKERRQ(ierr);
	pm_type = _BLOCK_;
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
	if(pm->type != pm_type) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes preconditioner matrix type used");

	// set matrix
	pc->pm = pm;

	// create preconditioner
	ierr = pc->Create(pc); CHKERRQ(ierr);

	// return preconditioner
	(*p_pc) = pc;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesSetup(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = pc->Setup(pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesDestroy(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = pc->Destroy(pc); CHKERRQ(ierr);
	ierr = PetscFree(pc);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
PetscErrorCode PCStokesBFCreate(PCStokes pc)
{
	PC          vpc;
	PCStokesBF *bf;
	JacRes     *jr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc)
{
	PCStokesBF *bf;

	PetscBool   flg;
	char        pname[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	bf = (PCStokesBF*)pc->data;

	// set factorization type
	ierr = PetscOptionsGetString(NULL, NULL,"-bf_type", pname, _str_len_, &flg); CHKERRQ(ierr);

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
		else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect block factorization type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : upper \n");

		bf->type = _UPPER_;
	}

	// set velocity solver type
	ierr = PetscOptionsGetString(NULL, NULL,"-bf_vs_type", pname, _str_len_, &flg); CHKERRQ(ierr);

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
		else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect velocity solver type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : user-defined\n");

		bf->vtype = _VEL_USER_;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesBFDestroy(PCStokes pc)
{
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
PetscErrorCode PCStokesBFSetup(PCStokes pc)
{
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
	PetscFunctionBeginUser;

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
PetscErrorCode PCStokesMGCreate(PCStokes pc)
{
	PCStokesMG *mg;
	JacRes     *jr;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
PetscErrorCode PCStokesMGDestroy(PCStokes pc)
{
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	mg = (PCStokesMG*)pc->data;

	ierr = MGDestroy(&mg->mg); CHKERRQ(ierr);
	ierr = PetscFree(mg);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesMGSetup(PCStokes pc)
{
	PCStokesMG *mg;
	PMatMono   *P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// acces context
	mg = (PCStokesMG*)pc->data;
	P  = (PMatMono*)  pc->pm->data;

	ierr = MGSetup(&mg->mg, P->A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y)
{
	PCStokes    pc;
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	mg = (PCStokesMG*)pc->data;

	// apply multigrid preconditioner
	ierr = PCApply(mg->mg.pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//............................. USER-DEFINED ................................
//---------------------------------------------------------------------------
PetscErrorCode PCStokesUserCreate(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesUser), (void**)&user); CHKERRQ(ierr);

	// store context
	pc->data = (void*)user;

	// create user-defined preconditioner
	// attach index sets in case fieldsplit preconditioner needs to be used
	// set additional options (whatever with -jp_ prefix)
	ierr = PCCreate(PETSC_COMM_WORLD, &user->pc); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(user->pc, "jp_");   CHKERRQ(ierr);
	ierr = PCStokesUserAttachIS(pc);              CHKERRQ(ierr);
	ierr = PCSetFromOptions(user->pc);            CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesUserAttachIS(PCStokes pc)
{
	PCStokesUser *user;
	JacRes       *jr;
	DOFIndex     *dof;
	PetscInt      st, lnv, lnp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	user = (PCStokesUser*)pc->data;
	jr   =  pc->pm->jr;
	dof  = &jr->fs->dof;
	st   =  dof->st;
	lnv  =  dof->lnv;
	lnp  =  dof->lnp;

	// create index sets
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnv, st,     1, &user->isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnp, st+lnv, 1, &user->isp); CHKERRQ(ierr);

	// this needs to be defined before index sets can be attached
	ierr = PCSetType(user->pc, PCFIELDSPLIT); CHKERRQ(ierr);

	// attach index sets
	ierr = PCFieldSplitSetIS(user->pc, "v", user->isv); CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(user->pc, "p", user->isp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesUserDestroy(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get context
	user = (PCStokesUser*)pc->data;

	// cleanup
	ierr = PCDestroy(&user->pc);  CHKERRQ(ierr);
	ierr = ISDestroy(&user->isv); CHKERRQ(ierr);
	ierr = ISDestroy(&user->isp); CHKERRQ(ierr);
	ierr = PetscFree(user);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesUserSetup(PCStokes pc)
{
	PetscBool    flg;
	PCStokesUser *user;
	PMatMono     *P;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	user = (PCStokesUser*)pc->data;
	P    = (PMatMono*)    pc->pm->data;

	// compute preconditioner
	ierr = PCSetOperators(user->pc, P->A, P->A);  CHKERRQ(ierr);
	ierr = PCSetUp(user->pc);                     CHKERRQ(ierr);

	// inspect preconditioner if requested
	ierr = PetscOptionsHasName(NULL, NULL, "-pc_view", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PCView(user->pc, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y)
{
	PCStokes      pc;
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	user = (PCStokesUser*)pc->data;

	// apply user-defined preconditioner
	ierr = PCApply(user->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

