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
#include "matData.h"
#include "matrix.h"
#include "matFree.h"
#include "multigrid.h"
#include "lsolve.h"
#include "JacRes.h"
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetFromOptions(PCData *pc)
{
	PetscBool mat_free;
	char      pc_type[_str_len_], bf_type[_str_len_];
	char      vs_type[_str_len_], sp_type[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set defaults
	sprintf(pc_type, "user");
	sprintf(bf_type, "upper");
	sprintf(vs_type, "user");
	sprintf(sp_type, "inv_eta");
	pc->pgamma = 1.0;

	// read options
	ierr = PetscOptionsHasName  (NULL, NULL, "-js_mat_free",   &mat_free);                CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-jp_type",       pc_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_type",       bf_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_vs_type",    vs_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_schur_type", sp_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-jp_pgamma",     &pc->pgamma,        NULL); CHKERRQ(ierr);

	if     (!strcmp(pc_type, "mg"))   pc->pc_type = _STOKES_MG_;
	else if(!strcmp(pc_type, "bf"))   pc->pc_type = _STOKES_BF_;
	else if(!strcmp(pc_type, "user")) pc->pc_type = _STOKES_USER_;
	else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes preconditioner type (jp_type): %s", pc_type);

	if     (!strcmp(bf_type, "upper")) pc->bf_type = _BF_UPPER_;
	else if(!strcmp(bf_type, "lower")) pc->bf_type = _BF_LOWER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect block factorization type (bf_type): %s", bf_type);

	if     (!strcmp(vs_type, "mg"))   pc->vs_type = _VEL_MG_;
	else if(!strcmp(vs_type, "user")) pc->vs_type = _VEL_USER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect velocity solver type (bf_vs_type): %s", vs_type);

	if     (!strcmp(sp_type, "inv_eta")) pc->sp_type = _SCHUR_INV_ETA_;
	else if(!strcmp(sp_type, "wbfbt"))   pc->sp_type = _SCHUR_WBFBT_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect Schur preconditioner type (bf_schur_type): %s", sp_type);

	// set matrix type
	if     (pc->pc_type == _STOKES_MG_)   pc->pm_type = _MONOLITHIC_;
	else if(pc->pc_type == _STOKES_BF_)   pc->pm_type = _BLOCK_;
	else if(pc->pc_type == _STOKES_USER_) pc->pm_type = _MONOLITHIC_;

	// set assembly flag
	if(mat_free) pc->ps_type =_PICARD_MAT_FREE_;
	else         pc->ps_type =_PICARD_ASSEMBLED_;

	// check errors
	if(pc->pgamma < 1.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Penalty parameter is less than unit (jp_pgamma)");

	if(pc->pc_type == _STOKES_MG_ && pc->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Coupled geometric multigrid is incompatible with matrix penalty (jp_type, jp_pgamma)");
	}

	if(pc->vs_type == _VEL_MG_ && pc->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Velocity geometric multigrid is incompatible with matrix penalty (bf_vs_type, jp_pgamma)");
	}

	if(pc->sp_type && pc->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "wBFBT preconditioner is incompatible with matrix penalty (bf_schur_type, jp_pgamma)");
	}

	// print parameters
	PetscPrintf(PETSC_COMM_WORLD, "Preconditioner parameters: \n");
	if     (pc->pm_type == _MONOLITHIC_)       PetscPrintf(PETSC_COMM_WORLD, "   Matrix type                : monolithic\n");
	else if(pc->pm_type == _BLOCK_)            PetscPrintf(PETSC_COMM_WORLD, "   Matrix type                : block\n");
	if     (pc->ps_type == _PICARD_ASSEMBLED_) PetscPrintf(PETSC_COMM_WORLD, "   Picard operator type       : assembled\n");
	else if(pc->ps_type == _PICARD_MAT_FREE_)  PetscPrintf(PETSC_COMM_WORLD, "   Picard operator type       : matrix-free\n");
	if     (pc->pc_type == _STOKES_MG_)        PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type        : coupled Galerkin geometric multigrid\n");
	else if(pc->pc_type == _STOKES_BF_)        PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type        : block factorization\n");
	else if(pc->pc_type == _STOKES_USER_)      PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type        : user-defined\n");
	if     (pc->pc_type == _STOKES_BF_)
	{
		if     (pc->bf_type == _BF_UPPER_)         PetscPrintf(PETSC_COMM_WORLD, "   Block factorization type   : upper \n");
		else if(pc->bf_type == _BF_LOWER_)         PetscPrintf(PETSC_COMM_WORLD, "   Block factorization type   : lower \n");
		if     (pc->vs_type == _VEL_MG_)           PetscPrintf(PETSC_COMM_WORLD, "   Velocity preconditioner    : Galerkin geometric multigrid\n");
		else if(pc->vs_type == _VEL_USER_)         PetscPrintf(PETSC_COMM_WORLD, "   Velocity preconditioner    : user-defined\n");
		if     (pc->sp_type == _SCHUR_INV_ETA_)    PetscPrintf(PETSC_COMM_WORLD, "   Schur preconditioner       : inverse viscosity\n");
		else if(pc->sp_type == _SCHUR_WBFBT_)      PetscPrintf(PETSC_COMM_WORLD, "   Schur preconditioner       : wBFBT\n");
	}
	if     (pc->pgamma > 1.0)                  PetscPrintf(PETSC_COMM_WORLD, "   Penalty parameter (pgamma) : %e\n", pc->pgamma);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataCreate(PCData *pc, JacRes *jr)
{
	PetscInt     buildwBFBT;
	PetscInt     buildBvv;
	MatData      *md;        // matrix assembly context

	PMatMono *Pm;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;


	ierr = PCDataSetFromOptions(pc); CHKERRQ(ierr);



	if(pc->pm_type == _MONOLITHIC_)
	{
		ierr = MatDataCreate(md, jr, _IDX_COUPLED_); CHKERRQ(ierr);

		ierr = PMatMonoCreate(Pm, md, pc->pgamma); CHKERRQ(ierr);
	}
	else if(pc->pm_type == _BLOCK_)
	{
		ierr = MatDataCreate(md, jr, _IDX_BLOCK_); CHKERRQ(ierr);

		if(pc->sp_type == _SCHUR_WBFBT_) buildwBFBT = 1;
		else                             buildwBFBT = 0;

		if     (pc->ps_type == _PICARD_MAT_FREE_) buildBvv = 0;
		else if(pc->pgamma > 1.0)                 buildBvv = 1;


		PMatBlock *Pb;
		PMatBlock *P;

		ierr = PMatBlockCreate(Pb, md, pc->pgamma, buildwBFBT, buildBvv);


	}









	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataDestroy(PCData *pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PCDataSetFromOptions(pc); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetup(PCData *pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetApply (PCData *pc, Mat P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetPicard(PCData *pc, Mat J)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(pc->ps_type == _PICARD_MAT_FREE_)
	{
		// ... matrix-free Picard operator
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)pc->md);                                   CHKERRQ(ierr);
	}
	else if(pc->ps_type == _PICARD_ASSEMBLED_)
	{
		// ... assembled Picard operator
//		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))pm->Picard); CHKERRQ(ierr);
//		ierr = MatShellSetContext(J, pm);                                      CHKERRQ(ierr);

		if(pc->pm_type == _MONOLITHIC_)
		{
			ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatMonoPicard); CHKERRQ(ierr);
			ierr = MatShellSetContext(J, (void*)pc->md);                               CHKERRQ(ierr);


		}
		else if(pc->pm_type == _BLOCK_)
		{
			ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatBlockPicard); CHKERRQ(ierr);
			ierr = MatShellSetContext(J, (void*)pc->md);                                CHKERRQ(ierr);

			PetscErrorCode PMatBlockPicard(Mat J, Vec x, Vec r);

		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

/*
PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm)
{
	//========================================================================
	// create Stokes preconditioner context
	//========================================================================

	PCStokes pc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;


	PetscErrorCode MatDataCreate(MatData *md, JacRes *jr, idxtype idxmod);


	// allocate space
	ierr = PetscMalloc(sizeof(PMatMono), (void**)&P); CHKERRQ(ierr);



	// allocate space
	ierr = PetscMalloc(sizeof(p_PCStokes), &pc); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pc, sizeof(p_PCStokes)); CHKERRQ(ierr);

	// read options
	ierr = PCStokesSetFromOptions(pc); CHKERRQ(ierr);


	// create context
//	ierr = MatDataCreate(&pm->md, jr);  CHKERRQ(ierr);


	if(pc->pctype == _STOKES_BF_)
	{
		// Block Factorization
		pc->Create  = PCStokesBFCreate;
		pc->Setup   = PCStokesBFSetup;
		pc->Destroy = PCStokesBFDestroy;
		pc->Apply   = PCStokesBFApply;
	}
	else if(pc->pctype == _STOKES_MG_)
	{
		// Galerkin multigrid
		pc->Create  = PCStokesMGCreate;
		pc->Setup   = PCStokesMGSetup;
		pc->Destroy = PCStokesMGDestroy;
		pc->Apply   = PCStokesMGApply;
	}
	else if(pc->pctype == _STOKES_USER_)
	{
		// user-defined
		pc->Create  = PCStokesUserCreate;
		pc->Setup   = PCStokesUserSetup;
		pc->Destroy = PCStokesUserDestroy;
		pc->Apply   = PCStokesUserApply;
	}

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

	// update context
//	ierr = MatDataSetup(pc->md, jr);  CHKERRQ(ierr);


	ierr = pc->Setup(pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCStokesDestroy(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = pc->Destroy(pc);     CHKERRQ(ierr);
	ierr = PetscFree(pc);       CHKERRQ(ierr);

	ierr = PetscFree(P);       CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
PetscErrorCode PCStokesBFCreate(PCStokes pc)
{
	PC          vpc;
	PCStokesBF *bf;
	PMatBlock  *P;
	MatData    *md;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;



	// allocate space
	ierr = PetscMalloc(sizeof(PMatBlock), (void**)&P); CHKERRQ(ierr);


	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesBF), (void**)&bf); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(bf, sizeof(PCStokesBF)); CHKERRQ(ierr);

	// store context
	pc->data = (void*)bf;

	// access context
	md = pc->pm->md;
	P  = (PMatBlock*)pc->pm->data;

	// create velocity solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &bf->vksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(bf->vksp,"vs_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(bf->vksp);            CHKERRQ(ierr);

	// create & set velocity multigrid preconditioner
	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGCreate(&bf->vmg, md);           CHKERRQ(ierr);
		ierr = KSPGetPC(bf->vksp, &vpc);         CHKERRQ(ierr);
		ierr = PCSetType(vpc, PCSHELL);          CHKERRQ(ierr);
		ierr = PCShellSetContext(vpc, &bf->vmg); CHKERRQ(ierr);
		ierr = PCShellSetApply(vpc, MGApply);    CHKERRQ(ierr);
	}

	// create & set pressure Schur complement solver
	if(pc->pm->buildwBFBT)
	{
		// create pressure solver
		ierr = KSPCreate(PETSC_COMM_WORLD, &bf->pksp); CHKERRQ(ierr);
		ierr = KSPSetDM(bf->pksp, P->DA_P);            CHKERRQ(ierr);
		ierr = KSPSetDMActive(bf->pksp, PETSC_FALSE);  CHKERRQ(ierr);
		ierr = KSPSetOptionsPrefix(bf->pksp,"ks_");    CHKERRQ(ierr);
		ierr = KSPSetFromOptions(bf->pksp);            CHKERRQ(ierr);
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

	if(pc->pm->buildwBFBT)
	{
		ierr = KSPDestroy(&bf->pksp);  CHKERRQ(ierr);
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

	if(pc->pm->buildwBFBT)
	{
		ierr = KSPSetOperators(bf->pksp, P->K, P->K); CHKERRQ(ierr);
		ierr = KSPSetUp(bf->pksp);                    CHKERRQ(ierr);
	}

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

	if(bf->ftype == _UPPER_)
	{
		//=======================
		// BLOCK UPPER TRIANGULAR
		//=======================

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->pm->buildwBFBT)
		{
			ierr = PCStokesBFBTApply(JP, P->rp, P->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
		else
		{
			ierr = MatMult(P->iS, P->rp, P->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}

		ierr = MatMult(P->Avp, P->xp, P->wv);    CHKERRQ(ierr); // wv = Avp*xp

		ierr = VecAXPY(P->rv, -1.0, P->wv);      CHKERRQ(ierr); // rv = rv - wv

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv
	}
	else if(bf->ftype == _LOWER_)
	{
		//=======================
		// BLOCK LOWER TRIANGULAR
		//=======================

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv

		ierr = MatMult(P->Apv, P->xv, P->wp);    CHKERRQ(ierr); // wp = Apv*xv

		ierr = VecAXPY(P->rp, -1.0, P->wp);      CHKERRQ(ierr); // rp = rp - wp

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->pm->buildwBFBT)
		{
			ierr = PCStokesBFBTApply(JP, P->rp, P->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
		else
		{
			ierr = MatMult(P->iS, P->rp, P->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
	}

	// compose approximate solution
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode wBFBTApply(wBFBTData *sp, Vec x, Vec y)
{
	//============================
	// wBFBT preconditioner action
	//============================

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// y   = -(S^⁻1)*x
	// S⁻1 =  (K^⁻1)*B*C*A*C*B^T*(K^⁻1)
	// K   =  B*C*B^T

	ierr = KSPSolve(bf->pksp, x, P->wp);   CHKERRQ(ierr); // wp = (K^⁻1)*x

	ierr = MatMult(P->Avp, P->wp, P->wv);  CHKERRQ(ierr); // wv = Avp*wp

	ierr = MatMult(P->C, P->wv, P->w);     CHKERRQ(ierr); // w = C*wv

	ierr = MatMult(P->Avv, P->w, P->wv);   CHKERRQ(ierr); // wv = Avv * w

	ierr = MatMult(P->C, P->wv, P->w);     CHKERRQ(ierr); // w = C*wv

	ierr = MatMult(P->Apv, P->w, P->wp);   CHKERRQ(ierr); // wp = Apv*w

	ierr = KSPSolve(bf->pksp, P->wp, y);   CHKERRQ(ierr); // y = (K^⁻1)*wp

	ierr = VecScale(y, -1.0);              CHKERRQ(ierr); // y = -y

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------
PetscErrorCode PCStokesMGCreate(PCStokes pc)
{
	PCStokesMG *mg;
	MatData    *md;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesMG), (void**)&mg); CHKERRQ(ierr);

	// store context
	pc->data = (void*)mg;

	// access context
	md = pc->pm->md;

	// create context
	ierr = MGCreate(&mg->mg, md); CHKERRQ(ierr);

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
	MatData      *md;
	DOFIndex     *dof;
	PetscInt      st, lnv, lnp;
	IS            isv, isp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesUser), (void**)&user); CHKERRQ(ierr);

	// store context
	pc->data = (void*)user;

	md   =  pc->pm->md;
	dof  = &md->fs->dof;
	st   =  dof->st;
	lnv  =  dof->lnv;
	lnp  =  dof->lnp;

	// create user-defined preconditioner
	// attach index sets in case fieldsplit preconditioner needs to be used
	// set additional options (whatever with -jp_ prefix)

	ierr = PCCreate(PETSC_COMM_WORLD, &user->pc); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(user->pc, "jp_");   CHKERRQ(ierr);

	// create index sets
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnv, st,     1, &isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnp, st+lnv, 1, &isp); CHKERRQ(ierr);

	// this needs to be defined before index sets can be attached
	ierr = PCSetType(user->pc, PCFIELDSPLIT); CHKERRQ(ierr);

	// attach index sets
	ierr = PCFieldSplitSetIS(user->pc, "v", isv); CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(user->pc, "p", isp); CHKERRQ(ierr);

	// destroy index sets
	ierr = ISDestroy(&isv); CHKERRQ(ierr);
	ierr = ISDestroy(&isp); CHKERRQ(ierr);

	// configure preconditioner
	ierr = PCSetFromOptions(user->pc); CHKERRQ(ierr);

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
		ierr = PCView(user->pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
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
*/
