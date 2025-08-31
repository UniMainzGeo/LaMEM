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
//......................   PRECONDITIONER ROUTINES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "matData.h"
#include "matrix.h"
#include "matFree.h"
#include "multigrid.h"
#include "lsolve.h"
#include "JacRes.h"
#include "tools.h"
#include "parsing.h"
//---------------------------------------------------------------------------
PetscErrorCode PCParamSetFromOptions(PCParam *p)
{
	PetscBool mat_free;
	char      pc_type   [_str_len_], bf_type[_str_len_];
	char      vs_type   [_str_len_], sp_type[_str_len_];
	char      vs_pc_type[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set defaults
	sprintf(pc_type,    "user");
	sprintf(bf_type,    "upper");
	sprintf(vs_type,    "user");
	sprintf(vs_pc_type, "general");
	sprintf(sp_type,    "inv_eta");
	p->pgamma = 1.0;

	// read options
	ierr = PetscOptionsHasName  (NULL, NULL, "-js_mat_free",   &mat_free);                   CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-jp_type",       pc_type,    _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_type",       bf_type,    _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_vs_type",    vs_type,    _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-vs_pc_type",    vs_pc_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_schur_type", sp_type,    _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-jp_pgamma",     &p->pgamma,            NULL);  CHKERRQ(ierr);

	if     (!strcmp(pc_type, "mg"))   p->pc_type = _STOKES_MG_;
	else if(!strcmp(pc_type, "bf"))   p->pc_type = _STOKES_BF_;
	else if(!strcmp(pc_type, "user")) p->pc_type = _STOKES_USER_;
	else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect preconditioner type (jp_type): %s", pc_type);

	if     (!strcmp(bf_type, "upper")) p->bf_type = _BF_UPPER_;
	else if(!strcmp(bf_type, "lower")) p->bf_type = _BF_LOWER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect block factorization type (bf_type): %s", bf_type);

	if     (!strcmp(vs_type, "mg"))   p->vs_type = _VEL_MG_;
	else if(!strcmp(vs_type, "user")) p->vs_type = _VEL_USER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect velocity solver type (bf_vs_type): %s", vs_type);

	if(!strcmp(vs_pc_type, "lu")) p->vu_type = _VEL_USER_DIRECT_;
	else                          p->vu_type = _VEL_USER_GENERAL_;

	if     (!strcmp(sp_type, "inv_eta")) p->sp_type = _SCHUR_INV_ETA_;
	else if(!strcmp(sp_type, "wbfbt"))   p->sp_type = _SCHUR_WBFBT_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect Schur preconditioner type (bf_schur_type): %s", sp_type);

	// set assembly flag
	if(mat_free) p->ps_type =_PICARD_MAT_FREE_;
	else         p->ps_type =_PICARD_ASSEMBLED_;

	// check errors
	if(p->pgamma < 1.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Penalty parameter is less than unit (jp_pgamma)");

	if(p->pc_type == _STOKES_MG_ && p->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Coupled geometric multigrid is incompatible with matrix penalty (jp_type, jp_pgamma)");
	}

	if(p->vs_type == _VEL_MG_ && p->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Velocity geometric multigrid is incompatible with matrix penalty (bf_vs_type, jp_pgamma)");
	}

	if(p->sp_type && p->pgamma != 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "wBFBT preconditioner is incompatible with matrix penalty (bf_schur_type, jp_pgamma)");
	}

	// print parameters
	PetscPrintf(PETSC_COMM_WORLD, "Preconditioner parameters: \n");
	if     (p->ps_type == _PICARD_ASSEMBLED_) PetscPrintf(PETSC_COMM_WORLD, "   Picard operator type          : assembled\n");
	else if(p->ps_type == _PICARD_MAT_FREE_)  PetscPrintf(PETSC_COMM_WORLD, "   Picard operator type          : matrix-free\n");
	if     (p->pc_type == _STOKES_MG_)        PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : coupled Galerkin geometric multigrid\n");
	else if(p->pc_type == _STOKES_BF_)        PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : block factorization\n");
	else if(p->pc_type == _STOKES_USER_)      PetscPrintf(PETSC_COMM_WORLD, "   Preconditioner type           : user-defined\n");
	if     (p->pc_type == _STOKES_BF_)
	{
		if     (p->bf_type == _BF_UPPER_)         PetscPrintf(PETSC_COMM_WORLD, "   Block factorization type      : upper \n");
		else if(p->bf_type == _BF_LOWER_)         PetscPrintf(PETSC_COMM_WORLD, "   Block factorization type      : lower \n");
		if     (p->vs_type == _VEL_MG_)           PetscPrintf(PETSC_COMM_WORLD, "   Velocity preconditioner       : Galerkin geometric multigrid\n");
		else if(p->vs_type == _VEL_USER_)         PetscPrintf(PETSC_COMM_WORLD, "   Velocity preconditioner       : user-defined\n");
		if     (p->sp_type == _SCHUR_INV_ETA_)    PetscPrintf(PETSC_COMM_WORLD, "   Schur preconditioner          : inverse viscosity\n");
		else if(p->sp_type == _SCHUR_WBFBT_)      PetscPrintf(PETSC_COMM_WORLD, "   Schur preconditioner          : wBFBT\n");
	}
	if     (p->pgamma > 1.0)                  PetscPrintf(PETSC_COMM_WORLD, "   Penalty parameter (pgamma)    : %e\n", p->pgamma);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataCreate(PCData *pc, JacRes *jr, Mat J, Mat P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	ierr = PCParamSetFromOptions(param); CHKERRQ(ierr);

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data;

		ierr = PetscMalloc(sizeof(PCDataMG), &data); CHKERRQ(ierr);
		ierr = PetscMemzero(data, sizeof(PCDataMG)); CHKERRQ(ierr);

		ierr = PCDataMGCreate(data, param, jr, J, P); CHKERRQ(ierr);

		pc->data = (void*)data;
	}
	else if(param->pc_type == _STOKES_BF_)
	{
		PCDataBF *data;

		ierr = PetscMalloc(sizeof(PCDataBF), &data); CHKERRQ(ierr);
		ierr = PetscMemzero(data, sizeof(PCDataBF)); CHKERRQ(ierr);

		ierr = PCDataBFCreate(data, param, jr, J, P); CHKERRQ(ierr);

		pc->data = (void*)data;
	}
	else if(param->pc_type == _STOKES_USER_)
	{
		PCDataUser *data;

		ierr = PetscMalloc(sizeof(PCDataUser), &data); CHKERRQ(ierr);
		ierr = PetscMemzero(data, sizeof(PCDataUser)); CHKERRQ(ierr);

		ierr = PCDataUserCreate(data, param, jr, J, P); CHKERRQ(ierr);

		pc->data = (void*)data;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataDestroy(PCData *pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data = (PCDataMG*)pc->data;

		ierr = PCDataMGDestroy(data); CHKERRQ(ierr);
		ierr = PetscFree      (data); CHKERRQ(ierr);

	}

	else if(param->pc_type == _STOKES_BF_)
	{
		PCDataBF *data = (PCDataBF*)pc->data;

		ierr = PCDataBFDestroy(data); CHKERRQ(ierr);
		ierr = PetscFree      (data); CHKERRQ(ierr);
	}

	else if(param->pc_type == _STOKES_USER_)
	{
		PCDataUser *data = (PCDataUser*)pc->data;

		ierr = PCDataUserDestroy(data); CHKERRQ(ierr);
		ierr = PetscFree        (data); CHKERRQ(ierr);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetup(PCData *pc, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		ierr = PCDataMGSetup((PCDataMG*)pc->data, jr); CHKERRQ(ierr);
	}
	else if(param->pc_type == _STOKES_BF_)
	{
		ierr = PCDataBFSetup((PCDataBF*)pc->data, jr); CHKERRQ(ierr);
	}
	else if(param->pc_type == _STOKES_USER_)
	{
		ierr = PCDataUserSetup((PCDataUser*)pc->data, jr); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGCreate(PCDataMG *pc, PCParam *param, JacRes *jr, Mat J, Mat P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	ierr = MatDataCreate(&pc->md, jr, _IDX_COUPLED_); CHKERRQ(ierr);

	// create matrix
	PMatMono *pm = &pc->pm;

	if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = PMatMonoCreate(pm, &pc->md, param->pgamma); CHKERRQ(ierr);
	}

	// create multigrid context
	ierr = MGCreate(&pc->mg, &pc->md, pm->A); CHKERRQ(ierr);

	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->md));                                CHKERRQ(ierr);
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatMonoPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->pm));                            CHKERRQ(ierr);
	}

	// set application operator
	ierr = MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataMGApply); CHKERRQ(ierr);
	ierr = MatShellSetContext(P, (void*)pc);                                  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGDestroy(PCDataMG *pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatDataDestroy(&pc->md); CHKERRQ(ierr);

	if(pc->param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = PMatMonoDestroy(&pc->pm); CHKERRQ(ierr);
	}

	ierr = MGDestroy(&pc->mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGSetup(PCDataMG *pc, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatDataSetup(&pc->md, jr); CHKERRQ(ierr);

	if(pc->param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = PMatMonoAssemble(&pc->pm); CHKERRQ(ierr);
	}

	ierr = MGSetup(&pc->mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGApply(Mat P, Vec r, Vec x)
{
	PCDataMG *pc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(P, (void**)&pc); CHKERRQ(ierr);

	MG *mg = &pc->mg;

	ierr = PCApply(mg->pc, r, x); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFCreate(PCDataBF *pc, PCParam *param, JacRes *jr, Mat J, Mat P)
{
	PC        vpc;
	PetscInt  buildwBFBT, buildBvv, set_null_space;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	ierr = MatDataCreate(&pc->md, jr, _IDX_BLOCK_); CHKERRQ(ierr);

	// set assembly flags
	if(param->sp_type == _SCHUR_WBFBT_) buildwBFBT = 1;
	else                                buildwBFBT = 0;

	if     (param->ps_type == _PICARD_MAT_FREE_) buildBvv = 0;
	else if(param->pgamma > 1.0)                 buildBvv = 1;

	// set null space flag
	if(param->vs_type == _VEL_USER_
	&& param->vu_type != _VEL_USER_DIRECT_) set_null_space = 1;
	else                                    set_null_space = 0;

	// create matrix
	PMatBlock *pm = &pc->pm;

	ierr = PMatBlockCreate(pm, &pc->md, param->pgamma, buildwBFBT, buildBvv, set_null_space); CHKERRQ(ierr);

	// create velocity solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &pc->vksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(pc->vksp,"vs_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(pc->vksp);            CHKERRQ(ierr);

	// create & set velocity multigrid preconditioner
	if(param->vs_type == _VEL_MG_)
	{
		ierr = MGCreate(&pc->vmg, &pc->md, pm->Avv); CHKERRQ(ierr);
		ierr = KSPGetPC(pc->vksp, &vpc);             CHKERRQ(ierr);
		ierr = PCSetType(vpc, PCSHELL);              CHKERRQ(ierr);
		ierr = PCShellSetContext(vpc, &pc->vmg);     CHKERRQ(ierr);
		ierr = PCShellSetApply(vpc, MGApply);        CHKERRQ(ierr);
	}

	// create & set pressure Schur complement solver
	if(param->sp_type == _SCHUR_WBFBT_)
	{
		// create pressure solver
		ierr = KSPCreate(PETSC_COMM_WORLD, &pc->pksp); CHKERRQ(ierr);
		ierr = KSPSetDM(pc->pksp, pm->wbfbt->DA_P);    CHKERRQ(ierr);
		ierr = KSPSetDMActive(pc->pksp, PETSC_FALSE);  CHKERRQ(ierr);
		ierr = KSPSetOptionsPrefix(pc->pksp,"ks_");    CHKERRQ(ierr);
		ierr = KSPSetFromOptions(pc->pksp);            CHKERRQ(ierr);
	}
	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->md));                                CHKERRQ(ierr);
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatBlockPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->pm));                             CHKERRQ(ierr);
	}

	// set application operator
	ierr = MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataBFApply); CHKERRQ(ierr);
	ierr = MatShellSetContext(P, (void*)pc);                                  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFDestroy(PCDataBF *pc)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatDataDestroy(&pc->md); CHKERRQ(ierr);

	ierr = PMatBlockDestroy(&pc->pm); CHKERRQ(ierr);

	PetscCall(ViewSolver(pc->vksp));

	ierr = KSPDestroy(&pc->vksp);  CHKERRQ(ierr);

	if(pc->param->vs_type == _VEL_MG_)
	{
		ierr = MGDestroy(&pc->vmg); CHKERRQ(ierr);
	}

	if(pc->param->sp_type == _SCHUR_WBFBT_)
	{
		PetscCall(ViewSolver(pc->pksp));

		ierr = KSPDestroy(&pc->pksp); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFSetup(PCDataBF *pc, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PMatBlock *pm = &pc->pm;

	ierr = MatDataSetup(&pc->md, jr); CHKERRQ(ierr);

	ierr = PMatBlockAssemble(pm); CHKERRQ(ierr);

	ierr = KSPSetOperators(pc->vksp, pm->Avv, pm->Avv); CHKERRQ(ierr);

	if(pc->param->vs_type == _VEL_MG_)
	{
		ierr = MGSetup(&pc->vmg); CHKERRQ(ierr);
	}

	ierr = KSPSetUp(pc->vksp); CHKERRQ(ierr);

	if(pc->param->sp_type == _SCHUR_WBFBT_)
	{
		ierr = KSPSetOperators(pc->pksp, pm->wbfbt->K, pm->wbfbt->K); CHKERRQ(ierr);
		ierr = KSPSetUp(pc->pksp);                                    CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFApply(Mat P, Vec r, Vec x)
{
	//======================================================================
	// r - residual vector      (input)
	// x - approximate solution (output)
	//======================================================================

	PCDataBF *pc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(P, (void**)&pc); CHKERRQ(ierr);

	PMatBlock *pm = &pc->pm;

	// extract residual blocks
	ierr = VecScatterBlockToMonolithic(pm->rv, pm->rp, r, SCATTER_REVERSE); CHKERRQ(ierr);

	if(pc->param->bf_type == _BF_UPPER_)
	{
		//=======================
		// BLOCK UPPER TRIANGULAR
		//=======================

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->param->sp_type == _SCHUR_WBFBT_)
		{
			ierr = PCDataBFBTApply(pc, pm->rp, pm->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
		else if(pc->param->sp_type == _SCHUR_INV_ETA_)
		{
			ierr = MatMult(pm->iS, pm->rp, pm->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}

		ierr = MatMult(pm->Avp, pm->xp, pm->wv);    CHKERRQ(ierr); // wv = Avp*xp

		ierr = VecAXPY(pm->rv, -1.0, pm->wv);      CHKERRQ(ierr); // rv = rv - wv

		ierr = KSPSolve(pc->vksp, pm->rv, pm->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv
	}

	else if(pc->param->bf_type == _BF_LOWER_)
	{
		//=======================
		// BLOCK LOWER TRIANGULAR
		//=======================

		ierr = KSPSolve(pc->vksp, pm->rv, pm->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv

		ierr = MatMult(pm->Apv, pm->xv, pm->wp);    CHKERRQ(ierr); // wp = Apv*xv

		ierr = VecAXPY(pm->rp, -1.0, pm->wp);      CHKERRQ(ierr); // rp = rp - wp

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->param->sp_type == _SCHUR_WBFBT_)
		{
			ierr = PCDataBFBTApply(pc, pm->rp, pm->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
		else if(pc->param->sp_type == _SCHUR_INV_ETA_)
		{
			ierr = MatMult(pm->iS, pm->rp, pm->xp); CHKERRQ(ierr); // xp = (S^-1)*rp
		}
	}

	// compose approximate solution
	ierr = VecScatterBlockToMonolithic(pm->xv, pm->xp, x, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFBTApply(PCDataBF *pc, Vec x, Vec y)
{
	//============================
	// wBFBT preconditioner action
	//============================

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PMatBlock *pm = &pc->pm;
	wBFBTData *sp = pm->wbfbt;

	// y   = -(S^⁻1)*x
	// S⁻1 =  (K^⁻1)*B*C*A*C*B^T*(K^⁻1)
	// K   =  B*C*B^T

	ierr = KSPSolve(pc->pksp, x, pm->wp);    CHKERRQ(ierr); // wp = (K^⁻1)*x

	ierr = MatMult(pm->Avp, pm->wp, pm->wv); CHKERRQ(ierr); // wv = Avp*wp

	ierr = MatMult(sp->C, pm->wv, sp->w);    CHKERRQ(ierr); // w = C*wv

	ierr = MatMult(pm->Avv, sp->w, pm->wv);  CHKERRQ(ierr); // wv = Avv * w

	ierr = MatMult(sp->C, pm->wv, sp->w);    CHKERRQ(ierr); // w = C*wv

	ierr = MatMult(pm->Apv, sp->w, pm->wp);  CHKERRQ(ierr); // wp = Apv*w

	ierr = KSPSolve(pc->pksp, pm->wp, y);    CHKERRQ(ierr); // y = (K^⁻1)*wp

	ierr = VecScale(y, -1.0);                CHKERRQ(ierr); // y = -y

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//............................. USER-DEFINED ................................
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserCreate(PCDataUser *pc, PCParam *param, JacRes *jr, Mat J, Mat P)
{
	// create user-defined preconditioner
	// attach index sets to activate fieldsplit
	// set additional options with -jp_ prefix

	MatData  *md;
	DOFIndex *dof;
	IS        isv, isp;
	PetscInt  set_null_space;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	ierr = MatDataCreate(&pc->md, jr, _IDX_COUPLED_); CHKERRQ(ierr);

	// access context
	md  = &pc->md;
	dof = &md->fs->dof;

	// create matrix
	ierr = PMatMonoCreate(&pc->pm, &pc->md, param->pgamma, set_null_space = 1); CHKERRQ(ierr);

	// create preconditioner context
	ierr = PCCreate(PETSC_COMM_WORLD, &pc->pc); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(pc->pc, "jp_");   CHKERRQ(ierr);

	// create index sets
	ierr = ISCreateStride(PETSC_COMM_WORLD, dof->lnv, dof->st,            1, &isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, dof->lnp, dof->st + dof->lnv, 1, &isp); CHKERRQ(ierr);

	// activate fieldsplit preconditioner
	ierr = PCSetType(pc->pc, PCFIELDSPLIT); CHKERRQ(ierr);

	// attach index sets
	ierr = PCFieldSplitSetIS(pc->pc, "v", isv); CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc->pc, "p", isp); CHKERRQ(ierr);

	// destroy index sets
	ierr = ISDestroy(&isv); CHKERRQ(ierr);
	ierr = ISDestroy(&isp); CHKERRQ(ierr);

	// configure preconditioner
	ierr = PCSetFromOptions(pc->pc); CHKERRQ(ierr);

	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->md));                                CHKERRQ(ierr);
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatMonoPicard); CHKERRQ(ierr);
		ierr = MatShellSetContext(J, (void*)(&pc->pm));                            CHKERRQ(ierr);
	}

	// set application operator
	ierr = MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataUserApply); CHKERRQ(ierr);
	ierr = MatShellSetContext(P, (void*)pc);                                    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserDestroy(PCDataUser *pc)
{
	PetscBool flg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// view preconditioner
	ierr = PetscOptionsHasName(NULL, NULL, "-pc_view", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PCView(pc->pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	ierr = MatDataDestroy(&pc->md); CHKERRQ(ierr);

	ierr = PMatMonoDestroy(&pc->pm); CHKERRQ(ierr);

	ierr = PCDestroy(&pc->pc);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserSetup(PCDataUser *pc, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PMatMono *P = &pc->pm;

	ierr = MatDataSetup(&pc->md, jr); CHKERRQ(ierr);

	ierr = PMatMonoAssemble(P); CHKERRQ(ierr);

	ierr = PCSetOperators(pc->pc, P->A, P->A); CHKERRQ(ierr);

	ierr = PCSetUp(pc->pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserApply(Mat P, Vec r, Vec x)
{
	PCDataUser *pc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(P, (void**)&pc); CHKERRQ(ierr);

	ierr = PCApply(pc->pc, r, x); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------

