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

	
	PetscFunctionBeginUser;

	// set defaults
	sprintf(pc_type,    "user");
	sprintf(bf_type,    "upper");
	sprintf(vs_type,    "user");
	sprintf(vs_pc_type, "general");
	sprintf(sp_type,    "inv_eta");
	p->pgamma = 1.0;

	// read options
	PetscCall(PetscOptionsHasName  (NULL, NULL, "-js_mat_free",   &mat_free));
	PetscCall(PetscOptionsGetString(NULL, NULL, "-jp_type",       pc_type,    _str_len_, NULL));
	PetscCall(PetscOptionsGetString(NULL, NULL, "-bf_type",       bf_type,    _str_len_, NULL));
	PetscCall(PetscOptionsGetString(NULL, NULL, "-bf_vs_type",    vs_type,    _str_len_, NULL));
	PetscCall(PetscOptionsGetString(NULL, NULL, "-vs_pc_type",    vs_pc_type, _str_len_, NULL));
	PetscCall(PetscOptionsGetString(NULL, NULL, "-bf_schur_type", sp_type,    _str_len_, NULL));
	PetscCall(PetscOptionsGetScalar(NULL, NULL, "-jp_pgamma",     &p->pgamma,            NULL));

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
	
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	PetscCall(PCParamSetFromOptions(param));

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data;

		PetscCall(PetscMalloc(sizeof(PCDataMG), &data));
		PetscCall(PetscMemzero(data, sizeof(PCDataMG)));

		PetscCall(PCDataMGCreate(data, param, jr, J, P));

		pc->data = (void*)data;
	}
	else if(param->pc_type == _STOKES_BF_)
	{
		PCDataBF *data;

		PetscCall(PetscMalloc(sizeof(PCDataBF), &data));
		PetscCall(PetscMemzero(data, sizeof(PCDataBF)));

		PetscCall(PCDataBFCreate(data, param, jr, J, P));

		pc->data = (void*)data;
	}
	else if(param->pc_type == _STOKES_USER_)
	{
		PCDataUser *data;

		PetscCall(PetscMalloc(sizeof(PCDataUser), &data));
		PetscCall(PetscMemzero(data, sizeof(PCDataUser)));

		PetscCall(PCDataUserCreate(data, param, jr, J, P));

		pc->data = (void*)data;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataDestroy(PCData *pc)
{
	
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data = (PCDataMG*)pc->data;

		PetscCall(PCDataMGDestroy(data));
		PetscCall(PetscFree      (data));

	}

	else if(param->pc_type == _STOKES_BF_)
	{
		PCDataBF *data = (PCDataBF*)pc->data;

		PetscCall(PCDataBFDestroy(data));
		PetscCall(PetscFree      (data));
	}

	else if(param->pc_type == _STOKES_USER_)
	{
		PCDataUser *data = (PCDataUser*)pc->data;

		PetscCall(PCDataUserDestroy(data));
		PetscCall(PetscFree        (data));
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetup(PCData *pc, JacRes *jr)
{
	
	PetscFunctionBeginUser;

	PCParam *param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		PetscCall(PCDataMGSetup((PCDataMG*)pc->data, jr));
	}
	else if(param->pc_type == _STOKES_BF_)
	{
		PetscCall(PCDataBFSetup((PCDataBF*)pc->data, jr));
	}
	else if(param->pc_type == _STOKES_USER_)
	{
		PetscCall(PCDataUserSetup((PCDataUser*)pc->data, jr));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGCreate(PCDataMG *pc, PCParam *param, JacRes *jr, Mat J, Mat P)
{
	
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	PetscCall(MatDataCreate(&pc->md, jr, _IDX_COUPLED_));

	// create matrix
	PMatMono *pm = &pc->pm;

	if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(PMatMonoCreate(pm, &pc->md, param->pgamma));
	}

	// create multigrid context
	PetscCall(MGCreate(&pc->mg, &pc->md, pm->A));

	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->md)));
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatMonoPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->pm)));
	}

	// set application operator
	PetscCall(MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataMGApply));
	PetscCall(MatShellSetContext(P, (void*)pc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGDestroy(PCDataMG *pc)
{
	
	PetscFunctionBeginUser;

	PetscCall(MatDataDestroy(&pc->md));

	if(pc->param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(PMatMonoDestroy(&pc->pm));
	}

	PetscCall(MGDestroy(&pc->mg));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGSetup(PCDataMG *pc, JacRes *jr)
{
	
	PetscFunctionBeginUser;

	PetscCall(MatDataSetup(&pc->md, jr));

	if(pc->param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(PMatMonoAssemble(&pc->pm));
	}

	PetscCall(MGSetup(&pc->mg));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGApply(Mat P, Vec r, Vec x)
{
	PCDataMG *pc;

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(P, (void**)&pc));

	MG *mg = &pc->mg;

	PetscCall(PCApply(mg->pc, r, x));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFCreate(PCDataBF *pc, PCParam *param, JacRes *jr, Mat J, Mat P)
{
	PC        vpc;
	PetscInt  buildwBFBT, buildBvv, set_null_space;

	
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	PetscCall(MatDataCreate(&pc->md, jr, _IDX_BLOCK_));

	// set assembly flags
	if(param->sp_type == _SCHUR_WBFBT_) buildwBFBT = 1;
	else                                buildwBFBT = 0;

	if     (param->ps_type == _PICARD_MAT_FREE_) buildBvv = 0;
	else if(param->pgamma > 1.0)                 buildBvv = 1;
	else                                         buildBvv = 0;

	// set null space flag
	if(param->vs_type == _VEL_USER_
	&& param->vu_type != _VEL_USER_DIRECT_) set_null_space = 1;
	else                                    set_null_space = 0;

	// create matrix
	PMatBlock *pm = &pc->pm;

	PetscCall(PMatBlockCreate(pm, &pc->md, param->pgamma, buildwBFBT, buildBvv, set_null_space));

	// create velocity solver
	PetscCall(KSPCreate(PETSC_COMM_WORLD, &pc->vksp));
	PetscCall(KSPSetOptionsPrefix(pc->vksp,"vs_"));
	PetscCall(KSPSetFromOptions(pc->vksp));

	// create & set velocity multigrid preconditioner
	if(param->vs_type == _VEL_MG_)
	{
		PetscCall(MGCreate(&pc->vmg, &pc->md, pm->Avv));
		PetscCall(KSPGetPC(pc->vksp, &vpc));
		PetscCall(PCSetType(vpc, PCSHELL));
		PetscCall(PCShellSetContext(vpc, &pc->vmg));
		PetscCall(PCShellSetApply(vpc, MGApply));
	}

	// create & set pressure Schur complement solver
	if(param->sp_type == _SCHUR_WBFBT_)
	{
		// create pressure solver
		PetscCall(KSPCreate(PETSC_COMM_WORLD, &pc->pksp));
		PetscCall(KSPSetDM(pc->pksp, pm->wbfbt->DA_P));
#if PETSC_VERSION_LT(3, 25, 0)
		PetscCall(KSPSetDMActive(pc->pksp,                   PETSC_FALSE));
#else
		PetscCall(KSPSetDMActive(pc->pksp, KSP_DMACTIVE_ALL, PETSC_FALSE));
#endif
		PetscCall(KSPSetOptionsPrefix(pc->pksp,"ks_"));
		PetscCall(KSPSetFromOptions(pc->pksp));
	}
	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->md)));
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatBlockPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->pm)));
	}

	// set application operator
	PetscCall(MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataBFApply));
	PetscCall(MatShellSetContext(P, (void*)pc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFDestroy(PCDataBF *pc)
{
	
	PetscFunctionBeginUser;

	PetscCall(MatDataDestroy(&pc->md));

	PetscCall(PMatBlockDestroy(&pc->pm));

	PetscCall(ViewSolver(pc->vksp));

	PetscCall(KSPDestroy(&pc->vksp));

	if(pc->param->vs_type == _VEL_MG_)
	{
		PetscCall(MGDestroy(&pc->vmg));
	}

	if(pc->param->sp_type == _SCHUR_WBFBT_)
	{
		PetscCall(ViewSolver(pc->pksp));

		PetscCall(KSPDestroy(&pc->pksp));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFSetup(PCDataBF *pc, JacRes *jr)
{
	
	PetscFunctionBeginUser;

	PMatBlock *pm = &pc->pm;

	PetscCall(MatDataSetup(&pc->md, jr));

	PetscCall(PMatBlockAssemble(pm));

	PetscCall(KSPSetOperators(pc->vksp, pm->Avv, pm->Avv));

	if(pc->param->vs_type == _VEL_MG_)
	{
		PetscCall(MGSetup(&pc->vmg));
	}

	PetscCall(KSPSetUp(pc->vksp));

	if(pc->param->sp_type == _SCHUR_WBFBT_)
	{
		PetscCall(KSPSetOperators(pc->pksp, pm->wbfbt->K, pm->wbfbt->K));
		PetscCall(KSPSetUp(pc->pksp));
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

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(P, (void**)&pc));

	PMatBlock *pm = &pc->pm;

	// extract residual blocks
	PetscCall(VecScatterBlockToMonolithic(pm->rv, pm->rp, r, SCATTER_REVERSE));

	if(pc->param->bf_type == _BF_UPPER_)
	{
		//=======================
		// BLOCK UPPER TRIANGULAR
		//=======================

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->param->sp_type == _SCHUR_WBFBT_)
		{
			PetscCall(PCDataBFBTApply(pc, pm->rp, pm->xp)); // xp = (S^-1)*rp
		}
		else if(pc->param->sp_type == _SCHUR_INV_ETA_)
		{
			PetscCall(MatMult(pm->iS, pm->rp, pm->xp)); // xp = (S^-1)*rp
		}

		PetscCall(MatMult(pm->Avp, pm->xp, pm->wv)); // wv = Avp*xp

		PetscCall(VecAXPY(pm->rv, -1.0, pm->wv)); // rv = rv - wv

		PetscCall(KSPSolve(pc->vksp, pm->rv, pm->xv)); // xv = (Avv^-1)*rv
	}

	else if(pc->param->bf_type == _BF_LOWER_)
	{
		//=======================
		// BLOCK LOWER TRIANGULAR
		//=======================

		PetscCall(KSPSolve(pc->vksp, pm->rv, pm->xv)); // xv = (Avv^-1)*rv

		PetscCall(MatMult(pm->Apv, pm->xv, pm->wp)); // wp = Apv*xv

		PetscCall(VecAXPY(pm->rp, -1.0, pm->wp)); // rp = rp - wp

		// Schur complement applies negative sign internally (no negative sign here)
		if(pc->param->sp_type == _SCHUR_WBFBT_)
		{
			PetscCall(PCDataBFBTApply(pc, pm->rp, pm->xp)); // xp = (S^-1)*rp
		}
		else if(pc->param->sp_type == _SCHUR_INV_ETA_)
		{
			PetscCall(MatMult(pm->iS, pm->rp, pm->xp)); // xp = (S^-1)*rp
		}
	}

	// compose approximate solution
	PetscCall(VecScatterBlockToMonolithic(pm->xv, pm->xp, x, SCATTER_FORWARD));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataBFBTApply(PCDataBF *pc, Vec x, Vec y)
{
	//============================
	// wBFBT preconditioner action
	//============================

	
	PetscFunctionBeginUser;

	PMatBlock *pm = &pc->pm;
	wBFBTData *sp = pm->wbfbt;

	// y   = -(S^⁻1)*x
	// S⁻1 =  (K^⁻1)*B*C*A*C*B^T*(K^⁻1)
	// K   =  B*C*B^T

	PetscCall(KSPSolve(pc->pksp, x, pm->wp)); // wp = (K^⁻1)*x

	PetscCall(MatMult(pm->Avp, pm->wp, pm->wv)); // wv = Avp*wp

	PetscCall(MatMult(sp->C, pm->wv, sp->w)); // w = C*wv

	PetscCall(MatMult(pm->Avv, sp->w, pm->wv)); // wv = Avv * w

	PetscCall(MatMult(sp->C, pm->wv, sp->w)); // w = C*wv

	PetscCall(MatMult(pm->Apv, sp->w, pm->wp)); // wp = Apv*w

	PetscCall(KSPSolve(pc->pksp, pm->wp, y)); // y = (K^⁻1)*wp

	PetscCall(VecScale(y, -1.0)); // y = -y

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

	
	PetscFunctionBeginUser;

	// set parameters
	pc->param = param;

	// create matrix evaluation context
	PetscCall(MatDataCreate(&pc->md, jr, _IDX_COUPLED_));

	// access context
	md  = &pc->md;
	dof = &md->fs->dof;

	// create matrix
	PetscCall(PMatMonoCreate(&pc->pm, &pc->md, param->pgamma, set_null_space = 1));

	// create preconditioner context
	PetscCall(PCCreate(PETSC_COMM_WORLD, &pc->pc));
	PetscCall(PCSetOptionsPrefix(pc->pc, "jp_"));

	// create index sets
	PetscCall(ISCreateStride(PETSC_COMM_WORLD, dof->lnv, dof->st,            1, &isv));
	PetscCall(ISCreateStride(PETSC_COMM_WORLD, dof->lnp, dof->st + dof->lnv, 1, &isp));

	// activate fieldsplit preconditioner
	PetscCall(PCSetType(pc->pc, PCFIELDSPLIT));

	// attach index sets
	PetscCall(PCFieldSplitSetIS(pc->pc, "v", isv));
	PetscCall(PCFieldSplitSetIS(pc->pc, "p", isp));

	// destroy index sets
	PetscCall(ISDestroy(&isv));
	PetscCall(ISDestroy(&isp));

	// configure preconditioner
	PetscCall(PCSetFromOptions(pc->pc));

	// set Picard operator
	if(param->ps_type == _PICARD_MAT_FREE_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))MatFreeApplyPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->md)));
	}
	else if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		PetscCall(MatShellSetOperation(J, MATOP_MULT, (void(*)(void))PMatMonoPicard));
		PetscCall(MatShellSetContext(J, (void*)(&pc->pm)));
	}

	// set application operator
	PetscCall(MatShellSetOperation(P, MATOP_MULT, (void(*)(void))PCDataUserApply));
	PetscCall(MatShellSetContext(P, (void*)pc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserDestroy(PCDataUser *pc)
{
	PetscBool flg;

	
	PetscFunctionBeginUser;

	// view preconditioner
	PetscCall(PetscOptionsHasName(NULL, NULL, "-pc_view", &flg));

	if(flg == PETSC_TRUE)
	{
		PetscCall(PCView(pc->pc, PETSC_VIEWER_STDOUT_WORLD));
	}

	PetscCall(MatDataDestroy(&pc->md));

	PetscCall(PMatMonoDestroy(&pc->pm));

	PetscCall(PCDestroy(&pc->pc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserSetup(PCDataUser *pc, JacRes *jr)
{
	
	PetscFunctionBeginUser;

	PMatMono *P = &pc->pm;

	PetscCall(MatDataSetup(&pc->md, jr));

	PetscCall(PMatMonoAssemble(P));

	PetscCall(PCSetOperators(pc->pc, P->A, P->A));

	PetscCall(PCSetUp(pc->pc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataUserApply(Mat P, Vec r, Vec x)
{
	PCDataUser *pc;

	
	PetscFunctionBeginUser;

	PetscCall(MatShellGetContext(P, (void**)&pc));

	PetscCall(PCApply(pc->pc, r, x));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------

