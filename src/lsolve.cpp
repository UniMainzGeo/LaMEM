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
//---------------------------------------------------------------------------
PetscErrorCode PCParamSetFromOptions(PCParam *p)
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
	p->pgamma = 1.0;

	// read options
	ierr = PetscOptionsHasName  (NULL, NULL, "-js_mat_free",   &mat_free);                CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-jp_type",       pc_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_type",       bf_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_vs_type",    vs_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL, "-bf_schur_type", sp_type, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-jp_pgamma",     &p->pgamma,        NULL); CHKERRQ(ierr);

	if     (!strcmp(pc_type, "mg"))   p->pc_type = _STOKES_MG_;
	else if(!strcmp(pc_type, "bf"))   p->pc_type = _STOKES_BF_;
	else if(!strcmp(pc_type, "user")) p->pc_type = _STOKES_USER_;
	else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes preconditioner type (jp_type): %s", pc_type);

	if     (!strcmp(bf_type, "upper")) p->bf_type = _BF_UPPER_;
	else if(!strcmp(bf_type, "lower")) p->bf_type = _BF_LOWER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect block factorization type (bf_type): %s", bf_type);

	if     (!strcmp(vs_type, "mg"))   p->vs_type = _VEL_MG_;
	else if(!strcmp(vs_type, "user")) p->vs_type = _VEL_USER_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Incorrect velocity solver type (bf_vs_type): %s", vs_type);

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
	PCParam *param;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	param = &pc->param;

	ierr = PCParamSetFromOptions(param); CHKERRQ(ierr);

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data;

		ierr = PetscMalloc(sizeof(PCDataMG), &data); CHKERRQ(ierr);
		ierr = PetscMemzero(data, sizeof(PCDataMG)); CHKERRQ(ierr);

		ierr = PCDataMGCreate(data, param, jr, J, P); CHKERRQ(ierr);

		pc->data = (void*)data;
	}
/*
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
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataDestroy(PCData *pc)
{
	PCParam *param;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		PCDataMG *data = (PCDataMG*)pc->data;

		ierr = PCDataMGDestroy(data); CHKERRQ(ierr);
		ierr = PetscFree      (data); CHKERRQ(ierr);

	}
/*
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
*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataSetup(PCData *pc, JacRes *jr)
{
	PCParam *param;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	param = &pc->param;

	if(param->pc_type == _STOKES_MG_)
	{
		ierr = PCDataMGSetup((PCDataMG*)pc->data, jr); CHKERRQ(ierr);
	}
/*
	else if(param->pc_type == _STOKES_BF_)
	{
		ierr = PCDataBFDestroy((PCDataBF*)pc->data, jr); CHKERRQ(ierr);
	}
	else if(param->pc_type == _STOKES_USER_)
	{
		ierr = PCDataUserDestroy((PCDataUser*)pc->data, jr); CHKERRQ(ierr);
	}
*/
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
	if(param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = PMatMonoCreate(&pc->pm, &pc->md, param->pgamma); CHKERRQ(ierr);
	}

	// create multigrid context
	ierr = MGCreate(&pc->mg, &pc->md, pc->pm.A); CHKERRQ(ierr);

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

	PMatMono *P = &pc->pm;

	ierr = MatDataSetup(&pc->md, jr); CHKERRQ(ierr);

	if(pc->param->ps_type == _PICARD_ASSEMBLED_)
	{
		ierr = PMatMonoAssemble(P); CHKERRQ(ierr);
	}

	ierr = MGSetup(&pc->mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PCDataMGApply(Mat P, Vec x, Vec y)
{
	PCDataMG *pc;
	MG       *mg;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatShellGetContext(P, (void**)&pc); CHKERRQ(ierr);

	mg = &pc->mg;

	ierr = PCApply(mg->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


//    USW. USF. .......



/*

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
