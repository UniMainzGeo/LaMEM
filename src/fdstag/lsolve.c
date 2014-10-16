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
// * get rid of AL (as this is a special case of block upper triangular with penalty)
// * get rid of fieldsplit + matnest (just implement block upper triangular)
// * add penalty terms to matrix types (block/monolithic)
// * matrix type (block/monolithic) and penalty should be a parameter
// * implement Picard Jacobian as a function defined in matrix structure
// * properly implement Picard Jacobian if penalty is added (as this is screwed up now)
// * ideally bf & mg should be PETSc preconditioners
// * preallocation for Stokes restriction/interpolation
// * restriction/interpolation for FDSTG Laplacian (temperature, Schur complement approximation)
// * adding ones to diagonal (PCSetUp) ???
// * near null space for Stokes
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesSetFromOptions"
PetscErrorCode PCStokesSetFromOptions(PCStokes pc)
{
	PetscBool found;
	char      pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetString(PETSC_NULL,"-pstokes", pname, MAX_NAME_LEN, &found); CHKERRQ(ierr);

	if(found == PETSC_TRUE)
	{
		if     (!strcmp(pname, "al")) pc->ptype = STOKES_AL;
		else if(!strcmp(pname, "bf")) pc->ptype = STOKES_BF;
		else if(!strcmp(pname, "mg")) pc->ptype = STOKES_MG;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#Incorrect Stokes preconditioner: %s \n", pname);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesCreate"
PetscErrorCode PCStokesCreate(PCStokes *p_pc, JacRes *jr)
{
	//========================================================================
	// create Stokes preconditioner context
	//========================================================================

	PCStokes pc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PCStokes), &pc); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pc, sizeof(p_PCStokes)); CHKERRQ(ierr);

	// set type
	ierr = PCStokesSetFromOptions(pc); CHKERRQ(ierr);

	// assign data
	pc->jr = jr;

	if(pc->ptype == STOKES_AL)
	{
		// augmented Lagrangian
		pc->Create  = &PCStokesALCreate;
		pc->Setup   = &PCStokesALSetup;
		pc->Destroy = &PCStokesALDestroy;
		pc->Apply   = &PCStokesALApply;
		pc->Picard  = &PCStokesALPicard;
	}
	else if(pc->ptype == STOKES_BF)
	{
		// Block Factorization
		pc->Create  = &PCStokesBFCreate;
		pc->Setup   = &PCStokesBFSetup;
		pc->Destroy = &PCStokesBFDestroy;
		pc->Apply   = &PCStokesBFApply;
		pc->Picard  = &PCStokesBFPicard;
	}
	else if(pc->ptype == STOKES_MG)
	{
		// Galerkin multigrid
		pc->Create  = &PCStokesMGCreate;
		pc->Setup   = &PCStokesMGSetup;
		pc->Destroy = &PCStokesMGDestroy;
		pc->Apply   = &PCStokesMGApply;
		pc->Picard  = &PCStokesMGPicard;
	}

	// create actual preconditioner
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
//....................... AUGMENTED LAGRANGIAN ..............................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesALCreate"
PetscErrorCode PCStokesALCreate(PCStokes pc)
{
	PCStokesAL  *al;
	JacRes      *jr;
	DOFIndex    *udof;
	PetscBool    flg;
	PetscInt     lnv, lnp;
	PetscScalar  pgamma;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesAL), (void**)&al); CHKERRQ(ierr);

	// store context
	pc->data = (void*)al;

	jr   = pc->jr;
	udof = &(jr->fs->udof); // uncoupled

	// get number of local velocity & pressure rows
	lnv = udof->numdof;
	lnp = udof->numdofp;

	// create matrices
	ierr = PMatCreateBlock(jr->fs, &al->P, &al->M); CHKERRQ(ierr);

	// read & store penalty parameter
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-pgamma", &pgamma, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) pgamma = 1e-3;
	al->pgamma = pgamma;

	// create augmented Lagrangian preconditioner solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &al->ksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(al->ksp,"al_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(al->ksp);            CHKERRQ(ierr);

	// create work vectors
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &al->f); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnp, PETSC_DETERMINE, &al->h); CHKERRQ(ierr);
	ierr = VecDuplicate(al->f, &al->u); CHKERRQ(ierr);
	ierr = VecDuplicate(al->h, &al->p); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesALDestroy"
PetscErrorCode PCStokesALDestroy(PCStokes pc)
{
	PCStokesAL *al;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	al = (PCStokesAL*)pc->data;

	ierr = BMatDestroy(&al->P);   CHKERRQ(ierr);
	ierr = MatDestroy (&al->M);   CHKERRQ(ierr);
	ierr = KSPDestroy (&al->ksp); CHKERRQ(ierr);
	ierr = VecDestroy (&al->f);   CHKERRQ(ierr);
	ierr = VecDestroy (&al->h);   CHKERRQ(ierr);
	ierr = VecDestroy (&al->u);   CHKERRQ(ierr);
	ierr = VecDestroy (&al->p);   CHKERRQ(ierr);

	ierr = PetscFree(al); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesALSetup"
PetscErrorCode PCStokesALSetup(PCStokes pc)
{
	PCStokesAL *al;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	al = (PCStokesAL*)pc->data;

	// assemble matrices
	ierr = PMatAssembleBlock(pc->jr, &al->P, al->M, al->pgamma); CHKERRQ(ierr);

	// tell to recompute preconditioner
	ierr = KSPSetOperators(al->ksp, al->P.Avv, al->P.Avv); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesALApply"
PetscErrorCode PCStokesALApply(Mat P, Vec r, Vec z)
{
	//======================================================================
	// r - residual vector      (input)
	// z - approximate solution (output)
	//======================================================================

	PCStokesAL *al;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access data
	ierr = MatShellGetContext(P, (void**)&al); CHKERRQ(ierr);

	// extract residual blocks
	ierr = VecScatterBlockToMonolithic(al->f, al->h, r, SCATTER_REVERSE); CHKERRQ(ierr);

	// compose right-hand-side vector: f = f - Avp*M*h
	ierr = MatMult(al->M, al->h, al->p);     CHKERRQ(ierr);
	ierr = MatMult(al->P.Avp, al->p, al->u); CHKERRQ(ierr);
	ierr = VecAXPY(al->f, -1.0, al->u);      CHKERRQ(ierr);

	// solve for velocity: u = Avv\f
	ierr = KSPSolve(al->ksp, al->f, al->u); CHKERRQ(ierr);

	// solve for pressure: p = M*(h - Apv*u)
	ierr = MatMult(al->P.Apv, al->u, al->p); CHKERRQ(ierr);
	ierr = VecAXPY(al->h, -1.0, al->p);      CHKERRQ(ierr);
	ierr = MatMult(al->M, al->h, al->p);     CHKERRQ(ierr);

	// compose approximate solution
	ierr = VecScatterBlockToMonolithic(al->u, al->p, z, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesALPicard"
PetscErrorCode PCStokesALPicard(Mat J, Vec x, Vec y)
{
	PCStokesAL *al;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(J, (void**)&al); CHKERRQ(ierr);

	// compute action of Picard Jacobian
	ierr = MatMult(al->P.A, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFCreate"
PetscErrorCode PCStokesBFCreate(PCStokes pc)
{
	PCStokesBF *bf;
	JacRes     *jr;
	KSP        *subksp;
	DOFIndex   *cdof, *udof;
	PetscInt    sid, lnv, lnp;
	PetscBool   flg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesBF), (void**)&bf); CHKERRQ(ierr);

	// store context
	pc->data = (void*)bf;

	jr   = pc->jr;
	cdof = &(jr->fs->cdof);
	udof = &(jr->fs->udof);

	// get number of local rows & starting index in global numbering
	sid = cdof->istart;
	lnv = udof->numdof;
	lnp = udof->numdofp;

	// create matrices
	ierr = PMatCreateBlock(jr->fs, &bf->P, &bf->M); CHKERRQ(ierr);

	// create index sets
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnv, sid,       1, &bf->isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnp, sid + lnv, 1, &bf->isp); CHKERRQ(ierr);

	// setup block factorization preconditioner
	ierr = PCCreate(PETSC_COMM_WORLD, &bf->pc);                                  CHKERRQ(ierr);
	ierr = PCSetOperators(bf->pc, bf->P.A, bf->P.A);                             CHKERRQ(ierr);
	ierr = PCSetType(bf->pc, PCFIELDSPLIT);                                      CHKERRQ(ierr);
	ierr = PCFieldSplitSetType(bf->pc, PC_COMPOSITE_SCHUR);                      CHKERRQ(ierr);
	ierr = PCFieldSplitSetSchurFactType(bf->pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER); CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(bf->pc, "u", bf->isv);                              CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(bf->pc, "p", bf->isp);                              CHKERRQ(ierr);
	ierr = PCFieldSplitSetSchurPre(bf->pc, PC_FIELDSPLIT_SCHUR_PRE_USER, bf->M); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(bf->pc, "bf_");                                    CHKERRQ(ierr);
	ierr = PCSetFromOptions(bf->pc);                                             CHKERRQ(ierr);

	// check whether Galerkin multigrid is requested for the velocity block
	ierr = PetscOptionsHasName(NULL, "-velgmg", &flg);

	// setup velocity multigrid
	if(flg == PETSC_TRUE)
	{
		// set type
		bf->vtype = VEL_MG;

		// retrieve velocity preconditioner handle
		ierr = PCSetUp(bf->pc);                              CHKERRQ(ierr);
		ierr = PCFieldSplitGetSubKSP(bf->pc, NULL, &subksp); CHKERRQ(ierr);
		ierr = KSPGetPC(subksp[0], &bf->vpc);                CHKERRQ(ierr);

		// create velocity multigrid preconditioner
		ierr = MGCtxCreate(&bf->vctx, jr->fs, jr->ubc, bf->vpc, IDXUNCOUPLED); CHKERRQ(ierr);
		ierr = PetscFree(subksp);                                              CHKERRQ(ierr);

		// MatSetNearNullSpace(J, nearNullSpace);
		// MatNullSpaceCreate
	}
	else
	{
		bf->vtype = VEL_USER;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFDestroy"
PetscErrorCode PCStokesBFDestroy(PCStokes pc)
{
	PetscBool   flg;
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	bf = (PCStokesBF*)pc->data;

	// view block factorization preconditioner if required
	ierr = PetscOptionsHasName(NULL, "-bf_pc_view", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PCView(bf->pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	ierr = BMatDestroy(&bf->P);   CHKERRQ(ierr);
	ierr = MatDestroy (&bf->M);   CHKERRQ(ierr);
	ierr = ISDestroy  (&bf->isv); CHKERRQ(ierr);
	ierr = ISDestroy  (&bf->isp); CHKERRQ(ierr);
	ierr = PCDestroy  (&bf->pc);  CHKERRQ(ierr);

	if(bf->vtype == VEL_MG)
	{
		ierr = MGCtxDestroy(&bf->vctx); CHKERRQ(ierr);
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

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	bf = (PCStokesBF*)pc->data;

	// assemble matrices
	ierr = PMatAssembleBlock(pc->jr, &bf->P, bf->M, 0); CHKERRQ(ierr);

	// tell to recompute preconditioner
	ierr = PCSetOperators(bf->pc, bf->P.A, bf->P.A); CHKERRQ(ierr);

	if(bf->vtype == VEL_MG)
	{
		ierr = MGCtxSetup(&bf->vctx, IDXUNCOUPLED);      CHKERRQ(ierr);
		ierr = MGCtxSetDiagOnLevels(&bf->vctx, bf->vpc); CHKERRQ(ierr);

		// dump matrices to MATLAB if requested
		ierr = MGCtxDumpMat(&bf->vctx, bf->vpc); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFApply"
PetscErrorCode PCStokesBFApply(Mat P, Vec x, Vec y)
{
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(P, (void**)&bf); CHKERRQ(ierr);

	// apply block factorization preconditioner
	ierr = PCApply(bf->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFPicard"
PetscErrorCode PCStokesBFPicard(Mat J, Vec x, Vec y)
{
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(J, (void**)&bf); CHKERRQ(ierr);

	// compute action of Picard Jacobian
	ierr = MatMult(bf->P.A, x, y); CHKERRQ(ierr);

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
	jr       = pc->jr;

	// create matrices
	ierr = PMatCreateMonolithic(jr->fs, &mg->P, &mg->M); CHKERRQ(ierr);

	// create pc object
	ierr = PCCreate(PETSC_COMM_WORLD, &mg->pc); CHKERRQ(ierr);

	// create velocity multigrid preconditioner
	ierr = MGCtxCreate(&mg->ctx, jr->fs, jr->cbc, mg->pc, IDXCOUPLED); CHKERRQ(ierr);

	// create work vector
	ierr = VecDuplicate(jr->gsol, &mg->w); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGDestroy"
PetscErrorCode PCStokesMGDestroy(PCStokes pc)
{
	PetscBool   flg;
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	mg = (PCStokesMG*)pc->data;

	// view preconditioner if required
	ierr = PetscOptionsHasName(NULL, "-gmg_pc_view", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		ierr = PCView(mg->pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	ierr = MatDestroy  (&mg->P);   CHKERRQ(ierr);
	ierr = MatDestroy  (&mg->M);   CHKERRQ(ierr);
	ierr = MGCtxDestroy(&mg->ctx); CHKERRQ(ierr);
	ierr = PCDestroy   (&mg->pc);  CHKERRQ(ierr);
	ierr = VecDestroy  (&mg->w);   CHKERRQ(ierr);

	ierr = PetscFree(mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGSetup"
PetscErrorCode PCStokesMGSetup(PCStokes pc)
{
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	mg = (PCStokesMG*)pc->data;

	// assemble monolithic matrix
	ierr = PMatAssembleMonolithic(pc->jr, mg->P, mg->M); CHKERRQ(ierr);

	// tell to recompute preconditioner
	ierr = PCSetOperators(mg->pc, mg->P, mg->P);  CHKERRQ(ierr);

	// setup multilevel structure
	ierr = MGCtxSetup(&mg->ctx, IDXCOUPLED);       CHKERRQ(ierr);
	ierr = MGCtxSetDiagOnLevels(&mg->ctx, mg->pc); CHKERRQ(ierr);

	// dump matrices to MATLAB if requested
	ierr = MGCtxDumpMat(&mg->ctx, mg->pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGApply"
PetscErrorCode PCStokesMGApply(Mat P, Vec x, Vec y)
{
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(P, (void**)&mg); CHKERRQ(ierr);

	// apply block factorization preconditioner
	ierr = PCApply(mg->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGPicard"
PetscErrorCode PCStokesMGPicard(Mat J, Vec x, Vec y)
{
	// actual operation is: y = J*x = P*x - M*x

	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(J, (void**)&mg); CHKERRQ(ierr);

	// compute action of preconditioner matrix
	ierr = MatMult(mg->P, x, y);     CHKERRQ(ierr);
	// compute compensation
	ierr = MatMult(mg->M, x, mg->w); CHKERRQ(ierr);
	// update result
	ierr = VecAXPY(y, -1.0, mg->w);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecScatterBlockToMonolithic"
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode)
{
	// scatter block vectors to monolithic format forward & reverse

	PetscInt     fs,  gs,  bs;
	PetscScalar *fp, *gp, *bp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get sizes of the blocks
	ierr = VecGetLocalSize(f, &fs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(g, &gs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(b, &bs); CHKERRQ(ierr);

	if(bs != fs+gs)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Block sizes don't match monolithic format\n");
	}

	// access vectors
	ierr = VecGetArray(f, &fp); CHKERRQ(ierr);
	ierr = VecGetArray(g, &gp); CHKERRQ(ierr);
	ierr = VecGetArray(b, &bp); CHKERRQ(ierr);

	if(mode == SCATTER_FORWARD)
	{
		// block-to-monolithic
		ierr = PetscMemcpy(bp,    fp, (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(bp+fs, gp, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);
	}
	if(mode == SCATTER_REVERSE)
	{
		// monolithic-to-block
		ierr = PetscMemcpy(fp, bp,    (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(gp, bp+fs, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);
	}

	// restore access
	ierr = VecRestoreArray(f, &fp); CHKERRQ(ierr);
	ierr = VecRestoreArray(g, &gp); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &bp); CHKERRQ(ierr);

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
