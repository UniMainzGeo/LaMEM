//---------------------------------------------------------------------------
//.......................   LINEAR SOLVER ROUTINES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "lsolve.h"
#include "matrix.h"
//---------------------------------------------------------------------------


	// compute Galerkin multigrid preconditioner
//	ierr = PCSetOperators(mg->pc, P, P, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
//	ierr = PCSetUp(mg->pc); CHKERRQ(ierr);

PetscErrorCode PCStokesALCreate(PCStokesALCtx *alctx)
{


	Mat         Avv, Avp; // velocity sub-matrices
	Mat         Apv, App; // pressure sub-matrices
	Vec         M;        // diagonal penalty matrix
	PetscScalar pgamma;   // penalty parameter
	KSP         ksp;      // internal Krylov solver context


}


// PetscErrorCode PCStokesALDestroy(PCStokesALCtx *alctx);

/*
PC       pc;
KSP     *subksp;
//	PetscInt n = 2;

PetscErrorCode ierr;
PetscFunctionBegin;

FDSTAG    *fs   = nlctx->fs;
BCCtx     *sbc  = nlctx->sbc;
BlockMat  *bmat = nlctx->bmat;

ierr = PCFieldSplitGetSubKSP(bmat->pc, NULL, &subksp); CHKERRQ(ierr);
ierr = KSPGetPC(subksp[0], &pc);                       CHKERRQ(ierr);
ierr = MGCtxCreate(mg, fs, pc, IDXUNCOUPLED);          CHKERRQ(ierr);
ierr = MGCtxSetup(mg, fs, sbc, IDXUNCOUPLED);          CHKERRQ(ierr);
ierr = MGCtxSetDiagOnLevels(mg, pc);                   CHKERRQ(ierr);

PetscFree(subksp);
*/

//---------------------------------------------------------------------------


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
	else if(mode == SCATTER_REVERSE)
	{
		// monolithic-to-block
		ierr = PetscMemcpy(fp, bp,    (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(gp, bp+fs, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);
	}
	else
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown scatter mode.\n");
	}

	// restore access
	ierr = VecRestoreArray(f, &fp); CHKERRQ(ierr);
	ierr = VecRestoreArray(g, &gp); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &bp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

/*
#undef __FUNCT__
#define __FUNCT__ "BlockMatCreate"
PetscErrorCode BlockMatCreate(BlockMat *bmat, FDSTAG *fs, Vec b)
{

	// Currently just implements monolithic Picard preconditioner for the Stokes block.
	// Later this routine should create all preconditioning structures, based on user choice.

	PetscInt lnv, lnp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

// ADHOC
	// create preconditioning matrix
//	ierr = PMatCreateMonolithic(fs, &bmat->A, NULL); CHKERRQ(ierr);
//	ierr = PMatCreateMonolithic(fs, &bmat->P, NULL); CHKERRQ(ierr);

	ierr = PMatCreateBlock(fs, bmat); CHKERRQ(ierr);

	// set local number of velocity & pressure points
	lnv = fs->nXFace + fs->nYFace + fs->nZFace;
	lnp = fs->nCells;

	// create velocity & pressure work vectors
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnv, PETSC_DETERMINE, &bmat->wv);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, lnp, PETSC_DETERMINE, &bmat->wp);

	ierr = ISCreateStride(PETSC_COMM_WORLD, lnv, fs->dofcoupl.istart,       1, &bmat->isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnp, fs->dofcoupl.istart + lnv, 1, &bmat->isp); CHKERRQ(ierr);

	// create vector scatter contexts (SCATTER_FORWARD for MonoliticToBlock)
	ierr = VecScatterCreate(b, bmat->isv, bmat->wv, NULL, &bmat->vsv); CHKERRQ(ierr);
	ierr = VecScatterCreate(b, bmat->isp, bmat->wp, NULL, &bmat->vsp); CHKERRQ(ierr);

    // clear matrix norms
//	bmat->nrmVV = 0.0; bmat->nrmVP = 0.0;
//	bmat->nrmPV = 0.0; bmat->nrmPP = 0.0;

 	// set default tolerances
//	bmat->rtolV=1e-8;
//	bmat->atolV=1e-16;
//	bmat->rtolP=1e-8;
//	bmat->atolP=1e-16;

	// read stop tolerances
//	ierr = PetscOptionsGetReal(PETSC_NULL, "-v_rtol", &bmat->rtolV, PETSC_NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetReal(PETSC_NULL, "-v_atol", &bmat->atolV, PETSC_NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetReal(PETSC_NULL, "-p_rtol", &bmat->rtolP, PETSC_NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetReal(PETSC_NULL, "-p_atol", &bmat->atolP, PETSC_NULL); CHKERRQ(ierr);

	// print norm type
//	PetscPrintf(PETSC_COMM_WORLD, "StokesResidual: Using L_inf \n");

	// print stop tolerances
//	ierr = PetscPrintf(PETSC_COMM_WORLD, "Stopping conditions: \n");
//	ierr = PetscPrintf(PETSC_COMM_WORLD, "Velocity: RTOL: %e	ATOL: %e\n", bmat->rtolV, bmat->atolV); CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD, "Pressure: RTOL: %e	ATOL: %e\n", bmat->rtolP, bmat->atolP); CHKERRQ(ierr);

	// setup fieldsplit preconditioner
//	ierr = PCCreate(PETSC_COMM_WORLD, &bmat->pc);                                         CHKERRQ(ierr);
//	ierr = PCSetOptionsPrefix(bmat->pc, "fsp_");                                          CHKERRQ(ierr);
//	ierr = PCSetType(bmat->pc, PCFIELDSPLIT);                                             CHKERRQ(ierr);
//	ierr = PCFieldSplitSetType(bmat->pc, PC_COMPOSITE_SCHUR);                             CHKERRQ(ierr);
//	ierr = PCFieldSplitSetSchurFactType(bmat->pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);        CHKERRQ(ierr);
//	ierr = PCFieldSplitSetIS(bmat->pc, PETSC_NULL, bmat->isv);                            CHKERRQ(ierr);
//	ierr = PCFieldSplitSetIS(bmat->pc, PETSC_NULL, bmat->isp);                            CHKERRQ(ierr);
//	ierr = PCFieldSplitSetSchurPre(bmat->pc, PC_FIELDSPLIT_SCHUR_PRE_USER, bmat->InvEta); CHKERRQ(ierr);
//	ierr = PCSetFromOptions(bmat->pc);                                                    CHKERRQ(ierr);


    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BlockMatDestroy"
PetscErrorCode BlockMatDestroy(BlockMat *bmat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatDestroy       (&bmat->Avv);      CHKERRQ(ierr);
	ierr = MatDestroy       (&bmat->Avp);      CHKERRQ(ierr);
	ierr = MatDestroy       (&bmat->Apv);      CHKERRQ(ierr);
	ierr = VecDestroy       (&bmat->kIM);      CHKERRQ(ierr);

	ierr = ISDestroy        (&bmat->isv);      CHKERRQ(ierr);
	ierr = ISDestroy        (&bmat->isp);      CHKERRQ(ierr);
	ierr = VecScatterDestroy(&bmat->vsv);      CHKERRQ(ierr);
	ierr = VecScatterDestroy(&bmat->vsp);      CHKERRQ(ierr);

	ierr = MatDestroy       (&bmat->A);        CHKERRQ(ierr);
	ierr = MatDestroy       (&bmat->P);        CHKERRQ(ierr);
	ierr = VecDestroy       (&bmat->wv);       CHKERRQ(ierr);
	ierr = VecDestroy       (&bmat->wp);       CHKERRQ(ierr);
	ierr = PCDestroy        (&bmat->pc);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BlockMatCompute"
PetscErrorCode BlockMatCompute(BlockMat *bmat, FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx)
{
	PetscBool   flg;
	PetscScalar pgamma;

	PetscErrorCode ierr;
	PetscFunctionBegin;

// ADHOC
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-ph_gamma", &pgamma, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) pgamma = 1e2;

	// compute preconditioning matrix
	ierr = PMatAssembleBlock(fs, bc, jrctx, bmat, pgamma); CHKERRQ(ierr);

//	ierr = PCSetOperators(bmat->pc, bmat->P, bmat->P); CHKERRQ(ierr);
//	ierr = PCSetUp(bmat->pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/*
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
#define __FUNCT__ "ApplyFieldSplit"
PetscErrorCode ApplyFieldSplit(PC pc, Vec x, Vec y)
{
	BlockMat *bmat;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PCShellGetContext(pc, (void**)&bmat); CHKERRQ(ierr);

	// apply fieldsplit preconditioner
	ierr = PCApply(bmat->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
*/
//---------------------------------------------------------------------------
/*
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
	ierr = PetscOptionsGetInt(PETSC_NULL, "-ph_maxit", &PH_MaxInnerIterations, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) PH_MaxInnerIterations = 50;
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
*/
//---------------------------------------------------------------------------


/*

 #undef __FUNCT__
#define __FUNCT__ "PCCreate_Jacobi"
PETSC_EXTERN PetscErrorCode PCCreate_Jacobi(PC pc)
{
  PC_Jacobi      *jac;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /
     Creates the private data structure for this preconditioner and
     attach it to the PC object.
  /
  ierr     = PetscNewLog(pc,&jac);CHKERRQ(ierr);
  pc->data = (void*)jac;

  /
     Initialize the pointers to vectors to ZERO; these will be used to store
     diagonal entries of the matrix for fast preconditioner application.
  /
  jac->diag      = 0;
  jac->diagsqrt  = 0;
  jac->userowmax = PETSC_FALSE;
  jac->userowsum = PETSC_FALSE;
  jac->useabs    = PETSC_FALSE;

  /
      Set the pointers for the functions that are provided above.
      Now when the user-level routines (such as PCApply(), PCDestroy(), etc.)
      are called, they will automatically call these functions.  Note we
      choose not to provide a couple of these functions since they are
      not needed.
  /
  pc->ops->apply               = PCApply_Jacobi;
  pc->ops->applytranspose      = PCApply_Jacobi;
  pc->ops->setup               = PCSetUp_Jacobi;
  pc->ops->reset               = PCReset_Jacobi;
  pc->ops->destroy             = PCDestroy_Jacobi;
  pc->ops->setfromoptions      = PCSetFromOptions_Jacobi;
  pc->ops->view                = 0;
  pc->ops->applyrichardson     = 0;
  pc->ops->applysymmetricleft  = PCApplySymmetricLeftOrRight_Jacobi;
  pc->ops->applysymmetricright = PCApplySymmetricLeftOrRight_Jacobi;

  ierr = PetscObjectComposeFunction((PetscObject)pc,"PCJacobiSetUseRowMax_C",PCJacobiSetUseRowMax_Jacobi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pc,"PCJacobiSetUseRowSum_C",PCJacobiSetUseRowSum_Jacobi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pc,"PCJacobiSetUseAbs_C",PCJacobiSetUseAbs_Jacobi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


 */
