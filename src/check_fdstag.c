/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   check_fdstag.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//.........   LaMEM - FDSTAG CANONICAL INTERFACE CHECKING ROUTINES   ........
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "matFree.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "tools.h"
#include "check_fdstag.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacTest"
PetscErrorCode JacTest(JacRes *jr, Vec diff)
{
	Mat          MF, MFFD;
	DOFIndex    *dof;
	BCCtx       *bc;
	Vec          z, mf, mffd;
	PetscScalar  nrm, *sol;
	PetscBool    set;
	PetscInt     i, num, *list, pmf, pmffd;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	dof = &(jr->fs->dof);
	bc  =  jr->bc;

	// read and set debugging flags
	pmf   = 0;
	pmffd = 0;

	ierr = PetscOptionsHasName(NULL, NULL, "-mf",   &set); CHKERRQ(ierr); if(set == PETSC_TRUE) pmf      = 1;
	ierr = PetscOptionsHasName(NULL, NULL, "-mffd", &set); CHKERRQ(ierr); if(set == PETSC_TRUE) pmffd    = 1;

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &MF); CHKERRQ(ierr);
	ierr = MatSetUp(MF);                              CHKERRQ(ierr);

	// create finite-difference Jacobian
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, &MFFD); CHKERRQ(ierr);
	ierr = MatSetUp(MFFD);                        CHKERRQ(ierr);

	// setup vectors
	ierr = VecDuplicate(jr->gres, &z);    CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &mf);   CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &mffd); CHKERRQ(ierr);

	//===================================
	// initialize vector to be multiplied
	//===================================

	ierr = VecCopy    (jr->gsol, z); CHKERRQ(ierr);
	ierr = VecGetArray(z, &sol);     CHKERRQ(ierr);

	// velocity
	num   = bc->vNumSPC;
	list  = bc->vSPCList;

	for(i = 0; i < num; i++) sol[list[i]] = 0.0;

	// pressure
	num   = bc->pNumSPC;
	list  = bc->pSPCList;

	for(i = 0; i < num; i++) sol[list[i]] = 0.0;

	ierr = VecRestoreArray(z, &sol); CHKERRQ(ierr);
	ierr = VecNormalize   (z, NULL); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"=======================================================\n");

	ierr = VecNorm(z, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the vector to be multiplied: %g\n", nrm);

	//===================================
	// MF Jacobian
	//===================================

	ierr = VecZeroEntries(mf);  CHKERRQ(ierr);

	// analytical Jacobian
	ierr = MatShellSetOperation(MF, MATOP_MULT, (void(*)(void))JacApplyJacobian); CHKERRQ(ierr);
	ierr = MatShellSetContext(MF, (void*)jr);                                     CHKERRQ(ierr);

	ierr = MatAssemblyBegin(MF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (MF, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatMult(MF, z, mf); CHKERRQ(ierr);

	ierr = VecNorm(mf, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the analytical Jacobian-vector product: %g\n", nrm);

	//===================================
	// MFFD Jacobian
	//===================================

	ierr = VecZeroEntries(mffd); CHKERRQ(ierr);

	// ... matrix-free finite-difference (MMFD)
	ierr = MatMFFDSetFunction(MFFD, FormResidualMFFD, (void*)jr); CHKERRQ(ierr);
	ierr = MatMFFDSetBase(MFFD, jr->gsol, NULL);                  CHKERRQ(ierr);

	ierr = MatAssemblyBegin(MFFD, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (MFFD, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatMult(MFFD, z, mffd); CHKERRQ(ierr);

	ierr = VecNorm(mffd, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the finite difference Jacobian-vector product: %g\n", nrm);

	//===================================
	// difference
	//===================================

	ierr =  VecWAXPY(diff, -1.0, mf, mffd); CHKERRQ(ierr);

	ierr = VecNorm(diff, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the difference vector: %g\n", nrm);

	PetscPrintf(PETSC_COMM_WORLD,"=======================================================\n");

	// ACHTUNG!!!
	if(pmf)   { ierr = VecCopy(mf,   diff); CHKERRQ(ierr); }
	if(pmffd) { ierr = VecCopy(mffd, diff); CHKERRQ(ierr); }

	ierr = MatDestroy(&MF);   CHKERRQ(ierr);
	ierr = MatDestroy(&MFFD); CHKERRQ(ierr);
	ierr = VecDestroy(&z);    CHKERRQ(ierr);
	ierr = VecDestroy(&mf);   CHKERRQ(ierr);
	ierr = VecDestroy(&mffd); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PicardTest"
PetscErrorCode PicardTest(NLSol *nl, Vec diff)
{
	Mat          J;
	DOFIndex    *dof;
	BCCtx       *bc;
	PMat         pm;
	JacRes      *jr;
	Vec          z, wm, wmf;
	PetscScalar  nrm, *sol;
	PetscInt     i, num, *list;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	pm  = nl->pc->pm;
	jr  = pm->jr;
	dof = &(jr->fs->dof);
	bc  =  jr->bc;

	// assemble preconditioner with restrictions
	ierr = PMatAssemble(pm); CHKERRQ(ierr);

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->ln, dof->ln,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &J); CHKERRQ(ierr);
	ierr = MatSetUp(J);                              CHKERRQ(ierr);

	// setup vectors
	ierr = VecDuplicate(jr->gres, &z);    CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &wm);   CHKERRQ(ierr);
	ierr = VecDuplicate(jr->gres, &wmf);  CHKERRQ(ierr);

	//===================================
	// initialize vector to be multiplied
	//===================================

	ierr = VecCopy    (jr->gsol, z); CHKERRQ(ierr);
	ierr = VecGetArray(z, &sol);     CHKERRQ(ierr);

	// velocity
	num   = bc->vNumSPC;
	list  = bc->vSPCList;

	for(i = 0; i < num; i++) sol[list[i]] = 0.0;

	// pressure
	num   = bc->pNumSPC;
	list  = bc->pSPCList;

	for(i = 0; i < num; i++) sol[list[i]] = 0.0;

	ierr = VecRestoreArray(z, &sol); CHKERRQ(ierr);
	ierr = VecNormalize   (z, NULL); CHKERRQ(ierr);

	//===================================

	ierr = VecZeroEntries(wm);  CHKERRQ(ierr);
	ierr = VecZeroEntries(wmf); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"=======================================================\n");

	ierr = VecNorm(z, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the vector to be multiplied: %g\n", nrm);

	// assembled Jacobian
	ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))pm->Picard); CHKERRQ(ierr);
	ierr = MatShellSetContext(J, pm->data);                                CHKERRQ(ierr);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatMult(J, z, wm); CHKERRQ(ierr);

	ierr = VecNorm(wm, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the assembled residual vector: %g\n", nrm);

	// matrix-free Jacobian
	ierr = MatShellSetOperation(J, MATOP_MULT, (void(*)(void))JacApplyPicard); CHKERRQ(ierr);
	ierr = MatShellSetContext(J, (void*)jr);                                   CHKERRQ(ierr);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatMult(J, z, wmf); CHKERRQ(ierr);

	ierr = VecNorm(wmf, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the matrix-free residual vector: %g\n", nrm);

	// difference
	ierr =  VecWAXPY(diff, -1.0, wm, wmf); CHKERRQ(ierr);

	ierr = VecNorm(diff, NORM_2, &nrm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Norm of the difference vector: %g\n", nrm);

	PetscPrintf(PETSC_COMM_WORLD,"=======================================================\n");

	ierr = MatDestroy(&J);    CHKERRQ(ierr);
	ierr = VecDestroy(&z);    CHKERRQ(ierr);
	ierr = VecDestroy(&wm);   CHKERRQ(ierr);
	ierr = VecDestroy(&wmf);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "DarcyPostProcess"
PetscErrorCode DarcyPostProcess(NLCtx *nlctx, UserCtx *user)
{
	FILE        *db;
	PetscBool   flg;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***vz, dx, dy, lflux, gflux, L, A, dp, eta, vf, pgrad, K;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context variables
	FDSTAG    *fs      = nlctx->fs;
	JacResCtx *jrctx   = nlctx->jrctx;
	Material_t *phases = jrctx->phases;

	// access z-velocity vector
	ierr = DMDAVecGetArray(fs->DA_Z, jrctx->lvz, &vz);  CHKERRQ(ierr);

	// compute local part of fluid volume flux [m^3/s]
	// approximate integral of abs(vz) over xy-plane at z=0 (outflux face)

	lflux = 0.0;

	//---------
	// Z-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// integrate over outflux face only
		if(k == 0)
		{
			// get local mesh sizes
			dx = SIZE_CELL(i, sx, fs->dsx);
			dy = SIZE_CELL(j, sy, fs->dsy);

			// update integral
			lflux += PetscAbsScalar(vz[k][j][i])*dx*dy;
		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_Z, jrctx->lvz, &vz);  CHKERRQ(ierr);

	// compute global flux
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lflux, &gflux, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gflux = lflux;
	}

	// get length of the specimen along the flow direction
	L = user->H;

	// get area of outflux face
	A = user->W*user->L;

	// get applied pressure difference
	ierr = PetscOptionsGetScalar(NULL, "-pgrad", &dp, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) dp = 1.0;

	// get fluid viscosity (fluid phase is #1)
	eta = 1.0/(2.0*phases[1].Bd);

	// ***

	// compute average fluid velocity (normalized by outlux area)
	vf = gflux/A;

	// compute pressure gradient (normalized by length along flow direction)
	pgrad = dp/L;

	// compute permeability
	K = vf*eta/pgrad;

	// ***

	// output to the screen and to the file
	PetscPrintf(PETSC_COMM_WORLD,"# ==============================================\n");
	PetscPrintf(PETSC_COMM_WORLD,"# EFFECTIVE PERMEABILITY CONSTANT: %E\n", K);
	PetscPrintf(PETSC_COMM_WORLD,"# ==============================================\n");

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		db = fopen("darcy.dat", "w");

		fprintf(db,"# ==============================================\n");
		fprintf(db,"# EFFECTIVE PERMEABILITY CONSTANT: %E\n", K);
		fprintf(db,"# ==============================================\n");

		fclose(db);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DoDarcyTests"
PetscErrorCode DoDarcyTests(NLCtx *nlctx, UserCtx *user)
{

	FDSTAG    *fs      = nlctx->fs;
	BCCtx     *sbc     = nlctx->sbc;
	JacResCtx *jrctx   = nlctx->jrctx;
	BlockMat  *bmat    = nlctx->bmat;
	PetscBool  isset;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsHasName(NULL, "-darcy_test", &isset); CHKERRQ(ierr);
	if(isset != PETSC_TRUE) PetscFunctionReturn(0);

	// assemble matrix & rhs
	ierr = VecSet(jrctx->gsol, 0.0); CHKERRQ(ierr);
	ierr = VecSet(jrctx->gres, 0.0); CHKERRQ(ierr);
	ierr = FDSTAGFormResidual(NULL, jrctx->gsol, jrctx->gres, nlctx); CHKERRQ(ierr);
	ierr = VecScale(jrctx->gres, -1.0); CHKERRQ(ierr);

	ierr = BlockMatCompute(bmat, fs, sbc, jrctx); CHKERRQ(ierr);

	ierr = PowellHestenes(bmat, jrctx->gres, jrctx->gsol); CHKERRQ(ierr);

	ierr = FDSTAGFormResidual(NULL, jrctx->gsol, jrctx->gres, nlctx); CHKERRQ(ierr);

	// compute & output permeability
	ierr = DarcyPostProcess(nlctx, user); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DoMGTests"
PetscErrorCode DoMGTests(NLCtx *nlctx, PVOut *pvout)
{

	PetscViewer binViewer;
	PetscBool   isset;
	Vec         InitGuess;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(pvout) pvout = NULL;

	JacResCtx *jrctx = nlctx->jrctx;

	ierr = VecDuplicate(jrctx->gsol, &InitGuess); CHKERRQ(ierr);
	ierr = VecZeroEntries(InitGuess);             CHKERRQ(ierr);

	// load initial guess
	ierr = PetscOptionsHasName(NULL, "-load_init", &isset); CHKERRQ(ierr);

	if(isset == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "=====================\n");
		PetscPrintf(PETSC_COMM_WORLD, "LOADING INITIAL GUESS\n");
		PetscPrintf(PETSC_COMM_WORLD, "=====================\n");

		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "init.bin", FILE_MODE_READ, &binViewer); CHKERRQ(ierr);

		ierr = VecLoad(InitGuess, binViewer);

		PetscViewerDestroy(&binViewer);

	}

	// coupled test
	ierr = PetscOptionsHasName(NULL, "-test_coupled", &isset); CHKERRQ(ierr);

	if(isset == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "======================\n");
		PetscPrintf(PETSC_COMM_WORLD, "COUPLED MULTIGRID TEST\n");
		PetscPrintf(PETSC_COMM_WORLD, "======================\n");

//		ierr = MGTest(nlctx, pvout, InitGuess, 0.0, 0);  CHKERRQ(ierr);

	}


	// uncoupled test
	ierr = PetscOptionsHasName(NULL, "-test_uncoupled", &isset); CHKERRQ(ierr);

	if(isset == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD, "========================\n");
		PetscPrintf(PETSC_COMM_WORLD, "UNCOUPLED MULTIGRID TEST\n");
		PetscPrintf(PETSC_COMM_WORLD, "========================\n");

//		ierr = FieldSplitTest(nlctx, pvout, InitGuess, 1.0, 1); CHKERRQ(ierr);

	}

	// write solution
	ierr = PetscOptionsHasName(NULL, "-dump_init", &isset); CHKERRQ(ierr);

	if(isset == PETSC_TRUE)
	{

		PetscPrintf(PETSC_COMM_WORLD, "================\n");
		PetscPrintf(PETSC_COMM_WORLD, "WRITING SOLUTION\n");
		PetscPrintf(PETSC_COMM_WORLD, "================\n");

		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "init.bin", FILE_MODE_WRITE, &binViewer); CHKERRQ(ierr);

		ierr = VecView(jrctx->gsol, binViewer);

		PetscViewerDestroy(&binViewer);

	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetLinRes"
PetscErrorCode GetLinRes(Mat A, Vec x, Vec rhs, Vec res)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatMult(A, x, res);     CHKERRQ(ierr);
	ierr = VecScale(res, -1.0);    CHKERRQ(ierr);
	ierr = VecAXPY(res, 1.0, rhs); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGOutputResTest"
PetscErrorCode FDSTAGOutputResTest(FDSTAG *fs, JacResCtx *jrctx, Vec f)
{

	PetscScalar *fx, *fy, *fz, *c, *res, *iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = VecGetArray(jrctx->gfx,  &fx);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gfy,  &fy);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gfz,  &fz);  CHKERRQ(ierr);
	ierr = VecGetArray(jrctx->gc,   &c);   CHKERRQ(ierr);
	ierr = VecGetArray(f, &res); CHKERRQ(ierr);

	// copy vectors component-wise
	iter = res;

	ierr  = PetscMemcpy(fx, iter, (size_t)fs->nXFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nXFace;

	ierr  = PetscMemcpy(fy, iter, (size_t)fs->nYFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nYFace;

	ierr  = PetscMemcpy(fz, iter, (size_t)fs->nZFace*sizeof(PetscScalar)); CHKERRQ(ierr);
	iter += fs->nZFace;

	ierr  = PetscMemcpy(c,  iter, (size_t)fs->nCells*sizeof(PetscScalar)); CHKERRQ(ierr);

	// restore access
	ierr = VecRestoreArray(jrctx->gfx,  &fx);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gfy,  &fy);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gfz,  &fz);  CHKERRQ(ierr);
	ierr = VecRestoreArray(jrctx->gc,   &c);   CHKERRQ(ierr);
	ierr = VecRestoreArray(f, &res); CHKERRQ(ierr);


	PetscScalar xnorm, ynorm, znorm ,cnorm;

	ierr = VecNorm(jrctx->gfx, NORM_2, &xnorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gfy, NORM_2, &ynorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gfz, NORM_2, &znorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gc,  NORM_2, &cnorm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "*** xnorm: %12.12e ynorm: %12.12e znorm: %12.12e cnorm: %12.12e***", xnorm, ynorm, znorm, cnorm);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BlockSetUpMGVelBlock"
PetscErrorCode BlockSetUpMGVelBlock(NLCtx *nlctx, MGCtx *mg)
{
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FieldSplitTest"
PetscErrorCode FieldSplitTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime)
{
	KSP        ksp;
	PC         pc;
	Vec        res;
	MGCtx      mg;
	FDSTAG    *fs      = nlctx->fs;
	JacResCtx *jrctx   = nlctx->jrctx;
	BCCtx     *cbc     = nlctx->cbc;
	BlockMat  *bmat    = nlctx->bmat;
	char      *DirName = NULL;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// CHECK FILEDSPLIT

	// assemble matrix & rhs
	ierr = VecSet(jrctx->gsol, 0.0); CHKERRQ(ierr);
	ierr = VecSet(jrctx->gres, 0.0); CHKERRQ(ierr);
	ierr = FDSTAGFormResidual(NULL, jrctx->gsol, jrctx->gres, nlctx); CHKERRQ(ierr);
	ierr = VecScale(jrctx->gres, -1.0);
	ierr = BlockMatCompute(bmat, fs, cbc, jrctx); CHKERRQ(ierr);

	// setup linear solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, bmat->A, bmat->P);
	ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

	// setup preconditioner
	ierr = KSPGetPC(ksp, &pc);                    CHKERRQ(ierr);
	ierr = PCSetType(pc, PCSHELL);                CHKERRQ(ierr);
	ierr = PCShellSetContext(pc, bmat);           CHKERRQ(ierr);
//	ierr = PCShellSetApply(pc, &ApplyFieldSplit); CHKERRQ(ierr);

	// setup multigrid
	ierr = BlockSetUpMGVelBlock(nlctx, &mg);      CHKERRQ(ierr);

	// load initial guess
	ierr = VecCopy(InitGuess, jrctx->gsol); CHKERRQ(ierr);

	// solve equations
	ierr = KSPSolve(ksp, jrctx->gres, jrctx->gsol); CHKERRQ(ierr);

	// copy results for output
	ierr = FDSTAGCopySol(fs, cbc, jrctx, jrctx->gsol); CHKERRQ(ierr);

	// compute & store linear residual
	ierr = VecDuplicate(jrctx->gres, &res);                   CHKERRQ(ierr);
	ierr = GetLinRes(bmat->A, jrctx->gsol, jrctx->gres, res); CHKERRQ(ierr);
	ierr = FDSTAGOutputResTest(fs, jrctx, res);                 CHKERRQ(ierr);

	// create directory
	asprintf(&DirName, "Timestep_%1.6lld",(LLD)itime);
	ierr = FDSTAGCreateOutputDirectory(DirName); CHKERRQ(ierr);

	// Paraview output FDSTAG fields
	ierr = PVOutWriteTimeStep(pvout, jrctx, time, itime); CHKERRQ(ierr);

	// clean up
	if(DirName) free(DirName);
	ierr = KSPDestroy(&ksp);  CHKERRQ(ierr);
	ierr = VecDestroy(&res);  CHKERRQ(ierr);
	ierr = MGCtxDestroy(&mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGTest"
PetscErrorCode MGTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime)
{
	KSP        ksp;
	PC         pc;
	Vec        res;
	MGCtx      mg;
	FDSTAG    *fs      = nlctx->fs;
	BCCtx     *cbc     = nlctx->cbc;
	JacResCtx *jrctx   = nlctx->jrctx;
	BlockMat  *bmat    = nlctx->bmat;
	char      *DirName = NULL;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// CHECK MULTIGRID

	// assemble matrix & rhs
	ierr = VecSet(jrctx->gsol, 0.0); CHKERRQ(ierr);
	ierr = VecSet(jrctx->gres, 0.0); CHKERRQ(ierr);
	ierr = FDSTAGFormResidual(NULL, jrctx->gsol, jrctx->gres, nlctx); CHKERRQ(ierr);
	ierr = VecScale(jrctx->gres, -1.0);
	ierr = BlockMatCompute(bmat, fs, cbc, jrctx); CHKERRQ(ierr);

	// setup linear solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, bmat->A, bmat->P);
	ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

	// setup COUPLED MG preconditioner
	ierr = KSPGetPC(ksp, &pc);                   CHKERRQ(ierr);
	ierr = MGCtxCreate(&mg, fs, pc, IDXCOUPLED); CHKERRQ(ierr);
	ierr = MGCtxSetup(&mg, fs, bc, IDXCOUPLED);  CHKERRQ(ierr);
	ierr = MGCtxSetDiagOnLevels(&mg, pc);        CHKERRQ(ierr);

	// load initial guess
	ierr = VecCopy(InitGuess, jrctx->gsol); CHKERRQ(ierr);

	// solve equations
	ierr = KSPSolve(ksp, jrctx->gres, jrctx->gsol); CHKERRQ(ierr);

	// copy results for output
	ierr = FDSTAGCopySol(fs, bc, jrctx, jrctx->gsol); CHKERRQ(ierr);

	// compute & store linear residual
	ierr = VecDuplicate(jrctx->gres, &res);                   CHKERRQ(ierr);
	ierr = GetLinRes(bmat->A, jrctx->gsol, jrctx->gres, res); CHKERRQ(ierr);
	ierr = FDSTAGOutputResTest(fs, jrctx, res);                 CHKERRQ(ierr);

	// create directory
	asprintf(&DirName, "Timestep_%1.6lld",(LLD)itime);
	ierr = FDSTAGCreateOutputDirectory(DirName); CHKERRQ(ierr);

	// Paraview output FDSTAG fields
	ierr = PVOutWriteTimeStep(pvout, jrctx, time, itime); CHKERRQ(ierr);

	// clean up
	if(DirName) free(DirName);
	ierr = KSPDestroy(&ksp);
	ierr = MGCtxDestroy(&mg);
	ierr = VecDestroy(&res); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar InterpolateLinear3D(PetscScalar cx, PetscScalar cy, PetscScalar cz,  BCValues BC)
{
	PetscScalar xb, xe, yb, ye, zb, ze;

	// get relative coordinates
	xe = (cx - BC.cxb)/(BC.cxe - BC.cxb); xb = 1.0 - xe;
	ye = (cy - BC.cyb)/(BC.cye - BC.cyb); yb = 1.0 - ye;
	ze = (cz - BC.czb)/(BC.cze - BC.czb); zb = 1.0 - ze;

	return
	BC.A[0]*xb*yb*zb +
	BC.A[1]*xe*yb*zb +
	BC.A[2]*xb*ye*zb +
	BC.A[3]*xe*ye*zb +
	BC.A[4]*xb*yb*ze +
	BC.A[5]*xe*yb*ze +
	BC.A[6]*xb*ye*ze +
	BC.A[7]*xe*ye*ze;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InitVelocityTest"
PetscErrorCode InitVelocityTest(
	JacRes      *jr,
	UserCtx     *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal)
{
	// Initialize velocity vectors to test strain rates.
	// Selects the velocity direction x, y, z (0, 1, 2) with vectDir,
	// and applies constant gradient in the gradDir direction
	// between the values begVal & endVal (in the positive direction of gradDir).
	// Initializes boundary ghost points in the tangential directions accordingly.

	DM          DA;
	BCValues    BC;
	Vec         lvec, gvec;
	PetscScalar ***larr,  ***garr, cx, cy, cz, bcx, bcy, bcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG *fs = jr->fs;

	// initialize maximal index in all directions
	mx = fs->dsx.tnods - 2;
	my = fs->dsy.tnods - 2;
	mz = fs->dsz.tnods - 2;

	// assign interpolation data
	BC.cxb = usr->x_left;  BC.cxe = BC.cxb + usr->W;
	BC.cyb = usr->y_front; BC.cye = BC.cyb + usr->L;
	BC.czb = usr->z_bot;   BC.cze = BC.czb + usr->H;

	if(gradDir == 0)
	{
		BC.A[0] = BC.A[2] = BC.A[4] = BC.A[6] = begVal;
		BC.A[1] = BC.A[3] = BC.A[5] = BC.A[7] = endVal;
	}
	if(gradDir == 1)
	{
		BC.A[0] = BC.A[1] = BC.A[4] = BC.A[5] = begVal;
		BC.A[2] = BC.A[3] = BC.A[6] = BC.A[7] = endVal;
	}
	if(gradDir == 2)
	{
		BC.A[0] = BC.A[1] = BC.A[2] = BC.A[3] = begVal;
		BC.A[4] = BC.A[5] = BC.A[6] = BC.A[7] = endVal;
	}

	// get DA & vectors
	if(vectDir == 0) { DA = fs->DA_X; lvec = jrctx->lvx; gvec = jrctx->gvx; }
	if(vectDir == 1) { DA = fs->DA_Y; lvec = jrctx->lvy; gvec = jrctx->gvy; }
	if(vectDir == 2) { DA = fs->DA_Z; lvec = jrctx->lvz; gvec = jrctx->gvz; }

	// access vectors
	ierr = DMDAVecGetArray(DA, lvec, &larr); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(DA, gvec, &garr); CHKERRQ(ierr);

	// get loop bounds
	if(vectDir == 0) GET_NODE_RANGE(nx, sx, fs->dsx) else GET_CELL_RANGE(nx, sx, fs->dsx)
	if(vectDir == 1) GET_NODE_RANGE(ny, sy, fs->dsy) else GET_CELL_RANGE(ny, sy, fs->dsy)
	if(vectDir == 2) GET_NODE_RANGE(nz, sz, fs->dsz) else GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// get coordinates
		if(vectDir == 0) cx = COORD_NODE(i, sx, fs->dsx); else cx = COORD_CELL(i, sx, fs->dsx);
		if(vectDir == 1) cy = COORD_NODE(j, sy, fs->dsy); else cy = COORD_CELL(j, sy, fs->dsy);
		if(vectDir == 2) cz = COORD_NODE(k, sz, fs->dsz); else cz = COORD_CELL(k, sz, fs->dsz);

		// initialize global array
		garr[k][j][i] = InterpolateLinear3D(cx, cy, cz, BC);

		// set boundary points in local array
		if(vectDir == 0)
		{
			if(j == 0)	{ bcy = COORD_CELL(j-1, sy, fs->dsy); larr[k][j-1][i] = InterpolateLinear3D(cx, bcy, cz,  BC); }
			if(j == my) { bcy = COORD_CELL(j+1, sy, fs->dsy); larr[k][j+1][i] = InterpolateLinear3D(cx, bcy, cz,  BC); }
			if(k == 0)	{ bcz = COORD_CELL(k-1, sz, fs->dsz); larr[k-1][j][i] = InterpolateLinear3D(cx, cy,  bcz, BC); }
			if(k == mz) { bcz = COORD_CELL(k+1, sz, fs->dsz); larr[k+1][j][i] = InterpolateLinear3D(cx, cy,  bcz, BC); }
		}

		if(vectDir == 1)
		{
			if(i == 0)  { bcx = COORD_CELL(i-1, sx, fs->dsx); larr[k][j][i-1] = InterpolateLinear3D(bcx, cy, cz,  BC); }
			if(i == mx) { bcx = COORD_CELL(i+1, sx, fs->dsx); larr[k][j][i+1] = InterpolateLinear3D(bcx, cy, cz,  BC); }
			if(k == 0)  { bcz = COORD_CELL(k-1, sz, fs->dsz); larr[k-1][j][i] = InterpolateLinear3D(cx,  cy, bcz, BC); }
			if(k == mz) { bcz = COORD_CELL(k+1, sz, fs->dsz); larr[k+1][j][i] = InterpolateLinear3D(cx,  cy, bcz, BC); }
		}

		if(vectDir == 2)
		{
			if(i == 0)  { bcx = COORD_CELL(i-1, sx, fs->dsx); larr[k][j][i-1] = InterpolateLinear3D(bcx, cy,  cz, BC); }
			if(i == mx) { bcx = COORD_CELL(i+1, sx, fs->dsx); larr[k][j][i+1] = InterpolateLinear3D(bcx, cy,  cz, BC); }
			if(j == 0)  { bcy = COORD_CELL(j-1, sy, fs->dsy); larr[k][j-1][i] = InterpolateLinear3D(cx,  bcy, cz, BC); }
			if(j == my) { bcy = COORD_CELL(j+1, sy, fs->dsy); larr[k][j+1][i] = InterpolateLinear3D(cx,  bcy, cz, BC); }
		}

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(DA, lvec, &larr); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(DA, gvec, &garr); CHKERRQ(ierr);

	// initialize internal ghost points
	GLOBAL_TO_LOCAL(DA, gvec, lvec)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCtxClearVelocity"
PetscErrorCode JacResCtxClearVelocity(JacRes *jr, PetscInt vectDir)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	if(vectDir == 0)
	{
		ierr = VecSet(jr->lvx, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jr->gvx, 0.0); CHKERRQ(ierr);
	}
	if(vectDir == 1)
	{
		ierr = VecSet(jr->lvy, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jr->gvy, 0.0); CHKERRQ(ierr);
	}
	if(vectDir == 2)
	{
		ierr = VecSet(jr->lvz, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jr->gvz, 0.0); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "StrainRateSingleComp"
PetscErrorCode StrainRateSingleComp(
	JacRes      *jr,
	UserCtx     *usr,
	PVOut       *pvout,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  ttime,
	PetscInt     itime)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	FDSTAG *fs = jr->fs;

	// initialize velocity
	ierr = InitVelocityTest(fs, jr, usr, vectDir, gradDir, -1.0, 1.0); CHKERRQ(ierr);

	// compute effective strain rate & invariant
	ierr = FDSTAGetEffStrainRate(fs, jr); CHKERRQ(ierr);

	// compute residual
	ierr = FDSTAGetResidual(fs, jr); CHKERRQ(ierr);

	// create output directory
	char *DirectoryName = NULL;
	asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
	ierr = FDSTAGCreateOutputDirectory(DirectoryName); CHKERRQ(ierr);
	if(DirectoryName) free(DirectoryName);

	// save output
	ierr = PVOutWriteTimeStep(pvout, jr, ttime, itime); CHKERRQ(ierr);

	// clear velocities
	ierr = JacResCtxClearVelocity(jr, vectDir); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "StrainRateInterpTest"
PetscErrorCode StrainRateInterpTest(
	JacRes      *jr,
	UserCtx     *usr,
	PVOut       *pvout)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// Initialize velocity field to generate single non-zero component
	// of the velocity gradient. Strain rate field and second invariant
	// should give obvious & predictable result in sequence & in parallel.

	FDSTAG *fs = jr->fs;

	// dvx/dx
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 0, 0, 1.0, 11); CHKERRQ(ierr);

	// dvx/dy
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 0, 1, 2.0, 12); CHKERRQ(ierr);

	// dvx/dz
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 0, 2, 3.0, 13); CHKERRQ(ierr);


	// dvy/dx
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 1, 0, 4.0, 21); CHKERRQ(ierr);

	// dvy/dy
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 1, 1, 5.0, 22); CHKERRQ(ierr);

	// dvy/dz
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 1, 2, 6.0, 23); CHKERRQ(ierr);


	// dvz/dx
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 2, 0, 7.0, 31); CHKERRQ(ierr);

	// dvz/dy
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 2, 1, 8.0, 32); CHKERRQ(ierr);

	// dvz/dz
	ierr = StrainRateSingleComp(fs, jr, usr, pvout, 2, 2, 9.0, 33); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
PetscScalar InterpolateLinear3D(PetscScalar cx, PetscScalar cy, PetscScalar cz,  BCValues BC)
{
	PetscScalar xb, xe, yb, ye, zb, ze;

	// get relative coordinates
	xe = (cx - BC.cxb)/(BC.cxe - BC.cxb); xb = 1.0 - xe;
	ye = (cy - BC.cyb)/(BC.cye - BC.cyb); yb = 1.0 - ye;
	ze = (cz - BC.czb)/(BC.cze - BC.czb); zb = 1.0 - ze;

	return
	BC.A[0]*xb*yb*zb +
	BC.A[1]*xe*yb*zb +
	BC.A[2]*xb*ye*zb +
	BC.A[3]*xe*ye*zb +
	BC.A[4]*xb*yb*ze +
	BC.A[5]*xe*yb*ze +
	BC.A[6]*xb*ye*ze +
	BC.A[7]*xe*ye*ze;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InitVelocityTest"
PetscErrorCode InitVelocityTest(
	JacRes      *jr,
	UserCtx     *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal,
	Tensor2RN    *L)
{
	// Initialize velocity vectors to test velocity gradients
	// Selects the velocity direction x, y, z (0, 1, 2) with vectDir,
	// and applies constant gradient in the gradDir direction
	// between the values begVal & endVal (in the positive direction of gradDir).
	// Initializes boundary ghost points in the tangential directions accordingly.

	DM          DA;
	BCValues    BC;
	Vec         lvec;
	PetscScalar ***larr, cx, cy, cz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscScalar l[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, d[3];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG *fs = jr->fs;

	// assign interpolation data
	BC.cxb = usr->x_left;  BC.cxe = BC.cxb + usr->W; d[0] = usr->W;
	BC.cyb = usr->y_front; BC.cye = BC.cyb + usr->L; d[1] = usr->L;
	BC.czb = usr->z_bot;   BC.cze = BC.czb + usr->H; d[2] = usr->H;

	if(gradDir == 0)
	{
		BC.A[0] = BC.A[2] = BC.A[4] = BC.A[6] = begVal;
		BC.A[1] = BC.A[3] = BC.A[5] = BC.A[7] = endVal;
	}
	if(gradDir == 1)
	{
		BC.A[0] = BC.A[1] = BC.A[4] = BC.A[5] = begVal;
		BC.A[2] = BC.A[3] = BC.A[6] = BC.A[7] = endVal;
	}
	if(gradDir == 2)
	{
		BC.A[0] = BC.A[1] = BC.A[2] = BC.A[3] = begVal;
		BC.A[4] = BC.A[5] = BC.A[6] = BC.A[7] = endVal;
	}

	// compute uniform gradient
	l[vectDir*3 + gradDir] = (endVal - begVal)/d[gradDir];

	L->xx = l[0]; L->xy = l[1]; L->xz = l[2];
	L->yx = l[3]; L->yy = l[4]; L->yz = l[5];
	L->zx = l[6]; L->zy = l[7]; L->zz = l[8];

	// get DA & vectors
	if(vectDir == 0) { DA = fs->DA_X; lvec = jr->lvx; }
	if(vectDir == 1) { DA = fs->DA_Y; lvec = jr->lvy; }
	if(vectDir == 2) { DA = fs->DA_Z; lvec = jr->lvz; }

	// access vectors
	ierr = DMDAVecGetArray(DA, lvec, &larr); CHKERRQ(ierr);

	// get loop bounds
	ierr = DMDAGetGhostCorners(DA, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get coordinates
		if(vectDir == 0) cx = fs->dsx.nbuff[(i-sx)]; else cx = fs->dsx.cbuff[(i-sx)];
		if(vectDir == 1) cy = fs->dsy.nbuff[(j-sy)]; else cy = fs->dsy.cbuff[(j-sy)];
		if(vectDir == 2) cz = fs->dsz.nbuff[(k-sz)]; else cz = fs->dsz.cbuff[(k-sz)];

		// initialize velocity array
		larr[k][j][i] = InterpolateLinear3D(cx, cy, cz, BC);
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(DA, lvec, &larr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ClearVelocity"
PetscErrorCode ClearVelocity(JacRes *jr)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	ierr = VecSet(jr->lvx, 0.0); CHKERRQ(ierr);
	ierr = VecSet(jr->lvy, 0.0); CHKERRQ(ierr);
	ierr = VecSet(jr->lvz, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VelGradSingleComp"
PetscErrorCode VelGradSingleComp(
	JacRes      *jr,
	UserCtx     *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal)
{
	FDSTAG      *fs;
	Tensor2RN   L, LRef;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar vel[3];
	PetscScalar ***lvx, ***lvy, ***lvz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// initialize velocity component
	ierr = InitVelocityTest(jr, usr, vectDir, gradDir, begVal, endVal, &LRef); CHKERRQ(ierr);

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &lvz); CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get gradient and velocities at cell center
		ierr = getGradientVel(fs, lvx, lvy, lvz, i, j, k, sx, sy, sz, &L, vel, NULL);

		if(!Tensor2RNCheckEq(&L, &LRef, 1e-12))
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Velocity gradient test failed");
		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VelGradTest"
PetscErrorCode VelGradTest(
	JacRes      *jr,
	UserCtx     *usr)
{
	PetscInt   comp;
	PetscBool found;

	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	PetscOptionsGetInt(NULL ,"-grad_comp", &comp, &found);

	if(found == PETSC_FALSE) PetscFunctionReturn(0);

	// Initialize velocity field to generate single non-zero component.

	if     (comp == 11)
	{
		// dvx/dx
		ierr = VelGradSingleComp(jr, usr, 0, 0, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 12)
	{
		// dvx/dy
		ierr = VelGradSingleComp(jr, usr, 0, 1, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 13)
	{
		// dvx/dz
		ierr = VelGradSingleComp(jr, usr, 0, 2, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 21)
	{
		// dvy/dx
		ierr = VelGradSingleComp(jr, usr, 1, 0, 0.0, 1.0); CHKERRQ(ierr);

	}
	else if(comp == 22)
	{
		// dvy/dy
		ierr = VelGradSingleComp(jr, usr, 1, 1, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 23)
	{
		// dvy/dz
		ierr = VelGradSingleComp(jr, usr, 1, 2, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 31)
	{
		// dvz/dx
		ierr = VelGradSingleComp(jr, usr, 2, 0, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 32)
	{
		// dvz/dy
		ierr = VelGradSingleComp(jr, usr, 2, 1, 0.0, 1.0); CHKERRQ(ierr);
	}
	else if(comp == 33)
	{
		// dvz/dz
		ierr = VelGradSingleComp(jr, usr, 2, 2, 0.0, 1.0); CHKERRQ(ierr);
	}
	else
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Wrong value for -grad_comp option: %lld\n", (LLD)comp);
	}

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
