//---------------------------------------------------------------------------
//.........   LaMEM - FDSTAG CANONICAL INTERFACE CHECKING ROUTINES   ........
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "ParaViewOutput.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "multigrid.h"
#include "check_fdstag.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DarcyPostProcess"
PetscErrorCode DarcyPostProcess(NLCtx *nlctx, UserContext *user)
{
	FILE        *db;
	PetscBool   flg;
	PetscMPIInt nproc, iproc;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***vz, dx, dy, lflux, gflux, L, A, dp, eta, vf, pgrad, K;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context variables
	FDSTAG    *fs      = nlctx->fs;
	JacResCtx *jrctx   = nlctx->jrctx;
	Material_t *phases = jrctx->phases;

	// get number of processors & rank
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &iproc); CHKERRQ(ierr);

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
	if(nproc != 1)
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
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-pgrad", &dp, &flg); CHKERRQ(ierr);
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

	if(iproc == 0)
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
PetscErrorCode DoDarcyTests(NLCtx *nlctx, UserContext *user)
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

/*
	PetscScalar xnorm, ynorm, znorm ,cnorm;

	ierr = VecNorm(jrctx->gfx, NORM_2, &xnorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gfy, NORM_2, &ynorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gfz, NORM_2, &znorm); CHKERRQ(ierr);
	ierr = VecNorm(jrctx->gc,  NORM_2, &cnorm); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "*** xnorm: %12.12e ynorm: %12.12e znorm: %12.12e cnorm: %12.12e***", xnorm, ynorm, znorm, cnorm);
*/
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
/*
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
	ierr = LaMEM_CreateOutputDirectory(DirName); CHKERRQ(ierr);

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
	ierr = LaMEM_CreateOutputDirectory(DirName); CHKERRQ(ierr);

	// Paraview output FDSTAG fields
	ierr = PVOutWriteTimeStep(pvout, jrctx, time, itime); CHKERRQ(ierr);

	// clean up
	if(DirName) free(DirName);
	ierr = KSPDestroy(&ksp);
	ierr = MGCtxDestroy(&mg);
	ierr = VecDestroy(&res); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
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
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
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
PetscErrorCode JacResCtxClearVelocity(JacResCtx *jrctx, PetscInt vectDir)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	if(vectDir == 0)
	{
		ierr = VecSet(jrctx->lvx, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jrctx->gvx, 0.0); CHKERRQ(ierr);
	}
	if(vectDir == 1)
	{
		ierr = VecSet(jrctx->lvy, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jrctx->gvy, 0.0); CHKERRQ(ierr);
	}
	if(vectDir == 2)
	{
		ierr = VecSet(jrctx->lvz, 0.0); CHKERRQ(ierr);
		ierr = VecSet(jrctx->gvz, 0.0); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "StrainRateSingleComp"
PetscErrorCode StrainRateSingleComp(
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
	PVOut       *pvout,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  ttime,
	PetscInt     itime)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// initialize velocity
	ierr = InitVelocityTest(fs, jrctx, usr, vectDir, gradDir, -1.0, 1.0); CHKERRQ(ierr);

	// compute effective strain rate & invariant
	ierr = FDSTAGetEffStrainRate(fs, jrctx); CHKERRQ(ierr);

	// compute residual
	ierr = FDSTAGetResidual(fs, jrctx); CHKERRQ(ierr);

	// create output directory
	char *DirectoryName = NULL;
	asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
	ierr = LaMEM_CreateOutputDirectory(DirectoryName); CHKERRQ(ierr);
	if(DirectoryName) free(DirectoryName);

	// save output
	ierr = PVOutWriteTimeStep(pvout, jrctx, ttime, itime); CHKERRQ(ierr);

	// clear velocities
	ierr = JacResCtxClearVelocity(jrctx, vectDir); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "StrainRateInterpTest"
PetscErrorCode StrainRateInterpTest(
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
	PVOut       *pvout)
{
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// Initialize velocity field to generate single non-zero component
	// of the velocity gradient. Strain rate field and second invariant
	// should give obvious & predictable result in sequence & in parallel.

	// dvx/dx
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 0, 0, 1.0, 11); CHKERRQ(ierr);

	// dvx/dy
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 0, 1, 2.0, 12); CHKERRQ(ierr);

	// dvx/dz
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 0, 2, 3.0, 13); CHKERRQ(ierr);


	// dvy/dx
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 1, 0, 4.0, 21); CHKERRQ(ierr);

	// dvy/dy
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 1, 1, 5.0, 22); CHKERRQ(ierr);

	// dvy/dz
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 1, 2, 6.0, 23); CHKERRQ(ierr);


	// dvz/dx
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 2, 0, 7.0, 31); CHKERRQ(ierr);

	// dvz/dy
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 2, 1, 8.0, 32); CHKERRQ(ierr);

	// dvz/dz
	ierr = StrainRateSingleComp(fs, jrctx, usr, pvout, 2, 2, 9.0, 33); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGetResMag"
PetscErrorCode FDSTAGetResMag(FDSTAG *fs, UserContext *usr, Vec f, Vec fmag)
{
	Field       ***pf;
	PetscScalar ***pfmag, fx, fy, fz, res, max, tol;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, im, jm, km, cnt;

	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// initialize
	tol =  1e-10;
	max = -1.0;
	im  = -1;
	jm  = -1;
	km  = -1;
	cnt =  0;

	ierr = DMDAVecGetArray(usr->DA_Vel,            f,    &pf);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(usr->FDSTAG.DA_CORNER,  fmag, &pfmag); CHKERRQ(ierr);

	GET_NODE_RANGE(nx, sx, fs->dsx);
	GET_NODE_RANGE(ny, sy, fs->dsy);
	GET_NODE_RANGE(nz, sz, fs->dsz);

	START_STD_LOOP
	{
		fx  = pf[k][j][i].Vx;
		fy  = pf[k][j][i].Vy;
		fz  = pf[k][j][i].Vz;
		res = sqrt(fx*fx + fy*fy + fz*fz);
		if(res > max)
		{	max = res;
			im  = i;
			jm  = j;
			km  = k;
		}
		if(res > tol) cnt++;
		pfmag[k][j][i] = res;
	}
	END_STD_LOOP

	PetscPrintf(PETSC_COMM_WORLD, "\n\n\n ");
	PetscPrintf(PETSC_COMM_WORLD, "***   MAXRES: %12.12e   I: %lld   J: %lld   K: %lld   ***", max, im, jm, km);
	PetscPrintf(PETSC_COMM_WORLD, "\n\n\n ");

	ierr = DMDAVecRestoreArray(usr->DA_Vel,            f,    &pf);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(usr->FDSTAG.DA_CORNER,  fmag, &pfmag); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCheckPhaseRatios"
PetscErrorCode FDSTAGCheckPhaseRatios(FDSTAG *fs, JacResCtx *jrctx)
{
	// check phase ratios in the FDSTAG data structures

	PetscScalar *phRat, tol, top, bot, sum;
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz, iter, numPhases;

	PetscFunctionBegin;

	// initialize
	numPhases = jrctx->numPhases;
	tol       = 1e-13;
	top       = 1.0 + tol;
	bot       = 1.0 - tol;
	//---------------
	// central points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svCell[iter++].svDev.phRat;
		sum = ArraySum(phRat, numPhases);
		if(sum < bot || sum > top) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "***   SUM: %12.12e   ***", sum);
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXYEdge[iter++].svDev.phRat;
		sum = ArraySum(phRat, numPhases);
		if(sum < bot || sum > top) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "***   SUM: %12.12e   ***", sum);
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXZEdge[iter++].svDev.phRat;
		sum = ArraySum(phRat, numPhases);
		if(sum < bot || sum > top) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "***   SUM: %12.12e   ***", sum);
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svYZEdge[iter++].svDev.phRat;
		sum = ArraySum(phRat, numPhases);
		if(sum < bot || sum > top) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "***   SUM: %12.12e   ***", sum);
	}
	END_STD_LOOP

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
PetscScalar ArraySum(PetscScalar *a, PetscInt n)
{
	PetscInt    i;
	PetscScalar sum = 0.0;
	for(i = 0; i < n; i++) sum += a[i];
	return sum;
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGDumpPhaseRatios"
PetscErrorCode FDSTAGDumpPhaseRatios(FDSTAG *fs, JacResCtx *jrctx)
{
	// check phase ratios in the FDSTAG data structures

	FILE        *db;
	PetscScalar *phRat, *buff, *ptr;
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, iter, numPhases, buffSz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	numPhases = jrctx->numPhases;

	buffSz = numPhases*(fs->nCells + fs->nXYEdg + fs->nXZEdg + fs->nYZEdg);
	ierr   = makeScalArray(&buff, NULL, buffSz);
	ptr    = buff;

	//---------------
	// central points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svCell[iter++].svDev.phRat;
		ArrayCopy(phRat, ptr, numPhases);
		ptr += numPhases;
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXYEdge[iter++].svDev.phRat;
		ArrayCopy(phRat, ptr, numPhases);
		ptr += numPhases;
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXZEdge[iter++].svDev.phRat;
		ArrayCopy(phRat, ptr, numPhases);
		ptr += numPhases;
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svYZEdge[iter++].svDev.phRat;
		ArrayCopy(phRat, ptr, numPhases);
		ptr += numPhases;
	}
	END_STD_LOOP

	// dump to file
	db = fopen("phRat.bin","w");
	fwrite(&buffSz, sizeof(PetscInt), 1, db);
	fwrite(buff, sizeof(PetscScalar), (size_t)buffSz, db);
	fclose(db);

	ierr = PetscFree(buff);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
void ArrayCopy(PetscScalar *a, PetscScalar *b, PetscInt n)
{
	PetscInt i;
	for(i = 0; i < n; i++) b[i] = a[i];
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGComparePhaseRatios"
PetscErrorCode FDSTAGComparePhaseRatios(FDSTAG *fs, JacResCtx *jrctx, PetscScalar rtol, PetscScalar atol)
{

	// check phase ratios in the FDSTAG data structures
	FILE        *db;
	PetscScalar *phRat, *phRatRef, *buff, *ptr;
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, iter, numPhases, buffSz, tx, ty, tz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	numPhases = jrctx->numPhases;

	// read reference from file
	db = fopen("phRat.bin","r");

	if(db == NULL) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "cannot open file phRat.bin");

	fread(&buffSz, sizeof(PetscInt), 1, db);

	ierr = makeScalArray(&buff, NULL, buffSz);

	fread(buff, sizeof(PetscScalar), (size_t)buffSz, db);

	fclose(db);

	ptr = buff;

	//---------------
	// central points
	//---------------
	iter = 0;
	GET_CELL_BOUNDS(nx, sx, fs->dsx)
	GET_CELL_BOUNDS(ny, sy, fs->dsy)
	GET_CELL_BOUNDS(nz, sz, fs->dsz)

	tx  = fs->dsx.tnods-1;
	ty  = fs->dsy.tnods-1;
	tz  = fs->dsz.tnods-1;

	START_STD_LOOP
	{
		phRat    = jrctx->svCell[iter++].svDev.phRat;
		phRatRef = getPtr(tx, ty, i, j, k, numPhases, ptr);
		ArrayCompare(phRat, phRatRef, numPhases, rtol, atol);

	}
	END_STD_LOOP

	ptr += numPhases*tx*ty*tz;

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_BOUNDS(nx, sx, fs->dsx)
	GET_NODE_BOUNDS(ny, sy, fs->dsy)
	GET_CELL_BOUNDS(nz, sz, fs->dsz)

	tx = fs->dsx.tnods;
	ty = fs->dsy.tnods;
	tz = fs->dsz.tnods-1;

	START_STD_LOOP
	{
		phRat    = jrctx->svXYEdge[iter++].svDev.phRat;
		phRatRef = getPtr(tx, ty, i, j, k, numPhases, ptr);
		ArrayCompare(phRat, phRatRef, numPhases, rtol, atol);

	}
	END_STD_LOOP

	ptr += numPhases*tx*ty*tz;

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_BOUNDS(nx, sx, fs->dsx)
	GET_CELL_BOUNDS(ny, sy, fs->dsy)
	GET_NODE_BOUNDS(nz, sz, fs->dsz)

	tx = fs->dsx.tnods;
	ty = fs->dsy.tnods-1;
	tz = fs->dsz.tnods;

	START_STD_LOOP
	{
		phRat    = jrctx->svXZEdge[iter++].svDev.phRat;
		phRatRef = getPtr(tx, ty, i, j, k, numPhases, ptr);
		ArrayCompare(phRat, phRatRef, numPhases, rtol, atol);

	}
	END_STD_LOOP

	ptr += numPhases*tx*ty*tz;

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_BOUNDS(nx, sx, fs->dsx)
	GET_NODE_BOUNDS(ny, sy, fs->dsy)
	GET_NODE_BOUNDS(nz, sz, fs->dsz)

	tx = fs->dsx.tnods-1;
	ty = fs->dsy.tnods;
	tz = fs->dsz.tnods;

	START_STD_LOOP
	{
		phRat    = jrctx->svYZEdge[iter++].svDev.phRat;
		phRatRef = getPtr(tx, ty, i, j, k, numPhases, ptr);
		ArrayCompare(phRat, phRatRef, numPhases, rtol, atol);

	}
	END_STD_LOOP

	ierr = PetscFree(buff); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
PetscScalar * getPtr(PetscInt nx, PetscInt ny, PetscInt i, PetscInt j, PetscInt k, PetscInt ndof, PetscScalar *a)
{
	return a + ndof*(k*nx*ny + j*nx + i);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "ArrayCompare"
PetscErrorCode ArrayCompare(PetscScalar *a, PetscScalar *b, PetscInt n, PetscScalar rtol, PetscScalar atol)
{
	PetscInt i;
	PetscFunctionBegin;
	for(i = 0; i < n; i++)
	{	if(!LAMEM_CHECKEQ(a[i], b[i], rtol, atol))
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "***   DIFF: %12.12e   ***", a[i]-b[i]);
		}
	}
	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGetResidualTest"
PetscErrorCode FDSTAGetResidualTest(FDSTAG *fs, JacResCtx *jrctx)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fx,  ***fy,  ***fz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear local residual vectors
	ierr = VecZeroEntries(jrctx->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->lfz); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_X, jrctx->lfx, &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jrctx->lfy, &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jrctx->lfz, &fz);  CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// momentum
		fx[k][j][i] += 1.0;   fx[k][j][i+1] -= 1.0;
		fy[k][j][i] += 1.0;   fy[k][j+1][i] -= 1.0;
		fz[k][j][i] += 1.0;   fz[k+1][j][i] -= 1.0;
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// momentum
		fx[k][j-1][i] += 1.0;   fx[k][j][i] -= 1.0;
		fy[k][j][i-1] += 1.0;   fy[k][j][i] -= 1.0;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points
	//-------------------------------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// momentum
		fx[k-1][j][i] += 1.0;   fx[k][j][i] -= 1.0;
		fz[k][j][i-1] += 1.0;   fz[k][j][i] -= 1.0;

	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points
	//-------------------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// update momentum residuals
		fy[k-1][j][i] += 1.0;   fy[k][j][i] -= 1.0;
		fz[k][j-1][i] += 1.0;   fz[k][j][i] -= 1.0;

	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_X, jrctx->lfx, &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jrctx->lfy, &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jrctx->lfz, &fz);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jrctx->lfx, jrctx->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jrctx->lfy, jrctx->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jrctx->lfz, jrctx->gfz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InitVec"
PetscErrorCode InitVec(DM da, Vec v)
{
	Field           ***va;
	PetscErrorCode  ierr;
	PetscInt        i, j, k, sx, sy, sz, nx, ny, nz, base;

	PetscFunctionBegin;

	ierr = VecZeroEntries(v); CHKERRQ(ierr);

	ierr = DMDAGetCorners(da, &sx,  &sy,  &sz,  &nx,  &ny,  &nz); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da, v, &va); CHKERRQ(ierr);

	START_STD_LOOP
	{
		base = (i+1)*(j+1) + (j+1)*(k+1) + (i+1)*(k+1) + (i+1)*(j+1)*(k+1);
		va[k][j][i].Vx = (PetscScalar)(base);
		va[k][j][i].Vy = (PetscScalar)(base+1);
		va[k][j][i].Vz = (PetscScalar)(base+2);
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(da, v, &va); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DumpVec"
PetscErrorCode DumpVec(DM da, Vec v)
{
	FILE            *db;
	Field           ***va;
	PetscScalar     *buff;
	PetscErrorCode  ierr;
	PetscInt        i, j, k, sx, sy, sz, nx, ny, nz, iter, buffSz;

	PetscFunctionBegin;

	iter = 0;

	ierr = DMDAGetCorners(da, &sx,  &sy,  &sz,  &nx,  &ny,  &nz); CHKERRQ(ierr);

	buffSz = 3*nx*ny*nz;

	ierr = makeScalArray(&buff, NULL, buffSz);

	ierr = DMDAVecGetArray(da, v, &va); CHKERRQ(ierr);

	START_STD_LOOP
	{
		buff[iter++] = va[k][j][i].Vx;
		buff[iter++] = va[k][j][i].Vy;
		buff[iter++] = va[k][j][i].Vz;
	}
	END_STD_LOOP

	// dump to file
	db = fopen("vec_test.bin","w");
	fwrite(&nx, sizeof(PetscInt), 1, db);
	fwrite(&ny, sizeof(PetscInt), 1, db);
	fwrite(&nz, sizeof(PetscInt), 1, db);
	fwrite(buff, sizeof(PetscScalar), (size_t)buffSz, db);
	fclose(db);

	ierr = PetscFree(buff); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da, v, &va); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "CountVec"
PetscErrorCode CountVec(DM da, Vec v, PetscScalar rtol, PetscScalar atol)
{
	FILE            *db;
	Field           ***va;
	PetscScalar     *buff, *vRef, pvec[3];
	PetscInt        i, j, k, sx, sy, sz, nx, ny, nz, tx, ty, tz, dcnt, cnt, buffSz;
	PetscErrorCode  ierr;
	PetscMPIInt     rank;
	MaxDiff         md;
	PetscFunctionBegin;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// read reference from file
	db = fopen("vec_test.bin","r");
	if(db == NULL) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "cannot open file vec_test.bin");
	fread(&tx, sizeof(PetscInt), 1, db);
	fread(&ty, sizeof(PetscInt), 1, db);
	fread(&tz, sizeof(PetscInt), 1, db);
	buffSz = 3*tx*ty*tz;
	ierr = makeScalArray(&buff, NULL, buffSz);
	fread(buff, sizeof(PetscScalar), (size_t)buffSz, db);
	fclose(db);

	// compare
	db = fopen("vec_comp.out","w");
	fprintf(db, "Velocity comparison database\n");
	cnt = 0;
	ierr = DMDAGetCorners(da, &sx,  &sy,  &sz,  &nx,  &ny,  &nz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, v, &va); CHKERRQ(ierr);
	START_STD_LOOP
	{
		vRef    = getPtr(tx, ty, i, j, k, 3, buff);
		pvec[0] = va[k][j][i].Vx;
		pvec[1] = va[k][j][i].Vy;
		pvec[2] = va[k][j][i].Vz;
		dcnt    = ArrayCount(vRef, pvec, 3, rtol, atol);
		cnt    += dcnt;
		if(dcnt)
		{	MDiff(vRef, pvec, &md, 3);
			fprintf(db, "rank=%lld  i=%lld  j=%lld  k=%lld  sec=%15.15e par=%15.15e, dif=%15.15e\n",(LLD)rank, (LLD)i, (LLD)j, (LLD)k, md.a, md.b, PetscAbsScalar(md.a-md.b));
			fflush(db);
		}
	}
	END_STD_LOOP
	fclose(db);

	ierr = PetscPrintf             (PETSC_COMM_WORLD, "\n\n\n"); CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf (PETSC_COMM_WORLD, "   *** rank = %lld   cnt = %5lld   ***\n", (LLD)rank, (LLD)cnt); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush  (PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = PetscPrintf             (PETSC_COMM_WORLD, "\n\n\n"); CHKERRQ(ierr);

	ierr = PetscFree(buff); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, v, &va); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
/*
PetscInt ArrayCount(PetscScalar *a, PetscScalar *b, PetscInt n, PetscScalar rtol, PetscScalar atol)
{
	PetscInt    i, cnt = 0;
	for(i = 0; i < n; i++)
	{
		if(!LAMEM_CHECKEQ(a[i], b[i], rtol, atol)) cnt++;
	}
	return cnt;
}
*/
//---------------------------------------------------------------------------
void MDiff(PetscScalar *a, PetscScalar *b, MaxDiff *md, PetscInt n)
{
	PetscInt    i;
	PetscScalar mdiff = DBL_MAX, diff;
	for(i = 0; i < n; i++)
	{	diff = PetscAbsScalar(a[i] - b[i]);
		if(diff && diff < mdiff)
		{	mdiff = diff;
			md->a = a[i];
			md->b = b[i];
		}
	}
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGDumpVelocities"
PetscErrorCode FDSTAGDumpVelocities(FDSTAG *fs, JacResCtx *jrctx)
{
	// dump global velocity vectors for comparison

	FILE        *db;
	PetscScalar *buff;
	PetscScalar ***vx,  ***vy,  ***vz;
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, iter, buffSz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	buffSz = fs->nXFace + fs->nYFace + fs->nZFace;
	ierr   = makeScalArray(&buff, NULL, buffSz);
	iter   = 0;

	ierr = DMDAVecGetArray(fs->DA_X, jrctx->gvx, &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jrctx->gvy, &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jrctx->gvz, &vz);  CHKERRQ(ierr);

	//---------------
	// x-faces
	//---------------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		 buff[iter++] = vx[k][j][i];
	}
	END_STD_LOOP

	//---------------
	// y - faces
	//---------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		 buff[iter++] = vy[k][j][i];
	}
	END_STD_LOOP

	//---------------
	// z - faces
	//---------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		 buff[iter++] = vz[k][j][i];
	}
	END_STD_LOOP

	// dump to file
	db = fopen("vel.bin","w");
	fwrite(&buffSz, sizeof(PetscInt), 1, db);
	fwrite(buff, sizeof(PetscScalar), (size_t)buffSz, db);
	fclose(db);

	ierr = PetscFree(buff);   CHKERRQ(ierr);

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X, jrctx->gvx, &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jrctx->gvy, &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jrctx->gvz, &vz);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCountVelocities"
PetscErrorCode FDSTAGCountVelocities(FDSTAG *fs, JacResCtx *jrctx, PetscScalar rtol, PetscScalar _atol)
{
	// compare global velocity vectors sequential & parallel versions

	FILE        *db;
	PetscScalar *buff, *ptr, *vRef, a, b;
	PetscScalar ***vx,  ***vy,  ***vz;
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, tx, ty, tz, buffSz;
	PetscInt     cx, cy, cz;
	PetscBool    nonzero;

	PetscMPIInt rank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	cx = 0;
	cy = 0;
	cz = 0;

	// read reference from file
	db = fopen("vel.bin","r");

	if(db == NULL) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "cannot open file vel.bin");

	fread(&buffSz, sizeof(PetscInt), 1, db);

	ierr = makeScalArray(&buff, NULL, buffSz);

	fread(buff, sizeof(PetscScalar), (size_t)buffSz, db);

	fclose(db);

	ptr = buff;

	// access velocities
	ierr = DMDAVecGetArray(fs->DA_X, jrctx->gvx, &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, jrctx->gvy, &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, jrctx->gvz, &vz);  CHKERRQ(ierr);

	nonzero = PETSC_FALSE;

	db = fopen("vec_comp.out","w");
	fprintf(db, "Velocity comparison database\n");

	//---------------
	// x-faces
	//---------------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	tx = fs->dsx.tnods;
	ty = fs->dsy.tnods-1;
	tz = fs->dsz.tnods-1;

	START_STD_LOOP
	{
		vRef = getPtr(tx, ty, i, j, k, 1, ptr);
		a    = (*vRef);
		b    = vx[k][j][i];
		if(a || b) nonzero = PETSC_TRUE;
		if(!LAMEM_CHECKEQ(a, b, rtol, _atol))
		{	fprintf(db, "x rank=%lld  i=%lld  j=%lld  k=%lld  sec=%15.15e par=%15.15e, dif=%15.15e\n",(LLD)rank, (LLD)i, (LLD)j, (LLD)k, a, b, PetscAbsScalar(a-b));
			fflush(db);
			cx++;
		}
	}
	END_STD_LOOP

	ptr += tx*ty*tz;

	//---------------
	// y - faces
	//---------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	tx = fs->dsx.tnods-1;
	ty = fs->dsy.tnods;
	tz = fs->dsz.tnods-1;

	START_STD_LOOP
	{
		vRef = getPtr(tx, ty, i, j, k, 1, ptr);
		a    = (*vRef);
		b    = vy[k][j][i];
		if(a || b) nonzero = PETSC_TRUE;
		if(!LAMEM_CHECKEQ(a, b, rtol, _atol))
		{	fprintf(db, "y rank=%lld  i=%lld  j=%lld  k=%lld  sec=%15.15e par=%15.15e, dif=%15.15e\n",(LLD)rank, (LLD)i, (LLD)j, (LLD)k, a, b, PetscAbsScalar(a-b));
			fflush(db);
			cy++;
		}
	}
	END_STD_LOOP

	ptr += tx*ty*tz;

	//---------------
	// z - faces
	//---------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	tx = fs->dsx.tnods-1;
	ty = fs->dsy.tnods-1;
//	tz = fs->dsz.tnods;

	START_STD_LOOP
	{
		vRef = getPtr(tx, ty, i, j, k, 1, ptr);
		a    = (*vRef);
		b    = vz[k][j][i];
		if(a || b) nonzero = PETSC_TRUE;
		if(!LAMEM_CHECKEQ(a, b, rtol, _atol))
		{	fprintf(db, "z rank=%lld  i=%lld  j=%lld  k=%lld  sec=%15.15e par=%15.15e, dif=%15.15e\n",(LLD)rank, (LLD)i, (LLD)j, (LLD)k, a, b, PetscAbsScalar(a-b));
			fflush(db);
			cz++;
		}
	}
	END_STD_LOOP

//	ptr += tx*ty*tz;

	fclose(db);

	if(nonzero == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD, "\n   *** ELEMENTS ARE NONZERO ***\n");

	ierr = PetscPrintf             (PETSC_COMM_WORLD, "\n\n\n");     CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf (PETSC_COMM_WORLD, "   *** rank = %lld   cx = %5lld   cy = %5lld   cz = %5lld   ***\n", (LLD)rank, (LLD)cx, (LLD)cy, (LLD)cz); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush  (PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
	ierr = PetscPrintf             (PETSC_COMM_WORLD, "\n\n\n");     CHKERRQ(ierr);


	ierr = PetscFree(buff);CHKERRQ(ierr);

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X, jrctx->gvx, &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, jrctx->gvy, &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, jrctx->gvz, &vz);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
