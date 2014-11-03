//---------------------------------------------------------------------------
//........................   TIME STEPPING SOLVER   .........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "Utils.h"
#include "tssolve.h"
//---------------------------------------------------------------------------
// * change time stepping for requested time interval, not integer number of steps
// * break-point files
// ...
//---------------------------------------------------------------------------
PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax)
{
	PetscScalar v, h, vmax, idt, idtmax;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, idx, ijk[3], jj, ln;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	idtmax = (*_idtmax);

	if(ds->h_uni != PETSC_TRUE)
	{
		// compute time step on variable spacing grid
		PetscScalar ***va;

		ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, gv, &va);                     CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get velocity
			v = va[k][j][i];

			// prepare node index buffer
			ijk[0] = i-sx;
			ijk[1] = j-sy;
			ijk[2] = k-sz;

			// anisotropic direction-dependent criterion
			if(v >= 0.0)  idx = ijk[dir];
			else          idx = ijk[dir]-1;

			// get mesh step
			h = ds->ncoor[idx+1] - ds->ncoor[idx];

			// get inverse time step (safe to compute)
			idt = v/h;

			// update maximum inverse time step
			if(idt > idtmax) idt = idtmax;
		}
		END_STD_LOOP

		ierr = DMDAVecRestoreArray(da, gv, &va); CHKERRQ(ierr);
	}
	else
	{
		// compute time step on uniform spacing grid
		PetscScalar *va;

		// get maximum local velocity
		ierr = VecGetLocalSize(gv, &ln); CHKERRQ(ierr);
		ierr = VecGetArray(gv, &va);     CHKERRQ(ierr);

		vmax = 0.0;
		for(jj = 0; jj < ln; jj++) { v = PetscAbsScalar(va[jj]); if(v > vmax) vmax = v;	}

		ierr = VecRestoreArray(gv, &va); CHKERRQ(ierr);

		// get inverse time step
		idt = vmax/ds->h;

		// update maximum inverse time step
		if(idt > idtmax) idt = idtmax;
	}

	// return result
	(*_idtmax) = idtmax;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolGetCourantStep"
PetscErrorCode TSSolGetCourantStep(TSSol *ts, JacRes *jr)
{
	//-------------------------------------
	// compute length of the next time step
	//-------------------------------------

	FDSTAG      *fs;
	PetscScalar dt, lidtmax, gidtmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs      = jr->fs;
	lidtmax = 0.0;

	// determine maximum local inverse time step
	ierr = getMaxInvStep1DLocal(&fs->dsx, fs->DA_X, jr->gvx, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsy, fs->DA_Y, jr->gvy, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsz, fs->DA_Z, jr->gvz, 0, &lidtmax); CHKERRQ(ierr);

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lidtmax, &gidtmax, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gidtmax = lidtmax;
	}

	// compute time step
	gidtmax /= ts->Cmax;

	if     (gidtmax < 1.0/ts->dtmax) dt = ts->dtmax;
	else if(gidtmax > 1.0/ts->dtmin) dt = ts->dtmin;
	else                             dt = 1.0/gidtmax;



	/*
		PetscScalar pdt;   // previous time step
		PetscScalar dt;    // current time step (to be defined)
		PetscScalar dtmin; // minimum time step
		PetscScalar dtmax; // maximum time step
		PetscScalar Cmax;  // dimensionless Courant number (should be {significantly} less than unit)
	*/


/*
	//--------------------------------
	// print velocity & time step info
	//--------------------------------

	if(user->Characteristic.Length > 1.0)
	{
		PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user->Characteristic.Velocity*user->Characteristic.cmYear);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [cm/year]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user->Characteristic.Velocity);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [m/s]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}

	if(user->Characteristic.Length > 1.0)
	{
		// Most likely a setup that runs in natural length scales
		PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ",user->time*user->Characteristic.Time/user->Characteristic.SecYear, user->dt*user->Characteristic.Time/user->Characteristic.SecYear);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [Years]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}
	else
	{
		// Lab time scale or non-dimensional units
		PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ", user->time*user->Characteristic.Time, user->dt*user->Characteristic.Time);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [s]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}

*/

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

/*
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * add bound checking for iterative solution vector in SNES
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolCreate"
PetscErrorCode NLSolCreate(NLSol *nl, PCStokes pc, SNES *p_snes)
{
	SNES            snes;
	KSP             ksp;
	PC              ipc;
	SNESLineSearch  ls;
	JacRes         *jr;
	DOFIndex       *dof;
	PetscBool       flg;

    PetscErrorCode ierr;
    PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(nl, sizeof(NLSol)); CHKERRQ(ierr);

	// store context
 	nl->pc = pc;

 	// access context
	jr  = pc->pm->jr;
	dof = &(jr->fs->cdof);

	// create matrix-free Jacobian operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->J); CHKERRQ(ierr);
	ierr = MatSetUp(nl->J);                              CHKERRQ(ierr);

	// create matrix-free Preconditioner operator
	ierr = MatCreateShell(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, NULL, &nl->P); CHKERRQ(ierr);
	ierr = MatSetUp(nl->P);                              CHKERRQ(ierr);

	// create finite-difference Jacobian
	ierr = MatCreateMFFD(PETSC_COMM_WORLD, dof->numdof, dof->numdof,
		PETSC_DETERMINE, PETSC_DETERMINE, &nl->MFFD); CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(nl->MFFD,"fd_");       CHKERRQ(ierr);
	ierr = MatSetFromOptions(nl->MFFD);               CHKERRQ(ierr);
	ierr = MatSetUp(nl->MFFD);                        CHKERRQ(ierr);

	// setup nonlinear solver
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                     CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESNEWTONLS);                         CHKERRQ(ierr);
	ierr = SNESGetLineSearch(snes, &ls);                            CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);          CHKERRQ(ierr);
	ierr = SNESSetFunction(snes, jr->gres, &FormResidual, nl);      CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes, nl->J, nl->P, &FormJacobian, nl);  CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);                                CHKERRQ(ierr);

	// setup linear solver & preconditioner
	ierr = SNESGetKSP(snes, &ksp);         CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"js_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &ipc);            CHKERRQ(ierr);
	ierr = PCSetType(ipc, PCMAT);          CHKERRQ(ierr);

//	ierr = KSPSetConvergenceTest(ksp, &KSPBlockStopTest, &bmat, NULL);CHKERRQ(ierr);
//	ierr = SNESSetConvergenceTest(snes, SNESBlockStopTest, &nlctx, NULL); CHKERRQ(ierr);

	// set Jacobian type & initial guess
	nl->jtype = _PICARD_;
	ierr = VecSet(jr->gsol, 0.0); CHKERRQ(ierr);

	// read number of Picard iterations
	ierr = PetscOptionsGetInt(PETSC_NULL, "-snes_npicard", &nl->nPicIt, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) nl->nPicIt = 5;

	// return solver
	(*p_snes) = snes;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "NLSolDestroy"
PetscErrorCode NLSolDestroy(NLSol *nl)
{
	PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = MatDestroy(&nl->J);    CHKERRQ(ierr);
	ierr = MatDestroy(&nl->P);    CHKERRQ(ierr);
	ierr = MatDestroy(&nl->MFFD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormResidual"
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx)
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear unused parameters
	if(snes) snes = NULL;

	// access context
	nl = (NLSol*)ctx;
	jr = nl->pc->pm->jr;

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute residual
	ierr = JacResGetResidual(jr); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = JacResCopyRes(jr, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx)
{
	// Compute FDSTAG Jacobian matrix and preconditioner

	NLSol    *nl;
	PCStokes pc;
	PMat     pm;
	JacRes   *jr;
	PetscInt it;

	// clear unused parameters
	if(Amat) Amat = NULL;
	if(Pmat) Pmat = NULL;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	nl = (NLSol*)ctx;
	pc = nl->pc;
	pm = pc->pm;
	jr = pm->jr;

	// setup preconditioner
	ierr = PMatAssemble(pm);                                                  CHKERRQ(ierr);
	ierr = PCStokesSetup(pc);                                                 CHKERRQ(ierr);
	ierr = MatShellSetOperation(nl->P, MATOP_MULT, (void(*)(void))pc->Apply); CHKERRQ(ierr);
	ierr = MatShellSetContext(nl->P, pc);                                     CHKERRQ(ierr);

	// switch Jacobian after fixed number of iterations
	ierr = SNESGetIterationNumber(snes, &it); CHKERRQ(ierr);
	if(it == nl->nPicIt) nl->jtype = _MFFD_;

	// setup Jacobian ...
	if(nl->jtype == _PICARD_)
	{
		// ... Picard
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))pm->Picard); CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, pm->data);                                CHKERRQ(ierr);

	}
	else if(nl->jtype == _MFFD_)
	{
		// ... matrix-free finite-difference (MMFD)
		ierr = MatMFFDSetFunction(nl->MFFD, (PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction, snes); CHKERRQ(ierr);
		ierr = MatMFFDSetBase(nl->MFFD, x, jr->gres);                                                      CHKERRQ(ierr);
		ierr = MatShellSetOperation(nl->J, MATOP_MULT, (void(*)(void))JacApplyMFFD);                        CHKERRQ(ierr);
		ierr = MatShellSetContext(nl->J, (void*)&nl->MFFD);                                                 CHKERRQ(ierr);
	}

	// assemble Jacobian & preconditioner
	ierr = MatAssemblyBegin(nl->P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (nl->P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(nl->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (nl->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacApplyMFFD"
PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y)
{
	Mat *FD;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(A, (void**)&FD); CHKERRQ(ierr);

	// compute Jacobian times vector product
	ierr = MatMult((*FD), x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SNESPrintConvergedReason"
PetscErrorCode SNESPrintConvergedReason(SNES snes)
{
	SNESConvergedReason reason;
	PetscInt            its;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = SNESGetIterationNumber(snes, &its);    CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(snes, &reason);  CHKERRQ(ierr);

    // CONVERGENCE

    if(reason == SNES_CONVERGED_FNORM_ABS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: ||F|| < atol \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_FNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: ||F|| < rtol*||F_initial|| \n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_SNORM_RELATIVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: Newton computed step size small; || delta x || < stol || x ||\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_ITS)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_CONVERGED_TR_DELTA)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Convergence Reason: SNES_CONVERGED_TR_DELTA\n"); CHKERRQ(ierr);
	}

    // DIVERGENCE

	else if(reason == SNES_DIVERGED_FUNCTION_DOMAIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the new x location passed the function is not in the domain of F\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FUNCTION_COUNT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: too many function evaluations\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINEAR_SOLVE)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the linear solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_FNORM_NAN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: residual norm is NAN\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_MAX_IT)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: maximum iterations reached\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LINE_SEARCH)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the line search failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_INNER)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: the inner solve failed\n"); CHKERRQ(ierr);
	}
	else if(reason == SNES_DIVERGED_LOCAL_MIN)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, " SNES Divergence Reason: || J^T b || is small, implies converged to local minimum of F\n"); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of iterations : %lld\n", (LLD)its);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/

