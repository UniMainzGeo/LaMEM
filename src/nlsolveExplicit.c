//---------------------------------------------------------------------------
//.....................   NONLINEAR SOLVER ROUTINES   .......................
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

//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "NLSolverExp"
PetscErrorCode NLSolverExp(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

//	ierr = FormMomentumResidual(&jr); CHKERRQ(ierr);
//	ierr = GetVelocities(&jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormMomentumResidual"
//PetscErrorCode FormMomentumResidual(JacRes *jr)
PetscErrorCode FormMomentumResidual(SNES snes, Vec x, Vec f, void *ctx)
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

    /*// x=jr->gsol;
    // f=jr->gres;*/

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x, _APPLY_SPC_); CHKERRQ(ierr);

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute residual
	ierr = JacResGetMomentumResidual(jr); CHKERRQ(ierr);

	// copy residuals to global vector
	ierr = JacResCopyMomentumRes(jr, f); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetVelocities"
PetscErrorCode GetVelocities(JacRes *jr, Vec x, Vec f)
{
	//  ... comment 

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;

	PetscScalar rho, dt;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***dxx;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// initialize maximum node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// access residual context variables
	dt        =  jr->ts.dt;     // time step

	// strain-rate component (also used as buffer vectors)
	Vec ldxx; // local (ghosted)

	// access work vectors

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_X,   jr->gfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->gfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->gfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);

	//-------------------------------
	// get density from central points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svBulk = &svCell->svBulk;
		rho   = svBulk->rho;   // effective density
		dxx[k][j][i]=rho;

	}
	END_STD_LOOP

	// communicate boundary values
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);


	//-------------------------------
	// side points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		rho=(dxx[k][j][i]-dxx[k][j][i-1])/2;
		vx[k][j][i] += fx[k][j][i]*dt/rho;
	}
	END_STD_LOOP
	
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	
	START_STD_LOOP
	{
		rho=(dxx[k][j][i]-dxx[k][j-1][i])/2;
		vy[k][j][i] += fy[k][j][i]*dt/rho;
	}
	END_STD_LOOP
	
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)
	
	START_STD_LOOP
	{
		rho=(dxx[k][j][i]-dxx[k-1][j][i])/2;
		vz[k][j][i] += fz[k][j][i]*dt/rho;
	}
	END_STD_LOOP


	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);

	//assemble global residuals from local contributions


	GLOBAL_TO_LOCAL(fs->DA_X, jr->gvx, jr->lvx)
	GLOBAL_TO_LOCAL(fs->DA_Y, jr->gvy, jr->lvy)
	GLOBAL_TO_LOCAL(fs->DA_Z, jr->gvz, jr->lvz)

	//LOCAL_TO_GLOBAL()


	// copy velocities to global vector
	//ierr = JacResCopyVel(jr, f, appSPC); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
