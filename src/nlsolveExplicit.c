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
#include "constEq.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormMomentumResidualAndTheta"
PetscErrorCode FormMomentumResidualAndTheta(Vec x, void *ctx)
{
	NLSol  *nl;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	nl = (NLSol*)ctx;
	jr = nl->pc->pm->jr;

	// copy solution from global to local vectors, enforce boundary constraints
	ierr = JacResCopySol(jr, x, _APPLY_SPC_); CHKERRQ(ierr);

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute momentum residual
	ierr = JacResGetMomentumResidualAndTheta(jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

//---------------------------------------------------------------------------
/*#undef __FUNCT__
#define __FUNCT__ "FormMomentumResidualAndTheta"
PetscErrorCode FormMomentumResidualAndTheta(SNES snes, Vec x, Vec K, void *ctx)
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

    // x=jr->gsol; f=jr->gres;

	// copy solution from global to local vectors, enforce boundary constraints

	ierr = JacResCopySol(jr, x, _APPLY_SPC_); CHKERRQ(ierr);

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute momentum residual and theta
	//ierr = JacResGetMomentumResidualAndTheta(jr); CHKERRQ(ierr);
	ierr = JacResGetMomentumResidual(jr); CHKERRQ(ierr);

	// copy residuals and theta to global vector
	//ierr = JacResCopyMomentumRes(jr, f); CHKERRQ(ierr);
	//ierr = JacResCopyRes(jr, f); CHKERRQ(ierr);
	//ierr = JacResCopyTheta(jr, gK); CHKERRQ(ierr);


	//ierr = JacResCopyK(jr, K); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}*/
//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetVelocities"
PetscErrorCode GetVelocities(JacRes *jr)
{

	//  ... comment 

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscScalar dt, rho_side, t2;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***rho;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access residual context variables
	dt        =  jr->ts.dt;     // time step

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &rho); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)
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
		rho[k][j][i]=svBulk->rho;   // effective density
		//PetscPrintf(PETSC_COMM_WORLD, "    rho  = %12.12e \n", rho);
	}
	END_STD_LOOP

	// communicate boundary values ??
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
		rho_side=(rho[k][j][i]+rho[k][j][i-1])/2;
		//PetscPrintf(PETSC_COMM_WORLD, "    velocityx  = %12.12e \n", vx[k][j][i]);
		vx[k][j][i] += fx[k][j][i]*dt/rho_side;
		/*if (fx[k][j][i] != 0) {
			PetscPrintf(PETSC_COMM_WORLD, "    k  = %i \n", k);
			PetscPrintf(PETSC_COMM_WORLD, "    j  = %i \n", j);
			PetscPrintf(PETSC_COMM_WORLD, "    i  = %i \n", i);
			PetscPrintf(PETSC_COMM_WORLD, "    fx  = %12.12e \n", fx[k][j][i]);
		}*/
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	
	START_STD_LOOP
	{
		rho_side=(rho[k][j][i]+rho[k][j-1][i])/2;
		vy[k][j][i] += fy[k][j][i]*dt/rho_side;
		/*if (fy[k][j][i] != 0) {
			PetscPrintf(PETSC_COMM_WORLD, "    k  = %i \n", k);
			PetscPrintf(PETSC_COMM_WORLD, "    j  = %i \n", j);
			PetscPrintf(PETSC_COMM_WORLD, "    i  = %i \n", i);
			PetscPrintf(PETSC_COMM_WORLD, "    fy  = %12.12e \n", fy[k][j][i]);
		}*/
	}
	END_STD_LOOP
	
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)
	
	// Check vp  ///to be removed  ////////////////////////////////////////////////////////
	if (vz[nz/2][ny/2][10+nx/2]!=0) {
		t2	  =  JacResGetTime(jr);
	}
	//then do W/(nel_x-1)*11/t2 to get aprox. vp
	////////////////////////////////////////////////////////////////////////////////////////

	START_STD_LOOP
	{
		rho_side=(rho[k][j][i]+rho[k-1][j][i])/2;
		vz[k][j][i] += fz[k][j][i]*dt/rho_side;

		/*if (fz[k][j][i] != 0) {
			PetscPrintf(PETSC_COMM_WORLD, "    k  = %i \n", k);
			PetscPrintf(PETSC_COMM_WORLD, "    j  = %i \n", j);
			PetscPrintf(PETSC_COMM_WORLD, "    i  = %i \n", i);
			PetscPrintf(PETSC_COMM_WORLD, "    fz  = %12.12e \n", fz[k][j][i]);
		}*/
	}
	END_STD_LOOP


	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &rho); CHKERRQ(ierr);
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




	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetPressure"
PetscErrorCode GetPressure(JacRes *jr)
{
	//  ... comment

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar IKdt, theta;
	PetscScalar ***up,  ***p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access work vectors

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr); //here is current pressure
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp,   &up);  CHKERRQ(ierr); //here is theta (dvx/dx + dvy/dy + dvz/dz)


	//-------------------------------
	// get pressure from central points
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
		IKdt  = svBulk->IKdt;  // inverse bulk viscosity
		//why no directly svBulk->theta?
		p[k][j][i] = p[k][j][i]-up[k][j][i]/IKdt;
	}
	END_STD_LOOP


	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN,   jr->lp,  &p);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN,   jr->gp,  &up);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_CEN, jr->lp, jr->gp)


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CheckTimeStep"
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user)
{
	// check time step as in Virieux, 1985

	PetscInt    i, numPhases;
	Material_t  *phases, *M;
	PetscScalar dx, dy, dz, dt, rho, shear, bulk, vp, stability;

	PetscFunctionBegin;

	if (user->ExplicitSolver == PETSC_TRUE)	{

		numPhases = jr->numPhases;
		phases    = jr->phases;
		dt        = user->dt;     // time step
		dx=(user->W)/(user->nel_x);
		dy=(user->L)/(user->nel_y);
		dz=(user->H)/(user->nel_z);

		if (user->ExplicitSolver == PETSC_TRUE)	{
			//  phases
			for(i = 0; i < jr->numPhases; i++)
			{
				M = &phases[i];
				rho=M->rho;
				shear=M->G;
				bulk=M->K;
				// P-wave velocity
				vp=sqrt((bulk+4.0/3.0*shear)/rho);
				stability = vp*dt*sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz));
				if ( stability >= 1) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Stability condition = %12.12e", stability);
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
