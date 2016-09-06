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
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "nlsolveExplicit.h"

//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetVelocities"
PetscErrorCode GetVelocities(JacRes *jr, UserCtx *user)
{

	// Calculates the velocity from the residual and the densities
	//  v_new = v_old - residual*dt/density
	// If there are absorbing boundaries, the velocity is dampened
	// v_new = v_new*damping_factor

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscScalar dt, rho_side;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***rho;

//	PetscScalar t2, max_vel;

	// Scaling computational density factor
	PetscScalar DensityFactor = user->DensityFactor;

	// To damp velocity in the absorbing boundaries
	PetscScalar damping; //, damping_x,damping_y,damping_z;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access residual context variables
	dt        =  jr->ts.dt;     // time step

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &rho); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)
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
		rho[k][j][i]=svBulk->rho*DensityFactor;   // effective density
		//PetscPrintf(PETSC_COMM_WORLD, "    rho  = %12.12e \n", rho[k][j][i]);
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &rho); CHKERRQ(ierr);

	// communicate boundary values ??
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldzz);

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &rho); CHKERRQ(ierr);

	//-------------------------------
	// side points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if (i==0) rho_side = rho[k][j][i];
		else rho_side = (rho[k][j][i]+rho[k][j][i-1])/2.0;

		/*damping_x = GetBoundaryDamping('x',user,i,j,k); // To damp velocity if we are in an absorbing boundary
		damping_y = GetBoundaryDamping('y',user,i,j,k);
		damping_z = GetBoundaryDamping('z',user,i,j,k);
		damping = damping_x*damping_y*damping_z;*/

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
		//PetscPrintf(PETSC_COMM_WORLD, "    damping[%i,%i,%i] = %12.12e \n", i,j,k, damping);
		vx[k][j][i] = (vx[k][j][i]-fx[k][j][i]*dt/rho_side)*damping;

		//PetscPrintf(PETSC_COMM_WORLD, "    vx[%i,%i,%i]  %12.12e \n", i,j,k,vx[k][j][i]);

		if (i==0)
		{
			rho_side = rho[k][j][i];
		}
		else
		{
			rho_side = (rho[k][j][i]+rho[k][j][i-1])/2.0;
		}
		//PetscPrintf(PETSC_COMM_WORLD, "    velocityx  = %12.12e \n", vx[k][j][i]);
		vx[k][j][i] -= fx[k][j][i]*dt/rho_side;

		//PetscPrintf(PETSC_COMM_WORLD, "    [k,j,i]  = [%i,%i,%i]  vx  = %12.12e fx =  %12.12e \n", k,j,i,vx[k][j][i], fx[k][j][i]);
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if (j==0) rho_side=rho[k][j][i];
		else rho_side=(rho[k][j][i]+rho[k][j-1][i])/2.0;

		/*damping_x = GetBoundaryDamping('x',user,i,j,k); // To damp velocity if we are in an absorbing boundary
		damping_y = GetBoundaryDamping('y',user,i,j,k);
		damping_z = GetBoundaryDamping('z',user,i,j,k);
		damping = damping_x*damping_y*damping_z;*/

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
		vy[k][j][i] = (vy[k][j][i]-fy[k][j][i]*dt/rho_side)*damping;

		//PetscPrintf(PETSC_COMM_WORLD, "    vy[%i,%i,%i]  %12.12e \n", i,j,k,vy[k][j][i]);
		if (j==0)
		{
			rho_side=rho[k][j][i];
		}
		else
		{
			rho_side=(rho[k][j][i]+rho[k][j-1][i])/2.0;
		}
		vy[k][j][i] -= fy[k][j][i]*dt/rho_side;
		//PetscPrintf(PETSC_COMM_WORLD, "    [k,j,i]  = [%i,%i,%i]  vy  = %12.12e fy =  %12.12e \n", k,j,i,vy[k][j][i], fy[k][j][i]);
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)


//max_vel=0.0;				////////////////////////////////////////////////////////////////////////////////////////

	START_STD_LOOP
	{
		if (k==0) rho_side=rho[k][j][i];
		else rho_side=(rho[k][j][i]+rho[k-1][j][i])/2.0;

		/*damping_x = GetBoundaryDamping('x',user,i,j,k); // To damp velocity if we are in an absorbing boundary
		damping_y = GetBoundaryDamping('y',user,i,j,k);
		damping_z = GetBoundaryDamping('z',user,i,j,k);
		damping = damping_x*damping_y*damping_z;*/

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
		vz[k][j][i] = (vz[k][j][i]-fz[k][j][i]*dt/rho_side)*damping;

		//PetscPrintf(PETSC_COMM_WORLD, "    vz[%i,%i,%i]  %12.12e \n", i,j,k,vz[k][j][i]);
		if (k==0)
		{
			rho_side=rho[k][j][i];
		}
		else
		{
			rho_side=(rho[k][j][i]+rho[k-1][j][i])/2.0;
		}
		vz[k][j][i] -= fz[k][j][i]*dt/rho_side;
		//PetscPrintf(PETSC_COMM_WORLD, "    [k,j,i]  = [%i,%i,%i]  vz  = %12.12e fz =  %12.12e \n", k,j,i,vz[k][j][i], fz[k][j][i]);


//if (vz[k][j][i] > max_vel) {////////////////////////////////////////////////////////////////////////////////////////
//	max_vel=vz[k][j][i];	////////////////////////////////////////////////////////////////////////////////////////
//}							////////////////////////////////////////////////////////////////////////////////////////

	}
	END_STD_LOOP


//PetscPrintf(PETSC_COMM_WORLD, "    max_vel = %12.12e \n", max_vel);////////////////////////////////////////////////////////////////////////////////////////

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &rho); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);


	//assemble global residuals from local contributions


	//GLOBAL_TO_LOCAL(fs->DA_X, jr->gvx, jr->lvx)
	//GLOBAL_TO_LOCAL(fs->DA_Y, jr->gvy, jr->lvy)
	//GLOBAL_TO_LOCAL(fs->DA_Z, jr->gvz, jr->lvz)

	//LOCAL_TO_GLOBAL()

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormMomentumResidualPressureAndVelocities"
PetscErrorCode FormMomentumResidualPressureAndVelocities(JacRes *jr, UserCtx *user)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// copy solution from global to local vectors, enforce boundary constraints
	//ierr = JacResCopySol(jr, x, _APPLY_SPC_); CHKERRQ(ierr);

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute momentum residual and pressure
	ierr = JacResGetMomentumResidualAndPressure(jr,user); CHKERRQ(ierr);

	ierr = GetVelocities(jr, user);	CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CheckElasticProperties"
PetscErrorCode CheckElasticProperties(JacRes *jr, UserCtx *user)
{
	// Check that elastic properties remain constant

	FDSTAG     *fs;
	SolVarCell *svCell;
//	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;

	PetscInt i, j, k, iter/*, numPhases*/;
	PetscInt nx, ny, nz, sx, sy, sz;
	Material_t  *phases, *M;
	PetscScalar /*dx, dy, dz,*/ dt, rho, shear, bulk;

	PetscFunctionBegin;

	fs = jr->fs;

//	numPhases = jr->numPhases;
	phases    = jr->phases;

	dt        = user->dt;     // time step

	if (user->ExplicitSolver == PETSC_TRUE)	{

		//-------------------------------
		// central points
		//-------------------------------
		iter = 0;
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{

			// access solution variables
			svCell = &jr->svCell[iter++];
			svDev  = &svCell->svDev;
			svBulk = &svCell->svBulk;

			rho   = svBulk->rho;
			bulk  = 1.0/svBulk->IKdt/dt;
			shear = 1.0/svDev->I2Gdt/2.0/dt;

			//PetscPrintf(PETSC_COMM_WORLD, "    rho, bulk, shear  = %12.12e, %12.12e, %12.12e \n", rho, bulk, shear);

			//  phases
			//for(i = 0; i < jr->numPhases; i++)
			//{
				M = &phases[0];
				if (M->rho != rho) 		PetscPrintf(PETSC_COMM_WORLD, "     ---------------------------- rho in file, rho in the cell  = %12.12e, %12.12e \n", M->rho, rho);
				if (M->G   != shear) 	PetscPrintf(PETSC_COMM_WORLD, "     ---------------------------- shear in file, shear in the cell  = %12.12e, %12.12e \n", M->G, shear);
				if (M->K   != bulk) 	PetscPrintf(PETSC_COMM_WORLD, "     ---------------------------- bulk in file, bulk in the cell  = %12.12e, %12.12e \n", M->K, bulk);
			//}
		}
		END_STD_LOOP
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CheckTimeStep"
PetscErrorCode CheckTimeStep(JacRes *jr, UserCtx *user)
{
	// check time step as in Virieux, 1985

	PetscInt    i /*,numPhases*/;
	Material_t  *phases, *M;
	PetscScalar dx, dy, dz, dt, rho, shear, bulk, vp, stability;

	PetscFunctionBegin;


//	numPhases = jr->numPhases;
	phases    = jr->phases;
	dt        = user->dt;     // time step
	dx = user->W/((PetscScalar)(user->nel_x));
	dy = user->L/((PetscScalar)(user->nel_y));
	dz = user->H/((PetscScalar)(user->nel_z));

	//if (user->ExplicitSolver == PETSC_TRUE)	{
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
			if ( stability > 1) {
				SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Stability condition = %12.12e", stability);
			}
		}
	//}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "UpdateHistoryFieldsAndGetAxialStressStrain"
PetscErrorCode UpdateHistoryFieldsAndGetAxialStressStrain(JacRes *jr, PetscScalar *axial_stress, PetscScalar *axial_strain)
{

	// Update svCell->hxx, yy, zz and svBulk->pn

	// Calculates the pressure from velocities and bulk
	//  p_new = p_old - (dv/dx + dv/dy + dv/dz)*K*dt



	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	SolVarEdge *svEdge;
	PetscScalar ***p;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscInt    iter;

	PetscScalar IKdt, theta;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);


	PetscScalar stress_xx, stress_yy, stress_zz;
	PetscScalar strain_xx, strain_yy, strain_zz;
	PetscScalar count;

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr); //here is current pressure
	//ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp,   &up);  CHKERRQ(ierr); //here is theta (dvx/dx + dvy/dy + dvz/dz)


	stress_xx = 0;
	stress_yy = 0;
	stress_zz = 0;
	strain_xx = 0;
	strain_yy = 0;
	strain_zz = 0;
	count = 0;

	//-------------------------------
	// central points
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

		svCell->hxx = svCell->sxx;
		svCell->hyy = svCell->syy;
		svCell->hzz = svCell->szz;

		svBulk->pn = p[k][j][i];

		// Calculate axial stress
		stress_xx = stress_xx + svCell->sxx - p[k][j][i];
		stress_yy = stress_yy + svCell->syy - p[k][j][i];
		stress_zz = stress_zz + svCell->szz - p[k][j][i];

		strain_xx = strain_xx + svCell->dxx;
		strain_yy = strain_yy + svCell->dyy;
		strain_zz = strain_zz + svCell->dzz;
		count = count + 1;

	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svEdge->h = svEdge->s;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svEdge->h = svEdge->s;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svEdge->h = svEdge->s;
	}
	END_STD_LOOP

	stress_xx = stress_xx/count;
	stress_yy = stress_yy/count;
	stress_zz = stress_zz/count;

	strain_xx = strain_xx/count;
	strain_yy = strain_yy/count;
	strain_zz = strain_zz/count;

	// COMPLETE COMMUNICATION HERE!

	//MPI_Reduce

	*axial_stress = sqrt(stress_xx*stress_xx + stress_yy*stress_yy + stress_zz*stress_zz);
	*axial_strain = sqrt(strain_xx*strain_xx + strain_yy*strain_yy + strain_zz*strain_zz);



	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ChangeTimeStep"
PetscErrorCode ChangeTimeStep(JacRes *jr, UserCtx *user)
{
	// from Virieux, 1985

	PetscInt    i, numPhases;
	Material_t  *phases, *M;
	PetscScalar dx, dy, dz, empty, dt, rho, shear, bulk, vp, stability, d_average, computational_density_factor;

	PetscFunctionBegin;

	computational_density_factor = user->DensityFactor;
	numPhases = jr->numPhases;
	phases    = jr->phases;
	empty  = 999999999;
	dt = empty;
	dx = user->W/((PetscScalar)(user->nel_x));
	dy = user->L/((PetscScalar)(user->nel_y));
	dz = user->H/((PetscScalar)(user->nel_z));

	//  phases
	for(i = 0; i < numPhases; i++)
	{
		M = &phases[i];
		rho=M->rho;
		shear=M->G;
		bulk=M->K;

		rho=rho*computational_density_factor;

		vp=sqrt((bulk+4.0/3.0*shear)/rho);				// P-wave velocity
		//d_average 	= sqrt(dx*dx + dz*dz + dy*dy);   	// average spacing
		d_average 	= sqrt(1/(dx*dx) + 1/(dz*dz) + 1/(dy*dy));

		if (dt == empty) {
			//if (d_average/vp < dt) {
			//dt = 1/vp/sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz));
				//dt = d_average/vp;
			dt = 1/vp/d_average;
			//}
		}
		else {
			//if (d_average/vp < dt) {
			if (1/vp/d_average < dt) {
				dt = d_average/vp;
			}
		}
	}
	user->dt = dt;
	jr->ts.dt = dt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------





//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "ShowValues"
PetscErrorCode ShowValues(JacRes *jr, UserCtx *user, PetscInt n)
{
	// Show the values of velocity, pressure, stress, strain,...
/*
	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;

	PetscInt i, j, k, iter, numPhases;
	PetscInt nx, ny, nz, sx, sy, sz;
	Material_t  *phases, *M;
	PetscScalar dx, dy, dz, dt, rho, shear, bulk, time;

	PetscScalar ***fx,  ***fy,  ***fz, ***gfx,  ***gfy,  ***gfz, ***vx,  ***vy,  ***vz, ***lvx,  ***lvy,  ***lvz;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxz, ***p, ***up;

	PetscErrorCode ierr;
	PetscFunctionBegin;

//	FILE *fseism;
//	fseism = user->Station.output_file;


	FILE *fseism;
	fseism = user->Station.output_file;


	fs = jr->fs;


	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp,   &up);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->gfx,  &gfx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->gfy,  &gfy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->gfz,  &gfz);  CHKERRQ(ierr);


	ierr = DMDAVecGetArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &lvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &lvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &lvz);  CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	//PetscPrintf(PETSC_COMM_WORLD, "    [%i,%i,%i] \n", nx,ny,nz);

	//PetscPrintf(PETSC_COMM_WORLD, "\n    %i ************************************** \n ",n);
	//fprintf(fseism, "%i ************************************** \n ",n);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		rho   = svBulk->rho;
		bulk  = 1.0/svBulk->IKdt/dt;
		shear = 1.0/svDev->I2Gdt/2.0/dt;

		//PetscPrintf(PETSC_COMM_WORLD, "    rho[%i,%i,%i]  = %12.12e \n", i,j,k, rho);

		if (k==1 ) //&& j==2 && i==2)
		{
			//fprintf(fseism, "%i %i %i\n", i,j,k);

			//PetscPrintf(PETSC_COMM_WORLD, "\n    %i - Values in [%i,%i,%i] ************************************** \n ",n, i,j,k);
			//PetscPrintf(PETSC_COMM_WORLD, "    rho, bulk, shear  = %12.12e, %12.12e, %12.12e \n", rho, bulk, shear);
			//PetscPrintf(PETSC_COMM_WORLD, "    p[%i,%i,%i]  = %12.12e \n", i,j,k, p[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    up[%i,%i,%i]  = %12.12e \n", i,j,k, up[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    lfx,y,z[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, fx[k][j][i],fy[k][j][i],fz[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    gfx,y,z[%i,%i,%i] = (%12.12e,%12.12e,%12.12e) \n", i,j,k, gfx[k][j][i],gfy[k][j][i],gfz[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    ldxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, dxx[k][j][i],dyy[k][j][i],dzz[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    svCell.hxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->hxx,svCell->hyy,svCell->hzz);
			//PetscPrintf(PETSC_COMM_WORLD, "    svCell.sxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->sxx,svCell->syy,svCell->szz);
			//PetscPrintf(PETSC_COMM_WORLD, "    svCell.dxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->dxx,svCell->dyy,svCell->dzz);
			//PetscPrintf(PETSC_COMM_WORLD, "    theta[%i,%i,%i]  = %12.12e \n", i,j,k, svBulk->theta);

		}
		//PetscPrintf(PETSC_COMM_WORLD, "    gfx,y,z[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, gfx[k][j][i],gfy[k][j][i],gfz[k][j][i]);
		//PetscPrintf(PETSC_COMM_WORLD, "    svBulk.theta[%i,%i,%i]  = %12.12e \n", i,j,k, svBulk->theta);
		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.sxx,syy,szz[%i,%i,%i]  = %12.12e, %12.12e, %12.12e \n", i,j,k, svCell->sxx, svCell->syy, svCell->szz);
		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.dxx[%i,%i,%i]  = %12.12e \n", i,j,k, svCell->dxx);
		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.sxx[%i,%i,%i]  = (%12.12e) \n", i,j,k,svCell->sxx);
<<<<<<< 45303369d0e2a4894286aa8ad0f8feefc46161a6

		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.hxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->hxx,svCell->hyy,svCell->hzz);
		//fprintf(fseism, "%i %i %i\n", i,j,k);
		//fprintf(fseism, "%12.12e\n", gfz[k][j][i]);
		//if ( p[k][j][i] != up[k][j][i])			fprintf(fseism, "%12.12e\n", p[k][j][i]);
		//if ( fx[k][j][i] != gfx[k][j][i])		 fprintf(fseism, "%12.12e\n", fx[k][j][i]);
		//PetscPrintf(PETSC_COMM_WORLD, "    p[%i,%i,%i]  = %12.12e \n", i,j,k, p[k][j][i]);
		//PetscPrintf(PETSC_COMM_WORLD, "    rho  = %12.12e \n", svBulk->rho);

=======
		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.hxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->hxx,svCell->hyy,svCell->hzz);
>>>>>>> First changes for density scaling
		//fprintf(fseism, "%i %i %i\n", i,j,k);
		//fprintf(fseism, "%12.12e\n", gfz[k][j][i]);
		//if ( p[k][j][i] != up[k][j][i])			fprintf(fseism, "%12.12e\n", p[k][j][i]);
		//if ( fx[k][j][i] != gfx[k][j][i])			fprintf(fseism, "%12.12e\n", fx[k][j][i]);

		//PetscPrintf(PETSC_COMM_WORLD, "    rho  = %12.12e \n", svBulk->rho);

		//fprintf(fseism, "%12.12e\n", svCell->szz);
	}
	END_STD_LOOP

	//-------------------------------
	// side points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		//if (k==0 && j==2 && i==2)
		//{
			//PetscPrintf(PETSC_COMM_WORLD, "    gvx, lvx [%i,%i,%i] = %12.12e, %12.12e\n", i,j,k, vx[k][j][i], lvx[k][j][i]);
			//fprintf(fseism, "    gvx, lvx [%i,%i,%i] = %12.12e, %12.12e\n", i,j,k, vx[k][j][i], lvx[k][j][i]);
		//if ( vx[k][j][i] != lvx[k][j][i]) fprintf(fseism, "%12.12e\n", vx[k][j][i]);
		//}
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		//if (k==0 && j==2 && i==2)
		//{
			//PetscPrintf(PETSC_COMM_WORLD, "    gvy, lvy [%i,%i,%i]= %12.12e, %12.12e\n", i,j,k, vy[k][j][i], lvy[k][j][i]);
		//fprintf(fseism, "%12.12e\n", lvy[k][j][i]);
		//}
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		//if (k==0 ) //&& j==2 && i==2)
		//{
			//PetscPrintf(PETSC_COMM_WORLD, "\n    gfz, lfz[%i,%i,%i]  = %12.12e, %12.12e) \n", i,j,k, gfz[k][j][i], fz[k][j][i]);
			//PetscPrintf(PETSC_COMM_WORLD, "    gvz, lvz [%i,%i,%i] = %12.12e, %12.12e\n", i,j,k, vz[k][j][i], lvz[k][j][i]);
		//fprintf(fseism, "    gvz, lvz [%i,%i,%i] = %12.12e, %12.12e\n", i,j,k, vz[k][j][i], lvz[k][j][i]);
		//fprintf(fseism, "    %12.12e\n", lvz[k][j][i]);
		//fprintf(fseism, "    [%i,%i,%i] %12.12e\n", i,j,k, vz[k][j][i]);
		//fprintf(fseism, "%12.12e\n", vz[k][j][i]);
		//}
	}
	END_STD_LOOP

	//PetscPrintf(PETSC_COMM_WORLD, "    **************************************************************************** \n \n");

	// restore vectors

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gp,   &up);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gfx,  &gfx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gfy,  &gfy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gfz,  &gfz);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &lvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &lvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &lvz);  CHKERRQ(ierr);

*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Apply source
#undef __FUNCT__
#define __FUNCT__ "GetCellCoordinatesSource"
PetscErrorCode GetCellCoordinatesSource(JacRes *jr)
{
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar iter;
	PetscScalar    xs[3], xe[3];
	FDSTAG     *fs;

	fs = jr->fs;

	if (jr->SeismicSource==PETSC_TRUE && (jr->SourceParams.source_type==POINT || jr->SourceParams.source_type==MOMENT) ) {
		//-------------------------------
		// central points
		//-------------------------------
		iter = 0;
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// get cell coordinates
			xs[0] = jr->fs->dsx.ncoor[i]; xe[0] = jr->fs->dsx.ncoor[i+1];
			xs[1] = jr->fs->dsy.ncoor[j]; xe[1] = jr->fs->dsy.ncoor[j+1];
			xs[2] = jr->fs->dsz.ncoor[k]; xe[2] = jr->fs->dsz.ncoor[k+1];

			if (jr->SourceParams.x > xs[0] && jr->SourceParams.x <= xe[0] && jr->SourceParams.y > xs[1] && jr->SourceParams.y <= xe[1] && jr->SourceParams.z > xs[2] && jr->SourceParams.z <= xe[2])
			{
				jr->SourceParams.i=i;
				jr->SourceParams.j=j;
				jr->SourceParams.k=k;
			}
		}
		END_STD_LOOP
	}
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetStressFromSource"
PetscErrorCode GetStressFromSource(JacRes *jr, UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz)
{
	PetscInt M, N, P;
	PetscScalar /*coor,*/ time;
	PetscScalar    xs[3], xe[3];

	time	  =  JacResGetTime(jr);

	//FILE *fseism;
	//fseism = user->Station.output_file;

	if (jr->SourceParams.source_type == PLANE)
	{

		//if (k==user->nel_z-1) //(k==1)
		if (i==user->nel_x-1)
		{
			//*sxx = 		jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0))); 	//jr->SourceParams.amplitude;
			//*syy =	- 	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2; 	//jr->SourceParams.amplitude/2.0;
			//*szz = 	-  	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;	//jr->SourceParams.amplitude/2.0;

			//*szz = 		jr->SourceParams.amplitude;
			//*sxx = 0.0;
			//*syy =	 	-jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
			//*szz = 0.0;

			//*szz=jr->SourceParams.amplitude;
			//*syy=0.0;
			//*sxx=0.0;

			*sxx = 	100.0;
			*syy =	50.0;;
			*szz = 	-50.0 ;

		if (k==user->nel_z-1) //(k==1)
		{
			*szz = 		jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0))); 	//jr->SourceParams.amplitude;
			*sxx =	- 	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2; 	//jr->SourceParams.amplitude/2.0;
			*syy = 	-  	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;	//jr->SourceParams.amplitude/2.0;
		}
	}
	else if (jr->SourceParams.source_type == COMPRES)
		{
			if (k==1)
			{
				*szz = 		jr->SourceParams.amplitude;
				*sxx =	- 	jr->SourceParams.amplitude/2.0;
				*syy = 	-  	jr->SourceParams.amplitude/2.0;
			}
			else if (k==user->nel_z-1)
			{
				*szz = 	-	jr->SourceParams.amplitude;
				*sxx =	- 	jr->SourceParams.amplitude/2.0;
				*syy = 	-  	jr->SourceParams.amplitude/2.0;
			}
		}
	else if (jr->SourceParams.source_type == POINT)
	{
		// get cell coordinates
		xs[0] = jr->fs->dsx.ncoor[i]; xe[0] = jr->fs->dsx.ncoor[i+1];
		xs[1] = jr->fs->dsy.ncoor[j]; xe[1] = jr->fs->dsy.ncoor[j+1];
		xs[2] = jr->fs->dsz.ncoor[k]; xe[2] = jr->fs->dsz.ncoor[k+1];

		if (jr->SourceParams.x > xs[0] && jr->SourceParams.x <= xe[0] && jr->SourceParams.y > xs[1] && jr->SourceParams.y <= xe[1] && jr->SourceParams.z > xs[2] && jr->SourceParams.z <= xe[2])
				{
					*sxx = jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
					*syy = jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;
					*szz = -jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;

					//*sxx = 	100.0;
					//*syy =	50.0;;
					//*szz = 	-50.0 ;
				}



		/*if (k==70 && j == 70 && i == 20)
		{
<<<<<<< 5c262a0cadb3ec043a0807d3e91d5f400eeed756
			*sxx = jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
			*szz = - jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
=======
			*sxx = *sxx*0 + jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
			// *syy = 0.0;
			*szz = *szz*0 + jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
>>>>>>> reduced number of warnings


			*sxx = 	100.0;
			*syy =	50.0;;
			*szz = 	-50.0 ;
		}*/
	}
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

// Get damping factor for absorbing boundary
#undef __FUNCT__
#define __FUNCT__ "GetBoundaryDamping"
PetscScalar GetBoundaryDamping( UserCtx *user, PetscInt i, PetscInt j, PetscInt k)
{
	PetscScalar Damping, A0;

	Damping = 1.0;

	if (user->AbsBoundaries == PETSC_TRUE)
	{
		A0 = 0.92;	// change by user->...?

		if (k<user->AB.NzL) {
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(k)/(user->AB.NzL-1))));
		}
		else if (k >(user->nel_z -user->AB.NzR)) {
			//Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_z-1 -k)/(user->AB.NzR-1))));
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_z -k)/(user->AB.NzR-1))));
		}
		if (j<user->AB.NyL) {
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(j)/(user->AB.NyL-1))));
		}
		else if (j >(user->nel_y -user->AB.NyR)) {
			//Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_y-1 -j)/(user->AB.NyR-1))));
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_y -j)/(user->AB.NyR-1))));
		}

		if (i<user->AB.NxL) {
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(i)/(user->AB.NxL-1))));
		}
		else if (i >(user->nel_x -user->AB.NxR)) {
			//Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_x-1 -i)/(user->AB.NxR-1))));
			Damping=Damping*(A0+(1/2)*(1-A0)*(1-cos(M_PI*(user->nel_x -i)/(user->AB.NxR-1))));
		}
	}
	PetscFunctionReturn(Damping);
}
