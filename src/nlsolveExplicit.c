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

	// Scaling computational density factor
	PetscScalar DensityFactor = user->DensityFactor;

	// To damp velocity in the absorbing boundaries
	PetscScalar damping;

	// To save velocity in a given point (seismic station)
	PetscScalar t;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access residual context variables
	dt    =  jr->ts.dt;     // time step
	t	  =  JacResGetTime(jr);

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

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
	//PetscPrintf(PETSC_COMM_WORLD, "    damping[%i,%i,%i] = %12.12e \n", i,j,k, damping);


		vx[k][j][i] = (vx[k][j][i]-fx[k][j][i]*dt/rho_side)*damping;

		if (jr->SeismicStation == PETSC_TRUE) { // To save velocity in a given point
			if (jr->StationParams.i==i && jr->StationParams.j==j && jr->StationParams.k==k) {
				//fprintf(jr->StationParams.output_file[0], "%12.12e %12.12e\n", t, (vx[k][j][i]+vx[k][j][i-1])/2.0);
				fprintf(jr->StationParams.output_file[0], "%12.12e %12.12e\n", t, vx[k][j][i]);
			}
		}
		//PetscPrintf(PETSC_COMM_WORLD, "    damping[%i,%i,%i] = %12.12e \n", i,j,k, damping);
		//PetscPrintf(PETSC_COMM_WORLD, "    vx[%i,%i,%i]  %12.12e \n", i,j,k,vx[k][j][i]);

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

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
		vy[k][j][i] = (vy[k][j][i]-fy[k][j][i]*dt/rho_side)*damping;


		if (jr->SeismicStation == PETSC_TRUE) { // To save velocity in a given point
			if (jr->StationParams.i==i && jr->StationParams.j==j && jr->StationParams.k==k) {
				//fprintf(jr->StationParams.output_file[1], "%12.12e %12.12e\n", t, (vy[k][j][i]+vy[k][j-1][i])/2.0);
				fprintf(jr->StationParams.output_file[1], "%12.12e %12.12e\n", t, vy[k][j][i]);
			}
		}
	}
	END_STD_LOOP

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if (k==0) rho_side=rho[k][j][i];
		else rho_side=(rho[k][j][i]+rho[k-1][j][i])/2.0;

		damping = GetBoundaryDamping(user,i,j,k); // To damp velocity if we are in an absorbing boundary
		vz[k][j][i] = (vz[k][j][i]-fz[k][j][i]*dt/rho_side)*damping;


		if (jr->SeismicStation == PETSC_TRUE) { // To save velocity in a given point
			if (jr->StationParams.i==i && jr->StationParams.j==j && jr->StationParams.k==k) {
				//fprintf(jr->StationParams.output_file[2], "%12.12e %12.12e\n", t, (vz[k][j][i]+vy[k-1][j][i])/2.0);
				fprintf(jr->StationParams.output_file[2], "%12.12e %12.12e\n", t, vz[k][j][i]);
			}
		}
	}
	END_STD_LOOP

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

	// apply pressure limit at the first visco-plastic timestep and iteration
	if(jr->ts.istep == 1 && jr->matLim.presLimAct == PETSC_TRUE)
	{
		jr->matLim.presLimFlg = PETSC_TRUE;
	}

	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);

	// compute momentum residual and pressure
	ierr = JacResGetMomentumResidualAndPressure(jr,user); CHKERRQ(ierr);

	ierr = GetVelocities(jr, user);	CHKERRQ(ierr);

	// deactivate pressure limit after it has been activated
	jr->matLim.presLimFlg = PETSC_FALSE;

	PetscFunctionReturn(0);

}
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "UpdateHistoryFieldsAndGetAxialStressStrain"
PetscErrorCode UpdateHistoryFieldsAndGetAxialStressStrain(JacRes *jr, PetscScalar *axial_stress)
{
	// Update svCell->hxx, yy, zz and svBulk->pn, and calculate second invariant of stresses DII

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	SolVarEdge *svEdge;
	PetscScalar ***p, ***JII_center, ***JII_xy, ***JII_xz, ***JII_yz, ***JII_corner;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,    &p);   		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx,  &JII_center); CHKERRQ(ierr);		// strain-rate component (used as buffer vectors)
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy,  &JII_xy); 	CHKERRQ(ierr);		// strain-rate component (used as buffer vectors)
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz,  &JII_xz); 	CHKERRQ(ierr);		// strain-rate component (used as buffer vectors)
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz,  &JII_yz); 	CHKERRQ(ierr);		// strain-rate component (used as buffer vectors)
	ierr = DMDAVecGetArray(fs->DA_COR, jr->lbcor, &JII_corner); CHKERRQ(ierr);

	PetscScalar stress_II, sum_stress_II, count, sum_count;

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

		// update historic fields
		svCell->hxx = svCell->sxx;
		svCell->hyy = svCell->syy;
		svCell->hzz = svCell->szz;
		svBulk->pn = p[k][j][i];

		// Finding DII
		JII_center[k][j][i] = 0.5 * (svCell->sxx*svCell->sxx + svCell->syy*svCell->syy + svCell->szz*svCell->szz);

	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &JII_center);  CHKERRQ(ierr);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	//ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &JII_center); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)

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

		// update historic field
		svEdge->h = svEdge->s;

		// Finding DII
		JII_xy[k][j][i] = svEdge->s*svEdge->s;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &JII_xy);  CHKERRQ(ierr);
	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);
	//ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &JII_xy); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)

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

		// update historic field
		svEdge->h = svEdge->s;

		// Finding DII
		JII_xz[k][j][i] = svEdge->s*svEdge->s;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &JII_xz);  CHKERRQ(ierr);
	LOCAL_TO_LOCAL(fs->DA_XZ, jr->ldxz);
	//ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &JII_xz); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)

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

		// update historic field
		svEdge->h = svEdge->s;

		// Finding DII
		JII_yz[k][j][i] = svEdge->s*svEdge->s;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &JII_yz);  CHKERRQ(ierr);
	LOCAL_TO_LOCAL(fs->DA_YZ, jr->ldyz);
	//ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &JII_yz); CHKERRQ(ierr);	// strain-rate component (used as buffer vectors)

	// Interpolate stress_II to the corners
	InterpFlags iflag;
	iflag.update = PETSC_TRUE;
	iflag.use_bound = PETSC_FALSE;
	ierr = VecSet(jr->lbcor, 0.0); CHKERRQ(ierr);
	InterpCenterCorner(fs,jr->ldxx,jr->lbcor,iflag);
	InterpXYEdgeCorner(fs,jr->ldxy,jr->lbcor,iflag);
	InterpXZEdgeCorner(fs,jr->ldxz,jr->lbcor,iflag);
	InterpYZEdgeCorner(fs,jr->ldyz,jr->lbcor,iflag);
	// store second invariant
	ierr = VecSqrtAbs(jr->lbcor); CHKERRQ(ierr);

	//-------------------------------
	// corner points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	stress_II = 0;
	sum_stress_II = 0;
	count = 0.0;
	sum_count = 0;

	START_STD_LOOP
	{
		stress_II += JII_corner[k][j][i];
		count += 1.0;
	}
	END_STD_LOOP

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&stress_II, &sum_stress_II, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
		ierr = MPI_Allreduce(&count, &sum_count, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		sum_stress_II = stress_II;
		sum_count = count;
	}

	*axial_stress = sum_stress_II/sum_count;


	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,    &p);   		CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx,  &JII_center); CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy,  &JII_xy); 	CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz,  &JII_xz); 	CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz,  &JII_yz); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, jr->lbcor, &JII_corner); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

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

	//dt        = user->dt;     // time step
	dt = jr -> ts.dt;

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

	PetscInt    i;
	Material_t  *phases, *M;
	PetscScalar dx, dy, dz, dt, rho, shear, bulk, vp, stability, computational_density_factor;

	PetscFunctionBegin;

	computational_density_factor = user->DensityFactor;

	phases    = jr->phases;
	dt        = user->dt;     // time step
	dx 		  = user->W/((PetscScalar)(user->nel_x));
	dy 		  = user->L/((PetscScalar)(user->nel_y));
	dz 		 = user->H/((PetscScalar)(user->nel_z));

	//  phases
	for(i = 0; i < jr->numPhases; i++)
	{
		M = &phases[i];
		rho=M->rho;
		shear=M->G;
		bulk=M->K;

		rho=rho*computational_density_factor;

		vp = sqrt((bulk+4.0/3.0*shear)/rho); 					// P-wave velocity
		stability = vp*dt*sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz));
		if ( stability >= 1) {
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Stability condition = %12.12e", stability);
		}
	}
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
	PetscScalar dx, dy, dz, dt_min, dt, rho, shear, bulk, vp;
	PetscScalar CFL, computational_density_factor;

	PetscFunctionBegin;

	computational_density_factor = user->DensityFactor;
	numPhases 	= 	jr->numPhases;
	phases    	= 	jr->phases;
	dt_min 		= 	PETSC_MAX_REAL;
	dx 			= 	user->W/((PetscScalar)(user->nel_x));
	dy 			= 	user->L/((PetscScalar)(user->nel_y));
	dz 			= 	user->H/((PetscScalar)(user->nel_z));
	CFL 		=	user->CFL;

	//  phases
	for(i = 0; i < numPhases; i++)
	{
		M 			= 	&phases[i];
		rho			=	M->rho;
		shear 		=	M->G;
		bulk 		=	M->K;

		rho 		=	rho*computational_density_factor;

		vp 			=	sqrt( (bulk+4.0/3.0*shear)/rho);				// P-wave velocity



		// Compute velocity such that the wave does not move more than CFL times a gridcell per dt
		dt 			= 	CFL*dx/vp;
		if (dt_min>dt){dt_min =dt;}
		dt 			= 	CFL*dy/vp;
		if (dt_min>dt){dt_min =dt;}
		dt 			= 	CFL*dz/vp;
		if (dt_min>dt){dt_min =dt;}
	}
	//user->dt 	= dt_min;
	jr->ts.dt 	= dt_min;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


/*//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetPressure"
PetscErrorCode GetPressure(JacRes *jr)
{
	// Calculates the pressure from velocities and bulk
	//  p_new = p_old - (dv/dx + dv/dy + dv/dz)*K*dt


	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar IKdt, theta;
	PetscScalar ***up,  ***p; //, pn;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access work vectors

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr); //here is current pressure
	//ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp,   &up);  CHKERRQ(ierr); //here is theta (dvx/dx + dvy/dy + dvz/dz)


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
		//pn    = svBulk->pn;    // pressure history

		////why no directly svBulk->theta?
		//p[k][j][i] -= up[k][j][i]/IKdt;
		p[k][j][i] -= svBulk->theta/IKdt; //theta=dv/dx + dv/dy + dv/dz
	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN,   jr->lp,  &p);   CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(fs->DA_CEN,   jr->gp,  &up);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_CEN, jr->lp, jr->gp)


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/

//
/* //-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetPressure2"
PetscErrorCode GetPressure2(JacRes *jr)
{
	//  ... comment

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar dt, K, IKdt, theta, dx, dy, dz;
	PetscScalar ***p;
	PetscScalar ***vx,  ***vy,  ***vz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;
	dt =  jr->ts.dt;     // time step

	// access work vectors

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr); //here is current pressure

	// access local (ghosted) velocity components
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);


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

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// just to try /////////////////////////////
		// volumetric strain rate
		theta = (vx[k][j][i+1] - vx[k][j][i])/dx + (vy[k][j+1][i] - vy[k][j][i])/dy + (vz[k+1][j][i] - vz[k][j][i])/dz;
		//PetscPrintf(PETSC_COMM_WORLD, "    [k,j,i]  = [%i,%i,%i]  theta in pressure2  = %12.12e  \n", k,j,i,theta);
		// K
		K=3.38e10;
		////////////////////////////////////////////

		// calculate new pressure
		//p[k][j][i] = p[k][j][i]-(svBulk->theta)/IKdt;
		p[k][j][i] = p[k][j][i]-(theta)*K*dt ;
	}
	END_STD_LOOP


	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &p);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_CEN, jr->lp, jr->gp)


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
 */

 /*
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetStress"
PetscErrorCode GetStress(JacRes *jr)
{
	//  ... comment

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	MatParLim  *matLim;
	PetscInt    iter, numPhases;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, i_src, j_src, k_src;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar dx, dy, dz;
	PetscScalar gx, gy, gz, tx, ty, tz, sxx, syy, szz, sxy, sxz, syz;
	PetscScalar J2Inv, theta, rho, Tc, pc, pShift, dt, fssa, *grav,time;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***up, ***T;
	PetscScalar eta_creep;
	PetscScalar mu, lamda;;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	mu 		= 	2.3400e+10; /////////////////////////////////////////////
	lamda 	=	1.8200e+10;



	// access residual context variables
	dt        =  jr->ts.dt;     // time step


	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	This is calculated in the residual
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

		// access current pressure
		//pc = p[k][j][i];

		// compute total Cauchy stresses
		//sxx = svCell->sxx - pc;
		//syy = svCell->syy - pc;
		//szz = svCell->szz - pc;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// stress
		svCell->sxx = svCell->sxx +( (lamda+2*mu )*(vx[k][j][i+1] - vx[k][j][i])/dx + lamda*( (vy[k][j+1][i] - vy[k][j][i])/dy + (vz[k+1][j][i] - vz[k][j][i])/dz))*dt;
		svCell->syy = svCell->sxx +( (lamda+2*mu )*(vy[k][j+1][i] - vy[k][j][i])/dy + lamda*( (vx[k][j][i+1] - vx[k][j][i])/dx + (vz[k+1][j][i] - vz[k][j][i])/dz))*dt;
		svCell->szz = svCell->szz +( (lamda+2*mu )*(vz[k+1][j][i] - vz[k][j][i])/dz + lamda*( (vx[k][j][i+1] - vx[k][j][i])/dx + (vy[k][j+1][i] - vy[k][j][i])/dy))*dt;


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
		svDev  = &svEdge->svDev;

		// access xy component of the Cauchy stress
		//sxy = svEdge->s;

		// get mesh steps for the backward and forward derivatives
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// stress
		svEdge->s = svEdge->s  + mu * ((vy[k][j][i] - vy[k][j][i-1])/dx + (vx[k][j][i] - vx[k][j-1][i])/dy)*dt;

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
		svDev  = &svEdge->svDev;

		// access xz component of the Cauchy stress
		//sxz = svEdge->s;

		// get mesh steps for the backward and forward derivatives
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// stress
		svEdge->s = svEdge->s  + mu * ((vz[k][j][i] - vz[k][j][i-1])/dx + (vx[k][j][i] - vx[k-1][j][i])/dz)*dt;

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
		svDev  = &svEdge->svDev;

		// access yz component of the Cauchy stress
		//syz = svEdge->s;

		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// stress
		svEdge->s = svEdge->s  + mu * ((vz[k][j][i] - vz[k][j-1][i])/dy + (vy[k][j][i] - vy[k-1][j][i])/dz)*dt;

	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	//GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gp, jr->lp)



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/



/*//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetMomentumResidual2"
PetscErrorCode JacResGetMomentumResidual2(JacRes *jr)
{
	// Compute residual of nonlinear momentum conservation

	// ...
	// equations, based on pre-computed components of effective
	// strain-rate tensor, current values of pressure and temperature.
	// Missing components of the second invariant of the effective strain-rate
	// tensor (squares of the corresponding strain rate components) are averaged
	// form the hosting nodes using arithmetic mean.
	// DII = (0.5*D_ij*D_ij)^0.5
	// NOTE: we interpolate and average D_ij*D_ij terms instead of D_ij

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	MatParLim  *matLim;
	PetscInt    iter, numPhases;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, i_src, j_src, k_src;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar gx, gy, gz, tx, ty, tz, sxx, syy, szz, sxy, sxz, syz;
	PetscScalar J2Inv, theta, rho, Tc, pc, pShift, dt, fssa, *grav,time;
	PetscScalar ***fx,  ***fy,  ***fz, ***vx,  ***vy,  ***vz, ***gc;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***up, ***T;
	PetscScalar eta_creep;
	PetscScalar depth;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;


//	PetscInt mcz = fs->dsz.tcels - 1;

	// initialize maximum node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// access residual context variables
	numPhases =  jr->numPhases; // number phases
	phases    =  jr->phases;    // phase parameters
	matLim    = &jr->matLim;    // phase parameters limiters
	dt        =  jr->ts.dt;     // time step
	fssa      =  jr->FSSA;      // density gradient penalty parameter
	grav      =  jr->grav;      // gravity acceleration
	pShift    =  jr->pShift;    // pressure shift

	// clear local residual vectors
	ierr = VecZeroEntries(jr->lfx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->lfz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->gc);  CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->gp,   &up);  CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

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

		//=================
		// SECOND INVARIANT
		//=================

		// access strain rates
		XX = dxx[k][j][i];
		YY = dyy[k][j][i];
		ZZ = dzz[k][j][i];

		// x-y plane, i-j indices
		XY1 = dxy[k][j][i];
		XY2 = dxy[k][j+1][i];
		XY3 = dxy[k][j][i+1];
		XY4 = dxy[k][j+1][i+1];

		// x-z plane, i-k indices
		XZ1 = dxz[k][j][i];
		XZ2 = dxz[k+1][j][i];
		XZ3 = dxz[k][j][i+1];
		XZ4 = dxz[k+1][j][i+1];

		// y-z plane, j-k indices
		YZ1 = dyz[k][j][i];
		YZ2 = dyz[k+1][j][i];
		YZ3 = dyz[k][j+1][i];
		YZ4 = dyz[k+1][j+1][i];

		// compute second invariant
		J2Inv = 0.5*(XX*XX + YY*YY + ZZ*ZZ) +
		0.25*(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25*(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
		0.25*(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);


		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// current temperature
		Tc = T[k][j][i];

		// access current pressure
		pc = p[k][j][i];

// compute depth below the free surface
depth = jr->avg_topo - COORD_CELL(k, sz, fs->dsz);

if(depth < 0.0) depth = 0.0;

// evaluate volumetric constitutive equations
ierr = VolConstEq(svBulk, numPhases, phases, svCell->phRat, matLim, depth, dt, pc-pShift , Tc); CHKERRQ(ierr);

// solve mass conservation equation 1/K dp/dt + (dvx/dx + dvy/dy + dvz/dz) = 0
//ierr = GetPressure2(jr); 	CHKERRQ(ierr);

// access current pressure
pc = p[k][j][i];


//		// access current pressure
//		pc = p[k][j][i];



//		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, numPhases, phases, svCell->phRat, matLim, dt, pc-pShift, Tc); CHKERRQ(ierr);

		// store creep viscosity
		svCell->eta_creep = eta_creep;


		// compute stress, plastic strain rate and shear heating term on cell
		ierr = GetStressCell(svCell, matLim, XX, YY, ZZ); CHKERRQ(ierr);


		// compute total Cauchy stresses
		sxx = svCell->sxx - pc;
		syy = svCell->syy - pc;
		szz = svCell->szz - pc;


		////////////////////////////////////
		// Add seismic source
		PetscScalar t0;
		PetscScalar alfa;
		PetscScalar amplitude;
		t0=1.0;
		alfa=40.0;
		amplitude=10.0;
		if (k==0)
		{
			szz = 0*szz + amplitude;
			sxx = 0*sxx -  amplitude/2.0;
			syy = 0*syy -  amplitude/2.0;
		}
		///////////////////////////////////

// compute depth below the free surface
//depth = jr->avg_topo - COORD_CELL(k, sz, fs->dsz);

//if(depth < 0.0) depth = 0.0;

// evaluate volumetric constitutive equations
//ierr = VolConstEq(svBulk, numPhases, phases, svCell->phRat, matLim, depth, dt, pc-pShift , Tc); CHKERRQ(ierr);

		// access
		theta = svBulk->theta; // volumetric strain rate
		rho   = svBulk->rho;   // effective density
		//IKdt  = svBulk->IKdt;  // inverse bulk viscosity
		//alpha = svBulk->alpha; // effective thermal expansion
		//pn    = svBulk->pn;    // pressure history
		//Tn    = svBulk->Tn;    // temperature history

		// compute gravity terms
		gx = rho*grav[0];
		gy = rho*grav[1];
		gz = rho*grav[2];

		// compute stabilization terms (lumped approximation)
		tx = -fssa*dt*gx;
		ty = -fssa*dt*gy;
		tz = -fssa*dt*gz;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// momentum
		fx[k][j][i] -= (sxx + vx[k][j][i]*tx)/bdx + gx/2.0;   fx[k][j][i+1] += (sxx + vx[k][j][i+1]*tx)/fdx - gx/2.0;
		fy[k][j][i] -= (syy + vy[k][j][i]*ty)/bdy + gy/2.0;   fy[k][j+1][i] += (syy + vy[k][j+1][i]*ty)/fdy - gy/2.0;
		fz[k][j][i] -= (szz + vz[k][j][i]*tz)/bdz + gz/2.0;   fz[k+1][j][i] += (szz + vz[k+1][j][i]*tz)/fdz - gz/2.0;

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
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// access strain rates
		XY = dxy[k][j][i];

		// x-y plane, i-j indices (i & j - bounded)
		XX1 = dxx[k][J1][I1];
		XX2 = dxx[k][J1][I2];
		XX3 = dxx[k][J2][I1];
		XX4 = dxx[k][J2][I2];

		// x-y plane, i-j indices (i & j - bounded)
		YY1 = dyy[k][J1][I1];
		YY2 = dyy[k][J1][I2];
		YY3 = dyy[k][J2][I1];
		YY4 = dyy[k][J2][I2];

		// x-y plane, i-j indices (i & j - bounded)
		ZZ1 = dzz[k][J1][I1];
		ZZ2 = dzz[k][J1][I2];
		ZZ3 = dzz[k][J2][I1];
		ZZ4 = dzz[k][J2][I2];

		// y-z plane j-k indices (j - bounded)
		XZ1 = dxz[k][J1][i];
		XZ2 = dxz[k+1][J1][i];
		XZ3 = dxz[k][J2][i];
		XZ4 = dxz[k+1][J2][i];

		// x-z plane i-k indices (i - bounded)
		YZ1 = dyz[k][j][I1];
		YZ2 = dyz[k+1][j][I1];
		YZ3 = dyz[k][j][I2];
		YZ4 = dyz[k+1][j][I2];

		// compute second invariant
		J2Inv = XY*XY +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
		0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (x-y plane, i-j indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k][j-1][i] + p[k][j-1][i-1]);

		// current temperature (x-y plane, i-j indices)
		Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k][j-1][i] + T[k][j-1][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, numPhases, phases, svEdge->phRat, matLim, dt, pc-pShift, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, XY); CHKERRQ(ierr);

		// access xy component of the Cauchy stress
		sxy = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);

		// momentum
		fx[k][j-1][i] -= sxy/bdy;   fx[k][j][i] += sxy/fdy;
		fy[k][j][i-1] -= sxy/bdx;   fy[k][j][i] += sxy/fdx;


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
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// access strain rates
		XZ = dxz[k][j][i];

		// x-z plane, i-k indices (i & k - bounded)
		XX1 = dxx[K1][j][I1];
		XX2 = dxx[K1][j][I2];
		XX3 = dxx[K2][j][I1];
		XX4 = dxx[K2][j][I2];

		// x-z plane, i-k indices (i & k - bounded)
		YY1 = dyy[K1][j][I1];
		YY2 = dyy[K1][j][I2];
		YY3 = dyy[K2][j][I1];
		YY4 = dyy[K2][j][I2];

		// x-z plane, i-k indices (i & k - bounded)
		ZZ1 = dzz[K1][j][I1];
		ZZ2 = dzz[K1][j][I2];
		ZZ3 = dzz[K2][j][I1];
		ZZ4 = dzz[K2][j][I2];

		// y-z plane, j-k indices (k - bounded)
		XY1 = dxy[K1][j][i];
		XY2 = dxy[K1][j+1][i];
		XY3 = dxy[K2][j][i];
		XY4 = dxy[K2][j+1][i];

		// xy plane, i-j indices (i - bounded)
		YZ1 = dyz[k][j][I1];
		YZ2 = dyz[k][j+1][I1];
		YZ3 = dyz[k][j][I2];
		YZ4 = dyz[k][j+1][I2];

		// compute second invariant
		J2Inv = XZ*XZ +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25 *(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (x-z plane, i-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j][i-1] + p[k-1][j][i] + p[k-1][j][i-1]);

		// current temperature (x-z plane, i-k indices)
		Tc = 0.25*(T[k][j][i] + T[k][j][i-1] + T[k-1][j][i] + T[k-1][j][i-1]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, numPhases, phases, svEdge->phRat, matLim, dt, pc-pShift, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, XZ); CHKERRQ(ierr);

		// access xz component of the Cauchy stress
		sxz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_CELL(i-1, sx, fs->dsx);   fdx = SIZE_CELL(i, sx, fs->dsx);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// momentum
		fx[k-1][j][i] -= sxz/bdz;   fx[k][j][i] += sxz/fdz;
		fz[k][j][i-1] -= sxz/bdx;   fz[k][j][i] += sxz/fdx;

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
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// access strain rates
		YZ = dyz[k][j][i];

		// y-z plane, j-k indices (j & k - bounded)
		XX1 = dxx[K1][J1][i];
		XX2 = dxx[K1][J2][i];
		XX3 = dxx[K2][J1][i];
		XX4 = dxx[K2][J2][i];

		// y-z plane, j-k indices (j & k - bounded)
		YY1 = dyy[K1][J1][i];
		YY2 = dyy[K1][J2][i];
		YY3 = dyy[K2][J1][i];
		YY4 = dyy[K2][J2][i];

		// y-z plane, j-k indices (j & k - bounded)
		ZZ1 = dzz[K1][J1][i];
		ZZ2 = dzz[K1][J2][i];
		ZZ3 = dzz[K2][J1][i];
		ZZ4 = dzz[K2][J2][i];

		// x-z plane, i-k indices (k -bounded)
		XY1 = dxy[K1][j][i];
		XY2 = dxy[K1][j][i+1];
		XY3 = dxy[K2][j][i];
		XY4 = dxy[K2][j][i+1];

		// x-y plane, i-j indices (j - bounded)
		XZ1 = dxz[k][J1][i];
		XZ2 = dxz[k][J1][i+1];
		XZ3 = dxz[k][J2][i];
		XZ4 = dxz[k][J2][i+1];

		// compute second invariant
		J2Inv = YZ*YZ +
		0.125*(XX1*XX1 + XX2*XX2 + XX3*XX3 + XX4*XX4) +
		0.125*(YY1*YY1 + YY2*YY2 + YY3*YY3 + YY4*YY4) +
		0.125*(ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4) +
		0.25 *(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
		0.25 *(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4);

		// store square root of second invariant
		svDev->DII = sqrt(J2Inv);

		//=======================
		// CONSTITUTIVE EQUATIONS
		//=======================

		// access current pressure (y-z plane, j-k indices)
		pc = 0.25*(p[k][j][i] + p[k][j-1][i] + p[k-1][j][i] + p[k-1][j-1][i]);

		// current temperature (y-z plane, j-k indices)
		Tc = 0.25*(T[k][j][i] + T[k][j-1][i] + T[k-1][j][i] + T[k-1][j-1][i]);

		// evaluate deviatoric constitutive equations
		ierr = DevConstEq(svDev, &eta_creep, numPhases, phases, svEdge->phRat, matLim, dt, pc-pShift, Tc); CHKERRQ(ierr);

		// compute stress, plastic strain rate and shear heating term on edge
		ierr = GetStressEdge(svEdge, matLim, YZ); CHKERRQ(ierr);

		// access yz component of the Cauchy stress
		syz = svEdge->s;

		//=========
		// RESIDUAL
		//=========

		// get mesh steps for the backward and forward derivatives
		bdy = SIZE_CELL(j-1, sy, fs->dsy);   fdy = SIZE_CELL(j, sy, fs->dsy);
		bdz = SIZE_CELL(k-1, sz, fs->dsz);   fdz = SIZE_CELL(k, sz, fs->dsz);

		// update momentum residuals
		fy[k-1][j][i] -= syz/bdz;   fy[k][j][i] += syz/fdz;
		fz[k][j-1][i] -= syz/bdy;   fz[k][j][i] += syz/fdy;

	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gc,   &gc);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->gp,   &up);   CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &dxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &dxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &dyz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lfx,  &fx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lfy,  &fy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lfz,  &fz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,  &vz);  CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gp, jr->lp)

	// assemble global residuals from local contributions
	LOCAL_TO_GLOBAL(fs->DA_X, jr->lfx, jr->gfx)
	LOCAL_TO_GLOBAL(fs->DA_Y, jr->lfy, jr->gfy)
	LOCAL_TO_GLOBAL(fs->DA_Z, jr->lfz, jr->gfz)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/

/*//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SolveEquationsWave"
PetscErrorCode SolveEquationsWave(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get average pressure near the top surface
	ierr = JacResGetPressShift(jr); CHKERRQ(ierr);

	// compute effective strain rate
	ierr = JacResGetEffStrainRate(jr); CHKERRQ(ierr);


	//// solve mass conservation equation 1/K dp/dt + (dvx/dx + dvy/dy + dvz/dz) = 0
	ierr = GetPressure2(jr); 	CHKERRQ(ierr);

	JacResGetMomentumResidual2(jr);

	// calculate velocities
	ierr = GetVelocities(jr);	CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/

/*//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PutSeismicSource"
PetscErrorCode PutSeismicSource(JacRes *jr, AdvCtx *actx, UserCtx *user)
{

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarBulk *svBulk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscInt     ii, jj, ID, I, J, K;
	Marker      *P;

	PetscScalar time;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscScalar t0;
	PetscScalar alfa;
	PetscScalar amplitude;

	fs = actx->fs;
	jr = actx->jr;

	// number of cells
	nx = fs->dsx.ncels;
	ny = fs->dsy.ncels;

	t0			=	1.0;
	alfa		=	40.0;
	amplitude	=	10.0;

	if (jr->SourceParams.source_type == POINT)
	{

//
	}
	else if (jr->SourceParams.source_type == PLANE)
	{
		time	  =  JacResGetTime(jr);

		// scan ALL markers
		for(jj = 0; jj < actx->nummark; jj++)
		{
			// access next marker
			P = &actx->markers[jj];

			// get consecutive index of the host cell
			ID = actx->cellnum[jj];

			// expand I, J, K cell indices
			GET_CELL_IJK(ID, I, J, K, nx, ny)

			if (K==0) {

				P->S.xx =	-  amplitude/2.0; // *exp(-alfa*((time-t0)*(time-t0)));
				P->S.yy = 	-  amplitude/2.0; // *exp(-alfa*((time-t0)*(time-t0)));
				P->S.zz =  	   amplitude; // *exp(-alfa*((time-t0)*(time-t0)));
			}
		}
	}

	PetscFunctionReturn(0);

}
*/

/*//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PrintStress"
PetscErrorCode PrintStress(JacRes *jr)
{
	//

	FDSTAG     *fs;
	SolVarCell *svCell;

	PetscInt i, j, k, iter;
	PetscInt nx, ny, nz, sx, sy, sz;
	PetscFunctionBegin;

	fs = jr->fs;


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



			PetscPrintf(PETSC_COMM_WORLD, "    svCell->hxx, svCell->sxx in [k,j,i]  = [%i,%i,%i] are  %12.12e, %12.12e \n", k,j,i,svCell->hxx, svCell->sxx);

		}
		END_STD_LOOP

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/


//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SaveVelocitiesForSeismicStation"
PetscErrorCode SaveVelocitiesForSeismicStation(JacRes *jr, UserCtx *user)
{

	//  ... comment

	FDSTAG     *fs;
//	PetscInt    iter;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, i_rec, j_rec, k_rec;

	PetscScalar /*dt,*/ t, vx_rec, vy_rec, vz_rec;
	PetscScalar ***vx,  ***vy,  ***vz;
	FILE      *fseism;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// access residual context variables

//	dt    =  jr->ts.dt;     // time step
	t	  =  JacResGetTime(jr);

	// file for seismic signals
	fseism = user->StationParams.output_file;


	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);

	//-------------------------------
	// side points
	//-------------------------------
//	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	// Save coordinates for station
	// station position (to improve getting from input file)
	i_rec=nx/2;
	j_rec=ny/2;
	k_rec=nz/2;


	START_STD_LOOP
	{
		if (i==i_rec && j==j_rec && k==k_rec) {
			vx_rec=vx[k][j][i]+vx[k][j][i-1]/2.0;
			break;
		}
	}
	END_STD_LOOP

//	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		if (i==i_rec && j==j_rec && k==k_rec) {
			vy_rec=vy[k][j][i]+vy[k][j-1][i]/2.0;
			break;
		}
	}
	END_STD_LOOP

//	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)
	START_STD_LOOP
	{
		if (i==i_rec && j==j_rec && k==k_rec) {
			vz_rec=vz[k][j][i]+vz[k-1][j][i]/2.0;
			break;
		}
	}
	END_STD_LOOP

	// restore vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->gvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->gvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->gvz,  &vz);  CHKERRQ(ierr);


	fprintf(fseism, "%12.12e %12.12e %12.12e %12.12e\n", t, vx_rec, vy_rec, vz_rec);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "ShowValues"
PetscErrorCode ShowValues(JacRes *jr, UserCtx *user, PetscInt n)
{
	// Show the values of velocity, pressure, stress, strain,...

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


	fs = jr->fs;
	dt = user->dt;

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
		//PetscPrintf(PETSC_COMM_WORLD, "    svCell.hxx,yy,zz[%i,%i,%i]  = (%12.12e,%12.12e,%12.12e) \n", i,j,k, svCell->hxx,svCell->hyy,svCell->hzz);
		//fprintf(fseism, "%i %i %i\n", i,j,k);
		//fprintf(fseism, "%12.12e\n", gfz[k][j][i]);
		//if ( p[k][j][i] != up[k][j][i])			fprintf(fseism, "%12.12e\n", p[k][j][i]);
		//if ( fx[k][j][i] != gfx[k][j][i])		 fprintf(fseism, "%12.12e\n", fx[k][j][i]);
		//PetscPrintf(PETSC_COMM_WORLD, "    p[%i,%i,%i]  = %12.12e \n", i,j,k, p[k][j][i]);
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Apply source
#undef __FUNCT__
#define __FUNCT__ "GetCellCoordinatesSourceAndSeismicStation"
PetscErrorCode GetCellCoordinatesSourceAndSeismicStation(JacRes *jr)
{
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz, ii, jj, kk, iii, jjj, kkk;
	PetscScalar    xs[3], xe[3];
	FDSTAG     *fs;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ii = -1;
	jj = -1;
	kk = -1;
	iii = -1;
	jjj = -1;
	kkk = -1;

	fs = jr->fs;

	if ((jr->SeismicSource==PETSC_TRUE && (jr->SourceParams.source_type==POINT || jr->SourceParams.source_type==MOMENT)) || jr->SeismicStation==PETSC_TRUE ) {
		//-------------------------------
		// central points
		//-------------------------------
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)
		START_STD_LOOP
		{
			PetscScalar x,y,z,dx,dy,dz;

			// get coordinate of center of cell 
			x  = COORD_CELL(i, sx, fs->dsx);
			y  = COORD_CELL(j, sy, fs->dsy);
			z  = COORD_CELL(k, sz, fs->dsz);

			// size of cell
			dx = SIZE_CELL(i, sx, fs->dsx);		
			dy = SIZE_CELL(j, sy, fs->dsy);
			dz = SIZE_CELL(k, sz, fs->dsz);

			// get cell coordinates
			xs[0] = x-dx/2; xe[0] = x+dx/2;
			xs[1] = y-dy/2; xe[1] = y+dy/2;
			xs[2] = z-dz/2; xe[2] = z+dz/2;
			
			if (jr->SeismicSource==PETSC_TRUE)
			{
				if (jr->SourceParams.x > xs[0] && jr->SourceParams.x <= xe[0] && jr->SourceParams.y > xs[1] && jr->SourceParams.y <= xe[1] && jr->SourceParams.z > xs[2] && jr->SourceParams.z <= xe[2])
				{
					ii=i;
					jj=j;
					kk=k;
				}
			}
			if (jr->SeismicStation==PETSC_TRUE)
			{
				if (jr->StationParams.x > xs[0] && jr->StationParams.x <= xe[0] && jr->StationParams.y > xs[1] && jr->StationParams.y <= xe[1] && jr->StationParams.z > xs[2] && jr->StationParams.z <= xe[2])
				{
					iii=i;
					jjj=j;
					kkk=k;
				}
			}
		}
		END_STD_LOOP
		
		if (jr->SeismicSource==PETSC_TRUE)
		{
			if(ISParallel(PETSC_COMM_WORLD))
			{
				ierr = MPI_Allreduce(&ii, &jr->SourceParams.i, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
				ierr = MPI_Allreduce(&jj, &jr->SourceParams.j, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
				ierr = MPI_Allreduce(&kk, &jr->SourceParams.k, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
			}
			else
			{
				jr->SourceParams.i = ii;
				jr->SourceParams.j = jj;
				jr->SourceParams.k = kk;
			}
		}
		if (jr->SeismicStation==PETSC_TRUE)
		{
			if(ISParallel(PETSC_COMM_WORLD))
			{
				ierr = MPI_Allreduce(&iii, &jr->StationParams.i, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
				ierr = MPI_Allreduce(&jjj, &jr->StationParams.j, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
				ierr = MPI_Allreduce(&kkk, &jr->StationParams.k, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
			}
			else
			{
				jr->StationParams.i = iii;
				jr->StationParams.j = jjj;
				jr->StationParams.k = kkk;
			}
		}
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Apply source
#undef __FUNCT__
#define __FUNCT__ "GetStressFromSource"
PetscErrorCode GetStressFromSource(JacRes *jr, UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz)
{
	PetscInt M, N, P;
	PetscScalar /*coor,*/ time;
	PetscScalar    xs[3], xe[3];

	time	  =  JacResGetTime(jr);

	if (jr->SourceParams.source_type == PLANE)
	{
		//if (k==user->nel_z-1) //(k==1)
		if (i==user->nel_x-1)
		{
			*sxx = 		jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
			*syy =	- 	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;
			*szz = 	-  	jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;
		}
	}
	else if (jr->SourceParams.source_type == COMPRES)
		{
			/*if (k==1)
			{
				*szz = 		jr->SourceParams.amplitude;
				*sxx =	- 	jr->SourceParams.amplitude/2.0;
				*syy = 	-  	jr->SourceParams.amplitude/2.0;
			}
			else */
			if (k==user->nel_z-1)
			{
				*szz = 	-	jr->SourceParams.amplitude;
				//*sxx =	- 	jr->SourceParams.amplitude/2.0;
				//*syy = 	-  	jr->SourceParams.amplitude/2.0;
			}
		}
	else if (jr->SourceParams.source_type == POINT)
	{

		// get cell coordinates
		//xs[0] = jr->fs->dsx.ncoor[i]; xe[0] = jr->fs->dsx.ncoor[i+1];
		//xs[1] = jr->fs->dsy.ncoor[j]; xe[1] = jr->fs->dsy.ncoor[j+1];
		//xs[2] = jr->fs->dsz.ncoor[k]; xe[2] = jr->fs->dsz.ncoor[k+1];
		//if (jr->SourceParams.x > xs[0] && jr->SourceParams.x <= xe[0] && jr->SourceParams.y > xs[1] && jr->SourceParams.y <= xe[1] && jr->SourceParams.z > xs[2] && jr->SourceParams.z <= xe[2])

		if (jr->SourceParams.i == i && jr->SourceParams.j == j && jr->SourceParams.k == k)
		{
			*sxx = jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)));
			*syy = jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;
			*szz = -jr->SourceParams.amplitude*exp(-jr->SourceParams.alfa*((time-jr->SourceParams.t0)*(time-jr->SourceParams.t0)))/2;
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



