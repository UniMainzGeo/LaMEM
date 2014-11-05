//---------------------------------------------------------------------------
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "bc.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * replace BC input specification consistently in the entire code
// * open box & Winkler (with tangential viscous friction)
// * tangential velocities
// * extend two-point constraint specification & (possibly) get rid bc-vectors
// * create bc-object only at fine level, coarse levels should have simple access
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCClear"
PetscErrorCode BCClear(BCCtx *bc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(bc, sizeof(BCCtx)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCreate"
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create boundary conditions vectors (velocity, pressure, temperature)
	ierr = DMCreateLocalVector(fs->DA_X,   &bc->bcvx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &bc->bcvy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &bc->bcvz);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcp);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcT);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCDestroy"
PetscErrorCode BCDestroy(BCCtx *bc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// destroy boundary conditions vectors (velocity, pressure, temperature)
	ierr = VecDestroy(&bc->bcvx);    CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcvy);    CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcvz);    CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcp);     CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcT);     CHKERRQ(ierr);

	// single-point constraints (combined)
	ierr = PetscFree(bc->SPCList);   CHKERRQ(ierr);
	ierr = PetscFree(bc->SPCVals);   CHKERRQ(ierr);

	// single-point constraints (pressure)
	ierr = PetscFree(bc->SPCListPres);   CHKERRQ(ierr);

	// two-point constraints
	ierr = PetscFree(bc->TPCList);      CHKERRQ(ierr);
	ierr = PetscFree(bc->TPCPrimeDOF);  CHKERRQ(ierr);
	ierr = PetscFree(bc->TPCVals);      CHKERRQ(ierr);
	ierr = PetscFree(bc->TPCLinComPar); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitBC"
PetscErrorCode FDSTAGInitBC(BCCtx *bc, FDSTAG *fs, idxtype idxmod)
{
	// initialize boundary conditions vectors
	PetscBool   flg;
	PetscScalar pgrad;
	PetscInt    mnx, mny, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    start, ln, numSPC, *SPCList;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcp;
	PetscScalar ***ivx,   ***ivy,   ***ivz,  ***ip;
	DOFIndex    *id;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetScalar(PETSC_NULL, "-pgrad", &pgrad, &flg); CHKERRQ(ierr);
	if(flg != PETSC_TRUE) pgrad = 1.0;

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mcz = fs->dsz.tcels - 1;

	if(idxmod == IDXCOUPLED)   id = &fs->dofcoupl;
	if(idxmod == IDXUNCOUPLED) id = &fs->dofsplit;

	// get total number of local matrix rows & global index of the first row
	start = id->istart;
	ln    = id->numdof;

	// allocate SPC arrays
	ierr = makeIntArray(&SPCList, NULL, ln); CHKERRQ(ierr);

	// mark all variables unconstrained
	ierr = VecSet(bc->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcp,  DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcT,  DBL_MAX); CHKERRQ(ierr);

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   id->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   id->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   id->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, id->ip,   &ip);   CHKERRQ(ierr);

	// set face-normal velocities to zero, count & store SPC information
	numSPC = 0.0;

	//=========================
	// SINGLE-POINT CONSTRAINTS
	//=========================

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// boundary x-normal points only
		if(i == 0 || i == mnx) { bcvx[k][j][i] = 0.0; SPCList[numSPC++] = start; }
		start++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// boundary y-normal points only
		if(j == 0 || j == mny) { bcvy[k][j][i] = 0.0; SPCList[numSPC++] = start; }
		start++;
	}
	END_STD_LOOP

	//======================
	// TWO-POINT CONSTRAINTS
	//======================

	//---------
	// X points
	//---------
	GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// all ghost points excluding z-tangential
		if(ivx[k][j][i] == -1 && k >= 0 && k <= mcz) bcvx[k][j][i] = 0.0;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_NODE_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// all ghost points excluding z-tangential
		if(ivy[k][j][i] == -1 && k >= 0 && k <= mcz) bcvy[k][j][i] = 0.0;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_NODE_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// all ghost points
		if(ivz[k][j][i] == -1) bcvz[k][j][i] = 0.0;
	}
	END_STD_LOOP

	//----------------
	// central points
	//---------------
	GET_CELL_RANGE_GHOST_ALL(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_ALL(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_ALL(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(ip[k][j][i] == -1)
		{
			// bottom ghost points (zero pressure)
			if(k < 0) bcp[k][j][i] = 0.0;

			// top ghost points (unit pressure)
			if(k > mcz) bcp[k][j][i] = pgrad;
		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   id->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   id->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   id->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, id->ip,   &ip);   CHKERRQ(ierr);

	// store constraints
	bc->numSPC  = numSPC;
	bc->SPCList = SPCList;

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCInit"
PetscErrorCode BCInit(BCCtx *bc, FDSTAG *fs, idxtype idxmod)
{
	// initialize boundary conditions vectors

	// *************************************************************
	// WARNING !!! AD-HOC FREE-SLIP BOX IS CURRENTLY ASSUMED    !!!
	// WARNING !!! PRESSURE CONSTRAINS ARE CURRENTLY NOT ALLOWED !!!
	// *************************************************************

//	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    start, ln, numSPC, *SPCList;
	PetscScalar *SPCVals;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz;
	DOFIndex    *dof;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

//	mcx = fs->dsx.tcels - 1;
//	mcy = fs->dsy.tcels - 1;
//	mcz = fs->dsz.tcels - 1;

	if(idxmod == IDXCOUPLED)   dof = &fs->cdof;
	if(idxmod == IDXUNCOUPLED) dof = &fs->udof;

	// get total number of local matrix rows & global index of the first row
	start = dof->istart;
	ln    = dof->numdof;

	// allocate SPC arrays
	ierr = makeIntArray (&SPCList, NULL, ln); CHKERRQ(ierr);
	ierr = makeScalArray(&SPCVals, NULL, ln); CHKERRQ(ierr);

	// mark all variables unconstrained
	ierr = VecSet(bc->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcp,  DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcT,  DBL_MAX); CHKERRQ(ierr);

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz); CHKERRQ(ierr);

	// set face-normal velocities to zero, count & store SPC information
	numSPC = 0.0;

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(i == 0 || i == mnx)
		{
			bcvx[k][j][i]   = 0.0;
			SPCList[numSPC] = start;
			SPCVals[numSPC] = 0.0;
			numSPC++;
		}
		start++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(j == 0 || j == mny)
		{
			bcvy[k][j][i]   = 0.0;
			SPCList[numSPC] = start;
			SPCVals[numSPC] = 0.0;
			numSPC++;
		}
		start++;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(k == 0 || k == mnz)
		{
			bcvz[k][j][i]   = 0.0;
			SPCList[numSPC] = start;
			SPCVals[numSPC] = 0.0;
			numSPC++;
		}
		start++;
	}
	END_STD_LOOP

	// add pushing BC
	//ierr = PBCGetIndices(bc, fs, bcvx, bcvy, SPCList, numSPC, start); CHKERRQ(ierr);

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz); CHKERRQ(ierr);

	// store constraints
	bc->numSPC  = numSPC;
	bc->SPCList = SPCList;
	bc->SPCVals = SPCVals;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PBCInit"
PetscErrorCode PBCInit(BCCtx *bc, UserContext *user)
{
	// MUST be called at the beginning of time step before setting boundary conditions
	// initialize pushing boundary conditions parameters

	PetscInt    i, ichange;
	PetscScalar time0, time1, vpush, omega;
	PetscScalar charTime, charNVel, charRVel;
	PetscScalar Vx, Vy, theta;
	PetscInt    c_adv, dir;

	PetscFunctionBegin;

	// set boundary conditions flag
	bc->pflag = PETSC_FALSE;

	// check if pushing option is activated
	if (user->AddPushing == 0) PetscFunctionReturn(0);

	// add pushing boundary conditions ONLY within the specified time interval - for that introduce a new flag
	if ((user->time>=user->Pushing.time[0]) && (user->time<=user->Pushing.time[user->Pushing.num_changes]))
	{
		// check which pushing stage
		ichange = 0;

		for (i = 0; i < user->Pushing.num_changes; i++)
		{
			if (user->time>=user->Pushing.time[i] && user->time<=user->Pushing.time[i+1]){
				ichange = i;
			}
		}
		user->Pushing.ind_change = ichange;

		// initialize parameters for the time step
		Vx    = 0.0;
		Vy    = 0.0;
		theta = user->Pushing.theta;

		if (user->Pushing.dir[ichange] == 0)
		{
			Vx    = cos(theta)*user->Pushing.V_push[ichange];
			Vy    = sin(theta)*user->Pushing.V_push[ichange];
		}

		if (user->Pushing.dir[ichange] == 1) Vx = user->Pushing.V_push[ichange];
		if (user->Pushing.dir[ichange] == 2) Vy = user->Pushing.V_push[ichange];

		// set boundary conditions parameters
		bc->pflag   = PETSC_TRUE;
		bc->vval[0] = Vx;
		bc->vval[1] = Vy;
		bc->theta   = theta;

		// prepare block coordinates
		bc->xbs[0] = user->Pushing.x_center_block - user->Pushing.L_block*0.5; //left
		bc->xbs[1] = user->Pushing.y_center_block - user->Pushing.W_block*0.5; //front
		bc->xbs[2] = user->Pushing.z_center_block - user->Pushing.H_block*0.5; //bottom
		bc->xbe[0] = user->Pushing.x_center_block + user->Pushing.L_block*0.5; //right
		bc->xbe[1] = user->Pushing.y_center_block + user->Pushing.W_block*0.5; //back
		bc->xbe[2] = user->Pushing.z_center_block + user->Pushing.H_block*0.5; //top

		// print useful info
		charTime = 1e-6/user->Characteristic.SecYear*user->Characteristic.Time;     // time in Ma
		charNVel = 1e2*user->Characteristic.SecYear*user->Characteristic.Velocity;  // velocity in cm/yr
		charRVel = 180/M_PI*user->Characteristic.SecYear/user->Characteristic.Time; // rotation velocity in deg/yr

		time0 = user->Pushing.time        [ichange  ]*charTime;
		time1 = user->Pushing.time        [ichange+1]*charTime;
		vpush = user->Pushing.V_push      [ichange  ]*charNVel;
		omega = user->Pushing.omega       [ichange  ]*charRVel;
		c_adv = user->Pushing.coord_advect[ichange  ];
		dir   = user->Pushing.dir         [ichange  ];

		PetscPrintf(PETSC_COMM_WORLD,"#  Pushing BC: Time interval = [%g - %g] Ma, V_push = %g cm/yr, Omega = %g deg/yr, mobile_block = %d, direction = %d\n",time0,time1,vpush,omega,c_adv,dir);
		PetscPrintf(PETSC_COMM_WORLD,"#  Pushing BC: Block center coordinates x = [%g], y = [%g], z = [%g]\n",user->Pushing.x_center_block*1000,user->Pushing.y_center_block*1000,user->Pushing.z_center_block*1000);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PBCGetIndices"
PetscErrorCode PBCGetIndices(BCCtx *bc, FDSTAG *fs, PetscScalar ***pbcvx, PetscScalar ***pbcvy, PetscInt *SPCListPush, PetscInt numSPCPush, PetscInt start)
{
	// get SPC information
	// initialize internal (pushing) boundary conditions vectors
	// only x, y velocities are constrained!! constraining vz makes a bad case - will not converge.

	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar xbl, xbr, ybf, ybb, zbb, zbt;
	PetscScalar x, y, z, rx, ry, rz;
	PetscScalar theta;

	PetscFunctionBegin;

	// check if pushing is activated in every time step
	if (!bc->pflag) PetscFunctionReturn(0);

	// prepare block coordinates
	xbl = bc->xbs[0]; //left
	xbr = bc->xbe[0]; //right
	ybf = bc->xbs[1]; //front
	ybb = bc->xbe[1]; //back
	zbb = bc->xbs[2]; //bottom
	zbt = bc->xbe[2]; //top

	theta = bc->theta; //rotation angle

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		x = fs->dsx.ncoor[i];
		y = fs->dsy.ccoor[j];
		z = fs->dsz.ccoor[k];

		// clockwise rotation with rotation matrix R = [ cos() sin() ; -sin() cos() ]
		rx = cos(theta)*x + sin(theta)*y;
		ry =-sin(theta)*x + cos(theta)*y;
		rz = z;

		// check for points inside the block
		if(rx > xbl && rx < xbr && ry > ybf && ry < ybb && rz > zbb && rz < zbt)
		{
			pbcvx[k][j][i] = bc->vval[0]; SPCListPush[numSPCPush++] = start;
		}
		start++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		x = fs->dsx.ccoor[i];
		y = fs->dsy.ncoor[j];
		z = fs->dsz.ccoor[k];

		// clockwise rotation with rotation matrix R = [ cos() sin() ; -sin() cos() ]
		rx = cos(theta)*x + sin(theta)*y;
		ry =-sin(theta)*x + cos(theta)*y;
		rz = z;

		// check for points inside the block
		if(rx > xbl && rx < xbr && ry > ybf && ry < ybb && rz > zbb && rz < zbt)
		{
			pbcvy[k][j][i] = bc->vval[1]; SPCListPush[numSPCPush++] = start;
		}
		start++;
	}
	END_STD_LOOP

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PBCAdvectBlock"
PetscErrorCode PBCAdvectBlock(UserContext *user)
{
	// Routine after solver!

	PetscScalar xc, yc , zc, Vx, Vy, dx, dy, dt, omega, theta, dtheta;
	PetscInt    advc, ichange;

	// check if pushing option is activated
	if (user->AddPushing == 0) PetscFunctionReturn(0);

	// check time interval
	if ((user->time>=user->Pushing.time[0]) && (user->time<=user->Pushing.time[user->Pushing.num_changes]))
	{
		// initialize variables
		ichange = user->Pushing.ind_change;
		advc    = user->Pushing.coord_advect[ichange];

		Vx = 0.0;
		Vy = 0.0;
		dy = 0.0;
		dx = 0.0;
		dt = user->dt;

		xc = user->Pushing.x_center_block;
		yc = user->Pushing.y_center_block;
		zc = user->Pushing.z_center_block;

		// block is rotated
		if (user->Pushing.dir[ichange] == 0)
		{
			omega = user->Pushing.omega[ichange];
			theta = user->Pushing.theta;
			Vx    = cos(theta)*user->Pushing.V_push[ichange];
			Vy    = sin(theta)*user->Pushing.V_push[ichange];

			// rotation
			dtheta = omega*dt;
			user->Pushing.theta = theta + dtheta;
		}

		// block moves in X-dir
		if (user->Pushing.dir[ichange] == 1) Vx = user->Pushing.V_push[ichange];

		// block moves in Y-dir
		if (user->Pushing.dir[ichange] == 2) Vy = user->Pushing.V_push[ichange];

		dx     = dt * Vx;
		dy     = dt * Vy;

		// advect the block coordinates
		//stationary block
		if (advc == 0)
		{
			user->Pushing.x_center_block = xc;
			user->Pushing.y_center_block = yc;
			user->Pushing.z_center_block = zc;
		}
		// advect the block
		else
		{
			user->Pushing.x_center_block = xc + dx;
			user->Pushing.y_center_block = yc + dy;
			user->Pushing.z_center_block = zc;
		}
	}

	PetscFunctionReturn(0);
}
