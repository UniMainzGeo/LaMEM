//---------------------------------------------------------------------------
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "bc.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCreate"
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(bc, sizeof(BCCtx)); CHKERRQ(ierr);

	// create boundary conditions vectors (velocity, pressure, temperature)
	ierr = DMCreateLocalVector(fs->DA_X,   &bc->bcvx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &bc->bcvy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &bc->bcvz);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcp);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcT);   CHKERRQ(ierr);

	// single-point constraints (combined)
	bc->numSPC  = 0;
	bc->SPCList = NULL;

	// single-point constraints (pressure)
	bc->numSPCPres  = 0;
	bc->SPCListPres = NULL;

	// two-point constraints
	bc->numTPC       = 0;
	bc->TPCList      = NULL;
	bc->TPCPrimeDOF  = NULL;
	bc->TPCVals      = NULL;
	bc->TPCLinComPar = NULL;

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
	// WARNING !!! PRESSUE CONSTRAINS ARE CURRENTLY NOT ALLOWED !!!
	// *************************************************************

//	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    start, ln, numSPC, *SPCList;
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
	ierr = makeIntArray(&SPCList, NULL, ln); CHKERRQ(ierr);

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
		if(j == 0 || j == mny) { bcvy[k][j][i] = 0.0; SPCList[numSPC++] = start; }
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
		if(k == 0 || k == mnz) { bcvz[k][j][i] = 0.0; SPCList[numSPC++] = start; }
		start++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz); CHKERRQ(ierr);

	// store constraints
	bc->numSPC  = numSPC;
	bc->SPCList = SPCList;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
