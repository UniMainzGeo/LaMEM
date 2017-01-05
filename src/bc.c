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
 **    filename:   bc.c
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
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * open box & Winkler (with tangential viscous friction)
// * tangential velocities
// * extend two-point constraint specification
//---------------------------------------------------------------------------
// Bezier block functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCBlockCreate"
PetscErrorCode BCBlockCreate(BCBlock *bcb, Scaling *scal, FB *fb)
{
	//	-npath - Number of path points of Bezier curve (end-points only!)
	//	-theta - Orientation angles at path points (counter-clockwise positive)
	//	-time  - Times at path points
	//	-path  - Bezier curve path & control points (6*npath-4 entries are expected)
	//	-npoly - Number of polygon vertices
	//	-poly  - Polygon x-y coordinates at initial time
	//	-bot   - Polygon bottom coordinate
	//	-top   - Polygon top coordinate

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bcb->npath = 2;
	bcb->npoly = 4;

	ierr = getIntParam   (fb, _OPTIONAL_, "npath", &bcb->npath, 1,              _max_path_points_); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "theta",  bcb->theta, bcb->npath,      scal->angle     ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "time",   bcb->time,  bcb->npath,      scal->time      ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "path",   bcb->path,  6*bcb->npath-4,  scal->length    ); CHKERRQ(ierr);

	ierr = getIntParam   (fb, _OPTIONAL_, "npoly", &bcb->npoly, 1,              _max_poly_points_); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "poly ",  bcb->poly,  2*bcb->npoly,    scal->length    ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "bot",   &bcb->bot,   1,               scal->length    ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "top",   &bcb->top,   1,               scal->length    ); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCBlockGetPosition"
PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *f, PetscScalar X[])
{
	// compute position along the path and rotation angle as a function of time

	PetscInt      i, n;
    PetscScalar   r, r2, r3, s, s2, s3;
	PetscScalar  *p1, *p2, *p3, *p4;
	PetscScalar  *path, *theta, *time;

	PetscFunctionBegin;

	n     = bcb->npath;
	path  = bcb->path;
	theta = bcb->theta;
	time  = bcb->time;

	// set flag
	(*f) = 1; if(t < time[0] || t > time[n-1]) { (*f) = 0; PetscFunctionReturn(0); }

	// find time interval
	for(i = 1; i < n-1; i++) { if(t < time[i]) break; } i--;

	// get path and control points
	p1 = path + 6*i;
	p2 = p1 + 2;
	p3 = p2 + 2;
	p4 = p3 + 2;

	// compute interpolation parameters
	r  = (t - time[i])/(time[i+1] - time[i]);
	r2 = r*r;
    r3 = r2*r;
    s  = 1.0 - r;
    s2 = s*s;
    s3 = s2*s;

	// interpolate Bezier path and rotation angle
	X[0] = s3*p1[0] + 3.0*s2*r*p2[0] + 3.0*s*r2*p3[0] + r3*p4[0];
    X[1] = s3*p1[1] + 3.0*s2*r*p2[1] + 3.0*s*r2*p3[1] + r3*p4[1];
    X[2] = s*theta[i] + r*theta[i+1];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCBlockGetPolygon"
PetscErrorCode BCBlockGetPolygon(BCBlock *bcb, PetscScalar Xb[], PetscScalar *cpoly)
{
	// compute current polygon coordinates

	PetscInt     i;
	PetscScalar *xa, *xb;
	PetscScalar  Xa[3], theta, costh, sinth;

	PetscFunctionBegin;

	// get initial polygon position
	Xa[0] = bcb->path[0];
	Xa[1] = bcb->path[1];
	Xa[2] = bcb->theta[0];

	// get rotation matrix
	theta = Xb[2] - Xa[2];
	costh = cos(theta);
	sinth = sin(theta);

	// compute current polygon coordinates
	for(i = 0; i < bcb->npoly; i++)
	{
		// get reference and current points
		xa = bcb->poly + 2*i;
		xb = cpoly     + 2*i;

		// rotate & displace
		RotDispPoint2D(Xa, Xb, costh, sinth, xa, xb);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Dropping boxes functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBoxReadCreate"
PetscErrorCode DBoxReadCreate(DBox *dbox, Scaling *scal, FB *fb)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	//========================
	// Dropping box parameters
	//========================

	ierr = getIntParam(fb, _OPTIONAL_, "dbox_num", &dbox->num, 1, _max_boxes_); CHKERRQ(ierr);

	if(dbox->num)
	{
		ierr = getScalarParam(fb, _REQUIRED_, "dbox_bounds",  dbox->bounds, 6*dbox->num, scal->time    ); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "dbox_zvel",   &dbox->zvel,   1,           scal->velocity); CHKERRQ(ierr);

	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// BCCtx functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCreate"
PetscErrorCode BCCreate(BCCtx *bc, FB *fb)
{
	Scaling     *scal;
	PetscInt     jj;
	PetscScalar  bz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = bc->scal;

	//=====================
	// VELOCITY CONSTRAINTS
	//=====================

	// horizontal background strain-rate parameters
	ierr = getIntParam   (fb, _OPTIONAL_, "ExxNumPeriods",  &bc->ExxNumPeriods,  1,                   _max_periods_    ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "ExxTimeDelims",   bc->ExxTimeDelims,  bc->ExxNumPeriods-1, scal->time       ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "ExxStrainRates",  bc->ExxStrainRates, bc->ExxNumPeriods,   scal->strain_rate); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "EyyNumPeriods",  &bc->EyyNumPeriods,  1,                   _max_periods_    ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "EyyTimeDelims",   bc->EyyTimeDelims,  bc->EyyNumPeriods-1, scal->time       ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "EyyStrainRates",  bc->EyyStrainRates, bc->EyyNumPeriods,   scal->strain_rate); CHKERRQ(ierr);

	// Bezier blocks
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<BCBlockStart>", "<BCBlockEnd>"); CHKERRQ(ierr);

	if(fb->nblocks)
	{
		// error checking
		if(fb->nblocks > _max_boxes_)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many Bezier blocks! found: %lld, max allowed: %lld", (LLD)fb->nblocks, (LLD)_max_boxes_);
		}

		// store actual number of Bezier blocks
		bc->nblocks = fb->nblocks;

		// read Bezier blocks
		for(jj = 0; jj < fb->nblocks; jj++)
		{
			ierr = BCBlockCreate(bc->blocks + jj, scal, fb); CHKERRQ(ierr);

			fb->blockID++;
		}
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// dropping boxes
	ierr = DBoxReadCreate(&bc->dbox, scal, fb); CHKERRQ(ierr);

	// boundary velocities
	ierr = getIntParam(fb, _OPTIONAL_, "bvel_face", &bc->face, 1, -1); CHKERRQ(ierr);

	if(bc->face)
	{
		ierr = getIntParam   (fb, _REQUIRED_, "bvel_phase", &bc->phase, 1, -1            ); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bvel_bot",   &bc->top,   1, scal->length  ); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bvel_top",   &bc->bot,   1, scal->length  ); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bvel_velin", &bc->velin, 1, scal->velocity); CHKERRQ(ierr);

		ierr = FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL); CHKERRQ(ierr);

		// compute outflow velocity
		// INTRODUCE CORRECTION FOR CELL SIZES
		// MUST BE MASS CONSERVATIVE IN DISCRETE SENSE

		bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);
	}

	// simple shear boundary condition
	ierr = getScalarParam(fb, _OPTIONAL_, "bc_gamma_xz", &bc->gamma_xz, 1, scal->strain_rate); CHKERRQ(ierr);

	// open boundary flag
	ierr = getIntParam(fb, _OPTIONAL_, "open_top_bound", &bc->top_open, 1, -1); CHKERRQ(ierr);

	// no-slip boundary condition mask
	ierr = getIntParam(fb, _OPTIONAL_, "noslip", bc->noslip, 6, -1); CHKERRQ(ierr);

	//========================
	// TEMPERATURE CONSTRAINTS
	//========================

	// read from options
	ierr = getScalarParam(fb, _OPTIONAL_, "temp_top", &bc->Tbot, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "temp_top", &bc->Ttop, 1, 1.0); CHKERRQ(ierr);

	// nondimensionalize temperature
	bc->Ttop = (bc->Ttop + scal->Tshift)/scal->temperature;
	bc->Tbot = (bc->Tbot + scal->Tshift)/scal->temperature;

	// allocate vectors
	ierr = BCCreateData(bc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCreateData"
PetscErrorCode BCCreateData(BCCtx *bc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	FDSTAG   *fs  =  bc->fs;
	DOFIndex *dof = &fs->dof;

	// create boundary conditions vectors (velocity, pressure, temperature)
	ierr = DMCreateLocalVector(fs->DA_X,   &bc->bcvx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &bc->bcvy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &bc->bcvz);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcp);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcT);   CHKERRQ(ierr);

	// SPC velocity-pressure
	ierr = makeIntArray (&bc->SPCList, NULL, dof->ln);   CHKERRQ(ierr);
	ierr = makeScalArray(&bc->SPCVals, NULL, dof->ln);   CHKERRQ(ierr);

	// SPC (temperature)
	ierr = makeIntArray (&bc->tSPCList, NULL, dof->lnp); CHKERRQ(ierr);
	ierr = makeScalArray(&bc->tSPCVals, NULL, dof->lnp); CHKERRQ(ierr);

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
	ierr = VecDestroy(&bc->bcvx); CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcvy); CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcvz); CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcp);  CHKERRQ(ierr);
	ierr = VecDestroy(&bc->bcT);  CHKERRQ(ierr);

	// SPC velocity-pressure
	ierr = PetscFree(bc->SPCList);  CHKERRQ(ierr);
	ierr = PetscFree(bc->SPCVals);  CHKERRQ(ierr);

	// SPC temperature
	ierr = PetscFree(bc->tSPCList); CHKERRQ(ierr);
	ierr = PetscFree(bc->tSPCVals); CHKERRQ(ierr);

	// two-point constraints
//	ierr = PetscFree(bc->TPCList);      CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCPrimeDOF);  CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCVals);      CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCLinComPar); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApply"
PetscErrorCode BCApply(BCCtx *bc, Vec x)
{
	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscScalar *SPCVals;
	PetscInt    i, ln, lnv, numSPC, vNumSPC, pNumSPC, *SPCList;

 	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs      = bc->fs;
	dof     = &fs->dof;
	ln      = dof->ln;
	lnv     = dof->lnv;
	SPCList = bc->SPCList;
	SPCVals = bc->SPCVals;

	// mark all variables unconstrained
	ierr = VecSet(bc->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcp,  DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(bc->bcT,  DBL_MAX); CHKERRQ(ierr);

	for(i = 0; i < ln; i++) SPCVals[i] = DBL_MAX;

	// apply boundary constraints
	ierr = BCApplyBound(bc); CHKERRQ(ierr);

	// apply Bezier block constraints
	ierr = BCApplyBezier(bc); CHKERRQ(ierr);

	// apply Boundary velocity
	ierr = BCApplyBoundVel(bc); CHKERRQ(ierr);

	// apply dropping boxes
	ierr = BCApplyDBox(bc); CHKERRQ(ierr);

	// apply simple shear BCs
	ierr = BCApplySimpleShear(bc); CHKERRQ(ierr);

	// exchange ghost point constraints
	// AVOID THIS BY SETTING CONSTRAINTS REDUNDANTLY
	// IN MULTIGRID ONLY REPEAT BC COARSENING WHEN THINGS CHANGE
	LOCAL_TO_LOCAL(fs->DA_X,   bc->bcvx)
	LOCAL_TO_LOCAL(fs->DA_Y,   bc->bcvy)
	LOCAL_TO_LOCAL(fs->DA_Z,   bc->bcvz)
	LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcp)
	LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcT)

	// form constraint lists
	numSPC = 0;

	// velocity
	for(i = 0; i < lnv; i++)
	{
		if(SPCVals[i] != DBL_MAX)
		{
			SPCList[numSPC] = i;
			SPCVals[numSPC] = SPCVals[i];
			numSPC++;
		}
	}
	vNumSPC = numSPC;

	// pressure
	for(i = lnv; i < ln; i++)
	{
		if(SPCVals[i] != DBL_MAX)
		{
			SPCList[numSPC] = i;
			SPCVals[numSPC] = SPCVals[i];
			numSPC++;
		}
	}
	pNumSPC = numSPC - vNumSPC;

	// set index (shift) type
	bc->stype = _GLOBAL_TO_LOCAL_;

	// store constraint lists
	bc->numSPC   = numSPC;
	bc->vNumSPC  = vNumSPC;
	bc->vSPCList = SPCList;
	bc->vSPCVals = SPCVals;
	bc->pNumSPC  = pNumSPC;
	bc->pSPCList = SPCList + vNumSPC;
	bc->pSPCVals = SPCVals + vNumSPC;

	// WARNING! currently primary temperature constraints are not implemented
	bc->tNumSPC = 0;

	// apply SPC to global solution vector
	ierr = BCApplySPC(bc, x); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplySPC"
PetscErrorCode BCApplySPC(BCCtx *bc, Vec x)
{
	// apply SPC to global solution vector

	PetscScalar *sol, *vals;
	PetscInt    i, num, *list;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecGetArray(x, &sol); CHKERRQ(ierr);

	//============================================
	// enforce single point constraints (velocity)
	//============================================

	num   = bc->vNumSPC;
	list  = bc->vSPCList;
	vals  = bc->vSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	//============================================
	// enforce single point constraints (pressure)
	//============================================

	num   = bc->pNumSPC;
	list  = bc->pSPCList;
	vals  = bc->pSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	ierr = VecRestoreArray(x, &sol); CHKERRQ(ierr);

 	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCShiftIndices"
PetscErrorCode BCShiftIndices(BCCtx *bc, ShiftType stype)
{
	FDSTAG   *fs;
	DOFIndex *dof;
	PetscInt i, vShift, pShift;

	PetscInt vNumSPC, pNumSPC, *vSPCList, *pSPCList;

	// error checking
	if(stype == bc->stype)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Cannot call same type of index shifting twice in a row");
	}

	// access context
	fs       = bc->fs;
	dof      = &fs->dof;
	vNumSPC  = bc->vNumSPC;
	vSPCList = bc->vSPCList;
	pNumSPC  = bc->pNumSPC;
	pSPCList = bc->pSPCList;

	// get local-to-global index shifts
	if(dof->idxmod == IDXCOUPLED)   { vShift = dof->st;  pShift = dof->st;             }
	if(dof->idxmod == IDXUNCOUPLED) { vShift = dof->stv; pShift = dof->stp - dof->lnv; }

	// shift constraint indices
	if(stype == _LOCAL_TO_GLOBAL_)
	{
		for(i = 0; i < vNumSPC; i++) vSPCList[i] += vShift;
		for(i = 0; i < pNumSPC; i++) pSPCList[i] += pShift;
	}
	else if(stype == _GLOBAL_TO_LOCAL_)
	{
		for(i = 0; i < vNumSPC; i++) vSPCList[i] -= vShift;
		for(i = 0; i < pNumSPC; i++) pSPCList[i] -= pShift;
	}

	// switch shit type
	bc->stype = stype;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Specific constraints
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyBound"
PetscErrorCode BCApplyBound(BCCtx *bc)
{
	// initialize boundary conditions vectors

	// *************************************************************
	//    NO-SLIP OR FREE-SLIP TANGENTIAL BOUNDARY VELOCITIES CAN BE SELECTED
	//    NORMAL VELOCITIES CAN BE SET VIA BACKGROUND STRAIN RATES
	//    PRESSURE CONSTRAINS ARE CURRENTLY NOT ALLOWED
	//    TEMPERATURE IS PRESCRIBED ON TOP AND BOTTOM BOUNDARIES
	//    ONLY ZERO BOUNDARY HEAT FLUXES ARE IMPLEMENTED
	//    TWO-POINT CONSTRAINTS MUST BE SET ON CROSS-PROCESSOR GHOST POINTS
	// *************************************************************

	FDSTAG      *fs;
	PetscScalar Tbot, Ttop;
	PetscScalar Exx, Eyy, Ezz;
	PetscScalar bx,  by,  bz;
	PetscScalar ex,  ey,  ez;
	PetscScalar vbx, vby, vbz;
	PetscScalar vex, vey, vez;
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter, top_open;
	PetscInt    nsLeft, nsRight, nsFront, nsBack, nsBottom, nsTop;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcT, *SPCVals;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = bc->fs;

	// set open boundary flag
	top_open = bc->top_open;

	// set no-slip flags
	nsLeft   		= bc->noslip[0];
	nsRight  		= bc->noslip[1];
	nsFront  		= bc->noslip[2];
	nsBack   		= bc->noslip[3];
	nsBottom 		= bc->noslip[4];
	nsTop    		= bc->noslip[5];

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// get current coordinates of the mesh boundaries
	ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	// get background strain rates
	ierr = BCGetBGStrainRates(bc, &Exx, &Eyy, &Ezz); CHKERRQ(ierr);

	// get boundary velocities
	// coordinate origin is assumed to be fixed
	// velocity is a product of strain rate and coordinate
	vbx = bx*Exx;   vex = ex*Exx;
	vby = by*Eyy;   vey = ey*Eyy;
	vbz = bz*Ezz;   vez = ez*Ezz;

	if(top_open)
	{
		vbz = 0.0;
		vez = 0.0;
	}

	// get boundary temperatures
	Tbot = bc->Tbot;
	Ttop = bc->Ttop;

	// access constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcT,  &bcT);  CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	//=========================================================================
	// SPC (normal velocities)
	//=========================================================================

	iter = 0;

	//------------------
	// X points SPC only
	//------------------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(i == 0)   { bcvx[k][j][i] = vbx; SPCVals[iter] = vbx; }
		if(i == mnx) { bcvx[k][j][i] = vex; SPCVals[iter] = vex; }
		iter++;
	}
	END_STD_LOOP

	//------------------
	// Y points SPC only
	//------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(j == 0)   { bcvy[k][j][i] = vby; SPCVals[iter] = vby; }
		if(j == mny) { bcvy[k][j][i] = vey; SPCVals[iter] = vey; }
		iter++;
	}
	END_STD_LOOP

	//------------------
	// Z points SPC only
	//------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(k == 0)                { bcvz[k][j][i] = vbz; SPCVals[iter] = vbz; }
		if(k == mnz && !top_open) { bcvz[k][j][i] = vez; SPCVals[iter] = vez; }
		iter++;
	}
	END_STD_LOOP

	//=========================================================================
	// TPC (no-slip boundary conditions)
	//=========================================================================

	//-----------------------------------------------------
	// X points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsFront || nsBack || nsBottom || nsTop)
	{
		GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(nsFront  && j == 0)   { bcvx[k][j-1][i] = 0.0; }
			if(nsBack   && j == mcy) { bcvx[k][j+1][i] = 0.0; }
			if(nsBottom && k == 0)   { bcvx[k-1][j][i] = 0.0; }
			if(nsTop    && k == mcz) { bcvx[k+1][j][i] = 0.0; }
		}
		END_STD_LOOP
	}

	//-----------------------------------------------------
	// Y points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsLeft || nsRight || nsBottom || nsTop)
	{

		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_NODE_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(nsLeft   && i == 0)   { bcvy[k][j][i-1] = 0.0; }
			if(nsRight  && i == mcx) { bcvy[k][j][i+1] = 0.0; }
			if(nsBottom && k == 0)   { bcvy[k-1][j][i] = 0.0; }
			if(nsTop    && k == mcz) { bcvy[k+1][j][i] = 0.0; }
		}
		END_STD_LOOP
	}

	//-----------------------------------------------------
	// Z points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsLeft || nsRight || nsFront || nsBack)
	{

		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_NODE_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(nsLeft  && i == 0)   { bcvz[k][j][i-1] = 0.0; }
			if(nsRight && i == mcx) { bcvz[k][j][i+1] = 0.0; }
			if(nsFront && j == 0)   { bcvz[k][j-1][i] = 0.0; }
			if(nsBack  && j == mcy) { bcvz[k][j+1][i] = 0.0; }
		}
		END_STD_LOOP
	}

	//-----------------------------------------------------
	// T points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(Tbot >= 0.0 || Ttop >= 0.0)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// only positive temperature!
			// negative will set zero-flux BC automatically
			if(Tbot >= 0.0 && k == 0)   { bcT[k-1][j][i] = Tbot; }
			if(Ttop >= 0.0 && k == mcz) { bcT[k+1][j][i] = Ttop; }
		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcT,  &bcT);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyBezier"
PetscErrorCode BCApplyBezier(BCCtx *bc)
{
	FDSTAG      *fs;
	BCBlock     *bcb;
	PetscInt    fbeg, fend, npoly, in;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter, ib;
	PetscScalar ***bcvx,  ***bcvy, *SPCVals;
	PetscScalar t, dt, theta, costh, sinth, atol, bot, top, vel;
	PetscScalar Xbeg[3], Xend[3], xbeg[3], xend[3], box[4], cpoly[2*_max_poly_points_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether constraint is activated
	if(!bc->nblocks) PetscFunctionReturn(0);

	// access context
	fs    =  bc->fs;
	t     =  bc->ts->time;
	dt    =  bc->ts->dt;

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	// loop over all bezier blocks
	for(ib = 0; ib < bc->nblocks; ib++)
	{
		bcb   =  bc->blocks + ib;
		bot   =  bcb->bot;
		top   =  bcb->top;
		npoly =  bcb->npoly;

		// get polygon positions in the beginning & end of the time step
		ierr = BCBlockGetPosition(bcb, t,    &fbeg, Xbeg); CHKERRQ(ierr);
		ierr = BCBlockGetPosition(bcb, t+dt, &fend, Xend); CHKERRQ(ierr);

		// check whether constraint applies to the current time step
		if(!fbeg || !fend) continue;

		// get current polygon geometry
		ierr = BCBlockGetPolygon(bcb, Xbeg, cpoly);

		// get bounding box
		polygon_box(&npoly, cpoly, 1e-12, &atol, box);

		// get time step rotation matrix
		theta = Xend[2] - Xbeg[2];
		costh = cos(theta);
		sinth = sin(theta);

		iter = 0;

		//---------
		// X points
		//---------
		GET_NODE_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// get node coordinates in the beginning of time step
			xbeg[0] = COORD_NODE(i, sx, fs->dsx);
			xbeg[1] = COORD_CELL(j, sy, fs->dsy);
			xbeg[2] = COORD_CELL(k, sz, fs->dsz);

			// perform point test
			if(xbeg[2] >= bot && xbeg[2] <= top)
			{
				in_polygon(1, xbeg, npoly, cpoly, box, atol, &in);

				// check whether point is inside polygon
				if(in)
				{
					// compute point position in the end of time step
					RotDispPoint2D(Xbeg, Xend, costh, sinth, xbeg, xend);

					// compute & set x-velocity
					vel = (xend[0] - xbeg[0])/dt;

					bcvx[k][j][i] = vel;
					SPCVals[iter] = vel;
				}
			}
			iter++;
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
			// get node coordinates in the beginning of time step
			xbeg[0] = COORD_CELL(i, sx, fs->dsx);
			xbeg[1] = COORD_NODE(j, sy, fs->dsy);
			xbeg[2] = COORD_CELL(k, sz, fs->dsz);

			// perform point test
			if(xbeg[2] >= bot && xbeg[2] <= top)
			{
				in_polygon(1, xbeg, npoly, cpoly, box, atol, &in);

				// check whether point is inside polygon
				if(in)
				{
					// compute point position in the end of time step
					RotDispPoint2D(Xbeg, Xend, costh, sinth, xbeg, xend);

					// compute & set y-velocity
					vel = (xend[1] - xbeg[1])/dt;

					bcvy[k][j][i] = vel;
					SPCVals[iter] = vel;
				}
			}
			iter++;
		}
		END_STD_LOOP
	}
	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyBoundVel"
PetscErrorCode BCApplyBoundVel(BCCtx *bc)
{
	FDSTAG      *fs;
	PetscInt    mnx, mny;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar ***bcvx,  ***bcvy, *SPCVals;
	PetscScalar z, bot, top, vel, velin, velout;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether constraint is activated
	if(!bc->face) PetscFunctionReturn(0);

	// access context
	fs     = bc->fs;
	bot    = bc->bot;
	top    = bc->top;
	velin  = bc->velin;
	velout = bc->velout;

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	iter = 0;

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	if(bc->face == 1 || bc->face == 2)
	{
		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			vel = 0.0;
			if(z <= top && z >= bot) vel = velin;
			if(z < bot)              vel = velout;

			if(bc->face == 1 && i == 0)   { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
			if(bc->face == 2 && i == mnx) { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
			iter++;
		}
		END_STD_LOOP
	}

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	if(bc->face == 3 || bc->face == 4)
	{
		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			vel = 0.0;
			if(z <= top && z >= bot) vel = velin;
			if(z < bot)              vel = velout;

			if(bc->face == 3 && j == 0)   { bcvy[k][j][i] = vel; SPCVals[iter] = vel; }
			if(bc->face == 4 && j == mny) { bcvy[k][j][i] = vel; SPCVals[iter] = vel; }
			iter++;
		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyDBox"
PetscErrorCode BCApplyDBox(BCCtx *bc)
{
	DBox        *dbox;
	FDSTAG      *fs;
	PetscInt    jj, i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar ***bcvz, *SPCVals, bounds[6*_max_boxes_], *pbounds;
	PetscScalar x, y, z, t, vz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs   = bc->fs;
	dbox = &bc->dbox;

	// check whether dropping box is activated
	if(!dbox->num) PetscFunctionReturn(0);

	// copy original coordinates
	ierr = PetscMemcpy(bounds, dbox->bounds, (size_t)(6*dbox->num)*sizeof(PetscScalar)); CHKERRQ(ierr);

	// integrate box positions
	t  = bc->ts->time;
	vz = dbox->zvel;

	for(jj = 0; jj < dbox->num; jj++)
	{
		pbounds     = bounds + 6*jj;
		pbounds[4] += t*vz;
		pbounds[5] += t*vz;
	}

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz); CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	//---------
	// Z points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	iter = fs->nXFace + fs->nYFace;

	START_STD_LOOP
	{
		// get node coordinates
		x = COORD_CELL(i, sx, fs->dsx);
		y = COORD_CELL(j, sy, fs->dsy);
		z = COORD_NODE(k, sz, fs->dsz);

		// check whether node is inside any of boxes
		for(jj = 0; jj < dbox->num; jj++)
		{
			pbounds = bounds + 6*jj;

			if(x >= pbounds[0] && x <= pbounds[1]
			&& y >= pbounds[2] && y <= pbounds[3]
			&& z >= pbounds[4] && z <= pbounds[5])
			{
				bcvz[k][j][i] = vz;
				SPCVals[iter] = vz;
				break;
			}
		}
		iter++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplySimpleShear"
PetscErrorCode BCApplySimpleShear(BCCtx *bc)
{
	FDSTAG      *fs;
	PetscInt    mnx, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar ***bcvx, *SPCVals;
	PetscScalar z, gamma_xz, vel;

	// This routine sets simple shear velocities with constant shear rate gamma.
	// This will imply inflow/outflow

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether constraint is activated
	if(bc->gamma_xz == 0.0) PetscFunctionReturn(0);

	// access context
	fs     		= bc->fs;
	gamma_xz    = bc->gamma_xz;

	// initialize maximal index
	mnx = fs->dsx.tnods - 1;
	mcz = fs->dsz.tcels - 1;

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	iter = 0;

	//----------------------------------------
	// normal x - velocities (SPC constraints)
	//----------------------------------------

	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		z   = COORD_CELL(k, sz, fs->dsz);
		vel = gamma_xz*z;

		if(i == 0)   { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
		if(i == mnx) { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
		iter++;
	}
	END_STD_LOOP

	//--------------------------------------------
	// tangential x - velocities (TPC constraints)
	//--------------------------------------------

	GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		z   = COORD_CELL(k, sz, fs->dsz);
		vel = gamma_xz*z;

		if(k == 0)   { bcvx[k-1][j][i] = vel; }
		if(k == mcz) { bcvx[k+1][j][i] = vel; }
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Service functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCGetBGStrainRates"
PetscErrorCode BCGetBGStrainRates(BCCtx *bc, PetscScalar *Exx_, PetscScalar *Eyy_, PetscScalar *Ezz_)
{
	// get current background strain rates

	PetscInt    jj;
	PetscScalar time, Exx, Eyy, Ezz;

	// initialize
	time = bc->ts->time;
	Exx  = 0.0;
	Eyy  = 0.0;
	Ezz  = 0.0;

	// x-direction background strain rate
	if(bc->ExxNumPeriods)
	{
		for(jj = 0; jj < bc->ExxNumPeriods-1; jj++)
		{
			if(time < bc->ExxTimeDelims[jj]) break;
		}

		Exx = bc->ExxStrainRates[jj];
	}

	// y-direction background strain rate
	if(bc->EyyNumPeriods)
	{
		for(jj = 0; jj < bc->EyyNumPeriods-1; jj++)
		{
			if(time < bc->EyyTimeDelims[jj]) break;
		}

		Eyy = bc->EyyStrainRates[jj];
	}

	// z-direction background strain rate
	Ezz = -(Exx+Eyy);

	// store result
	if(Exx_) (*Exx_) = Exx;
	if(Eyy_) (*Eyy_) = Eyy;
	if(Ezz_) (*Ezz_) = Ezz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCStretchGrid"
PetscErrorCode BCStretchGrid(BCCtx *bc)
{
	// apply background strain-rate "DWINDLAR" BC (Bob Shaw "Ship of Strangers")

	// Stretch grid with constant stretch factor about coordinate origin.
	// The origin point remains fixed, and the displacements of all points are
	// proportional to the distance from the origin (i.e. coordinate).
	// The origin (zero) point must remain within domain (checked at input).
	// Stretch factor is positive at extension, i.e.:
	// eps = (L_new-L_old)/L_old
	// L_new = L_old + eps*L_old
	// x_new = x_old + eps*x_old

	TSSol       *ts;
	FDSTAG      *fs;
	PetscScalar Exx, Eyy, Ezz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = bc->fs;
	ts = bc->ts;

	// get background strain rates
	ierr = BCGetBGStrainRates(bc, &Exx, &Eyy, &Ezz); CHKERRQ(ierr);

	// stretch grid
	if(Exx) { ierr = Discret1DStretch(&fs->dsx, Exx*ts->dt); CHKERRQ(ierr); }
	if(Eyy) { ierr = Discret1DStretch(&fs->dsy, Eyy*ts->dt); CHKERRQ(ierr); }
	if(Ezz) { ierr = Discret1DStretch(&fs->dsz, Ezz*ts->dt); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCOverridePhase"
PetscErrorCode BCOverridePhase(BCCtx *bc, PetscInt cellID, Marker *P)
{
	FDSTAG     *fs;
	PetscInt    i, j, k, M, N, mx, my, sx, sy;
	PetscScalar z;

	PetscFunctionBegin;

	fs = bc->fs;
	M  = fs->dsx.ncels;
	N  = fs->dsy.ncels;
	sx = fs->dsx.pstart;
	sy = fs->dsy.pstart;
	mx = fs->dsx.tcels-1;
	my = fs->dsy.tcels-1;
	z  = P->X[2];

	// expand i, j, k cell indices
	GET_CELL_IJK(cellID, i, j, k, M, N);

	if(((bc->face == 1 && i + sx == 0)
	||  (bc->face == 2 && i + sx == mx)
	||  (bc->face == 3 && j + sy == 0)
	||  (bc->face == 4 && j + sy == my))
	&&  (z >= bc->bot && z <= bc->top))
	{
		P->phase = bc->phase;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
