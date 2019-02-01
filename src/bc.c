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
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * replace BC input specification consistently in the entire code
// * open box & Winkler (with tangential viscous friction)
// * tangential velocities
// * extend two-point constraint specification & (possibly) get rid bc-vectors
// * create bc-object only at fine level, coarse levels should have simple access
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCBlockReadFromOptions"
PetscErrorCode BCBlockReadFromOptions(BCBlock *bcb, Scaling *scal)
{
	//	-bcb_npath - Number of path points of Bezier curve (end-points only!)
	//	-bcb_theta - Orientation angles at path points (counter-clockwise positive)
	//	-bcb_time  - Times at path points
	//	-bcb_path  - Bezier curve path & control points (6*npath-4 points are expected)
	//	-bcb_npoly - Number of polygon vertices
	//	-bcb_poly  - Polygon x-y coordinates at initial time
	//	-bcb_bot   - Polygon bottom coordinate
	//	-bcb_top   - Polygon top coordinate

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//==================
	// Bezier path curve
	//==================

	ierr = GetIntDataItemCheck("-bcb_npath", "Number of path points of Bezier curve",
		_NOT_FOUND_EXIT_, 1, &bcb->npath, 1, _max_path_points_); CHKERRQ(ierr);

	if(bcb->npath)
	{
		// angles
		ierr = GetScalDataItemCheckScale("-bcb_theta", "Orientation angles at path points",
			_NOT_FOUND_ERROR_, bcb->npath, bcb->theta, 0.0, 0.0, scal->angle); CHKERRQ(ierr);

		// times
		ierr = GetScalDataItemCheckScale("-bcb_time", "Times at path points",
			_NOT_FOUND_ERROR_, bcb->npath, bcb->time, 0.0, 0.0, scal->time); CHKERRQ(ierr);

		// path coordinates
		ierr = GetScalDataItemCheckScale("-bcb_path", "Bezier curve path & control points",
			_NOT_FOUND_ERROR_, 6*bcb->npath-4, bcb->path, 0.0, 0.0, scal->length); CHKERRQ(ierr);

		//========
		// polygon
		//========

		ierr = GetIntDataItemCheck("-bcb_npoly", "Number of polygon vertices",
			_NOT_FOUND_EXIT_, 1, &bcb->npoly, 1, _max_poly_points_); CHKERRQ(ierr);

		// polygon coordinates
		ierr = GetScalDataItemCheckScale("-bcb_poly", "Polygon coordinates",
			_NOT_FOUND_ERROR_, 2*bcb->npoly, bcb->poly, 0.0, 0.0, scal->length); CHKERRQ(ierr);

		// polygon bottom
		ierr = GetScalDataItemCheckScale("-bcb_bot", "Polygon bottom coordinate",
			_NOT_FOUND_ERROR_, 1, &bcb->bot, 0.0, 0.0, scal->length); CHKERRQ(ierr);

		// polygon top
		ierr = GetScalDataItemCheckScale("-bcb_top", "Polygon top coordinate",
			_NOT_FOUND_ERROR_, 1, &bcb->top, 0.0, 0.0, scal->length); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCBlockGetPosition"
PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *f, PetscScalar X[])
{
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

	// check whether bezier is applied
	if(bc->AddBezier != PETSC_TRUE) PetscFunctionReturn(0);

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
	for(ib = 0; ib < bc->nblo; ib++)
	{
		bcb   = &bc->blocks[ib];
		bot   =  bcb->bot;
		top   =  bcb->top;
		npoly =  bcb->npoly;

		// get polygon positions in the beginning & end of the time step
		ierr = BCBlockGetPosition(bcb, t,    &fbeg, Xbeg); CHKERRQ(ierr);
		ierr = BCBlockGetPosition(bcb, t+dt, &fend, Xend); CHKERRQ(ierr);

		// check whether constraint applies to the current time step
		if(!fbeg || !fend) continue;

		// print bezier block overview
		PetscPrintf(PETSC_COMM_WORLD,"Bezier block[%d]: npath = %d, npoly = %d\n",ib, bcb->npath, bcb->npoly);

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
#define __FUNCT__ "BCSetupBoundVel"
PetscErrorCode BCSetupBoundVel(BCCtx *bc, PetscScalar top)
{
	PetscScalar bz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set plate top
	bc->top = top;

	ierr = FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL); CHKERRQ(ierr);

	// compute outflow velocity
	bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);

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
#undef __FUNCT__
#define __FUNCT__ "DBoxReadFromOptions"
PetscErrorCode DBoxReadFromOptions(DBox *dbox, Scaling *scal)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	//========================
	// Dropping box parameters
	//========================

	ierr = GetIntDataItemCheck("-dbox_num", "Number of dropping boxes",
		_NOT_FOUND_EXIT_, 1, &dbox->num, 1, _max_boxes_); CHKERRQ(ierr);

	if(dbox->num)
	{
		// boundaries
		ierr = GetScalDataItemCheckScale("-dbox_bounds", "Dropping box boundaries",
			_NOT_FOUND_ERROR_, 6*dbox->num, dbox->bounds, 0.0, 0.0, scal->length); CHKERRQ(ierr);

		// vertical velocity
		ierr = GetScalDataItemCheckScale("-dbox_zvel", "Dropping box vertical velocity",
			_NOT_FOUND_ERROR_, 1, &dbox->zvel, 0.0, 0.0, scal->velocity); CHKERRQ(ierr);

	}

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
#define __FUNCT__ "BCApplySimpleShearVel"
PetscErrorCode BCApplySimpleShearVel(BCCtx *bc)
{
	FDSTAG      *fs;
	PetscInt    mnx, mny;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar ***bcvx,  ***bcvy, *SPCVals;
	PetscScalar z, gamma_xz, vel;

	/* This routine sets simple shear velocities, such that we have
	 * constant shear rate gamma.
	 * This will imply inflow/outflow
	 *
	 */

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether constraint is activated
	if(!bc->simpleshear) PetscFunctionReturn(0);

	// access context
	fs     		= bc->fs;
	gamma_xz    = bc->gamma_xz;

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

	START_STD_LOOP
	{
		z   = COORD_CELL(k, sz, fs->dsz);
		vel = gamma_xz*z;

		if(i == 0)   { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
		if(i == mnx) { bcvx[k][j][i] = vel; SPCVals[iter] = vel; }
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
		z   = COORD_CELL(k, sz, fs->dsz);
		vel = 0.0;
		if(j == 0)   { bcvy[k][j][i] = vel; SPCVals[iter] = vel; }
		if(j == mny) { bcvy[k][j][i] = vel; SPCVals[iter] = vel; }
		iter++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
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
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs, TSSol *ts, Scaling *scal)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create boundary conditions vectors (velocity, pressure, temperature)
	ierr = DMCreateLocalVector(fs->DA_X,   &bc->bcvx);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &bc->bcvy);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &bc->bcvz);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcp);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcT);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &bc->bcPl);  CHKERRQ(ierr); // liquidPressure/Darcy

	// SPC velocity-pressure
	ierr = makeIntArray (&bc->SPCList, NULL, fs->dof.ln);   CHKERRQ(ierr);
	ierr = makeScalArray(&bc->SPCVals, NULL, fs->dof.ln);   CHKERRQ(ierr);

	// SPC (temperature)
	ierr = makeIntArray (&bc->tSPCList, NULL, fs->dof.lnp); CHKERRQ(ierr);
	ierr = makeScalArray(&bc->tSPCVals, NULL, fs->dof.lnp); CHKERRQ(ierr);

	// SPC (LiquidPressure/Darcy)
	ierr = makeIntArray (&bc->Pl_SPCList, NULL, fs->dof.lnp); CHKERRQ(ierr);
	ierr = makeScalArray(&bc->Pl_SPCVals, NULL, fs->dof.lnp); CHKERRQ(ierr);

	bc->ExxAct = PETSC_FALSE;
	bc->EyyAct = PETSC_FALSE;
	bc->pbAct  = PETSC_FALSE;

	bc->fs   = fs;
	bc->ts   = ts;
	bc->scal = scal;

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
	ierr = VecDestroy(&bc->bcPl);  CHKERRQ(ierr);// liquidPressure/Darcy

	// SPC velocity-pressure
	ierr = PetscFree(bc->SPCList);  CHKERRQ(ierr);
	ierr = PetscFree(bc->SPCVals);  CHKERRQ(ierr);

	// SPC temperature
	ierr = PetscFree(bc->tSPCList); CHKERRQ(ierr);
	ierr = PetscFree(bc->tSPCVals); CHKERRQ(ierr);

	// SPC (LiquidPressure/Darcy)
	ierr = PetscFree(bc->Pl_SPCList); CHKERRQ(ierr);
	ierr = PetscFree(bc->Pl_SPCVals); CHKERRQ(ierr);

	// two-point constraints
//	ierr = PetscFree(bc->TPCList);      CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCPrimeDOF);  CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCVals);      CHKERRQ(ierr);
//	ierr = PetscFree(bc->TPCLinComPar); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCSetParam"
PetscErrorCode BCSetParam(BCCtx *bc, UserCtx *user)
{
	PetscFunctionBegin;

	bc->Tbot  = user->Temp_bottom;
	bc->Ttop  = user->Temp_top;

	if (user->Pl_bottom != 0.0 )
	{
		bc->Plbot  = user->Pl_bottom; //LiquidPressure/Darcy
	}
	bc->Pltop  = user->Pl_top;

	if(user->AddBezier)
	{
		bc->AddBezier = PETSC_TRUE;
		bc->blocks = user->blocks;
		bc->nblo   = user->nblo;
	}
	else
		bc->AddBezier = PETSC_FALSE;

	// Darcy
	//if(user->GradientStraintRate)
	//{
	//	bc->GradientStraintRate = PETSC_TRUE;
	//}
	//else
	//	bc->GradientStraintRate = PETSC_FALSE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCReadFromOptions"
PetscErrorCode BCReadFromOptions(BCCtx *bc)
{
	// set parameters from PETSc options

	PetscBool  set;
	Scaling   *scal;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = bc->scal;

	// x-direction background strain rate

	ierr = GetIntDataItemCheck("-ExxNumPeriods", "Number of Exx background strain rate periods",
		_NOT_FOUND_EXIT_, 1, &bc->ExxNumPeriods, 1, _max_periods_); CHKERRQ(ierr);

	if(bc->ExxNumPeriods)
	{
		bc->ExxAct = PETSC_TRUE;

		ierr = GetScalDataItemCheckScale("-ExxTimeDelims", "Exx background strain rate time delimiters",
			_NOT_FOUND_ERROR_, bc->ExxNumPeriods-1, bc->ExxTimeDelims, 0.0, 0.0, scal->time); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-ExxStrainRates", "Exx background strain rates",
			_NOT_FOUND_ERROR_, bc->ExxNumPeriods, bc->ExxStrainRates, 0.0, 0.0, scal->strain_rate); CHKERRQ(ierr);
	}

	// y-direction background strain rate

	ierr = GetIntDataItemCheck("-EyyNumPeriods", "Number of Eyy background strain rate periods",
		_NOT_FOUND_EXIT_, 1, &bc->EyyNumPeriods, 1, _max_periods_); CHKERRQ(ierr);

	if(bc->EyyNumPeriods)
	{
		bc->EyyAct = PETSC_TRUE;

		ierr = GetScalDataItemCheckScale("-EyyTimeDelims", "Eyy background strain rate time delimiters",
			_NOT_FOUND_ERROR_, bc->EyyNumPeriods-1, bc->EyyTimeDelims, 0.0, 0.0, scal->time); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-EyyStrainRates", "Eyy background strain rates",
			_NOT_FOUND_ERROR_, bc->EyyNumPeriods, bc->EyyStrainRates, 0.0, 0.0, scal->strain_rate); CHKERRQ(ierr);
	}

	// Bezier block
	//ierr = BCBlockReadFromOptions(&bc->blocks, scal); CHKERRQ(ierr);


	if(bc->AddBezier && (bc->ExxAct == PETSC_TRUE || bc->EyyAct == PETSC_TRUE))
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot combine background strain rate with moving block\n");
	}

	// boundary velocities
	ierr = GetIntDataItemCheck("-bvel_face", "Boundary velocity face identifier",
		_NOT_FOUND_EXIT_, 1, &bc->face, 1, 4); CHKERRQ(ierr);

	if(bc->face)
	{
		if(bc->face && (bc->AddBezier || bc->ExxAct == PETSC_TRUE || bc->EyyAct == PETSC_TRUE))
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot combine boundary velocity with either background strain rate or moving block\n");
		}

		ierr = GetIntDataItemCheck("-bvel_phase", "Boundary velocity phase",
			_NOT_FOUND_EXIT_, 1, &bc->phase, 0, 0); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-bvel_bot", "Boundary velocity bottom level",
			_NOT_FOUND_ERROR_, 1, &bc->bot, 0.0, 0.0, scal->length); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-bvel_velin", "Boundary velocity magnitude",
			_NOT_FOUND_ERROR_, 1, &bc->velin, 0.0, 0.0, scal->velocity); CHKERRQ(ierr);
	}

	// set open boundary flag
	ierr = PetscOptionsHasName(NULL, NULL, "-open_top_bound", &set); CHKERRQ(ierr); if(set == PETSC_TRUE) bc->top_open = 1;

	// set simple shear boundary condition
	ierr = PetscOptionsHasName(NULL, NULL, "-bc_simpleshear", &set); CHKERRQ(ierr); if(set == PETSC_TRUE) bc->simpleshear = 1;
	if(bc->simpleshear)
	{
		ierr = GetScalDataItemCheckScale("-bc_simpleshear_gamma_xz", "Simple shear strain rate",
					_NOT_FOUND_ERROR_, 1, &bc->gamma_xz, 0.0, 0.0, scal->time); CHKERRQ(ierr);
	}


	ierr = DBoxReadFromOptions(&bc->dbox, scal); CHKERRQ(ierr);

	// read no-slip boundary condition mask
	ierr = GetIntDataItemCheck("-noslip", "no-slip mask",
		_NOT_FOUND_EXIT_, 6, bc->noslip, 0, 1); CHKERRQ(ierr);

	// Darcy
	ierr = PetscOptionsHasName(NULL, NULL, "-gradientStraintRate", &set); CHKERRQ(ierr); if(set == PETSC_TRUE) bc->GradientStraintRate = PETSC_TRUE;

	PetscFunctionReturn(0);
}
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
	if(bc->ExxAct == PETSC_TRUE)
	{
		for(jj = 0; jj < bc->ExxNumPeriods-1; jj++)
		{
			if(time < bc->ExxTimeDelims[jj]) break;
		}

		Exx = bc->ExxStrainRates[jj];
	}

	// y-direction background strain rate
	if(bc->EyyAct == PETSC_TRUE)
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
	ierr = VecSet(bc->bcPl, DBL_MAX); CHKERRQ(ierr); //LiquidPressure/Darcy

	for(i = 0; i < ln; i++) SPCVals[i] = DBL_MAX;

	// apply boundary constraints
	ierr = BCApplyBound(bc); CHKERRQ(ierr);

	// compute pushing parameters
	//ierr = BCCompPush(bc); CHKERRQ(ierr);

	// apply pushing block constraints
	ierr = BCApplyPush(bc); CHKERRQ(ierr);

	// apply Bezier block constraints
	ierr = BCApplyBezier(bc); CHKERRQ(ierr);

	// apply Boundary velocity
	ierr = BCApplyBoundVel(bc); CHKERRQ(ierr);

	// apply dropping boxes
	ierr = BCApplyDBox(bc); CHKERRQ(ierr);

	// apply simple shear BCs
	ierr = BCApplySimpleShearVel(bc); CHKERRQ(ierr);


	// exchange ghost point constraints
	// AVOID THIS BY SETTING CONSTRAINTS REDUNDANTLY
	// IN MULTIGRID ONLY REPEAT BC COARSENING WHEN THINGS CHANGE
	LOCAL_TO_LOCAL(fs->DA_X,   bc->bcvx)
	LOCAL_TO_LOCAL(fs->DA_Y,   bc->bcvy)
	LOCAL_TO_LOCAL(fs->DA_Z,   bc->bcvz)
	LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcp)
	LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcT)
	LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcPl)

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

	// LiquidPressure/Darcy
	bc->Pl_NumSPC = 0;

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
	PetscScalar Plbot, Pltop;
	PetscScalar Exx, Eyy, Ezz;
	PetscScalar bx,  by,  bz;
	PetscScalar ex,  ey,  ez;
	PetscScalar vbx, vby, vbz;
	PetscScalar vex, vey, vez;
	PetscScalar gamma_XZ, z;
	PetscInt    mcx, mcy, mcz;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter, top_open;
	PetscInt    nsLeft, nsRight, nsFront, nsBack, nsBottom, nsTop;
	PetscInt 	simpleShear;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcT, ***bcPl, *SPCVals;
	PetscScalar fac;

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
	simpleShear	 	= bc->simpleshear;

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

	// get boundary liquid pressures/Darcy
	Plbot 	=	bc->Plbot;
	Pltop 	=	bc->Pltop;

	//Pltop = (0.6+1.0)*1000.0*(-10)*20000.0;

	// access constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcT,   &bcT);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcPl,  &bcPl); CHKERRQ(ierr); // LiquidPressure/Darcy

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
		//if(i == 0)   { bcvx[k][j][i] = fac*vbx; SPCVals[iter] = fac*vbx; }
		//if(i == mnx) { bcvx[k][j][i] = fac*vex; SPCVals[iter] = fac*vex; }

		// Darcy
		fac = 1.0;
		if (bc->GradientStraintRate == PETSC_TRUE) fac=(double)(mcz-k)/(double)mcz;
		if(i == 0)   { bcvx[k][j][i] = fac*vbx; SPCVals[iter] = fac*vbx; }
		if(i == mnx) { bcvx[k][j][i] = fac*vex; SPCVals[iter] = fac*vex; }
		///////////////////////////////////////////////////////////////////////////

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
		//if(j == 0)   { bcvy[k][j][i] = fac*vby; SPCVals[iter] = fac*vby; }
		//if(j == mny) { bcvy[k][j][i] = fac*vey; SPCVals[iter] = fac*vey; }

		// Darcy
		fac = 1.0;
		if (bc->GradientStraintRate == PETSC_TRUE) fac=(double)(mcz-k)/(double)mcz;
		if(j == 0)   { bcvy[k][j][i] = fac*vby; SPCVals[iter] = fac*vby; }
		if(j == mny) { bcvy[k][j][i] = fac*vey; SPCVals[iter] = fac*vey; }
		///////////////////////////////////////////////////////////////////////////

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
		//if(k == 0)                { bcvz[k][j][i] = vbz; SPCVals[iter] = vbz; }
		//if(k == mnz && !top_open) { bcvz[k][j][i] = vez; SPCVals[iter] = vez; }

		// Darcy
		fac = 1.0;
		if (bc->GradientStraintRate == PETSC_TRUE) fac=(double)(mcz-k)/(double)mcz;
		if(k == 0)                                                              { bcvz[k][j][i] = fac*vbz; SPCVals[iter] = fac*vbz; }
		if(k == mnz  && (!top_open || bc->GradientStraintRate == PETSC_TRUE) )  { bcvz[k][j][i] = fac*vez; SPCVals[iter] = fac*vez; }
		///////////////////////////////////////////////////////////////////////////

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
	if (simpleShear)
	{
		GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		gamma_XZ =  bc->gamma_xz;

		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			if(k == 0)   { bcvx[k-1][j][i] = z*gamma_XZ; }
			if(k == mcz) { bcvx[k+1][j][i] = z*gamma_XZ; }
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

	//-----------------------------------------------------
	// LiquidPressure/Darcy points
	//-----------------------------------------------------
	if(Plbot != 0.0 || Pltop != 0.0)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(k == 0)   {
				bcPl[k-1][j][i] = Plbot;    // bottom boundary, dirichlet value constant
			}
			if(k == mcz) {
				bcPl[k+1][j][i] = Pltop;
			}
		}
		END_STD_LOOP
	}


	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcT,  &bcT);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcPl, &bcPl); CHKERRQ(ierr); //LiquidPressure/Darcy

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCSetPush"
PetscErrorCode BCSetPush(BCCtx *bc, UserCtx *user)
{
	PetscInt ip;

	PetscFunctionBegin;

	if(user->AddPushing)
	{
		bc->pbAct = PETSC_TRUE;
		bc->pb    = user->Pushing;
		bc->nPblo = user->nPush;
		for(ip = 0; ip < bc->nPblo; ip++) bc->pbApp[ip] = 1;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCompPush"
PetscErrorCode BCCompPush(BCCtx *bc,PetscInt ip)
{
	// MUST be called at the beginning of time step before setting boundary conditions
	// compute pushing boundary conditions actual parameters
	PushParams    *pb;
	TSSol         *ts;
	Scaling       *scal;
	PetscInt      i, ichange;
	PetscScalar   Vx, Vy, theta;

	PetscFunctionBegin;

	// check if pushing option is activated
	if(!bc->pbApp[ip]) PetscFunctionReturn(0);

	// access
	pb   = &bc->pb[ip];

	// access contexts
	ts   = bc->ts;
	scal = bc->scal;

	// set boundary conditions flag
	bc->pbApp[ip] = 0;

	// add pushing boundary conditions ONLY within the specified time interval - for that introduce a new flag
	if(ts->time >= pb->time[0]
							&& ts->time <= pb->time[pb->num_changes])
	{
		// check which pushing stage
		ichange = 0;

		for(i = 0; i < pb->num_changes; i++)
		{
			if(ts->time >= pb->time[i] && ts->time <= pb->time[i+1])
			{
				ichange = i;
			}
		}

		pb->ind_change = ichange;

		// initialize parameters for the time step
		Vx    = 0.0;
		Vy    = 0.0;
		theta = pb->theta;

		if(pb->dir[ichange] == 0)
		{
			Vx = cos(theta)*pb->V_push[ichange];
			Vy = sin(theta)*pb->V_push[ichange];
		}

		if(pb->dir[ichange] == 1) Vx = pb->V_push[ichange];
		if(pb->dir[ichange] == 2) Vy = pb->V_push[ichange];

		// set boundary conditions parameters
		bc->pbApp[ip] = 1;
		bc->Vx     	  = Vx;
		bc->Vy     	  = Vy;
		bc->theta     = theta;

		PetscPrintf(PETSC_COMM_WORLD,"Pushing Block[%d]: Time interval=[%g-%g] %s, V_push=%g %s, Omega=%g %s, mobile_block=%lld, direction=%lld, theta=%g deg\n",ip,
				pb->time             [ichange  ]*scal->time,
				pb->time             [ichange+1]*scal->time,             scal->lbl_time,
				pb->V_push           [ichange  ]*scal->velocity,         scal->lbl_velocity,
				pb->omega            [ichange  ]*scal->angular_velocity, scal->lbl_angular_velocity,
				(LLD)pb->coord_advect[ichange  ],
				(LLD)pb->dir         [ichange  ],
				theta*scal->angle);

		PetscPrintf(PETSC_COMM_WORLD,"\t\t Block center coordinates [x, y, z]=[%g, %g, %g] %s\n",
				pb->x_center_block*scal->length,
				pb->y_center_block*scal->length,
				pb->z_center_block*scal->length, scal->lbl_length);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyPush"
PetscErrorCode BCApplyPush(BCCtx *bc)
{
	// initialize internal (pushing) boundary conditions vectors
	// only x, y velocities are constrained!!
	// constraining vz makes a bad case - will not converge.

	FDSTAG      *fs;
	PushParams  *pb;
	PetscScalar xc, yc, zc;
	PetscScalar	dx, dy, dz;
	PetscScalar px, py, pz;
	PetscScalar rx, ry, rz;
	PetscScalar costh, sinth;
	PetscScalar ***bcvx,  ***bcvy, *SPCVals;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter, ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether pushing is applied
	if(bc->pbAct != PETSC_TRUE) PetscFunctionReturn(0);

	// access
	fs    = bc->fs;

	// access velocity constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy); CHKERRQ(ierr);

	// access constraint arrays
	SPCVals = bc->SPCVals;

	// loop over for pushing block
	for(ip = 0; ip < bc->nPblo; ip++)
	{
		// compute pushing parameters
		pb   = &bc->pb[ip];
		ierr = BCCompPush(bc, ip); CHKERRQ(ierr);

		// prepare block coordinates, sizes & rotation angle parameters
		xc    = pb->x_center_block;
		yc    = pb->y_center_block;
		zc    = pb->z_center_block;
		dx    = pb->L_block/2.0;
		dy    = pb->W_block/2.0;
		dz    = pb->H_block/2.0;
		costh = cos(bc->theta);
		sinth = sin(bc->theta);

		iter = 0;

		//---------
		// X points
		//---------
		GET_NODE_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// get point coordinates in the block-centered system
			px = COORD_NODE(i, sx, fs->dsx) - xc;
			py = COORD_CELL(j, sy, fs->dsy) - yc;
			pz = COORD_CELL(k, sz, fs->dsz) - zc;

			// get point coordinates in the block-aligned system
			// rotation matrix R = [ cos() sin() ; -sin() cos() ]
			rx =  costh*px + sinth*py;
			ry = -sinth*px + costh*py;
			rz =  pz;

			// perform point test
			if(rx >= -dx && rx <= dx
					&& ry >= -dy && ry <= dy
					&& rz >= -dz && rz <= dz)
			{
				bcvx[k][j][i] = bc->Vx;
				SPCVals[iter] = bc->Vx;
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
			// get point coordinates in the block-centered system
			px = COORD_CELL(i, sx, fs->dsx) - xc;
			py = COORD_NODE(j, sy, fs->dsy) - yc;
			pz = COORD_CELL(k, sz, fs->dsz) - zc;

			// get point coordinates in the block-aligned system
			// rotation matrix R = [ cos() sin() ; -sin() cos() ]
			rx =  costh*px + sinth*py;
			ry = -sinth*px + costh*py;
			rz =  pz;

			// perform point test
			if(rx >= -dx && rx <= dx
					&& ry >= -dy && ry <= dy
					&& rz >= -dz && rz <= dz)
			{
				bcvy[k][j][i] = bc->Vy;
				SPCVals[iter] = bc->Vy;
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
#define __FUNCT__ "BCAdvectPush"
PetscErrorCode BCAdvectPush(BCCtx *bc)
{
	TSSol        *ts;
	PushParams   *pb;
	PetscInt     advc, ichange, ip;
	PetscScalar  xc, yc, zc, Vx, Vy, dx, dy, dt, omega, theta, dtheta;

	// check if pushing option is activated
	if(bc->pbAct != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	ts = bc->ts;

	for(ip = 0; ip < bc->nPblo; ip++)
	{
		pb = &bc->pb[ip];

		// check time interval
		if(ts->time >= pb->time[0]
								&& ts->time <= pb->time[pb->num_changes])
		{
			// initialize variables
			ichange = pb->ind_change;
			advc    = pb->coord_advect[ichange];

			Vx = 0.0;
			Vy = 0.0;
			dy = 0.0;
			dx = 0.0;
			dt = ts->dt;

			xc = pb->x_center_block;
			yc = pb->y_center_block;
			zc = pb->z_center_block;

			// block is rotated
			if(pb->dir[ichange] == 0)
			{
				omega = pb->omega[ichange];
				theta = pb->theta;
				Vx    = cos(theta)*pb->V_push[ichange];
				Vy    = sin(theta)*pb->V_push[ichange];

				// rotation
				dtheta = omega*dt;
				pb->theta = theta + dtheta;
			}

			// block moves in X-direction
			if(pb->dir[ichange] == 1) Vx = pb->V_push[ichange];

			// block moves in Y-direction
			if(pb->dir[ichange] == 2) Vy = pb->V_push[ichange];

			// get displacements
			dx = Vx*dt;
			dy = Vy*dt;

			if(advc == 0)
			{	// stationary block
				pb->x_center_block = xc;
				pb->y_center_block = yc;
				pb->z_center_block = zc;
			}
			else
			{	// moving block
				pb->x_center_block = xc + dx;
				pb->y_center_block = yc + dy;
				pb->z_center_block = zc;
			}
		}
	}

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
	if(Exx) { ierr = Discret1DStretch(&fs->dsx, &fs->msx, Exx*ts->dt); CHKERRQ(ierr); }
	if(Eyy) { ierr = Discret1DStretch(&fs->dsy, &fs->msy, Eyy*ts->dt); CHKERRQ(ierr); }
	if(Ezz) { ierr = Discret1DStretch(&fs->dsz, &fs->msz, Ezz*ts->dt); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
