/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This sofware was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   surf.c
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
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "Utils.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * stair-case type of free surface
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfClear"
PetscErrorCode FreeSurfClear(FreeSurf *surf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(surf, sizeof(FreeSurf)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfCreate"
PetscErrorCode FreeSurfCreate(FreeSurf *surf, JacRes *jr, UserCtx *user)
{
	FDSTAG         *fs;
	const PetscInt *lx, *ly;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// store context
	surf->jr = jr;

	ierr = FreeSurfReadFromOptions(surf, &jr->scal); CHKERRQ(ierr);

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	fs = jr->fs;

	// get grid partitioning in X & Y directions
	ierr = DMDAGetOwnershipRanges(fs->DA_COR, &lx, &ly, NULL); CHKERRQ(ierr);

	// create redundant free surface DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_BOX,
		fs->dsx.tnods, fs->dsy.tnods, fs->dsz.nproc,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, NULL, &surf->DA_SURF); CHKERRQ(ierr);

	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->ltopo);  CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(surf->DA_SURF, &surf->gtopo);  CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vx);     CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vy);     CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vz);     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(surf->DA_SURF, &surf->vpatch); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(surf->DA_SURF, &surf->vmerge); CHKERRQ(ierr);

	// set initial internal free surface level
	ierr = VecSet(surf->ltopo, surf->InitLevel); CHKERRQ(ierr);
	ierr = VecSet(surf->gtopo, surf->InitLevel); CHKERRQ(ierr);


	// Set topo rom file if a Topo file is specified in the input
	PetscPrintf(PETSC_COMM_WORLD, "FileName: %s\n",user->TopoFilename);
	if(strcmp(user->TopoFilename,"noTopoFileName")!=0)
	{
		ierr = FreeSurfSetTopoFromFile(surf,user);
		CHKERRQ(ierr);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfReadFromOptions"
PetscErrorCode FreeSurfReadFromOptions(FreeSurf *surf, Scaling *scal)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read output flags
	ierr = PetscOptionsGetBool  (NULL, "-surf_use",       &surf->UseFreeSurf, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, "-surf_level",     &surf->InitLevel,   NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (NULL, "-surf_air_phase", &surf->AirPhase,    NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, "-surf_max_angle", &surf->MaxAngle,    NULL); CHKERRQ(ierr);

	// nondimensionalize
	surf->InitLevel /= scal->length;
	surf->MaxAngle  /= scal->angle;

	// set average topography & flag level
	surf->avg_topo     = surf->InitLevel;
	surf->jr->avg_topo = surf->InitLevel;
	surf->flat         = PETSC_TRUE;

	//======================================
	// read erosion sedimentation parameters
	//======================================

	ierr = GetIntDataItemCheck("-ErosionModel", "Erosion model",
		_NOT_FOUND_EXIT_, 1, &surf->ErosionModel, 0, 1); CHKERRQ(ierr);

	ierr = GetIntDataItemCheck("-SedimentModel", "Sedimentation model",
		_NOT_FOUND_EXIT_, 1, &surf->SedimentModel, 0, 1); CHKERRQ(ierr);

	// read prescribed sedimentation rate model parameter
	if(surf->SedimentModel == 1)
	{
		ierr = GetIntDataItemCheck("-numLayers", "Number of sediment layers",
			_NOT_FOUND_ERROR_, 1, &surf->numLayers, 1, _max_layers_); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-timeDelims", "Sediment layers time delimiters",
			_NOT_FOUND_ERROR_, surf->numLayers-1, surf->timeDelims, 0.0, 0.0, scal->time); CHKERRQ(ierr);

		ierr = GetScalDataItemCheckScale("-sedRates", "Sedimentation Rates",
			_NOT_FOUND_ERROR_, surf->numLayers, surf->sedRates, 0.0, 0.0, scal->velocity); CHKERRQ(ierr);

		ierr = GetIntDataItemCheck("-sedPhases", "Sediment layers phase numbers",
			_NOT_FOUND_ERROR_, surf->numLayers, surf->sedPhases, 0, 0); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfDestroy"
PetscErrorCode FreeSurfDestroy(FreeSurf *surf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	ierr = DMDestroy (&surf->DA_SURF); CHKERRQ(ierr);
	ierr = VecDestroy(&surf->ltopo);   CHKERRQ(ierr);
	ierr = VecDestroy(&surf->gtopo);   CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vx);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vy);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vz);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vpatch);  CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vmerge);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAdvect"
PetscErrorCode FreeSurfAdvect(FreeSurf *surf)
{
	// advect topography of the free surface mesh

	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	jr = surf->jr;

	// get surface velocities
	ierr = FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz); CHKERRQ(ierr);

	// advect topography
	ierr = FreeSurfAdvectTopo(surf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetVelComp"
PetscErrorCode FreeSurfGetVelComp(
	FreeSurf *surf,
	PetscErrorCode (*interp)(FDSTAG *, Vec, Vec, InterpFlags),
	Vec vcomp_grid, Vec vcomp_surf)
{
	// project velocity component from grid faces on the free surface

	// WARNING! this function has a problem if surface is placed on top boundary
	// most likely FindPointInCell has an issue

	JacRes      *jr;
	FDSTAG      *fs;
	Discret1D   *dsz;
	InterpFlags iflag;
	PetscInt    i, j, nx, ny, sx, sy, sz, level, K;
	PetscScalar ***topo, ***vsurf, ***vgrid, *vpatch, *vmerge, z, w;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr    = surf->jr;
	fs    = jr->fs;
	dsz   = &fs->dsz;
	level = dsz->rank;

	// create column communicator
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

	// set interpolation flags
	iflag.update    = PETSC_FALSE;
	iflag.use_bound = PETSC_TRUE;

	// interpolate velocity component from grid faces to corners
	ierr = interp(fs, vcomp_grid, jr->lbcor, iflag); CHKERRQ(ierr);

	// load ghost values
	LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)

	// clear surface velocity patch vector
	ierr = VecZeroEntries(surf->vpatch); CHKERRQ(ierr);

	// access topograpy, grid and surface velocity
	ierr = DMDAVecGetArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// get topography
		z = topo[level][j][i];

		// check whether point belongs to domain
		if(z >= dsz->crdbeg && z < dsz->crdend)
		{
			// find containing cell
			K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, z);

			// get interpolation weight
			w = (z - dsz->ncoor[K])/(dsz->ncoor[K+1] - dsz->ncoor[K]);

			// interpolate velocity
			vsurf[level][j][i] = (1.0 - w)*vgrid[sz+K][j][i] + w*vgrid[sz+K+1][j][i];
		}
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_COR,    jr->lbcor,    &vgrid); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &vsurf); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo,  &topo);  CHKERRQ(ierr);

	// merge velocity patches
	// compute ghosted version of the velocity component
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(surf->vpatch, &vpatch); CHKERRQ(ierr);
		ierr = VecGetArray(surf->vmerge, &vmerge); CHKERRQ(ierr);

		ierr = MPI_Allreduce(vpatch, vmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);

		ierr = VecRestoreArray(surf->vpatch, &vpatch); CHKERRQ(ierr);
		ierr = VecRestoreArray(surf->vmerge, &vmerge); CHKERRQ(ierr);

		// compute ghosted version of the velocity component
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->vmerge, vcomp_surf);
	}
	else
	{
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->vpatch, vcomp_surf);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAdvectTopo"
PetscErrorCode FreeSurfAdvectTopo(FreeSurf *surf)
{
	// advect topography on the free surface mesh

	JacRes      *jr;
	FDSTAG      *fs;
	PetscInt    I, I1, I2, J, J1, J2;
	PetscInt    i, j, jj, found, nx, ny, sx, sy, L, mx, my;
	PetscScalar cx[13], cy[13], cz[13];
	PetscScalar X, X1, X2, Y, Y1, Y2, Z, Exx, Eyy, step, gtol, avg_topo;
	PetscScalar ***advect, ***topo, ***vx, ***vy, ***vz;

	// local search grid triangulation
	PetscInt tria [] =
	{
		// first layer
		4, 5,  12, // 0
		4, 12, 7,  // 1
		4, 7,  11, // 2
		4, 11, 3,  // 3
		4, 3,  9,  // 4
		4, 9,  1,  // 5
		4, 1,  10, // 6
		4, 10, 5,  // 7

		// second layer
		5, 8, 12,  // 8
		8, 7, 12,  // 9
		7, 6, 11,  // 10
		6, 3, 11,  // 11
		3, 0, 9,   // 12
		0, 1, 9,   // 13
		1, 2, 10,  // 14
		2, 5, 10   // 15
	};

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr   = surf->jr;
	fs   = jr->fs;
	step = jr->ts.dt;
	mx   = fs->dsx.tnods;
	my   = fs->dsy.tnods;
	L    = fs->dsz.rank;
	gtol = jr->gtol;

	// get current background strain rates
	ierr = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, NULL); CHKERRQ(ierr);

	// access surface topography and velocity
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo, &advect); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &topo);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx,    &vx);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy,    &vy);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz,    &vz);     CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// get node coordinates
		X  = COORD_NODE(i,   sx, fs->dsx);
		X1 = COORD_NODE(i-1, sx, fs->dsx);
		X2 = COORD_NODE(i+1, sx, fs->dsx);
		Y  = COORD_NODE(j,   sy, fs->dsy);
		Y1 = COORD_NODE(j-1, sy, fs->dsy);
		Y2 = COORD_NODE(j+1, sy, fs->dsy);

		// get node indices
		I  = i;
		I1 = I-1; if(I1 == -1) I1 = I;
		I2 = I+1; if(I2 == mx) I2 = I;
		J  = j;
		J1 = J-1; if(J1 == -1) J1 = J;
		J2 = J+1; if(J2 == my) J2 = J;

		// compute deformed grid x-coordinates
		cx[0]  = step*vx[L][J1][I1] + X1;
		cx[1]  = step*vx[L][J1][I ] + X;
		cx[2]  = step*vx[L][J1][I2] + X2;
		cx[3]  = step*vx[L][J ][I1] + X1;
		cx[4]  = step*vx[L][J ][I ] + X;
		cx[5]  = step*vx[L][J ][I2] + X2;
		cx[6]  = step*vx[L][J2][I1] + X1;
		cx[7]  = step*vx[L][J2][I ] + X;
		cx[8]  = step*vx[L][J2][I2] + X2;
		cx[9]  = (cx[0] + cx[1] + cx[3] + cx[4])/4.0;
		cx[10] = (cx[1] + cx[2] + cx[4] + cx[5])/4.0;
		cx[11] = (cx[3] + cx[4] + cx[6] + cx[7])/4.0;
		cx[12] = (cx[4] + cx[5] + cx[7] + cx[8])/4.0;

		// compute deformed grid y-coordinates
		cy[0]  = step*vy[L][J1][I1] + Y1;
		cy[1]  = step*vy[L][J1][I ] + Y1;
		cy[2]  = step*vy[L][J1][I2] + Y1;
		cy[3]  = step*vy[L][J ][I1] + Y;
		cy[4]  = step*vy[L][J ][I ] + Y;
		cy[5]  = step*vy[L][J ][I2] + Y;
		cy[6]  = step*vy[L][J2][I1] + Y2;
		cy[7]  = step*vy[L][J2][I ] + Y2;
		cy[8]  = step*vy[L][J2][I2] + Y2;
		cy[9]  = (cy[0] + cy[1] + cy[3] + cy[4])/4.0;
		cy[10] = (cy[1] + cy[2] + cy[4] + cy[5])/4.0;
		cy[11] = (cy[3] + cy[4] + cy[6] + cy[7])/4.0;
		cy[12] = (cy[4] + cy[5] + cy[7] + cy[8])/4.0;

		// compute deformed grid z-coordinates
		cz[0]  = step*vz[L][J1][I1] + topo[L][J1][I1];
		cz[1]  = step*vz[L][J1][I ] + topo[L][J1][I ];
		cz[2]  = step*vz[L][J1][I2] + topo[L][J1][I2];
		cz[3]  = step*vz[L][J ][I1] + topo[L][J ][I1];
		cz[4]  = step*vz[L][J ][I ] + topo[L][J ][I ];
		cz[5]  = step*vz[L][J ][I2] + topo[L][J ][I2];
		cz[6]  = step*vz[L][J2][I1] + topo[L][J2][I1];
		cz[7]  = step*vz[L][J2][I ] + topo[L][J2][I ];
		cz[8]  = step*vz[L][J2][I2] + topo[L][J2][I2];
		cz[9]  = (cz[0] + cz[1] + cz[3] + cz[4])/4.0;
		cz[10] = (cz[1] + cz[2] + cz[4] + cz[5])/4.0;
		cz[11] = (cz[3] + cz[4] + cz[6] + cz[7])/4.0;
		cz[12] = (cz[4] + cz[5] + cz[7] + cz[8])/4.0;

		// compute updated node position if background strain rate is defined
		X *= (1.0 + step*Exx);
		Y *= (1.0 + step*Eyy);

		// find point in the deformed grid, interpolate topography
		found = 0;

		for(jj = 0; jj < 16; jj++)
		{
			found = InterpolateTriangle(cx, cy, cz, tria + 3*jj, X, Y, gtol, &Z);

			if(found) break;
		}

		if(!found)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Problems with topography advection");
		}

		// store advected topography
		advect[L][J][I] = Z;

	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &advect); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &topo);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx,    &vx);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy,    &vy);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz,    &vz);     CHKERRQ(ierr);

	// compute ghosted version of the advected surface topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	// set flat flag
	surf->flat = PETSC_FALSE;

	// compute & store average topography
	ierr = VecSum(surf->gtopo, &avg_topo); CHKERRQ(ierr);
	avg_topo /= (PetscScalar)(fs->dsx.tnods*fs->dsy.tnods*fs->dsz.nproc);
	surf->avg_topo = avg_topo;
	jr  ->avg_topo = avg_topo;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetAirPhaseRatio"
PetscErrorCode FreeSurfGetAirPhaseRatio(FreeSurf *surf)
{
	// compute proper phase ratio of air phase

	JacRes      *jr;
	FDSTAG      *fs;
	PetscScalar cx[5], cy[5], cz[5];
	PetscScalar ***topo, *phRat, vcell, phRatAir, gtol, cf;
	PetscScalar xleft, xright, yfront, yback, zbot, ztop;
	PetscInt    L, jj, iter, numPhases, AirPhase;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	// cell triangulation
	PetscInt tria [] =
	{
		0, 1, 4, // 0
		1, 3, 4, // 1
		3, 2, 4, // 2
		2, 0, 4  // 3
	};

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	jr        = surf->jr;
	AirPhase  = surf->AirPhase;
	fs        = jr->fs;
	gtol      = jr->gtol;
	numPhases = jr->numPhases;
	L         = fs->dsz.rank;
	iter      = 0;

	// access surface topography
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &topo); CHKERRQ(ierr);

	// scan all local cells
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access phase ratio array
		phRat = jr->svCell[iter++].phRat;

		// get cell bounds
		xleft  = COORD_NODE(i,   sx, fs->dsx);
		xright = COORD_NODE(i+1, sx, fs->dsx);

		yfront = COORD_NODE(j,   sy, fs->dsy);
		yback  = COORD_NODE(j+1, sy, fs->dsy);

		zbot   = COORD_NODE(k,   sz, fs->dsz);
		ztop   = COORD_NODE(k+1, sz, fs->dsz);

		// get cell volume
		vcell  = (xright - xleft)*(yback - yfront)*(ztop - zbot);

		// setup coordinate arrays
		cx[0]  = xleft;
		cx[1]  = xright;
		cx[2]  = xleft;
		cx[3]  = xright;
		cx[4]  = (xleft + xright)/2.0;

		cy[0]  = yfront;
		cy[1]  = yfront;
		cy[2]  = yback;
		cy[3]  = yback;
		cy[4]  = (yfront + yback)/2.0;

		cz[0]  = topo[L][j  ][i  ];
		cz[1]  = topo[L][j  ][i+1];
		cz[2]  = topo[L][j+1][i  ];
		cz[3]  = topo[L][j+1][i+1];
		cz[4]  = (cz[0] + cz[1] + cz[2] + cz[3])/4.0;

		// compute actual air phase ratio in the cell
		phRatAir = 1.0;

		for(jj = 0; jj < 4; jj++)
		{
			phRatAir -= IntersectTriangularPrism(cx, cy, cz, tria + 3*jj, vcell, zbot, ztop, gtol);
		}

		// normalize cell phase ratio if necessary
		if(phRat[AirPhase] != 1.0)
		{
			// get scaling factor
			cf = (1.0 - phRatAir)/(1.0 - phRat[AirPhase]);

			// scale solid phases
			for(jj = 0; jj < numPhases; jj++)
			{
				if(jj != AirPhase) phRat[jj] *= cf;
			}

			// correct air phase
			phRat[AirPhase] = phRatAir;
		}

		// WARNING !!!
		// think what to do if(phRat[AirPhase] == 1.0 && phRatAir != 1.0)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAppErosion"
PetscErrorCode FreeSurfAppErosion(FreeSurf *surf)
{
	// Apply fast erosion to the internal free surface of the model

	Scaling * scal;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = &surf->jr->scal;

	if(surf->ErosionModel == 1)
	{
		// erase topography
		ierr = VecSet(surf->ltopo, surf->avg_topo); CHKERRQ(ierr);
		ierr = VecSet(surf->gtopo, surf->avg_topo); CHKERRQ(ierr);

		// set flag
		surf->flat = PETSC_TRUE;

		PetscPrintf(PETSC_COMM_WORLD, "Applying infinitely fast erosion to internal free surface. Average free surface height = %e %s\n",
			surf->avg_topo*scal->length, scal->lbl_length);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAppSedimentation"
PetscErrorCode FreeSurfAppSedimentation(FreeSurf *surf)
{

	// Apply sedimentation to the internal free surface.
	// Currently we only have the option to add a fixed sedimentation rate,
	// and in this routine we simply advect the internal free surface upwards with
	// this rate. In the future we can think about adding different sedimentation routines.

	JacRes      *jr;
	FDSTAG      *fs;
	PetscScalar ***topo;
	PetscScalar dt, time, rate, zbot, ztop, z, dz, avg_topo;
	PetscInt    L, jj, phase;
	PetscInt    i, j, nx, ny, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr   = surf->jr;
	fs   = jr->fs;
	dt   = jr->ts.dt;
	time = jr->ts.time;
	L    = fs->dsz.rank;

	// get z-coordinates of the top and bottom boundaries
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, &ztop); CHKERRQ(ierr);

	if(surf->SedimentModel == 1)
	{
		// determine sedimentation rate & phase number
		for(jj = 0; jj < surf->numLayers-1; jj++)
		{
			if(time < surf->timeDelims[jj]) break;
		}

		rate  = surf->sedRates [jj];
		phase = surf->sedPhases[jj];

		// store the phase that is being sedimented
		surf->phase = phase;

		// get incremental thickness of the sediments
		dz = rate*dt;

		// access topography
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

		// scan all free surface local points
		ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);

		START_PLANE_LOOP
		{
			// get topography
			z = topo[L][j][i];

			// uniformly advect
			z += dz;

			// check if internal free surface goes outside the model domain
			if(z > ztop) z = ztop;
			if(z < zbot) z = zbot;

			// store advected topography
			topo[L][j][i] = z;
		}
		END_PLANE_LOOP

		// restore access
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

		// compute ghosted version of the topography
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

		// compute & store average topography
		ierr = VecSum(surf->gtopo, &avg_topo); CHKERRQ(ierr);
		avg_topo /= (PetscScalar)(fs->dsx.tnods*fs->dsy.tnods*fs->dsz.nproc);
		surf->avg_topo = avg_topo;
		jr  ->avg_topo = avg_topo;

		// print info
		PetscPrintf(PETSC_COMM_WORLD, "Applying sedimentation to internal free surface. Phase that is currently being sedimented is %lld   \n",
			(LLD)phase);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfSetTopoFromFile"
PetscErrorCode FreeSurfSetTopoFromFile(
	FreeSurf *surf, UserCtx *user)
{
	
	JacRes      *jr;
	FDSTAG      *fs;
	Discret1D   *dsz;
	PetscInt    i, j, nx, ny, sx, sy, sz, level;
	PetscScalar ***topo;
	PetscScalar avg_topo;
	char         *LoadFileName;
	PetscScalar  xp,yp, Xc, Yc, xpL, ypL;
	PetscInt 	nxTopo, nyTopo;
	PetscInt     Fsize;
	PetscInt Ix,Iy;
	PetscScalar  DX,DY;
	
	PetscScalar  *Topo;
	PetscScalar  chLen;
	int          fd;
	PetscScalar  header[2],dim[2];

	PetscViewer  view_in;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr    = surf->jr;
	fs    = jr->fs;
	dsz   = &fs->dsz;
	level = dsz->rank;


	// characteristic length
	chLen  = jr->scal.length;

	// create column communicator
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

	// set interpolation flags
	//iflag.update    = PETSC_FALSE;
	//iflag.use_bound = PETSC_TRUE;

	// interpolate velocity component from grid faces to corners
	//ierr = interp(fs, vcomp_grid, jr->lbcor, iflag); CHKERRQ(ierr);

	// load ghost values
	LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)

	// access topography, grid and surface velocity
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);


	// ****** Read file ******
	// create filename
	asprintf(&LoadFileName, "./%s/%s",
	user->LoadInitialParticlesDirectory,
	user->TopoFilename);

	PetscPrintf(PETSC_COMM_WORLD," Loading topo redundantly from file: %s \n", LoadFileName);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, LoadFileName, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header & size of file
	ierr = PetscBinaryRead(fd, &header, 2, PETSC_SCALAR); CHKERRQ(ierr);
	Fsize = (PetscInt)(header[1])-2;

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)Fsize*sizeof(PetscScalar), &Topo); CHKERRQ(ierr);

	// read entire file
	ierr = PetscBinaryRead(fd, &dim, 2,     PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, Topo, Fsize, PETSC_SCALAR); CHKERRQ(ierr);

	// grid spacing
	DX = user->W/(dim[0] - 1.0);
	DY = user->L/(dim[1] - 1.0);

	nxTopo = (PetscInt)dim[0];
	nyTopo = (PetscInt)dim[1];

	START_PLANE_LOOP
	{
		
		xp = COORD_NODE(i,   sx, fs->dsx);
		yp = COORD_NODE(j,   sy, fs->dsy);
		// index of the lower left corner of the element (of the temperature grid) in which the particle is
		Ix = (PetscInt)floor((xp - user->x_left) /DX);
		Iy = (PetscInt)floor((yp - user->y_front)/DY);
		
		// Take care of boundaries
		if (Ix == nxTopo-1)
		{
			Ix =nxTopo-2;
		}
		if (Iy == nyTopo-1)
		{
			Iy = nyTopo-2;
		}
		
		// Coordinate of the first corner (lower left deepest)
		Xc = user->x_left + (PetscScalar)Ix*DX;
		Yc = user->y_front+ (PetscScalar)Iy*DY;
		
		
		// Local coordinate of the particle inside a temperature element
		// Using the bilinear element in Kwon and Bang, p.161
		xpL = ( (xp - Xc)/DX )*2-1;
		ypL = ( (yp - Yc)/DY )*2-1;
		//zpL = (zp - Zc)/DZ;
		
		// Interpolate value on the particle using trilinear shape functions
		topo[level][j][i] = (
		1.0/4.0 * (1.0-xpL) * (1.0-ypL)  * Topo[Iy     * nxTopo + Ix   ] +
		1.0/4.0 * (1.0+xpL) * (1.0-ypL)  * Topo[Iy     * nxTopo + Ix+1 ] +
		1.0/4.0 * (1.0+xpL) * (1.0+ypL)  * Topo[(Iy+1) * nxTopo + Ix+1 ] +
		1.0/4.0 * (1.0-xpL) * (1.0+ypL)  * Topo[(Iy+1) * nxTopo + Ix   ])/chLen;
		
	}
	END_PLANE_LOOP
	
	

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);


	// compute ghosted version of the advected surface topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	// set flat flag
	surf->flat = PETSC_FALSE;

	// compute & store average topography
	ierr = VecSum(surf->gtopo, &avg_topo); CHKERRQ(ierr);
	avg_topo /= (PetscScalar)(fs->dsx.tnods*fs->dsy.tnods*fs->dsz.nproc);
	surf->avg_topo = avg_topo;
	jr  ->avg_topo = avg_topo;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// SERVICE FUNCTIONS
//---------------------------------------------------------------------------
PetscInt InterpolateTriangle(
	PetscScalar *x,   // x-coordinates of triangle
	PetscScalar *y,   // y-coordinates of triangle
	PetscScalar *f,   // interpolated field
	PetscInt    *i,   // indices of triangle corners
	PetscScalar  xp,  // x-coordinate of point
	PetscScalar  yp,  // y-coordinate of point
	PetscScalar  tol, // relative tolerance
	PetscScalar *fp)  // field value in the point
{
	PetscScalar xa, xb, xc, ya, yb, yc, la, lb, lc, A, S;

	// function returns:
	// 1 - if point is inside, fp contains the interpolated value
	// 0 - otherwise

	// access coordinates
	xa = x[i[0]];
	xb = x[i[1]];
	xc = x[i[2]];
	ya = y[i[0]];
	yb = y[i[1]];
	yc = y[i[2]];

	// triangle p-b-c (a-vertex)
	la = GET_AREA_TRIANG(xp, xb, xc, yp, yb, yc);

	// triangle p-c-a (b-vertex)
	lb = GET_AREA_TRIANG(xp, xc, xa, yp, yc, ya);

	// triangle p-a-b (c-vertex)
	lc = GET_AREA_TRIANG(xp, xa, xb, yp, ya, yb);

	// triangle a-b-c (total area)
	A  = GET_AREA_TRIANG(xa, xb, xc, ya, yb, yc);

	// sum
	S = la + lb + lc;

	// perform point test
	if(S > A*(1.0 + tol)) return 0;

	// perform interpolation
	la /= S;
	lb /= S;
	lc /= S;

	(*fp) = la*f[i[0]] + lb*f[i[1]] + lc*f[i[2]];

	return 1;
}
//---------------------------------------------------------------------------
PetscScalar IntersectTriangularPrism(
	PetscScalar *x,     // x-coordinates of prism base
	PetscScalar *y,     // y-coordinates of prism base
	PetscScalar *z,     // z-coordinates of prism top surface
	PetscInt    *i,     // indices of base corners
	PetscScalar  vcell, // total volume of cell
	PetscScalar  bot,   // z-coordinate of bottom plane
	PetscScalar  top,   // z-coordinate of top plane
	PetscScalar  tol)   // relative tolerance
{
    // compute prism volume cut by top and bottom horizontal planes
	// relative to the total volume of the cell

	PetscScalar xa, xb, xc, ya, yb, yc, za, zb, zc;
	PetscScalar xab, xbc, xca, yab, ybc, yca, zab, zbc, zca;
	PetscScalar dh, w, vbot, vtop, zmin, zmax;

	// access coordinates
	xa = x[i[0]];
	xb = x[i[1]];
	xc = x[i[2]];
	ya = y[i[0]];
	yb = y[i[1]];
	yc = y[i[2]];
	za = z[i[0]];
	zb = z[i[1]];
	zc = z[i[2]];

	// get absolute tolerance
	dh = (top-bot)*tol;

	// get vertical coordinate range of the prism top surface
	zmin = za; if(zb < zmin) zmin = zb; if(zc < zmin) zmin = zc;
	zmax = za; if(zb > zmax) zmax = zb; if(zc > zmax) zmax = zc;

	// check for empty cell
	if(zmax <= bot)
	{
		return 0.0;
	}

	// check for filled cell
	if(zmin >= top)
	{
		return 0.25;
	}

	// get volume above bottom plane
	vbot = 0.0;

	// determine edge intersection points
	INTERSECT_EDGE(xa, ya, za, xb, yb, zb, xab, yab, zab, bot, dh); // edge a-b
	INTERSECT_EDGE(xb, yb, zb, xc, yc, zc, xbc, ybc, zbc, bot, dh); // edge b-c
	INTERSECT_EDGE(xc, yc, zc, xa, ya, za, xca, yca, zca, bot, dh); // edge c-a

	// compute volume
	vbot += GET_VOLUME_PRISM(xa,  xab, xca, ya,  yab, yca, za,  zab, zca, bot); // prism a--ab-ca
	vbot += GET_VOLUME_PRISM(xb,  xbc, xab, yb,  ybc, yab, zb,  zbc, zab, bot); // prism b--bc-ab
	vbot += GET_VOLUME_PRISM(xc,  xca, xbc, yc,  yca, ybc, zc,  zca, zbc, bot); // prism c--ca-bc
	vbot += GET_VOLUME_PRISM(xab, xbc, xca, yab, ybc, yca, zab, zbc, zca, bot); // prism ab-bc-ca

	// get volume above top plane
	vtop = 0.0;

	if(zmax > top)
	{
		// determine edge intersection points
		INTERSECT_EDGE(xa, ya, za, xb, yb, zb, xab, yab, zab, top, dh); // edge a-b
		INTERSECT_EDGE(xb, yb, zb, xc, yc, zc, xbc, ybc, zbc, top, dh); // edge b-c
		INTERSECT_EDGE(xc, yc, zc, xa, ya, za, xca, yca, zca, top, dh); // edge c-a

		// compute volume
		vtop += GET_VOLUME_PRISM(xa,  xab, xca, ya,  yab, yca, za,  zab, zca, top); // prism a--ab-ca
		vtop += GET_VOLUME_PRISM(xb,  xbc, xab, yb,  ybc, yab, zb,  zbc, zab, top); // prism b--bc-ab
		vtop += GET_VOLUME_PRISM(xc,  xca, xbc, yc,  yca, ybc, zc,  zca, zbc, top); // prism c--ca-bc
		vtop += GET_VOLUME_PRISM(xab, xbc, xca, yab, ybc, yca, zab, zbc, zca, top); // prism ab-bc-ca
	}

	// volume inside cell = volume above bottom plane - volume above top plane
	// NOTE! volume returned by the macro is two times larger than actual value
	return (vbot - vtop)/2.0/vcell;
}
//---------------------------------------------------------------------------
/*
{
	// TEST

	PetscScalar bot, top, gtol, v1, v2;

	PetscScalar x[] = { 1.0, 7.0, 2.0 };
	PetscScalar y[] = { 2.0, 1.0, 5.0 };
	PetscScalar z[] = { 1.0, 2.0, 3.0 };
	PetscInt    i[] = { 0,   1,   2   };

	gtol  = 1e-12;

	bot   = 1.0;
	top   = 1.5;

	v1 = IntersectTriangularPrism(x, y, z, i, 1.0, bot, top, gtol);

	bot   = 1.5;
	top   = 3.0;

	v2 = IntersectTriangularPrism(x, y, z, i, 1.0, bot, top, gtol);


	printf("\n\n\n v1: %f \n\n\n", v1);

	printf("\n\n\n v2: %f \n\n\n", v2);

	printf("\n\n\n sum: %f \n\n\n", v1+v2);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetSinusoidalPerturbation"
// There are some cases in which we set an initial sinusoidal perturbation on the free surface
PetscErrorCode SetSinusoidalPerturbation(PetscScalar SinusoidalFreeSurfaceAmplitude, UserCtx *user)
{
	PetscErrorCode ierr;
	PetscScalar    ***LocalSurfaceTopography, x,z, maxVec;
	PetscInt	   xs_Z, ys_Z, zs_Z, xm_Z, ym_Z, zm_Z, ix, iy;
	DM			   cda_SurfaceTopo;
	Vec			   gc_SurfaceTopo;
	DMDACoor3d	   ***coors_SurfaceTopo;
	// nondimensionalize
	SinusoidalFreeSurfaceAmplitude = SinusoidalFreeSurfaceAmplitude/user->Characteristic.Length;
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,	&cda_SurfaceTopo); 	                         CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography, &gc_SurfaceTopo); 	                     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo, &coors_SurfaceTopo); 	                     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,user->SurfaceTopography, &LocalSurfaceTopography); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z);               CHKERRQ(ierr);
	for (iy=ys_Z; iy<ys_Z+ym_Z; iy++)
	{	for(ix=xs_Z; ix<xs_Z+xm_Z; ix++)
		{	// Extract x,y,z coordinates of surface topography
			x = coors_SurfaceTopo[zs_Z][iy][ix].x;
//			y = coors_SurfaceTopo[zs_Z][iy][ix].y;
			z = LocalSurfaceTopography[zs_Z][iy][ix];
			z = z + cos(x/user->W*2*PETSC_PI)*SinusoidalFreeSurfaceAmplitude;
			// set topography
			LocalSurfaceTopography[zs_Z][iy][ix] = z;
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);
	VecMax(user->SurfaceTopography,PETSC_NULL, &maxVec);
	PetscPrintf(PETSC_COMM_WORLD,"max topo = %f ", maxVec*user->Characteristic.Length/1000.0);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// get partitioning of the free surface in the XY plane
// ierr = FreeSurfGetPartition(&fs->dsx, xs, dx, fs->dsx.h_min, &nx, &lx); CHKERRQ(ierr);
// ierr = FreeSurfGetPartition(&fs->dsy, ys, dy, fs->dsy.h_min, &ny, &ly); CHKERRQ(ierr);
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetPartition"
PetscErrorCode FreeSurfGetPartition(
	Discret1D    *ds,  // discretization
	PetscScalar   beg, // starting coordinate
	PetscScalar   len, // domain length
	PetscScalar   h,   // target mesh step
	PetscInt     *n,   // total number of nodes
	PetscInt    **l)   // free surface partitioning vector
{
	MPI_Comm    comm;
	PetscScalar step, tol, rtol;
	PetscInt    i, first, last, lnum, nnod, ncel, sum, *part;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	rtol = 1e-8;

	// compute total number of cells & nodes
	ncel = (PetscInt)ceil(len/h);
	nnod = ncel + 1;

	if(ds->nproc == 1)
	{
		// single processor case
		part = NULL;
	}
	else
	{
		// compute actual mesh step & geometric tolerance
		step = len/(PetscScalar)ncel;
		tol  = rtol*step;

		// get index of first local node
		if(ds->grprev != -1) first = (PetscInt)floor((ds->ncoor[0] - beg + tol)/step);
		else                 first = 0;

		// get index of last local node
		if(ds->grnext != -1) last = (PetscInt)floor((ds->ncoor[ds->ncels] - beg - tol)/step);
		else                 last = ncel;

		// get local number of nodes
		lnum = last - first + 1;

		// create column communicator
		ierr = Discret1DGetColumnComm(ds, &comm); CHKERRQ(ierr);

		// create partitioning vector
		ierr = makeIntArray(&part, NULL, ds->nproc); CHKERRQ(ierr);

		// gather partitioning vector
		ierr = MPI_Allgather(&lnum, 1, MPIU_INT, part, 1, MPIU_INT, comm); CHKERRQ(ierr);

		// free
		ierr = MPI_Comm_free(&comm); CHKERRQ(ierr);

		// checksum
		for(i = 0, sum = 0; i < ds->nproc; i++) sum += part[i];

		if(sum != nnod) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Inconsistent free surface partitioning");
	}

	// return partitioning
	(*n) = nnod;
	(*l) = part;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
