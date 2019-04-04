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
#include "surf.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "bc.h"
#include "phase.h"
#include "JacRes.h"
#include "interpolate.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * stair-case type of free surface
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfCreate"
PetscErrorCode FreeSurfCreate(FreeSurf *surf, FB *fb)
{
	Scaling  *scal;
	PetscInt  maxPhaseID;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	surf->phaseCorr   =  1;
	surf->AirPhase    = -1;

	// check whether free surface is activated
	ierr = getIntParam(fb, _OPTIONAL_, "surf_use", &surf->UseFreeSurf, 1,  1); CHKERRQ(ierr);

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	scal       = surf->jr->scal;
	maxPhaseID = surf->jr->dbm->numPhases-1;

	// read from options
	ierr = getIntParam   (fb, _OPTIONAL_, "surf_corr_phase",    &surf->phaseCorr,     1,  1);            CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "surf_level",         &surf->InitLevel,     1,  scal->length); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, "surf_air_phase",     &surf->AirPhase,      1,  maxPhaseID);   CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "surf_max_angle",     &surf->MaxAngle,      1,  scal->angle);  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "erosion_model",      &surf->ErosionModel,  1,  1);            CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "sediment_model",     &surf->SedimentModel, 1,  3);            CHKERRQ(ierr);

	if(surf->SedimentModel == 1 || surf->SedimentModel == 2 || surf->SedimentModel == 3 )
	{
		// sedimentation model parameters
		ierr = getIntParam   (fb, _REQUIRED_, "sed_num_layers",  &surf->numLayers,  1,                 _max_sed_layers_);  CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "sed_time_delims",  surf->timeDelims, surf->numLayers-1, scal->time);        CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "sed_rates",        surf->sedRates,   surf->numLayers,   scal->velocity);    CHKERRQ(ierr);
		ierr = getIntParam   (fb, _REQUIRED_, "sed_phases",       surf->sedPhases,  surf->numLayers,   maxPhaseID);        CHKERRQ(ierr);
	}

	if(surf->SedimentModel == 2)
	{
		// sedimentation model parameters
		ierr = getScalarParam(fb, _REQUIRED_, "marginO",          surf->marginO,    2,                 scal->length);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "marginE",          surf->marginE,    2,                 scal->length);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "marginE",          surf->marginE,    2,                 scal->length);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "hUp",             &surf->hUp,        1,                 scal->length);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "hDown",           &surf->hDown,      1,                 scal->length);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "dTrans",          &surf->dTrans,     1,                 scal->length);      CHKERRQ(ierr);
	}

	if (surf->SedimentModel == 3)
	{
		ierr = getScalarParam(fb, _REQUIRED_, "sed_rates2nd",        surf->sedRates2nd,   surf->numLayers,   scal->velocity);  CHKERRQ(ierr);
	}

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Free surface parameters: \n");

	PetscPrintf(PETSC_COMM_WORLD, "   Sticky air phase ID       : %lld \n",  (LLD)surf->AirPhase);
    PetscPrintf(PETSC_COMM_WORLD, "   Initial surface level     : %g %s \n",  surf->InitLevel*scal->length, scal->lbl_length);

	PetscPrintf(PETSC_COMM_WORLD, "   Erosion model             : ");
	if      (surf->ErosionModel == 0)  PetscPrintf(PETSC_COMM_WORLD, "none\n");
	else if (surf->ErosionModel == 1)  PetscPrintf(PETSC_COMM_WORLD, "infinitely fast\n");

	PetscPrintf(PETSC_COMM_WORLD, "   Sedimentation model       : ");
	if      (surf->SedimentModel == 0) PetscPrintf(PETSC_COMM_WORLD, "none\n");
	else if (surf->SedimentModel == 1) PetscPrintf(PETSC_COMM_WORLD, "prescribed rate\n");
	else if (surf->SedimentModel == 2) PetscPrintf(PETSC_COMM_WORLD, "directed sedimentation (continental margin) with prescribed rate\n");

	if(surf->numLayers) PetscPrintf(PETSC_COMM_WORLD, "   Number of sediment layers : %lld \n",  (LLD)surf->numLayers);
	if(surf->phaseCorr) PetscPrintf(PETSC_COMM_WORLD, "   Correct marker phases     @ \n");
	if(surf->MaxAngle)  PetscPrintf(PETSC_COMM_WORLD, "   Maximum surface slope     : %g %s\n",  surf->MaxAngle*scal->angle, scal->lbl_angle);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// create structures
	ierr = FreeSurfCreateData(surf); CHKERRQ(ierr);

	// set initial internal free surface level
	ierr = VecSet(surf->gtopo, surf->InitLevel); CHKERRQ(ierr);
	ierr = VecSet(surf->ltopo, surf->InitLevel); CHKERRQ(ierr);

	// initialize topography from file if provided
	ierr = FreeSurfSetTopoFromFile(surf, fb); CHKERRQ(ierr);

	// compute & store average topography
	ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfCreateData"
PetscErrorCode FreeSurfCreateData(FreeSurf *surf)
{
	FDSTAG         *fs;
	const PetscInt *lx, *ly;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = surf->jr->fs;

	// get grid partitioning in X & Y directions
	ierr = DMDAGetOwnershipRanges(fs->DA_COR, &lx, &ly, NULL); CHKERRQ(ierr);

	// create redundant free surface DMDA
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetAvgTopo"
PetscErrorCode FreeSurfGetAvgTopo(FreeSurf *surf)
{
	JacRes      *jr;
	FDSTAG      *fs;
	PetscScalar  avg_topo;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	jr = surf->jr;
	fs = jr->fs;

	// compute & set average topography
	ierr = VecSum(surf->gtopo, &avg_topo); CHKERRQ(ierr);

	avg_topo /= (PetscScalar)(fs->dsx.tnods*fs->dsy.tnods*fs->dsz.nproc);

	surf->avg_topo = avg_topo;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfReadRestart"
PetscErrorCode FreeSurfReadRestart(FreeSurf *surf, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// create structures
	ierr = FreeSurfCreateData(surf); CHKERRQ(ierr);

	// read topography vector
	ierr = VecReadRestart(surf->gtopo, fp); CHKERRQ(ierr);

	// get ghosted topography vector
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfWriteRestart"
PetscErrorCode FreeSurfWriteRestart(FreeSurf *surf, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// store topography vector
	ierr = VecWriteRestart(surf->gtopo, fp); CHKERRQ(ierr);

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
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

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
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	jr = surf->jr;

	// get surface velocities
	ierr = FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz); CHKERRQ(ierr);

	// advect topography
	ierr = FreeSurfAdvectTopo(surf); CHKERRQ(ierr);

	// smooth topography spikes
	ierr = FreeSurfSmoothMaxAngle(surf); CHKERRQ(ierr);

	// compute & store average topography
	ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

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
	InterpFlags iflags;
	PetscScalar bz, ez;
	PetscInt    i, j, nx, ny, sx, sy, sz, level, K;
	PetscScalar ***topo, ***vsurf, ***vgrid, *vpatch, *vmerge, z, w;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr    = surf->jr;
	fs    = jr->fs;
	dsz   = &fs->dsz;
	level = (PetscInt)dsz->rank;

	// get local coordinate bounds
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);

	// create column communicator
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

	// set interpolation flags
	iflags.update    = 0; // overwrite vectors
	iflags.use_bound = 1; // use boundary values

	// interpolate velocity component from grid faces to corners
	ierr = interp(fs, vcomp_grid, jr->lbcor, iflags); CHKERRQ(ierr);

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
		if(z >= bz && z < ez)
		{
			// find containing cell
			ierr = Discret1DFindPoint(&fs->dsz, z, K); CHKERRQ(ierr);

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
	PetscScalar X, X1, X2, Y, Y1, Y2, Z, Exx, Eyy, Rxx, Ryy, step, gtol;
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
	step = jr->ts->dt;
	mx   = fs->dsx.tnods;
	my   = fs->dsy.tnods;
	L    = (PetscInt)fs->dsz.rank;
	gtol = fs->gtol;

	// get current background strain rates
	ierr = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, NULL, &Rxx, &Ryy, NULL); CHKERRQ(ierr);

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
		X += step*Exx*(X - Rxx);
		Y += step*Eyy*(Y - Ryy);

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfSmoothMaxAngle"
PetscErrorCode FreeSurfSmoothMaxAngle(FreeSurf *surf)
{
	// smooth topography if maximum angle with horizon is exceeded

	JacRes      *jr;
	FDSTAG      *fs;
	Vec         cellTopo;
	PetscScalar ***ntopo, ***ctopo;
	PetscScalar tanMaxAng, zbot, dx, dy, h, t, tmax, cz[4], Ezz, Rzz, step;
	PetscInt    i, j, nx, ny, sx, sy, L, cnt, gcnt, I1, I2, J1, J2, mx, my;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether smoothing is activated
	if(surf->MaxAngle == 0.0) PetscFunctionReturn(0);

	// access context
	jr        = surf->jr;
	fs        = jr->fs;
	mx        = fs->dsx.tnods - 1;
	my        = fs->dsy.tnods - 1;
	L         = (PetscInt)fs->dsz.rank;
	step      = jr->ts->dt;
	tanMaxAng = PetscTanReal(surf->MaxAngle);

	// get global coordinate bounds
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, NULL); CHKERRQ(ierr);

	// get current background strain rates
	ierr = BCGetBGStrainRates(jr->bc, NULL, NULL, &Ezz, NULL, NULL, &Rzz); CHKERRQ(ierr);

	// update position of bottom boundary
	zbot += step*Ezz*(zbot - Rzz);

	// get cell topography vector
	ierr = DMGetLocalVector(jr->DA_CELL_2D, &cellTopo); CHKERRQ(ierr);

	ierr = VecZeroEntries(cellTopo); CHKERRQ(ierr);

	// access cell topography
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, cellTopo, &ctopo); CHKERRQ(ierr);

	// access surface topography (corner nodes)
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,  &ntopo); CHKERRQ(ierr);

	// scan all local cells
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	cnt = 0;

	START_PLANE_LOOP
	{
		// get horizontal cell sizes
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);

		// get topography at cell corners
		cz[0] = ntopo[L][j  ][i  ];
		cz[1] = ntopo[L][j  ][i+1];
		cz[2] = ntopo[L][j+1][i  ];
		cz[3] = ntopo[L][j+1][i+1];

		// estimate maximum free surface deviation from horizon
		t = PetscAbsScalar(cz[1] - cz[0])/dx;              tmax = t;
		t = PetscAbsScalar(cz[3] - cz[2])/dx; if(t > tmax) tmax = t;
		t = PetscAbsScalar(cz[2] - cz[0])/dy; if(t > tmax) tmax = t;
		t = PetscAbsScalar(cz[3] - cz[1])/dy; if(t > tmax) tmax = t;

		// get average cell height
		h = (cz[0] + cz[1] + cz[2] + cz[3])/4.0 - zbot;

		// mark (with negative value) and count local affected cells
		if(tmax > tanMaxAng) { h = -h; cnt++; }

		// store cell topography
		ctopo[L][j][i] = h;
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, cellTopo, &ctopo); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &ntopo); CHKERRQ(ierr);

	// count global affected cells
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&cnt, &gcnt, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gcnt = cnt;
	}

	// return if topography is within the limits
	if(!gcnt)
	{
		ierr = DMRestoreLocalVector(jr->DA_CELL_2D, &cellTopo); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// fill ghost points
	LOCAL_TO_LOCAL(jr->DA_CELL_2D, cellTopo)

	// access cell topography
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, cellTopo, &ctopo); CHKERRQ(ierr);

	// access surface topography (corner nodes)
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &ntopo); CHKERRQ(ierr);

	// scan all local nodes
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// get topography at neighboring cell centers
		cz[0] = ctopo[L][J1][I1];
		cz[1] = ctopo[L][J1][I2];
		cz[2] = ctopo[L][J2][I1];
		cz[3] = ctopo[L][J2][I2];

		// check whether nodal topography needs correction
		cnt = 0;

		if(cz[0] < 0.0) { cz[0] = -cz[0]; cnt++; }
		if(cz[1] < 0.0) { cz[1] = -cz[1]; cnt++; }
		if(cz[2] < 0.0) { cz[2] = -cz[2]; cnt++; }
		if(cz[3] < 0.0) { cz[3] = -cz[3]; cnt++; }

		// smooth nodal topography
		if(cnt)
		{
			ntopo[L][j][i] = (cz[0] + cz[1] + cz[2] + cz[3])/4.0 + zbot;
		}
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, cellTopo, &ctopo); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &ntopo); CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(jr->DA_CELL_2D, &cellTopo); CHKERRQ(ierr);

	// fill ghost points
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

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
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// check whether phase correction is activated
	if(!surf->phaseCorr) PetscFunctionReturn(0);

	// access context
	jr        = surf->jr;
	AirPhase  = surf->AirPhase;
	fs        = jr->fs;
	gtol      = fs->gtol;
	numPhases = jr->dbm->numPhases;
	L         = (PetscInt)fs->dsz.rank;
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

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	scal = surf->jr->scal;

	if(surf->ErosionModel == 1)
	{
		// erase topography
		ierr = VecSet(surf->gtopo, surf->avg_topo); CHKERRQ(ierr);
		ierr = VecSet(surf->ltopo, surf->avg_topo); CHKERRQ(ierr);

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
	PetscScalar dt, time, rate, zbot, ztop, z,zprop,zpropn, dz, dr, x,y;
	PetscScalar rsq, rsqn, t0, t0n, l[2], aE[2], aO[2],ln[2], aEn[2], aOn[2], b[2],c[2],d,dn;
	PetscInt    L, jj, phase;
	PetscInt    i, j, nx, ny, sx, sy, sz;
	PetscScalar   BoxWidth, ex,bx,dz_x,dz1,dz2,rate1,rate2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	jr   = surf->jr;
	fs   = jr->fs;
	dt   = jr->ts->dt;
	time = jr->ts->time;
	L    = (PetscInt)fs->dsz.rank;


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
		ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

		// print info
		PetscPrintf(PETSC_COMM_WORLD, "Applying sedimentation to internal free surface. Phase that is currently being sedimented is %lld   \n",
			(LLD)phase);
	}
	else if(surf->SedimentModel == 2) 
	{ 
		// sedimentation after Gemmer et al. 2004 - Moving Gaussian to mimic the sedimentation at a continental margin

		// determine sedimentation rate & phase number
		for(jj = 0; jj < surf->numLayers-1; jj++)
		{
			if(time < surf->timeDelims[jj]) break;
		}

		rate  = surf->sedRates [jj];
		phase = surf->sedPhases[jj];

		// store the phase that is being sedimented
		surf->phase = phase;

		// strike of current margin
		b[0] = surf->marginE[0] - surf->marginO[0];
		b[1] = surf->marginE[1] - surf->marginO[1];

		// right lateral normal vector of current margin
		c[0] = -b[1]/sqrt(b[0]*b[0]+b[1]*b[1]);
		c[1] = b[0]/sqrt(b[0]*b[0]+b[1]*b[1]);

		// lateral offset
		dr = rate*dt;

		// current margin 
		aO[0] = surf->marginO[0];  
		aO[1] = surf->marginO[1]; 
		aE[0] = surf->marginE[0];  
		aE[1] = surf->marginE[1]; 
		
		// new margin with offset computed with given rate
		aOn[0] = surf->marginO[0]-dr*c[0];  
		aOn[1] = surf->marginO[1]-dr*c[1]; 
		aEn[0] = surf->marginE[0]-dr*c[0];  
		aEn[1] = surf->marginE[1]-dr*c[1]; 

		// access topography
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

		// scan all free surface local points
		ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);
		
		START_PLANE_LOOP
		{

			x = COORD_NODE(i, sx, fs->dsx);
			y = COORD_NODE(j, sy, fs->dsy);

			// get topography
			z = topo[L][j][i];


			// compute distance to margin

			t0    = ((x-aO[0])*b[0] + (y-aO[1])*b[1]) / (b[0]*b[0]+b[1]*b[1]);
			t0n   = ((x-aOn[0])*b[0] + (y-aOn[1])*b[1]) / (b[0]*b[0]+b[1]*b[1]);
			l[0]  = aO[0] + t0 * b[0];
			l[1]  = aO[1] + t0 * b[1];
			ln[0] = aOn[0] + t0n * b[0];
			ln[1] = aOn[1] + t0n * b[1];
			rsq   = ((x-l[0]) *(x-l[0])+(y-l[1])*(y-l[1]));
			rsqn  = ((x-ln[0]) *(x-ln[0])+(y-ln[1])*(y-ln[1]));
			
			// identify side of margin (left lateral < 0)
			dn = (x-aOn[0])*(aEn[1]-aOn[1]) - (y-aOn[1])*(aEn[0]-aOn[0]);
			d  = (x-aO[0])*(aE[1]-aO[1]) - (y-aO[1])*(aE[0]-aO[0]);

			if(dn<0) 
			{

				zprop = surf->InitLevel+surf->hUp;
			}
			else if(d<0)
			{
				zpropn = surf->InitLevel+surf->hUp;
			}
			else
			{
				zprop  = surf->InitLevel + surf->hDown + (surf->hUp-surf->hDown) * exp(-rsq /surf->dTrans/surf->dTrans );
				zpropn = surf->InitLevel + surf->hDown + (surf->hUp-surf->hDown) * exp(-rsqn /surf->dTrans/surf->dTrans );
			}
			
			dz = zpropn-zprop;

			if (time == 0.0) z = zprop;
			if (z <= zprop) z += dz;
			
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
		ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

		// update margin
		surf->marginO[0] = aOn[0];
		surf->marginO[1] = aOn[1];
		surf->marginE[0] = aEn[0];
		surf->marginE[1] = aEn[1];

		// print info
		PetscPrintf(PETSC_COMM_WORLD, "Applying directed (cont. margin) sedimentation to internal free surface. Phase that is currently being sedimented is %lld   \n",
			(LLD)phase);
	}
	else if(surf->SedimentModel == 3)
	{
		// sedimentation after Gemmer et al. 2004 - Moving Gaussian to mimic the sedimentation at a continental margin

		// determine sedimentation rate & phase number
		for(jj = 0; jj < surf->numLayers-1; jj++)
		{
			if(time < surf->timeDelims[jj]) break;
		}

		rate1  = surf->sedRates [jj];
		rate2  = surf->sedRates2nd [jj];
		phase = surf->sedPhases[jj];

		// store the phase that is being sedimented
		surf->phase = phase;

		// lateral offset
		dz1 = rate1*dt;
		dz2 = rate2*dt;
		dz = dz1-dz2;

		// access topography
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

		// scan all free surface local points
		ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);
		
		// Get Global Box extent
		ierr = FDSTAGGetGlobalBox(fs, &bx, 0, 0, &ex,0, 0); CHKERRQ(ierr);
		BoxWidth = ex-bx;


		START_PLANE_LOOP
		{

			x = COORD_NODE(i, sx, fs->dsx);
			y = COORD_NODE(j, sy, fs->dsy);

			// get topography
			z = topo[L][j][i];


			dz_x = -dz/BoxWidth * x + dz1;

			// uniformly advect
			z += dz_x;

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
		ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

		// print info
		PetscPrintf(PETSC_COMM_WORLD, "Applying differential loading to internal free surface. Phase that is currently being sedimented is %lld   \n",
			(LLD)phase);
	}

	
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfSetTopoFromFile"
PetscErrorCode FreeSurfSetTopoFromFile(FreeSurf *surf, FB *fb)
{
	FDSTAG         *fs;
	int            fd;
	PetscLogDouble t;
	PetscViewer    view_in;
	char           filename[_str_len_];
	PetscInt 	   nxTopo, nyTopo, Ix, Iy, GridSize;
	PetscInt       i, j, nx, ny, sx, sy, sz, level;
	PetscScalar    ***topo, *Z, header, dim[2], start[2], spacing[2];
	PetscScalar    xp, yp, xpL, ypL, DX, DY, bx, by, ex, ey, leng, X1, Y1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "surf_topo_file", filename, NULL); CHKERRQ(ierr);

	// check whether file is provided
	if(!strlen(filename)) PetscFunctionReturn(0);

	PrintStart(&t, "Loading topography redundantly from", filename);

	// access context
	fs    = surf->jr->fs;
	level = fs->dsz.rank;
	leng  = surf->jr->scal->length;

	// open & read file
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd);                               CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header
	ierr = PetscBinaryRead(fd, &header, 1, PETSC_SCALAR); CHKERRQ(ierr);

	// read grid dimensions
	ierr = PetscBinaryRead(fd, &dim, 2,  PETSC_SCALAR); CHKERRQ(ierr);

	// read south-west corner coordinates
	ierr = PetscBinaryRead(fd, &start, 2,  PETSC_SCALAR); CHKERRQ(ierr);

	// read grid spacing
	ierr = PetscBinaryRead(fd, &spacing, 2,  PETSC_SCALAR); CHKERRQ(ierr);

	// compute grid size
	nxTopo = (PetscInt)dim[0];
	nyTopo = (PetscInt)dim[1];
	GridSize = nxTopo * nyTopo;

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)GridSize*sizeof(PetscScalar), &Z); CHKERRQ(ierr);
	
	// read entire file
	ierr = PetscBinaryRead(fd, Z, GridSize, PETSC_SCALAR); CHKERRQ(ierr);

	// destroy file handle
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);

	// get input topography grid spacing
	DX = spacing[0];
	DY = spacing[1];

	// get input topography south-west corner
	X1 = start[0];
	Y1 = start[1];

	// get mesh extents
	ierr = FDSTAGGetGlobalBox(fs, &bx, &by, 0, &ex, &ey, 0); CHKERRQ(ierr);

	// access topography vector
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);

	// check if input grid covers at least the entire LaMEM grid
	if(X1 > bx){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Topography input file does not cover western edge of the LaMEM box!");
	}

	if(X1+(nxTopo-1)*DX < ex){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Topography input file does not cover eastern edge of the LaMEM box!");
	}

	if(Y1 > by){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Topography input file does not cover southern edge of the LaMEM box!");
	}

	if(Y1+(nyTopo-1)*DY < ey){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Topography input file does not cover northern edge of the LaMEM box!");
	}


	// runs over all LaMEM nodes
	START_PLANE_LOOP
	{
		// get node coordinate
		xp = COORD_NODE(i, sx, fs->dsx);
		yp = COORD_NODE(j, sy, fs->dsy);

		// check in which element of the input grid the LaMEM node is
		Ix = (PetscInt)PetscFloorReal((xp-X1)/DX);
		Iy = (PetscInt)PetscFloorReal((yp-Y1)/DY);

		// in case the last LaMEM node is identical with the last input grid node
		if (Ix == nxTopo - 1) Ix = nxTopo - 2;
		if (Iy == nyTopo - 1) Iy = nyTopo - 2;

		// get relative coordinates of the LaMEM node in relation to SW node of inout grid element
		xpL = (PetscScalar)((xp-(X1+(Ix*DX)))/DX);
		ypL = (PetscScalar)((yp-(Y1+(Iy*DY)))/DY);

		// interpolate topography from input grid onto LaMEM nodes
		topo[level][j][i] = (
			(1.0-xpL) * (1.0-ypL) * Z[Ix   + Iy     * nxTopo] +
			(xpL)     * (1.0-ypL) * Z[Ix+1 + Iy     * nxTopo] +
			(xpL)     * (ypL)     * Z[Ix+1 + (Iy+1) * nxTopo] +
			(1.0-xpL) * (ypL)     * Z[Ix   + (Iy+1) * nxTopo])/leng;
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo);  CHKERRQ(ierr);

	// clear memory
	PetscFree(Z);

	// compute ghosted version of the advected surface topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	PrintDone(t);

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

	PetscScalar bot, top, g_tol, v1, v2;

	PetscScalar x[] = { 1.0, 7.0, 2.0 };
	PetscScalar y[] = { 2.0, 1.0, 5.0 };
	PetscScalar z[] = { 1.0, 2.0, 3.0 };
	PetscInt    i[] = { 0,   1,   2   };

	g_tol  = 1e-12;

	bot   = 1.0;
	top   = 1.5;

	v1 = IntersectTriangularPrism(x, y, z, i, 1.0, bot, top, g_tol);

	bot   = 1.5;
	top   = 3.0;

	v2 = IntersectTriangularPrism(x, y, z, i, 1.0, bot, top, g_tol);


	printf("\n\n\n v1: %f \n\n\n", v1);

	printf("\n\n\n v2: %f \n\n\n", v2);

	printf("\n\n\n sum: %f \n\n\n", v1+v2);

}
*/
//---------------------------------------------------------------------------

