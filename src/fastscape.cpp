#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <algorithm>

#include <stdio.h>
#include <string.h>

// LaMEM header file
#include "LaMEM.h"
#include "surf.h"
#include "interpolate.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "phase.h"
#include "JacRes.h"
#include "tools.h"
#include "Tensor.h"
#include "fastscape.h"
#include "paraViewOutSurf.h"
#include "paraViewOutBin.h"

//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreate(FastScapeLib *FSLib, FB *fb)
{
	PetscInt       maxPhaseID;
	Scaling        *scal;
	FreeSurf       *surf;

	// access context
	surf    = FSLib->surf;
	scal    = FSLib->scal;

	//=================================================================================
	// Load data from .dat file
	//=================================================================================
	// initialize
	// non uniform grid
	FSLib->non_uniform_grid  =   0;
	// 2D grid
	FSLib->fs2D              =   0;
	// extend range & nodes
	FSLib->extendedNodes     =   101;
	FSLib->extendedRange     =   100 * scal->length_fs;
	// refine times & load refined grid
	FSLib->refine            =   1;
	// max timestep
	FSLib->Max_dt            =   0.01 * scal->time_fs;
	// random noise
	FSLib->random_noise      =   1;
	// sedimentation
	FSLib->setMarine         =   0;
	// output information
	FSLib->surf_out_nstep    =   1;
	FSLib->vec_times         =   1;
	// total phases
	maxPhaseID = surf->jr->dbm->numPhases-1;
	// total steps
	FSLib->count_nsteps      = 0;
	// load information from .dat file
	// setup block access mode
	PetscCall(FBFindBlocks(fb, _REQUIRED_, "<FastScapeStart>", "<FastScapeEnd>"));

	if(fb->nblocks)
	{
		//-------------------------------
		// Grid information
		//-------------------------------
		// non uniform grid
		PetscCall(getIntParam(fb, _OPTIONAL_, "non_uniform_grid", &FSLib->non_uniform_grid, 1,  1)); // flag
		// 2D grid
		PetscCall(getIntParam(fb, _OPTIONAL_, "fs2D",             &FSLib->fs2D,             1,  1)); // flag
		if(FSLib->fs2D)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "extendedRange",    &FSLib->extendedRange,    1,  1/scal->length_fs)); // km (LaMEM) -> m (FastScape)

			if (!FSLib->non_uniform_grid)
			{
				PetscCall(getIntParam   (fb, _REQUIRED_, "extendedNodes",    &FSLib->extendedNodes,    1,  10000000)); // non-dimensional
				if(FSLib->extendedNodes <= 2)   SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "extendedNodes must be ≥ 2");
			}
		}
		// refined grid
		PetscCall(getIntParam   (fb, _OPTIONAL_, "fs_refine",        &FSLib->refine,           1,  100)); // non-dimensional

		//===============================
		// FastScape PARAMETER
		//===============================
		// dt & boundary condition
		PetscCall(getScalarParam(fb, _REQUIRED_, "Max_dt",            &FSLib->Max_dt,           1,  1 / scal->time_fs)); // Myr (LaMEM) ->yr (FastScape)

		// bottom-right-top-left; 0 = reflective, 1 = fixed height boundary; When two reflective boundaris face each other they become cyclic
		PetscCall(getStringParam(fb, _REQUIRED_, "topo_boundary",     FSLib->FS_BC,                 "1111"));

		if(strlen(FSLib->FS_BC) != 4)   SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "topo_boundary must be 4-character string (e.g., '1100')");

		// 1 -- boundary velocity == 0; 0 -- boundary velocity from LaMEM
		PetscCall(getStringParam(fb, _REQUIRED_, "vel_boundary",      FSLib->FS_VELBC,              "1111"));

		if(strlen(FSLib->FS_VELBC) != 4)    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "vel_boundary must be 4-character string (e.g., '1100')");

		// random noise
		PetscCall(getIntParam   (fb, _REQUIRED_, "random_noise",      &FSLib->random_noise,     1,  1));

		// sedimentation phase
		PetscCall(getIntParam   (fb, _REQUIRED_, "sed_phases",        &FSLib->sedPhases,        1,  maxPhaseID)); // non-dimensional

		surf->phase = FSLib->sedPhases;
		//-------------------------------
		// Erosion process
		//-------------------------------
		// kf, kd can be set as an array
		PetscCall(getScalarParam(fb, _REQUIRED_, "kf",                &FSLib->kf,               1,          1.0)); // m/yr
		PetscCall(getScalarParam(fb, _REQUIRED_, "kfsed",             &FSLib->kfsed,            1,          1.0)); // m/yr
		PetscCall(getScalarParam(fb, _REQUIRED_, "m",                 &FSLib->m,                1,          1.0)); // non-dimensional
		PetscCall(getScalarParam(fb, _REQUIRED_, "n",                 &FSLib->n,                1,          1.0)); // non-dimensional
		PetscCall(getScalarParam(fb, _REQUIRED_, "kd",                &FSLib->kd,               1,          1.0)); // m/yr
		PetscCall(getScalarParam(fb, _REQUIRED_, "kdsed",             &FSLib->kdsed,            1,          1.0)); // m/yr
		PetscCall(getScalarParam(fb, _REQUIRED_, "g",                 &FSLib->g,                1,          1.0)); // non-dimensional
		PetscCall(getScalarParam(fb, _REQUIRED_, "gsed",              &FSLib->gsed,             1,          1.0)); // non-dimensional
		PetscCall(getScalarParam(fb, _REQUIRED_, "p",                 &FSLib->p,                1,          1.0)); // non-dimensional

		//-------------------------------
		// Sedimentation process
		//-------------------------------
		PetscCall(getIntParam   (fb, _OPTIONAL_, "setMarine",         &FSLib->setMarine,        1,         1)); // flag
		if(FSLib->setMarine)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "sealevel",      &FSLib->sealevel,         1,         1 / scal->length_fs)); // m
			PetscCall(getScalarParam(fb, _REQUIRED_, "poroSilt",      &FSLib->poro_silt,        1,         1.0)); // non-dimensional
			PetscCall(getScalarParam(fb, _REQUIRED_, "poroSand",      &FSLib->poro_sand,        1,         1.0)); // non-dimensional
			PetscCall(getScalarParam(fb, _REQUIRED_, "zporoSilt",     &FSLib->zporo_silt,       1,         1 / scal->length_fs));
			PetscCall(getScalarParam(fb, _REQUIRED_, "zporoSand",     &FSLib->zporo_sand,       1,         1 / scal->length_fs));
			PetscCall(getScalarParam(fb, _REQUIRED_, "ratio",         &FSLib->ratio,            1,         1.0)); // non-dimensional
			PetscCall(getScalarParam(fb, _REQUIRED_, "Lsolve",        &FSLib->Lsolve,           1,         1 / scal->length_fs)); // m
			PetscCall(getScalarParam(fb, _REQUIRED_, "kdsSilt",       &FSLib->kds_silt,         1,         1.0)); // m2/yr
			PetscCall(getScalarParam(fb, _REQUIRED_, "kdsSand",       &FSLib->kds_sand,         1,         1.0)); // m2/yr
		}

		//-------------------------------
		// Output information
		//-------------------------------
		PetscCall(getIntParam   (fb, _OPTIONAL_, "surf_out_nstep",    &FSLib->surf_out_nstep,   1,        1e6)); // non-dimensional
		PetscCall(getScalarParam(fb, _OPTIONAL_, "vec_times",         &FSLib->vec_times,        1,        1.0)); // non-dimensional
	}

	PetscCall(FBFreeBlocks(fb));

	//=================================================================================
	// Load grid information from LaMEM
	//=================================================================================
	// load nx, ny, rangeX, rangeY, rangeZ

	PetscCall(FastScapeLoadGridInf(FSLib));

	if (ISRankZero(PETSC_COMM_WORLD))
	{
		PetscCall(FastScapeCreateSurfaceGrid(FSLib, 1));
	}

	// save nx & ny
	if(FSLib->fs2D)
	{
		if(FSLib->refine == 1)
		{
			FSLib->nx_solve    = (int)FSLib->extendedXNodes;
			FSLib->ny_solve    = (int)FSLib->extendedYNodes;
		}
		else
		{
			FSLib->nx_solve    = (int)FSLib->etRefineXNodes;
			FSLib->ny_solve    = (int)FSLib->etRefineYNodes;
		}
	}
	else
	{
		if(FSLib->refine == 1)
		{
			if(!FSLib->non_uniform_grid)
			{
				FSLib->nx_solve    = (int)FSLib->nx_LaMEM;
				FSLib->ny_solve    = (int)FSLib->ny_LaMEM;
			}
			else
			{
				FSLib->nx_solve    = (int)FSLib->msx_fs.nnodes_nug;
				FSLib->ny_solve    = (int)FSLib->msy_fs.nnodes_nug;
			}
		}
		else
		{
			FSLib->nx_solve    = (int)FSLib->nx_refine;
			FSLib->ny_solve    = (int)FSLib->ny_refine;
		}
	}

	FSLib->nodes_solve = FSLib->nx_solve * FSLib->ny_solve;

	//=================================================================================
	// Visualization
	//=================================================================================
	PetscPrintf(PETSC_COMM_WORLD, "FastScape parameters: \n");
	// LaMEM grid
	PetscPrintf(PETSC_COMM_WORLD, "    Original grid: \n");
	PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%" PetscInt_FMT ", %" PetscInt_FMT "]\n",     FSLib->nx_LaMEM, FSLib->ny_LaMEM);
	// non uniform grid
	if(FSLib->non_uniform_grid)
	{
		PetscPrintf(PETSC_COMM_WORLD, "    Non unifrom grid: \n");
		PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%" PetscInt_FMT ", %" PetscInt_FMT "]\n",     FSLib->msx_fs.nnodes_nug, FSLib->msy_fs.nnodes_nug);
	}
	// 2D extended grid
	if(FSLib->fs2D)
	{
		PetscPrintf(PETSC_COMM_WORLD, "    Extended  grid: \n");
		PetscPrintf(PETSC_COMM_WORLD, "    [rangeX,rangeY]       : [%g, %g] %s\n",   FSLib->extendedXRange / scal->length_fs, FSLib->extendedYRange / scal->length_fs, scal->lbl_length);
		PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%" PetscInt_FMT ", %" PetscInt_FMT "]\n",      FSLib->extendedXNodes, FSLib->extendedYNodes);
	}
	// refined grid
	if(FSLib->refine > 1)
	{
		PetscPrintf(PETSC_COMM_WORLD, "    Refined grid: \n");
		PetscPrintf(PETSC_COMM_WORLD, "    Refined times         : %" PetscInt_FMT "\n",       FSLib->refine);
		// 2D
		if(FSLib->fs2D) PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", FSLib->etRefineXNodes, FSLib->etRefineYNodes);
		// 3D
		else    PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", FSLib->nx_refine, FSLib->ny_refine);
	}

	// surface process parameter
	PetscPrintf(PETSC_COMM_WORLD, "  Surface process: \n");
	PetscPrintf(PETSC_COMM_WORLD, "    Max timestep          : %g %s\n",        FSLib->Max_dt / scal->time_fs, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "    Topography boundary   : %s\n",           FSLib->FS_BC);
	PetscPrintf(PETSC_COMM_WORLD, "    Velocity boundary     : %s\n",           FSLib->FS_VELBC);
	PetscPrintf(PETSC_COMM_WORLD, "    Sedimentation phase   : %" PetscInt_FMT "\n",           FSLib->sedPhases);
	PetscPrintf(PETSC_COMM_WORLD, "    SPL: \n");
	PetscPrintf(PETSC_COMM_WORLD, "      Kf                  : %g\n",           FSLib->kf);
	PetscPrintf(PETSC_COMM_WORLD, "      Kfsed               : %g\n",           FSLib->kfsed);
	PetscPrintf(PETSC_COMM_WORLD, "      m                   : %g\n",           FSLib->m);
	PetscPrintf(PETSC_COMM_WORLD, "      n                   : %g\n",           FSLib->n);
	PetscPrintf(PETSC_COMM_WORLD, "    Hillslope process: \n");
	PetscPrintf(PETSC_COMM_WORLD, "      Kd                  : %g\n",           FSLib->kd);
	PetscPrintf(PETSC_COMM_WORLD, "      Kdsed               : %g\n",           FSLib->kdsed);
	PetscPrintf(PETSC_COMM_WORLD, "      g                   : %g\n",           FSLib->g);
	PetscPrintf(PETSC_COMM_WORLD, "      gsed                : %g\n",           FSLib->gsed);
	PetscPrintf(PETSC_COMM_WORLD, "      p                   : %g\n",           FSLib->p);

	if(FSLib->setMarine)
	{
		PetscPrintf(PETSC_COMM_WORLD, "    Marine process: \n");
		PetscPrintf(PETSC_COMM_WORLD, "      sealevel            : %g %s\n",    FSLib->sealevel / scal->length_fs, scal->lbl_length);
		PetscPrintf(PETSC_COMM_WORLD, "      poro_silt           : %g\n",       FSLib->poro_silt);
		PetscPrintf(PETSC_COMM_WORLD, "      poro_sand           : %g\n",       FSLib->poro_sand);
		PetscPrintf(PETSC_COMM_WORLD, "      zporo_silt          : %g\n",       FSLib->zporo_silt);
		PetscPrintf(PETSC_COMM_WORLD, "      zporo_sand          : %g\n",       FSLib->zporo_sand);
		PetscPrintf(PETSC_COMM_WORLD, "      ratio               : %g\n",       FSLib->ratio);
		PetscPrintf(PETSC_COMM_WORLD, "      L                   : %g\n",       FSLib->Lsolve);
		PetscPrintf(PETSC_COMM_WORLD, "      kds_silt            : %g\n",       FSLib->kds_silt);
		PetscPrintf(PETSC_COMM_WORLD, "      kds_sand            : %g\n",       FSLib->kds_sand);
	}
	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateData(FastScapeLib *FSLib)
{
	PetscFunctionBeginUser;

	FreeSurf  *surf;
	PetscInt  nodes, nodes_nug, nodes_ori;

	// access context
	// Invalid pointer check
	surf      = FSLib->surf;

	// nodes
	nodes     = FSLib->nodes_solve;
	nodes_nug = FSLib->fsX.nodes * FSLib->fsY.nodes;
	nodes_ori = FSLib->nx_LaMEM * FSLib->ny_LaMEM;


	// for collection
	// vx & vy &vz
	PetscCall(DMCreateGlobalVector(surf->DA_SURF, &FSLib->vz_collect));
	PetscCall(DMCreateGlobalVector(surf->DA_SURF, &FSLib->vx_collect));
	PetscCall(DMCreateGlobalVector(surf->DA_SURF, &FSLib->vy_collect));

	if (ISRankZero(PETSC_COMM_WORLD))
	{
		// for create vector
#define FASTSCAPE_ROOT_VEC_CREATE(vec_ptr, size)                       \
            do                                                                  \
            {                                                                   \
                (vec_ptr) = nullptr;                                            \
                if (ISRankZero(PETSC_COMM_WORLD))                               \
                {                                                               \
                    PetscCall(VecCreateSeq(PETSC_COMM_SELF, (size), &(vec_ptr)));                                              \
                }                                                               \
            } while (0)

		// topography
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->gtopo_fs, nodes_ori);

		// for different grid
		// non_uniform_grid
		if(FSLib->non_uniform_grid)
		{
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->gtopo_nug, nodes_nug);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vx_nug, nodes_nug);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vy_nug, nodes_nug);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vz_nug, nodes_nug);
		}

		if(FSLib->fs2D)
		{
			// extended part
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->gtopo_extend, nodes);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vx_extend, nodes);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vy_extend, nodes);
			FASTSCAPE_ROOT_VEC_CREATE(FSLib->vz_extend, nodes);

			if(FSLib->refine > 1)
			{
				// extended part after refinement
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->gtopo_et_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vx_et_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vy_et_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vz_et_refine, nodes);
			}
		}
		else
		{
			if(FSLib->refine > 1)
			{
				// refined part
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->gtopo_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vx_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vy_refine, nodes);
				FASTSCAPE_ROOT_VEC_CREATE(FSLib->vz_refine, nodes);
			}
		}

		// FastScape solution
		// Create all vectors, as the part for creating output is later than the part for creating surf
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->silt_fraction, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->basement, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->total_erosion, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->drainage_area, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->erosion_rate, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->slope, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->curvature, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->chi, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->catchment, nodes);
		FASTSCAPE_ROOT_VEC_CREATE(FSLib->lake_depth, nodes);

#undef FASTSCAPE_ROOT_VEC_CREATE
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeLoadGridInf(FastScapeLib *FSLib)
{
	PetscFunctionBeginUser;

	// load global nx, ny, rangeX, rangeY, rangeZ
	FDSTAG   *fs;
	Scaling  *scal;
	PetscInt baseNodesX, baseNodesY;
	PetscScalar bx, by, bz, ex, ey, ez;

	fs    = FSLib->surf->jr->fs;
	scal  = fs->scal;

	// range X, Y, Z
	PetscCall(FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez));

	// km (LaMEM) in GEO
	FSLib->rangeX_begin = bx * scal->length;
	FSLib->rangeX_end   = ex * scal->length;
	FSLib->rangeY_begin = by * scal->length;
	FSLib->rangeY_end   = ey * scal->length;
	FSLib->rangeZ_begin = bz * scal->length;
	FSLib->rangeZ_end   = ez * scal->length;

	FSLib->rangeX       = (FSLib->rangeX_end - FSLib->rangeX_begin) * scal->length_fs; //(km) in LaMEM to (m) in FastScape (GEO)
	FSLib->rangeY       = (FSLib->rangeY_end - FSLib->rangeY_begin) * scal->length_fs;

	// original nx, ny
	FSLib->nx_LaMEM        = fs->dsx.tnods;
	FSLib->ny_LaMEM        = fs->dsy.tnods;
	FSLib->fsX.nodes    = FSLib->nx_LaMEM;
	FSLib->fsY.nodes    = FSLib->ny_LaMEM;

	// create new nodes & range
	// non uniform grid
	if(FSLib->non_uniform_grid)
	{
		// original grid
		// setting bias-flag, minimum grid spacing, grid nodes
		// x-direction
		PetscCall(FSLoadNonUniformGrid(&FSLib->msx_fs, FSLib->rangeX_end / scal->length, fs->scal));
		// y-direction
		PetscCall(FSLoadNonUniformGrid(&FSLib->msy_fs, FSLib->rangeY_end / scal->length, fs->scal));
		FSLib->fsX.nodes      = FSLib->msx_fs.nnodes_nug;
		FSLib->fsY.nodes      = FSLib->msy_fs.nnodes_nug;
	}

	// 3D refine
	if(FSLib->refine > 1 && FSLib->fs2D == 0)
	{
		baseNodesX = FSLib->non_uniform_grid ? FSLib->fsX.nodes : FSLib->nx_LaMEM;
		baseNodesY = FSLib->non_uniform_grid ? FSLib->fsY.nodes : FSLib->ny_LaMEM;

		// refined grid in FastScaoe
		FSLib->nx_refine = (baseNodesX - 1) * FSLib->refine + 1;
		FSLib->ny_refine = (baseNodesY - 1) * FSLib->refine + 1;
		FSLib->fsX.nodes_refine      = FSLib->nx_refine;
		FSLib->fsY.nodes_refine      = FSLib->ny_refine;
	}

	// 2D
	if(FSLib->fs2D)
	{
		// 2D No Refine
		// extend grid in FastScape
		// extend in rangeX
		if(FSLib->rangeX > FSLib->rangeY)
		{
			FSLib->extendedXRange = FSLib->rangeX;
			FSLib->extendedYRange = FSLib->extendedRange;
			FSLib->extendedXNodes = FSLib->non_uniform_grid ? FSLib->fsX.nodes : FSLib->nx_LaMEM;;
			FSLib->extendedYNodes = FSLib->non_uniform_grid ?
			                        (PetscInt)(FSLib->extendedYRange / scal->length_fs / FSLib->msx_fs.min_spacing) + 2 :
			                        FSLib->extendedNodes;
			FSLib->extendedX      = 0;
			FSLib->extendedY      = 1;
		}
		// extend in rangeY
		else
		{
			FSLib->extendedXRange = FSLib->extendedRange;
			FSLib->extendedYRange = FSLib->rangeY;
			FSLib->extendedXNodes = FSLib->non_uniform_grid ?
			                        (PetscInt)(FSLib->extendedXRange / scal->length_fs / FSLib->msy_fs.min_spacing) + 2 :
			                        FSLib->extendedNodes;
			FSLib->extendedYNodes = FSLib->non_uniform_grid ? FSLib->fsY.nodes : FSLib->ny_LaMEM;;
			FSLib->extendedX      = 1;
			FSLib->extendedY      = 0;
		}
		FSLib->fsX.nodes_extend      = FSLib->extendedXNodes;
		FSLib->fsY.nodes_extend      = FSLib->extendedYNodes;

		// extended grid after refinement in FastScape
		if(FSLib->refine > 1)
		{
			FSLib->etRefineXNodes = (FSLib->fsX.nodes_extend - 1) * FSLib->refine + 1;
			FSLib->etRefineYNodes = (FSLib->fsY.nodes_extend - 1) * FSLib->refine + 1;
			FSLib->fsX.nodes_refine      = FSLib->etRefineXNodes;
			FSLib->fsY.nodes_refine      = FSLib->etRefineYNodes;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FSLoadNonUniformGrid(MeshSeg1DFS *ms_fs, PetscScalar xend, Scaling *scal)
{
	PetscInt i;
	PetscScalar begSz, endSz, avgSz, bias;

	// initialize
	ms_fs->bias = 0;
	ms_fs->xstart[ms_fs->nsegs] = xend;

	// setting bias-flag, minimum grid spacing, grid nodes
	// bias == 1.0 ?
	for(i = 0; i < ms_fs->nsegs; i++)
	{
		if(ms_fs->biases[i] != 1.0)
		{
			ms_fs->bias = 1; break;
		}
	}

	// obtain the minimum grid spacing
	for(i = 0; i < ms_fs->nsegs; i++)
	{
		if(!ms_fs->bias)
		{
			ms_fs->grid_spacing_min[i]  = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (PetscScalar)(ms_fs->istart[i+1] - ms_fs->istart[i]);
			ms_fs->grid_spacing_max[i]  = ms_fs->grid_spacing_min[i];
		}
		else
		{
			// bias (last to first cell size ratio > 1 -> growing)
			bias  = ms_fs->biases[i];

			// average cell size
			avgSz = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (PetscScalar)(ms_fs->istart[i+1] - ms_fs->istart[i]);

			// cell size limits
			begSz = 2.0*avgSz/(1.0 + bias);
			endSz = bias*begSz;

			ms_fs->grid_spacing_min[i] = (begSz < endSz)? begSz : endSz;
			ms_fs->grid_spacing_max[i] = (begSz > endSz)? begSz : endSz;
		}
	}
	ms_fs->min_spacing = *min_element(ms_fs->grid_spacing_min, ms_fs->grid_spacing_min + ms_fs->nsegs);
	ms_fs->max_spacing = *max_element(ms_fs->grid_spacing_max, ms_fs->grid_spacing_max + ms_fs->nsegs);

	// get new spacing & nodes
	ms_fs->nnodes_nug  = (PetscInt)(floor( (ms_fs->xstart[ms_fs->nsegs] - ms_fs->xstart[0]) * scal->length / ms_fs->min_spacing ) ) + 2;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCopyMeshSeg1D(FastScapeLib *FSLib, MeshSeg1D *ms, const char *dir)
{
	PetscInt i;

	if(strcmp("x", dir) == 0)
	{
		FSLib->msx_fs.nsegs             = ms->nsegs;
		FSLib->msx_fs.tcels             = ms->tcels;

		for(i = 0; i < ms->nsegs; i++)
		{
			FSLib->msx_fs.istart[i]     = ms->istart[i];
			FSLib->msx_fs.xstart[i]     = ms->xstart[i];
			FSLib->msx_fs.biases[i]     = ms->biases[i];
		}
		FSLib->msx_fs.istart[ms->nsegs] = ms->istart[ms->nsegs];
	}
	else if(strcmp("y", dir) == 0)
	{
		FSLib->msy_fs.nsegs = ms->nsegs;
		FSLib->msy_fs.tcels             = ms->tcels;

		for(i = 0; i < ms->nsegs; i++)
		{
			FSLib->msy_fs.istart[i]     = ms->istart[i];
			FSLib->msy_fs.xstart[i]     = ms->xstart[i];
			FSLib->msy_fs.biases[i]     = ms->biases[i];
		}
		FSLib->msy_fs.istart[ms->nsegs] = ms->istart[ms->nsegs];
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void GenerateGridCoordinates(PetscScalar *coords, PetscInt numNodes, PetscScalar start, PetscScalar spacing)
{
	for (PetscInt i = 0; i < numNodes; i++)
	{
		coords[i] = start + spacing * (PetscScalar)i;
	}
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateSurfaceGrid(FastScapeLib *FSLib, PetscInt mode)
{
	//mode = 1, normal run (create), = 2, restart, = 0, update range
	FSGrid   *fsX;
	FSGrid   *fsY;
	Scaling  *scal;
	TSSol    *ts;
	JacRes   *jr;
	PetscInt i, step_fs;

	fsX     = &FSLib->fsX;
	fsY     = &FSLib->fsY;
	scal    = FSLib->scal;
	jr      = FSLib->jr;
	ts      = jr->ts;
	step_fs = ts->istep;

	fsX->dx = FSLib->rangeX / scal->length_fs / (PetscScalar)(fsX->nodes - 1);
	fsY->dx = FSLib->rangeY / scal->length_fs / (PetscScalar)(fsY->nodes - 1);

	// original fastscape grid
	if(mode == 1 || mode == 2)
	{
		PetscCall(PetscMalloc((size_t)(fsX->nodes) * sizeof(PetscScalar), &fsX->ncoor));
		PetscCall(PetscMalloc((size_t)(fsY->nodes) * sizeof(PetscScalar), &fsY->ncoor));
	}

	GenerateGridCoordinates(fsX->ncoor, fsX->nodes, FSLib->rangeX_begin, fsX->dx);
	GenerateGridCoordinates(fsY->ncoor, fsY->nodes, FSLib->rangeY_begin, fsY->dx);

	if(FSLib->fs2D == 0)
	{
		// 3D refined grid
		if(FSLib->refine > 1)
		{
			if(mode == 2 || step_fs == 0)
			{
				PetscCall(PetscMalloc((size_t)(fsX->nodes_refine) * sizeof(PetscScalar), &fsX->ncoor_refine));
				PetscCall(PetscMalloc((size_t)(fsY->nodes_refine) * sizeof(PetscScalar), &fsY->ncoor_refine));
			}

			fsX->dx_refine = FSLib->rangeX / scal->length_fs  / (PetscScalar)(fsX->nodes_refine - 1);
			fsY->dx_refine = FSLib->rangeY / scal->length_fs  / (PetscScalar)(fsY->nodes_refine - 1);

			GenerateGridCoordinates(fsX->ncoor_refine, fsX->nodes_refine, FSLib->rangeX_begin, fsX->dx_refine);
			GenerateGridCoordinates(fsY->ncoor_refine, fsY->nodes_refine, FSLib->rangeY_begin, fsY->dx_refine);
		}
	}
	else
	{
		// 2D grid
		if(mode == 2 || step_fs == 0)
		{
			PetscCall(PetscMalloc((size_t)(fsX->nodes_extend) * sizeof(PetscScalar), &fsX->ncoor_extend));
			PetscCall(PetscMalloc((size_t)(fsY->nodes_extend) * sizeof(PetscScalar), &fsY->ncoor_extend));
		}

		fsX->dx_extend = FSLib->extendedXRange / scal->length_fs  / (PetscScalar)(fsX->nodes_extend - 1);
		fsY->dx_extend = FSLib->extendedYRange / scal->length_fs  / (PetscScalar)(fsY->nodes_extend - 1);

		GenerateGridCoordinates(fsX->ncoor_extend, fsX->nodes_extend, FSLib->rangeX_begin, fsX->dx_extend);
		GenerateGridCoordinates(fsY->ncoor_extend, fsY->nodes_extend, FSLib->rangeY_begin, fsY->dx_extend);

		// 2D refined grid
		if(FSLib->refine > 1)
		{
			if(mode == 2 || step_fs == 0)
			{
				PetscCall(PetscMalloc((size_t)(fsX->nodes_refine) * sizeof(PetscScalar), &fsX->ncoor_refine));
				PetscCall(PetscMalloc((size_t)(fsY->nodes_refine) * sizeof(PetscScalar), &fsY->ncoor_refine));
			}

			fsX->dx_refine = FSLib->extendedXRange / scal->length_fs  / (PetscScalar)(fsX->nodes_refine - 1);
			fsY->dx_refine = FSLib->extendedYRange / scal->length_fs  / (PetscScalar)(fsY->nodes_refine - 1);

			GenerateGridCoordinates(fsX->ncoor_refine, fsX->nodes_refine, FSLib->rangeX_begin, fsX->dx_refine);
			GenerateGridCoordinates(fsY->ncoor_refine, fsY->nodes_refine, FSLib->rangeY_begin, fsY->dx_refine);
		}
	}

	if(FSLib->non_uniform_grid == 1)
	{
		// create a global coordinate for LaMEM
		if(FSLib->fs2D == 0)
		{
			if(mode == 2 || step_fs == 0)
			{
				PetscCall(PetscMalloc((size_t)((FSLib->nx_LaMEM + 1)) * sizeof(PetscScalar), &FSLib->ncoor_ori_x));
				PetscCall(PetscMalloc((size_t)((FSLib->ny_LaMEM + 1)) * sizeof(PetscScalar), &FSLib->ncoor_ori_y));
			}

			// x-direction
			for(i = 0; i < FSLib->msx_fs.nsegs; i++)
			{
				// generate nodal coordinates for the local part of the segment
				PetscCall(FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal));
			}
			// last nodes
			FSLib->ncoor_ori_x[FSLib->nx_LaMEM] = FSLib->rangeX_end;

			// y-direction
			for(i = 0; i < FSLib->msy_fs.nsegs; i++)
			{
				// generate nodal coordinates for the local part of the segment
				PetscCall(FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal));
			}
			// last nodes
			FSLib->ncoor_ori_y[FSLib->ny_LaMEM] = FSLib->rangeY_end;
		}
		else
		{
			if(FSLib->extendedX == 1)
			{
				if(mode == 2 || step_fs == 0 )  PetscCall(PetscMalloc((size_t)((FSLib->ny_LaMEM + 1)) * sizeof(PetscScalar), &FSLib->ncoor_ori_y));

				// y-direction
				for(i = 0; i < FSLib->msy_fs.nsegs; i++)
				{
					// generate nodal coordinates for the local part of the segment
					PetscCall(FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal));
				}
				// last nodes
				FSLib->ncoor_ori_y[FSLib->ny_LaMEM] = FSLib->rangeY_end;
			}
			else
			{
				if(mode == 2|| step_fs == 0)    PetscCall(PetscMalloc((size_t)((FSLib->nx_LaMEM + 1)) * sizeof(PetscScalar), &FSLib->ncoor_ori_x));

				// x-direction
				for(i = 0; i < FSLib->msx_fs.nsegs; i++)
				{
					// generate nodal coordinates for the local part of the segment
					PetscCall(FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal));
				}
				// last nodes
				FSLib->ncoor_ori_x[FSLib->nx_LaMEM] = FSLib->rangeX_end;
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ScalingFastScapeCreate(Scaling *scal)
{
	PetscFunctionBeginUser;
	PetscScalar km, yr, m;

	// read unit values
	km     = 1e3;
	m      = 1e2;
	yr     = 3600.0*24.0*365.25;

	if(scal->utype == _SI_)
	{
		// s (LaMEM) -> yr (FastScape)
		scal->time_fs             = 1/yr;                      snprintf(scal->lbl_time_fs, sizeof(scal->lbl_time_fs),        "[yr]");
		// m (LaMEM) -> m (FastSCape)
		scal->length_fs           = 1.0;                       snprintf(scal->lbl_length_fs, sizeof(scal->lbl_length_fs),    "[m]");
		// m/s (LaMEM) -> m/yr (FastScape)
		scal->velocity_fs         = 1/yr;                      snprintf(scal->lbl_velocity_fs, sizeof(scal->lbl_velocity_fs),   "[m/yr]");

		// output
		// m^2 --> m^2
		scal->area_fs             = 1.0;                       snprintf(scal->lbl_area_fs, sizeof(scal->lbl_area_fs),       "[m^2]");
		// m/yr --> m/yr
		scal->rate                = 1.0;                       snprintf(scal->lbl_rate, sizeof(scal->lbl_rate),        "[m/yr]");
		scal->degree              = 1.0;                       snprintf(scal->lbl_degree, sizeof(scal->lbl_degree),       "[°]");
	}
	else if(scal->utype == _GEO_)
	{
		// Myr (LaMEM) -> yr (FastScape)
		scal->time_fs             = 1e6;                       snprintf(scal->lbl_time_fs, sizeof(scal->lbl_time_fs),     "[yr]");
		// km (LaMEM) -> m (FastScape)
		scal->length_fs           = km;                        snprintf(scal->lbl_length_fs, sizeof(scal->lbl_length_fs),     "[m]");
		// cm/yr (LaMEM) -> m/yr (FastScape)
		scal->velocity_fs         = m;                         snprintf(scal->lbl_velocity_fs, sizeof(scal->lbl_velocity_fs),    "[m/yr]");

		// output
		// m^2 --> km^2
		scal->area_fs             = km * km;                 snprintf(scal->lbl_area_fs, sizeof(scal->lbl_area_fs),       "[km^2]");
		// m/yr --> km/yr
		scal->rate                = 1.0 * km;                    snprintf(scal->lbl_rate, sizeof(scal->lbl_rate),         "[km/yr]");
		scal->degree              = 1.0;                       snprintf(scal->lbl_degree, sizeof(scal->lbl_degree),      "[°]");
	}
	else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect unit type for FastScape");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfFastScapeCreate(FastScapeLib *FSLib, FB *fb)
{

	char filename[_str_len_];

	PetscFunctionBeginUser;

	// initialize
	// topography & pvd
	FSLib->outsurf_fs          = 1;
	FSLib->outpvd_fs           = 1;
	FSLib->out_topofs          = 1;

	// surface parameter
	FSLib->out_silt_fraction   = 1;
	FSLib->out_basement        = 1;
	FSLib->out_total_erosion   = 1;
	FSLib->out_drainage_area   = 1;
	FSLib->out_erosion_rate    = 1;
	FSLib->out_slope           = 1;
	FSLib->out_curvature       = 1;
	FSLib->out_chi             = 1;
	FSLib->out_catchment       = 1;
	FSLib->out_lake_depth      = 1;

	// check activation
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_fs",            &FSLib->outsurf_fs,        1, 1));

	// read
	PetscCall(getStringParam(fb, _OPTIONAL_, "out_file_name",          filename,               "output"));
	// pvd
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_fs_pvd",             &FSLib->outpvd_fs,         1, 1));

	// topography
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_topofs",        &FSLib->out_topofs,            1, 1));

	// surface parameter
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_silt_fraction", &FSLib->out_silt_fraction,     1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_basement",      &FSLib->out_basement,          1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_total_erosion", &FSLib->out_total_erosion,     1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_drainage_area", &FSLib->out_drainage_area,     1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_erosion_rate",  &FSLib->out_erosion_rate,      1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_slope",         &FSLib->out_slope,             1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_curvature",     &FSLib->out_curvature,         1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_chi",           &FSLib->out_chi,               1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_catchment",     &FSLib->out_catchment,         1, 1));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_surf_lake_depth",    &FSLib->out_lake_depth,        1, 1));

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "FastScape output parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Write .pvd fs file         : %s \n", FSLib->outpvd_fs    ? "yes" : "no");
	if(FSLib->outpvd_fs)
	{
		if(FSLib->out_topofs)        PetscPrintf(PETSC_COMM_WORLD, "     Topo                       @ \n");
		if(FSLib->out_silt_fraction) PetscPrintf(PETSC_COMM_WORLD, "     silt_fraction              @ \n");
		if(FSLib->out_basement)      PetscPrintf(PETSC_COMM_WORLD, "     basement                   @ \n");
		if(FSLib->out_drainage_area) PetscPrintf(PETSC_COMM_WORLD, "     drainage_area              @ \n");
		if(FSLib->out_erosion_rate)  PetscPrintf(PETSC_COMM_WORLD, "     erosion_rate               @ \n");
		if(FSLib->out_slope)         PetscPrintf(PETSC_COMM_WORLD, "     slope                      @ \n");
		if(FSLib->out_curvature)     PetscPrintf(PETSC_COMM_WORLD, "     curvature                  @ \n");
		if(FSLib->out_chi)           PetscPrintf(PETSC_COMM_WORLD, "     chi                        @ \n");
		if(FSLib->out_catchment)     PetscPrintf(PETSC_COMM_WORLD, "     catchment                  @ \n");
		if(FSLib->out_lake_depth)    PetscPrintf(PETSC_COMM_WORLD, "     lake_depth                 @ \n");
	}
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// set file name
	snprintf(FSLib->outfile_fs,   sizeof(FSLib->outfile_fs),     "%s_fs",        filename);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCopyVelocity(FastScapeLib *FSLib)
{
	JacRes *jr;
	FreeSurf *surf;

	surf = FSLib->surf;
	jr = surf->jr;

	PetscCall(FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx));
	PetscCall(FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy));
	PetscCall(FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateGlobalGrid(PetscScalar *ncoor, MeshSeg1DFS ms_fs, PetscInt iseg, Scaling *scal)
{
	PetscInt    i, M, sum = 0, istart;
	PetscScalar xstart, xclose, bias, avgSz, begSz, endSz, dx;

	PetscFunctionBeginUser;

	// start index
	istart = ms_fs.istart[iseg];

	// total number of cells
	M = ms_fs.istart[iseg+1] - ms_fs.istart[iseg];

	// starting & closing coordinates
	xstart = ms_fs.xstart[iseg] * scal->length;
	xclose = ms_fs.xstart[iseg+1] * scal->length;

	// bias (last to first cell size ratio > 1 -> growing)
	bias = ms_fs.biases[iseg];

	// average cell size
	avgSz = (xclose - xstart)/(PetscScalar)M;

	// uniform case
	if(bias == 1.0)
	{
		// generate coordinates of local nodes
		for(i = istart; i < M + istart + 1; i++)
		{
			ncoor[i] = xstart + (PetscScalar)(i - istart) * avgSz;
		}
	}
	// non-uniform case
	else
	{
		// cell size limits
		begSz = 2.0 * avgSz / (1.0 + bias);
		endSz = bias * begSz;

		// cell size increment (negative for bias < 1)
		dx = (endSz - begSz) / (PetscScalar)(M-1);

		// generate coordinates of local nodes
		for(i = istart; i < istart + M; i++)
		{
			ncoor[i] = xstart + (PetscScalar)(i - istart) * begSz + dx * (PetscScalar)sum ;
			sum += i - istart;
		}
	}

	// override last node coordinate
	ncoor[ ms_fs.istart[iseg+1] ] = xclose;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
GridIndex FindIndexForInterpolation(FastScapeLib* FSLib, PetscScalar x_coor, PetscScalar y_coor, PetscScalar dx, PetscScalar dy)
{
	GridIndex idx = {0, 0, 0, 0, false};
	PetscInt p, q;
	PetscInt find_indx = 0, find_indy = 0;
	PetscScalar x1, x2, y1, y2;
	bool found_x = false, found_y = false;

	// x-direction
	for(p = 0; p < FSLib->nx_LaMEM; )
	{
		x1 = FSLib->ncoor_ori_x[p];
		x2 = FSLib->ncoor_ori_x[p + 1];

		if((x_coor >= x1) && (x_coor <= x2))
		{
			idx.m = p;
			idx.mm = p + 1;

			if(FSLib->nx_LaMEM == idx.mm) idx.mm -= 1;

			found_x = true;
			break;
		}

		find_indx = (PetscInt)(floor((x_coor - x1) / dx));

		if(find_indx > 0) p += find_indx;
		else p++;
	}

	if(found_x)
	{
		// y_direction
		for(q = 0; q < FSLib->ny_LaMEM; )
		{
			y1 = FSLib->ncoor_ori_y[q];
			y2 = FSLib->ncoor_ori_y[q + 1];

			if((y_coor >= y1) && (y_coor <= y2))
			{
				idx.n = q;
				idx.nn = q + 1;

				if(FSLib->ny_LaMEM == idx.nn) idx.nn -= 1;

				found_y = true;
				break;
			}

			find_indy = (PetscInt)(floor((y_coor - y1) / dy));

			if(find_indy > 0) q += find_indy;
			else q++;
		}
	}

	idx.found = found_x && found_y;
	return idx;
}
//---------------------------------------------------------------------------
PetscScalar ReturnBiInterFunction(PetscScalar a, PetscScalar b, PetscScalar c, PetscScalar d,
                                  PetscScalar dx, PetscScalar dy, PetscScalar wtx, PetscScalar wty,
                                  PetscInt m, PetscInt mm, PetscInt n, PetscInt nn)
{
	// Bilinear interpolation
	if (m == mm && n == nn)  return a;
	if (m == mm)  return ((dy - wty) / dy) * a + (wty / dy) * c;
	if (n == nn)  return ((dx - wtx) / dx) * a + (wtx / dx) * b;
	return ((dy - wty) / dy) * ((dx - wtx) / dx) * a +
	       ((dy - wty) / dy) * (wtx / dx) * b +
	       (wty / dy) * ((dx - wtx) / dx) * c +
	       (wty / dy) * (wtx / dx) * d;
}
//---------------------------------------------------------------------------
PetscErrorCode InterpolationFor3DNonUniformGrid(FastScapeLib *FSLib, PetscScalar *value, PetscInt mode)
{
	FSGrid    *fsX;
	FSGrid    *fsY;
	PetscInt i, j, ind, ind_a, ind_b, ind_c, ind_d;
	PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
	PetscScalar *value_save;
	Vec target_vec;

	fsX     = &FSLib->fsX;
	fsY     = &FSLib->fsY;

	switch (mode)
	{
	case 1: target_vec = FSLib->gtopo_nug; break;
	case 2: target_vec = FSLib->vx_nug; break;
	case 3: target_vec = FSLib->vy_nug; break;
	case 4: target_vec = FSLib->vz_nug; break;
	default:
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
		        "Invalid interpolation mode %" PetscInt_FMT, mode);
	}
	PetscCall(VecGetArray(target_vec, &value_save));

	for(j = 0; j < fsY->nodes; j++)
	{
		for(i = 0; i < fsX->nodes; i++)
		{
			x_coor = fsX->ncoor[i];
			y_coor = fsY->ncoor[j];

			GridIndex indices = FindIndexForInterpolation(FSLib, x_coor, y_coor, FSLib->msx_fs.max_spacing, FSLib->msy_fs.max_spacing);

			if(indices.found)
			{
				// interpolate
				ind = j * fsX->nodes + i;
				ind_a = indices.n * FSLib->nx_LaMEM + indices.m;
				ind_b = indices.n * FSLib->nx_LaMEM + indices.mm;
				ind_c = indices.nn * FSLib->nx_LaMEM + indices.m;
				ind_d = indices.nn * FSLib->nx_LaMEM + indices.mm;

				// grid
				x1 = FSLib->ncoor_ori_x[indices.m];
				x2 = FSLib->ncoor_ori_x[indices.mm];
				y1 = FSLib->ncoor_ori_y[indices.n];
				y2 = FSLib->ncoor_ori_y[indices.nn];

				// bilinear interpolation
				dx    = x2 - x1;
				dy    = y2 - y1;
				wtx   = x_coor - x1;
				wty   = y_coor - y1;

				value_save[ind] = ReturnBiInterFunction(
				                      value[ind_a], value[ind_b], value[ind_c], value[ind_d],
				                      dx, dy, wtx, wty, indices.m, indices.mm, indices.n, indices.nn);
			}
		}
	}

	PetscCall(VecRestoreArray(target_vec, &value_save));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode InterpolationFor2DNonUniformGrid(FastScapeLib *FSLib, PetscScalar *value_in, PetscScalar *value_out)
{
	FSGrid *fsGrid;
	PetscReal *coords_ori, *coords_ext;
	PetscInt nodes_ori, nodes_ext;
	PetscFunctionBegin;
	PetscInt i;
	PetscReal x_min, x_max;

	// interpolation direction
	if (FSLib->extendedX)
	{
		fsGrid = &FSLib->fsY;
		coords_ori = FSLib->ncoor_ori_y;
		nodes_ori = FSLib->ny_LaMEM;
		coords_ext = fsGrid->ncoor_extend;
		nodes_ext = fsGrid->nodes_extend;
	}
	else
	{
		fsGrid = &FSLib->fsX;
		coords_ori = FSLib->ncoor_ori_x;
		nodes_ori = FSLib->nx_LaMEM;
		coords_ext = fsGrid->ncoor_extend;
		nodes_ext = fsGrid->nodes_extend;
	}

	// original range
	x_min = coords_ori[0];
	x_max = coords_ori[nodes_ori-1];

	for(i = 0; i < nodes_ext; i++)
	{
		PetscReal x_coor = coords_ext[i];
		PetscInt idx_left = -1;
		PetscInt idx_right = -1;

		// boundary point
		if (x_coor <= x_min)
		{
			value_out[i] = value_in[0];
			continue;
		}
		else if (x_coor >= x_max)
		{
			value_out[i] = value_in[nodes_ori-1];
			continue;
		}

		// find nearest nodes
		PetscInt low = 0;
		PetscInt high = nodes_ori - 1;

		while (low <= high)
		{
			PetscInt mid = low + (high - low) / 2;

			if (coords_ori[mid] <= x_coor)
			{
				idx_left = mid;
				low = mid + 1;
			}
			else
			{
				high = mid - 1;
			}
		}

		if (idx_left >= 0 && idx_left < (nodes_ori - 1))
		{
			idx_right = idx_left + 1;
			PetscReal x1 = coords_ori[idx_left];
			PetscReal x2 = coords_ori[idx_right];
			PetscReal dx = x2 - x1;

			// linear interploation
			PetscReal weight = (x_coor - x1) / dx;
			value_out[i] = (1.0 - weight) * value_in[idx_left] + weight * value_in[idx_right];
		}
		else
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,
			        "Interpolation failed at point %" PetscInt_FMT " (coord=%.4f) between nodes [%.4f, %.4f]",
			        i, x_coor, x_min, x_max);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode GatherVariableFromLaMEM(FastScapeLib *FSLib, PetscScalar *topo_alloc, PetscScalar *vx_alloc, PetscScalar *vy_alloc, PetscScalar *vz_alloc, PetscInt step_fs)
{
	PetscInt L, sx, sy, sz, nx, ny, nz, tproc, i, j;
	PetscMPIInt rankZ_id;
	PetscScalar ***topo;
	PetscScalar ***vz, ***vx, ***vy;
	PetscScalar ***vz_collect, ***vx_collect, ***vy_collect;
	PetscScalar *topo_collect = nullptr;
	PetscScalar *vz_LaMEM = nullptr, *vx_LaMEM = nullptr, *vy_LaMEM = nullptr;
	FDSTAG      *fs;
	FreeSurf    *surf;
	VecScatter ctx_topo = nullptr;
	VecScatter ctx_vx   = nullptr;
	VecScatter ctx_vy   = nullptr;
	VecScatter ctx_vz   = nullptr;

	surf   = FSLib->surf;
	fs     = surf->jr->fs;

	// local process info
	L    = (PetscInt)fs->dsz.rank;
	PetscCall(DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz));

	// Gather topography and velocity
	// topography
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo));
	if(step_fs == 0)
	{
		PetscCall(VecScatterCreateToZero(surf->gtopo, &ctx_topo, &FSLib->gtopo_collect));
		PetscCall(VecScatterBegin(ctx_topo, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD));
		PetscCall(VecScatterEnd(ctx_topo, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD));
	}

	// velocity
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->vx, &vx));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->vy, &vy));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->vz, &vz));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect));

	// load local value to global varibles
	START_PLANE_LOOP
	{
		vx_collect[L][j][i] = vx[L][j][i]; // (km)
		vy_collect[L][j][i] = vy[L][j][i]; // (km)
		vz_collect[L][j][i] = vz[L][j][i]; // (km)
	}
	END_PLANE_LOOP

	// note: if the vec is a local vec, it will create a local output, which isn't correct
	// A global vector is needed
	PetscCall(VecScatterCreateToZero(FSLib->vx_collect, &ctx_vx, &FSLib->vx_LaMEM));
	PetscCall(VecScatterBegin(ctx_vx, FSLib->vx_collect, FSLib->vx_LaMEM, INSERT_VALUES, SCATTER_FORWARD));
	PetscCall(VecScatterEnd(ctx_vx, FSLib->vx_collect, FSLib->vx_LaMEM, INSERT_VALUES, SCATTER_FORWARD));
	PetscCall(VecScatterCreateToZero(FSLib->vy_collect, &ctx_vy, &FSLib->vy_LaMEM));
	PetscCall(VecScatterBegin(ctx_vy, FSLib->vy_collect, FSLib->vy_LaMEM, INSERT_VALUES, SCATTER_FORWARD));
	PetscCall(VecScatterEnd(ctx_vy, FSLib->vy_collect, FSLib->vy_LaMEM, INSERT_VALUES, SCATTER_FORWARD));
	PetscCall(VecScatterCreateToZero(FSLib->vz_collect, &ctx_vz, &FSLib->vz_LaMEM));
	PetscCall(VecScatterBegin(ctx_vz, FSLib->vz_collect, FSLib->vz_LaMEM, INSERT_VALUES, SCATTER_FORWARD));
	PetscCall(VecScatterEnd(ctx_vz, FSLib->vz_collect, FSLib->vz_LaMEM, INSERT_VALUES, SCATTER_FORWARD));

	// reallocate
	tproc = fs->dsx.nproc * fs->dsy.nproc * fs->dsz.nproc;
	rankZ_id = NO_NEED;

	if (fs->dsz.rank == 0)  PetscCall(PetscMPIIntCast(fs->dsx.nproc * fs->dsy.rank + fs->dsx.rank, &rankZ_id));

	// local process info
	ProcInfo local_info = {sx, sy, sz, nx, ny, nz, rankZ_id};
	ProcInfo *global_infos = nullptr;

	// allocate memory
	if (ISRankZero(PETSC_COMM_WORLD))
	{
		PetscCall(PetscMalloc((size_t)(tproc) * sizeof(ProcInfo), &global_infos));
		if (!global_infos)  SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_MEM, "Memory allocation failed");
	}

	// collect information for the surface grid
	PetscCall(MPI_Gather(&local_info, sizeof(ProcInfo), MPI_BYTE, global_infos, sizeof(ProcInfo), MPI_BYTE, 0, PETSC_COMM_WORLD));
	if (ISRankZero(PETSC_COMM_WORLD))
	{
		if(step_fs == 0)    PetscCall(VecGetArray(FSLib->gtopo_collect, &topo_collect));

		PetscCall(VecGetArray(FSLib->vx_LaMEM,  &vx_LaMEM));
		PetscCall(VecGetArray(FSLib->vy_LaMEM,  &vy_LaMEM));
		PetscCall(VecGetArray(FSLib->vz_LaMEM,  &vz_LaMEM));

		// Sort the process information in logical order
		std::vector<ProcInfo> valid_infos;

		for (i = 0; i < tproc; i++)
		{
			if (global_infos[i].rankZ_id != NO_NEED)    valid_infos.push_back(global_infos[i]);
		}

		std::sort(valid_infos.begin(), valid_infos.end(), [](const ProcInfo &a, const ProcInfo &b)
		{
			if (a.sy != b.sy) return a.sy < b.sy;
			return a.sx < b.sx;
		});

		// put the information into target array
		PetscInt global_index = 0;
		for (const auto &pi : valid_infos)
		{
			for (PetscInt y = pi.sy; y < pi.sy + pi.ny; y++)
			{
				for (PetscInt x = pi.sx; x < pi.sx + pi.nx; x++)
				{
					PetscInt target_idx = y * FSLib->nx_LaMEM + x;

					if (step_fs == 0)   topo_alloc[target_idx] = topo_collect[global_index];
					vx_alloc[target_idx] = vx_LaMEM[global_index];
					vy_alloc[target_idx] = vy_LaMEM[global_index];
					vz_alloc[target_idx] = vz_LaMEM[global_index];

					global_index++;
				}
			}
		}

		PetscFree(global_infos);

		if(step_fs == 0)
		{
			PetscCall(VecRestoreArray(FSLib->gtopo_collect, &topo_collect));
			PetscCall(VecDestroy(&FSLib->gtopo_collect));
		}
		PetscCall(VecRestoreArray(FSLib->vz_LaMEM, &vz_LaMEM));
		PetscCall(VecRestoreArray(FSLib->vx_LaMEM, &vx_LaMEM));
		PetscCall(VecRestoreArray(FSLib->vy_LaMEM, &vy_LaMEM));
		PetscCall(VecDestroy(&FSLib->vx_LaMEM));
		PetscCall(VecDestroy(&FSLib->vy_LaMEM));
		PetscCall(VecDestroy(&FSLib->vz_LaMEM));
	}

	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->vx, &vx));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->vy, &vy));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->vz, &vz));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect));

	PetscCall(VecScatterDestroy(&ctx_topo));
	PetscCall(VecScatterDestroy(&ctx_vx));
	PetscCall(VecScatterDestroy(&ctx_vy));
	PetscCall(VecScatterDestroy(&ctx_vz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BilinearInterpolate(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_refine,
                                   PetscScalar *data_refine_pass, Scaling *scal, PetscInt corMode, PetscInt nx_refine, PetscInt ny_refine)
{
	//PetscInt interpolationMode = 1;
	// linear interpolation
	// refine = 1; no nodes; 2, add a node between two original nearest nodes; 3, add two nodes between two nearest original nodes;
	// corMode: 1 -- topography, (km) in LaMEM to (m) in FastScape; 2 -- velocity, (cm/yr) in LaMEM to (m/yr) in FastScape
	PetscScalar unit_factor = 1.0;

	if (corMode == 1)
	{
		unit_factor = scal->length * scal->length_fs;
	}
	else if (corMode == 2)
	{
		unit_factor = scal->velocity / scal->velocity_fs;
	}

	const PetscInt nx_source = (nx_refine - 1) / FSLib->refine + 1;
	const PetscInt ny_source = (ny_refine - 1) / FSLib->refine + 1;

	for (PetscInt j = 0; j < ny_refine; ++j)
	{
		const PetscInt j0 = j / FSLib->refine;
		const PetscInt j1 = PetscMin(j0 + 1, ny_source - 1);
		const PetscScalar ty =
		    (PetscScalar)(j % FSLib->refine) / (PetscScalar)(FSLib->refine);

		for (PetscInt i = 0; i < nx_refine; ++i)
		{
			const PetscInt i0 = i / FSLib->refine;
			const PetscInt i1 = PetscMin(i0 + 1, nx_source - 1);
			const PetscScalar tx =
			    (PetscScalar)(i % FSLib->refine) / (PetscScalar)(FSLib->refine);

			const PetscScalar q00 = data[j0 * nx_source + i0] * unit_factor;
			const PetscScalar q10 = data[j0 * nx_source + i1] * unit_factor;
			const PetscScalar q01 = data[j1 * nx_source + i0] * unit_factor;
			const PetscScalar q11 = data[j1 * nx_source + i1] * unit_factor;

			const PetscInt out = j * nx_refine + i;

			data_refine[out] =
			    (1.0 - tx) * (1.0 - ty) * q00 +
			    tx         * (1.0 - ty) * q10 +
			    (1.0 - tx) * ty         * q01 +
			    tx         * ty         * q11;

			data_refine_pass[out] = data_refine[out];
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Extended2D(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_extended, PetscScalar *data_extend_pass, Scaling *scal, PetscInt corMode, PetscBool isRefine)
{
	PetscScalar sum, unit_factor    = 1.0;
	PetscInt i, j, idx, linearIdx;
	PetscScalar *data_aver_ori = nullptr;
	PetscScalar *data_aver     = nullptr;
	bool extendX        = (FSLib->extendedX == 1);
	PetscInt origDim    = extendX ? FSLib->ny_LaMEM : FSLib->nx_LaMEM;
	PetscInt targetDim  = extendX ? FSLib->extendedYNodes : FSLib->extendedXNodes;
	PetscInt otherDim   = extendX ? FSLib->nx_LaMEM : FSLib->ny_LaMEM;

	// scaling factor
	if(!isRefine)
	{
		switch (corMode)
		{
		case 1: // topography: (LaMEM) → m (FastScape)
			unit_factor = scal->length * scal->length_fs;
			break;
		case 2: // velocity：cm/yr (LaMEM) → m/yr (FastScape)
			unit_factor = scal->velocity / scal->velocity_fs;
			break;
		default: // no transtition
			unit_factor = 1.0;
		}
	}

	PetscCall(PetscMalloc((size_t)(origDim) * sizeof(PetscScalar), &data_aver_ori));
	PetscCall(PetscMalloc((size_t)(targetDim) * sizeof(PetscScalar), &data_aver));

	// calculate averange values with scaling
	for (i = 0; i < origDim; i++)
	{
		sum = 0.0;

		for (j = 0; j < otherDim; j++)
		{
			idx = extendX ? i * FSLib->nx_LaMEM + j : j * FSLib->nx_LaMEM + i;
			sum += data[idx] * unit_factor;
		}

		data_aver_ori[i] = sum / (PetscScalar)(otherDim);
	}

	// non-uniform grid
	if (FSLib->non_uniform_grid)
	{
		PetscCall(InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver));
	}
	else
	{
		for (i = 0; i < origDim; i++)
		{
			data_aver[i] = data_aver_ori[i];
		}
	}

	for (j = 0; j < FSLib->extendedYNodes; j++)
	{
		for (i = 0; i < FSLib->extendedXNodes; i++)
		{
			linearIdx = j * FSLib->extendedXNodes + i;

			if (extendX)   data_extended[linearIdx] = data_aver[j];
			else   data_extended[linearIdx] = data_aver[i];

			if(!isRefine)   data_extend_pass[linearIdx] = data_extended[linearIdx];
		}
	}

	PetscCall(PetscFree(data_aver_ori));
	PetscCall(PetscFree(data_aver));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeRun(FastScapeLib *FSLib)
{
	/* Unit
	FastScape:
	model range (rangeX, rangeY): m
	timestep: yr
	velocity: m/yr
	topography: m

	LaMEM: when using Geo unit
	model range : km
	timestep after scalling: Myr
	velocity: after scalling: cm/yr
	topography: km
	*/

	// Apply surface process to the internal free surface of the model
	// free surface cases only
	FreeSurf *surf;
	surf = FSLib->surf;
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	PetscScalar dt, dt_scal, dt_fs, time_fs;
	PetscInt i, j, tnodes, tnodes_ori, L, sx, sy, sz, nx, ny, nz, step_fs;
	PetscScalar *topo_fs = nullptr, *topo_solve = nullptr, *topo_solve_refine = nullptr, *topo_alloc = nullptr, *gtopo_rank0 = nullptr;
	PetscScalar *vx_alloc = nullptr, *vy_alloc = nullptr, *vz_alloc = nullptr;
	PetscScalar *vx_solve = nullptr, *vy_solve = nullptr, *vz_solve = nullptr;
	PetscScalar *vx_solve_refine = nullptr, *vy_solve_refine = nullptr, *vz_solve_refine = nullptr;
	PetscScalar ***topo;

	// load global nx, ny, dt, time, rangeX, rangeY
	// nx, ny, dt, time
	FDSTAG   *fs;
	TSSol    *ts;
	Scaling  *scal;
	JacRes   *jr;

	jr      = surf->jr;
	fs      = jr->fs;
	ts      = jr->ts;
	scal    = ts->scal;

	// load time
	dt      = ts->dt;
	dt_scal = dt * scal->time;
	dt_fs   = dt_scal * scal->time_fs; // (Myr) in LaMEM to (yr) in FastScape (GEO)
	time_fs = ts->time * scal->time + dt_scal; // time after finishing surface processes
	step_fs = ts->istep;

	// Gather topography and velocity
	tnodes_ori    = FSLib->nx_LaMEM * FSLib->ny_LaMEM;

	PetscCall(PetscMalloc((size_t)(tnodes_ori) * sizeof(PetscScalar), &topo_fs));
	PetscCall(PetscMalloc((size_t)(tnodes_ori) * sizeof(PetscScalar), &topo_alloc));
	PetscCall(PetscMalloc((size_t)(tnodes_ori) * sizeof(PetscScalar), &vx_alloc));
	PetscCall(PetscMalloc((size_t)(tnodes_ori) * sizeof(PetscScalar), &vy_alloc));
	PetscCall(PetscMalloc((size_t)(tnodes_ori) * sizeof(PetscScalar), &vz_alloc));

	PetscCall(GatherVariableFromLaMEM(FSLib, topo_alloc, vx_alloc, vy_alloc, vz_alloc, step_fs));

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		PetscScalar dt_max = FSLib->Max_dt; // Maximum step length, if dt_LaMEM is larger than this, use this
		PetscScalar dt_n = 0; //dt_residual
		PetscScalar quotient = dt_fs/dt_max;
		PetscInt nsteps = (PetscInt)(floor( dt_fs/dt_max ));
		PetscScalar *topo_pass = nullptr, *vx_pass = nullptr, *vy_pass = nullptr, *vz_pass = nullptr;

		// store the phase that is being sedimented
		surf->phase = FSLib->sedPhases;

		tnodes      = FSLib->nodes_solve;

		PetscCall(PetscMalloc((size_t)(tnodes) * sizeof(PetscScalar), &topo_pass));
		PetscCall(PetscMalloc((size_t)(tnodes) * sizeof(PetscScalar), &vx_pass));
		PetscCall(PetscMalloc((size_t)(tnodes) * sizeof(PetscScalar), &vy_pass));
		PetscCall(PetscMalloc((size_t)(tnodes) * sizeof(PetscScalar), &vz_pass));

		// timestep
		if(nsteps < 1)
		{
			nsteps = 1;
			dt_max = dt_fs;
		}
		else if( (nsteps >= 1) && ((PetscScalar)(nsteps) != quotient) )
		{
			nsteps = nsteps + 1;
			dt_n   = dt_fs - (PetscScalar)(nsteps - 1) * dt_max;
		}

		FSLib->count_nsteps += nsteps;

		// run fastscape when using 3D geodynamic model
		if(FSLib->fs2D == 0)
		{
			// for non uniform grid
			if(FSLib->non_uniform_grid )
			{
				if(step_fs == 0)
				{
					PetscCall(InterpolationFor3DNonUniformGrid(FSLib, topo_alloc, 1));
				}
				PetscCall(InterpolationFor3DNonUniformGrid(FSLib, vx_alloc, 2));
				PetscCall(InterpolationFor3DNonUniformGrid(FSLib, vy_alloc, 3));
				PetscCall(InterpolationFor3DNonUniformGrid(FSLib, vz_alloc, 4));
			}

			// don't apply refinement
			if(FSLib->refine == 1)
			{
				if(FSLib->non_uniform_grid == 0)
				{
					PetscCall(VecGetArray(FSLib->gtopo_fs, &gtopo_rank0));

					for(i = 0; i < tnodes; i++)
					{
						if (0 == step_fs)   topo_pass[i] = topo_alloc[i] * scal->length * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
						else    topo_pass[i] = gtopo_rank0[i] * scal->length_fs;
						vx_pass[i]   = vx_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
						vy_pass[i]   = vy_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
						vz_pass[i]   = vz_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)

					}

					PetscCall(VecRestoreArray(FSLib->gtopo_fs, &gtopo_rank0));
				}
				else
				{
					PetscCall(VecGetArray(FSLib->gtopo_nug, &topo_solve));
					PetscCall(VecGetArray(FSLib->vx_nug, &vx_solve));
					PetscCall(VecGetArray(FSLib->vy_nug, &vy_solve));
					PetscCall(VecGetArray(FSLib->vz_nug, &vz_solve));

					for(i = 0; i < tnodes; i++)
					{
						topo_pass[i] = topo_solve[i] * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
						vx_pass[i]   = vx_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
						vy_pass[i]   = vy_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
						vz_pass[i]   = vz_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
					}

					PetscCall(VecRestoreArray(FSLib->gtopo_nug, &topo_solve));
					PetscCall(VecRestoreArray(FSLib->vx_nug, &vx_solve));
					PetscCall(VecRestoreArray(FSLib->vy_nug, &vy_solve));
					PetscCall(VecRestoreArray(FSLib->vz_nug, &vz_solve));
				}
			}
			// apply refinement
			else
			{
				printf("Refined times                    : %" PetscInt_FMT "\n", FSLib->refine);
				printf("Refined grid cells [nx, ny]      : [%" PetscInt_FMT ", %" PetscInt_FMT "] \n", FSLib->nx_refine, FSLib->ny_refine);

				PetscCall(VecGetArray(FSLib->gtopo_refine, &topo_solve_refine));
				PetscCall(VecGetArray(FSLib->vx_refine, &vx_solve_refine));
				PetscCall(VecGetArray(FSLib->vy_refine, &vy_solve_refine));
				PetscCall(VecGetArray(FSLib->vz_refine, &vz_solve_refine));

				if(FSLib->non_uniform_grid == 0)
				{
					// load velocity and topography field
					if(step_fs == 0)
					{
						PetscCall(BilinearInterpolate(FSLib, topo_alloc, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{
							topo_pass[i] = topo_solve_refine[i] * scal->length_fs;

						}
					}
					PetscCall(BilinearInterpolate(FSLib, vx_alloc, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));
					PetscCall(BilinearInterpolate(FSLib, vy_alloc, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));
					PetscCall(BilinearInterpolate(FSLib, vz_alloc, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));

				}
				else
				{
					PetscCall(VecGetArray(FSLib->gtopo_nug, &topo_solve));
					PetscCall(VecGetArray(FSLib->vz_nug, &vz_solve));
					PetscCall(VecGetArray(FSLib->vx_nug, &vx_solve));
					PetscCall(VecGetArray(FSLib->vy_nug, &vy_solve));

					// load velocity and topography field
					if(step_fs == 0)
					{
						PetscCall(BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{

							topo_pass[i] = topo_solve_refine[i] * scal->length_fs;

						}
					}

					PetscCall(BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));
					PetscCall(BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));
					PetscCall(BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine));

					PetscCall(VecRestoreArray(FSLib->gtopo_nug, &topo_solve));
					PetscCall(VecRestoreArray(FSLib->vz_nug, &vz_solve));
					PetscCall(VecRestoreArray(FSLib->vx_nug, &vx_solve));
					PetscCall(VecRestoreArray(FSLib->vy_nug, &vy_solve));
				}

				PetscCall(VecRestoreArray(FSLib->gtopo_refine, &topo_solve_refine));
				PetscCall(VecRestoreArray(FSLib->vx_refine, &vx_solve_refine));
				PetscCall(VecRestoreArray(FSLib->vy_refine, &vy_solve_refine));
				PetscCall(VecRestoreArray(FSLib->vz_refine, &vz_solve_refine));
			}
		}
		// run fastscape when using 2D geodynamic model
		else
		{
			PetscCall(VecGetArray(FSLib->gtopo_extend, &topo_solve));
			PetscCall(VecGetArray(FSLib->vx_extend, &vx_solve));
			PetscCall(VecGetArray(FSLib->vy_extend, &vy_solve));
			PetscCall(VecGetArray(FSLib->vz_extend, &vz_solve));

			if(FSLib->non_uniform_grid == 0)
			{
				// don't apply refinement
				if(FSLib->refine == 1)
				{
					// load velocity and topography field
					// step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
					// topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
					if(step_fs == 0)
					{
						PetscCall(Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1, PETSC_FALSE));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{
							topo_pass[i] = topo_solve[i] * scal->length_fs;
						}
					}

					// velocity
					// vx, vy
					if(FSLib->extendedX == 1)
					{
						PetscCall(Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2, PETSC_FALSE));
						std::fill(vx_pass, vx_pass + tnodes, 0);
					}
					else
					{
						PetscCall(Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2, PETSC_FALSE));
						std::fill(vy_pass, vy_pass + tnodes, 0);
					}
					// vz
					PetscCall(Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2, PETSC_FALSE));

				}
				// apply refinement
				else
				{
					// load velocity and topography field
					PetscCall(VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine));
					PetscCall(VecGetArray(FSLib->vx_et_refine, &vx_solve_refine));
					PetscCall(VecGetArray(FSLib->vy_et_refine, &vy_solve_refine));
					PetscCall(VecGetArray(FSLib->vz_et_refine, &vz_solve_refine));

					if(step_fs == 0)
					{
						PetscCall(Extended2D(FSLib, topo_alloc, topo_solve, NULL, NULL, 0, PETSC_TRUE));
						PetscCall(BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{
							topo_pass[i] = topo_solve_refine[i] * scal->length_fs;
						}
					}

					if(FSLib->extendedX == 1)
					{
						PetscCall(Extended2D(FSLib, vy_alloc, vy_solve, NULL, NULL, 0, PETSC_TRUE));
						std::fill(vx_solve, vx_solve + tnodes, 0);
					}
					else
					{
						PetscCall(Extended2D(FSLib, vx_alloc, vx_solve, NULL, NULL, 0, PETSC_TRUE));
						std::fill(vy_solve, vy_solve + tnodes, 0);
					}

					PetscCall(Extended2D(FSLib, vz_alloc, vz_solve, NULL, NULL, 0, PETSC_TRUE));

					PetscCall(BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					PetscCall(BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					PetscCall(BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));

					PetscCall(VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine));

				}
			}
			else
			{
				// don't apply refinement
				if(FSLib->refine == 1)
				{
					// load velocity and topography field
					// step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
					// topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
					if(step_fs == 0)
					{
						PetscCall(Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1, PETSC_FALSE));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{
							topo_pass[i] = topo_solve[i] * scal->length_fs;
						}
					}

					// velocity
					// vx, vy
					if(FSLib->extendedX == 1)
					{
						PetscCall(Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2, PETSC_FALSE));
						std::fill(vx_pass, vx_pass + tnodes, 0);
					}
					else
					{
						PetscCall(Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2, PETSC_FALSE));
						std::fill(vy_pass, vy_pass + tnodes, 0);
					}
					// vz
					PetscCall(Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2, PETSC_FALSE));
				}
				// apply refinement
				else
				{
					// load velocity and topography field
					PetscCall(VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine));
					PetscCall(VecGetArray(FSLib->vx_et_refine, &vx_solve_refine));
					PetscCall(VecGetArray(FSLib->vy_et_refine, &vy_solve_refine));
					PetscCall(VecGetArray(FSLib->vz_et_refine, &vz_solve_refine));

					if(step_fs == 0)
					{
						PetscCall(Extended2D(FSLib, topo_alloc, topo_solve, NULL, NULL, 0, PETSC_TRUE));
						PetscCall(BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					}
					else
					{
						for(i = 0; i < tnodes; i++)
						{
							topo_pass[i] = topo_solve_refine[i] * scal->length_fs;

						}
					}

					if(FSLib->extendedX == 1)
					{
						PetscCall(Extended2D(FSLib, vy_alloc, vy_solve, NULL, NULL, 0, PETSC_TRUE));
						std::fill(vx_solve, vx_solve + tnodes, 0);
					}
					else
					{
						PetscCall(Extended2D(FSLib, vx_alloc, vx_solve, NULL, NULL, 0, PETSC_TRUE));
						std::fill(vy_solve, vy_solve + tnodes, 0);
					}

					PetscCall(Extended2D(FSLib, vz_alloc, vz_solve, NULL, NULL, 0, PETSC_TRUE));

					PetscCall(BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					PetscCall(BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));
					PetscCall(BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes));

					PetscCall(VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine));
					PetscCall(VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine));
				}
			}

			PetscCall(VecRestoreArray(FSLib->gtopo_extend, &topo_solve));
			PetscCall(VecRestoreArray(FSLib->vx_extend, &vx_solve));
			PetscCall(VecRestoreArray(FSLib->vy_extend, &vy_solve));
			PetscCall(VecRestoreArray(FSLib->vz_extend, &vz_solve));
		}
		// run FastScape
		PetscCall( FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, FSLib->count_nsteps, (PetscScalar)step_fs, vx_pass, vy_pass, vz_pass, topo_pass));

		PetscCall(VecGetArray(FSLib->gtopo_fs, &gtopo_rank0));

		// pass data to original LaMEM grid
		if(FSLib->fs2D == 0)
		{
			PetscCall(PassValue3D(FSLib, topo_pass, gtopo_rank0));
		}
		else
		{
			PetscCall(PassValue2D(FSLib, topo_pass, gtopo_rank0));
		}

		PetscCall(PetscArraycpy(topo_fs, gtopo_rank0, tnodes_ori));
		PetscCall(VecRestoreArray(FSLib->gtopo_fs, &gtopo_rank0));

		// save data in a new grid used by fastscape
		if(step_fs == 0 || ((step_fs + 1) % FSLib->surf_out_nstep) == 0 )
		{
			PetscCall(FastScapeSave(FSLib, step_fs, time_fs));
		}

		PetscCall(PetscFree(topo_pass));
		PetscCall(PetscFree(vx_pass));
		PetscCall(PetscFree(vy_pass));
		PetscCall(PetscFree(vz_pass));
	}

	// Broadcast
	if(ISParallel(PETSC_COMM_WORLD))
	{
		PetscCall(MPI_Bcast(topo_fs, (PetscMPIInt)tnodes_ori, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD));
	}

	L    = (PetscInt)fs->dsz.rank;
	PetscCall(DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->gtopo, &topo));

	// Save topography in different rank
	START_PLANE_LOOP
	{
		topo[L][j][i] = topo_fs[j * FSLib->nx_LaMEM + i] / scal->length; // (km)
	}
	END_PLANE_LOOP

	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo));

	// compute ghosted version of the topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	// compute & store average topography
	PetscCall(FreeSurfGetAvgTopo(surf));

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscCall(PetscFree(topo_alloc));
	PetscCall(PetscFree(vx_alloc));
	PetscCall(PetscFree(vy_alloc));
	PetscCall(PetscFree(vz_alloc));
	PetscCall(PetscFree(topo_fs));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue2D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
	// initialize
	FSGrid  *fsX = &FSLib->fsX;
	FSGrid  *fsY = &FSLib->fsY;
	Scaling *scal = FSLib->scal;
	PetscScalar *topo_extend, *topo_aver, *topo_aver_ori, *topo_et_refine;
	PetscInt mainNodes, secNodes, laMEMNodes, i, j, ind, ind2, idx, sec, n, nn, countY = 0, countX = 0;
	PetscScalar sum, coord, begin, delta, coord0, coord1;
	// scaling function
	auto convertValue = [&](PetscScalar value) -> PetscScalar
	{
		PetscScalar converted = value / scal->length_fs;
		if (converted > FSLib->rangeZ_end) return FSLib->rangeZ_end;
		if (converted < FSLib->rangeZ_begin) return FSLib->rangeZ_begin;
		return converted;
	};

	// 1D linear interplotation
	auto linearInterp = [](PetscScalar coord, PetscScalar coord0, PetscScalar coord1,
	                       PetscScalar val0, PetscScalar val1) -> PetscScalar
	{
		PetscScalar weight = (coord - coord0) / (coord1 - coord0);
		return (1 - weight) * val0 + weight * val1;
	};

	PetscCall(VecGetArray(FSLib->gtopo_extend, &topo_extend));

	mainNodes = FSLib->extendedX ? FSLib->extendedYNodes : FSLib->extendedXNodes;
	secNodes = FSLib->extendedX ? FSLib->extendedXNodes : FSLib->extendedYNodes;
	laMEMNodes = FSLib->extendedX ? FSLib->ny_LaMEM : FSLib->nx_LaMEM;

	PetscCall(PetscMalloc1((size_t)(mainNodes), &topo_aver));
	PetscCall(PetscMalloc1((size_t)(laMEMNodes), &topo_aver_ori));

	// refine
	if (FSLib->refine > 1)
	{
		PetscCall(VecGetArray(FSLib->gtopo_et_refine, &topo_et_refine));

		for (j = 0; j < FSLib->etRefineYNodes; j++)
		{
			for (i = 0; i < FSLib->etRefineXNodes; i++)
			{
				ind = j * FSLib->etRefineXNodes + i;
				topo_et_refine[ind] = convertValue(topo_pass_f[ind]);
			}
		}

		countY = 0;
		for (j = 0; j < FSLib->etRefineYNodes; j += FSLib->refine)
		{
			countX = 0;
			for (i = 0; i < FSLib->etRefineXNodes; i += FSLib->refine)
			{
				ind = j * FSLib->etRefineXNodes + i;
				ind2 = countY * FSLib->extendedXNodes + countX;
				topo_extend[ind2] = topo_et_refine[ind];
				countX++;
			}
			countY++;
		}
		PetscCall(VecRestoreArray(FSLib->gtopo_et_refine, &topo_et_refine));
	}
	else
	{
		// no refine
		for (j = 0; j < secNodes; j++)
		{
			for (i = 0; i < mainNodes; i++)
			{
				ind = FSLib->extendedX ? (j * mainNodes + i) : (i * secNodes + j);
				topo_extend[ind] = convertValue(topo_pass_f[ind]);
			}
		}
	}

	for (idx = 0; idx < mainNodes; idx++)
	{
		sum = 0.0;
		for (sec = 0; sec < secNodes; sec++)
		{
			ind = FSLib->extendedX ? (idx * secNodes + sec) : (sec * mainNodes + idx);
			sum += topo_extend[ind];
		}

		topo_aver[idx] = sum / (PetscScalar)(secNodes);
	}

	// non-uniform grid
	if (FSLib->non_uniform_grid)
	{
		for (idx = 0; idx < laMEMNodes; idx++)
		{
			coord = FSLib->extendedX ? FSLib->ncoor_ori_y[idx] : FSLib->ncoor_ori_x[idx];
			begin = FSLib->extendedX ? FSLib->ncoor_ori_y[0] : FSLib->ncoor_ori_x[0];
			delta = FSLib->extendedX ? fsY->dx : fsX->dx;

			n = (PetscInt)(PetscFloorReal((coord - begin) / delta));
			nn = n + 1;

			if (nn >= mainNodes) nn = n;
			if (n == nn)  topo_aver_ori[idx] = topo_aver[n];
			else
			{
				coord0 = begin + (PetscScalar)n * delta;
				coord1 = begin + (PetscScalar)nn * delta;
				topo_aver_ori[idx] = linearInterp(coord, coord0, coord1, topo_aver[n], topo_aver[nn]);
			}
		}
	}
	else    PetscArraycpy(topo_aver_ori, topo_aver, laMEMNodes); // uniform grid

	for (j = 0; j < FSLib->ny_LaMEM; j++)
	{
		for (i = 0; i < FSLib->nx_LaMEM; i++)
		{
			ind = j * FSLib->nx_LaMEM + i;
			topo_fs[ind] = FSLib->extendedX ? topo_aver_ori[j] : topo_aver_ori[i];
		}
	}

	PetscCall(VecRestoreArray(FSLib->gtopo_extend, &topo_extend));
	PetscCall(PetscFree(topo_aver));
	PetscCall(PetscFree(topo_aver_ori));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue3D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
	FSGrid  *fsX;
	FSGrid  *fsY;
	Scaling *scal;
	PetscInt m, n, mm, nn, ind_a, ind_b, ind_c, ind_d, refine, i, j, ind, ind2, countX = 0, countY = 0;
	PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x_begin = 0.0, y_begin = 0.0;
	PetscScalar *topo_refine = nullptr, *topo_nug = nullptr;

	fsX     = &FSLib->fsX;
	fsY     = &FSLib->fsY;
	scal    = FSLib->scal;
	refine  = FSLib->refine;

	// for scaling
	auto convertValue = [&](PetscScalar value) -> PetscScalar
	{
		PetscScalar converted = value / scal->length_fs;
		if (FSLib->rangeZ_end < converted) return FSLib->rangeZ_end;
		if (FSLib->rangeZ_begin > converted) return FSLib->rangeZ_begin;
		return converted;
	};

	if(FSLib->non_uniform_grid == 1)
	{
		x_begin = FSLib->ncoor_ori_x[0];
		y_begin = FSLib->ncoor_ori_y[0];
	}

	if(FSLib->non_uniform_grid == 0)
	{
		if(refine == 1)
		{
			for(j = 0; j < FSLib->ny_LaMEM; j++)
			{
				for(i = 0; i < FSLib->nx_LaMEM; i++)
				{
					ind = j * FSLib->nx_LaMEM + i;
					topo_fs[ind] = convertValue(topo_pass_f[ind]);
				}
			}
		}
		else
		{
			PetscCall(VecGetArray(FSLib->gtopo_refine, &topo_refine));

			for(j = 0; j < FSLib->ny_refine; j++)
			{
				for (i = 0; i < FSLib->nx_refine; i++)
				{
					ind = j * FSLib->nx_refine + i;
					topo_refine[ind] = convertValue(topo_pass_f[ind]);
				}
			}

			for(j = 0; j < FSLib->ny_refine; j += refine)
			{
				for(i = 0; i < FSLib->nx_refine; i += refine)
				{
					ind = j * FSLib->nx_refine + i;
					ind2 = countY * FSLib->nx_LaMEM + countX;
					topo_fs[ind2] = topo_refine[ind];

					if (i == FSLib->nx_refine - 1) {   countX = 0;     countY++;}
					else    countX++;
				}
			}
			PetscCall(VecRestoreArray(FSLib->gtopo_refine,  &topo_refine));
		}
	}
	else
	{
		PetscCall(VecGetArray(FSLib->gtopo_nug, &topo_nug));

		if(refine == 1)
		{
			// save value
			for(j = 0; j < fsY->nodes; j++)
			{
				for(i = 0; i < fsX->nodes; i++)
				{
					ind = j * fsX->nodes + i;
					topo_nug[ind] =  convertValue(topo_pass_f[ind]);
				}
			}

			// interploate
			for(j = 0; j < FSLib->ny_LaMEM; j++)
			{
				for(i = 0; i < FSLib->nx_LaMEM; i++)
				{
					x_coor = FSLib->ncoor_ori_x[i];
					y_coor = FSLib->ncoor_ori_y[j];

					dx     = fsX->dx;
					dy     = fsY->dx;

					// get nearest four index
					// x-direction
					m  = (PetscInt)(floor( (x_coor - x_begin) / dx ));
					mm = (m + 1 < fsX->nodes) ? m + 1 : m;

					// y-direction
					n  = (PetscInt)(floor( (y_coor - y_begin) / dy ));
					nn = (n + 1 < fsY->nodes) ? n + 1 : n;

					// interpolate
					ind   = j  * FSLib->nx_LaMEM + i;
					ind_a = n  * fsX->nodes   + m;
					ind_b = n  * fsX->nodes   + mm;
					ind_c = nn * fsX->nodes   + m;
					ind_d = nn * fsX->nodes   + mm;
					// bilinear interpolation
					wtx   = x_coor - fsX->ncoor[m];
					wty   = y_coor - fsY->ncoor[n];

					topo_fs[ind] = ReturnBiInterFunction(
					                   topo_nug[ind_a], topo_nug[ind_b], topo_nug[ind_c], topo_nug[ind_d],
					                   dx, dy, wtx, wty, m, mm, n, nn);
				}
			}
		}
		else
		{
			PetscCall(VecGetArray(FSLib->gtopo_refine, &topo_refine));

			// save value
			for(j = 0; j < FSLib->ny_refine; j++)
			{
				for (i = 0; i < FSLib->nx_refine; i++)
				{
					ind = j * FSLib->nx_refine + i;
					topo_refine[ind] = convertValue(topo_pass_f[ind]);
				}
			}

			for(j = 0; j < FSLib->ny_refine; j += refine)
			{
				for(i = 0; i < FSLib->nx_refine; i += refine)
				{
					ind = j * FSLib->nx_refine + i;
					ind2 = countY * fsX->nodes + countX;
					topo_nug[ind2] = topo_refine[ind];

					if (i == FSLib->nx_refine - 1) { countX = 0; countY++;}
					else    countX++;
				}
			}

			// interploate
			for(j = 0; j < FSLib->ny_LaMEM; j++)
			{
				for(i = 0; i < FSLib->nx_LaMEM; i++)
				{
					x_coor = FSLib->ncoor_ori_x[i];
					y_coor = FSLib->ncoor_ori_y[j];

					dx     = fsX->dx;
					dy     = fsY->dx;

					// get nearest four index
					// x-direction
					m  = (PetscInt)(floor( (x_coor - x_begin) / dx ));
					mm = (m + 1 < fsX->nodes) ? m + 1 : m;

					// y-direction
					n  = (PetscInt)(floor( (y_coor - y_begin) / dy ));
					nn = (n + 1 < fsY->nodes) ? n + 1 : n;

					// interpolate
					ind   = j  * FSLib->nx_LaMEM + i;
					ind_a = n  * fsX->nodes   + m;
					ind_b = n  * fsX->nodes   + mm;
					ind_c = nn * fsX->nodes   + m;
					ind_d = nn * fsX->nodes   + mm;
					// bilinear interpolation
					wtx   = x_coor - fsX->ncoor[m];
					wty   = y_coor - fsY->ncoor[n];

					topo_fs[ind] = ReturnBiInterFunction(
					                   topo_nug[ind_a], topo_nug[ind_b], topo_nug[ind_c], topo_nug[ind_d],
					                   dx, dy, wtx, wty, m, mm, n, nn);
				}
			}
			PetscCall(VecRestoreArray(FSLib->gtopo_refine,  &topo_refine));
		}
		PetscCall(VecRestoreArray(FSLib->gtopo_nug, &topo_nug));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteVTSFS(FastScapeLib *FSLib, const char *dirName, PetscScalar *topo, PetscInt mode)
{
	FILE      *fp;
	Scaling   *scal;
	FreeSurf  *surf;
	size_t    offset = 0;
	PVSurf    *pvsurf;
	PetscFunctionBeginUser;

	PetscInt i, nx, ny, numFields;
	PetscScalar *silt_fraction = nullptr, *basement = nullptr, *total_erosion = nullptr;
	PetscScalar *drainage_area = nullptr, *erosion_rate = nullptr, *slope = nullptr;
	PetscScalar *curvature = nullptr, *chi = nullptr, *catchment = nullptr, *lake_depth = nullptr;

	// access context
	pvsurf  = FSLib->pvsurf;
	surf    = pvsurf->surf;
	scal    = surf->jr->scal;

	// only processor 0 run the code
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	fp = NULL;

	// initialize
	nx = FSLib->nx_solve;
	ny = FSLib->ny_solve;

	char fname[PETSC_MAX_PATH_LEN];

	const size_t gridSize = (size_t)(nx * ny);

	// create file name
	PetscCall(PetscSNPrintf(fname, sizeof(fname), "%s/%s_p0.vts", dirName, FSLib->outfile_fs));

	// 打开文件
	fp = fopen(fname, "wb");
	if (fp == NULL)     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open file %s", fname);

	WriteXMLHeader(fp, "StructuredGrid");

	// grid info
	fprintf(fp, "\t<StructuredGrid WholeExtent=\"%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT " 1 1\">\n",
	        (PetscInt)1, nx, (PetscInt)1, ny);
	fprintf(fp, "\t\t<Piece Extent=\"%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT " 1 1\">\n",
	        (PetscInt)1, nx, (PetscInt)1, ny);
	fprintf(fp, "\t\t\t<CellData/>\n");

	// offset & nodes
	fprintf(fp, "\t\t<Points>\n");
	fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" "
	            "format=\"appended\" offset=\"%zu\"/>\n", offset);
	offset += sizeof(uint64_t) + sizeof(float) * gridSize * 3;
	fprintf(fp, "\t\t</Points>\n");

	// point data
	fprintf(fp, "\t\t<PointData>\n");

	// output info
	typedef struct
	{
		PetscInt *flag;
		const char *name;
		const char *unit;
	} OutputField;

	OutputField fields[] =
	{
		{&FSLib->out_topofs, "topoFs", scal->lbl_length},
		{&FSLib->out_silt_fraction, "silt_fraction", scal->lbl_unit},
		{&FSLib->out_basement, "basement", scal->lbl_length},
		{&FSLib->out_total_erosion, "total_erosion", scal->lbl_length},
		{&FSLib->out_drainage_area, "drainage_area", scal->lbl_area_fs},
		{&FSLib->out_erosion_rate, "erosion_rate", scal->lbl_rate},
		{&FSLib->out_slope, "slope", scal->lbl_degree},
		{&FSLib->out_curvature, "curvature", scal->lbl_unit},
		{&FSLib->out_chi, "chi", scal->lbl_length},
		{&FSLib->out_catchment, "catchment", scal->lbl_area_fs},
		{&FSLib->out_lake_depth, "lake_depth", scal->lbl_length}
	};

	numFields = sizeof(fields) / sizeof(fields[0]);

	for (i = 0; i < numFields; i++)
	{
		if (*(fields[i].flag))
		{
			fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"%s %s\" "
			            "NumberOfComponents=\"1\" format=\"appended\" offset=\"%zu\"/>\n",
			        fields[i].name, fields[i].unit, offset);
			offset += sizeof(uint64_t) + sizeof(float) * gridSize;
		}
	}

	fprintf(fp, "\t\t</PointData>\n");
	fprintf(fp, "\t\t</Piece>\n");
	fprintf(fp, "\t</StructuredGrid>\n");

	fprintf(fp, "\t<AppendedData encoding=\"raw\">\n_");
	// write point coordinates
	// allocate output buffer
	PetscCall(PetscMalloc((size_t)(_max_num_comp_surf_ * nx * ny)*sizeof(float), &FSLib->buff_fs));

	PetscCall(PVSurfWriteCoordFS (FSLib, fp, topo, mode));

	// topography
	if(FSLib->out_topofs)
	{
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, topo, 1));
	}
	// silt fraction
	if(FSLib->out_silt_fraction)
	{
		PetscCall(VecGetArray(FSLib->silt_fraction, &silt_fraction));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, silt_fraction, 2));
		PetscCall(VecRestoreArray(FSLib->silt_fraction, &silt_fraction));
	}
	// basement
	if(FSLib->out_basement)
	{
		PetscCall(VecGetArray(FSLib->basement, &basement));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, basement, 3));
		PetscCall(VecRestoreArray(FSLib->basement, &basement));
	}
	// total_erosion
	if(FSLib->out_total_erosion)
	{
		PetscCall(VecGetArray(FSLib->total_erosion, &total_erosion));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, total_erosion, 4));
		PetscCall(VecRestoreArray(FSLib->total_erosion, &total_erosion));
	}
	// drainage_area
	if(FSLib->out_drainage_area)
	{
		PetscCall(VecGetArray(FSLib->drainage_area, &drainage_area));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, drainage_area, 5));
		PetscCall(VecRestoreArray(FSLib->drainage_area, &drainage_area));
	}
	// erosion_rate
	if(FSLib->out_erosion_rate)
	{
		PetscCall(VecGetArray(FSLib->erosion_rate, &erosion_rate));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, erosion_rate, 6));
		PetscCall(VecRestoreArray(FSLib->erosion_rate, &erosion_rate));
	}
	// slope
	if(FSLib->out_slope)
	{
		PetscCall(VecGetArray(FSLib->slope, &slope));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, slope, 7));
		PetscCall(VecRestoreArray(FSLib->slope, &slope));
	}
	// curvature
	if(FSLib->out_curvature)
	{
		PetscCall(VecGetArray(FSLib->curvature, &curvature));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, curvature, 8));
		PetscCall(VecRestoreArray(FSLib->curvature, &curvature));
	}
	// chi
	if(FSLib->out_chi)
	{
		PetscCall(VecGetArray(FSLib->chi, &chi));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, chi, 9));
		PetscCall(VecRestoreArray(FSLib->chi, &chi));
	}
	// catchment
	if(FSLib->out_catchment)
	{
		PetscCall(VecGetArray(FSLib->catchment, &catchment));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, catchment, 10));
		PetscCall(VecRestoreArray(FSLib->catchment, &catchment));
	}
	// lake_depth
	if(FSLib->out_lake_depth)
	{
		PetscCall(VecGetArray(FSLib->lake_depth, &lake_depth));
		PetscCall(PVSurfWriteInfFS  (FSLib, fp, lake_depth, 11));
		PetscCall(VecRestoreArray(FSLib->lake_depth, &lake_depth));
	}
	// close appended data section and file
	fprintf(fp, "\n\t</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");

	// close file
	fclose(fp);

	PetscCall(PetscFree(FSLib->buff_fs));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteCoordFS(FastScapeLib *FSLib, FILE *fp, PetscScalar *topo, PetscInt mode)
{
	float       *buff;
	PetscInt    i, j, ind, cn, nx, ny;
	FSGrid  *fsX;
	FSGrid  *fsY;
	PetscFunctionBeginUser;

	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	fsX = &FSLib->fsX;
	fsY = &FSLib->fsY;

	// initialize
	buff   = FSLib->buff_fs;

	nx = FSLib->nx_solve;
	ny = FSLib->ny_solve;

	cn     = 0;

	for(j = 0; j < ny; j++)
	{
		for(i = 0; i < nx; i++)
		{
			ind = j * nx + i;

			if(mode == 0)
			{
				// store node coordinates
				buff[cn++] = (float)(fsX->ncoor[i]);
				buff[cn++] = (float)(fsY->ncoor[j]);
				buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m
			}
			else if(mode == 2)
			{
				// store node coordinates
				buff[cn++] = (float)(fsX->ncoor_extend[i]);
				buff[cn++] = (float)(fsY->ncoor_extend[j]);
				buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m
			}
			else if(mode == 1 || mode == 3)
			{
				buff[cn++] = (float)(fsX->ncoor_refine[i]);
				buff[cn++] = (float)(fsY->ncoor_refine[j]);
				buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m
			}
		}
	}

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteInfFS(FastScapeLib *FSLib, FILE *fp, PetscScalar *Inf, PetscInt InfMode)
{
	PetscFunctionBeginUser;
	float       *buff;
	PetscInt    ind, cn, nx, ny;
	Scaling *scal;

	scal = FSLib->scal;

	// initialize
	buff   = FSLib->buff_fs;

	nx = FSLib->nx_solve;
	ny = FSLib->ny_solve;

	cn   = 0;

	for (PetscInt j = 0; j < ny; j++)
	{
		for (PetscInt i = 0; i < nx; i++)
		{
			ind = j * nx + i;

			switch (InfMode)
			{
			case 1:  // topography
			case 2:  // silt fraction
			case 7:  // slope
			case 8:  // curvature
				buff[cn++] = (float)(Inf[ind]);
				break;
			case 3:  // basement
			case 4:  // total_erosion
			case 6:  // erosion
			case 9:  // chi
			case 11: // lake_depth
				buff[cn++] = (float)(Inf[ind] / scal->length_fs);
				break;
			case 5:  // drainage_area
			case 10: // catchment
				buff[cn++] = (float)(Inf[ind] / scal->area_fs);
				break;
			}
		}
	}

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode UpdatePVDFileFS(
    const char *dirName, const char *outfile, const char *ext,
    long int *offset, PetscScalar ttime, PetscInt outpvd, PetscInt step)
{
	PetscFunctionBeginUser;

	FILE        *fp;
	char        *fname;

	// check whether pvd is requested
	if(!outpvd) PetscFunctionReturn(0);

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// open outfile.pvd file (write or update mode)
	asprintf(&fname, "%s.pvd", outfile);
	if(step == 1) fp = fopen(fname,"wb");
	else       fp = fopen(fname,"r+b");

	if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);

	free(fname);

	if(step == 1)
	{
		// write header
		WriteXMLHeader(fp, "Collection");

		// open time step collection
		fprintf(fp,"<Collection>\n");
	}
	else
	{
		// put the file pointer on the next entry
		PetscCall(fseek(fp, (*offset), SEEK_SET));
	}

	// add entry to .pvd file
	fprintf(fp,"\t<DataSet timestep=\"%1.6e\" file=\"%s/%s_%s\"/>\n",
	        ttime, dirName, outfile, ext);

	// store current position in the file
	(*offset) = ftell(fp);

	// close time step collection
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode SavePvtsFS(FastScapeLib *FSLib, PetscScalar ttime, PetscInt step, const char *dirName, PetscScalar *topo)
{
	PetscFunctionBeginUser;
	PetscInt mode = 0;

	// check activation
	if(!FSLib->outsurf_fs) PetscFunctionReturn(0);

	// update .pvd file if necessary
	PetscCall(UpdatePVDFileFS(dirName, FSLib->outfile_fs, "p0.vts", &FSLib->offset_fs, ttime, FSLib->outpvd_fs, step));

	// write sub-domain data .vts files
	if (!FSLib->fs2D)   mode = (FSLib->refine == 1) ? 0 : 1;
	else    mode = (FSLib->refine == 1) ? 2 : 3;

	PetscCall(PVSurfWriteVTSFS(FSLib, dirName, topo, mode));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeSave(FastScapeLib *FSLib, PetscInt step_fs, PetscScalar time_fs)
{
	PetscInt status;
	char *dirName = NULL;
	PetscScalar *saveArray = nullptr;
	Vec saveVec = nullptr;

	// create directory(encode current time & steo number)
	// update time stamp and counter
	step_fs++;
	PetscCall(PetscMalloc1(256, &dirName));
	PetscCall(PetscSNPrintf(dirName, 256, "Timestep_%1.8" PetscInt_FMT "_%1.8e", step_fs, time_fs));

	// create output directory
#ifdef _WIN32
	// call this on windows machines
	status = mkdir(dirName);
#else
	// standard access pattern drwxr-xr-x
	status = mkdir(dirName, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
	if(status && errno != EEXIST)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to create directory %s", dirName);
	}

	// only saved in processor 0
	if(ISRankZero(PETSC_COMM_WORLD) && FSLib->outsurf_fs)
	{
		if (!FSLib->fs2D)
		{
			// 3D
			saveVec = (FSLib->refine == 1) ?
			          (FSLib->non_uniform_grid ? FSLib->gtopo_nug : FSLib->gtopo_fs) :
			          FSLib->gtopo_refine;
		}
		else
		{
			// 2D
			saveVec = (FSLib->refine == 1) ?
			          FSLib->gtopo_extend :
			          FSLib->gtopo_et_refine;
		}

		// save data
		PetscCall(VecGetArray(saveVec, &saveArray));
		PetscCall(SavePvtsFS(FSLib, time_fs, step_fs, dirName, saveArray));
		PetscCall(VecRestoreArray(saveVec, &saveArray));
	}

	// free resource
	PetscCall(PetscFree(dirName));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeInitialize(FastScapeLib *FSLib, PetscScalar *topo_pass, PetscInt restart)
{
	PetscScalar rangeX, rangeY;
	PetscScalar *topo_random = nullptr;
	PetscInt ind;
	int ibc_int;
	char *endptr;

	// random noise (uniform distribution)
	mt19937 generator;
	uniform_int_distribution<int> distribution(1, 10000);

	fastscape_init_();
	fastscape_set_nx_ny_(&FSLib->nx_solve, &FSLib->ny_solve);

	// allocate memory
	fastscape_setup_();

	// set model dimensions
	rangeX = (0 == FSLib->fs2D) ? FSLib->rangeX : FSLib->extendedXRange;
	rangeY = (0 == FSLib->fs2D) ? FSLib->rangeY : FSLib->extendedYRange;
	fastscape_set_xl_yl_(&rangeX, &rangeY);

	// set topography boundary conditions
	ibc_int = (int)(strtol(FSLib->FS_BC, &endptr, 10));
	fastscape_set_bc_(&ibc_int);

	// random noise
	if (FSLib->random_noise && restart == 0)
	{
		PetscCall(PetscMalloc((size_t)(FSLib->nodes_solve) * sizeof(PetscScalar), &topo_random));
		for (ind = 0; ind < FSLib->nodes_solve; ind++)
		{
			topo_random[ind] = distribution(generator)/10000.0;
			topo_pass[ind]   += topo_random[ind];
		}
		PetscCall(PetscFree(topo_random));
	}

	// initialize topography
	fastscape_init_h_(topo_pass);

	if(restart == 1)
	{
		PetscScalar *silt_fraction = nullptr, *basement = nullptr;

		PetscCall(VecGetArray(FSLib->basement, &basement));
		fastscape_set_basement_(basement);
		PetscCall(VecRestoreArray(FSLib->basement, &basement));

		if(FSLib->setMarine == 1)
		{
			PetscCall(VecGetArray(FSLib->silt_fraction, &silt_fraction));
			fastscape_init_f_(silt_fraction);
			PetscCall(VecRestoreArray(FSLib->silt_fraction, &silt_fraction));
		}
	}

	printf("initialize fastscape process\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeFortranCppAdvc(FastScapeLib *FSLib, PetscScalar dt_max, PetscScalar dt_n, PetscInt nstep,
                                       PetscScalar step_fs, PetscScalar *vx_pass, PetscScalar *vy_pass, PetscScalar *vz_pass, PetscScalar *topo_pass)
{
	int istep;
	PetscInt ind, i, j;
	PetscScalar *kf = nullptr, *kd = nullptr;
	PetscScalar *silt_fraction = nullptr, *basement = nullptr, *total_erosion = nullptr;
	PetscScalar *drainage_area = nullptr, *erosion_rate = nullptr, *slope = nullptr;
	PetscScalar *curvature = nullptr, *chi = nullptr, *catchment = nullptr, *lake_depth = nullptr;


	std::vector<PetscScalar*> output_arrays;
	std::vector<Vec> output_vecs;
	std::vector<PetscInt> output_flags;
	std::vector<void (*)(PetscScalar*)> output_functions;

	PetscCall(PetscMalloc((size_t)(FSLib->nodes_solve) * sizeof(PetscScalar), &kf));
	PetscCall(PetscMalloc((size_t)(FSLib->nodes_solve) * sizeof(PetscScalar), &kd));

	// output
	auto manage_output = [&](PetscInt flag, Vec vec, PetscScalar** array, void (*copy_func)(PetscScalar*)) -> PetscErrorCode
	{
		if (flag)
	{
		PetscCall(VecGetArray(vec, array));
			output_flags.push_back(flag);
			output_arrays.push_back(*array);
			output_vecs.push_back(vec);
			output_functions.push_back(copy_func);
		}
		return 0;
	};

	PetscCall(manage_output(FSLib->out_silt_fraction, FSLib->silt_fraction, &silt_fraction, fastscape_copy_f_));
	PetscCall(manage_output(FSLib->out_basement, FSLib->basement, &basement, fastscape_copy_basement_));
	PetscCall(manage_output(FSLib->out_total_erosion, FSLib->total_erosion, &total_erosion, fastscape_copy_total_erosion_));
	PetscCall(manage_output(FSLib->out_drainage_area, FSLib->drainage_area, &drainage_area, fastscape_copy_drainage_area_));
	PetscCall(manage_output(FSLib->out_erosion_rate, FSLib->erosion_rate, &erosion_rate, fastscape_copy_erosion_rate_));
	PetscCall(manage_output(FSLib->out_slope, FSLib->slope, &slope, fastscape_copy_slope_));
	PetscCall(manage_output(FSLib->out_curvature, FSLib->curvature, &curvature, fastscape_copy_curvature_));
	PetscCall(manage_output(FSLib->out_chi, FSLib->chi, &chi, fastscape_copy_chi_));
	PetscCall(manage_output(FSLib->out_catchment, FSLib->catchment, &catchment, fastscape_copy_catchment_));
	PetscCall(manage_output(FSLib->out_lake_depth, FSLib->lake_depth, &lake_depth, fastscape_copy_lake_depth_));

	// initialize FastScape
	if(step_fs == 0)
	{
		FastScapeInitialize(FSLib, topo_pass, 0);
	}
	// restart
	else if(FSLib->restart == 1)
	{
		FastScapeInitialize(FSLib, topo_pass, 1);
		FSLib->restart = 0;
	}

	// kf & kd
	std::fill_n(kf, FSLib->nodes_solve, FSLib->kf);
	std::fill_n(kd, FSLib->nodes_solve, FSLib->kd);

	// velocity boundary conditions
	auto set_velocity_boundary = [&](PetscInt i, PetscInt j, PetscInt ind)
	{
		// bottom
		if (FSLib->FS_VELBC[0] == '1' && j == 0)
		{
			vx_pass[ind] = 0;
			vy_pass[ind] = 0;
			vz_pass[ind] = 0;
		}
		// right
		if (FSLib->FS_VELBC[1] == '1' && i == (FSLib->nx_solve - 1))
		{
			vx_pass[ind] = 0;
			vy_pass[ind] = 0;
			vz_pass[ind] = 0;
		}
		// top
		if (FSLib->FS_VELBC[2] == '1' && j == (FSLib->ny_solve - 1))
		{
			vx_pass[ind] = 0;
			vy_pass[ind] = 0;
			vz_pass[ind] = 0;
		}
		// left
		if (FSLib->FS_VELBC[3] == '1' && i == 0)
		{
			vx_pass[ind] = 0;
			vy_pass[ind] = 0;
			vz_pass[ind] = 0;
		}
	};

	for(j = 0; j < FSLib->ny_solve; j++)
	{
		for(i = 0; i < FSLib->nx_solve; i++)
		{
			ind = j * FSLib->nx_solve + i;
			set_velocity_boundary(i, j, ind);
		}
	}

	// set velocity
	fastscape_set_u_(vz_pass);
	fastscape_set_v_(vx_pass, vy_pass);

	// set time step
	fastscape_set_dt_(&dt_max);

	// reset topography
	if(step_fs > 0)
	{
		fastscape_set_h_(topo_pass);
	}

	// set erosional parameters
	fastscape_set_erosional_parameters_(kf, &FSLib->kfsed, &FSLib->m, &FSLib->n, kd, &FSLib->kdsed, &FSLib->g, &FSLib->gsed, &FSLib->p);

	// set marine transport parameters
	if (FSLib->setMarine)
	{
		fastscape_set_marine_parameters_(&FSLib->sealevel, &FSLib->poro_silt, &FSLib->poro_sand, &FSLib->zporo_silt,
		                                 &FSLib->zporo_sand, &FSLib->ratio, &FSLib->Lsolve, &FSLib->kds_silt, &FSLib->kds_sand);
	}

	// set number of time steps and initialize counter istep
	fastscape_get_step_(&istep);

	// loop on time stepping
	do
	{
		// time step in the last step
		if (istep == nstep - 1 && dt_n > 0)
		{
			fastscape_set_dt_(&dt_n);
		}

		// execute step
		fastscape_execute_step_();
		// get value of time step counter
		fastscape_get_step_(&istep);

		if(istep == nstep)
		{
			// output
			// topography
			fastscape_copy_h_(topo_pass);
			// others
			for (size_t i = 0; i < output_arrays.size(); i++)
			{
				if (output_flags[i]) output_functions[i](output_arrays[i]);
			}
		}
	}
	while (istep < nstep);

	// output timing
	// fastscape_debug_();

	// end FastScape run
	// fastscape_destroy_();

	// save data
	for (size_t i = 0; i < output_arrays.size(); i++)
	{
		if(output_flags[i])
		{
			PetscCall(VecRestoreArray(output_vecs[i], &output_arrays[i]));
		}
	}

	PetscCall(PetscFree(kf));
	PetscCall(PetscFree(kd));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeReadRestart(FastScapeLib *FSLib, FILE *fp)
{
	PetscFunctionBeginUser;

	FSLib->restart = 1;
	FSLib->count_nsteps = 0;

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		PetscCall(FastScapeCreateSurfaceGrid(FSLib, 2));

		if(FSLib->fs2D == 0 && FSLib->refine == 1)
		{
			// read topography vector
			PetscCall(VecReadRestart(FSLib->gtopo_fs, fp));
		}
		if(FSLib->non_uniform_grid == 1)
		{
			// read topography vector
			PetscCall(VecReadRestart(FSLib->gtopo_nug, fp));
		}
		if(FSLib->fs2D == 0 && FSLib->refine > 1)
		{
			// read topography vector
			PetscCall(VecReadRestart(FSLib->gtopo_refine, fp));
		}
		else if(FSLib->fs2D == 1 && FSLib->refine == 1)
		{
			// read topography vector
			PetscCall(VecReadRestart(FSLib->gtopo_extend, fp));
		}
		else if(FSLib->fs2D == 1 && FSLib->refine > 1)
		{
			// read topography vector
			PetscCall(VecReadRestart(FSLib->gtopo_et_refine, fp));
		}

		// FastScape part
		PetscCall(VecReadRestart(FSLib->basement, fp));

		if(FSLib->setMarine == 1)
		{
			PetscCall(VecReadRestart(FSLib->silt_fraction, fp));
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeWriteRestart(FastScapeLib *FSLib, FILE *fp)
{
	PetscFunctionBeginUser;

	if(ISRankZero(PETSC_COMM_WORLD))
	{

		if(FSLib->fs2D == 0 && FSLib->refine == 1)
		{
			// store topography vector
			PetscCall(VecWriteRestart(FSLib->gtopo_fs, fp));
		}
		if(FSLib->non_uniform_grid == 1)
		{
			// store topography vector
			PetscCall(VecWriteRestart(FSLib->gtopo_nug, fp));
		}
		if(FSLib->fs2D == 0 && FSLib->refine > 1)
		{
			// store topography vector
			PetscCall(VecWriteRestart(FSLib->gtopo_refine, fp));
		}
		else if( (FSLib->fs2D == 1 && FSLib->refine == 1) )
		{
			// store topography vector
			PetscCall(VecWriteRestart(FSLib->gtopo_extend, fp));
		}
		else if(FSLib->fs2D == 1 && FSLib->refine > 1)
		{
			// store topography vector
			PetscCall(VecWriteRestart(FSLib->gtopo_et_refine, fp));
		}

		// FastScape part
		PetscCall(VecWriteRestart(FSLib->basement, fp));

		if(FSLib->setMarine)
		{
			PetscCall(VecWriteRestart(FSLib->silt_fraction, fp));
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeDestroy(FastScapeLib *FSLib)
{
	PetscFunctionBeginUser;

	PetscCall(VecDestroy(&FSLib->vx_collect));
	PetscCall(VecDestroy(&FSLib->vy_collect));
	PetscCall(VecDestroy(&FSLib->vz_collect));

	// release fastscape
	if (ISRankZero(PETSC_COMM_WORLD))
	{
		fastscape_destroy_();

		PetscCall(VecDestroy(&FSLib->gtopo_fs));

		if(FSLib->non_uniform_grid)
		{
			PetscCall(VecDestroy(&FSLib->gtopo_nug));
			PetscCall(VecDestroy(&FSLib->vx_nug));
			PetscCall(VecDestroy(&FSLib->vy_nug));
			PetscCall(VecDestroy(&FSLib->vz_nug));
		}

		if(FSLib->fs2D)
		{
			// extended part
			PetscCall(VecDestroy(&FSLib->gtopo_extend));
			PetscCall(VecDestroy(&FSLib->vx_extend));
			PetscCall(VecDestroy(&FSLib->vy_extend));
			PetscCall(VecDestroy(&FSLib->vz_extend));

			if(FSLib->refine > 1)
			{
				// extended part after refinement
				PetscCall(VecDestroy(&FSLib->gtopo_et_refine));
				PetscCall(VecDestroy(&FSLib->vx_et_refine));
				PetscCall(VecDestroy(&FSLib->vy_et_refine));
				PetscCall(VecDestroy(&FSLib->vz_et_refine));
			}
		}
		else
		{
			if(FSLib->refine > 1)
			{
				// refined part
				PetscCall(VecDestroy(&FSLib->gtopo_refine));
				PetscCall(VecDestroy(&FSLib->vx_refine));
				PetscCall(VecDestroy(&FSLib->vy_refine));
				PetscCall(VecDestroy(&FSLib->vz_refine));
			}
		}

		// FastScape solution
		PetscCall(VecDestroy(&FSLib->silt_fraction));
		PetscCall(VecDestroy(&FSLib->basement));
		PetscCall(VecDestroy(&FSLib->total_erosion));
		PetscCall(VecDestroy(&FSLib->drainage_area));
		PetscCall(VecDestroy(&FSLib->erosion_rate));
		PetscCall(VecDestroy(&FSLib->slope));
		PetscCall(VecDestroy(&FSLib->curvature));
		PetscCall(VecDestroy(&FSLib->chi));
		PetscCall(VecDestroy(&FSLib->catchment));
		PetscCall(VecDestroy(&FSLib->lake_depth));

		PetscCall(PetscFree(FSLib->fsX.ncoor));
		PetscCall(PetscFree(FSLib->fsY.ncoor));

		if(FSLib->non_uniform_grid)
		{
			if(FSLib->extendedX)
			{
				PetscCall(PetscFree(FSLib->ncoor_ori_y));
			}
			else
			{
				PetscCall(PetscFree(FSLib->ncoor_ori_x));
			}
		}

		if(FSLib->fs2D)
		{
			PetscCall(PetscFree(FSLib->fsX.ncoor_extend));
			PetscCall(PetscFree(FSLib->fsY.ncoor_extend));
		}

		if(FSLib->refine > 1)
		{
			PetscCall(PetscFree(FSLib->fsX.ncoor_refine));
			PetscCall(PetscFree(FSLib->fsY.ncoor_refine));
		}
	}

	PetscFunctionReturn(0);
}
