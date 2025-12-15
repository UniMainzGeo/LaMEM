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
#include "scaling.h"
#include "JacRes.h"
#include "tssolve.h"
#include "advect.h"
#include "interpolate.h"
#include "fastscape.h"
#include "paraViewOutSurf.h"
#include "paraViewOutBin.h"

//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreate(FastScapeLib *FSLib, FB *fb)
{ 
    PetscErrorCode ierr;
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
    maxPhaseID = FSLib->surf->jr->dbm->numPhases-1;
 
    // load information from .dat file
    // setup block access mode                                                                                                                                      
    ierr = FBFindBlocks(fb, _REQUIRED_, "<FastScapeStart>", "<FastScapeEnd>");      CHKERRQ(ierr);

    if(fb->nblocks)
    {
        //-------------------------------
        // Grid information
        //-------------------------------
        // non uniform grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "non_uniform_grid", &FSLib->non_uniform_grid, 1,  1);                      CHKERRQ(ierr); // flag 
        // 2D grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "fs2D",             &FSLib->fs2D,             1,  1);                      CHKERRQ(ierr); // flag
        if(FSLib->fs2D )
        {
            ierr     = getScalarParam(fb, _REQUIRED_, "extendedRange",    &FSLib->extendedRange,    1,  1/scal->length_fs);      CHKERRQ(ierr); // km (LaMEM) -> m (FastScape)
            
            if (!FSLib->non_uniform_grid)
            {
                ierr = getIntParam   (fb, _REQUIRED_, "extendedNodes",    &FSLib->extendedNodes,    1,  10000000);               CHKERRQ(ierr); // non-dimensional
                if(2 >= FSLib->extendedNodes)   SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "extendedNodes must be ≥ 2");
            }
        }    
        // refined grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "fs_refine",        &FSLib->refine,           1,  100);                    CHKERRQ(ierr); // non-dimensional
       
        //===============================                                                                                                                                            
        // FastScape PARAMETER                                                                                                               
        //===============================
        // dt & boundary condition
        ierr = getScalarParam(fb, _REQUIRED_, "Max_dt",            &FSLib->Max_dt,           1,  1 / scal->time_fs);    CHKERRQ(ierr); // Myr (LaMEM) ->yr (FastScape)
      
        // bottom-right-top-left; 0 = reflective, 1 = fixed height boundary; When two reflective boundaris face each other they become cyclic
        ierr = getStringParam(fb, _REQUIRED_, "topo_boundary",     FSLib->FS_BC,                 "1111");               CHKERRQ(ierr); 
        
        if(strlen(FSLib->FS_BC) != 4)   SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "topo_boundary must be 4-character string (e.g., '1100')");

        // 1 -- boundary velocity == 0; 0 -- boundary velocity from LaMEM
        ierr = getStringParam(fb, _REQUIRED_, "vel_boundary",      FSLib->FS_VELBC,              "1111");               CHKERRQ(ierr); 
        
        if(strlen(FSLib->FS_VELBC) != 4)    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "vel_boundary must be 4-character string (e.g., '1100')");

        // random noise
        ierr = getIntParam   (fb, _REQUIRED_, "random_noise",      &FSLib->random_noise,     1,  1);                    CHKERRQ(ierr); 

        // sedimentation phase
        ierr = getIntParam   (fb, _REQUIRED_, "sed_phases",        &FSLib->sedPhases,        1,  maxPhaseID);           CHKERRQ(ierr); // non-dimensional

        surf->phase = FSLib->sedPhases;
        //-------------------------------
        // Erosion process
        //-------------------------------
        // kf, kd can be set as an array
        ierr = getScalarParam(fb, _REQUIRED_, "kf",                &FSLib->kf,               1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "kfsed",             &FSLib->kfsed,            1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "m",                 &FSLib->m,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "n",                 &FSLib->n,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "kd",                &FSLib->kd,               1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "kdsed",             &FSLib->kdsed,            1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "g",                 &FSLib->g,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "gsed",              &FSLib->gsed,             1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "p",                 &FSLib->p,                1,          1.0);         CHKERRQ(ierr); // non-dimensional

        //-------------------------------
        // Sedimentation process
        //-------------------------------    
        ierr = getIntParam   (fb, _OPTIONAL_, "setMarine",         &FSLib->setMarine,        1,         1);                     CHKERRQ(ierr); // flag
        if(FSLib->setMarine)
        {
            ierr = getScalarParam(fb, _REQUIRED_, "sealevel",      &FSLib->sealevel,         1,         1 / scal->length_fs);   CHKERRQ(ierr); // m
            ierr = getScalarParam(fb, _REQUIRED_, "poroSilt",      &FSLib->poro_silt,        1,         1.0);                   CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "poroSand",      &FSLib->poro_sand,        1,         1.0);                   CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "zporoSilt",     &FSLib->zporo_silt,       1,         1.0);                   CHKERRQ(ierr); 
            ierr = getScalarParam(fb, _REQUIRED_, "zporoSand",     &FSLib->zporo_sand,       1,         1.0);                   CHKERRQ(ierr); 
            ierr = getScalarParam(fb, _REQUIRED_, "ratio",         &FSLib->ratio,            1,         1.0);                   CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "Lsolve",        &FSLib->Lsolve,           1,         1.0);                   CHKERRQ(ierr); // m
            ierr = getScalarParam(fb, _REQUIRED_, "kdsSilt",       &FSLib->kds_silt,         1,         1.0);                   CHKERRQ(ierr); // m2/yr
            ierr = getScalarParam(fb, _REQUIRED_, "kdsSand",       &FSLib->kds_sand,         1,         1.0);                   CHKERRQ(ierr); // m2/yr
        }

        //-------------------------------
        // Output information
        //-------------------------------
        ierr = getIntParam   (fb, _OPTIONAL_, "surf_out_nstep",    &FSLib->surf_out_nstep,   1,        1e6);           CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _OPTIONAL_, "vec_times",         &FSLib->vec_times,        1,        1.0);           CHKERRQ(ierr); // non-dimensional
    }

    ierr = FBFreeBlocks(fb);       CHKERRQ(ierr);

    //=================================================================================
    // Load grid information from LaMEM
    //=================================================================================
    // load nx, ny, rangeX, rangeY, rangeZ
    ierr = FastScapeLoadGridInf(FSLib);           CHKERRQ(ierr);
    ierr = FastScapeCreateSurfaceGrid(FSLib, 1);  CHKERRQ(ierr);

    // save nx & ny
    if(FSLib->fs2D)
    {
        if(1 == FSLib->refine)
        {
            FSLib->nx_solve    = FSLib->extendedXNodes;
            FSLib->ny_solve    = FSLib->extendedYNodes;
        }
        else
        {
            FSLib->nx_solve    = FSLib->etRefineXNodes;
            FSLib->ny_solve    = FSLib->etRefineYNodes;
        }
    }
    else
    {
        if(1 == FSLib->refine)
        {
            if(!FSLib->non_uniform_grid)
            {
                FSLib->nx_solve    = FSLib->nx_LaMEM;
                FSLib->ny_solve    = FSLib->ny_LaMEM;
            }
            else
            {
                FSLib->nx_solve    = FSLib->msx_fs.nnodes_nug;
                FSLib->ny_solve    = FSLib->msy_fs.nnodes_nug;                
            }
        }
        else
        {
            FSLib->nx_solve    = FSLib->nx_refine;
            FSLib->ny_solve    = FSLib->ny_refine;
        }        
    }
    
    FSLib->nodes_solve = FSLib->nx_solve * FSLib->ny_solve;

    //=================================================================================
    // Visualization
    //=================================================================================
    PetscPrintf(PETSC_COMM_WORLD, "FastScape parameters: \n");
    // LaMEM grid
    PetscPrintf(PETSC_COMM_WORLD, "    Original grid: \n");  
    PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",     FSLib->nx_LaMEM, FSLib->ny_LaMEM);        
    // non uniform grid
    if(FSLib->non_uniform_grid )
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Non unifrom grid: \n");  
        PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",     FSLib->msx_fs.nnodes_nug, FSLib->msy_fs.nnodes_nug);         
    }
    // 2D extended grid
    if(FSLib->fs2D )  
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Extended  grid: \n");
        PetscPrintf(PETSC_COMM_WORLD, "    [rangeX,rangeY]       : [%g, %g] %s\n",   FSLib->extendedXRange / scal->length_fs, FSLib->extendedYRange / scal->length_fs, scal->lbl_length);
        PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",      FSLib->extendedXNodes, FSLib->extendedYNodes);
    }
    // refined grid
    if( 1 < FSLib->refine ) 
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Refined grid: \n");
        PetscPrintf(PETSC_COMM_WORLD, "    Refined times         : %d\n",       FSLib->refine);   
        // 2D
        if(FSLib->fs2D) PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n", FSLib->etRefineXNodes, FSLib->etRefineYNodes);        
        // 3D
        else    PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n", FSLib->nx_refine, FSLib->ny_refine);  
    }

    // surface process parameter
    PetscPrintf(PETSC_COMM_WORLD, "  Surface process: \n");  
    PetscPrintf(PETSC_COMM_WORLD, "    Max timestep          : %g %s\n",        FSLib->Max_dt / scal->time_fs, scal->lbl_time);
    PetscPrintf(PETSC_COMM_WORLD, "    Topography boundary   : %s\n",           FSLib->FS_BC);
    PetscPrintf(PETSC_COMM_WORLD, "    Velocity boundary     : %s\n",           FSLib->FS_VELBC);
    PetscPrintf(PETSC_COMM_WORLD, "    Sedimentation phase   : %d\n",           FSLib->sedPhases);    
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
    PetscErrorCode ierr;
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

    // for create vector   
    #define FASTSCAPE_VEC_CREATE(vec_ptr, size) \
      ierr = VecCreate(PETSC_COMM_WORLD, &(vec_ptr));  CHKERRQ(ierr);\
      ierr = VecSetSizes((vec_ptr), (size), (size));  CHKERRQ(ierr);\
      ierr = VecSetFromOptions((vec_ptr)); CHKERRQ(ierr);\

    // for collection
    // vx & vy &vz
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vz_collect);   CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vx_collect);   CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vy_collect);   CHKERRQ(ierr);

    // topography 
    FASTSCAPE_VEC_CREATE(FSLib->gtopo_fs, nodes_ori);

    // for different grid
    // non_uniform_grid
    if(FSLib->non_uniform_grid)
    {
        FASTSCAPE_VEC_CREATE(FSLib->gtopo_nug, nodes_nug);
        FASTSCAPE_VEC_CREATE(FSLib->vx_nug, nodes_nug);
        FASTSCAPE_VEC_CREATE(FSLib->vy_nug, nodes_nug);
        FASTSCAPE_VEC_CREATE(FSLib->vz_nug, nodes_nug);
    }

    if(FSLib->fs2D)
    {
        // extended part
        FASTSCAPE_VEC_CREATE(FSLib->gtopo_extend, nodes);
        FASTSCAPE_VEC_CREATE(FSLib->vx_extend, nodes);
        FASTSCAPE_VEC_CREATE(FSLib->vy_extend, nodes);
        FASTSCAPE_VEC_CREATE(FSLib->vz_extend, nodes);

        if(1 < FSLib->refine)
        {
            // extended part after refinement
            FASTSCAPE_VEC_CREATE(FSLib->gtopo_et_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vx_et_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vy_et_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vz_et_refine, nodes);
        }
    }
    else
    {
        if(1 < FSLib->refine)
        {
            // refined part
            FASTSCAPE_VEC_CREATE(FSLib->gtopo_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vx_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vy_refine, nodes);
            FASTSCAPE_VEC_CREATE(FSLib->vz_refine, nodes);
        }
    }

    // FastScape solution
    // Create all vectors, as the part for creating output is later than the part for creating surf
    FASTSCAPE_VEC_CREATE(FSLib->silt_fraction, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->basement, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->total_erosion, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->drainage_area, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->erosion_rate, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->slope, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->curvature, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->chi, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->catchment, nodes);
    FASTSCAPE_VEC_CREATE(FSLib->lake_depth, nodes);
    
    #undef FASTSCAPE_VEC_CREATE

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeLoadGridInf(FastScapeLib *FSLib)
{

    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // load global nx, ny, rangeX, rangeY, rangeZ
    FDSTAG   *fs;    
    Scaling  *scal;
    PetscInt baseNodesX, baseNodesY;
    PetscScalar bx, by, bz, ex, ey, ez;

    fs    = FSLib->surf->jr->fs;
    scal  = fs->scal;

    // range X, Y, Z   
    ierr  = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

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
        ierr = FSLoadNonUniformGrid(&FSLib->msx_fs, FSLib->rangeX_end / scal->length, fs->scal);  CHKERRQ(ierr);
        // y-direction
        ierr = FSLoadNonUniformGrid(&FSLib->msy_fs, FSLib->rangeY_end / scal->length, fs->scal);  CHKERRQ(ierr);
        FSLib->fsX.nodes      = FSLib->msx_fs.nnodes_nug;
        FSLib->fsY.nodes      = FSLib->msy_fs.nnodes_nug;
    }
    
    // 3D refine
    if(1 < FSLib->refine && 0 == FSLib->fs2D)
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
        if(1 < FSLib->refine)
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
            ms_fs->grid_spacing_min[i]  = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (ms_fs->istart[i+1] - ms_fs->istart[i]);
            ms_fs->grid_spacing_max[i]  = ms_fs->grid_spacing_min[i];
        }
        else
        {      
            // bias (last to first cell size ratio > 1 -> growing)
            bias  = ms_fs->biases[i];
        
            // average cell size
            avgSz = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (ms_fs->istart[i+1] - ms_fs->istart[i]);

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
    ms_fs->nnodes_nug  = floor( (ms_fs->xstart[ms_fs->nsegs] - ms_fs->xstart[0]) * scal->length / ms_fs->min_spacing ) + 2; 

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void GenerateGridCoordinates(PetscScalar *coords, PetscInt numNodes, PetscScalar start, PetscScalar spacing)
{
    for (PetscInt i = 0; i < numNodes; i++) 
    {
        coords[i] = start + spacing * i;
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
    PetscErrorCode ierr;

    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;
    jr      = FSLib->jr;
    ts      = jr->ts;
    step_fs = ts->istep;

    fsX->dx = FSLib->rangeX / scal->length_fs / (fsX->nodes - 1);
    fsY->dx = FSLib->rangeY / scal->length_fs / (fsY->nodes - 1);

    // original fastscape grid
    if(2 == mode || 1 == mode)
    {
        ierr = PetscMalloc(fsX->nodes * sizeof(PetscScalar), &fsX->ncoor); CHKERRQ(ierr);
        ierr = PetscMalloc(fsY->nodes * sizeof(PetscScalar), &fsY->ncoor); CHKERRQ(ierr); 
    }

    GenerateGridCoordinates(fsX->ncoor, fsX->nodes, FSLib->rangeX_begin, fsX->dx);
    GenerateGridCoordinates(fsY->ncoor, fsY->nodes, FSLib->rangeY_begin, fsY->dx);

    if(0 == FSLib->fs2D)
    {
        // 3D refined grid
        if( 1 < FSLib->refine)
        {
            if(2 == mode || 0 == step_fs)
            {
                ierr = PetscMalloc(fsX->nodes_refine * sizeof(PetscScalar), &fsX->ncoor_refine); CHKERRQ(ierr);
                ierr = PetscMalloc(fsY->nodes_refine * sizeof(PetscScalar), &fsY->ncoor_refine); CHKERRQ(ierr);        
            }

            fsX->dx_refine = FSLib->rangeX / scal->length_fs  / (fsX->nodes_refine - 1);
            fsY->dx_refine = FSLib->rangeY / scal->length_fs  / (fsY->nodes_refine - 1);

            GenerateGridCoordinates(fsX->ncoor_refine, fsX->nodes_refine, FSLib->rangeX_begin, fsX->dx_refine);
            GenerateGridCoordinates(fsY->ncoor_refine, fsY->nodes_refine, FSLib->rangeY_begin, fsY->dx_refine); 
        }
    }
    else
    {
        // 2D grid 
        if(2 == mode || 0 == step_fs)
        {
            ierr = PetscMalloc(fsX->nodes_extend * sizeof(PetscScalar), &fsX->ncoor_extend); CHKERRQ(ierr);
            ierr = PetscMalloc(fsY->nodes_extend * sizeof(PetscScalar), &fsY->ncoor_extend); CHKERRQ(ierr);
        }

        fsX->dx_extend = FSLib->extendedXRange / scal->length_fs  / (fsX->nodes_extend - 1);
        fsY->dx_extend = FSLib->extendedYRange / scal->length_fs  / (fsY->nodes_extend - 1); 

        GenerateGridCoordinates(fsX->ncoor_extend, fsX->nodes_extend, FSLib->rangeX_begin, fsX->dx_extend);
        GenerateGridCoordinates(fsY->ncoor_extend, fsY->nodes_extend, FSLib->rangeY_begin, fsY->dx_extend);

        // 2D refined grid 
        if(1 < FSLib->refine)
        {
            if(2 == mode || 0 == step_fs)
            {
                ierr = PetscMalloc(fsX->nodes_refine * sizeof(PetscScalar), &fsX->ncoor_refine); CHKERRQ(ierr);
                ierr = PetscMalloc(fsY->nodes_refine * sizeof(PetscScalar), &fsY->ncoor_refine); CHKERRQ(ierr);      
            }

            fsX->dx_refine = FSLib->extendedXRange / scal->length_fs  / (fsX->nodes_refine - 1);
            fsY->dx_refine = FSLib->extendedYRange / scal->length_fs  / (fsY->nodes_refine - 1);

            GenerateGridCoordinates(fsX->ncoor_refine, fsX->nodes_refine, FSLib->rangeX_begin, fsX->dx_refine);
            GenerateGridCoordinates(fsY->ncoor_refine, fsY->nodes_refine, FSLib->rangeY_begin, fsY->dx_refine);       
        }
    }

    if(1 == FSLib->non_uniform_grid)
    {
        // create a global coordinate for LaMEM
        if(0 == FSLib->fs2D)
        {
            if(2 == mode || 0 == step_fs)
            {
                ierr = PetscMalloc((FSLib->nx_LaMEM + 1) * sizeof(PetscScalar), &FSLib->ncoor_ori_x); CHKERRQ(ierr);
                ierr = PetscMalloc((FSLib->ny_LaMEM + 1) * sizeof(PetscScalar), &FSLib->ncoor_ori_y); CHKERRQ(ierr); 
            }

            // x-direction
            for(i = 0; i < FSLib->msx_fs.nsegs; i++)
            {   
                // generate nodal coordinates for the local part of the segment
                ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal); CHKERRQ(ierr);
            } 
            // last nodes
            FSLib->ncoor_ori_x[FSLib->nx_LaMEM] = FSLib->rangeX_end;

            // y-direction
            for(i = 0; i < FSLib->msy_fs.nsegs; i++)
            {   
                // generate nodal coordinates for the local part of the segment
                ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal); CHKERRQ(ierr);
            } 
            // last nodes
            FSLib->ncoor_ori_y[FSLib->ny_LaMEM] = FSLib->rangeY_end;   
        }
        else
        {
            if(1 == FSLib->extendedX)
            {
                if(2 == mode || 0 == step_fs)
                {
                    ierr = PetscMalloc((FSLib->ny_LaMEM + 1) * sizeof(PetscScalar), &FSLib->ncoor_ori_y); CHKERRQ(ierr);
                }

                // y-direction
                for(i = 0; i < FSLib->msy_fs.nsegs; i++)
                {   
                    // generate nodal coordinates for the local part of the segment
                    ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal); CHKERRQ(ierr);
                } 
                // last nodes
                FSLib->ncoor_ori_y[FSLib->ny_LaMEM] = FSLib->rangeY_end;                  
            }
            else
            {
                if(2 == mode || 0 == step_fs)
                {
                    ierr = PetscMalloc((FSLib->nx_LaMEM + 1) * sizeof(PetscScalar), &FSLib->ncoor_ori_x); CHKERRQ(ierr);
                }

                // x-direction
                for(i = 0; i < FSLib->msx_fs.nsegs; i++)
                {   
                    // generate nodal coordinates for the local part of the segment
                    ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal); CHKERRQ(ierr);
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

    if( _SI_ == scal->utype)
    {
        // s (LaMEM) -> yr (FastScape)
        scal->time_fs             = 1/yr;                      sprintf(scal->lbl_time_fs,         "[yr]");  
        // m (LaMEM) -> m (FastSCape)
        scal->length_fs           = 1.0;                       sprintf(scal->lbl_length_fs,       "[m]");   
        // m/s (LaMEM) -> m/yr (FastScape)
        scal->velocity_fs         = 1/yr;                      sprintf(scal->lbl_velocity_fs,     "[m/yr]"); 

        // output 
        // m^2 --> m^2
        scal->area_fs             = 1.0;                       sprintf(scal->lbl_area_fs,            "[m^2]"); 
        // m/yr --> m/yr
        scal->rate                = 1.0;                       sprintf(scal->lbl_rate,            "[m/yr]"); 
        scal->degree              = 1.0;                       sprintf(scal->lbl_degree,          "[°]");   
    }
    else if( _GEO_ == scal->utype)
    {
        // Myr (LaMEM) -> yr (FastScape)
        scal->time_fs             = 1e6;                       sprintf(scal->lbl_time_fs,         "[yr]");  
        // km (LaMEM) -> m (FastScape) 
        scal->length_fs           = km;                        sprintf(scal->lbl_length_fs,       "[m]");   
        // cm/yr (LaMEM) -> m/yr (FastScape)
        scal->velocity_fs         = m;                         sprintf(scal->lbl_velocity_fs,     "[m/yr]"); 

        // output 
        // m^2 --> km^2
        scal->area_fs             = km * km;                 sprintf(scal->lbl_area_fs,         "[km^2]"); 
        // m/yr --> km/yr
        scal->rate                = 1.0 * km;                    sprintf(scal->lbl_rate,            "[km/yr]"); 
        scal->degree              = 1.0;                       sprintf(scal->lbl_degree,          "[°]");   
    }
    else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect unit type for FastScape");

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfFastScapeCreate(FastScapeLib *FSLib, FB *fb)
{

    char filename[_str_len_];

    PetscErrorCode ierr;
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
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_fs",            &FSLib->outsurf_fs,        1, 1); CHKERRQ(ierr); 
       
    // read
    ierr = getStringParam(fb, _OPTIONAL_, "out_file_name",          filename,               "output"); CHKERRQ(ierr);
    // pvd
    ierr = getIntParam   (fb, _OPTIONAL_, "out_fs_pvd",             &FSLib->outpvd_fs,         1, 1); CHKERRQ(ierr); 

    // topography
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_topofs",        &FSLib->out_topofs,            1, 1); CHKERRQ(ierr);

    // surface parameter
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_silt_fraction", &FSLib->out_silt_fraction,     1, 1); CHKERRQ(ierr);    
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_basement",      &FSLib->out_basement,          1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_total_erosion", &FSLib->out_total_erosion,     1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_drainage_area", &FSLib->out_drainage_area,     1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_erosion_rate",  &FSLib->out_erosion_rate,      1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_slope",         &FSLib->out_slope,             1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_curvature",     &FSLib->out_curvature,         1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_chi",           &FSLib->out_chi,               1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_catchment",     &FSLib->out_catchment,         1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_lake_depth",    &FSLib->out_lake_depth,        1, 1); CHKERRQ(ierr);  

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
    sprintf(FSLib->outfile_fs,        "%s_fs",        filename);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCopyVelocity(FastScapeLib *FSLib)
{

    PetscErrorCode ierr;
    JacRes *jr;
    FreeSurf *surf;
    
    surf = FSLib->surf;
    jr = surf->jr;

    ierr = FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx); CHKERRQ(ierr);
    ierr = FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy); CHKERRQ(ierr);
    ierr = FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz); CHKERRQ(ierr);

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
    if(1.0 == bias)
    {
        // generate coordinates of local nodes
        for(i = istart; i < M + istart + 1; i++)
        {   
            ncoor[i] = xstart + (PetscScalar)( (i - istart) * avgSz );
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
			ncoor[i] = xstart + (i - istart) * begSz + dx * sum ;
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

        if((x1 <= x_coor) && (x2 >= x_coor))
        {
            idx.m = p;
            idx.mm = p + 1;

            if(FSLib->nx_LaMEM == idx.mm) idx.mm -= 1; 

            found_x = true;
            break;
        }

        find_indx = floor((x_coor - x1) / dx);

        if(0 < find_indx) p += find_indx; 
        else p++; 
    }

    if(found_x)
    {
        // y_direction
        for(q = 0; q < FSLib->ny_LaMEM; )
        {
            y1 = FSLib->ncoor_ori_y[q];
            y2 = FSLib->ncoor_ori_y[q + 1];

            if((y1 <= y_coor) && (y2 >= y_coor))
            {
                idx.n = q;
                idx.nn = q + 1;

                if(FSLib->ny_LaMEM == idx.nn) idx.nn -= 1;

                found_y = true;
                break;
            }

            find_indy = floor((y_coor - y1) / dy);

            if(0 < find_indy) q += find_indy;
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
    PetscErrorCode ierr; 
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
    }
    ierr = VecGetArray(target_vec, &value_save);          CHKERRQ(ierr); 

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
    
    ierr = VecRestoreArray(target_vec, &value_save);          CHKERRQ(ierr); 
    
    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode InterpolationFor2DNonUniformGrid(FastScapeLib *FSLib, PetscScalar *value_in, PetscScalar *value_out)
{
    FSGrid *fsGrid;    
    PetscReal *coords_ori, *coords_ext;    
    PetscInt nodes_ori, nodes_ext;    
    PetscFunctionBegin;
    MeshSeg1DFS ms;
    PetscInt i, j, m, n, mm, nn, ind, ind_a, ind_b, p, q, find_m = 0, find_n = 0, find_indx = 0, find_indy = 0, idx1, idx2, find_ind;
    PetscScalar dx, dy, wtx, wty, x1, x2, y1, y2, x_coor, y_coor, coord1, coord2;

    // 选择插值方向    
    if (FSLib->extendedX) 
    {
        fsGrid = &FSLib->fsY;        
        coords_ori = FSLib->ncoor_ori_y;
        nodes_ori = FSLib->ny_LaMEM;
        coords_ext = fsGrid->ncoor_extend;
        nodes_ext = fsGrid->nodes_extend;
        ms        = FSLib->msy_fs;
    } 
    else 
    {
        fsGrid = &FSLib->fsX;
        coords_ori = FSLib->ncoor_ori_x;
        nodes_ori = FSLib->nx_LaMEM;
        coords_ext = fsGrid->ncoor_extend;
        nodes_ext = fsGrid->nodes_extend;
        ms        = FSLib->msx_fs;
    }

    // interpolate after calculating mean value
        for(i = 0; i < nodes_ext; i++)
        {
            x_coor = coords_ext[i];

            // get nearest 2 index
            for(p = 0; p < nodes_ori; )
            {
                x1 = coords_ori[p];
                x2 = coords_ori[p + 1];

                if( (x1 <= x_coor) && (x2 >= x_coor) )
                {
                    m      = p;
                    mm     = p + 1;
                    if( nodes_ori == mm )  mm -= 1;
                    find_m = 1;                   
                    break;
                }   
                   
                find_indx = floor( (x_coor - x1) / ms.max_spacing );
                if( 0 < find_indx)
                {
                    p += find_indx;
                }
                else
                {
                    p++;
                }
            }

            if( 1 == find_m )
            {
                // interpolate
                ind = i;    
                ind_a = m;
                ind_b = mm;                          
                // bilinear interpolation
                dx    = x2 - x1;
                wtx   = x_coor - x1;
                
                if(m == mm)
                {
                    value_out[ind] = value_in[ind_a];
                }
                else
                {
                    value_out[ind] = ( (dx - wtx) / dx ) * value_in[ind_a]+ ( wtx / dx ) * value_in[ind_b];
                }
            } 

            find_m = 0;
        }

    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode GatherVariableFromLaMEM(FastScapeLib *FSLib, PetscScalar *topo_alloc, PetscScalar *vx_alloc, PetscScalar *vy_alloc, PetscScalar *vz_alloc, PetscInt step_fs)
{
    PetscErrorCode ierr;
    PetscInt L, sx, sy, sz, nx, ny, nz, tproc, rankZ_id, i, j;
    PetscScalar ***topo;
    PetscScalar ***vz, ***vx, ***vy;
    PetscScalar ***vz_collect, ***vx_collect, ***vy_collect;
    PetscScalar *topo_collect = PETSC_NULL;
    PetscScalar *vz_LaMEM = PETSC_NULL, *vx_LaMEM = PETSC_NULL, *vy_LaMEM = PETSC_NULL;
    FDSTAG      *fs;
    FreeSurf    *surf;

    surf   = FSLib->surf;
    fs     = surf->jr->fs;

    // local process info
    L    = (PetscInt)fs->dsz.rank;
    ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz);    CHKERRQ(ierr);

    // Gather topography and velocity
    // topography
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);         CHKERRQ(ierr);
    if(0 == step_fs)
    {
        ierr = VecScatterCreateToZero(surf->gtopo, &FSLib->ctx, &FSLib->gtopo_collect);                        CHKERRQ(ierr); 
        ierr = VecScatterBegin(FSLib->ctx, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(FSLib->ctx, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);
    }

    // velocity
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx, &vx);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect);     CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy, &vy);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect);     CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz, &vz);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect);     CHKERRQ(ierr);

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
    ierr = VecScatterCreateToZero(FSLib->vx_collect, &FSLib->ctx, &FSLib->vx_LaMEM);                            CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vx_collect, FSLib->vx_LaMEM, INSERT_VALUES, SCATTER_FORWARD);     CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vx_collect, FSLib->vx_LaMEM, INSERT_VALUES, SCATTER_FORWARD);       CHKERRQ(ierr);
    ierr = VecScatterCreateToZero(FSLib->vy_collect, &FSLib->ctx, &FSLib->vy_LaMEM);                            CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vy_collect, FSLib->vy_LaMEM, INSERT_VALUES, SCATTER_FORWARD);     CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vy_collect, FSLib->vy_LaMEM, INSERT_VALUES, SCATTER_FORWARD);       CHKERRQ(ierr);
    ierr = VecScatterCreateToZero(FSLib->vz_collect, &FSLib->ctx, &FSLib->vz_LaMEM);                            CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vz_collect, FSLib->vz_LaMEM, INSERT_VALUES, SCATTER_FORWARD);     CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vz_collect, FSLib->vz_LaMEM, INSERT_VALUES, SCATTER_FORWARD);       CHKERRQ(ierr);

    // reallocate
    tproc = fs->dsx.nproc * fs->dsy.nproc * fs->dsz.nproc;
    rankZ_id = NO_NEED;    

    if (0 == fs->dsz.rank)  rankZ_id = fs->dsx.nproc * fs->dsy.rank + fs->dsx.rank;

    // local process info  
    ProcInfo local_info = {sx, sy, sz, nx, ny, nz, rankZ_id};    
    ProcInfo *global_infos = nullptr;    

    // allocate memory    
    if (ISRankZero(PETSC_COMM_WORLD)) 
    {        
        ierr = PetscMalloc(tproc * sizeof(ProcInfo), &global_infos); CHKERRQ(ierr);
        if (!global_infos)  SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_MEM, "Memory allocation failed");        
    } 

    // collect information for the surface grid   
    ierr = MPI_Gather(&local_info, sizeof(ProcInfo), MPI_BYTE, global_infos, sizeof(ProcInfo), MPI_BYTE, 0, PETSC_COMM_WORLD);    CHKERRQ(ierr);
    if (ISRankZero(PETSC_COMM_WORLD)) 
    {        
        if(0 == step_fs)
        {
            ierr = VecGetArray(FSLib->gtopo_collect, &topo_collect);  CHKERRQ(ierr);
        }
        ierr = VecGetArray(FSLib->vx_LaMEM,  &vx_LaMEM);      CHKERRQ(ierr);
        ierr = VecGetArray(FSLib->vy_LaMEM,  &vy_LaMEM);      CHKERRQ(ierr);        
        ierr = VecGetArray(FSLib->vz_LaMEM,  &vz_LaMEM);      CHKERRQ(ierr);

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

                    if (0 == step_fs)   topo_alloc[target_idx] = topo_collect[global_index];
                    vx_alloc[target_idx] = vx_LaMEM[global_index];
                    vy_alloc[target_idx] = vy_LaMEM[global_index];
                    vz_alloc[target_idx] = vz_LaMEM[global_index];

                    global_index++;
                }
            }
        }

        PetscFree(global_infos);
        if(0 == step_fs)
        {
            ierr = VecRestoreArray(FSLib->gtopo_collect, &topo_collect);      CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(FSLib->vz_LaMEM, &vz_LaMEM);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->vx_LaMEM, &vx_LaMEM);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->vy_LaMEM, &vy_LaMEM);       CHKERRQ(ierr);
    }

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo);                CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx, &vx);                     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy, &vy);                     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect);    CHKERRQ(ierr);         
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz, &vz);                     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect);    CHKERRQ(ierr);    
    ierr = VecScatterDestroy(&FSLib->ctx);                                        CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeStretchGrid(FastScapeLib *FSLib)
{
    PetscErrorCode ierr;
    PetscScalar Exx, Eyy, Ezz, Rxx, Ryy, Rzz, step;
    Scaling  *scal;
    JacRes   *jr;

    jr    = FSLib->jr;
    scal  = FSLib->scal;
    step  = jr->ts->dt;

    // update model range
    // range X, Y, Z   
    ierr  = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, &Ezz, NULL, NULL, NULL, &Rxx, &Ryy, &Rzz); CHKERRQ(ierr);

    auto update_range = [&](PetscScalar& pos, PetscScalar strain_rate, PetscScalar rot_center) 
    {           
        pos += step * strain_rate * (pos / scal->length - rot_center) * scal->length;
    };

    update_range(FSLib->rangeX_begin, Exx, Rxx);
    update_range(FSLib->rangeX_end,   Exx, Rxx);
    update_range(FSLib->rangeY_begin, Eyy, Ryy);
    update_range(FSLib->rangeY_end,   Eyy, Ryy);
    update_range(FSLib->rangeZ_begin, Ezz, Rzz);
    update_range(FSLib->rangeZ_end,   Ezz, Rzz);

    FSLib->rangeX       = (FSLib->rangeX_end - FSLib->rangeX_begin) * scal->length_fs; //(km) in LaMEM to (m) in FastScape (GEO)
    FSLib->rangeY       = (FSLib->rangeY_end - FSLib->rangeY_begin) * scal->length_fs;

    if(FSLib->fs2D)
    {
        // 2D No Refine
        // extend grid in FastScape
        // extend in rangeX
        if(FSLib->extendedX)
        {
            if(Eyy) FSLib->extendedXRange += (step * Exx * (FSLib->extendedXRange / scal->length - Rxx)) * scal->length * scal->length_fs;
            else    FSLib->extendedXRange += FSLib->rangeX_begin * scal->length_fs;         
            FSLib->extendedYRange = FSLib->rangeY;
        }
        else
        {
            if(Exx) FSLib->extendedYRange += (step * Eyy * (FSLib->extendedYRange / scal->length - Ryy)) * scal->length * scal->length_fs;
            else    FSLib->extendedYRange += FSLib->rangeY_begin * scal->length_fs;                
            FSLib->extendedXRange = FSLib->rangeX;
        }
    }

    if(FSLib->non_uniform_grid)
    { 
        auto update_coords = [&](PetscScalar* coords, int count, PetscScalar strain_rate, PetscScalar rot_center) 
        {            
            for (int idx = 0; idx < count; idx++) 
            {                
                coords[idx] += step * strain_rate * (coords[idx] - rot_center);            
            }        
        };

        if(!FSLib->fs2D)
        {
            update_coords(FSLib->msx_fs.xstart, FSLib->msx_fs.nsegs + 1, Exx, Rxx);
            update_coords(FSLib->msy_fs.xstart, FSLib->msy_fs.nsegs + 1, Eyy, Ryy);
        }   
        else
        {
            if(FSLib->extendedX)  update_coords(FSLib->msy_fs.xstart, FSLib->msy_fs.nsegs + 1, Eyy, Ryy); // y-direction
            else  update_coords(FSLib->msx_fs.xstart, FSLib->msx_fs.nsegs + 1, Exx, Rxx); // x-direction
        }
    }

    // update coordinate
    ierr = FastScapeCreateSurfaceGrid(FSLib, 0); CHKERRQ(ierr);

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
    PetscInt i, j, ind, ind_a, ind_b, ind_aa, ind_bb, countX1, countX2, countY, tnodes_refine;
    PetscScalar distance_ax, distance_ay, unit_factor;

    tnodes_refine = nx_refine * ny_refine;
    countX1       = 0;
    countX2       = 0;
    countY        = 0;
    ind_a         = countX1 * FSLib->refine;
    ind_b         = (countX1 + 1) * FSLib->refine;
    ind_aa        = countX2;
    ind_bb        = countX2 + 1;

    // scaling factor
    switch (corMode) 
    {        
        case 1: // topography：km (LaMEM) → m (FastScape)
            unit_factor = scal->length * scal->length_fs;
            break;
        case 2: // velocity：cm/yr (LaMEM) → m/yr (FastScape)
            unit_factor = scal->velocity / scal->velocity_fs;
            break;
        default: // no transition
            unit_factor = 1.0;
    }

    // interpolate in x-direction
    for(j = 0; j < ny_refine; j += FSLib->refine) 
    {   
        for(i = 0; i < nx_refine; i++) 
        {
            ind = j * nx_refine + i;

            if(0 == ind % nx_refine) 
            { 
                countX1 = 0;   

                if(0 == ind)
                {
                    data_refine[ind_a] = data[ind_aa] * unit_factor;
                    data_refine[ind_b] = data[ind_bb] * unit_factor;                                        
                }
                else
                {
                    ind_a  = countX1 * FSLib->refine + j * nx_refine;
                    ind_b  = (countX1+1) * FSLib->refine + j * nx_refine; 
                    
                    countX2++; 
                    
                    ind_aa = countX2;
                    ind_bb = countX2 + 1;
 
                    data_refine[ind_a] = data[ind_aa] * unit_factor;
                    data_refine[ind_b] = data[ind_bb] * unit_factor;  
                }

                data_refine_pass[ind_a] = data_refine[ind_a];
                data_refine_pass[ind_b] = data_refine[ind_b]; 
            }
            else
            {
                if(ind == ind_b) 
                {
                    if(ind == tnodes_refine-1) 
                    {
                        data_refine[ind_b] = data[ind_bb] * unit_factor;
                        goto skip;
                    }
                    else
                    {
                        countX1++;
                        countX2++;

                        ind_a  = countX1 * FSLib->refine + j * nx_refine;
                        ind_b  = (countX1+1) * FSLib->refine + j * nx_refine;
                        ind_aa = countX2;
                        ind_bb = countX2 + 1;

                        data_refine[ind_b] = data[ind_bb] * unit_factor;  
                    }

                    data_refine_pass[ind_b]   = data_refine[ind_b];
                }
                else
                {
                    distance_ax           = ((ind - ind_a) % nx_refine) * 1.0 / FSLib->refine;   
                    data_refine[ind]      = data_refine[ind_a] * (1 - distance_ax)  + data_refine[ind_b] * distance_ax;               
                    data_refine_pass[ind] = data_refine[ind];
                }
            }                   
        }
    }
    skip: // printf("interpolation in x-direction done\n");     

    // interpolate in y-direction
    for(j = 0; j < ny_refine; j++) 
    {   
        if(0 == j % FSLib->refine) 
        {
            if( 0 != j) countY += FSLib->refine;    
            continue;
        }

        for(i = 0; i < nx_refine; i++) 
        {
            ind   = j * nx_refine + i;
            ind_a = countY * nx_refine + i;
            ind_b = (countY + FSLib->refine) * nx_refine + i;

            distance_ay           = (j - countY) * 1.0 / FSLib->refine;  
            data_refine[ind]      = data_refine[ind_a] * (1 - distance_ay)  + data_refine[ind_b] * distance_ay;                   
            data_refine_pass[ind] = data_refine[ind]; 
        }         
    }   

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Extended2D(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_extended, PetscScalar *data_extend_pass, Scaling *scal, PetscInt corMode, PetscBool isRefine)
{
    PetscErrorCode ierr;
    PetscScalar sum, unit_factor    = 1.0;
    PetscInt i, j, idx, linearIdx; 
    PetscScalar *data_aver_ori = PETSC_NULL;
    PetscScalar *data_aver     = PETSC_NULL;
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

    ierr = PetscMalloc(origDim * sizeof(PetscScalar), &data_aver_ori);  CHKERRQ(ierr);
    ierr = PetscMalloc(targetDim * sizeof(PetscScalar), &data_aver);    CHKERRQ(ierr);

    // calculate averange values with scaling
    for (i = 0; i < origDim; i++) 
    {
        sum = 0.0;
        
        for (j = 0; j < otherDim; j++) 
        {
            idx = extendX ? i * FSLib->nx_LaMEM + j : j * FSLib->nx_LaMEM + i;
            sum += data[idx] * unit_factor;
        }

        data_aver_ori[i] = sum / otherDim;
    }

    // non-uniform grid
    if (FSLib->non_uniform_grid) 
    {
        ierr = InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver); CHKERRQ(ierr);
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

    ierr = PetscFree(data_aver_ori); CHKERRQ(ierr);
    ierr = PetscFree(data_aver); CHKERRQ(ierr);

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

    PetscErrorCode ierr;
    PetscScalar dt, dt_scal, dt_fs, time_fs, Exx, Eyy, Rxx, Ryy;
    PetscInt i, j, tnodes, tnodes_ori, L, sx, sy, sz, nx, ny, nz, step_fs;
    PetscScalar *topo_fs = PETSC_NULL, *topo_solve = PETSC_NULL, *topo_solve_refine = PETSC_NULL, *topo_alloc = PETSC_NULL; 
    PetscScalar *vx_alloc = PETSC_NULL, *vy_alloc = PETSC_NULL, *vz_alloc = PETSC_NULL;
    PetscScalar *vx_solve = PETSC_NULL, *vy_solve = PETSC_NULL, *vz_solve = PETSC_NULL;
    PetscScalar *vx_solve_refine = PETSC_NULL, *vy_solve_refine = PETSC_NULL, *vz_solve_refine = PETSC_NULL; 
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

    ierr = PetscMalloc(tnodes_ori * sizeof(PetscScalar), &topo_fs); CHKERRQ(ierr);
    ierr = PetscMalloc(tnodes_ori * sizeof(PetscScalar), &topo_alloc); CHKERRQ(ierr);
    ierr = PetscMalloc(tnodes_ori * sizeof(PetscScalar), &vx_alloc);   CHKERRQ(ierr);
    ierr = PetscMalloc(tnodes_ori * sizeof(PetscScalar), &vy_alloc);   CHKERRQ(ierr);
    ierr = PetscMalloc(tnodes_ori * sizeof(PetscScalar), &vz_alloc);   CHKERRQ(ierr);  

    ierr = GatherVariableFromLaMEM(FSLib, topo_alloc, vx_alloc, vy_alloc, vz_alloc, step_fs);  CHKERRQ(ierr);        

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        PetscScalar dt_max = FSLib->Max_dt; // Maximum step length, if dt_LaMEM is larger than this, use this
        PetscScalar dt_n = 0; //dt_residual
        PetscScalar quotient = dt_fs/dt_max;
        PetscInt nsteps = floor( dt_fs/dt_max );
        PetscScalar *topo_pass = PETSC_NULL, *vx_pass = PETSC_NULL, *vy_pass = PETSC_NULL, *vz_pass = PETSC_NULL; 

        // store the phase that is being sedimented
        surf->phase = FSLib->sedPhases;

        tnodes      = FSLib->nodes_solve;  

        ierr = PetscMalloc(tnodes * sizeof(PetscScalar), &topo_pass); CHKERRQ(ierr);
        ierr = PetscMalloc(tnodes * sizeof(PetscScalar), &vx_pass);   CHKERRQ(ierr);
        ierr = PetscMalloc(tnodes * sizeof(PetscScalar), &vy_pass);   CHKERRQ(ierr);
        ierr = PetscMalloc(tnodes * sizeof(PetscScalar), &vz_pass);   CHKERRQ(ierr);             

        // timestep
        if(1 > nsteps) 
        {
            nsteps = 1;
            dt_max = dt_fs;
        }
        else if( (1 <= nsteps) && (quotient != nsteps) )
        {
            nsteps = nsteps + 1;
            dt_n   = dt_fs - (nsteps - 1) * dt_max;
        }

        // get background strain rates
        ierr = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, NULL, NULL, NULL, NULL, &Rxx, &Ryy, NULL);  CHKERRQ(ierr);

        if( Exx || Eyy )
        {
            ierr = FastScapeStretchGrid(FSLib);       CHKERRQ(ierr);  
        }

        // run fastscape when using 3D geodynamic model
        if(0 == FSLib->fs2D)
        {
            // for non uniform grid
            if(FSLib->non_uniform_grid )
            {
                if( 0 == step_fs)
                {
                    ierr = InterpolationFor3DNonUniformGrid(FSLib, topo_alloc, 1);              CHKERRQ(ierr);  
                }
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vx_alloc, 2);                    CHKERRQ(ierr);             
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vy_alloc, 3);                    CHKERRQ(ierr);  
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vz_alloc, 4);                    CHKERRQ(ierr);  
            }               

            // don't apply refinement
            if( 1 == FSLib->refine)
            {
                if(0 == FSLib->non_uniform_grid)
                {
                    ierr = VecGetArray(FSLib->gtopo_fs, &topo_fs);                        CHKERRQ(ierr);

                    for(i = 0; i < tnodes; i++) 
                    {   
                            if (0 == step_fs)   topo_pass[i] = topo_alloc[i] * scal->length * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
                            else    topo_pass[i] = topo_fs[i] * scal->length_fs;
                            vx_pass[i]   = vx_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                            vy_pass[i]   = vy_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)                            
                            vz_pass[i]   = vz_alloc[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                         
                    }

                    ierr = VecRestoreArray(FSLib->gtopo_fs, &topo_fs);                    CHKERRQ(ierr);                     
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);                    
                    ierr = VecGetArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                    
                    for(i = 0; i < tnodes; i++) 
                    {   
                            topo_pass[i] = topo_solve[i] * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
                            vx_pass[i]   = vx_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                            vy_pass[i]   = vy_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)  
                            vz_pass[i]   = vz_solve[i] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                    }

                    ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);                    
                    ierr = VecRestoreArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                }                      
            }                           
            // apply refinement
            else 
            {
                printf("Refined times                    : %d\n", FSLib->refine);
                printf("Refined grid cells [nx, ny]      : [%d, %d] \n", FSLib->nx_refine, FSLib->ny_refine);

                ierr = VecGetArray(FSLib->gtopo_refine, &topo_solve_refine);  CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vx_refine, &vx_solve_refine);       CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vy_refine, &vy_solve_refine);       CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vz_refine, &vz_solve_refine);       CHKERRQ(ierr);
    
                if(0 == FSLib->non_uniform_grid)
                {
                    // load velocity and topography field
                    if(0 == step_fs)
                    {
                        ierr = BilinearInterpolate(FSLib, topo_alloc, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(i = 0; i < tnodes; i++)
                        {
                                topo_pass[i] = topo_solve_refine[i] * scal->length_fs;
                            
                        }
                    }
                    ierr = BilinearInterpolate(FSLib, vx_alloc, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_alloc, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_alloc, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
    
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);
                    
                    // load velocity and topography field
                    if(0 == step_fs)
                    {
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(i = 0; i < tnodes; i++)
                        {

                                topo_pass[i] = topo_solve_refine[i] * scal->length_fs;
                            
                        }
                    }

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                                        
                    ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr); 
                }

                ierr = VecRestoreArray(FSLib->gtopo_refine, &topo_solve_refine);  CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vx_refine, &vx_solve_refine);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vy_refine, &vy_solve_refine);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vz_refine, &vz_solve_refine);       CHKERRQ(ierr);
            }
        }
        // run fastscape when using 2D geodynamic model
        else
        {
            ierr = VecGetArray(FSLib->gtopo_extend, &topo_solve);  CHKERRQ(ierr);
            ierr = VecGetArray(FSLib->vx_extend, &vx_solve);       CHKERRQ(ierr);
            ierr = VecGetArray(FSLib->vy_extend, &vy_solve);       CHKERRQ(ierr); 
            ierr = VecGetArray(FSLib->vz_extend, &vz_solve);       CHKERRQ(ierr); 

            if(0 == FSLib->non_uniform_grid)
            {
                // don't apply refinement
                if(1 == FSLib->refine)
                {
                    // load velocity and topography field
                    // step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
                    // topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1, PETSC_FALSE);  CHKERRQ(ierr);
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
                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2, PETSC_FALSE);        CHKERRQ(ierr);
                        std::fill(vx_pass, vx_pass + tnodes, 0);
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2, PETSC_FALSE);        CHKERRQ(ierr);    
                        std::fill(vy_pass, vy_pass + tnodes, 0);
                    }
                    // vz
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2, PETSC_FALSE);            CHKERRQ(ierr);

                }
                // apply refinement
                else
                {
                    // load velocity and topography field
                    ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine);   CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_et_refine, &vx_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_et_refine, &vy_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_et_refine, &vz_solve_refine);        CHKERRQ(ierr);

                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, NULL, NULL, 0, PETSC_TRUE);        CHKERRQ(ierr);
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(i = 0; i < tnodes; i++)
                        {
                                topo_pass[i] = topo_solve_refine[i] * scal->length_fs;
                        }
                    }

                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, NULL, NULL, 0, PETSC_TRUE);           CHKERRQ(ierr);
                        std::fill(vx_solve, vx_solve + tnodes, 0);
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, NULL, NULL, 0, PETSC_TRUE);           CHKERRQ(ierr);   
                        std::fill(vy_solve, vy_solve + tnodes, 0);
                    }
                    
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, NULL, NULL, 0, PETSC_TRUE);               CHKERRQ(ierr);

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);

                    ierr = VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine);       CHKERRQ(ierr);            

                }
            }
            else
            {
                // don't apply refinement
                if(1 == FSLib->refine)
                {
                    // load velocity and topography field
                    // step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
                    // topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1, PETSC_FALSE); CHKERRQ(ierr);
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
                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2, PETSC_FALSE);      CHKERRQ(ierr);
                        std::fill(vx_pass, vx_pass + tnodes, 0);
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2, PETSC_FALSE);      CHKERRQ(ierr);    
                        std::fill(vy_pass, vy_pass + tnodes, 0);
                    }
                    // vz
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2, PETSC_FALSE);          CHKERRQ(ierr);
                }
                // apply refinement
                else
                {
                    // load velocity and topography field
                    ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine);   CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_et_refine, &vx_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_et_refine, &vy_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_et_refine, &vz_solve_refine);        CHKERRQ(ierr);

                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, NULL, NULL, 0, PETSC_TRUE);     CHKERRQ(ierr);
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(i = 0; i < tnodes; i++) 
                        {
                                topo_pass[i] = topo_solve_refine[i] * scal->length_fs;
                         
                        }
                    }

                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, NULL, NULL, 0, PETSC_TRUE);        CHKERRQ(ierr);
                        std::fill(vx_solve, vx_solve + tnodes, 0);
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, NULL, NULL, 0, PETSC_TRUE);        CHKERRQ(ierr);   
                        std::fill(vy_solve, vy_solve + tnodes, 0);
                    }
                    
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, NULL, NULL, 0, PETSC_TRUE);            CHKERRQ(ierr);

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);

                    ierr = VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine);       CHKERRQ(ierr);
                }
            }

            ierr = VecRestoreArray(FSLib->gtopo_extend, &topo_solve);  CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vx_extend, &vx_solve);       CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vy_extend, &vy_solve);       CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vz_extend, &vz_solve);       CHKERRQ(ierr);       
        }
     
        // run FastScape
        ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
        
        ierr = VecGetArray(FSLib->gtopo_fs, &topo_fs);                 CHKERRQ(ierr);

        // pass data to original LaMEM grid
        if(0 == FSLib->fs2D)
        {
            ierr = PassValue3D(FSLib, topo_pass, topo_fs);             CHKERRQ(ierr); 
        }
        else
        {
            ierr = PassValue2D(FSLib, topo_pass, topo_fs);             CHKERRQ(ierr); 
        }

        // save data in a new grid used by fastscape
        if(0 == step_fs || 0 == (step_fs + 1) % FSLib->surf_out_nstep)
        {
            ierr = FastScapeSave(FSLib, step_fs, time_fs);             CHKERRQ(ierr);
        }

        ierr = PetscFree(topo_pass); CHKERRQ(ierr);
        ierr = PetscFree(vx_pass);   CHKERRQ(ierr);
        ierr = PetscFree(vy_pass);   CHKERRQ(ierr);
        ierr = PetscFree(vz_pass);   CHKERRQ(ierr);
    }
       
    // Broadcast      
    if(ISParallel(PETSC_COMM_WORLD))
    {
        ierr = MPI_Bcast(topo_fs, (PetscMPIInt)tnodes_ori, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD);   CHKERRQ(ierr);    
    }

    L    = (PetscInt)fs->dsz.rank;
    ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo, &topo);       CHKERRQ(ierr);

    // Save topography in different rank
    START_PLANE_LOOP
    {
        topo[L][j][i] = topo_fs[j * FSLib->nx_LaMEM + i] / scal->length; // (km)
    }
    END_PLANE_LOOP  

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo);  CHKERRQ(ierr);
    ierr = VecRestoreArray(FSLib->gtopo_fs, &topo_fs);              CHKERRQ(ierr);

    // compute ghosted version of the topography
    GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

    // compute & store average topography
    ierr = FreeSurfGetAvgTopo(surf);                                CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
    
    ierr = PetscFree(topo_alloc); CHKERRQ(ierr);
    ierr = PetscFree(vx_alloc);   CHKERRQ(ierr);
    ierr = PetscFree(vy_alloc);   CHKERRQ(ierr);
    ierr = PetscFree(vz_alloc);   CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue2D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
    PetscErrorCode ierr;
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

    ierr = VecGetArray(FSLib->gtopo_extend, &topo_extend); CHKERRQ(ierr);

    mainNodes = FSLib->extendedX ? FSLib->extendedYNodes : FSLib->extendedXNodes;
    secNodes = FSLib->extendedX ? FSLib->extendedXNodes : FSLib->extendedYNodes;
    laMEMNodes = FSLib->extendedX ? FSLib->ny_LaMEM : FSLib->nx_LaMEM;

    ierr = PetscMalloc1(mainNodes, &topo_aver); CHKERRQ(ierr);
    ierr = PetscMalloc1(laMEMNodes, &topo_aver_ori); CHKERRQ(ierr);

    // refine
    if (FSLib->refine > 1) 
    {
        ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_et_refine); CHKERRQ(ierr);

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
        ierr = VecRestoreArray(FSLib->gtopo_et_refine, &topo_et_refine); CHKERRQ(ierr);
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

        topo_aver[idx] = sum / secNodes;
    }

    // non-uniform grid
    if (FSLib->non_uniform_grid)
    {
        for (idx = 0; idx < laMEMNodes; idx++) 
        {
            coord = FSLib->extendedX ? FSLib->ncoor_ori_y[idx] : FSLib->ncoor_ori_x[idx];
            begin = FSLib->extendedX ? FSLib->ncoor_ori_y[0] : FSLib->ncoor_ori_x[0];
            delta = FSLib->extendedX ? fsY->dx : fsX->dx;

            n = PetscFloorReal((coord - begin) / delta);
            nn = n + 1;

            if (nn >= mainNodes) nn = n;
            if (n == nn)  topo_aver_ori[idx] = topo_aver[n];
            else 
            {
                coord0 = begin + n * delta;
                coord1 = begin + nn * delta;
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

    ierr = VecRestoreArray(FSLib->gtopo_extend, &topo_extend); CHKERRQ(ierr);
    ierr = PetscFree(topo_aver); CHKERRQ(ierr);
    ierr = PetscFree(topo_aver_ori); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue3D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
    PetscErrorCode ierr;
    FSGrid  *fsX;
    FSGrid  *fsY;
    Scaling *scal;
    PetscInt m, n, mm, nn, ind_a, ind_b, ind_c, ind_d, refine, i, j, ind, ind2, countX = 0, countY = 0;
    PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x_begin = 0.0, y_begin = 0.0;
    PetscScalar *topo_refine = PETSC_NULL, *topo_nug = PETSC_NULL;

    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;
    refine  = FSLib->refine;

    // for scaling    
    auto convertValue = [&](PetscScalar value) -> PetscScalar 
    {
        PetscScalar converted = value / scal->length_fs;
        if (converted > FSLib->rangeZ_end) return FSLib->rangeZ_end;
        if (converted < FSLib->rangeZ_begin) return FSLib->rangeZ_begin;
        return converted;
    };

    if(1 == FSLib->non_uniform_grid)
    {
        x_begin = FSLib->ncoor_ori_x[0];
        y_begin = FSLib->ncoor_ori_y[0];
    }

    if(0 == FSLib->non_uniform_grid)
    {
        if(1 == refine)
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
            ierr = VecGetArray(FSLib->gtopo_refine, &topo_refine);         CHKERRQ(ierr);

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

                    if (i == FSLib->nx_refine - refine) {   countX = 0;     countY++;} 
                    else    countX++;           
                }
            }
            ierr = VecRestoreArray(FSLib->gtopo_refine,  &topo_refine);    CHKERRQ(ierr);
        } 
    }
    else
    {
        ierr = VecGetArray(FSLib->gtopo_nug, &topo_nug);    CHKERRQ(ierr);

        if(1 == refine)
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
                    m  = floor( (x_coor - x_begin) / dx );
                    mm = (m + 1 < fsX->nodes) ? m + 1 : m;

                    // y-direction
                    n  = floor( (y_coor - y_begin) / dy );
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
            ierr = VecGetArray(FSLib->gtopo_refine, &topo_refine);         CHKERRQ(ierr);

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

                    if (i == FSLib->nx_refine - refine) { countX = 0; countY++;} 
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
                    m  = floor( (x_coor - x_begin) / dx );
                    mm = (m + 1 < fsX->nodes) ? m + 1 : m;

                    // y-direction
                    n  = floor( (y_coor - y_begin) / dy );
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
            ierr = VecRestoreArray(FSLib->gtopo_refine,  &topo_refine);    CHKERRQ(ierr);
        } 
        ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_nug);    CHKERRQ(ierr);
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
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PetscInt i, nx, ny, numFields;
    PetscScalar *silt_fraction = PETSC_NULL, *basement = PETSC_NULL, *total_erosion = PETSC_NULL;
    PetscScalar *drainage_area = PETSC_NULL, *erosion_rate = PETSC_NULL, *slope = PETSC_NULL;
    PetscScalar *curvature = PETSC_NULL, *chi = PETSC_NULL, *catchment = PETSC_NULL, *lake_depth = PETSC_NULL;

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
    ierr = PetscSNPrintf(fname, sizeof(fname), "%s/%s_p0.vts", dirName, FSLib->outfile_fs); CHKERRQ(ierr);

    // 打开文件
    fp = fopen(fname, "wb");
    if (fp == NULL)     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open file %s", fname);

    WriteXMLHeader(fp, "StructuredGrid");

    // grid info
    fprintf(fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld 1 1\">\n",
            (LLD)1, (LLD)nx, (LLD)1, (LLD)ny);
    fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld 1 1\">\n",
            (LLD)1, (LLD)nx, (LLD)1, (LLD)ny);
    fprintf(fp, "\t\t\t<CellData/>\n");  

    // offset & nodes
    fprintf(fp, "\t\t<Points>\n");
    fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset);
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
                "NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
                fields[i].name, fields[i].unit, (LLD)offset);
            offset += sizeof(uint64_t) + sizeof(float) * gridSize;
        }
    }

    fprintf(fp, "\t\t</PointData>\n");
    fprintf(fp, "\t\t</Piece>\n");
    fprintf(fp, "\t</StructuredGrid>\n");

    fprintf(fp, "\t<AppendedData encoding=\"raw\">\n_");
    // write point coordinates
    // allocate output buffer
    ierr = PetscMalloc((size_t)(_max_num_comp_surf_ * nx * ny)*sizeof(float), &FSLib->buff_fs);    CHKERRQ(ierr);

    ierr = PVSurfWriteCoordFS (FSLib, fp, topo, mode); CHKERRQ(ierr);

    // topography
    if(FSLib->out_topofs)
    {
        ierr = PVSurfWriteInfFS  (FSLib, fp, topo, 1);                      CHKERRQ(ierr);
    }
    // silt fraction
    if(FSLib->out_silt_fraction) 
    {
        ierr = VecGetArray(FSLib->silt_fraction, &silt_fraction);       CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, silt_fraction, 2);         CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->silt_fraction, &silt_fraction);   CHKERRQ(ierr);
    }
    // basement
    if(FSLib->out_basement) 
    {
        ierr = VecGetArray(FSLib->basement, &basement);                 CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, basement, 3);              CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->basement, &basement);             CHKERRQ(ierr);
    }
    // total_erosion
    if(FSLib->out_total_erosion) 
    {
        ierr = VecGetArray(FSLib->total_erosion, &total_erosion);       CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, total_erosion, 4);         CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->total_erosion, &total_erosion);   CHKERRQ(ierr);
    }
    // drainage_area
    if(FSLib->out_drainage_area) 
    {
        ierr = VecGetArray(FSLib->drainage_area, &drainage_area);       CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, drainage_area, 5);         CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->drainage_area, &drainage_area);   CHKERRQ(ierr);
    }
    // erosion_rate
    if(FSLib->out_erosion_rate) 
    {
        ierr = VecGetArray(FSLib->erosion_rate, &erosion_rate);         CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, erosion_rate, 6);          CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->erosion_rate, &erosion_rate);     CHKERRQ(ierr);
    }
    // slope
    if(FSLib->out_slope) 
    {
        ierr = VecGetArray(FSLib->slope, &slope);                       CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, slope, 7);                 CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->slope, &slope);                   CHKERRQ(ierr);
    }
    // curvature
    if(FSLib->out_curvature) 
    {
        ierr = VecGetArray(FSLib->curvature, &curvature);               CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, curvature, 8);             CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->curvature, &curvature);           CHKERRQ(ierr);
    }
    // chi
    if(FSLib->out_chi) 
    {
        ierr = VecGetArray(FSLib->chi, &chi);                           CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, chi, 9);                   CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->chi, &chi);                       CHKERRQ(ierr);
    }
    // catchment
    if(FSLib->out_catchment) 
    {
        ierr = VecGetArray(FSLib->catchment, &catchment);               CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, catchment, 10);            CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->catchment, &catchment);           CHKERRQ(ierr);
    }
    // lake_depth
    if(FSLib->out_lake_depth) 
    {
        ierr = VecGetArray(FSLib->lake_depth, &lake_depth);             CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, lake_depth, 11);           CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->lake_depth, &lake_depth);         CHKERRQ(ierr);
    }
    // close appended data section and file
    fprintf(fp, "\n\t</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // close file
    fclose(fp);

    ierr = PetscFree(FSLib->buff_fs);    CHKERRQ(ierr);

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

            if( 0 == mode )
            {
                // store node coordinates
                buff[cn++] = (float)(fsX->ncoor[i]); 
                buff[cn++] = (float)(fsY->ncoor[j]);
                buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m
            }
            else if( 2 == mode )
            {
                // store node coordinates
                buff[cn++] = (float)(fsX->ncoor_extend[i]); 
                buff[cn++] = (float)(fsY->ncoor_extend[j]);
                buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m                
            }
            else if(1 == mode || 3 == mode)
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
    PetscErrorCode ierr;
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
    free(fname);

    if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);

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
        ierr = fseek(fp, (*offset), SEEK_SET); CHKERRQ(ierr);
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
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PetscInt mode = 0;

    // check activation
    if(!FSLib->outsurf_fs) PetscFunctionReturn(0);

    // update .pvd file if necessary
    ierr = UpdatePVDFileFS(dirName, FSLib->outfile_fs, "p0.vts", &FSLib->offset_fs, ttime, FSLib->outpvd_fs, step); CHKERRQ(ierr);

    // write sub-domain data .vts files
    if (!FSLib->fs2D)   mode = (FSLib->refine == 1) ? 0 : 1;
    else    mode = (FSLib->refine == 1) ? 2 : 3;

    ierr = PVSurfWriteVTSFS(FSLib, dirName, topo, mode);           CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeSave(FastScapeLib *FSLib, PetscInt step_fs, PetscScalar time_fs)
{
    PetscErrorCode ierr;    
    PetscInt status;    
    char *dirName = NULL;
    PetscScalar *saveArray = PETSC_NULL;        
    Vec saveVec = PETSC_NULL;

    // create directory(encode current time & steo number)    
    // update time stamp and counter
    step_fs++;
    ierr = PetscMalloc1(256, &dirName); CHKERRQ(ierr);
    ierr = PetscSNPrintf(dirName, 256, "Timestep_%1.8lld_%1.8e", (LLD)step_fs, time_fs); CHKERRQ(ierr);

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
        ierr = VecGetArray(saveVec, &saveArray); CHKERRQ(ierr);
        ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, saveArray); CHKERRQ(ierr);
        ierr = VecRestoreArray(saveVec, &saveArray); CHKERRQ(ierr);   
    }

    // free resource
    ierr = PetscFree(dirName); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeFortranCppAdvc(FastScapeLib *FSLib, PetscScalar dt_max, PetscScalar dt_n, PetscInt nstep, 
    PetscScalar step_fs, PetscScalar *vx_pass, PetscScalar *vy_pass, PetscScalar *vz_pass, PetscScalar *topo_pass)
{
    PetscErrorCode ierr;
    PetscInt istep, ind, i, j, idx;
    PetscScalar rangeX, rangeY;
    PetscScalar *topo_random = PETSC_NULL, *kf = PETSC_NULL, *kd = PETSC_NULL;
    PetscScalar *silt_fraction = PETSC_NULL, *basement = PETSC_NULL, *total_erosion = PETSC_NULL;
    PetscScalar *drainage_area = PETSC_NULL, *erosion_rate = PETSC_NULL, *slope = PETSC_NULL;
    PetscScalar *curvature = PETSC_NULL, *chi = PETSC_NULL, *catchment = PETSC_NULL, *lake_depth = PETSC_NULL;
    char *endptr; 
    PetscInt ibc_int;
    std::vector<PetscScalar*> output_arrays;
    std::vector<Vec> output_vecs;
    std::vector<PetscInt> output_flags;
    std::vector<void (*)(PetscScalar*)> output_functions;

    // random noise (uniform distribution)
    mt19937 generator;
    uniform_int_distribution<int> distribution(1, 10000);

    ierr = PetscMalloc(FSLib->nodes_solve * sizeof(PetscScalar), &kf); CHKERRQ(ierr);
    ierr = PetscMalloc(FSLib->nodes_solve * sizeof(PetscScalar), &kd); CHKERRQ(ierr);
 
    // output 
    auto manage_output = [&](int flag, Vec vec, PetscScalar** array, void (*copy_func)(PetscScalar*)) -> PetscErrorCode 
    {
        if (flag) 
        {
            ierr = VecGetArray(vec, array); CHKERRQ(ierr);
            output_flags.push_back(flag);            
            output_arrays.push_back(*array);
            output_vecs.push_back(vec);
            output_functions.push_back(copy_func);
        }
        return 0;
    };
    manage_output(FSLib->out_silt_fraction, FSLib->silt_fraction, &silt_fraction, fastscape_copy_f_);
    manage_output(FSLib->out_basement, FSLib->basement, &basement, fastscape_copy_basement_);
    manage_output(FSLib->out_total_erosion, FSLib->total_erosion, &total_erosion, fastscape_copy_total_erosion_);
    manage_output(FSLib->out_drainage_area, FSLib->drainage_area, &drainage_area, fastscape_copy_drainage_area_);
    manage_output(FSLib->out_erosion_rate, FSLib->erosion_rate, &erosion_rate, fastscape_copy_erosion_rate_);
    manage_output(FSLib->out_slope, FSLib->slope, &slope, fastscape_copy_slope_);
    manage_output(FSLib->out_curvature, FSLib->curvature, &curvature, fastscape_copy_curvature_);  
    manage_output(FSLib->out_chi, FSLib->chi, &chi, fastscape_copy_chi_);
    manage_output(FSLib->out_catchment, FSLib->catchment, &catchment, fastscape_copy_catchment_);
    manage_output(FSLib->out_lake_depth, FSLib->lake_depth, &lake_depth, fastscape_copy_lake_depth_);

    // initialize FastScape
    fastscape_init_();
    fastscape_set_nx_ny_(&FSLib->nx_solve, &FSLib->ny_solve);  

    // allocate memory
    fastscape_setup_();

    // set model dimensions
    rangeX = (0 == FSLib->fs2D) ? FSLib->rangeX : FSLib->extendedXRange;     
    rangeY = (0 == FSLib->fs2D) ? FSLib->rangeY : FSLib->extendedYRange;    
    fastscape_set_xl_yl_(&rangeX, &rangeY);

    // set time step
    fastscape_set_dt_(&dt_max);

    // set random initial topography & erosional parameters
    if(0 == step_fs) 
    {
        // random noise
        if (FSLib->random_noise) 
        {
            ierr = PetscMalloc(FSLib->nodes_solve * sizeof(PetscScalar), &topo_random); CHKERRQ(ierr);
            for (ind = 0; ind < FSLib->nodes_solve; ind++) 
            {
                topo_random[ind] = distribution(generator)/10000.0;
                topo_pass[ind]   += topo_random[ind];
            }
        }
        // kf & kd
        std::fill_n(kf, FSLib->nodes_solve, FSLib->kf);
        std::fill_n(kd, FSLib->nodes_solve, FSLib->kd);
    }
    else
    {
        // kf & kd
        std::fill_n(kf, FSLib->nodes_solve, FSLib->kf);
        std::fill_n(kd, FSLib->nodes_solve, FSLib->kd);        
    }

    // velocity boundary conditions
    auto set_velocity_boundary = [&]() 
    {        
        std::vector<std::function<bool()>> bc_conditions = 
        {
            [&](){ return 0 == j; },
            [&](){ return FSLib->nx_solve - 1 == i; },
            [&](){ return FSLib->ny_solve - 1 == j; },
            [&](){ return 0 == i; }
        };

        for (idx = 0; idx < 4; idx++)
        {
            if ('1' == FSLib->FS_VELBC[idx] && bc_conditions[idx]()) 
            {
                std::fill_n(vx_pass, FSLib->nodes_solve, 0);
                std::fill_n(vy_pass, FSLib->nodes_solve, 0);
                std::fill_n(vz_pass, FSLib->nodes_solve, 0);
            }    
        }
    };
    set_velocity_boundary();

    // set velocity
    fastscape_set_u_(vz_pass);
    fastscape_set_v_(vx_pass, vy_pass);

    // set topography boundary conditions
    ibc_int = strtol(FSLib->FS_BC, &endptr, 10);
    fastscape_set_bc_(&ibc_int);

    // initialize topography
    fastscape_init_h_(topo_pass);

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
        if (nstep - 1 == istep && dt_n > 0) 
        {
            fastscape_set_dt_(&dt_n);
        }

        // execute step
        fastscape_execute_step_();
        // get value of time step counter
        fastscape_get_step_(&istep);

        if(nstep == istep)
        { 
            // output
            // topography
            fastscape_copy_h_(topo_pass);
            // others
            for (i = 0; i < output_arrays.size(); i++) 
            {
                if (output_flags[i]) output_functions[i](output_arrays[i]);
            }           
        }      
    } while (istep < nstep);

    // output timing
    fastscape_debug_();

    // end FastScape run
    fastscape_destroy_();

    // clear
    for (i = 0; i < output_arrays.size(); i++) 
    {
        if(output_flags[i])  
        {
            ierr = VecRestoreArray(output_vecs[i], &output_arrays[i]); CHKERRQ(ierr);
        }
    }

    ierr = PetscFree(kf); CHKERRQ(ierr);
    ierr = PetscFree(kd); CHKERRQ(ierr);

    if (FSLib->random_noise) 
    {
        ierr = PetscFree(topo_random); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeReadRestart(FastScapeLib *FSLib, FILE *fp)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        ierr =  FastScapeCreateSurfaceGrid(FSLib, 2); CHKERRQ(ierr); 

        if( 0 == FSLib->fs2D && 1 == FSLib->refine  )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_fs, fp); CHKERRQ(ierr);
        }
        if( 1 == FSLib->non_uniform_grid)
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_nug, fp); CHKERRQ(ierr);            
        }
        if( 0 == FSLib->fs2D && FSLib->refine > 1  )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_refine, fp); CHKERRQ(ierr);
        }
        else if( 1 == FSLib->fs2D && 1 == FSLib->refine )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_extend, fp); CHKERRQ(ierr);        
        }
        else if( 1 == FSLib->fs2D && FSLib->refine > 1)
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_et_refine, fp); CHKERRQ(ierr);          
        }
    }   

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeWriteRestart(FastScapeLib *FSLib, FILE *fp)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    if(ISRankZero(PETSC_COMM_WORLD))
    {

        if( 0 == FSLib->fs2D && 1 == FSLib->refine  )
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_fs, fp); CHKERRQ(ierr);
        }
        if(1 == FSLib->non_uniform_grid)
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_nug, fp); CHKERRQ(ierr);            
        }
        if( 0 == FSLib->fs2D && FSLib->refine > 1)
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_refine, fp); CHKERRQ(ierr);
        }
        else if( (1 == FSLib->fs2D && 1 == FSLib->refine)  )
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_extend, fp); CHKERRQ(ierr);        
        }
        else if( 1 == FSLib->fs2D && FSLib->refine > 1)
        {
            // store topography vector 
            ierr = VecWriteRestart(FSLib->gtopo_et_refine, fp); CHKERRQ(ierr);         
        }
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeDestroy(FastScapeLib *FSLib)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    ierr = VecDestroy(&FSLib->vx_LaMEM); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_LaMEM); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vz_LaMEM); CHKERRQ(ierr);  
    ierr = VecDestroy(&FSLib->gtopo_fs); CHKERRQ(ierr);  
    ierr = VecDestroy(&FSLib->vx_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vz_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->gtopo_collect); CHKERRQ(ierr);      

    if(FSLib->non_uniform_grid)
    {
        ierr = VecDestroy(&FSLib->gtopo_nug); CHKERRQ(ierr);         
        ierr = VecDestroy(&FSLib->vx_nug); CHKERRQ(ierr);   
        ierr = VecDestroy(&FSLib->vy_nug); CHKERRQ(ierr);   
        ierr = VecDestroy(&FSLib->vz_nug); CHKERRQ(ierr);   
    }

    if(FSLib->fs2D)
    {
        // extended part
        ierr = VecDestroy(&FSLib->gtopo_extend); CHKERRQ(ierr);  
        ierr = VecDestroy(&FSLib->vx_extend); CHKERRQ(ierr); 
        ierr = VecDestroy(&FSLib->vy_extend); CHKERRQ(ierr); 
        ierr = VecDestroy(&FSLib->vz_extend); CHKERRQ(ierr); 

        if(1 < FSLib->refine)
        {
            // extended part after refinement
            ierr = VecDestroy(&FSLib->gtopo_refine); CHKERRQ(ierr);   
            ierr = VecDestroy(&FSLib->vx_refine); CHKERRQ(ierr);   
            ierr = VecDestroy(&FSLib->vy_refine); CHKERRQ(ierr);  
            ierr = VecDestroy(&FSLib->vz_refine); CHKERRQ(ierr);  
        }
    }
    else
    {
        if(1 < FSLib->refine)
        {
            // refined part
            ierr = VecDestroy(&FSLib->gtopo_et_refine); CHKERRQ(ierr); 
            ierr = VecDestroy(&FSLib->vx_et_refine); CHKERRQ(ierr);
            ierr = VecDestroy(&FSLib->vy_et_refine); CHKERRQ(ierr);
            ierr = VecDestroy(&FSLib->vz_et_refine); CHKERRQ(ierr);
        }
    }

    // FastScape solution
    ierr = VecDestroy(&FSLib->silt_fraction); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->basement); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->total_erosion); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->drainage_area); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->erosion_rate); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->slope); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->curvature); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->chi); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->catchment); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->lake_depth); CHKERRQ(ierr);

    PetscFree(FSLib->buff_fs);

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        free(FSLib->fsX.ncoor);
        free(FSLib->fsY.ncoor);
        FSLib->fsX.ncoor = PETSC_NULL;
        FSLib->fsY.ncoor = PETSC_NULL;

        if(FSLib->non_uniform_grid)
        {
            if(FSLib->extendedX)
            {
                free(FSLib->ncoor_ori_y);
                FSLib->ncoor_ori_y = PETSC_NULL;
            }
            else
            {
                free(FSLib->ncoor_ori_x);
                FSLib->ncoor_ori_x = PETSC_NULL;
            }       
        }

        if(FSLib->fs2D)
        {
            free(FSLib->fsX.ncoor_extend);
            free(FSLib->fsY.ncoor_extend);     
            FSLib->fsX.ncoor_extend = PETSC_NULL;
            FSLib->fsY.ncoor_extend = PETSC_NULL;
        }

        if(1 < FSLib->refine)
        {
            free(FSLib->fsX.ncoor_refine);
            free(FSLib->fsY.ncoor_refine);
            FSLib->fsX.ncoor_refine = PETSC_NULL;
            FSLib->fsY.ncoor_refine = PETSC_NULL;
        }
    }

    PetscFunctionReturn(0);
}
