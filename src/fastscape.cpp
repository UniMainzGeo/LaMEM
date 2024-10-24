#include <iostream>
#include <fstream>

#include <stdio.h>
#include <string.h>

#include "xtensor.hpp"
#include "xtensor/xtensor.hpp" 
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/eroders/diffusion_adi.hpp"
#include "fastscapelib/eroders/spl.hpp"

// LaMEM header file
#include "LaMEM.h"
#include "surf.h"
#include "scaling.h"
#include "JacRes.h"
#include "tssolve.h"
#include "advect.h"
#include "interpolate.h"
#include "fastscape.h"

namespace fs = fastscapelib;

PetscErrorCode fastscape(FreeSurf *surf, AdvCtx *actx)
{
    // Apply erosion to the internal free surface of the model
    PetscErrorCode ierr;

    // free surface cases only
    if(!surf->UseFreeSurf) PetscFunctionReturn(0);

    // load global nx, ny, dt, time, rangeX, rangeY
    PetscScalar dt, time, rangeX_begin, rangeY_begin, rangeX_end, rangeY_end, rangeZ_begin, rangeZ_end;
    PetscScalar rangeX, rangeY, bx, by, bz, ex, ey, ez, chLen, level, cf;
    PetscInt nx_fs, ny_fs, ind;

    // nx, ny, dt, time
    JacRes      *jr;
    FDSTAG      *fs;

    PetscPrintf(PETSC_COMM_WORLD,"\n============= receive Vz, Topography, grid, model_range from LaMEM ==============\n");
    
    jr   = surf->jr;
    fs = jr->fs;

    nx_fs = fs->dsx.tnods;
    ny_fs = fs->dsy.tnods;

    dt   = jr->ts->dt;
    time = jr->ts->time;
    PetscPrintf(PETSC_COMM_WORLD,"\n nx: %d, ny: %d \n", nx_fs, ny_fs);
    PetscPrintf(PETSC_COMM_WORLD,"\n dt_LaMEM %f, time: %f \n", dt, time);
    
    chLen = fs->scal->length;
    ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

    rangeX_begin = bx*chLen;
    rangeY_begin = by*chLen;
    rangeZ_begin = bz*chLen;   
    rangeX_end = ex*chLen;
    rangeY_end = ey*chLen;
    rangeZ_end = ez*chLen;

    rangeX = rangeX_end - rangeX_begin;
    rangeY = rangeY_end - rangeY_begin;

    PetscPrintf(PETSC_COMM_WORLD, "x_begin: %f, x_end: %f, y_begin: %f, y_end: %f; \n",
            rangeX_begin, rangeX_end, rangeY_begin, rangeY_end);

    // Gather topography and velocity
    PetscInt    i, j, tnodes;
    PetscScalar *topo_fs=PETSC_NULL;
    PetscScalar *vz_fs=PETSC_NULL;

    tnodes = nx_fs*ny_fs;

    topo_fs = (PetscScalar *)malloc(tnodes*sizeof(PetscScalar)); 
    vz_fs = (PetscScalar *)malloc(tnodes*sizeof(PetscScalar)); 

    // topography
    ierr = VecScatterCreateToZero(surf->gtopo, &surf->ctx, &surf->gtopo_fs); CHKERRQ(ierr); 
    ierr = VecScatterBegin(surf->ctx, surf->gtopo, surf->gtopo_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(surf->ctx, surf->gtopo, surf->gtopo_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    //velocity
    ierr = VecScatterCreateToZero(surf->vz, &surf->ctx, &surf->vz_fs); CHKERRQ(ierr); 
    ierr = VecScatterBegin(surf->ctx, surf->vz, surf->vz_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(surf->ctx, surf->vz, surf->vz_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);


    // Computing topography in rank 0
    if(ISRankZero(PETSC_COMM_WORLD))
    {
        ierr = VecGetArray(surf->gtopo_fs,  &topo_fs);  CHKERRQ(ierr);
        ierr = VecGetArray(surf->vz_fs,  &vz_fs);  CHKERRQ(ierr);
/*
        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                printf("topo_fs[0][%d][%d]: %f      ",j,i,topo_fs[j*ny_fs+i]); // (km)
                printf("vz_fs[0][%d][%d]: %f      ",j,i,vz_fs[j*ny_fs+i]); // (km)
            } 
        }
        */
        
        // high level that begin erosion // wating for finish
        level = 1;
    /* 
        jj = surf->numErPhs-1
        level = surf->erLevels[jj];
    */
        // Setting sedimental phase // wating for finish
     //   jj = surf->numErPhs-1;
        surf->phase = 3;
        printf("\n ============= receive from LaMEM done ====================\n");

        // raster grid and boundary conditions
        fs::raster_boundary_status bs(fs::node_status::fixed_value);

        auto grid = fs::raster_grid<>::from_length({ nx_fs, ny_fs }, { rangeX, rangeY }, bs);

        // flow graph with single direction flow routing
        fs::flow_graph<fs::raster_grid<>> flow_graph(
            grid, { fs::single_flow_router(), fs::mst_sink_resolver() });

        // Setup eroders
        auto spl_eroder = fs::make_spl_eroder(flow_graph, 1e-4, 0.4, 1, 1e-3);
        auto diffusion_eroder = fs::diffusion_adi_eroder(grid, 0.01);

        // initial topographic surface elevation (flat surface + random perturbations)
        xt::xarray<double> init_elevation = xt::random::rand<double>(grid.shape());

        // A loop judgment, only the first step the first loop terrain plus random noise, nothing else
        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
                init_elevation[ind] = init_elevation[ind] + topo_fs[ind]*1e3; // km(LaMEM) to m(FastScape)
         //       printf("init_elevation[%d]:%f      ",ind,init_elevation[ind]);
            }
        }

        xt::xarray<double> elevation = init_elevation;

        // init drainage area and temp arrays
        xt::xarray<double> drainage_area(elevation.shape());
        xt::xarray<double> uplifted_elevation(elevation.shape());

        // uplift rate
        xt::xarray<double> uplift_rate(elevation.shape(), 1e-3);  // uplift rate (m/yr)
        
        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
                if(!(PetscInt)fs->dsz.rank)
                {
                    uplift_rate[ind] = vz_fs[ind]*cf*1e2; // cm/yr(LaMEM) to m/yr(FastScape)
                } 
                else
                {
                    uplift_rate[ind] = vz_fs[ind]*1e2; // cm/yr(LaMEM) to m/yr(FastScape)
                }
    //            printf("uplift_rate[%d]: %f      ",ind,uplift_rate[ind]); // cm/yr(LaMEM)
            }
        }

        auto row_bounds = xt::view(uplift_rate, xt::keep(0, -1), xt::all());
        row_bounds = 0.0;
        auto col_bounds = xt::view(uplift_rate, xt::all(), xt::keep(0, -1));
        col_bounds = 0.0;

        //
        // run model
        //
        double dt_max = 10; // Maximum step length, if LaMEM calculates a time step larger than this, use this
        double dt_n; 
        int nsteps = dt/dt_max;
        // remain supplement (about reminder)
        if(nsteps < 1) {
            nsteps = 1;
            dt_max = dt;
        }
     //   if(isnan(dt)) dt = 1; // vy,vx do not converge, dt is assigned 1
        printf("nsteps:%d;  dt: %f\n",nsteps,dt_max);

        if(nsteps > 1) dt_n=dt-nsteps*dt_max;

        int error_fs =0 ;
        // get size of box

        for (int step = 0; step < nsteps; step++)
        {
            // apply uplift
            uplifted_elevation = elevation + dt_max * uplift_rate;

            // flow routing
            flow_graph.update_routes(uplifted_elevation);

            // flow accumulation (drainage area)
            flow_graph.accumulate(drainage_area, 1.0);

            // apply channel erosion then hillslope diffusion
            auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, dt_max);
            auto diff_erosion = diffusion_eroder.erode(uplifted_elevation - spl_erosion, dt_max);

            // update topography
            elevation = uplifted_elevation - spl_erosion - diff_erosion;

            error_fs ++;
            if(error_fs - step > 100){
                printf("FastScape Error");
                break;
            }
        }

    //    std::cout<<"mean final elevation:"<< xt::mean(elevation)<<std::endl;
        // Correction
        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
            
                if(elevation[ind] > rangeZ_end) elevation[ind] = rangeZ_end;
                if(elevation[ind] < rangeZ_begin) elevation[ind] = rangeZ_begin;

                topo_fs[ind] = elevation[ind]/1e3; // m(FastScape) to km(LaMEM) 
         //       printf("topo_new[%d][%d][%d]: %f      ",L,j,i,topo[L][j][i]); // (km)
         //       printf("init_elevation[%d]:%f      ",ind,init_elevation[ind]);
            }
        }
        

    }

    // Broadcast      
    if(ISParallel(PETSC_COMM_WORLD))
    {
        ierr = MPI_Bcast(topo_fs, (PetscMPIInt)tnodes, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    // Save topography in different ranks
    PetscInt    L, sx, sy, sz, nx, ny, nz;
    PetscScalar ***topo;

    FDSTAG      *fs_actx;

    fs_actx = actx->fs;

    sx = fs_actx->dsx.pstart; nx = fs_actx->dsx.ncels;
    sy = fs_actx->dsy.pstart; ny = fs_actx->dsy.ncels;
    sz = fs_actx->dsz.pstart; nz = fs_actx->dsz.ncels;

    L    = (PetscInt)fs->dsz.rank;
//    printf("\n L: %d \n", L);
    
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);


    START_PLANE_LOOP
    {
//        printf("topo_old[%d][%d][%d] : %f      ", L, j, i, topo[L][j][i]);
        topo[L][j][i] = topo_fs[j*ny_fs+i]; // (km)
//        printf("topo_new[%d][%d][%d] : %f      ", L, j, i, topo[L][j][i]);
    }
    END_PLANE_LOOP  

    ierr = VecRestoreArray(surf->gtopo_fs,  &topo_fs);  CHKERRQ(ierr);
    
    ierr = VecScatterDestroy(&surf->ctx); CHKERRQ(ierr);
    ierr = VecDestroy(&surf->gtopo_fs);   CHKERRQ(ierr);
    ierr = VecDestroy(&surf->vz_fs);   CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"\n============= FastScape Done =================\n");

    PetscFunctionReturn(0);
}





