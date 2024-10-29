#include <iostream>
#include <fstream>

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

PetscErrorCode fastscape(FreeSurf *surf, AdvCtx *actx)
{
    // Apply erosion to the internal free surface of the model
    PetscErrorCode ierr;

    // free surface cases only
    if(!surf->UseFreeSurf) PetscFunctionReturn(0);

    // load global nx, ny, dt, time, rangeX, rangeY
    PetscScalar dt, time, rangeX_begin, rangeY_begin, rangeX_end, rangeY_end, rangeZ_begin, rangeZ_end;
    PetscScalar rangeX, rangeY, bx, by, bz, ex, ey, ez, chLen, level, cf;
    PetscInt nx_fs, ny_fs, ind, ind2;

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

    rangeX = (rangeX_end - rangeX_begin)*1e3;
    rangeY = (rangeY_end - rangeY_begin)*1e3;

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

        double vz_pass[tnodes], topo_pass[tnodes];

        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
                vz_pass[ind]=vz_fs[ind]; // (km)
                topo_pass[ind]=topo_fs[ind]; // (km)
    //            printf("topo_fs[0][%d][%d]: %f      \n",j,i,topo_fs[j*ny_fs+i]); // (km)
    //            printf("vz_fs[0][%d][%d]: %f      \n",j,i,vz_fs[j*ny_fs+i]); // (km)
            } 
        }        
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

        double dt_max = 10; // 最大步长,如果LaMEM计算出来时间步长比这个大，就用这个
        double dt_n=0.0; 
        int nsteps = dt/dt_max;
        // remain supplement (about reminder)
        if(nsteps < 1) {
            nsteps = 1;
            dt_max = dt;
        }
        else
        {
            nsteps = nsteps + 1;
            dt_n=dt-(nsteps-1)*dt_max;
        }
        printf("nsteps:%d;  dt: %f, dt_n: %f\n",nsteps,dt_max,dt_n);
        
        printf("\n ============= receive from LaMEM done ====================\n");

        double* topo_pass_f = NULL;

        topo_pass_f = fastscapeFortran(&nx_fs,&ny_fs,&rangeX,&rangeY,&dt_max,&dt_n,&nsteps,vz_pass,topo_pass);
  /*      
        for(j = 1; j < ny_fs+1; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
         //       init_elevation[ind] = init_elevation[ind] + topo_fs[ind]*1e3; // km(LaMEM) to m(FastScape)
                printf("elevation[0][%d][%d]: %f      \n",j,i,topo_pass_f[ind]); // (km)
            }
        }
*/
        clearArray();

        // Correction
        for(j = 0; j < ny_fs; j++) 
        {   
            for(i = 0; i < nx_fs; i++) 
            {
                ind = j*ny_fs+i;
                ind2 = (j+1)*ny_fs+i;

                topo_fs[ind] = topo_pass_f[ind2]/1e3; // m(FastScape) to km(LaMEM) 

                if(topo_fs[ind] > rangeZ_end) topo_fs[ind] = rangeZ_end;
                if(topo_fs[ind] < rangeZ_begin) topo_fs[ind] = rangeZ_begin;

         //       printf("topo_new[%d][%d][%d]: %f      ",L,j,i,topo[L][j][i]); // (km)
         //       printf("init_elevation[%d]:%f      ",ind,init_elevation[ind]);
            }
        }
        topo_pass_f = NULL;
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





