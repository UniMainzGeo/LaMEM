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
#include "paraViewOutSurf.h"
#include "paraViewOutBin.h"

PetscErrorCode FastScapeCreate(FastScapeLib *FSLib, FB *fb)
{
    PetscErrorCode ierr;
    PetscInt maxPhaseID;

    // initialize
	FSLib->refine	  =  1;

	maxPhaseID = FSLib->surf->jr->dbm->numPhases-1;

    // kf, kd, kfsed, kdsed, can be set as an array
	ierr = getScalarParam(fb, _REQUIRED_, "Max_dt",          &FSLib->Max_dt,     1,                 1.0);      				CHKERRQ(ierr); // m/yr
 //	ierr = getScalarParam(fb, _REQUIRED_, "SPL_kf",          &FSLib->kf,   	     1,                 1.0);      				CHKERRQ(ierr); // m/yr
	ierr = getScalarParam(fb, _REQUIRED_, "SPL_kfsed",       &FSLib->kfsed,      1,                 1.0);      				CHKERRQ(ierr); // m/yr
	ierr = getScalarParam(fb, _REQUIRED_, "SPL_m",           &FSLib->m,    	     1,                 1.0);      				CHKERRQ(ierr); // non-dimensional
	ierr = getScalarParam(fb, _REQUIRED_, "SPL_n",           &FSLib->n,          1,                 1.0);     				CHKERRQ(ierr); // non-dimensional
//	ierr = getScalarParam(fb, _REQUIRED_, "hillslope_kd",    &FSLib->kd,         1,                 1.0);      				CHKERRQ(ierr); // m/yr
	ierr = getScalarParam(fb, _REQUIRED_, "hillslope_kdsed", &FSLib->kdsed,      1,                 1.0);      				CHKERRQ(ierr); // m/yr
	ierr = getScalarParam(fb, _REQUIRED_, "hillslope_g",     &FSLib->g,          1,                 1.0);   				CHKERRQ(ierr); // non-dimensional
	ierr = getScalarParam(fb, _REQUIRED_, "hillslope_gsed",  &FSLib->gsed,       1,                 1.0);     			 	CHKERRQ(ierr); // non-dimensional
	ierr = getScalarParam(fb, _REQUIRED_, "multiFlow_p",     &FSLib->p,     	 1,                 1.0);    				CHKERRQ(ierr); // non-dimensional
	ierr = getIntParam   (fb, _REQUIRED_, "fs_refine",       &FSLib->refine,     1,                 100);  				 	CHKERRQ(ierr); // non-dimensional
	ierr = getIntParam   (fb, _REQUIRED_, "sed_phases",      &FSLib->sedPhases,  1,   		 maxPhaseID);                   CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, "fs_bc",           &FSLib->FS_BC,      1,   		       1111);                   CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, "fs_bgphase",      &FSLib->bgphase,    1,   		 maxPhaseID);                   CHKERRQ(ierr);
    ierr = getIntParam   (fb, _OPTIONAL_, "setMarine",       &FSLib->setMarine, 1,                   1);                    CHKERRQ(ierr);

    if(FSLib->setMarine == 1)
    {
        ierr = getScalarParam(fb, _REQUIRED_, "marine_sealevel",  &FSLib->sealevel,   1,       1.0);                      CHKERRQ(ierr); // m
        ierr = getScalarParam(fb, _REQUIRED_, "marine_poroSilt",  &FSLib->poro_silt,  1,       1.0);                      CHKERRQ(ierr); 
        ierr = getScalarParam(fb, _REQUIRED_, "marine_poroSand",  &FSLib->poro_sand,  1,       1.0);                      CHKERRQ(ierr); 
        ierr = getScalarParam(fb, _REQUIRED_, "marine_zporoSilt", &FSLib->zporo_silt, 1,       1.0);                      CHKERRQ(ierr); 
        ierr = getScalarParam(fb, _REQUIRED_, "marine_zporoSand", &FSLib->zporo_sand, 1,       1.0);                      CHKERRQ(ierr); 
        ierr = getScalarParam(fb, _REQUIRED_, "marine_ratio",     &FSLib->ratio,      1,       1.0);                      CHKERRQ(ierr); 
        ierr = getScalarParam(fb, _REQUIRED_, "marine_L",         &FSLib->Lsolve,     1,       1.0);                      CHKERRQ(ierr);  // m
        ierr = getScalarParam(fb, _REQUIRED_, "marine_kdsSilt",   &FSLib->kds_silt,   1,       1.0);                      CHKERRQ(ierr);  // m2/yr
        ierr = getScalarParam(fb, _REQUIRED_, "marine_kdsSand",   &FSLib->kds_sand,   1,       1.0);                      CHKERRQ(ierr);  // m2/yr
    }

	PetscPrintf(PETSC_COMM_WORLD, "FastScape parameters: \n");	
	PetscPrintf(PETSC_COMM_WORLD, "   Max timestep         : %g\n", FSLib->Max_dt);			
	PetscPrintf(PETSC_COMM_WORLD, "   SPL: \n");	
	PetscPrintf(PETSC_COMM_WORLD, "      Kf                : %g\n", FSLib->kf);
    PetscPrintf(PETSC_COMM_WORLD, "      Kfsed             : %g\n", FSLib->kfsed);
	PetscPrintf(PETSC_COMM_WORLD, "      m                 : %g\n", FSLib->m);	
	PetscPrintf(PETSC_COMM_WORLD, "      n                 : %g\n", FSLib->n);
	PetscPrintf(PETSC_COMM_WORLD, "   Hillslope process: \n");	
	PetscPrintf(PETSC_COMM_WORLD, "      Kd                : %g\n", FSLib->kd);
	PetscPrintf(PETSC_COMM_WORLD, "      Kdsed             : %g\n", FSLib->kdsed);
	PetscPrintf(PETSC_COMM_WORLD, "      g                 : %g\n", FSLib->g);	
	PetscPrintf(PETSC_COMM_WORLD, "      gsed              : %g\n", FSLib->gsed);
	PetscPrintf(PETSC_COMM_WORLD, "      p                 : %g\n", FSLib->p);	
	PetscPrintf(PETSC_COMM_WORLD, "   Sedimentation phase  : %d\n", FSLib->sedPhases);	
	PetscPrintf(PETSC_COMM_WORLD, "   Background phase     : %d\n", FSLib->bgphase);	 
	if(FSLib->refine > 1) PetscPrintf(PETSC_COMM_WORLD, "   refined times        : %d\n", FSLib->refine);	
    if(FSLib->setMarine == 1)
    {
        PetscPrintf(PETSC_COMM_WORLD, "   Marine process:\n"); 
        PetscPrintf(PETSC_COMM_WORLD, "      sealevel          : %g\n", FSLib->sealevel);
        PetscPrintf(PETSC_COMM_WORLD, "      poro_silt         : %g\n", FSLib->poro_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      poro_sand         : %g\n", FSLib->poro_sand);
        PetscPrintf(PETSC_COMM_WORLD, "      zporo_silt        : %g\n", FSLib->zporo_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      zporo_sand        : %g\n", FSLib->zporo_sand);
        PetscPrintf(PETSC_COMM_WORLD, "      ratio             : %g\n", FSLib->ratio);
        PetscPrintf(PETSC_COMM_WORLD, "      L                 : %g\n", FSLib->Lsolve);
        PetscPrintf(PETSC_COMM_WORLD, "      kds_silt          : %g\n", FSLib->kds_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      kds_sand          : %g\n", FSLib->kds_sand);        
    }
	PetscPrintf(PETSC_COMM_WORLD, "   Boundary condition   : %d\n", FSLib->FS_BC);	
	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

    PetscFunctionReturn(0);
}

PetscErrorCode savePvtsFS(PVSurf *pvsurf, FastScapeLib *FSLib, PetscScalar ttime, PetscInt step, const char *dirName, PetscScalar *topo)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // check activation
    if(!pvsurf->outsurf) PetscFunctionReturn(0);

    // update .pvd file if necessary
    ierr = UpdatePVDFileRefine(dirName, pvsurf->outfile_refine, "p0.vts", &pvsurf->offset_refine, ttime, pvsurf->outpvd, step); CHKERRQ(ierr);

    // write sub-domain data .vts files
   ierr = PVSurfWriteVTSRefine(pvsurf, dirName, FSLib->nx_refine, FSLib->ny_refine, FSLib->rangeX, FSLib->rangeY, topo); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} 

PetscErrorCode bilinearInterpolate(FastScapeLib *FSLib, PetscScalar data[], PetscScalar data_refine[], Scaling *scal, PetscInt corMode)
{
    //PetscInt interpolationMode = 1;
    // linear interpolation
    // refine = 1; no nodes; 2, add a node between two original nearest nodes; 3, add two nodes between two nearest original nodes;
    // corMode: 1 -- topography, (km) in LaMEM to (m) in FastScape; 2 -- velocity, (cm/yr) in LaMEM to (m/yr) in FastScape

    PetscInt i, j, ind, ind_a, ind_b, ind_aa, ind_bb, countX1, countX2, countY, tnodes_refine;
    PetscScalar distance_ax, distance_ay;

    tnodes_refine = FSLib->nx_refine * FSLib->ny_refine;
    countX1 = 0;
    countX2 = 0;
    countY = 0;
    ind_a = countX1 * FSLib->refine;
    ind_b = (countX1 + 1) * FSLib->refine;
    ind_aa = countX2;
    ind_bb = countX2 + 1;

    // interpolate in x-direction
    for(j = 0; j < FSLib->ny_refine; j += FSLib->refine) 
    {   
        for(i = 0; i < FSLib->nx_refine; i++) 
        {
            ind = j * FSLib->nx_refine+i;

            if(0 == ind%FSLib->nx_refine) 
            {
                countX1 = 0;   
                if(0 == ind)
                {
                    if(1 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * 1e3;
                        data_refine[ind_b] = data[ind_bb] * 1e3;
                    }
                    else if(2 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * scal->velocity/1e2;
                        data_refine[ind_b] = data[ind_bb] * scal->velocity/1e2; 
                    }
                    else
                    {
                        data_refine[ind_a] = data[ind_aa];
                        data_refine[ind_b] = data[ind_bb];                          
                    }
                }
                else
                {
                    ind_a = countX1 * FSLib->refine + j * FSLib->nx_refine;
                    ind_b = (countX1+1) * FSLib->refine + j * FSLib->nx_refine; 
                    countX2++; 
                    ind_aa = countX2;
                    ind_bb = countX2 + 1;
 
                    if(1 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * 1e3;
                        data_refine[ind_b] = data[ind_bb] * 1e3;                       
                    }
                    else if(2 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * scal->velocity/1e2;
                        data_refine[ind_b] = data[ind_bb] * scal->velocity/1e2; 
                    }
                    else            
                    {
                        data_refine[ind_a] = data[ind_aa];
                        data_refine[ind_b] = data[ind_bb];                        
                    }        

                }
            }
            else
            {
                if(ind == ind_b) 
                {
                    if(ind == tnodes_refine-1) 
                    {
                        if(1 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * 1e3;                           
                        }
                        else if(2 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * scal->velocity/1e2; 
                        }
                        else 
                        {
                            data_refine[ind_b] = data[ind_bb];
                        }
                        

                        goto skip;
                    }
                    else
                    {
                        countX1++;
                        countX2++;
                        ind_a = countX1 * FSLib->refine + j * FSLib->nx_refine;
                        ind_b = (countX1+1) * FSLib->refine + j * FSLib->nx_refine;
                        ind_aa = countX2;
                        ind_bb = countX2 + 1;

                        if(1 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * 1e3;                       
                        }
                        else if(2 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * scal->velocity / 1e2; 
                        }
                        else  
                        {
                            data_refine[ind_b] = data[ind_bb];
                        }
                    }
                }
                else
                {
                    distance_ax = ((ind-ind_a)%FSLib->nx_refine) * 1.0 / FSLib->refine;   
                    data_refine[ind] = data_refine[ind_a] * (1 - distance_ax)  + data_refine[ind_b] * distance_ax;               
                }
            }                   
        }
    }
    skip: // printf("interpolation in x-direction done\n");     

    // interpolate in y-direction
    for(j = 0; j < FSLib->ny_refine; j++) 
    {   
        if(0 == j%FSLib->refine) 
        {
            if( 0 != j) countY += FSLib->refine;    
            continue;
        }

        for(i = 0; i < FSLib->nx_refine; i++) 
        {

            ind = j * FSLib->nx_refine + i;
            ind_a = countY * FSLib->nx_refine + i;
            ind_b = (countY+FSLib->refine) * FSLib->nx_refine + i;

            distance_ay = (j - countY)*1.0/FSLib->refine;  
            data_refine[ind] = data_refine[ind_a] * (1 - distance_ay)  + data_refine[ind_b] * distance_ay;                   
        }         
    }   
  //  printf("interpolation in y-direction done\n");

    PetscFunctionReturn(0);
}

PetscErrorCode fastscape(FastScapeLib *FSLib)
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
    PetscScalar dt, dt_scal, dt_fs, rangeX_begin, rangeY_begin, rangeX_end, rangeY_end, rangeZ_begin, rangeZ_end;
    PetscScalar bx, by, bz, ex, ey, ez, chLen, time_fs, topo_f, vz_f;
    PetscInt nx_fs, ny_fs, ind, ind2, i, j, k, tnodes, L, sx, sy, sz, nx, ny, nz, step_fs, tproc, cproc, rank_id;
    PetscScalar *topo_fs=PETSC_NULL;
    PetscScalar *vz_fs=PETSC_NULL;
    PetscScalar ***topo;
    PetscScalar ***vz;
    PetscScalar ***vz_collect;
    char           *dirName;

    // load global nx, ny, dt, time, rangeX, rangeY
    // nx, ny, dt, time
    JacRes      *jr;
    FDSTAG      *fs;
    TSSol       *ts;
    Scaling     *scal;
	PVSurf *pvsurf;
    DBMat  *dbm;
    Material_t *m;

   // PetscPrintf(PETSC_COMM_WORLD,"\n------------- receive information from LaMEM -------------\n");
    
    jr   = surf->jr;
    fs = jr->fs;
    ts = jr->ts;
    scal     = ts->scal;
    pvsurf = FSLib->pvsurf;
    dbm    = FSLib->dbm;
    m      = dbm->phases;

    nx_fs = fs->dsx.tnods;
    ny_fs = fs->dsy.tnods;

    dt   = ts->dt;
    dt_scal = dt * scal->time;
    dt_fs = dt_scal * 1e6; // (Myr) in LaMEM to (yr) in FastScape
    time_fs = ts->time * scal->time + dt_scal; // time after finishing surface processes
    step_fs = ts->istep;

 //   PetscPrintf(PETSC_COMM_WORLD,"\nnx_fs: %d, ny_fs: %d \n", nx_fs, ny_fs);
 //   PetscPrintf(PETSC_COMM_WORLD," dt_LaMEM: %f, dt_scal: %f (Myr), dt_fs: %f (yr) \n", dt, dt_scal, dt_fs); 
    
    chLen = fs->scal->length;
    ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

    rangeX_begin = bx*chLen;
    rangeY_begin = by*chLen;
    rangeZ_begin = bz*chLen;   
    rangeX_end = ex*chLen;
    rangeY_end = ey*chLen;
    rangeZ_end = ez*chLen;

    FSLib->rangeX = (rangeX_end - rangeX_begin)*1e3; //(km) in LaMEM to (m) in FastScape
    FSLib->rangeY = (rangeY_end - rangeY_begin)*1e3;

//    PetscPrintf(PETSC_COMM_WORLD, "x_begin: %g, x_end: %g, y_begin: %g, y_end: %g (km); \n",
 //           rangeX_begin, rangeX_end, rangeY_begin, rangeY_end);

    PetscPrintf(PETSC_COMM_WORLD, "Lower coordinate bounds [lx, ly] : [%g, %g] [km]\n", rangeX_begin, rangeY_begin);
    PetscPrintf(PETSC_COMM_WORLD, "Upper coordinate bounds [ux, uy] : [%g, %g] [km]\n", rangeX_end, rangeY_end);

    // Gather topography and velocity
    tnodes = nx_fs * ny_fs;
    topo_fs = (PetscScalar *)malloc(tnodes * sizeof(PetscScalar)); 
 
    // topography
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
    
    ierr = VecScatterCreateToZero(surf->gtopo, &FSLib->ctx, &FSLib->gtopo_fs); CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, surf->gtopo, FSLib->gtopo_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, surf->gtopo, FSLib->gtopo_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    // velocity
    L    = (PetscInt)fs->dsz.rank;
    ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz,    &vz);     CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vz_collect,    &vz_collect);     CHKERRQ(ierr);


    START_PLANE_LOOP
    {
        vz_collect[L][j][i] = vz[L][j][i]; // (km)
 //       printf("topo_ori[%d][%d][%d]:%f\n", L, j, i, topo[L][j][i]);
  //      printf("vz_ori[%d][%d][%d]:%f\n", L, j, i, vz[L][j][i]);
    }
    END_PLANE_LOOP  

    // note: if the vec is a local vec, it will create a local output, which isn't correct
    // A global vector is needed
    ierr = VecScatterCreateToZero(FSLib->vz_collect, &FSLib->ctx, &FSLib->vz_fs); CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vz_collect, FSLib->vz_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vz_collect, FSLib->vz_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    // reallocate
    tproc = fs->dsx.nproc * fs->dsy.nproc * fs->dsz.nproc;
    //printf("total processor: %d\n", tproc);

    PetscInt para_info[7] = {sx, sy, sz, nx, ny, nz, 0};

    MPI_Comm_rank(PETSC_COMM_WORLD, &cproc);

    topo_f = topo[L][sy][sx];
    vz_f = vz[L][sy][sx];

    // send message to processor 0
    if(!ISRankZero(PETSC_COMM_WORLD))
    {
        if(0 == fs->dsz.rank)
        {
            // current order
            if(1 == fs->dsx.nproc)
            {
                rank_id = fs->dsy.rank;
            }
            else
            { 
                rank_id = fs->dsy.nproc * fs->dsy.rank + fs->dsx.rank;
            }
        }
        else
        {
            rank_id = 1000000000;
        }

        para_info[6] = rank_id;
 
        MPI_Send(para_info, 7, MPIU_INT, 0, 0, PETSC_COMM_WORLD);
        MPI_Send(&topo_f, 1, MPIU_SCALAR, 0, 1, PETSC_COMM_WORLD);
        MPI_Send(&vz_f, 1, MPIU_SCALAR, 0, 2, PETSC_COMM_WORLD); 

 //       printf("\nsend: topo_first: %f vz_first: %f, current rank: %d\n", topo_f, vz_f, cproc); 
               
 //       printf("send:  sx:%d, sy%d, sz:%d, nx:%d, ny:%d, nz:%d, id:%d, current rank: %d\n", 
 //           para_info[0], para_info[1], para_info[2], para_info[3], para_info[4], para_info[5], 
  //          para_info[6], cproc);
    }

    // Computing topography in rank 0
    if(ISRankZero(PETSC_COMM_WORLD))
    {
        PetscInt para_info_rec[tproc][7];
        PetscScalar topo_first[tproc], vz_first[tproc], topo_alloc[tnodes], vz_alloc[tnodes];
        PetscScalar* topo_pass_f = NULL;
        PetscScalar dt_max = FSLib->Max_dt; // Maximum step length, if dt_LaMEM is larger than this, use this
        PetscScalar dt_n=0.0; 
        PetscInt nsteps = dt_fs/dt_max;
        PetscInt tnodes_refine;

        PetscInt countI = 0, ind_alloc = 0, countJ;

        for(i = 0; i < tproc; i++)
        {
            if(0 == i)
            {
                for(j = 0; j < 7; j++)
                {
                    para_info_rec[i][j] = para_info[j];
                }

                topo_first[i] = topo[L][0][0];
                vz_first[i] = vz[L][0][0];
 //               printf("\nrank 0: topo_first: %f vz_first: %f\n", topo_first[i], vz_first[i]);
  //              printf("\nrank 0:  sx:%d, sy:%d, sz:%d, nx:%d, ny:%d, nz:%d, id:%d, current rank: %d\n", 
   //                 para_info_rec[i][0], para_info_rec[i][1], para_info_rec[i][2], para_info_rec[i][3], 
   //                 para_info_rec[i][4], para_info_rec[i][5], para_info_rec[i][6], i);
            }
            else
            {
                MPI_Recv(para_info_rec[i], 7, MPIU_INT, i, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&topo_first[i], 1, MPIU_SCALAR, i, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&vz_first[i], 1, MPIU_SCALAR, i, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
   //              printf("receive: topo_first: %f vz_first: %f, current rank: %d\n", topo_first[i], vz_first[i], i);                       
    //            printf("receive:  sx:%d, sy:%d, sz:%d, nx:%d, ny:%d, nz:%d, id:%d, current rank: %d\n", 
     //               para_info_rec[i][0], para_info_rec[i][1], para_info_rec[i][2], para_info_rec[i][3], 
    //                para_info_rec[i][4], para_info_rec[i][5], para_info_rec[i][6], i);
            }
        }

        // reallocate
        // All the data of a process is put into a new array in a whole, and the index may not be correct

        ierr = VecGetArray(FSLib->gtopo_fs,  &topo_fs);  CHKERRQ(ierr);
        ierr = VecGetArray(FSLib->vz_fs,  &vz_fs);  CHKERRQ(ierr);

        do
        {
            for(i = 0; i < tproc; i++)
            {
                if((topo_fs[ind_alloc] == topo_first[i]) && (vz_fs[ind_alloc] == vz_first[i]) && (para_info_rec[i][6] != 1000000000))
                {

                    countJ = 0;          

                    for(j = para_info_rec[i][1]; j < para_info_rec[i][1] + para_info_rec[i][4]; j++)
                    {
                        for(k = para_info_rec[i][0]; k < para_info_rec[i][0] + para_info_rec[i][3]; k++)
                        {
                            ind = j * nx_fs + k;
                            topo_alloc[ind] = topo_fs[countI]; 
                            vz_alloc[ind] = vz_fs[countI];
                            countI++;
                            countJ++;
//                            printf("topo_alloc[%d]:%f, ind_alloc: %d, countJ: %d\n",ind,topo_alloc[ind], ind_alloc, countJ);                          
                        }
                    }

                    ind_alloc += countJ; 
                }
            }
        } while (ind_alloc < tnodes);

/*
        for(j = 0; j < ny_fs; j++)
        {
            for(i = 0; i < nx_fs; i++)
            {
                ind = j * nx_fs + i;
                printf("topo_alloc[%d][%d]:%f\n",j,i,topo_alloc[ind]);
 //               printf("vz_alloc[%d][%d]:%f\n",j,i,vz_alloc[ind]);

            }
        }
*/
        // refine
        FSLib->nx_refine = (nx_fs - 1) * FSLib->refine + 1;
        FSLib->ny_refine = (ny_fs - 1) * FSLib->refine + 1;
        tnodes_refine = FSLib->nx_refine * FSLib->ny_refine;

        PetscScalar vz_pass[tnodes], topo_pass[tnodes], topo_refine[tnodes_refine], vz_refine[tnodes_refine];

        if(FSLib->refine > 1)    
        {
            printf("Refined times                    : %d\n", FSLib->refine);
       //     printf("nx_refine: %d, ny_refine: %d, tnodes_refine: %d \n", nx_refine, ny_refine, tnodes_refine);
            printf("Refined grid cells [nx, ny]      : [%d, %d] \n", FSLib->nx_refine, FSLib->ny_refine);

            bilinearInterpolate(FSLib, topo_alloc, topo_refine, scal, 1);
            bilinearInterpolate(FSLib, vz_alloc, vz_refine, scal, 2);


     //       printf("End refinement\n");
        }
        else
        {
            for(j = 0; j < ny_fs; j++) 
            {   
                for(i = 0; i < nx_fs; i++) 
                {
                    ind = j * nx_fs + i;
                    topo_pass[ind] = topo_alloc[ind] * 1e3; // (km) in LaMEM to (m) in FastScape
                    vz_pass[ind] = vz_alloc[ind] * scal->velocity / 1e2; // (cm/yr) in LaMEM to (m/yr) in FastScape
                } 
            }      
        }

        if(nsteps < 1) 
        {
            nsteps = 1;
            dt_max = dt_fs;
        }
        else
        {
            nsteps = nsteps + 1;
            dt_n = dt_fs-(nsteps-1)*dt_max;
        }
        
//        printf("nsteps:%d;  dt: %f, dt_n: %f\n",nsteps,dt_max,dt_n);
 //       printf("\n --------------- Begin FastScape ---------------\n");

        // store the phase that is being sedimented
        surf->phase = FSLib->sedPhases;

    //    printf("surf phase: %d\n", surf->phase);

        // get different kd, kf in different phases
        // background
        m = m + FSLib->bgphase;
        FSLib->kd = m->kd;
        FSLib->kf = m->kf;

     //   printf("kf: %g, kd: %g\n", FSLib->kf, FSLib->kd);

        // other phases (waiting)

        // needed a struct transport, too many parameters
        if(FSLib->refine > 1)
        {
            topo_pass_f = fastscapeFortran(&FSLib->nx_refine,&FSLib->ny_refine,&FSLib->rangeX,&FSLib->rangeY,&dt_max,&dt_n,&nsteps,
                vz_refine,topo_refine, &FSLib->kf, &FSLib->kfsed, &FSLib->m, &FSLib->n, &FSLib->kd, &FSLib->kdsed, &FSLib->g, 
                &FSLib->gsed, &FSLib->p, &FSLib->setMarine, &FSLib->sealevel, &FSLib->poro_silt, &FSLib->poro_sand, &FSLib->zporo_silt, &FSLib->zporo_sand,
                &FSLib->ratio, &FSLib->Lsolve, &FSLib->kds_silt, &FSLib->kds_sand, &FSLib->FS_BC);
        }
        else
        {
            topo_pass_f = fastscapeFortran(&nx_fs,&ny_fs,&FSLib->rangeX,&FSLib->rangeY,&dt_max,&dt_n,&nsteps,vz_pass,topo_pass,
                &FSLib->kf, &FSLib->kfsed, &FSLib->m, &FSLib->n, &FSLib->kd, &FSLib->kdsed, &FSLib->g, &FSLib->gsed, &FSLib->p, 
                &FSLib->setMarine, &FSLib->sealevel, &FSLib->poro_silt, &FSLib->poro_sand, &FSLib->zporo_silt, &FSLib->zporo_sand, &FSLib->ratio, 
                &FSLib->Lsolve, &FSLib->kds_silt, &FSLib->kds_sand, &FSLib->FS_BC);
        }

        // save high-resolution topography in processor 0
        // create directory(encode current time & steo number)
        // Correction
        if(FSLib->refine > 1)
        {
            PetscInt status, countX, countY;
            countX = 0; 
            countY = 0;

            if(pvsurf->outsurf_refine)
            {
                // update time stamp and counter
                step_fs++;

                asprintf(&dirName, "Timestep_%1.8lld_%1.8e", (LLD)step_fs, time_fs);

                // create output directory
                #ifdef _WIN32
                        // call this on windows machines
                        status = mkdir(name);
                #else
                        // standard access pattern drwxr-xr-x
                        status = mkdir(dirName, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
                #endif 
                        if(status && errno != EEXIST)
                        {
                            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to create directory %s", dirName);
                        }
            }

            for(j = 0; j < FSLib->ny_refine; j += FSLib->refine) 
            {   
                for(i = 0; i < FSLib->nx_refine; i += FSLib->refine) 
                {
                    ind = j * FSLib->nx_refine+i;  

                    if(0 == ind%FSLib->nx_refine) 
                    {
                        countX = 0;
                        if(0 != ind) 
                        {
                            countY++; 
                        }
                    }

                    ind2 = countY * nx_fs + countX;

                    topo_fs[ind2] = topo_pass_f[ind]/1e3; // m(FastScape) to km(LaMEM)

                    if(topo_fs[ind2] > rangeZ_end) topo_fs[ind2] = rangeZ_end;
                    if(topo_fs[ind2] < rangeZ_begin) topo_fs[ind2] = rangeZ_begin;

                    countX++;                    
                }
            }

            // save output in a high-resolution grid
            if(pvsurf->outsurf_refine)
            {
                ierr = savePvtsFS(pvsurf, FSLib, time_fs, step_fs, dirName, topo_pass_f); CHKERRQ(ierr);
            }
        } 
        else
        {
            for(j = 0; j < ny_fs; j++) 
            {   
                for(i = 0; i < nx_fs; i++) 
                {
                    ind = j * nx_fs+i;

                    topo_fs[ind] = topo_pass_f[ind]/1e3; // m(FastScape) to km(LaMEM) 
     
                    if(topo_fs[ind] > rangeZ_end) topo_fs[ind] = rangeZ_end;
                    if(topo_fs[ind] < rangeZ_begin) topo_fs[ind] = rangeZ_begin;
               }
            }
        }
        topo_pass_f = NULL;
        clearArray();
    }

    // Broadcast      
    if(ISParallel(PETSC_COMM_WORLD))
    {
        ierr = MPI_Bcast(topo_fs, (PetscMPIInt)tnodes, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);
//        ierr = MPI_Bcast(topo_alloc, (PetscMPIInt)tnodes, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);        
    }

    // Save topography in different ranks
    START_PLANE_LOOP
    {
        topo[L][j][i] = topo_fs[j * nx_fs + i]; // (km)
    }
    END_PLANE_LOOP  

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz,  &vz);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vz_collect,  &vz_collect);  CHKERRQ(ierr);

    // compute ghosted version of the topography
    GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

    // compute & store average topography
    ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

    ierr = VecScatterDestroy(&FSLib->ctx); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->gtopo_fs);   CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->vz_fs);   CHKERRQ(ierr);

//    PetscPrintf(PETSC_COMM_WORLD,"\n---------------- FastScape Done ---------------\n");
    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

    PetscFunctionReturn(0);
}

