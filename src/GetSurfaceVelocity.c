/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de
 
GetSurfaceVelocity.c

Created on: 09.02.2012
Author: tobibaumann

This sourcefile includes following functions:
 - 	GetMisfit_HVel
 -	GetMisfit_Gravity

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "GetSurfaceVelocity.h"
#include "ParaViewOutput.h"
#include "Utils.h"

#undef __FUNCT__
#define __FUNCT__ "GetSurfaceVelocity"
PetscErrorCode GetSurfaceVelocity(UserContext *user,DM da_nodes,Vec gvec_Vel)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	//  --- general ---
	PetscErrorCode		ierr;
	PetscMPIInt			rank,surfrank;
	MPI_Comm			SURF_COMM;
	PetscInt			SurfProc;


	//  --- data ---
	PetscInt 			xs,ys,zs,xm,ym,zm,nnode_x,nnode_y,nnode_z;
	PetscInt			ix,iy,iz;

	Field				***field_velocity;

	DMDALocalInfo 		linfo;

	Vec 				lvec_Vel,lvec_Vxall,lvec_Vyall;
	Vec 				gvec_Vxall,gvec_Vyall,gvec_Qualall;
	PetscScalar         *array_Vxall,*array_Vyall;
	PetscScalar 		global_VelMisfit;


	// 	--- file i/o ---
	PetscViewer			view_out;
	char				filename_out[PETSC_MAX_PATH_LEN];


//  --- Code ----------------------------------------------------------------------------------------------------------


// 1. === get local process information ===

	// MPI communication
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// get dimensions
	ierr = DMDAGetLocalInfo(da_nodes,&linfo);

	// if process owns velocity-nodes at the model surface
	if (linfo.zs+linfo.zm == linfo.mz) SurfProc = 1; else SurfProc = 0;

	// Create communicator for processes owning surface nodes
	MPI_Comm_split(PETSC_COMM_WORLD,SurfProc,rank,&SURF_COMM);
	MPI_Comm_rank(SURF_COMM,&surfrank);


// 2. === extract surface velocity field ===

	// --- PETSC_COMM_WORLD ---
	// get indices and number of nodes
	ierr = DMDAGetInfo(da_nodes, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_nodes, &xs, &ys, &zs, &xm, &ym, &zm); 				CHKERRQ(ierr);

	// get local vectors and arrays
	ierr = DMGetLocalVector(da_nodes,&lvec_Vel); 								CHKERRQ(ierr);
	ierr = VecZeroEntries(lvec_Vel);											CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da_nodes,gvec_Vel, INSERT_VALUES, lvec_Vel); 	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_nodes,gvec_Vel, INSERT_VALUES, lvec_Vel); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_nodes, lvec_Vel,&field_velocity); 					CHKERRQ(ierr);
	// --- PETSC_COMM_WORLD END---


	if(SurfProc){
		ierr = VecCreateSeq(PETSC_COMM_SELF,linfo.mx*linfo.my,&lvec_Vxall);CHKERRQ(ierr);
		ierr = VecCreateSeq(PETSC_COMM_SELF,linfo.mx*linfo.my,&lvec_Vyall);CHKERRQ(ierr);

		ierr = VecSet(lvec_Vxall,0.0); CHKERRQ(ierr);
		ierr = VecSet(lvec_Vyall,0.0); CHKERRQ(ierr);

		ierr = VecGetArray(lvec_Vxall,&array_Vxall);CHKERRQ(ierr);
		ierr = VecGetArray(lvec_Vyall,&array_Vyall);CHKERRQ(ierr);
	}

	// --- PETSC_COMM_WORLD ---
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){

				// extract SURFACE VALUES
				if (iz== (nnode_z-1)){

					array_Vxall[iy*(linfo.mx-1)+ix+iy] 	= field_velocity[iz][iy][ix].Vx;
					array_Vyall[iy*(linfo.mx-1)+ix+iy] 	= field_velocity[iz][iy][ix].Vy;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(da_nodes,lvec_Vel,&field_velocity); 	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da_nodes,&lvec_Vel); 				CHKERRQ(ierr);
	// --- PETSC_COMM_WORLD END---


	if(SurfProc){

		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,linfo.mx*linfo.my,array_Vxall,&lvec_Vxall);CHKERRQ(ierr);			ierr = VecAssemblyBegin(lvec_Vxall); CHKERRQ(ierr);		ierr = VecAssemblyEnd(lvec_Vxall); CHKERRQ(ierr);
		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,linfo.mx*linfo.my,array_Vyall,&lvec_Vyall);CHKERRQ(ierr);			ierr = VecAssemblyBegin(lvec_Vyall); CHKERRQ(ierr);		ierr = VecAssemblyEnd(lvec_Vyall); CHKERRQ(ierr);

		ierr = VecRestoreArray(lvec_Vxall,&array_Vxall);CHKERRQ(ierr);
		ierr = VecRestoreArray(lvec_Vyall,&array_Vyall);CHKERRQ(ierr);

	}


//	--- gathering process ---------------------------------------------------------------------------------------------

	if(SurfProc){

	// store all in all local vectors
	ierr = Sum_Vectors(SURF_COMM ,&lvec_Vxall,linfo.mx*linfo.my);							CHKERRQ(ierr);
	ierr = Sum_Vectors(SURF_COMM ,&lvec_Vyall,linfo.mx*linfo.my);							CHKERRQ(ierr);

	// scatter sequential vector lvec_x -> distributed vector gvec_x
	ierr = VecCreateMPI(SURF_COMM,PETSC_DECIDE,linfo.mx*linfo.my,&gvec_Vxall);				CHKERRQ(ierr);
	ierr = VecCreateMPI(SURF_COMM,PETSC_DECIDE,linfo.mx*linfo.my,&gvec_Vyall);				CHKERRQ(ierr);
	ierr = VecCreateMPI(SURF_COMM,PETSC_DECIDE,linfo.mx*linfo.my,&gvec_Qualall);			CHKERRQ(ierr);
	ierr = VecSet(gvec_Qualall,1); 															CHKERRQ(ierr);

	ierr = VecSeq2VecMpi(surfrank,lvec_Vxall,&gvec_Vxall);									CHKERRQ(ierr);
	ierr = VecSeq2VecMpi(surfrank,lvec_Vyall,&gvec_Vyall);									CHKERRQ(ierr);

	// destroy sequential vectors
	ierr = VecDestroy(&lvec_Vxall); 														CHKERRQ(ierr);
	ierr = VecDestroy(&lvec_Vyall); 														CHKERRQ(ierr);

	}

// 	--- save binaries -------------------------------------------------------------------------------------------------

	if((SurfProc) && (user->SurfVelField.SaveRef==1)){

		sprintf(filename_out,"SurfaceVelocity_%s_REF_%lld.dat",user->OutputFile,(LLD)user->Optimisation.mpi_group_id);
		if(surfrank == 0){
			ierr = PetscPrintf(PETSC_COMM_SELF,"#  -- save reference data -- \n");				CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"#     save file: %s \n",filename_out);		CHKERRQ(ierr);
		}

		ierr = PetscViewerBinaryOpen(SURF_COMM,filename_out , FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
		ierr = VecView(gvec_Vxall,			view_out); 										CHKERRQ(ierr);
		ierr = VecView(gvec_Vyall,			view_out); 										CHKERRQ(ierr);
		ierr = VecView(gvec_Qualall,		view_out); 										CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out); 												CHKERRQ(ierr);
	}

	// wait for all nodes
	MPI_Barrier(PETSC_COMM_WORLD);


// 	--- get misfit ----------------------------------------------------------------------------------------------------
	if(user->Optimisation.GetIt==1){
		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- get surface velocity misfit -----------\n");


		if(SurfProc){
			GetMisfit_SurfaceVelocity(user,gvec_Vxall,gvec_Vyall,SURF_COMM);
		}

		// wait for all nodes
		MPI_Barrier(PETSC_COMM_WORLD);

		// send misfit to global root node and broadcast to all
		if (surfrank==0) global_VelMisfit = user->Optimisation.MisfitSurfVel; else global_VelMisfit = 0.0;

		//PetscPrintf(PETSC_COMM_SELF,"--- DEBUG rank %lld misift %g\n",(LLD)rank,global_VelMisfit);
		ierr = MPI_Allreduce(&global_VelMisfit,&user->Optimisation.MisfitSurfVel,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);	CHKERRQ(ierr);
		//PetscPrintf(PETSC_COMM_SELF,"--- DEBUG rank %lld misift %g\n",(LLD)rank,global_VelMisfit);

	}
// 	--- clean up ------------------------------------------------------------------------------------------------------

	// --- Clean Vectors ---
	if(SurfProc){
		ierr = VecDestroy(&gvec_Vxall); 		CHKERRQ(ierr);
		ierr = VecDestroy(&gvec_Vyall); 		CHKERRQ(ierr);
		ierr = VecDestroy(&gvec_Qualall); 		CHKERRQ(ierr);
	}

	// --- Clean communicators
	MPI_Comm_free(&SURF_COMM);

	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "VecSeq2VecMpi"
PetscErrorCode VecSeq2VecMpi(PetscMPIInt rank,Vec lvec,Vec *gvec)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
	PetscScalar 		*array_l;
	PetscInt			i,size_lvec,*array_num;

//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	ierr = VecGetSize(lvec,&size_lvec);CHKERRQ(ierr);


	if(rank==0){
		ierr = VecGetArray(lvec,&array_l);CHKERRQ(ierr);
		ierr = PetscMalloc((size_t)size_lvec*sizeof(PetscInt),&array_num);CHKERRQ(ierr);

		for (i=0;i<size_lvec;i++){
			array_num[i] = i;
		}

		ierr = VecSetValues(*gvec,size_lvec,array_num,array_l,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecRestoreArray(lvec,&array_l);	CHKERRQ(ierr);
		ierr = PetscFree(array_num); 			CHKERRQ(ierr);
	}


	ierr = VecAssemblyBegin(*gvec);  		CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*gvec);  			CHKERRQ(ierr);



	PetscFunctionReturn(0);
}



//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "Sum_Vectors"
PetscErrorCode Sum_Vectors(MPI_Comm G_COMM ,Vec *lvec_Vall,PetscInt vec_length)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	Vec					gvec_Vall;
	VecScatter			scattercontext;
	IS              	isglobal,islocal;
	PetscErrorCode		ierr;

//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	// create global vector of the same size as local vectors
	ierr = VecCreateMPI(G_COMM,PETSC_DECIDE,vec_length,&gvec_Vall);								CHKERRQ(ierr);
	ierr = VecSet(gvec_Vall,0.0);																CHKERRQ(ierr);
	ierr = VecAssemblyBegin(gvec_Vall);  														CHKERRQ(ierr);
	ierr = VecAssemblyEnd(gvec_Vall);  															CHKERRQ(ierr);

	// create indices
	ierr = ISCreateStride(PETSC_COMM_SELF,vec_length,0,1,&islocal);								CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_SELF,vec_length,0,1,&isglobal);							CHKERRQ(ierr);

	// gather process
	ierr = VecScatterCreate(*lvec_Vall,islocal,gvec_Vall,isglobal,&scattercontext);				CHKERRQ(ierr);
	ierr = VecScatterBegin(scattercontext,*lvec_Vall,gvec_Vall,ADD_VALUES,SCATTER_FORWARD);		CHKERRQ(ierr);
	ierr = VecScatterEnd(scattercontext  ,*lvec_Vall,gvec_Vall,ADD_VALUES,SCATTER_FORWARD);		CHKERRQ(ierr);
	ierr = VecScatterDestroy(&scattercontext);													CHKERRQ(ierr);

	// back scatter to local vector
	ierr = VecScatterCreate(gvec_Vall,isglobal,*lvec_Vall,islocal,&scattercontext);				CHKERRQ(ierr);
	ierr = VecScatterBegin(scattercontext,gvec_Vall,*lvec_Vall,INSERT_VALUES,SCATTER_FORWARD);	CHKERRQ(ierr);
	ierr = VecScatterEnd(scattercontext,gvec_Vall,*lvec_Vall,INSERT_VALUES,SCATTER_FORWARD);	CHKERRQ(ierr);
	ierr = VecScatterDestroy(&scattercontext);													CHKERRQ(ierr);

	// Destroy indices
	ierr = ISDestroy(&isglobal);																CHKERRQ(ierr);
	ierr = ISDestroy(&islocal);																	CHKERRQ(ierr);

	// Destroy global vector
	ierr = VecDestroy(&gvec_Vall);																CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetMisfit_SurfaceVelocity"
PetscErrorCode GetMisfit_SurfaceVelocity(UserContext *user,Vec gvec_Vx, Vec gvec_Vy,MPI_Comm SURF_COMM)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	// --- general ---
	PetscErrorCode		ierr;

	// --- file i/o ---
	PetscViewer			view_in;
	char				filename_in[PETSC_MAX_PATH_LEN];

	//  --- misfit calculation ---
	Vec 				gvec_VxRef,gvec_VyRef,gvec_Vqual;
	PetscScalar 		err_Vx,err_Vy,max_VxRefAbs,max_VyRefAbs;

	PetscInt			size_VxRef,i,num,global_num, NumOfPoints;
	PetscScalar			*garray_VxRef,*garray_Vx;
	PetscScalar			*garray_VyRef,*garray_Vy;
	PetscScalar 		*garray_Vqual;


//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

//	0.) --- initialize variables ---

	max_VxRefAbs 			= 0.0;
	max_VyRefAbs 			= 0.0;
	err_Vx 					= 0.0;
	err_Vy 					= 0.0;



// 1.) --- load reference vectors ---

	// file name for reference data
	//sprintf(filename_in,"SurfaceVelocity_%s_REF_%lld.dat",user->OutputFile,(LLD)user->Optimisation.mpi_group_id);
	sprintf(filename_in,"%s_%lld.dat",user->SurfVelField.RefDatFile2load,(LLD)user->Optimisation.mpi_group_id);

					/* --- DEBUG ---
					// load data you just saved
						PetscPrintf(PETSC_COMM_WORLD,"--- DEBUG: load debug data: \n");
						sprintf(filename_in,"TEST_parallel.dat");
						sprintf(filename_in,"SurfaceVelocity_%s_REF_%lld.dat",user->OutputFile,(LLD)user->mpi_group_id);
						//sprintf(filename_in,"TTT.dat");
					// --- DEBUG end ---*/

	// initialize & zero vectors to store data
	ierr = VecDuplicate(gvec_Vx,&gvec_VxRef);											CHKERRQ(ierr);
	ierr = VecDuplicate(gvec_Vy,&gvec_VyRef);											CHKERRQ(ierr);
	ierr = VecDuplicate(gvec_Vy,&gvec_Vqual);											CHKERRQ(ierr);

	ierr = VecZeroEntries(gvec_VxRef);													CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_VyRef);													CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_Vqual);													CHKERRQ(ierr);


					/* --- DEBUG ---
						PetscPrintf(PETSC_COMM_SELF,"(OPEN)==> rank %lld: , mz: %lld, zs+zm: %lld; filename: %s\n",(LLD)rank,(LLD)linfo.mz,(LLD)(linfo.zs+linfo.zm), filename_in);
					// --- DEBUG end---*/

	// load data
	ierr = PetscViewerBinaryOpen(SURF_COMM,filename_in,FILE_MODE_READ,&view_in);		CHKERRQ(ierr);
	ierr = PetscPrintf(SURF_COMM,"#  -- load reference data --\n");						CHKERRQ(ierr);
	ierr = PetscPrintf(SURF_COMM,"#     load file: %s\n",filename_in);					CHKERRQ(ierr);
	ierr = VecLoad(gvec_VxRef,view_in);													CHKERRQ(ierr);
	ierr = VecLoad(gvec_VyRef,view_in);													CHKERRQ(ierr);
	ierr = VecLoad(gvec_Vqual,view_in);													CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in);												CHKERRQ(ierr);

					/* --- DEBUG ---
						PetscPrintf(PETSC_COMM_SELF,"Vxref \n");
						VecView(gvec_VxRef, PETSC_VIEWER_STDOUT_WORLD);
						PetscPrintf(PETSC_COMM_SELF,"Vx \n");
						VecView(gvec_Vx, PETSC_VIEWER_STDOUT_WORLD);

						PetscPrintf(PETSC_COMM_SELF,"Vyref \n");
						VecView(gvec_VyRef, PETSC_VIEWER_DRAW_WORLD);
						PetscPrintf(PETSC_COMM_SELF,"Vy \n");
						VecView(gvec_Vy, PETSC_VIEWER_DRAW_WORLD);


						PetscPrintf(PETSC_COMM_SELF,"Vqual \n");
						VecView(gvec_Vqual, PETSC_VIEWER_STDOUT_WORLD);
					// --- DEBUG end ---*/

// 2.) --- check whether reference velocities contain NANs ---
	/* There are a couple of possibilities why a vel-node contain NAN:
	 *  - it's a boundary node which was set to NAN to prevent the total error being influenced by the boundary conditions
	 *  - the vel-node was set to NAN on purpose since there is no vel data available
	 *  - it's NAN because interpolation failed because it lacked for data while creating a the ref velocities
	 */

	ierr = VecGetLocalSize(gvec_VxRef,&size_VxRef);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_VxRef	,&garray_VxRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vx		,&garray_Vx);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_VyRef	,&garray_VyRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vy		,&garray_Vy);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vqual 	,&garray_Vqual);									CHKERRQ(ierr);


	num  = 0;
	ierr =  VecGetSize(gvec_VxRef,&NumOfPoints);										CHKERRQ(ierr);

	for (i=0; i<size_VxRef; i++) {

					/* --- DEBUG ---
						if(rank==8) PetscPrintf(PETSC_COMM_SELF,"--- DEBUG : local %lld:  larray_VxRef[%lld]=  %g \n",(LLD)rank,(LLD)i,larray_VxRef[i]);
						if(rank==8) PetscPrintf(PETSC_COMM_SELF,"--- DEBUG : local %lld:  larray_Vx[%lld]=  %g \n",(LLD)rank, (LLD)i,larray_Vx[i]);
					// --- DEBUG end ---*/

		if(isnan(garray_VxRef[i])){

			if(isnan(garray_VyRef[i])==0) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VyRef is a number while VxRef is NaN !! ==> CHECK input data\n");

			// set NAN-components to zero so that Petsc vector routines work fine
			// also set the corresponding grid velocity to zero so that total error isn't effected
			garray_VxRef[i] 	= 0.0;
			garray_VyRef[i] 	= 0.0;
			garray_Vx[i] 		= 0.0;
			garray_Vy[i] 		= 0.0;
			garray_Vqual[i] 	= 0.0;

						/* --- DEBUG ---
							PetscPrintf(PETSC_COMM_SELF,"> rank: %lld; garray_Vxref[%lld]= %g \n",(LLD)rank,(LLD)i,garray_VxRef[i]);
						// --- DEBUG end ---*/

			// if velocity at this node is NAN, also the total number of points needs to be corrected.
			num = num + 1;

		}
		else{
			if(isnan(garray_VyRef[i])) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VyRef is NaN while VxRef is not! ==> CHECK input data\n");
					/* --- DEBUG ---
						PetscPrintf(PETSC_COMM_SELF,"--- DEBUG larray_VxRef[%lld]= %g ,larray_Vx[%lld]= %g \n",(LLD)i,larray_VxRef[i],(LLD)i,larray_Vx[i]);
						PetscPrintf(PETSC_COMM_SELF,"--- DEBUG larray_VyRef[%lld]= %g ,larray_Vy[%lld]= %g \n",(LLD)i,larray_VyRef[i],(LLD)i,larray_Vy[i]);
					// --- DEBUG end ---*/
		}
	}

	ierr = MPI_Allreduce(&num,&global_num,1,MPIU_INT,MPI_SUM,SURF_COMM);					CHKERRQ(ierr);
	NumOfPoints = NumOfPoints - global_num;

	ierr = VecRestoreArray(gvec_Vqual	,&garray_Vqual);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_Vx		,&garray_Vx);									CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_VxRef	,&garray_VxRef);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_Vy		,&garray_Vy);									CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_VyRef	,&garray_VyRef);								CHKERRQ(ierr);

// 3.) --- Normalization --

	// global max(|VxRef|), global max(|VyRef|)
		ierr = VecAbsMax(gvec_VxRef,&max_VxRefAbs);											CHKERRQ(ierr);
		ierr = VecAbsMax(gvec_VyRef,&max_VyRefAbs);											CHKERRQ(ierr);

						/* --- DEBUG---
							PetscPrintf(SURF_COMM,"--- DEBUG vxref max: %g \n",max_VxRefAbs);
							PetscPrintf(SURF_COMM,"--- DEBUG vyref max: %g \n",max_VyRefAbs);
						// --- DEBUG end ---*/


		// Normalize with the maximum of the reference "true" signal
		ierr = VecScale(gvec_VxRef,1/max_VxRefAbs);											CHKERRQ(ierr);
		ierr = VecScale(gvec_Vx   ,1/max_VxRefAbs);											CHKERRQ(ierr);
		ierr = VecScale(gvec_VyRef,1/max_VyRefAbs);											CHKERRQ(ierr);
		ierr = VecScale(gvec_Vy   ,1/max_VyRefAbs);											CHKERRQ(ierr);

						/* --- DEBUG---
							 PetscPrintf(SURF_COMM,"==> Normed: VxRef\n");
							 if (surfrank ==0 ) VecView(lvec_VxRef, PETSC_VIEWER_STDOUT_SELF);
							 PetscPrintf(SURF_COMM,"==> Normed: Vx\n");
							 if (surfrank ==0 ) VecView(lvec_Vx, PETSC_VIEWER_STDOUT_SELF);
						// --- DEBUG end ---*/


// 4.) --- calculate error ---

		//  lvec_VxRef = (VxRef-Vx) ,		lvec_VyRef = (VyRef-Vy)
		ierr = VecAXPBY(gvec_VxRef,1.0,-1.0	,gvec_Vx);										CHKERRQ(ierr);
		ierr = VecAXPBY(gvec_VyRef,1.0,-1.0	,gvec_Vy);										CHKERRQ(ierr);

						/* --- DEBUG---
							 PetscPrintf(SURF_COMM,"==> Error: \n");
							 if (rank ==0 ) VecView(lvec_VxRef, PETSC_VIEWER_STDOUT_SELF);
							 if (rank ==0 ) VecView(lvec_VxRef, PETSC_VIEWER_STDOUT_SELF);
						// --- DEBUG end ---*/


		//lvec_VxRef = vx.^2,  				lvec_VyRef = vy.^2
		ierr = VecPointwiseMult(gvec_VxRef,gvec_VxRef,gvec_VxRef);							CHKERRQ(ierr);
		ierr = VecPointwiseMult(gvec_VyRef,gvec_VyRef,gvec_VyRef);							CHKERRQ(ierr);

		// lvec_VxRef = vx.^2 * quality, 	lvec_VyRef = vy.^2 * quality
		ierr = VecPointwiseMult(gvec_VxRef,gvec_VxRef,gvec_Vqual);							CHKERRQ(ierr);
		ierr = VecPointwiseMult(gvec_VyRef,gvec_VyRef,gvec_Vqual);							CHKERRQ(ierr);

		// lvec_VxRef = sum(vx.^2 *quality), lvec_VyRef = sum(vyx.^2 *quality)
		ierr = VecSum(gvec_VxRef,&err_Vx);													CHKERRQ(ierr);
		ierr = VecSum(gvec_VyRef,&err_Vy);													CHKERRQ(ierr);


	// RMS = sqrt(1/N * sum(x_i.^2)  )
	user->Optimisation.MisfitSurfVel = sqrt( (err_Vx + err_Vy)/(2*NumOfPoints) );


					//* --- DEBUG---
						PetscPrintf(SURF_COMM,"#     Global err_Vx: %g \n",err_Vx);
						PetscPrintf(SURF_COMM,"#     Global err_Vy: %g \n",err_Vy);
						PetscPrintf(SURF_COMM,"#     NumberofPoints: %lld \n",(LLD)NumOfPoints);
					// --- DEBUG end ---*/

					/* --- DEBUG ---
						PetscPrintf(SURF_COMM,"--- DEBUG ==> velocity misfit: %g \n",user->misfit_Velocity);
					// --- DEBUG end ---*/

	// 5.) --- clean reference vectors ---
	ierr = VecDestroy(&gvec_VxRef); 													CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_VyRef); 													CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_Vqual); 													CHKERRQ(ierr);



	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "VecAbsMax"
PetscErrorCode VecAbsMax(Vec gvec,PetscScalar *max)
{

	PetscErrorCode	ierr;
	PetscInt		ind;
	Vec 			gvec_abs;


	PetscFunctionBegin;

	ierr = VecDuplicate(gvec,&gvec_abs); 			CHKERRQ(ierr);
	ierr = VecCopy(gvec,gvec_abs);					CHKERRQ(ierr);

	ierr = VecAbs(gvec_abs);						CHKERRQ(ierr);
	ierr = VecMax(gvec_abs,&ind,max);				CHKERRQ(ierr);

	ierr = VecDestroy(&gvec_abs); 			 		CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/**********************************************************************************************************************
 *
 *
 *	functions that handle the internal free surface ...
 *
 *
 *
 **********************************************************************************************************************/

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetSurfaceVelocity_InternalFreeSurface"
PetscErrorCode GetSurfaceVelocity_InternalFreeSurface(UserContext *user,PetscInt itime)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
	Vec					gvec_Vx,gvec_Vy,gvec_Vz,gvec_Topo,gvec_Quali;
	MPI_Comm			PLANE_COMM;
	PetscViewer			view_out;
	DM					da_plane;
	char 				*FileName,*DirectoryName;
	PetscInt			gp=0;
	PetscMPIInt			planerank=0;
	PetscScalar			global_SumAbsVel=0.0,global_SumSqrsVel=0.0;

//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;


	user->Optimisation.SumAbsVel  = 0.0;
	user->Optimisation.SumSqrsVel = 0.0;


	// Get DMDA2d for including gridpoint 0
	ierr = DMDACreate2dFrom3d(user->DA_SurfaceTopography,gp,DMDA_Z,&da_plane,&PLANE_COMM);		CHKERRQ(ierr);

	// Get associated Vecs
	ierr = DMExtractGlobalVec2dFromGlobalVec3d(user->DA_SurfaceTopography,user->SurfaceTopography_Vx,gp,da_plane,PLANE_COMM,&gvec_Vx);	CHKERRQ(ierr);
	ierr = DMExtractGlobalVec2dFromGlobalVec3d(user->DA_SurfaceTopography,user->SurfaceTopography_Vy,gp,da_plane,PLANE_COMM,&gvec_Vy);	CHKERRQ(ierr);
	ierr = DMExtractGlobalVec2dFromGlobalVec3d(user->DA_SurfaceTopography,user->SurfaceTopography_Vz,gp,da_plane,PLANE_COMM,&gvec_Vz);	CHKERRQ(ierr);
	ierr = DMExtractGlobalVec2dFromGlobalVec3d(user->DA_SurfaceTopography,user->SurfaceTopography   ,gp,da_plane,PLANE_COMM,&gvec_Topo);	CHKERRQ(ierr);

	// Create dummy quali vector
	if(PLANE_COMM!=MPI_COMM_NULL){
		ierr = VecDuplicate(gvec_Vx,&gvec_Quali);												CHKERRQ(ierr);
		ierr = VecSet(gvec_Quali,1.0); 															CHKERRQ(ierr);
	}

	// Save binaries
	if(user->SurfVelField.SaveRef==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save reference data -- \n");					CHKERRQ(ierr);
		asprintf(&DirectoryName, "ReferenceData_%1.6lld",(LLD)itime);
		ierr = LaMEMCreateOutputDirectory(DirectoryName); 										CHKERRQ(ierr);
		asprintf( &FileName,"%s/REF_FreeSurfaceVel.bin",DirectoryName);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     save file: %s \n",FileName);					CHKERRQ(ierr);

		if(PLANE_COMM!=MPI_COMM_NULL){
		ierr = PetscViewerBinaryOpen(PLANE_COMM,FileName, FILE_MODE_WRITE, &view_out);			CHKERRQ(ierr);
		ierr = VecView(gvec_Vx,view_out);														CHKERRQ(ierr);
		ierr = VecView(gvec_Vy,view_out);														CHKERRQ(ierr);
		ierr = VecView(gvec_Vz,view_out);														CHKERRQ(ierr);
	  //ierr = VecView(gvec_Topo,view_out);														CHKERRQ(ierr);
		ierr = VecView(gvec_Quali,view_out);													CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out);													CHKERRQ(ierr);
		}
		free(FileName);
		free(DirectoryName);
	}

	// Get misfit
	if(user->Optimisation.GetIt==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");	CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#-- get surface velocity misfit -----------\n");	CHKERRQ(ierr);

		if(PLANE_COMM!=MPI_COMM_NULL){
		ierr = GetMisfit_SurfaceVelocity_InternalFreeSurface(user,gvec_Vx,gvec_Vy,gvec_Vz,PLANE_COMM);CHKERRQ(ierr);
		}

		// send misfit to all nodes � only nodes of PLANE_COMM know about the misfit
		planerank = 1;
		MPI_Comm_rank(PLANE_COMM,&planerank);

		if (planerank!=0){// only root node of PLANE_COMM has global_VelMisfit !=0.0
			global_SumAbsVel = 0.0;
			global_SumSqrsVel = 0.0;
		}
		else{
			global_SumAbsVel = user->Optimisation.SumAbsVel;
			global_SumSqrsVel = user->Optimisation.SumSqrsVel;
		}


		ierr = MPI_Allreduce(&global_SumAbsVel,&user->Optimisation.SumAbsVel,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
		ierr = MPI_Allreduce(&global_SumSqrsVel,&user->Optimisation.SumSqrsVel,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	}

	// Destroy Vecs
	ierr = DMDestroyGlobalVec2dFromGlobalVec3d(gvec_Vx,PLANE_COMM);								CHKERRQ(ierr);
	ierr = DMDestroyGlobalVec2dFromGlobalVec3d(gvec_Vy,PLANE_COMM);								CHKERRQ(ierr);
	ierr = DMDestroyGlobalVec2dFromGlobalVec3d(gvec_Vz,PLANE_COMM);								CHKERRQ(ierr);
	ierr = DMDestroyGlobalVec2dFromGlobalVec3d(gvec_Topo,PLANE_COMM);							CHKERRQ(ierr);
	ierr = DMDestroyGlobalVec2dFromGlobalVec3d(gvec_Quali,PLANE_COMM);							CHKERRQ(ierr);

	// Destroy DMDA2d + PLANE_COMM
	ierr = DMDADestroy2dFrom3d(da_plane,PLANE_COMM);											CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetMisfit_SurfaceVelocity_InternalFreeSurface"
PetscErrorCode GetMisfit_SurfaceVelocity_InternalFreeSurface(UserContext *user,Vec gvec_Vx, Vec gvec_Vy,Vec gvec_Vz,MPI_Comm PLANE_COMM)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	// --- general ---
	PetscErrorCode		ierr;

	// --- file i/o ---
	PetscViewer			view_in;
	char				filename_in[PETSC_MAX_PATH_LEN];

	//  --- misfit calculation ---
	Vec 				gvec_VxRef,gvec_VyRef,gvec_VzRef,gvec_Vqual;
	PetscScalar 		err_Vx,err_Vy,err_Vz;
	PetscScalar         l1_Vx,l1_Vy,l1_Vz;

	PetscInt			size_VxRef,i,num,global_num, NumOfPoints;
	PetscScalar			*garray_VxRef,*garray_Vx;
	PetscScalar			*garray_VyRef,*garray_Vy;
	PetscScalar			*garray_VzRef,*garray_Vz;
	PetscScalar 		*garray_Vqual;


//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

//	0.) --- initialize variables ---
	err_Vx 					= 0.0;
	err_Vy 					= 0.0;
	err_Vz 					= 0.0;
	l1_Vx                   = 0.0;
	l1_Vy                   = 0.0;
	l1_Vz                   = 0.0;

// 1.) --- load reference vectors ---
	// file name for reference data
	//sprintf(filename_in,"%s_%lld.dat",user->SurfVelField.RefDatFile2load,(LLD)user->Optimisation.mpi_group_id);
	sprintf(filename_in,"%s",user->SurfVelField.RefDatFile2load);

	// initialize & zero vectors to store data
	ierr = VecDuplicate(gvec_Vx,&gvec_VxRef);											CHKERRQ(ierr);
	ierr = VecDuplicate(gvec_Vy,&gvec_VyRef);											CHKERRQ(ierr);
	ierr = VecDuplicate(gvec_Vz,&gvec_VzRef);											CHKERRQ(ierr);
	ierr = VecDuplicate(gvec_Vx,&gvec_Vqual);											CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_VxRef);													CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_VyRef);													CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_VzRef);													CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_Vqual);													CHKERRQ(ierr);

	// load data
	ierr = PetscViewerBinaryOpen(PLANE_COMM,filename_in,FILE_MODE_READ,&view_in);		CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#  -- load reference data --\n");					CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     load file: %s\n",filename_in);					CHKERRQ(ierr);
	ierr = VecLoad(gvec_VxRef,view_in);													CHKERRQ(ierr);
	ierr = VecLoad(gvec_VyRef,view_in);													CHKERRQ(ierr);
	ierr = VecLoad(gvec_VzRef,view_in);													CHKERRQ(ierr);
	ierr = VecLoad(gvec_Vqual,view_in);													CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in);												CHKERRQ(ierr);

// 2.) --- check whether reference velocities contain NANs ---
	/* There are a couple of possibilities why a vel-node contain NAN:
	 *  - it's a boundary node which was set to NAN to prevent the total error being influenced by the boundary conditions
	 *  - the vel-node was set to NAN on purpose since there is no vel data available
	 *  - it's NAN because interpolation failed because it lacked for data while creating a the ref velocities
	 */

	ierr = VecGetLocalSize(gvec_VxRef,&size_VxRef);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_VxRef	,&garray_VxRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vx		,&garray_Vx);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_VyRef	,&garray_VyRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vy		,&garray_Vy);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_VzRef	,&garray_VzRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vz		,&garray_Vz);										CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Vqual 	,&garray_Vqual);									CHKERRQ(ierr);

	num  = 0;
	ierr =  VecGetSize(gvec_VxRef,&NumOfPoints);										CHKERRQ(ierr);

	for (i=0; i<size_VxRef; i++) {
		if(isnan(garray_VxRef[i])){
			if(isnan(garray_VyRef[i])==0) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VyRef is a number while VxRef is NaN !! ==> CHECK input data\n");
			if(isnan(garray_VzRef[i])==0) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VzRef is a number while VxRef is NaN !! ==> CHECK input data\n");

			// set NAN-components to zero so that Petsc vector routines work fine
			// also set the corresponding grid velocity to zero so that total error isn't effected
			garray_VxRef[i] 	= 0.0;
			garray_VyRef[i] 	= 0.0;
			garray_VzRef[i] 	= 0.0;
			garray_Vx[i] 		= 0.0;
			garray_Vy[i] 		= 0.0;
			garray_Vz[i] 		= 0.0;
			garray_Vqual[i] 	= 0.0;

			// if velocity at this node is NAN, also the total number of points needs to be corrected.
			num = num + 1;
		}
		else{
			if(isnan(garray_VyRef[i])) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VyRef is NaN while VxRef is not! ==> CHECK input data\n");
			if(isnan(garray_VzRef[i])) PetscPrintf(PETSC_COMM_SELF,"--- ERROR: VzRef is NaN while VxRef is not! ==> CHECK input data\n");
		}
	}

	ierr = MPI_Allreduce(&num,&global_num,1,MPIU_INT,MPI_SUM,PLANE_COMM);				CHKERRQ(ierr);
	NumOfPoints = NumOfPoints - global_num;
	user->Optimisation.NSurfVel = (PetscScalar) NumOfPoints*3;

	ierr = VecRestoreArray(gvec_Vqual	,&garray_Vqual);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_Vx		,&garray_Vx);									CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_VxRef	,&garray_VxRef);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_Vy		,&garray_Vy);									CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_VyRef	,&garray_VyRef);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_Vz		,&garray_Vz);									CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_VzRef	,&garray_VzRef);								CHKERRQ(ierr);

// 3.) --- Normalization --

	// global max(|VxRef|), global max(|VyRef|)
//	ierr = VecAbsMax(gvec_VxRef,&max_VxRefAbs);											CHKERRQ(ierr);
//	ierr = VecAbsMax(gvec_VyRef,&max_VyRefAbs);											CHKERRQ(ierr);
//	ierr = VecAbsMax(gvec_VzRef,&max_VzRefAbs);											CHKERRQ(ierr);

	// Normalize with the maximum of the reference "true" signal
//	ierr = VecScale(gvec_VxRef,1/max_VxRefAbs);											CHKERRQ(ierr);
//	ierr = VecScale(gvec_Vx   ,1/max_VxRefAbs);											CHKERRQ(ierr);
//	ierr = VecScale(gvec_VyRef,1/max_VyRefAbs);											CHKERRQ(ierr);
//	ierr = VecScale(gvec_Vy   ,1/max_VyRefAbs);											CHKERRQ(ierr);
//	ierr = VecScale(gvec_VzRef,1/max_VzRefAbs);											CHKERRQ(ierr);
//	ierr = VecScale(gvec_Vz   ,1/max_VzRefAbs);											CHKERRQ(ierr);
	ierr = VecScale(gvec_VxRef,1/(user->SurfVelField.VxStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_Vx   ,1/(user->SurfVelField.VxStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_VyRef,1/(user->SurfVelField.VyStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_Vy   ,1/(user->SurfVelField.VyStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_VzRef,1/(user->SurfVelField.VzStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_Vz   ,1/(user->SurfVelField.VzStdDev));						CHKERRQ(ierr);


// 4.) --- calculate error ---

	//  lvec_VxRef = (VxRef-Vx),
	ierr = VecAXPBY(gvec_VxRef,1.0,-1.0	,gvec_Vx);										CHKERRQ(ierr);
	ierr = VecAXPBY(gvec_VyRef,1.0,-1.0	,gvec_Vy);										CHKERRQ(ierr);
	ierr = VecAXPBY(gvec_VzRef,1.0,-1.0	,gvec_Vz);										CHKERRQ(ierr);

	// l1 norm of lvec_VRef
	ierr = VecNorm(gvec_VxRef,NORM_1,&l1_Vx);									        CHKERRQ(ierr);
	ierr = VecNorm(gvec_VyRef,NORM_1,&l1_Vy);									        CHKERRQ(ierr);
	ierr = VecNorm(gvec_VzRef,NORM_1,&l1_Vz);									        CHKERRQ(ierr);

	user->Optimisation.SumAbsVel = l1_Vx+l1_Vy+l1_Vz;

	//lvec_VxRef = vx.^2,
	ierr = VecPointwiseMult(gvec_VxRef,gvec_VxRef,gvec_VxRef);							CHKERRQ(ierr);
	ierr = VecPointwiseMult(gvec_VyRef,gvec_VyRef,gvec_VyRef);							CHKERRQ(ierr);
	ierr = VecPointwiseMult(gvec_VzRef,gvec_VzRef,gvec_VzRef);							CHKERRQ(ierr);

	// lvec_VxRef = vx.^2 * quality,
	ierr = VecPointwiseMult(gvec_VxRef,gvec_VxRef,gvec_Vqual);							CHKERRQ(ierr);
	ierr = VecPointwiseMult(gvec_VyRef,gvec_VyRef,gvec_Vqual);							CHKERRQ(ierr);
	ierr = VecPointwiseMult(gvec_VzRef,gvec_VzRef,gvec_Vqual);							CHKERRQ(ierr);

	// lvec_VxRef = sum(vx.^2 *quality),
	ierr = VecSum(gvec_VxRef,&err_Vx);													CHKERRQ(ierr);
	ierr = VecSum(gvec_VyRef,&err_Vy);													CHKERRQ(ierr);
	ierr = VecSum(gvec_VzRef,&err_Vz);													CHKERRQ(ierr);

	user->Optimisation.SumSqrsVel = err_Vx+err_Vy+err_Vz;


	// RMS = sqrt(1/N * sum(x_i.^2)  )
	//*global_VelMisfit = sqrt( (err_Vx + err_Vy + err_Vz)/(3*NumOfPoints) );

	ierr = PetscPrintf(PLANE_COMM,"#     Global err_Vx: %g \n",err_Vx);					CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     Global err_Vy: %g \n",err_Vy);					CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     Global err_Vz: %g \n",err_Vz);					CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     NumberofPoints: %lld \n",(LLD)NumOfPoints);	CHKERRQ(ierr);

// 5.) --- clean reference vectors ---
	ierr = VecDestroy(&gvec_VxRef); 													CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_VyRef); 													CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_VzRef); 													CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_Vqual); 													CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
