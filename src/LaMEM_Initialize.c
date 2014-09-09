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

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "LaMEM_Initialize.h"
#include "Utils.h"
#include "Mesh.h"
#include "Assembly_FDSTAG.h"

/*==============================================================================================================================
		Driver routine to initialize LaMEM - checks whether we use FDSTAG or FEM, and calls the relevant subroutines
*/
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_StokesSolver"
PetscErrorCode LaMEM_Initialize_StokesSolver( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C)
{
	PetscErrorCode 	ierr;

	if (vpt_element_type == DAVP_FDSTAG) {
		/* Finite difference method */
		 ierr = LaMEM_Initialize_StokesSolver_FDSTAG(user, vpt_element_type, C); CHKERRQ(ierr);
	} else if (vpt_element_type == DAVP_Q1P0) {
		ierr = LaMEM_Initialize_StokesSolver_FEM(user, vpt_element_type, C); CHKERRQ(ierr);
	} else if (vpt_element_type == DAVP_Q1Q1) {
		ierr = LaMEM_Initialize_StokesSolver_FEM(user, vpt_element_type, C); CHKERRQ(ierr);
	} else if ( (vpt_element_type == DAVP_Q2PM1L) || (vpt_element_type == DAVP_Q2PM1G) ) {
		ierr = LaMEM_Initialize_StokesSolver_FEM(user, vpt_element_type, C); CHKERRQ(ierr);
	} else {
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Unknown element type provided; vpt_element_type=%lld",(LLD)vpt_element_type);
	}
	/* Finite element method
	 * NOTE: At a later stage, we should split this routine into more subroutines, depending on
	 * the type of element that is employed.
	 */
	PetscFunctionReturn(0);
}
//==============================================================================================================================
// Initialize matrices required for the Stokes solver and for defining Material properties, in case of using FEM
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_StokesSolver_FEM"
PetscErrorCode LaMEM_Initialize_StokesSolver_FEM( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C)
{
	PetscMPIInt     size;
	PetscErrorCode 	ierr;
	PetscInt 		num_cpu_x, num_cpu_y, num_cpu_z;
	PetscInt 		nel_x=0, nel_y=0, nel_z=0, cpu_x, cpu_y, cpu_z;
	PetscInt		nel_x_local, nel_y_local, nel_z_local, mx, my, mz;
	PetscInt		nnode_x, nnode_y, nnode_z;
	PetscInt		vel_sdof, pres_sdof, fac;
	PetscInt		nvel_rowEnd, 	nvel_rowStart, nvel_rows_Local;
	PetscInt		npres_rowStart, npres_rowEnd, npres_rows_Local;
	PetscInt		*lx, *ly, *lz, dof;
	PetscInt        ndof_mat,i;
	const PetscInt 	 *llx, *lly, *llz;
	PetscInt		*llx_q, *lly_q, *llz_q;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	// Set number of elements in every direction
	if ( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) )
	{	nel_x =	 user->nnode_x-1;
		nel_y =	 user->nnode_y-1;
		nel_z =  user->nnode_z-1;
	}
	else if ((vpt_element_type == DAVP_Q2PM1L) || (vpt_element_type == DAVP_Q2PM1G))
	{	// check number of nodes is sane for Q2 elements
		if(((user->nnode_x-1))%2)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Number of nodes specified in x-dir is incompatible with Q2 elements" );
		if(((user->nnode_y-1))%2)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Number of nodes specified in y-dir is incompatible with Q2 elements" );
		if(((user->nnode_z-1))%2)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Number of nodes specified in z-dir is incompatible with Q2 elements" );
		nel_x = (user->nnode_x-1)/2; 
		nel_y = (user->nnode_y-1)/2; 
		nel_z = (user->nnode_z-1)/2;
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Unknown element type!");
	}

	// This DMDA is used to distribute the elements over the processors
	// NOTE: This is a base DMDA in LaMEM (uses PETSC_DECIDE to compute partitioning)
	// We typically let petsc decide how to distribute the grid. If you do want to override this, do it
	// by adding e.g. -da_proc_z 1 to the command-line (and NOT the PETSC default -da_processors_z 1 !!)
	num_cpu_x = PETSC_DECIDE;
	num_cpu_y = PETSC_DECIDE;
	num_cpu_z = PETSC_DECIDE;
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_x",		&num_cpu_x	, PETSC_NULL);  		//	fix # of processors in x-direction
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_y",		&num_cpu_y	, PETSC_NULL);  		//	fix # of processors in y-direction
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_z",		&num_cpu_z	, PETSC_NULL);  		//	fix # of processors in z-direction
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x,nel_y,nel_z,num_cpu_x,num_cpu_y,num_cpu_z,1,0,0,0,0,&user->DA_Processors); CHKERRQ(ierr);

	// Compute number of elements on the local processor
	ierr = DMDAGetInfo(user->DA_Processors,0,0,0,0, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); 			CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,0,0,0,&nel_x_local,&nel_y_local,&nel_z_local);	CHKERRQ(ierr);

	// compute number "degrees-of-freedom" per "node" of material DMDA
	// product of number material parameters and number of quadrature points
	ndof_mat = NumMaterialPropsElem*(C->ngp_vel);

	// The material properties should only be created locally
	//=====================================================================================
	// Why do we need DMDA anyway? A local array of material data structure would be better
	//=====================================================================================
	ierr = DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x_local,nel_y_local,nel_z_local,1,1,1,ndof_mat,0,0,0,0,&user->DA_Materials); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Materials,0,"Materials"); CHKERRQ(ierr);

	// global vector with local DMDA, that is indeed interesting
	ierr = DMCreateGlobalVector(user->DA_Materials, &user->Materials); CHKERRQ(ierr);
	ierr = VecSet(user->Materials, 0.0); CHKERRQ(ierr);

	// Compute the precise number of quadrature points in every direction
	ierr = DMDAGetOwnershipRanges(user->DA_Processors,&llx, &lly, &llz); CHKERRQ(ierr);
	ierr = makeIntArray(&llx_q, llx, cpu_x); CHKERRQ(ierr);
	ierr = makeIntArray(&lly_q, lly, cpu_y); CHKERRQ(ierr);
	ierr = makeIntArray(&llz_q, llz, cpu_z); CHKERRQ(ierr);

	for (i=0; i<cpu_x; i++){
		llx_q[i] = llx_q[i]*C->nintp_1D;
	}
	for (i=0; i<cpu_y; i++){
		lly_q[i] = lly_q[i]*C->nintp_1D;
	}
	for (i=0; i<cpu_z; i++){
		llz_q[i] = llz_q[i]*C->nintp_1D;
	}

	// Create the quadrature-point DMDA - this is solely created to facilitate saving output data into Paraview VTK format
	if (user->VTKOutputFiles==1){

		ierr = 	DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x*C->nintp_1D,nel_y*C->nintp_1D,nel_z*C->nintp_1D,cpu_x,cpu_y,cpu_z,1,1, llx_q, lly_q, llz_q,&user->DA_Quadrature); CHKERRQ(ierr);
		ierr = 	DMDASetUniformCoordinates(user->DA_Quadrature,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L, user->z_bot, user->z_bot + user->H);	CHKERRQ(ierr);

	}
	PetscFree(llx_q);
	PetscFree(lly_q);
	PetscFree(llz_q);


	ierr = DMDAGetInfo(user->DA_Processors, 0, &mx, &my, &mz, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Fine Grid parameters, nel=[%lld,%lld,%lld], cpu=[%lld,%lld,%lld] \n",(LLD)nel_x,(LLD)nel_y,(LLD)nel_z,(LLD)cpu_x,(LLD)cpu_y,(LLD)cpu_z);

	// Allocate arrays that specify the distribution of nodes on every cpu
	ierr = PetscMalloc((size_t)cpu_x*sizeof(PetscInt), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)cpu_y*sizeof(PetscInt), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)cpu_z*sizeof(PetscInt), &lz); CHKERRQ(ierr);
	ierr = ComputeNodesPerProcessor(user->DA_Processors,lx,ly,lz,vpt_element_type); CHKERRQ(ierr);

	nnode_x     = user->nnode_x;		
	nnode_y     = user->nnode_y;	
	nnode_z     = user->nnode_z;
	user->nel_x = nel_x;
	user->nel_y = nel_y;
	user->nel_z = nel_z;
	
	user->finest_nnode_x = nnode_x;
	user->finest_nnode_y = nnode_y;
	user->finest_nnode_z = nnode_z;
	user->finest_nelx    = nel_x;
	user->finest_nely    = nel_y;
	user->finest_nelz    = nel_z;
	
	// correct for periodicity
	if (user->BC.LeftBound==3){
		nnode_x	 		= nnode_x-1;
		lx[cpu_x-1] 	= lx[cpu_x-1]-1;
	}
	if (user->BC.FrontBound==3){
		nnode_y	 		= nnode_y-1;
		ly[cpu_y-1] 	= ly[cpu_y-1]-1;
	}
	// Create Velocity DMDA
	dof=3;
	if ( (C->ElementType==2) ){
		// Q1P0/Q1Q1 element
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof,1,lx,ly,lz,&user->DA_Vel); CHKERRQ(ierr);
	}
	else {
		// Q2Pm1 element (seems to require one more ghost node then FEM - NOTE: THIS IS SLOOOOWWW; check it out later)
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof,2,lx,ly,lz,&user->DA_Vel); CHKERRQ(ierr);
	}
	ierr = DMDASetFieldName(user->DA_Vel,0,"Vx"); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Vel,1,"Vy"); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Vel,2,"Vz"); CHKERRQ(ierr);

	// set coarsening options to be used with geometric multigrid
	ierr = DMDASetInterpolationType(user->DA_Vel, DMDA_Q1); CHKERRQ(ierr);
	ierr = DMDASetRefinementFactor(user->DA_Vel, user->refinex, user->refiney, user->refinez);  CHKERRQ(ierr);

	// Create Pressure DMDA
	if ((vpt_element_type == DAVP_Q1Q1) ){
		// Continuous pressure
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,	&user->DA_Pres); CHKERRQ(ierr);
	}
	else {
		// Discontinuous pressure
		ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,nel_x,nel_y,nel_z,cpu_x,cpu_y,cpu_z,C->npres,0,0,0,0, &user->DA_Pres); CHKERRQ(ierr);
	}
	ierr = DMDASetFieldName(user->DA_Pres,0,"P"); CHKERRQ(ierr);

	// Create stiffness matrices
	// The nomenclature for the global stiffness matrixes is:
	// | VV_MAT VP_MAT | | V |    |Rhs |
	// |        	   | |	 |  = |    |
	// | PV_MAT PP_MAT | | P |    | 0  |
	
	// Create a DMDA that holds the surface (and later: bottom) topography.
	// Note: the FDSTAG method has something similar, but the surface topography can in that case be inside the model ("sticky air approach")
	if ( (C->ElementType==2) ){
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,0,&user->DA_SurfaceTopography); CHKERRQ(ierr);
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,0,&user->DA_BottomTopography); CHKERRQ(ierr);
	}
	else{
		// note: we have to add 2 more in z-direction for the stencil
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z*2,cpu_x,cpu_y,cpu_z,1,2,lx,ly,0,&user->DA_SurfaceTopography); CHKERRQ(ierr);
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z*2,cpu_x,cpu_y,cpu_z,1,2,lx,ly,0,&user->DA_BottomTopography); CHKERRQ(ierr);
	}
	ierr = DMDASetFieldName(user->DA_SurfaceTopography,0,"SurfaceTopography"); CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(user->DA_SurfaceTopography,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L,	user->z_bot, user->z_bot+ user->H);	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&user->SurfaceTopography); 	CHKERRQ(ierr);
	ierr = VecSet(user->SurfaceTopography,-100 + user->z_bot+ user->H);	CHKERRQ(ierr);


	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vx);	CHKERRQ(ierr);	// will store Vx velocity on free surface
	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vy);	CHKERRQ(ierr);	// will store Vy velocity on free surface
	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vz);	CHKERRQ(ierr);	// will store Vz velocity on free surface


	ierr = DMDASetFieldName(user->DA_BottomTopography,0,"BottomTopography"); CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(user->DA_BottomTopography,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L, user->z_bot,	user->z_bot	+ user->H);	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_BottomTopography, &user->BottomTopography); CHKERRQ(ierr);
	ierr = VecSet(user->BottomTopography, user->z_bot); CHKERRQ(ierr);

	// Create velocity stiffness matrix
	ierr = DMCreateMatrix(user->DA_Vel,	&user->VV_MAT); CHKERRQ(ierr); // get matrix that matches the DA

	// Create PP matrix (pressure mass matrix)
	ierr = DMCreateMatrix(user->DA_Pres, &user->PP_MAT); CHKERRQ(ierr);	// get matrix that matches the DA
	ierr = MatGetSize(user->VV_MAT, &vel_sdof , PETSC_NULL); CHKERRQ(ierr);	    // global # of velocity unknowns
	ierr = MatGetSize(user->PP_MAT, &pres_sdof, PETSC_NULL); CHKERRQ(ierr);	    // global # of pressure unknowns

	// Create VP matrix, with the correct parallel layout (dictated by the VV & PP matrixes).
	ierr = MatGetOwnershipRange(user->VV_MAT,&nvel_rowStart,&nvel_rowEnd); CHKERRQ(ierr);
	nvel_rows_Local = nvel_rowEnd-nvel_rowStart;
	ierr = MatGetOwnershipRange(user->PP_MAT,&npres_rowStart,&npres_rowEnd); CHKERRQ(ierr);
	npres_rows_Local = npres_rowEnd-npres_rowStart;
	ierr = MatCreate(PETSC_COMM_WORLD,&user->VP_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->VP_MAT, nvel_rows_Local, npres_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType(user->VP_MAT,	MATAIJ); CHKERRQ(ierr);
	if (size==1){
		ierr = MatSeqAIJSetPreallocation(user->VP_MAT,C->npres*8,PETSC_NULL); CHKERRQ(ierr); // very coarse preallocation
	}
	else {
		ierr = MatMPIAIJSetPreallocation(user->VP_MAT,C->npres*8,PETSC_NULL,C->npres*8,PETSC_NULL); CHKERRQ(ierr); // very coarse preallocation
	}
	// Create PV matrix, with the correct parallel layout (dictated by the VV & PP matrixes).
	ierr = MatCreate(PETSC_COMM_WORLD, &user->PV_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->PV_MAT, npres_rows_Local, nvel_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType(user->PV_MAT,	MATAIJ); CHKERRQ(ierr);
	// This part is extremely crude, particularly for Q1Q1 elements - please improve in the future (uses too much memory)
	if ((vpt_element_type == DAVP_Q1Q1) ){	fac =	4;	}
	else{									fac =	1;	}
	if (size==1){
		ierr = MatSeqAIJSetPreallocation(user->PV_MAT,C->edof*fac,PETSC_NULL); CHKERRQ(ierr); 				  		// very coarse preallocation
	}
	else{
		ierr = MatMPIAIJSetPreallocation(user->PV_MAT,C->edof*fac,PETSC_NULL,C->edof*fac,PETSC_NULL); CHKERRQ(ierr); // very coarse preallocation
	}
	// Store cpu-distribution to create TEMP DA
	user->cpu_x = cpu_x;
	user->cpu_y = cpu_y;
	user->cpu_z = cpu_z;
	

	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	// Initialize preconditioning matrix
	ierr = MatDuplicate(user->PP_MAT, MAT_COPY_VALUES, &user->approx_S ); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//==============================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_Preallocation_A11FDSTAG"
PetscErrorCode FDSTAG_Preallocation_A11FDSTAG(DM daproc,DM da,DM da_pres,Mat A, PetscInt nvel_rowStart, PetscInt nvel_rowEnd)
{
	PetscErrorCode 	ierr;
	PetscInt        vx,vy,vz;
	MatStencil      row,row_vx,row_vy,row_vz, col,col_vx[15],col_vy[15],col_vz[15];
	PetscInt		n,i,j,k,xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt        xm,ym,zm,xs,ys,zs;
	PetscInt		iel_x, iel_y, iel_z, nnode_x,nnode_y,nnode_z;
	PetscInt        rowidx,vgidx[19];
	PetscInt        start_row,end_row,start_col,end_col,local_rowidx;
	PetscInt        *d_nnz,*od_nnz;

	PetscFunctionBegin;

	ierr = DMDAGetCorners(daproc,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);  // # of nodes 	 in all directions

	start_row = start_col = nvel_rowStart;
	end_row   = end_col   = nvel_rowEnd;

	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&d_nnz);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&od_nnz);CHKERRQ(ierr);
	for (i=0; i<(end_row-start_row); i++) {
		d_nnz[i] = od_nnz[i] = 0;
	}

	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables

				// x-momentum
				vx = 0;
				if ((i>0) && (i<(nnode_x-1)) ){
					row_vx.i = i;   row_vx.j = j; row_vx.k = k; row_vx.c = 0;	// 1th force balance equation is located @ Vx(i) nodes

					// [Sxx(i)-Sxx(i-1)]/dx
					col.i = 	i+1; col.j = j; col.k = k; col.c = 0;				col_vx[vx++] = col;			//  Vx(i+1,j,k)
					col.i = 	i  ; col.j = j; col.k = k; col.c = 0;				col_vx[vx++] = col;			//  Vx(i  ,j,k)
					col.i = 	i-1; col.j = j; col.k = k; col.c = 0;				col_vx[vx++] = col;;			//  Vx(i-1,j,k)

					// [Sxy(j)-Sxy(j-1)]/dy
					// Sxy(j)/dy
					col.i = 	i;   col.j = j+1; col.k = k; col.c = 0;			col_vx[vx++] = col;			//  Vx(i,j+1)
				//	col.i = 	i;   col.j = j;   col.k = k; col.c = 0;			col_vx[4] = col; v_vv[4] = v;			//  Vx(i,j  )
					col.i =   i  ; col.j = j+1; col.k = k; col.c = 1;			col_vx[vx++] = col;			//  Vy(i  ,j+1)
					col.i =	  i-1; col.j = j+1; col.k = k; col.c = 1;			col_vx[vx++] = col;			//  Vy(i-1,j+1)

					// -Sxy(j-1)/dy
				//	col.i = 	i; col.j = j;   col.k = k; col.c = 0;			col_vx[7] = col; v_vv[7] = v;			//  Vx(i,j  )
					col.i = 	i;   col.j = j-1; col.k = k; col.c = 0;			col_vx[vx++] = col;			//  Vx(i,j-1)
					col.i =   i  ; col.j = j  ; col.k = k; col.c = 1;			col_vx[vx++] = col;			//  Vy(i  ,j  )
					col.i =   i-1; col.j = j  ; col.k = k; col.c = 1;			col_vx[vx++]= col;			//  Vy(i-1,j  )

					// [Sxz(k)-Sxz(k-1)]/dz
					// Sxz(k)/dz
					col.i = 	i;   col.j = j; col.k = k+1; col.c = 0;			col_vx[vx++] = col;	//  Vx(i,k+1)
				//	col.i = 	i;   col.j = j; col.k = k;   col.c = 0;			col_vx[12] = col; v_vv[12] = v;	//  Vx(i,k  )
					col.i =   i  ; col.j = j; col.k = k+1; col.c = 2;			col_vx[vx++] = col;	//  Vz(i  ,k+1)
					col.i =	  i-1; col.j = j; col.k = k+1; col.c = 2;			col_vx[vx++] = col;	//  Vz(i-1,k+1)

					// -Sxz(k-1)/dz
				//	col.i = 	i;   col.j = j; col.k = k;   col.c = 0;			col_vx[15] = col; v_vv[15] = v;	//  Vx(i,k  )
					col.i = 	i;   col.j = j; col.k = k-1; col.c = 0;			col_vx[vx++] = col;	//  Vx(i,k-1)
					col.i =   i  ; col.j = j; col.k = k;   col.c = 2;			col_vx[vx++] = col;	//  Vz(i  ,k  )
					col.i =   i-1; col.j = j; col.k = k;   col.c = 2;			col_vx[vx++] = col;	//  Vz(i-1,k  )

					if(vx!=15) {
						SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"nnz count for vx is wrong");
					}

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row_vx.i,row_vx.j,row_vx.k,row_vx.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,15,col_vx,0,PETSC_NULL,vgidx,PETSC_NULL);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<15; n++) {
						if (vgidx[n]<start_col) {
							od_nnz[local_rowidx]++;
						} else if (vgidx[n]>end_col) {
							od_nnz[local_rowidx]++;
						} else {
							d_nnz[local_rowidx]++;
						}
					}

				}

				// y-momentum
				vy = 0;
				if ((j>0) && (j<(nnode_y-1))){
					row_vy.i = i;   row_vy.j = j; row_vy.k = k; row_vy.c = 1;	// 2nd force balance equation is located @ Vy(i,j,k) nodes

					// [Syy(j)-Syy(j-1)]/dy
					col.i = 	i; col.j = j+1; col.k = k; col.c = 1;	        col_vy[vy++] = col;			//  Vy(i,j+1,k)
					col.i = 	i; col.j = j  ; col.k = k; col.c = 1;	        col_vy[vy++] = col;			//  Vy(i,j  ,k)
					col.i = 	i; col.j = j-1; col.k = k; col.c = 1;	        col_vy[vy++] = col;			//  Vy(i,j-1,k)

					// [Sxy(i+1)-Sxy(i)]/dx
					// Sxy(i+1)/dx
					col.i =   i+1; col.j = j;   col.k = k; col.c = 1;	     col_vy[vy++] = col;			//  Vy(i+1,j)
					//col.i = 	i;   col.j = j;   col.k = k; col.c = 1;	   col_vy[vy++] = col;			//  Vy(i  ,j  )
					col.i =   i+1; col.j = j-1; col.k = k; col.c = 0;   	 col_vy[vy++] = col;			//  Vx(i+1,j-1)
					col.i =   i+1; col.j = j  ; col.k = k; col.c = 0;	     col_vy[vy++] = col;			//  Vx(i+1,j  )

					// -Sxy(i)/dx
					//col.i = 	i;   col.j = j;   col.k = k; col.c = 1;	     col_vy[vy++] = col;			//  Vy(i  ,j  )
					col.i =   i-1; col.j = j;   col.k = k; col.c = 1;	     col_vy[vy++] = col;			//  Vy(i-1,j)
					col.i =	  i  ; col.j = j-1; col.k = k; col.c = 0;	     col_vy[vy++] = col;			//  Vx(i  ,j-1)
					col.i =   i  ; col.j = j  ; col.k = k; col.c = 0;	     col_vy[vy++]= col;			//  Vx(i  ,j  )

					// [Syz(k+1)-Syz(k)]/dz
					// Syz(k+1)/dz
					col.i = i; col.j = j;   col.k = k+1; col.c = 1;	       col_vy[vy++] = col;	//  Vy(i,j,k+1)
					//col.i = i; col.j = j;   col.k = k;   col.c = 1;	       col_vy[vy++] = col;	//  Vy(i,j,k  )
					col.i = i; col.j = j  ; col.k = k+1; col.c = 2;	       col_vy[vy++] = col;	//  Vz(i  ,j,k+1)
					col.i =	i; col.j = j-1; col.k = k+1; col.c = 2;	       col_vy[vy++] = col;	//  Vz(i-1,j,k+1)

					// -Syz(k)/dz
					//col.i = i; col.j = j;   col.k = k;   col.c = 1;	       col_vy[vy++] = col;	//  Vy(i,j,k  )
					col.i = i; col.j = j;   col.k = k-1; col.c = 1;	       col_vy[vy++] = col;	//  Vy(i,j,k-1)
					col.i = i; col.j = j  ; col.k = k;   col.c = 2;	       col_vy[vy++] = col;	//  Vz(i  ,j,k  )
					col.i = i; col.j = j-1; col.k = k;   col.c = 2;	       col_vy[vy++] = col;	//  Vz(i-1,j,k  )

					if(vy!=15) {
						SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"nnz count for vy is wrong");
					}

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row_vy.i,row_vy.j,row_vy.k,row_vy.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,15,col_vy,0,PETSC_NULL,vgidx,PETSC_NULL);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<15; n++) {
						if (vgidx[n]<start_col) {
							od_nnz[local_rowidx]++;
						} else if (vgidx[n]>end_col) {
							od_nnz[local_rowidx]++;
						} else {
							d_nnz[local_rowidx]++;
						}
					}
				}

				// z-momentum
				vz = 0;
				if ( (k>0) && (k<(nnode_z-1)) ){

					row_vz.i = i;   row_vz.j = j; row_vz.k = k; row_vz.c = 2;		// 3rd force balance equation is located @ Vz(i,j,k) nodes

					// [Szz(k)-Szz(k-1)]/dz
					col.i = 	i; col.j = j; col.k = k+1; col.c = 2;	  col_vz[vz++] = col;			//  Vz(i,j,k+1)
					col.i = 	i; col.j = j; col.k = k  ; col.c = 2;	  col_vz[vz++] = col;			//  Vz(i,j,k)
					col.i = 	i; col.j = j; col.k = k-1; col.c = 2;	  col_vz[vz++] = col;			//  Vz(i,j,k-1)

					// [Sxz(i+1)-Sxz(i)]/dx
					// Sxz(i+1)/dx
					col.i = 	i+1; col.j = j; col.k = k;   col.c = 2;	  col_vz[vz++] = col;			//  Vz(i+1,j,k)
					//col.i = 	i; 	 col.j = j; col.k = k;   col.c = 2;	  col_vz[vz++] = col;			//  Vz(i  ,j,k)
					col.i =   i+1; col.j = j; col.k = k  ; col.c = 0;	  col_vz[vz++] = col;			//  Vx(i+1,j, k  )
					col.i =	  i+1; col.j = j; col.k = k-1; col.c = 0;	  col_vz[vz++] = col;			//  Vx(i+1,j, k-1)

					// -Sxz(i)/dx
					//col.i = 	i; 	 col.j = j; col.k = k;   col.c = 2;	  col_vz[vz++] = col;			//  Vz(i  ,j,k)
					col.i = 	i-1; col.j = j; col.k = k;   col.c = 2;	  col_vz[vz++] = col;			//  Vz(i-1,j,k)
					col.i =   i  ; col.j = j; col.k = k;   col.c = 0;	  col_vz[vz++] = col;			//  Vx(i  ,j, k  )
					col.i =   i  ; col.j = j; col.k = k-1; col.c = 0;	  col_vz[vz++]= col;			//  Vx(i  ,j, k-1)

					// [Syz(j+1)-Syz(j)]/dy
					// Syz(j+1)/dy
					col.i = i; col.j = j+1; col.k = k;   col.c = 2;	  col_vz[vz++] = col;	//  Vz(i,j+1,k)
					//col.i = i; col.j = j  ; col.k = k;   col.c = 2;	  col_vz[vz++] = col;	//  Vz(i,j  ,k)
					col.i = i; col.j = j+1; col.k = k  ; col.c = 1;   col_vz[vz++] = col;	//  Vy(i,j+1,k  )
					col.i =	i; col.j = j+1; col.k = k-1; col.c = 1;	  col_vz[vz++] = col;	//  Vy(i,j+1,k-1)

					// -Syz(j)/dy
					//col.i = i; col.j = j  ; col.k = k;   col.c = 2;	  col_vz[vz++] = col;	//  Vz(i,j  ,k)
					col.i = i; col.j = j-1; col.k = k;   col.c = 2;	  col_vz[vz++] = col;	//  Vz(i,j-1,k)
					col.i = i; col.j = j;   col.k = k  ; col.c = 1;	  col_vz[vz++] = col;	//  Vy(i,j  ,k  )
					col.i = i; col.j = j;   col.k = k-1; col.c = 1;	  col_vz[vz++] = col;	//  Vy(i,j  ,k-1)


					if(vz!=15) {
						SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"nnz count for vz is wrong");
					}

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row_vz.i,row_vz.j,row_vz.k,row_vz.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,15,col_vz,0,PETSC_NULL,vgidx,PETSC_NULL);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<15; n++) {
						if (vgidx[n]<start_col) {
							od_nnz[local_rowidx]++;
						} else if (vgidx[n]>end_col) {
							od_nnz[local_rowidx]++;
						} else {
							d_nnz[local_rowidx]++;
						}
					}

				}


			}
		}
	}


	// Set Dirichlet BCs for velocity normal to the side boundaries
	// NOTE: THERE SHOULD BE A PETSC ROUTINE TO DO THIS!

	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for (i=xs; i<xs+xm; i++){

				if (((i==0) || (i==nnode_x-1)) && (j<nnode_y-1) && (k<nnode_z-1)){
					/* Left or right boundary [Vx=constant] */
					row.i = i;   row.j = j; row.k = k; row.c = 0;
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					local_rowidx = rowidx - start_row;
					d_nnz[local_rowidx] = 1;
				}
				if ( ((j==0) || (j==nnode_y-1)) && (i<nnode_x-1) && (k<nnode_z-1)){
					/* Front or back boundary [Vy=constant]*/
					row.i = i;   row.j = j; row.k = k; row.c = 1;
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					local_rowidx = rowidx - start_row;
					d_nnz[local_rowidx] = 1;
				}
				if (((k==0) || (k==nnode_z-1)) && (i<nnode_x-1) && (j<nnode_y-1) ){
					/* Lower or upper boundary [Vz=constant]*/
					row.i = i;   row.j = j; row.k = k; row.c = 2;
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					local_rowidx = rowidx - start_row;
					d_nnz[local_rowidx] = 1;
				}

			}
		}
	}



	// The FDSTAG formulation, in combination with the DA, results in some velocity points that are not used in the discretization.
	// Yet, the matrix should not contain any rows that are zero.
	// Therefore, we have to set the diagonal of the matrix to 1 for these velocity points
	ierr = DMDAGetInfo(da,	0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);  // # of nodes 	 in all directions
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	/* Set additional Vx points */
	for (j=ys; j<ys+ym; j++){
		for	(i=xs; i<xs+xm; i++){
			if (zs+zm==nnode_z){  // the local processor borders this boundary


				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row.i = i;   row.j = j; row.k = k; row.c = 0;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}


	for (k=zs; k<zs+zm; k++){
		for(i=xs; i<xs+xm; i++){
			if (ys+ym==nnode_y){
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row.i = i;   row.j = j; row.k = k; row.c = 0;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}

	/* Set additional Vy points */
	for (j=ys; j<ys+ym; j++){
		for(i=xs; i<xs+xm; i++){
			if (zs+zm==nnode_z){
				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row.i = i;   row.j = j; row.k = k; row.c = 1;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}
	for (k=zs; k<zs+zm; k++){
		for(j=ys; j<ys+ym; j++){
			if (xs+xm==nnode_x){
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row.i = i;   row.j = j; row.k = k; row.c = 1;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}

	/* Set additional Vz points */
	for (k=zs; k<zs+zm; k++){
		for(i=xs; i<xs+xm; i++){
			if (ys+ym==nnode_y){
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row.i = i;   row.j = j; row.k = k; row.c = 2;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}
	for (k=zs; k<zs+zm; k++){
		for(j=ys; j<ys+ym; j++){
			if (xs+xm==nnode_x){
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row.i = i;   row.j = j; row.k = k; row.c = 2;
				ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
				local_rowidx = rowidx - start_row;
				d_nnz[local_rowidx] = 1;
			}
		}
	}

	for (i=0; i<(end_row-start_row); i++) {
		if ( (d_nnz[i]+od_nnz[i]) == 0 ) {
			printf("nnz[%lld] %lld , o_nnz[%lld] %lld , sum %lld \n", (LLD)(i+start_row), (LLD)d_nnz[i], (LLD)(i+start_row), (LLD)od_nnz[i], (LLD)(d_nnz[i]+od_nnz[i]) );
			SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"total nnz <= 0  : i=%lld nnz = %lld, onnz = %lld \n", (LLD)i,(LLD)d_nnz[i],(LLD)od_nnz[i] );
		}
	}

	ierr = MatSeqAIJSetPreallocation(A,0,d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A,0,d_nnz,0,od_nnz); CHKERRQ(ierr);
	// set memory option
	ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);

	ierr = PetscFree(d_nnz);CHKERRQ(ierr);
	ierr = PetscFree(od_nnz);CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_Preallocation_A12FDSTAG"
PetscErrorCode FDSTAG_Preallocation_A12FDSTAG(DM daproc,DM da,DM da_pres,Mat A,
		PetscInt nvel_rowStart,  PetscInt nvel_rowEnd,
		PetscInt npres_rowStart, PetscInt npres_rowEnd)
{
	PetscErrorCode 	ierr;
	MatStencil      row,col,col_vp[2];
	PetscInt		n,i,j,k,xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt        xm,ym,zm,xs,ys,zs;
	PetscInt		iel_x, iel_y, iel_z, nnode_x,nnode_y,nnode_z;
	PetscInt        rowidx,pgidx[2];
	PetscInt        start_row,end_row,start_col,end_col,local_rowidx;
	PetscInt        *d_nnz,*od_nnz;

	PetscFunctionBegin;

	ierr = DMDAGetCorners(daproc,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,      			      &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);  // # of nodes 	 in all directions

	start_row = nvel_rowStart;
	end_row   = nvel_rowEnd;
	start_col = npres_rowStart;
	end_col   = npres_rowEnd;

	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&d_nnz);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&od_nnz);CHKERRQ(ierr);
	for (i=0; i<(end_row-start_row); i++) {
		d_nnz[i] = od_nnz[i] = 0;
	}

	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables

				// gradient - dp/dx
				if ((i>0) && (i<(nnode_x-1)) ){
					row.i = i;   row.j = j; row.k = k; row.c = 0;		// 1th force balance equation is located @ Vx(i) nodes
					// -[ P(i) - P(i-1) ]/dx_P_west		- Pressure equations
					col.i = i-1; col.j = j; col.k = k; col.c = 0;  col_vp[0] = col; 	//  P(i-1)/dx
					col.i = i  ; col.j = j; col.k = k; col.c = 0;  col_vp[1] = col;	// -P(i  )/dx

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,0,PETSC_NULL,2,col_vp,PETSC_NULL,pgidx);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<2; n++) {
						if        (pgidx[n]<start_col) {  od_nnz[local_rowidx]++;
						} else if (pgidx[n]>end_col)   {  od_nnz[local_rowidx]++;
						} else {                          d_nnz[local_rowidx]++;
						}
					}
				}

				// gradient - dp/dy
				if ((j>0) && (j<(nnode_y-1))){
					row.i = i;   row.j = j; row.k = k; row.c = 1;	 // 2nd force balance equation is located @ Vy(i,j,k) nodes
					// -[ P(j) - P(j-1) ]/dy_P_south		- Pressure equations
					col.i = i; col.j = j-1; col.k = k; col.c = 0;  col_vp[0] = col;  //  P(j-1)/dy
					col.i = i; col.j = j  ; col.k = k; col.c = 0;  col_vp[1] = col;	// -P(j  )/dy

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,0,PETSC_NULL,2,col_vp,PETSC_NULL,pgidx);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<2; n++) {
						if        (pgidx[n]<start_col) {  od_nnz[local_rowidx]++;
						} else if (pgidx[n]>end_col)   {  od_nnz[local_rowidx]++;
						} else {                          d_nnz[local_rowidx]++;
						}
					}

				}

				// gradient - dp/dz
				if ( (k>0) && (k<(nnode_z-1)) ){
					row.i = i;   row.j = j; row.k = k; row.c = 2;	 // 3rd force balance equation is located @ Vz(i,j,k) nodes
					// -[ P(k) - P(k-1) ]/dz_P_bottom		- Pressure equations
					col.i = i; col.j = j; col.k = k-1; col.c = 0;  col_vp[0] = col;	//  P(k-1)/dz
					col.i = i; col.j = j; col.k = k  ; col.c = 0;  col_vp[1] = col;	// -P(k  )/dz

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da,row.i,row.j,row.k,row.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,0,PETSC_NULL,2,col_vp,PETSC_NULL,pgidx);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<2; n++) {
						if        (pgidx[n]<start_col) {  od_nnz[local_rowidx]++;
						} else if (pgidx[n]>end_col)   {  od_nnz[local_rowidx]++;
						} else {                          d_nnz[local_rowidx]++;
						}
					}

				}

			}
		}
	}

	ierr = MatSeqAIJSetPreallocation(A,0,d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A,0,d_nnz,0,od_nnz); CHKERRQ(ierr);
	// set memory option
	ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);


	ierr = PetscFree(d_nnz);CHKERRQ(ierr);
	ierr = PetscFree(od_nnz);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_Preallocation_A21FDSTAG"
PetscErrorCode FDSTAG_Preallocation_A21FDSTAG(DM daproc,DM da,DM da_pres,Mat A,
		PetscInt nvel_rowStart,  PetscInt nvel_rowEnd,
		PetscInt npres_rowStart, PetscInt npres_rowEnd)
{
	PetscErrorCode 	ierr;
	PetscInt        peqn;
	MatStencil      col,row_p,col_v[6];
	PetscInt		n,i,j,k,xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt        xm,ym,zm,xs,ys,zs;
	PetscInt		iel_x, iel_y, iel_z, nnode_x,nnode_y,nnode_z;
	PetscInt        rowidx,vgidx[19];
	PetscInt        start_row,end_row,start_col,end_col,local_rowidx;
	PetscInt        *d_nnz,*od_nnz;

	PetscFunctionBegin;

	ierr = DMDAGetCorners(daproc,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,      			      &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);  // # of nodes 	 in all directions

	start_row = npres_rowStart;
	end_row   = npres_rowEnd;
	start_col = nvel_rowStart;
	end_col   = nvel_rowEnd;

	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&d_nnz);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(end_row-start_row),&od_nnz);CHKERRQ(ierr);
	for (i=0; i<(end_row-start_row); i++) {
		d_nnz[i] = od_nnz[i] = 0;
	}


	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables

				// incompressibility
				peqn = 0;
				if ( (i<(nnode_x-1)) && (j<(nnode_y-1)) && (k<(nnode_z-1))  ){
					row_p.i = i;   row_p.j = j; row_p.k = k; row_p.c = 0;	// 1th force balance equation is located @ Vx(i) nodes

					// (Vx(i+1)-Vx(i))/dx
					col.i = 	i+1; col.j = j; col.k = k; col.c = 0;	 col_v[peqn++] = col;	// 	Vx(i+1)/dx
					col.i = 	i  ; col.j = j; col.k = k; col.c = 0;	 col_v[peqn++] = col; 	// 	-Vx(i)/dx

					// (Vy(j+1)-Vy(j))/dy
					col.i = 	i; col.j = j+1; col.k = k; col.c = 1;	 col_v[peqn++] = col;	//   Vy(i+1)/dy
					col.i = 	i; col.j = j  ; col.k = k; col.c = 1;	 col_v[peqn++] = col;	//  -Vy(i)/dy


					// (Vz(k+1)-Vz(k))/dz
					col.i = 	i; col.j = j; col.k = k+1; col.c = 2;	 col_v[peqn++] = col;	//   Vz(i+1)/dz
					col.i = 	i; col.j = j; col.k = k  ; col.c = 2;	 col_v[peqn++] = col;	//  -Vz(i)/dz

					if(peqn!=6) {
						SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"nnz count for continuity is wrong");
					}

					/* convert values to global indices */
					ierr = DAGetGlobalIndex(da_pres,row_p.i,row_p.j,row_p.k,row_p.c,&rowidx);CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,6,col_v,0,PETSC_NULL,vgidx,PETSC_NULL);CHKERRQ(ierr);

					/* offset into local indices in nnz[], i.e. convert rowidx to local index  */
					local_rowidx = rowidx - start_row;

					/* count */
					for (n=0; n<6; n++) {
						if (vgidx[n]<start_col) {
							od_nnz[local_rowidx]++;
						} else if (vgidx[n]>end_col) {
							od_nnz[local_rowidx]++;
						} else {
							d_nnz[local_rowidx]++;
						}
					}

				}

			}
		}
	}

	 for (i=0; i<(end_row-start_row); i++) {
		 if ( (d_nnz[i]+od_nnz[i]) == 0 ) {
			 printf("nnz[%lld] %lld , o_nnz[%lld] %lld , sum %lld \n", (LLD)(i+start_row),(LLD)d_nnz[i], (LLD)(i+start_row),(LLD)od_nnz[i], (LLD)(d_nnz[i]+od_nnz[i]) );
			 SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"total nnz <= 0  : i=%lld nnz = %lld, onnz = %lld \n", (LLD)i,(LLD)d_nnz[i],(LLD)od_nnz[i] );
		 }
	 }

	ierr = MatSeqAIJSetPreallocation(A,0,d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A,0,d_nnz,0,od_nnz); CHKERRQ(ierr);
	// set memory option
	ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);

	ierr = PetscFree(d_nnz);CHKERRQ(ierr);
	ierr = PetscFree(od_nnz);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//============================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "SetSinusoidalPerturbation"
// There are some cases in which we set an initial sinusoidal perturbation on the free surface
PetscErrorCode SetSinusoidalPerturbation(PetscScalar SinusoidalFreeSurfaceAmplitude, UserContext *user)
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
//============================================================================================================================
// Initialize matrices required for the Stokes solver and for defining Material properties, in case of using FDSTAG
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_StokesSolver_FDSTAG"
PetscErrorCode LaMEM_Initialize_StokesSolver_FDSTAG( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C)
{
	//==================
	// e c l e c t i c a
	//==================
	PetscErrorCode 	ierr;
	PetscBool	    flg;
	PetscInt		num_cpu_x, num_cpu_y, num_cpu_z;
	PetscInt 		i, nel_x, nel_y, nel_z, cpu_x, cpu_y, cpu_z;
	PetscInt		nel_x_local, nel_y_local, nel_z_local, mx, my, mz;
	PetscInt		nnode_x, nnode_y, nnode_z;
	PetscInt		vel_sdof, pres_sdof;
	PetscInt		nvel_rowEnd, 	nvel_rowStart, nvel_rows_Local;
	PetscInt		npres_rowStart, npres_rowEnd, npres_rows_Local;
	PetscInt		*lx, *ly, *lz, dof,ii;
	PetscScalar     SinusoidalFreeSurfaceAmplitude;
	PetscScalar 	    z_surface;
	Vec             tmp_u,tmp_p;
	const PetscInt 	 *llx, *lly, *llz;
	PetscInt		*llx_q, *lly_q, *llz_q;


	// Define the number of elements in every direction based on the number of nodes
	nel_x = user->nnode_x-1;
	nel_y =	user->nnode_y-1;
	nel_z = user->nnode_z-1;

	// This DMDA is used to distribute the elements over the processors
	// We typically let petsc decide how to distribute the grid. If you do want to override this, do it
	// by adding e.g. -da_proc_z 1 to the command-line (and NOT the PETSC default -da_processors_z 1 !!)
	num_cpu_x = PETSC_DECIDE;
	num_cpu_y = PETSC_DECIDE;
	num_cpu_z = PETSC_DECIDE;
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_x",		&num_cpu_x	, PETSC_NULL);  		//	fix # of processors in x-direction
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_y",		&num_cpu_y	, PETSC_NULL);  		//	fix # of processors in y-direction
	PetscOptionsGetInt(PETSC_NULL ,"-da_proc_z",		&num_cpu_z	, PETSC_NULL);  		//	fix # of processors in z-direction

	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x,nel_y,nel_z,num_cpu_x,num_cpu_y,num_cpu_z,1,0,0,0,0,&user->DA_Processors); CHKERRQ(ierr);

	// Compute number of elements on the local processor
	ierr = DMDAGetInfo(user->DA_Processors,0,0,0,0, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,0,0,0,&nel_x_local,&nel_y_local,&nel_z_local); CHKERRQ(ierr);

	// The material properties should only be created locally
	ierr = DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x_local,nel_y_local,nel_z_local,1,1,1,NumMaterialPropsElem*C->ngp_vel,0,0,0,0,&user->DA_Materials); CHKERRQ(ierr);


	// for FDSTAG formulation, we need to define viscosity at corners and at centers.
	// In addition, the stress tensor should be defined at this location, but density should be defined @ different location.

	ierr = DMDASetFieldName(user->DA_Materials,0,"Materials"); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_Materials,&user->Materials); CHKERRQ(ierr);
	ierr = VecSet(user->Materials, 0.0); CHKERRQ(ierr);

	ierr = DMDAGetOwnershipRanges(user->DA_Processors,&llx, &lly, &llz); CHKERRQ(ierr);
	ierr = makeIntArray(&llx_q, llx, cpu_x); CHKERRQ(ierr);
	ierr = makeIntArray(&lly_q, lly, cpu_y); CHKERRQ(ierr);
	ierr = makeIntArray(&llz_q, llz, cpu_z); CHKERRQ(ierr);

	for (ii=0; ii<cpu_x; ii++){
		llx_q[ii] = llx_q[ii]*C->nintp_1D;
	}
	for (ii=0; ii<cpu_y; ii++){
		lly_q[ii] = lly_q[ii]*C->nintp_1D;
	}
	for (ii=0; ii<cpu_z; ii++){
		llz_q[ii] = llz_q[ii]*C->nintp_1D;
	}

	// Create the quadrature-point DMDA - this is solely created to facilitate saving output data into Paraview VTK format
	if (user->VTKOutputFiles==1){

		ierr = 	DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x*C->nintp_1D,nel_y*C->nintp_1D,nel_z*C->nintp_1D,cpu_x,cpu_y,cpu_z,1,1, llx_q, lly_q, llz_q,&user->DA_Quadrature); CHKERRQ(ierr);
		ierr = 	DMDASetUniformCoordinates(user->DA_Quadrature,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L, user->z_bot, user->z_bot + user->H);	CHKERRQ(ierr);

	}
	PetscFree(llx_q);
	PetscFree(lly_q);
	PetscFree(llz_q);

	ierr = DMDAGetInfo(user->DA_Processors, 0, &mx, &my, &mz, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Fine Grid parameters, nel=[%lld,%lld,%lld], cpu=[%lld,%lld,%lld] \n",(LLD)nel_x,(LLD)nel_y,(LLD)nel_z,(LLD)cpu_x,(LLD)cpu_y,(LLD)cpu_z);

	// Allocate arrays that specify the distribution of nodes on every cpu
	ierr = PetscMalloc((size_t)cpu_x*sizeof(PetscInt), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)cpu_y*sizeof(PetscInt), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)cpu_z*sizeof(PetscInt), &lz); CHKERRQ(ierr);
	ierr = ComputeNodesPerProcessor(user->DA_Processors,lx,ly,lz,vpt_element_type); CHKERRQ(ierr);

	nel_x = user->nnode_x-1;
	nel_y =	user->nnode_y-1;
	nel_z = user->nnode_z-1;

	nnode_x     = user->nnode_x;
	nnode_y     = user->nnode_y;
	nnode_z     = user->nnode_z;
	user->nel_x	= nel_x;
	user->nel_y	= nel_y;
	user->nel_z	= nel_z;

	user->finest_nnode_x = nnode_x;
	user->finest_nnode_y = nnode_y;
	user->finest_nnode_z = nnode_z;
	user->finest_nelx    = nel_x;
	user->finest_nely    = nel_y;
	user->finest_nelz    = nel_z;

	// correct for periodicity
	if (user->BC.LeftBound==3){
		nnode_x	 		= nnode_x-1;
		lx[cpu_x-1] 	= lx[cpu_x-1]-1;
	}
	if (user->BC.FrontBound==3){
		nnode_y	 		= nnode_y-1;
		ly[cpu_y-1] 	= ly[cpu_y-1]-1;
	}

	// Create Velocity DMDA for FDSTAG
	dof=3;

	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof,1,lx,ly,lz,&user->DA_Vel); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Vel,0,"Vx"); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Vel,1,"Vy"); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_Vel,2,"Vz"); CHKERRQ(ierr);

	// Create Pressure DMDA for FDSTAG. As pressure is continuous we need ghost points
	//======================================================
	// *** ??? WOW! continuous pressure? is this FDSTAG? ***
	//======================================================

// a bit more careful, please!
	lx[cpu_x-1]--; ly[cpu_y-1]--; lz[cpu_z-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x, user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nel_x,nel_y,nel_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,	&user->DA_Pres); CHKERRQ(ierr);
	lx[cpu_x-1]++; ly[cpu_y-1]++; lz[cpu_z-1]++;
//	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x, user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nel_x,nel_y,nel_z,cpu_x,cpu_y,cpu_z,1,1,0,0,0,	&user->DA_Pres); CHKERRQ(ierr);

	ierr = DMDASetFieldName(user->DA_Pres,0,"P"); CHKERRQ(ierr);

	//=======================================================================================================
	// for FDSTAG formulation, we need to define viscosity at corners and at centers.
	// In addition, the stress tensor should be defined at the center of the control volume.
	// Yet, density should be defined at a different location.

	// An additional complication is that the FDSTAG requires knowledge of neighboring control volumes,
	// which might be located on different processors. This thus requires creation of a DMDA with BOX-stencil
	// and stencil width 1.
	// Moreover, we are not allowed to create a DMDA with
	// too many DOF per element, as then the global DMDA vector might have more then 2^31 entries.
	// Therefore, we create separate vectors for each field.

	// At best would be local arrays of material data structure
	//=======================================================================================================

	// Create a DMDA that holds the surface topography.
	// Note: for FDSTAG method the surface topography can be inside the model ("sticky air approach")

	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,0,&user->DA_SurfaceTopography); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_SurfaceTopography,0,"SurfaceTopography"); CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(user->DA_SurfaceTopography,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L,	user->z_bot, user->z_bot+ user->H);	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&user->SurfaceTopography); 	CHKERRQ(ierr);
	z_surface = user->z_bot	+ user->H;
	if (user->ErosionParameters.UseInternalFreeSurface == 1){
		// We use an internal free surface
		z_surface = user->ErosionParameters.InitialFreeSurfaceHeight;
	}
	// initialize (internal) surface topography
	ierr = VecSet(user->SurfaceTopography,z_surface); CHKERRQ(ierr);

	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vx);	CHKERRQ(ierr);	// will store Vx velocity on free surface
	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vy);	CHKERRQ(ierr);	// will store Vy velocity on free surface
	ierr = VecDuplicate(user->SurfaceTopography,&user->SurfaceTopography_Vz);	CHKERRQ(ierr);	// will store Vz velocity on free surface


	// set sinusoidal perturbation if necessary
	PetscOptionsGetReal(PETSC_NULL ,"-Initial_SinusoidalFreeSurfaceAmplitude", 	&SinusoidalFreeSurfaceAmplitude, &flg);
	if(flg) ierr = SetSinusoidalPerturbation(SinusoidalFreeSurfaceAmplitude, user); CHKERRQ(ierr);

	// Initialize bottom topography

	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,cpu_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,0,&user->DA_BottomTopography); CHKERRQ(ierr);
	ierr = DMDASetFieldName(user->DA_BottomTopography,0,"BottomTopography"); CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(user->DA_BottomTopography,  user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L,user->z_bot, user->z_bot+ user->H);	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_BottomTopography,&user->BottomTopography);	CHKERRQ(ierr);
	ierr = VecSet(user->BottomTopography,user->z_bot); CHKERRQ(ierr);

	// Center values DMDA
// a bit more careful, please!
	lx[cpu_x-1]--; ly[cpu_y-1]--; lz[cpu_z-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,nel_x,nel_y,nel_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz, &user->FDSTAG.DA_CENTER); 	 CHKERRQ(ierr);
	lx[cpu_x-1]++; ly[cpu_y-1]++; lz[cpu_z-1]++;
//	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_STAR,nel_x,nel_y,nel_z,cpu_x,cpu_y,cpu_z,1,1,0,0,0, &user->FDSTAG.DA_CENTER); 	 CHKERRQ(ierr);

	// Create a vector for each of the history variables we need
	for (i=0;  i<user->num_phases; i++){
		// Phase proportions
		ierr = DMCreateGlobalVector(user->FDSTAG.DA_CENTER,&user->FDSTAG.Center_PhaseProportions[i]); 	CHKERRQ(ierr);
		ierr = VecSet(user->FDSTAG.Center_PhaseProportions[i], 	0.0);									CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(user->FDSTAG.DA_CENTER,&user->FDSTAG.Center_Temperature); 				CHKERRQ(ierr);	// Temperature
	ierr = VecSet(user->FDSTAG.Center_Temperature, 	0.0);												CHKERRQ(ierr);
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_Pressure); 				CHKERRQ(ierr);	// Pressure
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_NumParticles); 			CHKERRQ(ierr);	// # of particles
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_Strain); 					CHKERRQ(ierr);	// Strain
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_PlasticStrain); 			CHKERRQ(ierr);	// Plastic strain
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_T2nd); 					CHKERRQ(ierr);	// 2nd invariant of deviatoric stress tensor
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_E2nd); 					CHKERRQ(ierr);	// 2nd invariant of strain rate tensor
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_EffectiveViscosity); 		CHKERRQ(ierr);	// Effective viscosity
	ierr = VecDuplicate(user->FDSTAG.Center_Temperature,&user->FDSTAG.Center_Density); 					CHKERRQ(ierr);	// Density

	// Corner values DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,&user->FDSTAG.DA_CORNER); CHKERRQ(ierr);

	// Create a vector for each of the history variables we need
	for (i=0;  i<user->num_phases; i++)
	{	// Phase proportions
		ierr = DMCreateGlobalVector(user->FDSTAG.DA_CORNER,&user->FDSTAG.Corner_PhaseProportions[i]); 	CHKERRQ(ierr);
		ierr = VecSet(user->FDSTAG.Corner_PhaseProportions[i], 	0.0);									CHKERRQ(ierr);
		ierr = DMCreateLocalVector(user->FDSTAG.DA_CORNER,&user->FDSTAG.Corner_PhaseProportions_local[i]);	CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(user->FDSTAG.DA_CORNER,&user->FDSTAG.Corner_Temperature); 	    	CHKERRQ(ierr);	// Temperature
	ierr = VecSet(user->FDSTAG.Corner_Temperature, 0.0);											CHKERRQ(ierr);
	ierr = VecDuplicate(user->FDSTAG.Corner_Temperature,&user->FDSTAG.Corner_Pressure); 			CHKERRQ(ierr);	// Pressure
	ierr = VecDuplicate(user->FDSTAG.Corner_Temperature,&user->FDSTAG.Corner_Density); 				CHKERRQ(ierr);	// Density
	ierr = VecDuplicate(user->FDSTAG.Corner_Temperature,&user->FDSTAG.Corner_HeatCapacity); 		CHKERRQ(ierr);	// Heat capacity
	ierr = VecDuplicate(user->FDSTAG.Corner_Temperature,&user->FDSTAG.Corner_Conductivity); 		CHKERRQ(ierr);	// conductivity
	ierr = VecDuplicate(user->FDSTAG.Corner_Temperature,&user->FDSTAG.Corner_RadioactiveHeat); 		CHKERRQ(ierr);	// radioactive heat production

	// In addition to viscosity at center points we also need viscosity at the vertices (to compute Sxz, Sxy and Sxz).
	// Corner values are similar to the velocity DA; yet the velocity DMDA has 3 dof's which we don't want here.
	// For this reason, an additional DMDA is constructed with only 1 DOF.
	// For some reason (that is not yet well understood), we need to use a stencilwidth of 2 here

	// Viscosity at Sxy points // DMDA_STENCIL_STAR ??? stencil width = 2 ???
// a bit more careful, please!
	lz[cpu_z-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nnode_x,nnode_y,nel_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,&user->FDSTAG.DA_XY_POINTS);
	lz[cpu_z-1]++;
//	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,nnode_x,nnode_y,nel_z,cpu_x,cpu_y,cpu_z,1,2,0,0,0,&user->FDSTAG.DA_XY_POINTS);

	for (i=0;  i<user->num_phases; i++){
		// Phase proportions
		ierr = DMCreateGlobalVector(user->FDSTAG.DA_XY_POINTS,&user->FDSTAG.XYPoints_PhaseProportions[i]); 	CHKERRQ(ierr);
		ierr = VecSet(user->FDSTAG.XYPoints_PhaseProportions[i], 	0.0);									CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(user->FDSTAG.DA_XY_POINTS,	&user->FDSTAG.XYPoints_Temperature); 			CHKERRQ(ierr);	// Temperature
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_Pressure); 				CHKERRQ(ierr);	// Pressure
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_Strain); 				CHKERRQ(ierr);	// Strain
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_PlasticStrain); 			CHKERRQ(ierr);	// Plastic strain
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_T2nd); 					CHKERRQ(ierr);	// 2nd invariant of deviatoric stress tensor
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_E2nd); 					CHKERRQ(ierr);	// 2nd invariant of strain rate tensor
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_EffectiveViscosity); 	CHKERRQ(ierr);	// Effective viscosity
	ierr = VecDuplicate(user->FDSTAG.XYPoints_Temperature,	&user->FDSTAG.XYPoints_Density); 				CHKERRQ(ierr);	// Density

	// Viscosity at Sxz points // DMDA_STENCIL_STAR ??? stencil width = 2 ???
// a bit more careful, please!
	ly[cpu_y-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nnode_x,nel_y,nnode_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,&user->FDSTAG.DA_XZ_POINTS);
	ly[cpu_y-1]++;
//	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,nnode_x,nel_y,nnode_z,cpu_x,cpu_y,cpu_z,1,2,0,0,0,&user->FDSTAG.DA_XZ_POINTS);

	for (i=0;  i<user->num_phases; i++){
		// Phase proportions
		ierr = DMCreateGlobalVector(user->FDSTAG.DA_XZ_POINTS,&user->FDSTAG.XZPoints_PhaseProportions[i]); 	CHKERRQ(ierr);
		ierr = VecSet(user->FDSTAG.XZPoints_PhaseProportions[i], 0.0);									CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(user->FDSTAG.DA_XZ_POINTS,	&user->FDSTAG.XZPoints_Temperature); 			CHKERRQ(ierr);	// Temperature
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_Pressure); 				CHKERRQ(ierr);	// Pressure
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_Strain); 				CHKERRQ(ierr);	// Strain
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_PlasticStrain); 			CHKERRQ(ierr);	// Plastic strain
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_T2nd); 					CHKERRQ(ierr);	// 2nd invariant of deviatoric stress tensor
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_E2nd); 					CHKERRQ(ierr);	// 2nd invariant of strain rate tensor
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_EffectiveViscosity); 	CHKERRQ(ierr);	// Effective viscosity
	ierr = VecDuplicate(user->FDSTAG.XZPoints_Temperature,	&user->FDSTAG.XZPoints_Density); 				CHKERRQ(ierr);	// Density

	// Viscosity at Syz points // DMDA_STENCIL_STAR ??? stencil width = 2 ???
// a bit more careful, please!
	lx[cpu_x-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nel_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,1,1,lx,ly,lz,&user->FDSTAG.DA_YZ_POINTS);
	lx[cpu_x-1]++;
//	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,nel_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,1,2,0,0,0,&user->FDSTAG.DA_YZ_POINTS);

	for (i=0;  i<user->num_phases; i++){
		// Phase proportions
		ierr = DMCreateGlobalVector(user->FDSTAG.DA_YZ_POINTS,&user->FDSTAG.YZPoints_PhaseProportions[i]); 	CHKERRQ(ierr);
		ierr = VecSet(user->FDSTAG.YZPoints_PhaseProportions[i], 	0.0);									CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(user->FDSTAG.DA_YZ_POINTS,	&user->FDSTAG.YZPoints_Temperature); 			CHKERRQ(ierr);	// Temperature
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_Pressure); 				CHKERRQ(ierr);	// Pressure
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_Strain); 				CHKERRQ(ierr);	// Strain
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_PlasticStrain); 			CHKERRQ(ierr);	// Plastic strain
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_T2nd); 					CHKERRQ(ierr);	// 2nd invariant of deviatoric stress tensor
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_E2nd); 					CHKERRQ(ierr);	// 2nd invariant of strain rate tensor
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_EffectiveViscosity); 	CHKERRQ(ierr);	// Effective viscosity
	ierr = VecDuplicate(user->FDSTAG.YZPoints_Temperature,	&user->FDSTAG.YZPoints_Density); 				CHKERRQ(ierr);	// Density

	// Set coordinates of DMDA (for plotting later)
	ierr = SetUniformCoordinates_FDSTAG(user); CHKERRQ(ierr);

	// Create stiffness matrices
	// The nomenclature for the global stiffness matrices is:
	// | VV_MAT VP_MAT | | V |    |Rhs |
	// |               | |   |  = |    |
	// | PV_MAT PP_MAT | | P |    | 0  |

	ierr = DMGetGlobalVector(user->DA_Vel,&tmp_u);CHKERRQ(ierr);
	ierr = DMGetGlobalVector(user->DA_Pres,&tmp_p);CHKERRQ(ierr);
	ierr = VecGetSize(tmp_u, &vel_sdof); CHKERRQ(ierr);	 // global # of velocity unknowns
	ierr = VecGetSize(tmp_p, &pres_sdof); CHKERRQ(ierr); // global # of pressure unknowns
	ierr = VecGetOwnershipRange(tmp_u,&nvel_rowStart,&nvel_rowEnd); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(tmp_p,&npres_rowStart,&npres_rowEnd); CHKERRQ(ierr);
	nvel_rows_Local  = nvel_rowEnd - nvel_rowStart;
	npres_rows_Local = npres_rowEnd - npres_rowStart;
	ierr = DMRestoreGlobalVector(user->DA_Vel,&tmp_u);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(user->DA_Pres,&tmp_p);CHKERRQ(ierr);

	//=============================================================================================================
	// Create velocity stiffness matrix (A_11)
	ierr = MatCreate(PETSC_COMM_WORLD,&user->VV_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->VV_MAT, nvel_rows_Local, nvel_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);

	if(user->StokesSolver != 1)
	{	// Not good for Powel-Hesteness iterations
		// Matrix duplication complains about block sizes
		ierr = MatSetBlockSize( user->VV_MAT, 3 ); CHKERRQ(ierr);
	}

	ierr = MatSetType(user->VV_MAT,	MATAIJ); CHKERRQ(ierr);


	// constant preallocation for A_11
	ierr = MatSeqAIJSetPreallocation(user->VV_MAT, 15, PETSC_NULL); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(user->VV_MAT, 15, PETSC_NULL, 15, PETSC_NULL); CHKERRQ(ierr);
	ierr = MatSetOption(user->VV_MAT, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

	// variable preallocation for A_11
	//ierr = FDSTAG_Preallocation_A11FDSTAG( user->DA_Processors, user->DA_Vel,user->DA_Pres, user->VV_MAT,
	//		nvel_rowStart, nvel_rowEnd ); CHKERRQ(ierr);

	//=============================================================================================================
	// Create PP matrix (A_22), aka pressure mass matrix (this only needs to be a diagonal matrix)
	ierr = MatCreate(PETSC_COMM_WORLD,&user->PP_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->PP_MAT, npres_rows_Local, npres_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType(user->PP_MAT,	MATAIJ); CHKERRQ(ierr);

	// constant preallocation for A_22
	ierr = MatSeqAIJSetPreallocation(user->PP_MAT,1,PETSC_NULL); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(user->PP_MAT,1,PETSC_NULL,0,PETSC_NULL); CHKERRQ(ierr);
	ierr = MatSetOption(user->PP_MAT, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

	// assemble A_22 now as MatDuplicate only functions on assembled matrices
	for (i=npres_rowStart; i<npres_rowEnd; i++) {
		ierr = MatSetValue(user->PP_MAT,i,i,0.0,INSERT_VALUES);CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(user->PP_MAT,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(user->PP_MAT,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	//=============================================================================================================
	// Create VP matrix (A_12), with the correct parallel layout (dictated by the VV & PP matrices).
	ierr = MatCreate(PETSC_COMM_WORLD,&user->VP_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->VP_MAT, nvel_rows_Local, npres_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType(user->VP_MAT,	MATAIJ); CHKERRQ(ierr);

	// constant preallocation for A_12
	ierr = MatSeqAIJSetPreallocation(user->VP_MAT,2,PETSC_NULL); CHKERRQ(ierr); // dp/dx, dp/dy, dp/dz computed using 2 points
	ierr = MatMPIAIJSetPreallocation(user->VP_MAT,2,PETSC_NULL,2,PETSC_NULL); CHKERRQ(ierr); // coarse preallocation
	ierr = MatSetOption(user->VP_MAT, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

	// variable preallocation for A_12
	//ierr = FDSTAG_Preallocation_A12FDSTAG( user->DA_Processors, user->DA_Vel,user->DA_Pres, user->VP_MAT,
	//		nvel_rowStart,  nvel_rowEnd,
	//		npres_rowStart, npres_rowEnd); CHKERRQ(ierr);
	//=============================================================================================================
	// Create PV matrix (A_21), with the correct parallel layout (dictated by the VV & PP matrixes).
	ierr = MatCreate(PETSC_COMM_WORLD, &user->PV_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(user->PV_MAT, npres_rows_Local, nvel_rows_Local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType(user->PV_MAT,	MATAIJ); CHKERRQ(ierr);

	// constant preallocation for A_21
	ierr = MatSeqAIJSetPreallocation(user->PV_MAT,6,PETSC_NULL); CHKERRQ(ierr);  // div(u) computed using 6 points
	ierr = MatMPIAIJSetPreallocation(user->PV_MAT,6,PETSC_NULL,6,PETSC_NULL); CHKERRQ(ierr); // coarse preallocation
	ierr = MatSetOption(user->PV_MAT, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

	// variable preallocation for A_21
//	ierr = FDSTAG_Preallocation_A21FDSTAG( user->DA_Processors, user->DA_Vel, user->DA_Pres, user->PV_MAT,
//			nvel_rowStart,  nvel_rowEnd,
//			npres_rowStart, npres_rowEnd); CHKERRQ(ierr);
	//=============================================================================================================

	// Store cpu-distribution to create TEMP DA
	user->cpu_x = cpu_x;
	user->cpu_y = cpu_y;
	user->cpu_z = cpu_z;

	ierr = MatDuplicate( user->PP_MAT, MAT_COPY_VALUES, &user->approx_S ); CHKERRQ(ierr);

	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==============================================================================================================================
// Initialize matrixes required for the Temperature solver - Driver routines that calls the appropriate routine depending on
// whether FDSTAG/FEM is being used.
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_TemperatureSolver"
PetscErrorCode LaMEM_Initialize_TemperatureSolver( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C)
{
	PetscErrorCode 	ierr;
	if (vpt_element_type == DAVP_FDSTAG ){
		// Finite Difference
		ierr = LaMEM_Initialize_TemperatureSolver_FDSTAG( user, vpt_element_type); CHKERRQ(ierr);
	}
	else{
		// FEM
		ierr = LaMEM_Initialize_TemperatureSolver_FEM( user, vpt_element_type, C); CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}
//==============================================================================================================================
// Initialize matrixes required for the Temperature solver in case of FEM
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_TemperatureSolver_FEM"
PetscErrorCode LaMEM_Initialize_TemperatureSolver_FEM( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C)
{
	PetscMPIInt     size;
	PetscErrorCode 	ierr;
	PetscInt 		dof_temp, nnode_x, nnode_y, nnode_z, cpu_x, cpu_y, cpu_z;
	PetscInt		*lx, *ly, *lz;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	ierr = PetscMalloc((size_t)user->cpu_x*sizeof(PetscInt), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)user->cpu_y*sizeof(PetscInt), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)user->cpu_z*sizeof(PetscInt), &lz); CHKERRQ(ierr);
	ierr = ComputeNodesPerProcessor(user->DA_Processors,lx,ly,lz,vpt_element_type); CHKERRQ(ierr);

	nnode_x = 	user->finest_nnode_x;
	nnode_y = 	user->finest_nnode_y;
	nnode_z = 	user->finest_nnode_z;
	cpu_x 	=	user->cpu_x;
	cpu_y 	=	user->cpu_y;
	cpu_z 	=	user->cpu_z;

	// correct for periodicity
	if (user->BC.LeftBound==3){
		nnode_x	 		= nnode_x-1;
		lx[cpu_x-1] 	= lx[cpu_x-1]-1;
	}
	if (user->BC.FrontBound==3){
		nnode_y	 	= nnode_y-1;
		ly[cpu_y-1]	= ly[cpu_y-1]-1;
	}

	dof_temp = 1;
	if ((C->ElementType==2) ){
		// Q1P0 or Q1Q1
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof_temp,1,lx,ly,lz,&user->DA_Temp); CHKERRQ(ierr);
	}
	else {
		// Q2Pm1
		ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof_temp,2,lx,ly,lz,&user->DA_Temp); CHKERRQ(ierr);
	}
	ierr = DMDASetFieldName(	user->DA_Temp,0,"Temp"); CHKERRQ(ierr);

	ierr = DMCreateMatrix(user->DA_Temp,&user->TEMP_MAT); CHKERRQ(ierr);

	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==============================================================================================================================*/

/*==============================================================================================================================
Initialize matrixes required for the Temperature solver in case of using FDSTAG
*/
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_TemperatureSolver_FDSTAG"
PetscErrorCode LaMEM_Initialize_TemperatureSolver_FDSTAG( UserContext *user, DAVPElementType vpt_element_type)
{
	PetscMPIInt     size;
	PetscErrorCode 	ierr;
	PetscInt 		dof_temp, nnode_x, nnode_y, nnode_z, cpu_x, cpu_y, cpu_z;
	PetscInt		*lx, *ly, *lz;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	ierr = PetscMalloc((size_t)user->cpu_x*sizeof(PetscInt), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)user->cpu_y*sizeof(PetscInt), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)user->cpu_z*sizeof(PetscInt), &lz); CHKERRQ(ierr);
	ierr = ComputeNodesPerProcessor(user->DA_Processors,lx,ly,lz,vpt_element_type); CHKERRQ(ierr);

	nnode_x = 	user->finest_nnode_x;
	nnode_y = 	user->finest_nnode_y;
	nnode_z = 	user->finest_nnode_z;
	cpu_x 	=	user->cpu_x;
	cpu_y 	=	user->cpu_y;
	cpu_z 	=	user->cpu_z;

	// correct for periodicity
	if (user->BC.LeftBound==3){
		nnode_x	 		= nnode_x-1;
		lx[cpu_x-1] 	= lx[cpu_x-1]-1;
	}
	if (user->BC.FrontBound==3){
		nnode_y	 	= nnode_y-1;
		ly[cpu_y-1]	= ly[cpu_y-1]-1;
	}

	dof_temp = 1;
	// Q1P0, Q1Q1 or FDSTAG
	ierr = DMDACreate3d(PETSC_COMM_WORLD,user->BC.BCType_x,user->BC.BCType_y,user->BC.BCType_z,DMDA_STENCIL_BOX,nnode_x,nnode_y,nnode_z,cpu_x,cpu_y,cpu_z,dof_temp,1,lx,ly,lz,&user->DA_Temp); CHKERRQ(ierr);
	ierr = DMDASetFieldName(	user->DA_Temp,0,"Temp"); CHKERRQ(ierr);


	/* if using FDSTAG, create global vectors that hold info about Cp, k and Q @ nodes */
	ierr = VecSet(user->FDSTAG.Corner_HeatCapacity, 	0.0); 				CHKERRQ(ierr);
	ierr = VecSet(user->FDSTAG.Corner_Conductivity, 	0.0); 				CHKERRQ(ierr);
	ierr = VecSet(user->FDSTAG.Corner_RadioactiveHeat, 	0.0); 				CHKERRQ(ierr);

	/* Get matrix */
	ierr = DMCreateMatrix(user->DA_Temp,&user->TEMP_MAT); CHKERRQ(ierr);

	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==============================================================================================================================*/

/*==============================================================================================================================*/
/*
Check some of the input variable to ensure that they are sane
*/
#undef __FUNCT__
#define __FUNCT__ "LaMEM_Initialize_CheckInput"
PetscErrorCode LaMEM_Initialize_CheckInput(void)
{
	PetscMPIInt size;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	if(size>MaxNumCPU) {
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP,"The maximum allowable cpus LaMEM can be run on is %lld. Please update MaxNumCPU in LaMEM.h", (LLD)MaxNumCPU );
	}


	PetscFunctionReturn(0);
}
/*==============================================================================================================================*/




/*==========================================================================================================*/
/* Initialize an internal erosion surface on rank zero, which will contain a high-resolution DMDA
 * that typically has a higher resolution than the numerical grid for the 3D deformation.
 */
#undef __FUNCT__
#define __FUNCT__ "InitializeInternalErosionSurfaceOnRankZero"
PetscErrorCode InitializeInternalErosionSurfaceOnRankZero(UserContext *user )
{
	PetscErrorCode	 	ierr;
	PetscMPIInt			rank;


	PetscInt			Nx_erosion, Ny_erosion;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	if ((rank==0) && (user->ErosionParameters.ErosionModel==2) ){
		/* Only perform this routine if we selected the FD erosion model, we are on rank 0 AND we have an internal free surface
		 * as the code currently does not work with the FEM part of LaMEM
		 */
		PetscScalar 	z_surface;
		PetscRandom		rctx;
		PetscScalar		**Topography;
		PetscInt 		xs,ys,xm,ym,ix,iy, Nx,Ny;

		/* Catch potential errors */
		if (user->ErosionParameters.UseInternalFreeSurface==0){
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "The finite difference-based erosion code can only be used in combination with an internal free surface!");
		}

		/* Resolution of erosion DA */
		Nx_erosion = user->ErosionParameters.FE_ErosionCode.ResolutionFactorX*(user->finest_nnode_x-1) + 1;
		Ny_erosion = user->ErosionParameters.FE_ErosionCode.ResolutionFactorY*(user->finest_nnode_y-1) + 1;

		/* Initialize the DA */
		ierr =  DMDACreate2d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR, Nx_erosion,Ny_erosion,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,&user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode);

		/* Set uniform coordinates */
		ierr = DMDASetUniformCoordinates(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,0,0); CHKERRQ(ierr);

		/* Get vector related to this DA */
		ierr = DMCreateGlobalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&user->ErosionParameters.FE_ErosionCode.ErosionSurface);

		/* Set the eroded surface initially to the same level as the internal free surface */
		z_surface = user->z_bot	+ user->H;
		if (user->ErosionParameters.UseInternalFreeSurface == 1){
			// We use an internal free surface
			z_surface = user->ErosionParameters.InitialFreeSurfaceHeight;
		}

		// initialize (internal) surface topography, which will be in NONDIMENSIONAL units to be consistent with the rest of the code.
		// The erosion code can potentially be in dimensional units.
		ierr = VecSet(user->ErosionParameters.FE_ErosionCode.ErosionSurface,z_surface); CHKERRQ(ierr);

		/* Add random noise on the surface if requested */
		ierr = 	PetscRandomCreate(PETSC_COMM_SELF, &rctx); 	CHKERRQ(ierr);
		ierr = 	PetscRandomSetFromOptions(rctx); 			CHKERRQ(ierr);

		ierr = 	DMDAGetInfo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,PETSC_NULL,&Nx,&Ny,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
		ierr = 	DMDAGetCorners(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&xs,&ys, PETSC_NULL, &xm,&ym, PETSC_NULL); CHKERRQ(ierr); // required for adding noise
		ierr = 	DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,user->ErosionParameters.FE_ErosionCode.ErosionSurface,	&Topography); CHKERRQ(ierr);
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){
				PetscScalar	cf_rand, fac;

				/* Add random noise */
				ierr 				= 	PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				cf_rand 			= 	cf_rand-0.5;
				Topography[iy][ix] 	=	Topography[iy][ix] + cf_rand*user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m/user->Characteristic.Length;	// add random noise in ND units

				/* Add inclination */
				//fac					=	((PetscScalar) ix)/((PetscScalar) Nx);
				fac					=	((PetscScalar) iy)/((PetscScalar) Ny);

				Topography[iy][ix] 	=	Topography[iy][ix] - (1-fac)*user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m/user->Characteristic.Length;	// Uplift in x-direction

			}
		}
		ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,user->ErosionParameters.FE_ErosionCode.ErosionSurface,	&Topography); CHKERRQ(ierr);
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);         // Destroy random context

		PetscPrintf(PETSC_COMM_SELF,"# Finished initializing FD erosion code\n");

		// check that the Topography was set correctly
		//VecView(user->ErosionParameters.FE_ErosionCode.ErosionSurface, PETSC_VIEWER_STDOUT_SELF);

	}




PetscFunctionReturn(0);
}
/*==========================================================================================================*/



