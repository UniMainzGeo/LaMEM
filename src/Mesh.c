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

Mesh.c, contains the following subroutines:

InterpolateMeshes				-	Interpolate meshes from fine to coarse grid
AdvectGrid						-	Advect grid in lagrangian fashion, were required
ComputeNodesPerProcessor		-	Compute number of nodes for each processor, given number of elements
SaveInitialMesh					-
ReadMeshFromFile				-
GenerateMeshFine				-	Generate the FEM mesh
SetBoundaryConditionsRHS		-	Set BC values to RHS vector
ComputeGlobalProperties			-	Compute global properties such as Vrms at every timestep
ApplyErosion					-	Apply one of several erosion models to the surface of the mesh

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Mesh.h"
#include "Utils.h"
#include "LaMEM_FE_ErosionCode.h"


#undef __FUNCT__
#define __FUNCT__ "DAUpdatedGhostedCoordinates"
PetscErrorCode DAUpdatedGhostedCoordinates( DM da )
{
	PetscErrorCode ierr;
	Vec da_coordinates, gcoords;
	DM _dac;

	ierr = DMGetCoordinateDM(da,&_dac);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
	ierr = DMGetCoordinates( da, &da_coordinates );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(_dac,da_coordinates,INSERT_VALUES,gcoords);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(_dac,da_coordinates,INSERT_VALUES,gcoords);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DASetCoordinatesFromLocalVector"
/*@
 DASetCoordinatesFromLocalVector - Sets into the DA's global coordinate vector
 the coordinate value stored in a local vector (ghost nodes are ignored).
 The local vector has the same parallel layout as the one obtained from
 calling DMGetCoordinatesLocal().

 Not Collective

 Input Parameter:
 +  da - the distributed array
 -  local_coords - local coordinate vector

 Note:
 The coordinates may include those for all ghost points, however ghost values
 will be ignored when the local coords are inserted into the global coordinate vector.

 The user is responsible for destroying the vector local_coords.

 Level: intermediate

 .keywords: distributed array, get, corners, nodes, local indices, coordinates

 .seealso: DMDAGetGhostCorners(), DMGetCoordinates(), DAGetGhostCoordinates(), DMDASetUniformCoordinates(), DMGetCoordinateDM()
 @*/
PetscErrorCode DASetCoordinatesFromLocalVector( DM da, Vec local_coords )
{
	PetscErrorCode ierr;
	Vec da_coordinates;
	DM dac;

	ierr = DMGetCoordinateDM(da,&dac);CHKERRQ(ierr);

	/* create new vector to store new coords */
	//ierr = DMCreateGlobalVector( dac, &global_coords );CHKERRQ(ierr);
	
	/* scatter new local coords into global_coords */
	//ierr = DMLocalToGlobal( dac, local_coords, INSERT_VALUES, global_coords );CHKERRQ(ierr);
	
	/* scatter new existing coords into global_coords */
	ierr = DMGetCoordinates( da, &da_coordinates );
	ierr = DMLocalToGlobalBegin( dac, local_coords, INSERT_VALUES, da_coordinates );CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  ( dac, local_coords, INSERT_VALUES, da_coordinates );CHKERRQ(ierr);
	
	//	ierr = DMSetCoordinates( da, da_coordinates );CHKERRQ(ierr); /* do not free global coords */
	
	ierr = DAUpdatedGhostedCoordinates(da);CHKERRQ(ierr);
	
#if(  (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==1) )
	//	ierr = VecDestroy( &da_coordinates );
#endif
	//	ierr = DMDestroy(&dac);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


double round_double( double d, PetscInt places ) {
	return round(d * pow(10, (double) places)) / pow(10,(double) places);
}

PetscErrorCode DARoundCoordinates( DM da, const PetscInt places )
{
	PetscErrorCode ierr;
	PetscInt xs,ys,zs,xm,ym,zm, i,j,k;
	DM cda;
	Vec global_coords;
	DMDACoor3d ***coords;

	ierr = DMGetCoordinates( da, &global_coords ); CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);

	ierr = DMDAVecGetArray( cda, global_coords, &coords ); CHKERRQ(ierr);
	ierr = DMDAGetCorners( da,&xs,&ys,&zs,&xm,&ym,&zm ); CHKERRQ(ierr);
	for( k=zs; k<zs+zm; k++ ) {
		for( j=ys; j<ys+ym; j++ ) {
			for( i=xs; i<xs+xm; i++ ) {
				double txi,tyi,tzi;

				txi = round_double( coords[k][j][i].x, places );
				tyi = round_double( coords[k][j][i].y, places );
				tzi = round_double( coords[k][j][i].z, places );

				coords[k][j][i].x = (double)txi;
				coords[k][j][i].y = (double)tyi;
				coords[k][j][i].z = (double)tzi;

			}}}


	ierr = DMDAVecRestoreArray( cda, global_coords, &coords ); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy( global_coords ); CHKERRQ(ierr);

	/* propagate rounding into ghost arrays */
	ierr = DAUpdatedGhostedCoordinates( da ); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Advect the grid in a Lagrangian manner with given velocity and timestep */
#undef __FUNCT__
#define __FUNCT__ "AdvectGrid"
PetscErrorCode AdvectGrid( DM da, Vec Velocity, DM da_coord, UserContext *user )
{
	PetscMPIInt rank;
	PetscErrorCode ierr;
	DM				cda;
	DMDACoor3d		***coors;
	Vec				local_Vel, gc;
	PetscInt		xs,ys,zs,xm,ym, zm, ix,iy,iz, nnode_x, nnode_y, nnode_z, dummy;
	PetscScalar		Vx,Vy,Vz;
	Field 			***velocity;
	PetscBool		FastErosion;
	//PetscScalar 	maxZ, minZ;
	PetscBool flg,ApplyKinematicCondition_FromElement_Z;


	ierr = DACompareStructures( da, da_coord, &flg ); CHKERRQ(ierr);
	if( flg == PETSC_FALSE ) {
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "DA'S MUST MATCH");
	}

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	ierr = DMGetCoordinateDM(da_coord,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da_coord,    &gc); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da_coord,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// make sure to use Ghost corners!!!
	ierr = DMDAGetInfo(da_coord, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);

	/* Copy global velocity solution to local processor, including ghostpoints */
	ierr = DMGetLocalVector(da,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_Vel,	&velocity); CHKERRQ(ierr);


	ierr = PetscOptionsGetBool( PETSC_NULL, "-FastErosion", 							&flg  	, &FastErosion); 																								CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL , "-ApplyKinematicCondition_FromElement_Z",	&dummy	, &ApplyKinematicCondition_FromElement_Z);	CHKERRQ(ierr);

	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				Vx = 	velocity[iz][iy][ix].Vx;
				Vy = 	velocity[iz][iy][ix].Vy;
				Vz = 	velocity[iz][iy][ix].Vz;

				/*Special case if we only push part of the domain.
				 * Note: this worked, but was found to be not very useful for detachment folding which is why I commented it again
				 * BK, June 2011
				 * */
				//if ( ((ApplyKinematicCondition_FromElement_Z) && (ix==0) ) ||  ((ApplyKinematicCondition_FromElement_Z) && (ix==nnode_x-1) ) ){
				//	Vx =	-user->BC.Exx*coors[iz][iy][ix].x;
				//}



				/* Advection of the grid is done if
				 * 	(1) we have a fully lagrangian grid
				 *  (2) we have an ALE grid, and we're close to the surface
				 */
				if ( ((iz> (nnode_z-user->NumSurfaceNodes)) && (user->GridAdvectionMethod==2)) ||
						(user->GridAdvectionMethod==1) ){

					if ( (FastErosion)   && (iz==nnode_z-1) ){
						coors[iz][iy][ix].x = coors[iz][iy][ix].x + user->dt*Vx;
						coors[iz][iy][ix].y = coors[iz][iy][ix].y + user->dt*Vx;
						coors[iz][iy][ix].z = coors[iz][iy][ix].z;

					}
					else{
						coors[iz][iy][ix].x = coors[iz][iy][ix].x + user->dt*Vx;
						coors[iz][iy][ix].y = coors[iz][iy][ix].y + user->dt*Vy;
						coors[iz][iy][ix].z = coors[iz][iy][ix].z + user->dt*Vz;
					}

				}



			}
		}
	}

	///* For debugging: get max & min height of free surface*/
	//maxZ = 0;
	//for (iy=ys; iy<ys+ym; iy++){
	//	for(ix=xs; ix<xs+xm; ix++){
	//		maxZ	= PetscMax(coors[zs+zm-1][iy][ix].z,maxZ);
	//	}
	//}

	ierr = DMDAVecRestoreArray(da,local_Vel,&velocity); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&local_Vel); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);
	//ierr = DMDestroy( cda ); CHKERRQ(ierr);

	ierr = DASetCoordinatesFromLocalVector( da_coord, gc ); CHKERRQ(ierr);
	//ierr = VecDestroy( gc ); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the number of nodes on each processor, given a DM that defines the number of
 * elements in each direction on each processor. */
#undef __FUNCT__
#define __FUNCT__ "ComputeNodesPerProcessor"
PetscErrorCode ComputeNodesPerProcessor( DM da, PetscInt *lx, PetscInt *ly,PetscInt *lz, const DAVPElementType element_type )
{
	PetscMPIInt rank, size;
	PetscErrorCode ierr;
	PetscInt 	i;
	PetscInt	cpu_x, cpu_y, cpu_z;
	PetscInt	mx, my, mz;
	PetscInt 	sumx,sumy,sumz;
	const PetscInt *_lx,*_ly,*_lz;

	ierr = DMDAGetInfo(da, 0, &mx, &my, &mz, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	ierr = DMDAGetOwnershipRanges( da, &_lx,&_ly,&_lz); CHKERRQ(ierr);

	if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) || (element_type==DAVP_FDSTAG)) {
		//
		for( i=0; i<cpu_x; i++ ) {
			lx[i] = _lx[i];
		}
		lx[cpu_x-1] = _lx[cpu_x-1] + 1;

		//
		for( i=0; i<cpu_y; i++ ) {
			ly[i] = _ly[i];
		}
		ly[cpu_y-1] = _ly[cpu_y-1] + 1;

		//
		for( i=0; i<cpu_z; i++ ) {
			lz[i] = _lz[i];
		}
		lz[cpu_z-1] = _lz[cpu_z-1] + 1;

	}
	else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ) {
		//
		for( i=0; i<cpu_x; i++ ) {
			lx[i] = 2*_lx[i];
		}
		lx[cpu_x-1] = 2*_lx[cpu_x-1] + 1;

		//
		for( i=0; i<cpu_y; i++ ) {
			ly[i] = 2*_ly[i];
		}
		ly[cpu_y-1] = 2*_ly[cpu_y-1] + 1;

		//
		for( i=0; i<cpu_z; i++ ) {
			lz[i] = 2*_lz[i];
		}
		lz[cpu_z-1] = 2*_lz[cpu_z-1] + 1;

	}
	else{
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown element type!");
	}


	sumx = sumy = sumz = 0;
	for (i=0; i<cpu_x; i++){ sumx += lx[i]; }
	for (i=0; i<cpu_y; i++){ sumy += ly[i]; }
	for (i=0; i<cpu_z; i++){ sumz += lz[i]; }


	if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ) {
		if( sumx != 2*mx+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in X: %lld", (LLD)sumx );
		}
		if( sumy != 2*my+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in Y: %lld", (LLD)sumy );
		}
		if( sumz != 2*mz+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in Z: %lld", (LLD)sumz );
		}
	}

	if( (element_type==DAVP_Q1P0) ||  (element_type==DAVP_Q1Q1) || (element_type==DAVP_FDSTAG)) {
		if( sumx != mx+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in X: %lld", (LLD)sumx );
		}
		if( sumy != my+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in Y: %lld", (LLD)sumy );
		}
		if( sumz != mz+1 ) {
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "*** Velocity node count exceeds global dimension in Z: %lld", (LLD)sumz );
		}
	}
	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Save initial mesh - mainly used to adapt the mesh for crustal thickness */
#undef __FUNCT__
#define __FUNCT__ "SaveInitialMesh"
PetscErrorCode SaveInitialMesh( UserContext *user, DM da_nodes, const char *FileName)
{
	PetscMPIInt rank, size;
	PetscErrorCode 		ierr;
	PetscInt			mx,my,mz;
	PetscInt			xs,ys,zs,xm,ym,zm,ix,iy,iz,num;
	Vec					information, coord_x, coord_y, coord_z, coord, global;
	DM					cda;
	DMDACoor3d			***coords;
	PetscScalar			*coords_x, *coords_y, *coords_z;
	char				SaveFileName[PETSC_MAX_PATH_LEN];
	PetscViewer			view_out;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	/* Create output directory */
	if (rank==0){
		num = mkdir("InitialMesh", S_IRWXU);
	}

	// All other procs should wait
	MPI_Barrier(PETSC_COMM_WORLD);



	sprintf(SaveFileName,"./InitialMesh/%s.%lld.out",FileName,(LLD)rank);  // construct the filename


	ierr = DMDAGetInfo(da_nodes, 0, &mx,    &my, &mz, 0,0,0,0,0,0,0,0,0);CHKERRQ(ierr); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_nodes,      &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Mesh size: %lld %lld %lld \n",(LLD)mx,(LLD)my,(LLD)mz);

	/* Write info about the grid to a vector */
	ierr = VecCreate(PETSC_COMM_SELF,&information); CHKERRQ(ierr);
	ierr = VecSetSizes(information,PETSC_DECIDE,40); CHKERRQ(ierr);
	ierr = VecSetFromOptions(information); CHKERRQ(ierr);
	ierr = VecSetValue(information,0 ,	(PetscScalar)(mx), 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,1 ,	(PetscScalar)(my), 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,2 ,	(PetscScalar)(mz), 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,3 ,	(PetscScalar)(size), 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,4,	(PetscScalar)(xs)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,5,	(PetscScalar)(ys)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,6,	(PetscScalar)(zs)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,7,	(PetscScalar)(xm)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,8,	(PetscScalar)(ym)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,9,	(PetscScalar)(zm)	,				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,10,	 user->Characteristic.Length,	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,11,	 user->Characteristic.km	,	INSERT_VALUES); CHKERRQ(ierr);

	ierr = VecAssemblyBegin(information); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(information); CHKERRQ(ierr);

	/* Create local vectors that contains coordinates and velocity @ nodal points */
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_x); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_y); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_z); CHKERRQ(ierr);

	/* Get the coordinates from the DM and put them into vectors */
	ierr = DMGetCoordinateDM(da_nodes,&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da_nodes,&coord); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&coords); CHKERRQ(ierr);
	num = 0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){
				coords_x[num]   = coords[iz][iy][ix].x*user->Characteristic.Length;
				coords_y[num]   = coords[iz][iy][ix].y*user->Characteristic.Length;
				coords_z[num]   = coords[iz][iy][ix].z*user->Characteristic.Length;
				num			  	= num+1;
			}
		}
	}

	ierr = DMDAVecRestoreArray(cda,coord,&coords); CHKERRQ(ierr);
	//ierr = VecDestroy(coord); CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_nodes,&global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(cda,coord,INSERT_VALUES,global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (cda,coord,INSERT_VALUES,global); CHKERRQ(ierr);
	//	ierr = VecDestroy(&global); CHKERRQ(ierr);
	//	ierr = DMDestroy(&cda); CHKERRQ(ierr);
	
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm,coords_x,&coord_x); CHKERRQ(ierr);		VecAssemblyBegin(coord_x); VecAssemblyEnd(coord_x);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm,coords_y,&coord_y); CHKERRQ(ierr);		VecAssemblyBegin(coord_y); VecAssemblyEnd(coord_y);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm,coords_z,&coord_z); CHKERRQ(ierr);		VecAssemblyBegin(coord_z); VecAssemblyEnd(coord_z);

	// Write the actual file
//	ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_SELF,SaveFileName,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileName,FILE_MODE_WRITE,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(view_out, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

	ierr = VecView(information, 		view_out); CHKERRQ(ierr);
	ierr = VecView(coord_x, 			view_out); CHKERRQ(ierr); 		// save coordinate array
	ierr = VecView(coord_y, 			view_out); CHKERRQ(ierr); 		// save coordinate array
	ierr = VecView(coord_z, 			view_out); CHKERRQ(ierr); 		// save coordinate array

	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Saved Initial Mesh to file %s \n",SaveFileName);

	// Destroy vectors
	ierr = VecDestroy(&information); CHKERRQ(ierr);
	ierr = VecDestroy(&coord_x); CHKERRQ(ierr);
	ierr = VecDestroy(&coord_y); CHKERRQ(ierr);
	ierr = VecDestroy(&coord_z); CHKERRQ(ierr);
	ierr = PetscFree(coords_x); CHKERRQ(ierr);
	ierr = PetscFree(coords_y); CHKERRQ(ierr);
	ierr = PetscFree(coords_z); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Read a mesh from file */
#undef __FUNCT__
#define __FUNCT__ "ReadMeshFromFile"
PetscErrorCode ReadMeshFromFile( DM da, UserContext *user )
{
	PetscErrorCode 	ierr;
	DMDACoor3d 		***coords;
	DM				cda;
	Vec				gc, global, coord_x, coord_y, coord_z, information;
	PetscInt		xs,ys,zs,xm,ym,zm, nx,ny,nz,ix,iy,iz, ind[1], num;
	PetscScalar		n_read[1],  *coord_x_in, *coord_y_in, *coord_z_in;
	PetscViewer		view_in;

	ierr = DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions

	/* Read the input coordinates on cpu 0 */


	/* Load the file */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, user->InitialMeshFileName,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&information);CHKERRQ(ierr);
	ierr = VecLoad(information,view_in); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&coord_x);CHKERRQ(ierr);
	ierr = VecLoad(coord_x,view_in); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&coord_y);CHKERRQ(ierr);
	ierr = VecLoad(coord_y,view_in); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&coord_z);CHKERRQ(ierr);
	ierr = VecLoad(coord_z,view_in); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);	//close binary
	MPI_Barrier(PETSC_COMM_WORLD);
	

	/* Check that input grid size is conform with current grid size */
	ind[0] = 0; VecGetValues(information,1, ind, n_read);
	if (  ((PetscInt) n_read[0]) != nx  ){
		PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in x-direction !! %lld %lld",(LLD)n_read[0], (LLD)nx);
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	ind[0] = 1; VecGetValues(information,1, ind, n_read);
	if (  ((PetscInt) n_read[0]) != ny  ){
		PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in y-direction !! %lld %lld",(LLD)n_read[0], (LLD)ny);
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	ind[0] = 2; VecGetValues(information,1, ind, n_read);
	if (  ((PetscInt) n_read[0]) != nz  ){
		PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in z-direction !! %lld %lld",(LLD)n_read[0], (LLD)nz);
		MPI_Abort(PETSC_COMM_WORLD,1);
	}

	/* Create 3D arrays with file coordinates */
	ierr = VecGetArray(coord_x, &coord_x_in); CHKERRQ(ierr);
	ierr = VecGetArray(coord_y, &coord_y_in); CHKERRQ(ierr);
	ierr = VecGetArray(coord_z, &coord_z_in); CHKERRQ(ierr);


	/* Set the local portion of the coordinate vector */
	ierr = DMGetCoordinateDM(da,			&cda); 			CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,		&gc); 			CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,			&coords); 		CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				num = ix*(ny*nz) + iy*nz + iz;		// number, conform with the way the vector was created in MATLAB

				coords[iz][iy][ix].x = coord_x_in[num]/user->Characteristic.Length;
				coords[iz][iy][ix].y = coord_y_in[num]/user->Characteristic.Length;
				coords[iz][iy][ix].z = coord_z_in[num]/user->Characteristic.Length;

			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
	
	//	ierr = VecDestroy(&global); 								CHKERRQ(ierr);
	//	ierr = VecDestroy(&gc); 									CHKERRQ(ierr);
	//	ierr = DMDestroy(&cda); 									CHKERRQ(ierr);
	
	ierr = VecRestoreArray(coord_x,&coord_x_in); CHKERRQ(ierr);
	ierr = VecRestoreArray(coord_y,&coord_y_in); CHKERRQ(ierr);
	ierr = VecRestoreArray(coord_z,&coord_z_in); CHKERRQ(ierr);


	PetscPrintf(PETSC_COMM_WORLD,"Reading initial mesh from file: %s \n", user->InitialMeshFileName);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Generate the numerical FEM grid */
#undef __FUNCT__
#define __FUNCT__ "GenerateMeshFine"
PetscErrorCode GenerateMeshFine( DM da, UserContext *user )
{
	/* This routine generates a mesh, and adapts in the vertical direction.
	 * The horizontal coordinates, however, should be computed with DMDASetUniformCoordinates
	 * as it is otherwise quite tricky to take periodic BCs into account (DMDASetUniformCoordinates does it for us)
	 *
	 */
	PetscMPIInt rank, size;
	PetscErrorCode 	ierr;
	DMDACoor3d 		***coords;
	DM				cda;
	Vec				gc, global;
	PetscInt		xs,ys,zs,xm,ym,zm, nx,ny,nz,ix,iy,iz, i, indbot, indlay;
	PetscInt		*NodesDistributionGlobal;
	PetscScalar		dz, dzbot, dzlay;
	PetscRandom		rctx;

	//	ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);
	//	PetscFunctionReturn(0);

	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); 	CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx); 			CHKERRQ(ierr);

	ierr = DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);

	/* Send info about the distribution of the fine grid on the various CPUs to all processors */
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
	if( user->NodesDistributionCPUsFineGrid == PETSC_NULL ) {
		ierr = PetscMalloc( (size_t)size*6*sizeof(PetscInt),  	&user->NodesDistributionCPUsFineGrid ); CHKERRQ(ierr);
	}

	ierr = PetscMalloc( (size_t)size*6*sizeof(PetscInt),  	&NodesDistributionGlobal ); CHKERRQ(ierr);

	for (i=0; i<size*6; i++){ user->NodesDistributionCPUsFineGrid[i]=0; };
	user->NodesDistributionCPUsFineGrid[rank*6+0] = xs;		user->NodesDistributionCPUsFineGrid[rank*6+1] = ys;
	user->NodesDistributionCPUsFineGrid[rank*6+2] = zs;  	user->NodesDistributionCPUsFineGrid[rank*6+3] = xm;
	user->NodesDistributionCPUsFineGrid[rank*6+4] = ym;  	user->NodesDistributionCPUsFineGrid[rank*6+5] = zm;
	MPI_Reduce(user->NodesDistributionCPUsFineGrid,NodesDistributionGlobal, size*6, MPIU_INT ,MPI_SUM, 0	,PETSC_COMM_WORLD);
	MPI_Bcast(NodesDistributionGlobal, size*6, MPIU_INT, 0, PETSC_COMM_WORLD);
	for (i=0; i<size*6; i++){
		user->NodesDistributionCPUsFineGrid[i] =  NodesDistributionGlobal[i];
	}

	if (user->Setup.Model==0){
		// diapir setup
		PetscScalar 	fac, cf_rand;
		PetscInt		AdjustMeshToInterface;

		AdjustMeshToInterface = 1;
		if ((user->amplNoise==0) && (user->ampl2D==0)){
			AdjustMeshToInterface = 0;			// do not adjust the mesh in this case; important for an FDSTAG setup!
		}


		dz = user->H/((double) (nz-1));
		ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);


		ierr = DMGetCoordinateDM(da,&cda); 						CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gc); 				CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,&coords); 					CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);
		for(iz=zs; iz<zs+zm; iz++){
			for(iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){
					fac 	=	PetscAbs((coords[iz][iy][ix].z- (user->Hinterface*user->H ) ) /dz);

					if ((fac<0.5) && (iz>0) && (iz<nz-1) && (AdjustMeshToInterface==1)) {
						// do not deform the upper and lower boundary of the computational domain

						ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
						cf_rand = (cf_rand-0.5);

						coords[iz][iy][ix].z = user->Hinterface*user->H + cf_rand*user->amplNoise + user->ampl2D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)) );

					}

				}
			}
		}
		ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
		ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		//		ierr = VecDestroy(&gc); 									CHKERRQ(ierr);
		//		ierr = VecDestroy(&global); 								CHKERRQ(ierr);
		//		ierr = DMDestroy(&cda); 									CHKERRQ(ierr
	}

	else if (user->Setup.Model==5){		// single layer folding setup
		PetscInt nzbot, nzlay, nztop;

		dz = user->H/((double) (nz-1));
		ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);
		//PetscPrintf(PETSC_COMM_WORLD,"dx: (%g), dy: (%g), dz: (%g) \n",dx,dy,dz);

		// Get info from the command line:
		PetscOptionsGetInt(PETSC_NULL,"-nzbot"      ,	&nzbot 		, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL,"-nzlay"      ,	&nzlay 		, PETSC_NULL);
		PetscOptionsGetInt(PETSC_NULL,"-nztop"      ,	&nztop 		, PETSC_NULL);


		if (( nzbot + nzlay + nztop) != user->nnode_z){
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"nzbot + nzlay + nztop should be equal to nnode_z and is not; change the input file!");
		}


		dzbot = ((user->H-1.0)/2.0)/((double)(nzbot-1));
		dzlay = 1.0/((double)(nzlay-1));

		indbot = nzbot - 1;
		indlay = indbot + nzlay - 1;

		user->Setup.ind_fold_bot = indbot;
		user->Setup.ind_fold_top = indlay;

		ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
		for(iz=zs; iz<zs+zm; iz++){
			for(iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){
					coords[iz][iy][ix].z = ((double)(iz))*dz + user->z_bot;

					if (iz< indbot)				   { coords[iz][iy][ix].z =  ((double)(iz))*dzbot + user->z_bot;}
					if((iz<=indlay) && (iz>=indbot)){ coords[iz][iy][ix].z = ((double)(iz)-(double)(indbot))*dzlay + user->z_bot+((user->H-1.0)/2.0);}
					if (iz>user->Setup.ind_fold_top){ coords[iz][iy][ix].z = user->z_bot + user->H -  (user->nnode_z - iz -1)*dz; }


					if((user->ampl3D)==(0.0)){
						if(iz==indbot){ coords[iz][iy][ix].z -= user->ampl2D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)) );}
						if(iz==indlay){ coords[iz][iy][ix].z -= user->ampl2D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)));}
					}
					else{
						if(iz==indbot){ coords[iz][iy][ix].z -= user->ampl3D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)))*cos(2.0*M_PI*coords[iz][iy][ix].y/((double)(user->L)));}
						if(iz==indlay){ coords[iz][iy][ix].z -= user->ampl3D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)))*cos(2.0*M_PI*coords[iz][iy][ix].y/((double)(user->L)));}
					}

				}
			}
		}
		ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
		ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		//		ierr = VecDestroy(&gc); 									CHKERRQ(ierr);
		//		ierr = VecDestroy(&global); 								CHKERRQ(ierr);
		//		ierr = DMDestroy(&cda); 									CHKERRQ(ierr);	
		}

	else if (user->Setup.Model==7){

		dz = user->H/((double) (nz-1));
		ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);
		//		PetscPrintf(PETSC_COMM_WORLD,"dx: (%g), dy: (%g), dz: (%g) \n",dx,dy,dz);

		ierr = DMGetCoordinateDM(da,&cda); 						CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gc); 				CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,&coords); 					CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);
		for(iz=zs; iz<zs+zm; iz++){
			for(iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){
					coords[iz][iy][ix].z = ((double)(iz))*dz + user->z_bot;
				}
			}
		}
		ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
		ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		//		ierr = VecDestroy(&gc); 									CHKERRQ(ierr);
		//		ierr = VecDestroy(&global); 								CHKERRQ(ierr);
		//		ierr = DMDestroy(&cda); 									CHKERRQ(ierr);
	}

	else if (user->Setup.Model==9){
		PetscScalar 	fac, fac0, cf_rand, dif_z;
		PetscInt 		iz_Salt, iz_Base, PertType;
		PetscInt 		nnode_x,nnode_y, nnode_z;


		/* Multilayer folding setup - we need to explicitly resolve the salt/overburden interface and add random noise on that interface */
		//PetscPrintf(PETSC_COMM_WORLD," Initiating mesh for multilayer setup with detachment/overburden interface at z=%g and noise amplitude=%g \n",user->Hinterface, user->amplNoise);


		dz = user->H/((double) (nz-1));

		ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);

		ierr = DMGetCoordinateDM(da,&cda); 						CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gc); 				CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,&coords); 					CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);
		ierr = DMDAGetInfo(user->DA_Vel, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);

		//detect iz_Salt
		iz_Salt=iz_Base=0;
		for(iz=zs; iz<zs+zm; iz++){
			for(iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){
					fac0 	=	PetscAbs((coords[iz][iy][ix].z-user->Hinterface)/dz);
					if (fac0<0.5){
						iz_Salt = iz;
					}
					if (fac0>0.5 && (coords[iz][iy][ix].z < dz)){
						iz_Base = iz;
					}

				}
			}
		}

		PertType = 99;
		PetscOptionsGetInt(PETSC_NULL,"-Perturbation"		, &PertType	, PETSC_NULL);

		for(iz=zs; iz<zs+zm; iz++){
			for(iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){
					fac 	=	PetscAbs((coords[iz][iy][ix].z-user->Hinterface)/dz);
					//If Perturbation Type = 1, all the mesh above the salt/overburden interface is modified with the cos perturbation
					//
					if (PertType == 1){
						if ((fac<0.5) || ((fac > 0.5) && (coords[iz][iy][ix].z > (user->Hinterface)))){
							if ((user->amplNoise) != (0.0)) {
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand = (cf_rand-0.5);
								//PetscPrintf(PETSC_COMM_WORLD," rand=%g \n",rand);
								coords[iz][iy][ix].z = coords[iz][iy][ix].z - cf_rand*user->amplNoise;
							}
							//Modifies the Hinterface by adding a cos perturbation with amplitude defined either in 2D or 3D (ampl2D or ampl3D).
							else if((user->ampl2D)!=(0.0)){
								coords[iz][iy][ix].z = coords[iz][iy][ix].z - user->ampl2D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)));
							}
							else if((user->ampl3D)!=(0.0)){
								coords[iz][iy][ix].z = coords[iz][iy][ix].z - user->ampl3D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)))*cos(2.0*M_PI*coords[iz][iy][ix].y/((double)(user->L)));
							}
						}
					}
					//If Perturbation Type = 0, the cos perturbation is applied to the salt/overburden interface. Care must be taken not to set a too high Amplitude in this case
					// otherwise the mesh will be too deformed
					//							if (PertType == 0){
					else if (PertType == 0 || PertType == 99){

						if (fac<0.5) {
							if ((user->amplNoise) != (0.0)) {
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand = (cf_rand-0.5);
								//PetscPrintf(PETSC_COMM_WORLD," rand=%g \n",rand);
								coords[iz][iy][ix].z = user->Hinterface + cf_rand*user->amplNoise;
							}
							//Modifies the Hinterface by adding a cos perturbation with amplitude defined either in 2D or 3D (ampl2D or ampl3D).
							else if((user->ampl2D)!=(0.0)){
								coords[iz][iy][ix].z = user->Hinterface - user->ampl2D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W)));
							}
							else if((user->ampl3D)!=(0.0)){
								coords[iz][iy][ix].z = user->Hinterface + user->ampl3D*cos(2.0*M_PI*coords[iz][iy][ix].x/((double)(user->W))) *cos(2.0*M_PI*coords[iz][iy][ix].y/((double)(user->L)));
							}
						}
					}
				}
			}
		}

		// If Perturbation Type = 1, the mesh bellow the salt/overburden interface must be modified and regularized.
		if (PertType ==1){
			for(iz=zs; iz<zs+zm; iz++){
				for(iy=ys; iy<ys+ym; iy++){
					for(ix=xs; ix<xs+xm; ix++){
						fac 	=	PetscAbs((coords[iz][iy][ix].z-user->Hinterface)/dz);
						// detect the nodes that were below the salt/overburden interface using the iz level
						if ((fac>0.5) && iz < iz_Salt){
							if((user->ampl2D)!=(0.0)){
								dif_z =  (coords[iz_Salt][iy][ix].z-coords[iz_Base][iy][ix].z)/((double) (iz_Salt));
								coords[iz][iy][ix].z = ((double) (iz))*dif_z + coords[iz_Base][iy][ix].z;
							}
						}
					}
				}
			}
		}

		ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
		ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		//		ierr = VecDestroy(&global); 								CHKERRQ(ierr);
		//		ierr = VecDestroy(&gc); 									CHKERRQ(ierr);
		//		ierr = DMDestroy(&cda); 									CHKERRQ(ierr);
	}

	else{

		//	PetscPrintf(PETSC_COMM_WORLD," Setting uniform mesh coordinates \n");

		dz =     2.0/((double) (nz-1));
		ierr = DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H); CHKERRQ(ierr);





		//================================ DEBUGGING SWISS CROSS GRID FOR FDSTAG
#if 0
		{
			PetscScalar	dx, dy, dz;

			PetscPrintf(PETSC_COMM_WORLD," DEFORMING GRID FOR DEBUGGING !!!! \n");
			dx =     user->W/((double) (nx-1));
			dy =     user->L/((double) (ny-1));
			dz =     user->H/((double) (nz-1));


			ierr = DMGetCoordinateDM(da,&cda); 						CHKERRQ(ierr);
			ierr = DMGetCoordinatesLocal(da,&gc); 				CHKERRQ(ierr);
			ierr = DMDAVecGetArray(cda,gc,&coords); 					CHKERRQ(ierr);
			ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 	CHKERRQ(ierr);
			for(iz=zs; iz<zs+zm; iz++){
				for(iy=ys; iy<ys+ym; iy++){
					for(ix=xs; ix<xs+xm; ix++){
						coords[iz][iy][ix].x = ((double)(ix))*dx + user->x_left;
						//coords[iz][iy][ix].y = ((double)(iy))*dy + user->y_front;
						coords[iz][iy][ix].z = ((double)(iz))*dz + user->z_bot;


						coords[iz][iy][ix].x = ((double)(ix))*dx*((double)(ix))*dx + user->x_left;
						coords[iz][iy][ix].y = ((double)(iy))*dy*((double)(iy))*dy + user->y_front;
						coords[iz][iy][ix].z = ((double)(iz))*dz*((double)(iz))*dz + user->z_bot;


					}
				}
			}
			ierr = DMDAVecRestoreArray(cda,gc,&coords); 				CHKERRQ(ierr);
			ierr = DMGetCoordinates(da,&global); 					CHKERRQ(ierr);
			ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); 	CHKERRQ(ierr);
		}
#endif
		//================================




#if 0
		/* This part was buggy and has thus been deactivated */

		/* Create 1D vector with z-coordinates */
		ierr = PetscMalloc( (nz)*sizeof(PetscScalar), 	&zvec);
		ierr = PetscMalloc( (nz)*sizeof(PetscScalar), 	&zvec_sign);
		for (iz=0; iz<nz; iz++){
			zvec[iz] = -1.0 + ((PetscScalar) iz)*dz;
			if (zvec[iz]<0){zvec_sign[iz] 	=  -1;	}
			else{			zvec_sign[iz]	=	1;	}
		};		// regularly spaced grid between [-1..1]
		if (user->Setup.Model==1){ // single_layer folding setup. Refined grid around fold
			for (iz=0; iz<nz; iz++){	zvec[iz] = pow(PetscAbs(zvec[iz]),2); zvec[iz]=zvec[iz]*zvec_sign[iz];};				// refine around center
		}


		for (iz=0; iz<nz; iz++){	zvec[iz] = zvec[iz]*(user->H)/2.0 + (user->H)/2.0 + user->z_bot; };	// rescale to orig. length


		if ( (user->Setup.Model==4)  || (user->Setup.Model==9) ){ // refined grid towards the top if we have a detachment folding setup

			/*
			dz =     1.0/((double) (nz-1));
			for (iz=0; iz<nz; iz++){
				zvec[iz] = -1 + ((PetscScalar) iz)*dz;
			};		// regularly spaced grid between [0..1]

			//for (iz=0; iz<nz; iz++){	zvec[iz] = -1*pow(PetscAbs(zvec[iz]),2.0) ; };	// refine towards top

			// rescale back to real coordinates
			for (iz=0; iz<nz; iz++){	zvec[iz] = zvec[iz]*(user->H)/1.0 + (user->H)/1.0 + user->z_bot; };			// rescale to orig. length
			 */
		}


		/* Find interface height and match grid in z-direction to this interface */
		diff = 100; diff1=100;
		step = 1;
		if (__ELEMENT_TYPE__ == ELEMENT_Q2P1){
			step = 2;		// Q2 elements should make the element EDGES aligned with the
		}
		for (iz=0; iz<nz; iz=iz+step){

			if ( (user->Setup.Model==0) || (user->Setup.Model==4) ){ // diapir or multilayer det. folding
				if (PetscAbs(zvec[iz]-user->Setup.Diapir_Hi)<diff){
					diff 			= PetscAbs(zvec[iz]-user->Setup.Diapir_Hi);
					ind_Hi_diapir 	= iz;
				}
			}
			if (user->Setup.Model==1){ // single_layer folding setup
				if (PetscAbs(zvec[iz]-0.5*(-user->Setup.SingleFold_H))<diff){
					diff 			= PetscAbs(zvec[iz]-0.5*(-user->Setup.SingleFold_H));
					ind_fold_bot 	= iz;
				}
				if (PetscAbs(zvec[iz]-0.5*(user->Setup.SingleFold_H))<diff1){
					diff1 			= PetscAbs(zvec[iz]-0.5*(user->Setup.SingleFold_H));
					ind_fold_top 	= iz;
				}
			}

		};

		if ( (user->Setup.Model==0) || (user->Setup.Model==4) ){
			PetscPrintf(PETSC_COMM_WORLD,"ind_Hi_diapir=%lld \n",(LLD)ind_Hi_diapir);
			zvec[ind_Hi_diapir]	=	user->Setup.Diapir_Hi;
		}
		else if  (user->Setup.Model==1){
			zvec[ind_fold_bot]	=	0.5*(-user->Setup.SingleFold_H);
			zvec[ind_fold_top]	=	0.5*( user->Setup.SingleFold_H);
		}


		/* The Zagros multilayer folding setup has adjustable layer thickness (from command line options) */
		if (user->Setup.Model==9){
			{
				PetscScalar LayerHeight;
				PetscInt layernumber, ind_layer;

				for (layernumber=1; layernumber<5; layernumber++){
					// Extract info from command line
					LayerHeight = PETSC_NULL;
					if (layernumber==1){
						PetscOptionsGetReal(PETSC_NULL,"-Layer1_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}
					else if (layernumber==2){
						PetscOptionsGetReal(PETSC_NULL,"-Layer2_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}
					else if (layernumber==3){
						PetscOptionsGetReal(PETSC_NULL,"-Layer3_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}
					else if (layernumber==4){
						PetscOptionsGetReal(PETSC_NULL,"-Layer4_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}
					else if (layernumber==5){
						PetscOptionsGetReal(PETSC_NULL,"-Layer5_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}
					else if (layernumber==6){
						PetscOptionsGetReal(PETSC_NULL,"-Layer6_bottom"      ,	&LayerHeight 		, PETSC_NULL);
					}

					if (LayerHeight != PETSC_NULL){

						PetscPrintf(PETSC_COMM_WORLD,"# Layer %lld has lower interface at height %g \n", (LLD)layernumber, LayerHeight);
						diff = 100;
						for (iz=0; iz<nz; iz=iz+step){
							if (PetscAbs(zvec[iz]-LayerHeight)<diff){
								diff 		= PetscAbs(zvec[iz]-LayerHeight);
								ind_layer 	= iz;
							}
						};
						zvec[ind_layer] = LayerHeight;

						if (layernumber==1){
							// In this setup, the first layer indicates the location of the salt-overburden interface
							ind_Hi_diapir 			= ind_layer;
							user->Setup.Diapir_Hi 	= LayerHeight;
						}
					}
				}
			}
		}

		/* Align internal nodes in case of Q2 element*/
		if (__ELEMENT_TYPE__ == ELEMENT_Q2P1){
			for (iz=1; iz<nz-1; iz=iz+step){
				zvec[iz] = (zvec[iz-1] + zvec[iz+1])/2.0;
			}
		}


		if ((user->Setup.Model==0) || (user->Setup.Model==4) || (user->Setup.Model==9)){  	// diapir or multilayer detachment folding setup

			if (user->Setup.Model==9){
				PetscPrintf(PETSC_COMM_WORLD,"# Zagros multilayer detachment folding setup with adjustable layer thickness (from command line) \n");
			}
			else if (user->Setup.Model==4){
				PetscPrintf(PETSC_COMM_WORLD,"# Multilayer detachment folding setup with fixed layers \n");
			}
			else if (user->Setup.Model==0){
				PetscPrintf(PETSC_COMM_WORLD,"# Diapir two-layer setup \n");
			}

			//zvec[ind_Hi_diapir] = user->Setup.Diapir_Hi + user->z_bot;
			user->Setup.ind_Hi_diapir = ind_Hi_diapir;
		}
		if (user->Setup.Model==1){ 		// single layer folding setup
			PetscPrintf(PETSC_COMM_WORLD,"# Single layer folding setup \n");

			zvec[ind_fold_bot]  = -0.5*(user->Setup.SingleFold_H);
			zvec[ind_fold_top]  =  0.5*(user->Setup.SingleFold_H);
			user->Setup.ind_fold_bot = ind_fold_bot;
			user->Setup.ind_fold_top = ind_fold_top;
		}


		ierr = DMGetCoordinateDM(da,			&cda); CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,		&gc); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,			&coords); CHKERRQ(ierr);
		ierr = DMDAGetCorners(da,&xs_noghost,&ys_noghost, PETSC_NULL, &xm_noghost,&ym_noghost, PETSC_NULL); CHKERRQ(ierr); // required for adding noise
		ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

		for (iz=zs; iz<zs+zm; iz++){
			for (iy=ys; iy<ys+ym; iy++){
				for(ix=xs; ix<xs+xm; ix++){

					//					coords[iz][iy][ix].z = ((double) (iz))*dz + user->z_bot;
					coords[iz][iy][ix].z = zvec[iz];

					// Deform the grid in a sinusoidal manner

					fac	=	0;
					if ((user->Setup.Model==0) || (user->Setup.Model==4) || (user->Setup.Model==9)){		// diapir or det. fold
						if  (iz<ind_Hi_diapir)	{ 	fac = (double)(iz)/(double)(ind_Hi_diapir);	    	}
						if  (iz==ind_Hi_diapir)	{ 	fac = 1;  };
						if  (iz>ind_Hi_diapir)  {	fac = (double)(nz-iz-1)/((double)(nz-ind_Hi_diapir-1));	}
					}
					else if (user->Setup.Model==1){	// single-layer fold
						if  (iz <ind_fold_bot)  { 	fac = (double)(iz)/(double)(ind_fold_bot);	    	}
						if  (iz >ind_fold_top)  {	fac = (double)(nz-iz)/((double)(nz-ind_fold_top));	}
						if  (iz==ind_fold_bot)	{   fac = 1;											}
						if  (iz==ind_fold_top)	{   fac = 1;											}
					}
					if (iz==0 || (iz==nz-1)){	fac=0; }
					//	if (iz==0){	fac=0; }


					coords[iz][iy][ix].z = zvec[iz];

					coords[iz][iy][ix].z = zvec[iz]
					                            +	  fac*user->ampl2D*cos(2*coords[iz][iy][ix].x/(user->W)*M_PI)
					+	0*fac*user->ampl2D*cos(2*coords[iz][iy][ix].y/(user->L)*M_PI)
					+	  fac*user->ampl3D*cos(2*coords[iz][iy][ix].x/(user->W)*M_PI)*cos(2*coords[iz][iy][ix].y/(user->L)*M_PI);

					//		coords[iz][iy][ix].x = coords[iz][iy][ix].x
					//			+	fac*user->ampl2D*cos(2*coords[iz][iy][ix].y/(user->L)*M_PI);

					//	coords[iz][iy][ix].y = coords[iz][iy][ix].y
					//		+	fac*user->ampl2D*cos(2*coords[iz][iy][ix].x/(user->W)*M_PI);


					if ((PetscAbs(user->amplNoise)>0.0) && (fac==1)){
						ierr = PetscRandomGetValueReal(rctx, &rand); CHKERRQ(ierr);
						rand=rand-0.5;

						// only set random noise in points that are NOT ghost points!
						if ((iy>ys_noghost) && (iy<ys_noghost+ym_noghost) && (ix>xs_noghost) && (ix<xs_noghost+xm_noghost)){
							coords[iz][iy][ix].z = coords[iz][iy][ix].z + rand*(user->amplNoise);
						}
					}

					if (iz == nz-1){

						// Add random noise to free surface in pts that are NOT ghost points
						if ((iy>ys_noghost) && (iy<ys_noghost+ym_noghost) && (ix>xs_noghost) && (ix<xs_noghost+xm_noghost)){
							ierr = PetscRandomGetValueReal(rctx, &rand); CHKERRQ(ierr);
							rand=rand-0.5;
							coords[iz][iy][ix].z = coords[iz][iy][ix].z + rand*(user->SurfaceNoiseAmplitude);

						}

						// Put an inclined slope to the top surface
						coords[iz][iy][ix].z = coords[iz][iy][ix].z + (coords[iz][iy][ix].x-user->x_left)*tan(user->SurfaceAngle/180*M_PI );

					}


				}
			}
		}

		ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);

		/*=========================================================================================
		 * Ensure that also ghost coordinates are send to all other procs.
		 * This is particularly an issue if random noise is set as initial perturbation.
		 */
		ierr = VecDuplicate( gc, &new_gc_coords ); CHKERRQ(ierr);
		ierr = VecCopy( gc, new_gc_coords ); CHKERRQ(ierr);

		//ierr = VecDestroy(gc); CHKERRQ(ierr);
		//ierr = DMDestroy(cda); CHKERRQ(ierr);

		ierr = DASetGhostedCoordinates(da,new_gc_coords); CHKERRQ(ierr);
		ierr = VecDestroy(&new_gc_coords); CHKERRQ(ierr);
		/*=========================================================================================*/

		// DEBUGGING
		/* Not updated crap
		{
			PetscMPIInt rank;
					Vec my_gc;
			DMDACoor3d ***coords;

					MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

			DMGetCoordinateDM(da,&cda);
					DMGetCoordinatesLocal(da,&my_gc);

				DMDAVecGetArray(cda,my_gc,                   &coords);
				DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
			if(rank==1) {
					for (iz=zs; iz<zs+zm; iz++){
						for (iy=ys; iy<ys+ym; iy++){
							for(ix=xs; ix<xs+xm; ix++){
					printf("[%lld]:  (%lld,%lld,%lld) : %lf,%lf,%lf \n", rank, ix,iy,iz, coords[iz][iy][ix].x, coords[iz][iy][ix].y,coords[iz][iy][ix].z );

				}}}
			}
			DMDAVecRestoreArray(cda,my_gc,&coords);
			DMDestroy(cda);
			VecDestroy(my_gc);
		}
		 */
		ierr = PetscFree(zvec); CHKERRQ(ierr);
		ierr = PetscFree(zvec_sign); CHKERRQ(ierr);

#endif
	}

	if( user->NodesDistributionCPUsFineGrid == PETSC_NULL ) {
		ierr = PetscFree(user->NodesDistributionCPUsFineGrid ); CHKERRQ(ierr);
	}

	ierr = PetscFree(NodesDistributionGlobal); CHKERRQ(ierr);
	if ( user->NodesDistributionCPUsFineGrid != PETSC_NULL ) {
		ierr = PetscFree( user->NodesDistributionCPUsFineGrid ); CHKERRQ(ierr);
	}


	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);         // Destroy random context



	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Advect and update internal free surface
 *
 * In case we use a 'sticky-air' approach with an internal free surface, we need to advect the free surface
 * at each timestep, and interpolate it back to the regular computational grid, prior to doing things like
 * erosion or sedimentation to the internal free surface.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "AdvectAndUpdate_InternalFreeSurface"
PetscErrorCode AdvectAndUpdate_InternalFreeSurface( UserContext *user, Vec Velocity )
{
	PetscMPIInt     rank;
	PetscErrorCode	ierr;
	PetscScalar		***LocalSurfaceTopography;
	DM				cda, cda_SurfaceTopo;
	DMDACoor3d		***coors, ***coors_SurfaceTopo;
	Vec 			gc_SurfaceTopo, gc;
	Vec 			local_Vel;
	PetscScalar		*Local_Surface_Array_x,	*Local_Surface_Array_y,	*Local_Surface_Array_z;
	PetscScalar		*Local_Surface_Array_x_summed,*Local_Surface_Array_y_summed,*Local_Surface_Array_z_summed;
	PetscScalar		*Local_Surface_Array_Vx,	*Local_Surface_Array_Vy,	*Local_Surface_Array_Vz;
	PetscScalar		*Local_Surface_Array_Vx_summed,*Local_Surface_Array_Vy_summed,*Local_Surface_Array_Vz_summed;
	PetscScalar		z_min, z_max, gmin[3], gmax[3];
	Field 			***velocity;
	PetscInt		xs,ys,zs,xm,ym,zm,xs_Z,ys_Z,zs_Z,xm_Z,ym_Z,zm_Z,ix,iy,num,nnode_x,nnode_y, NumProcs_z, zs_g, zm_g;
	MPI_Comm 		comm;
	Vec				Surface_Deformed_Global_x,	Surface_Deformed_Global_y,	Surface_Deformed_Global_z;
	Vec				Surface_Deformed_Local_x,	Surface_Deformed_Local_y,	Surface_Deformed_Local_z;
	PetscScalar		***Surface_Deformed_x,		***Surface_Deformed_y,		***Surface_Deformed_z;
	PetscScalar		***SurfaceVelocity_Vx,		***SurfaceVelocity_Vy,		***SurfaceVelocity_Vz;
	PetscScalar	 	xmin_advected, xmax_advected, ymin_advected, ymax_advected;

	PetscLogDouble	cputime_start, cputime_end;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscTime(&cputime_start);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"# Advecting internal free surface ... ");

	/* The outline of this routine is as follows:
	 *
	 * 1) Create local x,y,z arrays that contain the local non-ghosted portions of the free surface mesh.
	 * 2) If part of a free surface happens to be on the current processor, advect it.
	 * 	  If not, put all coordinates to zero.
	 * 3) Broadcast the vectors with x,y,z coordinates to all processors in the z-direction. [we need a list of this]
	 * 4) Redundantly, all processors interpolate the deformed free surface back onto the regular grid specified by the DA.
	 * 5) Put back the topo information into the global DM vector.
	 *
	 * We make a few assumptions in the routine, namely:
	 * 	-	the x,y coordinates do not change in the vertical direction (with other words: it will work for an FDSTAG
	 * 		formulation or for an eulerian FE formulation, but not in combination with a lagrangian mesh)
	 */
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,				&gc_SurfaceTopo			); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);


	/* Fetch velocity & coordinates of the velocity grid */
	ierr = DMCreateLocalVector(user->DA_Vel,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->DA_Vel, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_Vel,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_Vel, local_Vel,	&velocity); CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(user->DA_Vel, &cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_Vel, &gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coors); CHKERRQ(ierr);

	/* Extract minimum and maximum z-coordinates of the non-ghosted local portion of the velocity grid */
	ierr = DMDAGetCorners(user->DA_Vel,&xs,&ys,&zs,PETSC_NULL,PETSC_NULL,&zm); CHKERRQ(ierr);
	z_min = coors[zs   ][ys][xs].z;
	if ((zs+zm)< (user->finest_nnode_z)){
		z_max = coors[zs+zm][ys][xs].z;
	}
	else{
		z_max = coors[zs+zm-1][ys][xs].z;
	}

	ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);	// surface topography
	ierr = DMDAGetInfo(user->DA_SurfaceTopography,0,&nnode_x,&nnode_y,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);	// surface topography

	ierr = DMDAGetGhostCorners(user->DA_Vel,PETSC_NULL,PETSC_NULL,&zs_g,PETSC_NULL,PETSC_NULL,&zm_g); CHKERRQ(ierr);

	DMDAGetBoundingBox(user->DA_Vel,gmin,gmax);	// global bounding box
/*
	x_min = gmin[0];	x_max = gmax[0];
	y_min = gmin[1];	y_max = gmax[1];
*/


	/* Computed the boundaries of the FDSTAG grid AFTER advection by the background strain rate */
	xmin_advected 	=	user->x_left 			;
	xmax_advected 	=	user->x_left 	+ user->W;
	ymin_advected 	=	user->y_front			;
	ymax_advected 	=	user->y_front 	+ user->L;
	xmin_advected	=	xmin_advected - user->dt*xmin_advected*user->BC.Exx;
	xmax_advected	=	xmax_advected - user->dt*xmax_advected*user->BC.Exx;
	ymin_advected	=	ymin_advected - user->dt*ymin_advected*user->BC.Eyy;
	ymax_advected	=	ymax_advected - user->dt*ymax_advected*user->BC.Eyy;





	/* Create local Vec's to store x,y,z data */
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_x);			CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_y);			CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_z);			CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vx);			CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vy);			CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vz);			CHKERRQ(ierr);


	/* 1) Extract local x,y,z coordinates of surface DA
	 * 2) Check if they can be advected or not [are we on the correct PROC for this?]
	 *    Advect the ones that can be advected, set the others to zero
	 */
	num = 0;
	for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
		for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){

			PetscScalar x,y,z, fac, Vx, Vy, Vz;
			PetscInt 	iz;

			/* Extract x,y,z coordinates of surface topography */
			x = coors_SurfaceTopo[zs_Z][iy][ix].x;
			y = coors_SurfaceTopo[zs_Z][iy][ix].y;
			z = LocalSurfaceTopography[zs_Z][iy][ix];	// topography


			/* Can we advect? */
			if ((z >= z_min) && (z < z_max)){
				/* Yes -> interpolate velocity and advect */
				PetscScalar dz ;

				/* Find the correct z-level [works for variable z-spacing] */
				iz 	= zs_g;
				while ((coors[iz][ys][xs].z < z)){
					iz=iz+1;
				}
				if (iz>zs_g){iz = iz-1;}


				dz  = 	(coors[iz+1][ys][xs].z-coors[iz  ][ys][xs].z);

				fac	=	1.0-(z-coors[iz][ys][xs].z)/dz;		// linear interpolation factor

				// Interpolate velocity in the z-direction
				Vx 	=	fac*velocity[iz][iy][ix].Vx + (1-fac)*velocity[iz+1][iy][ix].Vx;
				Vy 	=	fac*velocity[iz][iy][ix].Vy + (1-fac)*velocity[iz+1][iy][ix].Vy;
				Vz 	=	fac*velocity[iz][iy][ix].Vz + (1-fac)*velocity[iz+1][iy][ix].Vz;

				// Advect
				x 	=	x +	Vx*user->dt;
				y 	=	y +	Vy*user->dt;
				z 	=	z +	Vz*user->dt;

				if (z>z_max){	z = z_max;}	// above bounding box
				if (z<z_min){	z = z_min;}	// below bounding box


			}
			else{
				/* No, can't advect-> set values to zero */
				if (z >= user->z_bot + user->H){
					z = user->z_bot + user->H;
				}
				else{
					z=0;
				}
				x=0; 	y=0;
				Vx=0; 	Vy=0;	Vz=0;

		//		ierr = PetscPrintf(PETSC_COMM_SELF," rank %i [xold,yold,zold]=[%f,%f,%f] [x,y,z]=[%f,%f,%f] z_min=%f, z_max=%f  delta_zmax=%f dz=%f dt=%f ix,iy=[%i,%i]\n",rank, coors_SurfaceTopo[zs_Z][iy][ix].x, coors_SurfaceTopo[zs_Z][iy][ix].y,LocalSurfaceTopography[zs_Z][iy][ix],x,y,z,z_min,z_max, user->dt*z_max*(user->BC.Exx + user->BC.Eyy), Vz*user->dt,user->dt,ix,iy);

			}



			/* Add (advected or zeroed) coords to local 1D array [used later with MPI] */
			Local_Surface_Array_x[num] = x;
			Local_Surface_Array_y[num] = y;
			Local_Surface_Array_z[num] = z;

			Local_Surface_Array_Vx[num] = Vx;
			Local_Surface_Array_Vy[num] = Vy;
			Local_Surface_Array_Vz[num] = Vz;


			num=num+1;
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);


	/* 3) Broadcast the vectors with the (advected) internal free surface coordinates to the PROCS above and below
	 * 	For this to work we:
	 * 		- Find all the PROCS above and below
	 * 		- Broadcast only to those procs
	 *
	 * 		- Create a global Vec to store the x,y,z coordinates of the deformed surface, such that we later have access to the ghost values.
	 */

	/* Create a new MPI sub-communicator that contains all PROCS within the current collumn */
	DAGetProcessorSubset_VerticalDirection(user->DA_Vel,&NumProcs_z, &comm);

	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_x_summed);	CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_y_summed);	CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_z_summed);	CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vx_summed);	CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vy_summed);	CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm_Z*ym_Z)*sizeof(PetscScalar),&Local_Surface_Array_Vz_summed);	CHKERRQ(ierr);

	if (NumProcs_z>1){
		/* Send data to other processors in current collumn */
		ierr = MPI_Allreduce(Local_Surface_Array_x,	Local_Surface_Array_x_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);
		ierr = MPI_Allreduce(Local_Surface_Array_y,	Local_Surface_Array_y_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);
		ierr = MPI_Allreduce(Local_Surface_Array_z,	Local_Surface_Array_z_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);

		ierr = MPI_Allreduce(Local_Surface_Array_Vx,	Local_Surface_Array_Vx_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);
		ierr = MPI_Allreduce(Local_Surface_Array_Vy,	Local_Surface_Array_Vy_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);
		ierr = MPI_Allreduce(Local_Surface_Array_Vz,	Local_Surface_Array_Vz_summed,	xm_Z*ym_Z, MPIU_SCALAR, MPIU_SUM, comm); CHKERRQ(ierr);
	}
	MPI_Comm_free(&comm);	// free communicator along collumns

	/* Create global vector & update ghost points of the deformed surface grid */
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&Surface_Deformed_Global_x); 	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&Surface_Deformed_Global_y); 	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&Surface_Deformed_Global_z); 	CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->DA_SurfaceTopography,		&Surface_Deformed_Local_x); 	CHKERRQ(ierr);
	ierr = DMGetLocalVector(user->DA_SurfaceTopography,		&Surface_Deformed_Local_y); 	CHKERRQ(ierr);
	ierr = DMGetLocalVector(user->DA_SurfaceTopography,		&Surface_Deformed_Local_z); 	CHKERRQ(ierr);

	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_x,			&Surface_Deformed_x); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_y,			&Surface_Deformed_y); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_z,			&Surface_Deformed_z); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,			user->SurfaceTopography_Vx,			&SurfaceVelocity_Vx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,			user->SurfaceTopography_Vy,			&SurfaceVelocity_Vy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,			user->SurfaceTopography_Vz,			&SurfaceVelocity_Vz); CHKERRQ(ierr);




	num=0;
	for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
		for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){
			PetscScalar 	eps=1e-6*(xmax_advected-xmin_advected);

			if (NumProcs_z>1){
				Surface_Deformed_x[zs_Z][iy][ix] 	=	Local_Surface_Array_x_summed[num];
				Surface_Deformed_y[zs_Z][iy][ix] 	=	Local_Surface_Array_y_summed[num];
				Surface_Deformed_z[zs_Z][iy][ix] 	=	Local_Surface_Array_z_summed[num];

				SurfaceVelocity_Vx[zs_Z][iy][ix]	=	Local_Surface_Array_Vx_summed[num];
				SurfaceVelocity_Vy[zs_Z][iy][ix]	=	Local_Surface_Array_Vy_summed[num];
				SurfaceVelocity_Vz[zs_Z][iy][ix]	=	Local_Surface_Array_Vz_summed[num];
			}
			else{
				Surface_Deformed_x[zs_Z][iy][ix] 	=	Local_Surface_Array_x[num];
				Surface_Deformed_y[zs_Z][iy][ix] 	=	Local_Surface_Array_y[num];
				Surface_Deformed_z[zs_Z][iy][ix] 	=	Local_Surface_Array_z[num];

				SurfaceVelocity_Vx[zs_Z][iy][ix]	=	Local_Surface_Array_Vx[num];
				SurfaceVelocity_Vy[zs_Z][iy][ix]	=	Local_Surface_Array_Vy[num];
				SurfaceVelocity_Vz[zs_Z][iy][ix]	=	Local_Surface_Array_Vz[num];

			}

			/* Make sure that the deformed surface has xy boundaries that are always slightly outside the computational domain - interpolation
			 * back to a regular grid will result in sane results
			 * */
			if ((ix==0)){// && (Surface_Deformed_x[zs_Z][iy][ix]>x_min)){					// point at left boundary is inside box
				Surface_Deformed_x[zs_Z][iy][ix] = xmin_advected-eps; 	// move point to outside regular mesh
			}
			if ((ix==(nnode_x-1))){// && (Surface_Deformed_x[zs_Z][iy][ix]<x_max)){			// point at right boundary is inside box
				Surface_Deformed_x[zs_Z][iy][ix] = xmax_advected+eps;		// move point to outside regular mesh
			}
			if ((iy==0)){// && (Surface_Deformed_y[zs_Z][iy][ix]>y_min)){					// point at left boundary is inside box
				Surface_Deformed_y[zs_Z][iy][ix] = ymin_advected-eps; 	// move point to outside regular mesh
			}

			if ((iy==(nnode_y-1))){// && (Surface_Deformed_y[zs_Z][iy][ix]<y_max)){		// point at right boundary is inside box
				Surface_Deformed_y[zs_Z][iy][ix] = ymax_advected+eps;	// move point to outside regular mesh
			}


			num 									=	num+1;
		}
	}

	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_x,			&Surface_Deformed_x); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_y,			&Surface_Deformed_y); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		Surface_Deformed_Local_z,			&Surface_Deformed_z); CHKERRQ(ierr);

	/* Restore */
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		user->SurfaceTopography_Vx,			&SurfaceVelocity_Vx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		user->SurfaceTopography_Vy,			&SurfaceVelocity_Vy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,		user->SurfaceTopography_Vz,			&SurfaceVelocity_Vz); CHKERRQ(ierr);



	/* No need for this anymore */
	ierr 	= 	PetscFree(Local_Surface_Array_x); 			CHKERRQ(ierr);
	ierr 	= 	PetscFree(Local_Surface_Array_y); 			CHKERRQ(ierr);
	ierr 	= 	PetscFree(Local_Surface_Array_z); 			CHKERRQ(ierr);

	ierr 	= 	PetscFree(Local_Surface_Array_x_summed); 	CHKERRQ(ierr);
	ierr 	=	PetscFree(Local_Surface_Array_y_summed); 	CHKERRQ(ierr);
	ierr 	=	PetscFree(Local_Surface_Array_z_summed); 	CHKERRQ(ierr);

	ierr 	= 	PetscFree(Local_Surface_Array_Vx); 			CHKERRQ(ierr);
	ierr 	= 	PetscFree(Local_Surface_Array_Vy); 			CHKERRQ(ierr);
	ierr 	= 	PetscFree(Local_Surface_Array_Vz); 			CHKERRQ(ierr);

	ierr 	= 	PetscFree(Local_Surface_Array_Vx_summed); 	CHKERRQ(ierr);
	ierr 	=	PetscFree(Local_Surface_Array_Vy_summed); 	CHKERRQ(ierr);
	ierr 	=	PetscFree(Local_Surface_Array_Vz_summed); 	CHKERRQ(ierr);


	/* Set back coordinates of surface_topography grid (as we might redefine them in the next step to correct for BG deformation rate)*/
	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);

//	ierr 	= 	VecDestroy(gc_SurfaceTopo); 				CHKERRQ(ierr);
//	ierr 	=	DMDestroy(cda_SurfaceTopo); 				CHKERRQ(ierr);

	// Copy local portions to global vectors - at this stage, the deformed surface grid is stored in Global vectors
	ierr	=	DMLocalToGlobalBegin(user->DA_SurfaceTopography,Surface_Deformed_Local_x,INSERT_VALUES,Surface_Deformed_Global_x);	CHKERRQ(ierr);
	ierr	=	DMLocalToGlobalEnd  (user->DA_SurfaceTopography,Surface_Deformed_Local_x,INSERT_VALUES,Surface_Deformed_Global_x);	CHKERRQ(ierr);
	
	ierr	=	DMLocalToGlobalBegin(user->DA_SurfaceTopography,Surface_Deformed_Local_y,INSERT_VALUES,Surface_Deformed_Global_y);	CHKERRQ(ierr);
	ierr	=	DMLocalToGlobalEnd  (user->DA_SurfaceTopography,Surface_Deformed_Local_y,INSERT_VALUES,Surface_Deformed_Global_y);	CHKERRQ(ierr);

	ierr	=	DMLocalToGlobalBegin(user->DA_SurfaceTopography,Surface_Deformed_Local_z,INSERT_VALUES,Surface_Deformed_Global_z);	CHKERRQ(ierr);
	ierr	=	DMLocalToGlobalEnd  (user->DA_SurfaceTopography,Surface_Deformed_Local_z,INSERT_VALUES,Surface_Deformed_Global_z);	CHKERRQ(ierr);


	/* 4) Interpolate the data of the deformed free surface grid back on to regular (FDSTAG or FE) grid
	 *
	 */

	/* - Recreate the x,y coordinates in case the FDSTAG mesh is being deformed with a constant BG strainrate
	 * - Extend the surface of the deformed grid at the boundaries such that it always is larger than the FDSTAG grid
	 * - Loop over all surface elements of the deformed surface including ghost points.
	 * - Divide each square in 2 triangles.
	 * 	- Loop over the regular grid and detect if a regular point is inside the triangle.
	 * 	- If it is interpolate the height of the surface
	 * - Set back the regular grid points to the global vector
	 *
	 */
	{
		PetscScalar xmin,xmax,ymin,ymax,zmin,zmax;


		xmin 	=	user->x_left 			;
		xmax 	=	user->x_left 	+ user->W;
		ymin 	=	user->y_front			;
		ymax 	=	user->y_front 	+ user->L;
		zmin 	=	user->z_bot				;
		zmax 	=	user->z_bot	 	+ user->H;

		// advect borders of grid based on BG strain rate
		xmin	=	xmin - user->dt*xmin*user->BC.Exx;
		xmax	=	xmax - user->dt*xmax*user->BC.Exx;
		ymin	=	ymin - user->dt*ymin*user->BC.Eyy;
		ymax	=	ymax - user->dt*ymax*user->BC.Eyy;
		zmin	=	zmin + user->dt*zmin*(user->BC.Exx + user->BC.Eyy);
		zmax	=	zmax + user->dt*zmax*(user->BC.Exx + user->BC.Eyy);

		// set new uniform coordinates
		ierr = DMDASetUniformCoordinates(user->DA_SurfaceTopography, xmin,xmax,ymin,ymax,zmin,zmax);	CHKERRQ(ierr);


		ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,					&gc_SurfaceTopo			); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);

	}



	/* Copy global coords to local vector, including ghost values */
	ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, Surface_Deformed_Global_x, INSERT_VALUES, Surface_Deformed_Local_x); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   Surface_Deformed_Global_x, INSERT_VALUES, Surface_Deformed_Local_x); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, 		Surface_Deformed_Local_x,	&Surface_Deformed_x); CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, Surface_Deformed_Global_y, INSERT_VALUES, Surface_Deformed_Local_y); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   Surface_Deformed_Global_y, INSERT_VALUES, Surface_Deformed_Local_y); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, 		Surface_Deformed_Local_y,	&Surface_Deformed_y); CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, Surface_Deformed_Global_z, INSERT_VALUES, Surface_Deformed_Local_z); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   Surface_Deformed_Global_z, INSERT_VALUES, Surface_Deformed_Local_z); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, 		Surface_Deformed_Local_z,	&Surface_Deformed_z); CHKERRQ(ierr);

	ierr = DMDAGetGhostCorners(user->DA_Vel,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);


	num = 0;
	/* Loop over deformed surface [including ghost points]*/
	/* Loop over regular surface */
	{
		PetscInt iix,iiy;
		for (iiy=ys_Z; iiy<ys_Z+ym_Z; iiy++){
			for(iix=xs_Z; iix<xs_Z+xm_Z; iix++){
				LocalSurfaceTopography[zs_Z][iiy][iix] = user->z_bot	 	+ user->H;
			}
		}
	}

	for (iy=ys; iy<ys+ym-1; iy++){
		for(ix=xs; ix<xs+xm-1; ix++){
			PetscInt		iix,iiy;
			PetscScalar		p_x[5],  p_y[5],  p_z[5];

			/* Extract a contour of the linear element and use an isinside routine to find out if the point is inside the element */
			p_x[0] = Surface_Deformed_x[zs_Z][iy  ][ix  ];	p_y[0] = Surface_Deformed_y[zs_Z][iy  ][ix  ]; p_z[0] = Surface_Deformed_z[zs_Z][iy  ][ix  ];
			p_x[1] = Surface_Deformed_x[zs_Z][iy  ][ix+1];	p_y[1] = Surface_Deformed_y[zs_Z][iy  ][ix+1]; p_z[1] = Surface_Deformed_z[zs_Z][iy  ][ix+1];
			p_x[2] = Surface_Deformed_x[zs_Z][iy+1][ix+1];	p_y[2] = Surface_Deformed_y[zs_Z][iy+1][ix+1]; p_z[2] = Surface_Deformed_z[zs_Z][iy+1][ix+1];
			p_x[3] = Surface_Deformed_x[zs_Z][iy+1][ix  ];	p_y[3] = Surface_Deformed_y[zs_Z][iy+1][ix  ]; p_z[3] = Surface_Deformed_z[zs_Z][iy+1][ix  ];
			p_x[4] = p_x[0]; 								p_y[4] = p_y[0]; 							   p_z[4] = p_z[0];

			/* Loop over regular surface */
			for (iiy=ys_Z; iiy<ys_Z+ym_Z; iiy++){
				for(iix=xs_Z; iix<xs_Z+xm_Z; iix++){
					PetscScalar 	x,y,z_interpolated;
					PetscBool		inside;

					/* Regular point */
					x = coors_SurfaceTopo[zs_Z][iiy][iix].x;
					y = coors_SurfaceTopo[zs_Z][iiy][iix].y;

					/* Is point inside square polygon? */
					inside = pnpoly(4, p_x, p_y, x,y);
					if (inside){
						z_interpolated =  InterpolateWithin2DLinearElement(p_x, p_y, p_z, x, y);		// interpolate surface coordinates, using linear shape function


						LocalSurfaceTopography[zs_Z][iiy][iix] = z_interpolated;


					//	PetscPrintf(PETSC_COMM_WORLD,"rank=%i; Point [%f,%f] is inside element z_interpolated=%f Surface_Deformed_x=[%f %f %f %f]; Surface_Deformed_y=[%f %f %f %f]; Surface_Deformed_z=[%f %f %f %f]\n",rank,x,y,z_interpolated,p_x[0],p_x[1],p_x[2],p_x[3],p_y[0],p_y[1],p_y[2],p_y[3],p_z[0],p_z[1],p_z[2],p_z[3]);

						coors_SurfaceTopo[zs_Z][iiy][iix].z = z_interpolated;
					}


				}

			}
		}

	}
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_Deformed_Local_x,	&Surface_Deformed_x); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_Deformed_Local_y,	&Surface_Deformed_y); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_Deformed_Local_z,	&Surface_Deformed_z); CHKERRQ(ierr);


	ierr 	= 	DMDAVecRestoreArray(user->DA_Vel,local_Vel,&velocity); CHKERRQ(ierr);
//	ierr 	= 	DMRestoreLocalVector(user->DA_Vel,&local_Vel); 		CHKERRQ(ierr);
	ierr 	= 	VecDestroy(&local_Vel); 		CHKERRQ(ierr);


	ierr 	= 	DMDAVecRestoreArray(cda,gc,&coors); 	CHKERRQ(ierr);

//	ierr 	= 	VecDestroy(gc); 					CHKERRQ(ierr);
//	ierr 	= 	DMDestroy( cda ); 					CHKERRQ(ierr);


	/* Set back coordinates */
	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);

	ierr 	= 	DASetCoordinatesFromLocalVector( user->DA_SurfaceTopography, gc_SurfaceTopo ); 						CHKERRQ(ierr);

	/* Cleaning up */
//	ierr 	= 	VecDestroy(gc_SurfaceTopo); 				CHKERRQ(ierr);
//	ierr 	=	DMDestroy(cda_SurfaceTopo); 				CHKERRQ(ierr);

	DMRestoreLocalVector(user->DA_SurfaceTopography, &Surface_Deformed_Local_x);
	DMRestoreLocalVector(user->DA_SurfaceTopography, &Surface_Deformed_Local_y);
	DMRestoreLocalVector(user->DA_SurfaceTopography, &Surface_Deformed_Local_z);

	VecDestroy(&Surface_Deformed_Global_x);
	VecDestroy(&Surface_Deformed_Global_y);
	VecDestroy(&Surface_Deformed_Global_z);

	/* Print time */
	PetscTime(&cputime_end);

	ierr = PetscPrintf(PETSC_COMM_WORLD," [%f s] \n",cputime_end-cputime_start);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Update surface and bottom topography, in case we use a deformable FE mesh and the topography coincides
 * with the top & bottom of the model domain.
 *
 * Every processor will get info about the surface and the bottom, which is used to
 *  (1) Remesh, based on bottom and surface topography
 *  (2) Postprocessing, which makes it easier to visualize and track topography
 *  (3) Perform 'realistic' erosion models.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "UpdateSurfaceAndBottomTopography_FEM"
PetscErrorCode UpdateSurfaceAndBottomTopography_FEM( UserContext *user )
{
	PetscMPIInt     rank;
	PetscErrorCode 	ierr;
	DM				cda, cda_SurfaceTopo, cda_BottomTopo;
	DMDACoor3d		***coors, ***coors_SurfaceTopo, ***coors_BottomTopo;
	Vec 			gc, gc_SurfaceTopo, gc_BottomTopo;
	PetscInt 		xs,ys,zs,xm,ym,zm, nnode_x,nnode_y, nnode_z;
	PetscInt		ix,iy,iz,iiz, zs_Z, zm_Z, xs_Z, xm_Z, ys_Z, ym_Z;
	Vec 			Surface_FIELD, Bottom_FIELD;
	PetscScalar 	***SURFACE_ZCOORD, ***BOTTOM_ZCOORD;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	/* Extract coordinates of the local portion of the velocity grid ----------*/
	ierr = DMGetCoordinateDM(user->DA_Vel,			&cda	); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_Vel,	&gc		); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,									&coors	); 	CHKERRQ(ierr);

	/* Extract coordinates of the local portion of the surface grid ----------*/
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,			&cda_SurfaceTopo	); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,		&gc_SurfaceTopo		); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,			&coors_SurfaceTopo	); 	CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(user->DA_BottomTopography,				&cda_BottomTopo		); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_BottomTopography,		&gc_BottomTopo		); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_BottomTopo,gc_BottomTopo,				&coors_BottomTopo	); 	CHKERRQ(ierr);


	ierr = DMDAGetInfo(user->DA_Vel, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);

	ierr = DMDAGetGhostCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(user->DA_Vel,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	// The parallel layout of the SurfaceTopography grid MUST be the same as that of the Vel grid (not updated crap)
	//PetscPrintf(PETSC_COMM_SELF," rank %lld, Corners: VEL: [%lld,%lld,%lld]-[%lld,%lld,%lld]  SURFACE: [%lld,%lld,%lld]-[%lld,%lld,%lld]\n",rank, xs,ys,zs,xs+xm-1,ys+ym-1,zs+zm-1, xs_Z,ys_Z,zs_Z,xs_Z+xm_Z-1,ys_Z+ym_Z-1,zs_Z+zm_Z);
	if ( (xm_Z != xm) || (ym_Z != ym) ){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"Parallel layout of SurfaceTopography and Vel arrays are not the same. Fix the bug!");
	}


	/* Extract z-coordinate */
	ierr = DMGetLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, Surface_FIELD, &SURFACE_ZCOORD );CHKERRQ(ierr);

	ierr = DMGetLocalVector( user->DA_BottomTopography, &Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_BottomTopography, Bottom_FIELD, &BOTTOM_ZCOORD );CHKERRQ(ierr);


	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				if (user->ErosionParameters.UseInternalFreeSurface==0){
					// surface topo
					if (iz== (nnode_z-1)){

						/* Update all x,y-coordinates of uppermost layer */
						for (iiz=zs_Z; iiz<zs_Z+zm_Z; iiz++){
							coors_SurfaceTopo[iiz][iy][ix].x = coors[nnode_z-1][iy][ix].x;
							coors_SurfaceTopo[iiz][iy][ix].y = coors[nnode_z-1][iy][ix].y;
							coors_SurfaceTopo[iiz][iy][ix].z = coors[nnode_z-1][iy][ix].z;
						}

						/* Update all z-coordinate values */
						for (iiz=zs_Z; iiz<zs_Z+zm_Z; iiz++){
							SURFACE_ZCOORD[iiz][iy][ix] = coors[nnode_z-1][iy][ix].z;
						}

					}
				}

				// bottom topo
				if (iz==0){

					/* Update all x,y-coordinates of uppermost layer */
					for (iiz=zs_Z; iiz<zs_Z+zm_Z; iiz++){
						coors_BottomTopo[iiz][iy][ix].x = coors[0][iy][ix].x;
						coors_BottomTopo[iiz][iy][ix].y = coors[0][iy][ix].y;
						coors_BottomTopo[iiz][iy][ix].z = coors[0][iy][ix].z;
					}

					/* Update all z-coordinate values */
					for (iiz=zs_Z; iiz<zs_Z+zm_Z; iiz++){
						BOTTOM_ZCOORD[iiz][iy][ix] = coors[0][iy][ix].z;
					}

				}

			}
		}
	}

	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_FIELD,&SURFACE_ZCOORD ); 						CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalBegin(user->DA_SurfaceTopography,Surface_FIELD,INSERT_VALUES,user->SurfaceTopography);	CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalEnd  (user->DA_SurfaceTopography,Surface_FIELD,INSERT_VALUES,user->SurfaceTopography);	CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );									CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(user->DA_BottomTopography, Bottom_FIELD,&BOTTOM_ZCOORD ); 						CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalBegin(user->DA_BottomTopography,Bottom_FIELD,INSERT_VALUES,user->BottomTopography);		CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalEnd  (user->DA_BottomTopography,Bottom_FIELD,INSERT_VALUES,user->BottomTopography);		CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_BottomTopography, &Bottom_FIELD );									CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(cda,gc,	&coors); 																CHKERRQ(ierr);
	//ierr 	= 	DMDestroy(cda); 																					CHKERRQ(ierr);
	//ierr 	= 	VecDestroy(gc); 																					CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
	ierr 	= 	DMDAVecRestoreArray(cda_BottomTopo,gc_BottomTopo,		&coors_BottomTopo); 							CHKERRQ(ierr);

	ierr 	= 	DASetCoordinatesFromLocalVector( user->DA_SurfaceTopography, gc_SurfaceTopo ); 						CHKERRQ(ierr);
	//ierr 	= 	VecDestroy(gc_SurfaceTopo); CHKERRQ(ierr);
	//ierr 	=	DMDestroy(cda_SurfaceTopo); CHKERRQ(ierr);

	ierr 	= 	DASetCoordinatesFromLocalVector( user->DA_BottomTopography, gc_BottomTopo ); 						CHKERRQ(ierr);
	//ierr 	= 	VecDestroy(gc_BottomTopo); CHKERRQ(ierr);
	//ierr 	=	DMDestroy(cda_BottomTopo); CHKERRQ(ierr);



	/* Scatter surface and bottom topography data to all 'internal' nodes 												*/
	if (user->ErosionParameters.UseInternalFreeSurface==0){
		ierr = ScatterTopographyData_to_DA_Surface(user->DA_SurfaceTopography, user->SurfaceTopography, 1); CHKERRQ(ierr);	// surface topo
	}
	ierr = ScatterTopographyData_to_DA_Surface(user->DA_BottomTopography,  user->BottomTopography,  0); CHKERRQ(ierr);	// bottom topo



#if 0
	///////// DEBUGGING: WRITE SOME INFO TO SCREEN

	ierr = DMGetLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, Surface_FIELD, &SURFACE_ZCOORD );CHKERRQ(ierr);

	ierr = DMGetLocalVector( user->DA_BottomTopography, &Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_BottomTopography, Bottom_FIELD, &BOTTOM_ZCOORD );CHKERRQ(ierr);


	ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);
	for (iz=zs_Z; iz<zs_Z+zm_Z; iz++){

		for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
			for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){

				PetscPrintf(PETSC_COMM_SELF,"rank %lld, [iz][iy][ix]=[%lld][%lld][%lld] SURFACE ZCOORD=%f; BOTTOM ZCOORD=%f \n",(LLD)rank,(LLD)iz,(LLD)iy,(LLD)ix,SURFACE_ZCOORD[iz][iy][ix], BOTTOM_ZCOORD[iz][iy][ix]);
			}
		}
	}


	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_FIELD,&SURFACE_ZCOORD ); 							CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );										CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(user->DA_BottomTopography, Bottom_FIELD,&BOTTOM_ZCOORD ); 							CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_BottomTopography, &Bottom_FIELD );										CHKERRQ(ierr);

#endif




	PetscFunctionReturn(0);
}
/*==========================================================================================================*/





/*==========================================================================================================*/
/* Ensures that info on the surface or bottom topography is correctly transmitted to every
 * PROC. This is particularly an issue if we have a large number of PROCs in the z-direction.
 * We do this by creating an index set, which finds the index of the points at the free surface and copies them
 *
 * This routine can do this for either the surface or the bottom data.
 *
 * SurfaceTopography=0 -> bottom topography
 * SurfaceTopography=1 -> surface topography
 *
 */
#undef __FUNCT__
#define __FUNCT__ "ScatterTopographyData_to_DA_Surface"
PetscErrorCode ScatterTopographyData_to_DA_Surface( DM da, Vec GlobalTopography, PetscInt SurfaceTopography)
{

	PetscErrorCode			ierr;
	Vec 					Topography_NaturalOrdering, Local_Topography;
	Vec 					gc_Topo, FIELD;
	DMDACoor3d				***coors_Topo;
	DM 					cda;
	PetscInt				*idx, num;
	PetscScalar			***LocalTopography, ***ZCOORD;
	IS						is;
	VecScatter 			newctx;
	PetscInt				ix,iy,iz,nx,ny,nz, xs_Z, ys_Z, zs_Z, xm_Z, ym_Z, zm_Z;

	ierr = DMDAGetCorners(da,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &nx,&ny, &nz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);


	// create global vector in natural ordering that holds the z-coordinates
	DMDACreateNaturalVector(da,&Topography_NaturalOrdering);
	DMDAGlobalToNaturalBegin(da,	GlobalTopography, INSERT_VALUES, Topography_NaturalOrdering);
	DMDAGlobalToNaturalEnd(da,	GlobalTopography, INSERT_VALUES, Topography_NaturalOrdering);

	VecCreate(PETSC_COMM_SELF,&Local_Topography);
	VecSetSizes(Local_Topography,xm_Z*ym_Z*zm_Z,xm_Z*ym_Z*zm_Z);
	VecSetType(Local_Topography, VECSEQ);

	// dynamically allocate PetscInt array of size xm_Z*ym_Z;
	idx = (PetscInt *)malloc(sizeof(PetscInt)*(size_t)(xm_Z*ym_Z*zm_Z));

	// convert i,j,k into a single number
	num=0;
	for (iz=zs_Z; iz<zs_Z+zm_Z; iz++){
		for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
			for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){
				if (SurfaceTopography==1){
					// we are scattering data of the surface topography
					idx[num] 	= 	( (nz-1)*(nx*ny) + iy*nx + ix) + 0;		// index of point @ free surface but with current [ix,iy] coordinates
					// equation for that is: index  = totdof*(k*(nx*ny) + j*nx + i) + dof;
				}
				else if (SurfaceTopography==0){
					// we are dealing with the bottom topography
					idx[num] 	= 	( (0)*(nx*ny) + iy*nx + ix) + 0;		// index of point @ free surface but with current [ix,iy] coordinates

				}
				else {
					// unknown
					SETERRQ1( PETSC_COMM_SELF, PETSC_ERR_SUP, "ScatterTopographyData_to_DA_Surface: Unknown option SurfaceTopography; should be either 0 or 1 and is %lld",(LLD)SurfaceTopography);
				}
				num					=	num+1;
			}
		}
	}



	// Create Index set
	ISCreateGeneral(PETSC_COMM_WORLD,xm_Z*ym_Z*zm_Z,idx,PETSC_USE_POINTER,&is);		// try PETSC_OWN_POINTER if it doesn't work; replaces ISCreateGeneralWithArray of earlier PETSC versions

	// Scatter data from global vector to local vector
	VecScatterCreate(Topography_NaturalOrdering,is, Local_Topography, 0, &newctx);
	VecScatterBegin(newctx,Topography_NaturalOrdering,Local_Topography, INSERT_VALUES,  SCATTER_FORWARD);
	VecScatterEnd  (newctx,Topography_NaturalOrdering,Local_Topography, INSERT_VALUES,  SCATTER_FORWARD);


	/* Send the data from the local portion to the global vector */
	ierr = DMGetCoordinateDM(da,				&cda	); 									CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,			&gc_Topo		); 							CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc_Topo,	&coors_Topo	); 										CHKERRQ(ierr);

	ierr = DMGetLocalVector( da, &FIELD );											CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, GlobalTopography, 	INSERT_VALUES, FIELD );		CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, GlobalTopography, 	INSERT_VALUES, FIELD );		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, FIELD, &ZCOORD );												CHKERRQ(ierr);

	ierr = DMDAGetCorners(da,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);

	ierr = VecGetArray3d(Local_Topography,zm_Z,ym_Z,xm_Z,0,0,0,&LocalTopography); CHKERRQ(ierr);
	for (iz=zs_Z; iz<zs_Z+zm_Z; iz++){
		for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
			for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){
				ZCOORD[iz][iy][ix] = LocalTopography[iz-zs_Z][iy-ys_Z][ix-xs_Z];
			}
		}
	}

	ierr = VecRestoreArray3d(Local_Topography,xm_Z,ym_Z,zm_Z,0,0,0,&LocalTopography); CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(cda, gc_Topo,	&coors_Topo); 					CHKERRQ(ierr);
	//ierr 	= 	VecDestroy(gc_Topo); 											CHKERRQ(ierr);
	//ierr 	= 	DMDestroy(cda); 												CHKERRQ(ierr);
	ierr 	= 	DMDAVecRestoreArray(da, FIELD,&ZCOORD ); 							CHKERRQ(ierr);

	ierr 	=	DMLocalToGlobalBegin(da,FIELD,INSERT_VALUES,GlobalTopography);		CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalEnd(da,FIELD,INSERT_VALUES,GlobalTopography);		CHKERRQ(ierr);

	ierr 	= 	DMRestoreLocalVector(da, &FIELD );								CHKERRQ(ierr);

	// Cleanup
	free(idx);
	ISDestroy(&is);
	VecDestroy(&Topography_NaturalOrdering);
	VecDestroy(&Local_Topography);
	VecScatterDestroy(&newctx);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Apply sedimentation to an internal free surface of the model.
 * Currently we only have the option to add a fixed sedimentation rate, and in this routine we simply advect
 * the internal free surface upwards with this rate.
 * In the future we can think about adding different sedimentation routines.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "ApplySedimentationToInternalFreeSurface"
PetscErrorCode ApplySedimentationToInternalFreeSurface( UserContext *user)
{
	PetscErrorCode  ierr;
	PetscInt		ix, iy, ys_Z, ym_Z, xs_Z, xm_Z, zs_Z, zm_Z;
	PetscScalar 	***LocalSurfaceTopography;
	DM				cda_SurfaceTopo;
	DMDACoor3d		***coors_SurfaceTopo;
	Vec				gc_SurfaceTopo;

	if (user->ErosionParameters.SedimentationModel>0){	// only perform this routine if we actually want to do sedimentation


		/* Developer note:
		 * 	Each processor has a local portion of the free surface.
		 * 	As we only advect it upwards with a constant sedimentation rate it is sufficient to update the local portion of the free surface.
		 * 	Once we decide to implement more complicated sedimentation routines this might have to be changed.
		 *
		 */

		/* Create global vector & update ghost points of the deformed surface grid */
		ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,				&gc_SurfaceTopo			); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);



		/* 1) Extract local x,y,z coordinates of surface DA
		 * 2) Advect them with the sedimentation velocity
		 */
		ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); CHKERRQ(ierr);	// surface topography
		for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
			for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){

				PetscScalar z, SedimentationVelocity, zmin, zmax;
				PetscInt 	numPhases, PhaseSedimented;

				/* Extract x,y,z coordinates of surface topography */
				//x 	= 	coors_SurfaceTopo[zs_Z][iy][ix].x;
				//y 	= 	coors_SurfaceTopo[zs_Z][iy][ix].y;
				z 		= 	LocalSurfaceTopography[zs_Z][iy][ix];	// topography
				zmin 	=	user->z_bot				;
				zmax 	=	user->z_bot	 	+ user->H;


				if (user->ErosionParameters.SedimentationModel==1){
					// constant sedimentation velocity throughout the model

					SedimentationVelocity = user->ErosionParameters.SedimentationRate;

					z 	= 	z	+ SedimentationVelocity*user->dt;		// update internal free surface coordinates according to sedimentation rate.

					// check if internal free surface goes outside the model domain
					if 		(z>zmax){	z = zmax;	}
					else if (z<zmin){	z = zmin;	}

					LocalSurfaceTopography[zs_Z][iy][ix] = z;		// set back to local surface topography


					// Compute which phase we should be sedimenting.
					numPhases 		= user->ErosionParameters.PhaseLastSedimentedLayer + 1 - user->ErosionParameters.PhaseFirstSedimentedLayer;
					PhaseSedimented = (PetscInt) (user->time/user->ErosionParameters.SedimentLayerThicknessYears);
					PhaseSedimented = PhaseSedimented- (PetscInt) (PhaseSedimented/numPhases)*numPhases + user->ErosionParameters.PhaseFirstSedimentedLayer;	// this ensures that we 'loop' the phases

					user->ErosionParameters.PhaseSedimented = PhaseSedimented;		//store the phase that is being sedimented
				}
				else{
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown sedimentation model");
				}


			}
		}




		/* Info */
		if (user->ErosionParameters.SedimentationModel==1){
			PetscPrintf(PETSC_COMM_WORLD,"Applying sedimentation to internal free surface of the model. Phase that is currently being sedimented is %i   \n", user->ErosionParameters.PhaseSedimented);
		}


		/* Cleaning up */
		ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
		ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
		ierr 	= 	DASetCoordinatesFromLocalVector( user->DA_SurfaceTopography, gc_SurfaceTopo ); 						CHKERRQ(ierr);


	}	// only do it if we wants sedimentation in our code


	PetscFunctionReturn(0);

}
/*==========================================================================================================*/




/*==========================================================================================================*/
/* Apply fast erosion to the internal free surface of the model
 *
 */
#undef __FUNCT__
#define __FUNCT__ "ApplyInfinitelyFastErosionToInternalFreeSurface"
PetscErrorCode ApplyInfinitelyFastErosionToInternalFreeSurface( UserContext *user)
{
	PetscErrorCode  ierr;
	PetscInt		size;
	PetscScalar 	AverageHeight;

	if (user->ErosionParameters.ErosionModel==1){	// only perform this routine if we actually want to do infinitely fast erosion

		/* Compute average height of internal free surface */
		ierr 			=	VecNorm(user->SurfaceTopography,NORM_1, &AverageHeight);	CHKERRQ(ierr);
		ierr			=	VecGetSize(user->SurfaceTopography,&size);					CHKERRQ(ierr);
		AverageHeight 	= 	AverageHeight/((PetscScalar) size);

		ierr 			= 	VecSet(user->SurfaceTopography,AverageHeight);				CHKERRQ(ierr);		// set all values to the average height

		PetscPrintf(PETSC_COMM_WORLD,"Applying infinitely fast erosion to internal free surface of the model. Average free surface height = %f m \n", AverageHeight*user->Characteristic.Length);

	}	// only do it if we want infinitely fast erosion in our code


	PetscFunctionReturn(0);

}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Apply a serial erosion code to the internal free surface
 *
 */
#undef __FUNCT__
#define __FUNCT__ "ApplySerialErosionCodeToInternalFreeSurface"
PetscErrorCode ApplySerialErosionCodeToInternalFreeSurface(UserContext *user )
{

	PetscErrorCode  ierr;
	PetscMPIInt     rank;
	PetscInt		M,N,P;
	Vec				InternalFreeSurfaceTopography, FreeSurface_Vx, FreeSurface_Vy, FreeSurface_Vz;
	DM				DMDA_InternalFreeSurfaceOnRankZero;
    PetscLogDouble  cputime_start0, cputime_end0;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscPrintf(PETSC_COMM_WORLD,"Applying a serial erosion code to the internal free surface *********** \n");

	/* Create a DA on rank 0 that has the same size as the global internal free surface but is 2D (3rd dimension not required, as they are copies of each other) */
	ierr = DMDAGetInfo(user->DA_SurfaceTopography,PETSC_NULL,&M,&N,&P,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
	if (rank==0){
		ierr = DMDACreate2d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,M, N,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL, &DMDA_InternalFreeSurfaceOnRankZero);
		ierr = DMDASetUniformCoordinates(DMDA_InternalFreeSurfaceOnRankZero,user->x_left,user->x_left + user->W,user->y_front,user->y_front+user->L,user->z_bot,user->z_bot+user->H);	// set coordinates; will be overwritten in CopyInternalFreeSurfaceToRankZero
		ierr = DMCreateLocalVector(DMDA_InternalFreeSurfaceOnRankZero, &InternalFreeSurfaceTopography);	// 2D vector that holds topography
		ierr = DMCreateLocalVector(DMDA_InternalFreeSurfaceOnRankZero, &FreeSurface_Vx);	// 2D vector that holds topography
		ierr = DMCreateLocalVector(DMDA_InternalFreeSurfaceOnRankZero, &FreeSurface_Vy);	// 2D vector that holds topography
		ierr = DMCreateLocalVector(DMDA_InternalFreeSurfaceOnRankZero, &FreeSurface_Vz);	// 2D vector that holds topography

	}

    PetscTime(&cputime_start0);

	/* Step 1: collect the internal free surface on processor 0 on a 2D DMDA */
	ierr = CopyInternalFreeSurfaceToRankZero(user, DMDA_InternalFreeSurfaceOnRankZero, InternalFreeSurfaceTopography, FreeSurface_Vx, FreeSurface_Vy, FreeSurface_Vz );	CHKERRQ(ierr);
    PetscTime(&cputime_end0);
    PetscPrintf(PETSC_COMM_WORLD,"  FE Erosion code: CopyInternalFreeSurfaceToRankZero took %f s \n",cputime_end0-cputime_start0);


	if (rank==0){

		/* Step 2: 		Call an erosion code on processor 0*/

		/* Step 2.1:  	Advect the (high-resolution) 'erosion' free surface with the tectonic velocities (given on a coarser resolution) */
	    PetscTime(&cputime_start0);
		ierr = AdvectInternalErosionSurfaceOnRankZeroWithTectonicVelocitiesAndReinterpolateToRegularMesh(user, DMDA_InternalFreeSurfaceOnRankZero, FreeSurface_Vx, FreeSurface_Vy, FreeSurface_Vz); CHKERRQ(ierr);
  
        
	    PetscTime(&cputime_end0);
	    PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: AdvectInternalErosionSurfaceOnRankZeroWithTectonicVelocitiesAndReinterpolateToRegularMesh took %f s \n",cputime_end0-cputime_start0);

		/* Step 2.2:  	Call the external erosion code and perform erosion sub-timesteps until we reach dt */
//		ierr = FE_ErosionCode_TectonicTimestep( user, DMDA_InternalFreeSurfaceOnRankZero, InternalFreeSurfaceTopography); 					CHKERRQ(ierr);
	    PetscTime(&cputime_start0);
	    ierr = FE_ErosionCode_TectonicTimestep( user); 					CHKERRQ(ierr);
	    PetscTime(&cputime_end0);
	    PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: FE_ErosionCode_TectonicTimestep took %f s \n",cputime_end0-cputime_start0);

		/* Step 2.3:	Interpolate the high-resolution 'erosion' free surface to the lower-resolution internal free surface all on rank 0*/
	    PetscTime(&cputime_start0);
	   	ierr = InterpolateErosionSurfaceToInternalFreeSurface(user, DMDA_InternalFreeSurfaceOnRankZero, InternalFreeSurfaceTopography);		CHKERRQ(ierr);
	    PetscTime(&cputime_end0);
	    PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: InterpolateErosionSurfaceToInternalFreeSurface took %f s \n",cputime_end0-cputime_start0);

	}

	/* Step 3: interpolate results of the erosion code back onto the (coarser) internal free surface & distribute to all processors	*/
	PetscTime(&cputime_start0);
	ierr =  CopyErodedSurfaceFromRankZeroToAllOtherRanks(user, DMDA_InternalFreeSurfaceOnRankZero, InternalFreeSurfaceTopography ); CHKERRQ(ierr);
    PetscTime(&cputime_end0);
    PetscPrintf(PETSC_COMM_WORLD,"  FE Erosion code: CopyErodedSurfaceFromRankZeroToAllOtherRanks took %f s \n",cputime_end0-cputime_start0);

	/* Cleaning up */
	if (rank==0){
		DMDestroy(&DMDA_InternalFreeSurfaceOnRankZero);
		VecDestroy(&InternalFreeSurfaceTopography);
		VecDestroy(&FreeSurface_Vx);
		VecDestroy(&FreeSurface_Vy);
		VecDestroy(&FreeSurface_Vz);

	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Interpolate the high-resolution erosion free surface back to the lower-resolution internal free surface
 * that is employed for the mechanics.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "InterpolateErosionSurfaceToInternalFreeSurface"
PetscErrorCode InterpolateErosionSurfaceToInternalFreeSurface(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography)
{
	PetscErrorCode 	ierr;
	PetscInt 		ResolutionFactorX,ResolutionFactorY, nx_fine, ny_fine, nx_coarse, ny_coarse, ix, iy,i,j;
	PetscScalar		**Topography_coarse, **Topography_fine;

	ResolutionFactorX = user->ErosionParameters.FE_ErosionCode.ResolutionFactorX;
	ResolutionFactorY = user->ErosionParameters.FE_ErosionCode.ResolutionFactorY;


	/* Interpolate tectonic velocity from a coarse mesh (based on the Internal Free Surface), to a fine mesh (on which ErosionTopography is defined) */
	ierr = DMDAGetInfo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,PETSC_NULL,&nx_fine,&ny_fine,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
	ierr = DMDAGetInfo(DMDA_InternalFreeSurfaceOnRankZero,						PETSC_NULL,&nx_coarse,&ny_coarse,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);


	/* Extract coarse and fine resolution topography */

	VecSet(InternalFreeSurfaceTopography,0.0);

	ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		InternalFreeSurfaceTopography, 			&Topography_coarse); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Topography_fine   ); 	CHKERRQ(ierr);

	/* Loop over coarse grid coordinates */
	for (iy=0; iy< ny_coarse; iy++){
		for(ix=0; ix< nx_coarse; ix++){
			PetscInt 		TopographyMethod, iix, iiy;
			PetscScalar		Topography_Value;

			TopographyMethod = 2;	// 0 - injection; 1-averaging


			if 		(TopographyMethod==0){
				// Inject value
				iix 				= 	ResolutionFactorX*ix;
				iiy 				= 	ResolutionFactorY*iy;

				Topography_Value 	=	Topography_fine[iiy][iix];

			}
			else if (TopographyMethod==1){
				PetscScalar num;

				// average as long as we're not at the boundary
				if ((ix>0) && (iy>0) && (ix<(nx_coarse-1)) && (iy<(ny_coarse-1))){

					// Average topography
					Topography_Value   	=   0.0;
					num 				=	0.0;
					for (i=0; i<ResolutionFactorX; i++){
						for (j=0; j<ResolutionFactorY; j++){
							iix 				= 	ResolutionFactorX*ix + i;
							iiy 				= 	ResolutionFactorY*iy + j;

							Topography_Value 	=	Topography_Value + Topography_fine[iiy][iix];

							num = num + 1.0;
						}
					}
					Topography_Value = Topography_Value/num;	// average topography

				}
				else{
					iix 				= 	ResolutionFactorX*ix;
					iiy 				= 	ResolutionFactorY*iy;

					Topography_Value 	=	Topography_fine[iiy][iix];
				}

			}
                        
            else if (TopographyMethod==2){
				PetscScalar num;
                Topography_Value   	=   0.0;
                num 				=	0.0;
                for (i=-ResolutionFactorX; i<ResolutionFactorX+1; i++){
                    for (j=-ResolutionFactorY; j<ResolutionFactorY+1; j++){
                        
                        iix 				= 	ResolutionFactorX*ix + i;
                        iiy 				= 	ResolutionFactorY*iy + j;
                        
                        if ( iix<0 || iix>((nx_coarse-1)*ResolutionFactorX) ) {continue;}
                        if ( iiy<0 || iiy>((ny_coarse-1)*ResolutionFactorY) ) {continue;}
                        
                        Topography_Value 	=	Topography_Value + Topography_fine[iiy][iix];
                        
                        num = num + 1.0;
                    }
                }
                
                Topography_Value = Topography_Value/num;	// average topography
            }

                        
        
			// Set value in coarse grid
			Topography_coarse[iy][ix] = Topography_Value;

		}
	}



	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Topography_fine   ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		InternalFreeSurfaceTopography, 			&Topography_coarse); 	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Advects the internal erosion surface with tectonic velocities (typically given at a coarser mesh).
 * The advected surface will typically not be a regular mesh anymore, so we interpolate it to a regular mesh.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "AdvectInternalErosionSurfaceOnRankZeroWithTectonicVelocitiesAndReinterpolateToRegularMesh"
PetscErrorCode AdvectInternalErosionSurfaceOnRankZeroWithTectonicVelocitiesAndReinterpolateToRegularMesh(UserContext *user,
		DM DMDA_InternalFreeSurfaceOnRankZero, Vec FreeSurface_Vx, Vec FreeSurface_Vy, Vec FreeSurface_Vz )
{
	PetscErrorCode  ierr;
	PetscInt		xs,ys,xm,ym,ix,iy, ResolutionFactorX,ResolutionFactorY, iix, iiy, nx_coarse, ny_coarse, nx_fine, ny_fine, i, j;
	Vec 			Topography_Vx, Topography_Vy, Topography_Vz, Topography_Advected_x, Topography_Advected_y, Topography_Advected_z, gc;
	PetscScalar		**Vx_fine,   **Vy_fine,   **Vz_fine, **Topo_x, **Topo_y, **Topo_z, **Erosion_Topo;
	PetscScalar		**Vx_coarse, **Vy_coarse, **Vz_coarse;
	DMDACoor2d       **coors;
	DM 				cda;
	PetscScalar	 	xmin_advected, xmax_advected, ymin_advected, ymax_advected;


	ResolutionFactorX = user->ErosionParameters.FE_ErosionCode.ResolutionFactorX;
	ResolutionFactorY = user->ErosionParameters.FE_ErosionCode.ResolutionFactorY;


	/* Interpolate tectonic velocity from a coarse mesh (based on the Internal Free Surface), to a fine mesh (on which ErosionTopography is defined) */
	ierr = DMDAGetCorners(DMDA_InternalFreeSurfaceOnRankZero,&xs,&ys, PETSC_NULL, &xm,&ym, PETSC_NULL); CHKERRQ(ierr); // get size of coarse grid (that stores tectonic velocities)

	ierr = DMDAGetInfo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,PETSC_NULL,&nx_fine,&ny_fine,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
	ierr = DMDAGetInfo(DMDA_InternalFreeSurfaceOnRankZero,						PETSC_NULL,&nx_coarse,&ny_coarse,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);


	/* Computed the boundaries of the FDSTAG grid AFTER advection by the background strain rate */
	xmin_advected 	=	user->x_left 			;
	xmax_advected 	=	user->x_left 	+ user->W;
	ymin_advected 	=	user->y_front			;
	ymax_advected 	=	user->y_front 	+ user->L;
	xmin_advected	=	xmin_advected - user->dt*xmin_advected*user->BC.Exx;
	xmax_advected	=	xmax_advected - user->dt*xmax_advected*user->BC.Exx;
	ymin_advected	=	ymin_advected - user->dt*ymin_advected*user->BC.Eyy;
	ymax_advected	=	ymax_advected - user->dt*ymax_advected*user->BC.Eyy;


	/* Extract coarse resolution velocities */
	ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vx, 			&Vx_coarse ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vy, 			&Vy_coarse ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vz, 			&Vz_coarse ); 	CHKERRQ(ierr);

	/* Extract fine resolution velocities */
	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Vx);	// fine grid Vx velocity
	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Vy);	// fine grid Vy velocity
	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Vz);	// fine grid Vz velocity

	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vx, 			&Vx_fine   ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vy, 			&Vy_fine   ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vz, 			&Vz_fine   ); 	CHKERRQ(ierr);

	/* Set coarse grid values in fine grid [every ResolutionFactor cells] */
	for (iy=0; iy<ny_coarse-1; iy++){
		for(ix=0; ix<nx_coarse-1; ix++){
			PetscScalar		p_x[5],  p_y[5],  Vx[4],  Vy[4],  Vz[4];

			/* Interpolate velocity on fine grid from the one on the coarse grid */
			p_x[0] = 0;			p_y[0] = 0;
			p_x[1] = 0;			p_y[1] = 1;
			p_x[2] = 1;			p_y[2] = 1;
			p_x[3] = 1;			p_y[3] = 0;
			p_x[4] = p_x[0]; 	p_y[4] = p_y[0];

			Vx[0] = Vx_coarse[iy  ][ix  ];	Vx[1] = Vx_coarse[iy+1][ix  ];	Vx[2] = Vx_coarse[iy+1][ix+1];	Vx[3] = Vx_coarse[iy  ][ix+1];
			Vy[0] = Vy_coarse[iy  ][ix  ];	Vy[1] = Vy_coarse[iy+1][ix  ];	Vy[2] = Vy_coarse[iy+1][ix+1];	Vy[3] = Vy_coarse[iy  ][ix+1];
			Vz[0] = Vz_coarse[iy  ][ix  ];	Vz[1] = Vz_coarse[iy+1][ix  ];	Vz[2] = Vz_coarse[iy+1][ix+1];	Vz[3] = Vz_coarse[iy  ][ix+1];

			for (i=0; i<=ResolutionFactorX; i++){
				for (j=0; j<=ResolutionFactorY; j++){
					PetscScalar fac_x, fac_y, Vz_val, Vx_val, Vy_val;


					fac_x = ((PetscScalar) i)/((PetscScalar) ResolutionFactorX );
					fac_y = ((PetscScalar) j)/((PetscScalar) ResolutionFactorY );

					// interpolate as if it were a linear shape function
					Vx_val 	= 	InterpolateWithin2DLinearElement(p_x, p_y, Vx, fac_x, fac_y);		// interpolate surface coordinates, using linear shape function
					Vy_val 	=  	InterpolateWithin2DLinearElement(p_x, p_y, Vy, fac_x, fac_y);		// interpolate surface coordinates, using linear shape function
					Vz_val 	=  	InterpolateWithin2DLinearElement(p_x, p_y, Vz, fac_x, fac_y);		// interpolate surface coordinates, using linear shape function


					iix 	= 	ResolutionFactorX*ix + i;
					iiy 	= 	ResolutionFactorY*iy + j;


					Vx_fine[iiy][iix] = Vx_val;
					Vy_fine[iiy][iix] = Vy_val;
					Vz_fine[iiy][iix] = Vz_val;

				}
			}



		}
	}


	ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vx, 			&Vx_coarse ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vy, 			&Vy_coarse ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vz, 			&Vz_coarse ); 	CHKERRQ(ierr);


	/* Perform the actual advection of the erosion surface  --- */
	DMGetCoordinateDM(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&cda);
	DMGetCoordinatesLocal(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&gc);
	DMDAVecGetArray(cda,gc,&coors);

	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Advected_x);	// advected x coords
	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Advected_y);	// advected y coords
	ierr = DMCreateLocalVector(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, &Topography_Advected_z);	// advected z coords

	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_x, 			&Topo_x   ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_y, 			&Topo_y   ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_z, 			&Topo_z   ); 	CHKERRQ(ierr);

	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Erosion_Topo   ); 	CHKERRQ(ierr);


	for (iy=0; iy<ny_fine; iy++){
		for(ix=0; ix<nx_fine; ix++){
			PetscScalar  eps=1e-6*(xmax_advected-xmin_advected);

			/* Advect topography */
			Topo_x[iy][ix] = coors[iy][ix].x 		+ Vx_fine[iy][ix]* user->dt;
			Topo_y[iy][ix] = coors[iy][ix].y 		+ Vy_fine[iy][ix]* user->dt;
			Topo_z[iy][ix] = Erosion_Topo[iy][ix] 	+ Vz_fine[iy][ix]* user->dt;

			/* Extend grid slightly at the boundaries */
			if (iy==0){
				//Topo_y[iy][ix] = Topo_y[iy][ix] - eps;
				Topo_y[iy][ix] = ymin_advected - eps;

			}
			if (iy==ny_fine-1){
				//Topo_y[iy][ix] = Topo_y[iy][ix] + eps;
				Topo_y[iy][ix] = ymax_advected + eps;
			}
			if (ix==0){
				//Topo_x[iy][ix] = Topo_x[iy][ix] - eps;
				Topo_x[iy][ix] = xmin_advected - eps;

			}
			if (ix==nx_fine-1){
				//Topo_x[iy][ix] = Topo_x[iy][ix] + eps;
				Topo_x[iy][ix] = xmax_advected + eps;
			}

		}
	}


	DMDAVecRestoreArray(cda,gc,&coors);


	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vx, 			&Vx_fine   ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vy, 			&Vy_fine   ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Vz, 			&Vz_fine   ); 	CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Erosion_Topo   ); 	CHKERRQ(ierr);
	/* End of advection --------------------------------------- */


	/* Interpolate it back on a regular grid ------------------ */

	// Remesh the erosion surface to a regular grid taking the advection of the grid by the BG strainrate into account
	ierr = DMDASetUniformCoordinates(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,xmin_advected,xmax_advected,ymin_advected,ymax_advected,0,0);

	DMGetCoordinateDM(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&cda);
	DMGetCoordinatesLocal(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&gc);
	DMDAVecGetArray(cda,gc,&coors);

	VecSet(user->ErosionParameters.FE_ErosionCode.ErosionSurface,0.0);
	ierr = DMDAVecGetArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Erosion_Topo   ); 	CHKERRQ(ierr);

	/* Interpolate back to regular grid */
	for (iy=0; iy<ny_fine-1; iy++){
		for(ix=0; ix<nx_fine-1; ix++){

			PetscScalar		p_x[5],  p_y[5],  p_z[5];
			PetscInt 		iix_min, iix_max, iiy_min, iiy_max; //FoundPolygon;

			/* Extract a contour of the linear element and use an isinside routine to find out if the point is inside the element */
			p_x[0] = Topo_x[iy  ][ix  ];	p_y[0] = Topo_y[iy  ][ix  ]; p_z[0] = Topo_z[iy  ][ix  ];
			p_x[1] = Topo_x[iy  ][ix+1];	p_y[1] = Topo_y[iy  ][ix+1]; p_z[1] = Topo_z[iy  ][ix+1];
			p_x[2] = Topo_x[iy+1][ix+1];	p_y[2] = Topo_y[iy+1][ix+1]; p_z[2] = Topo_z[iy+1][ix+1];
			p_x[3] = Topo_x[iy+1][ix  ];	p_y[3] = Topo_y[iy+1][ix  ]; p_z[3] = Topo_z[iy+1][ix  ];
			p_x[4] = p_x[0]; 				p_y[4] = p_y[0]; 			 p_z[4] = p_z[0];


			/* Loop over regular surface */
			iix_min = ix-20;	iix_min = PetscMax(iix_min, 0		);		// assumes that points can not move too much in one timestep
			iix_max = ix+20;	iix_max = PetscMin(iix_max, nx_fine	);
			iiy_min = iy-20;	iiy_min = PetscMax(iiy_min, 0		);
			iiy_max = iy+20;	iiy_max = PetscMin(iiy_max, ny_fine	);

//			iix_min =0; iix_max=nx_fine;
//			iiy_min =0; iiy_max=ny_fine;

//			FoundPolygon=0;
			for (iiy=iiy_min; iiy<iiy_max; iiy++){
					for(iix=iix_min; iix<iix_max; iix++){

					PetscScalar 	x,y,z_interpolated;
					PetscBool		inside;

					/* Regular point */
					x = coors[iiy][iix].x;
					y = coors[iiy][iix].y;

					/* Is point inside square polygon? */
					inside = pnpoly(4, p_x, p_y, x,y);

					/*
					// debugging
					if (iiy==ny_fine-1 && iix==10){
						PetscInt T;
						if (inside){
							T=1;
						}
						else{
							T=0;
						}

						PetscPrintf(PETSC_COMM_WORLD," x,y=[%f,%f] p_x=[%f,%f,%f,%f] p_y=[%f,%f,%f,%f] ymin/max_advected = %f %f, T=%i\n",x,y,p_x[0],p_x[1],p_x[2],p_x[3],p_y[0],p_y[1],p_y[2],p_y[3],ymin_advected, ymax_advected, T);
					}
					*/

					if (inside){
						z_interpolated 			=  	InterpolateWithin2DLinearElement(p_x, p_y, p_z, x, y);		// interpolate surface coordinates, using linear shape function
						Erosion_Topo[iiy][iix] 	= 	z_interpolated;
//						FoundPolygon 			=	1;
					}


				}
}
		//	if (FoundPolygon==0){
		//		PetscPrintf(PETSC_COMM_WORLD,"Ooops, did not find surrounding polygon for ix=%i, iy=%i \n",ix,iy);
		//	}


		}
	}
	DMDAVecRestoreArray(cda,gc,&coors);


	/* End of interpolation to regular grid ------------------- */

	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		user->ErosionParameters.FE_ErosionCode.ErosionSurface, 		&Erosion_Topo   ); 	CHKERRQ(ierr);


	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_x, 			&Topo_x   ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_y, 			&Topo_y   ); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, 		Topography_Advected_z, 			&Topo_z   ); 	CHKERRQ(ierr);


	/* Cleanup */
	ierr = VecDestroy(&Topography_Vx); CHKERRQ(ierr);
	ierr = VecDestroy(&Topography_Vy); CHKERRQ(ierr);
	ierr = VecDestroy(&Topography_Vz); CHKERRQ(ierr);

	ierr = VecDestroy(&Topography_Advected_x); CHKERRQ(ierr);
	ierr = VecDestroy(&Topography_Advected_y); CHKERRQ(ierr);
	ierr = VecDestroy(&Topography_Advected_z); CHKERRQ(ierr);



PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================
 * As the function name suggests, this routine copies the eroded internal free surface from rank zero (where it is on a 2D DA) back
 * to all other ranks on a 3D array.
 *
 */

PetscErrorCode CopyErodedSurfaceFromRankZeroToAllOtherRanks(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography )
{

	PetscErrorCode  ierr;
	PetscMPIInt     rank;
	PetscInt		M,N,P;
	DM				DMDA_InternalFreeSurface_3D;
	Vec 			FreeSurfaceAtRankZero, natural_ordering;
	VecScatter		ctx;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	/* Info on the 3D DMDA that contains the internal free surface */
	ierr = DMDAGetInfo(user->DA_SurfaceTopography,PETSC_NULL,&M,&N,&P,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);

	/* Create a 3D DMDA on rank zero to which we will copy the values */
	if (rank==0){
		ierr = DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,M, N,P, 1, 1,1,1,1,PETSC_NULL,PETSC_NULL, PETSC_NULL, &DMDA_InternalFreeSurface_3D);
	}

	/* Create scatter context */
	VecScatterCreateToZero(user->SurfaceTopography,&ctx,&FreeSurfaceAtRankZero);						// create a SEQ vector on rank 0 from which we'll scatter the values. On all other procs it's zero

	if (rank==0){
		/* Update vector from 2D DA to 3D array on processor 0 */

		PetscInt 		i,j,k;
		PetscScalar 	***FreeSurface3D, **FreeSurface2D;

		// arrays that contain free surface in the 3D and 2D arrays
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceAtRankZero, 					&FreeSurface3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		InternalFreeSurfaceTopography, 			&FreeSurface2D ); 	CHKERRQ(ierr);

		for (i=0; i<M; i++){
			for (j=0; j<N; j++){

				for (k=0; k<P; k++){

					// Update topography in 2D DA
					FreeSurface3D[k][j][i] = FreeSurface2D[j][i];
				}

			}
		}

		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 	InternalFreeSurfaceTopography, 	&FreeSurface2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurface_3D, 		FreeSurfaceAtRankZero, 			&FreeSurface3D ); 	CHKERRQ(ierr);


	}

	/* Scatter the vector from proc 0 to all other procs */
	DMDACreateNaturalVector(user->DA_SurfaceTopography, &natural_ordering);		// natural ordering vector


	VecScatterBegin(ctx,FreeSurfaceAtRankZero,natural_ordering,INSERT_VALUES,SCATTER_REVERSE);
	VecScatterEnd(ctx,FreeSurfaceAtRankZero,natural_ordering,INSERT_VALUES,SCATTER_REVERSE);

	/* Transfer from natural to DMDA ordering */
	DMDANaturalToGlobalBegin(user->DA_SurfaceTopography, natural_ordering, INSERT_VALUES, user->SurfaceTopography);
	DMDANaturalToGlobalEnd(user->DA_SurfaceTopography,   natural_ordering, INSERT_VALUES, user->SurfaceTopography);
	VecDestroy(&natural_ordering);



	/* Cleaning up */
	if (rank==0){
		ierr = DMDestroy(&DMDA_InternalFreeSurface_3D);	CHKERRQ(ierr);
	}
	VecScatterDestroy(&ctx);
	VecDestroy(&FreeSurfaceAtRankZero);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* Copy internal free surface to a 2D DMDA on rank 0, including the correct coordinates.
 * This simplifies calling a serial erosion code.
 *
 * NOTE: this algorithm is inherently NOT scalable. The main issue is the lack of a MPI-based scalable erosion code.
 *
 */
PetscErrorCode CopyInternalFreeSurfaceToRankZero(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography,
		Vec FreeSurface_Vx, Vec FreeSurface_Vy, Vec FreeSurface_Vz )
{
	PetscMPIInt     rank;
	PetscErrorCode	ierr;
	PetscInt		M,N,P,i,j;
	PetscScalar		***FreeSurface3D, **FreeSurface2D, ***FreeSurfaceVx_3D, **FreeSurfaceVx_2D, ***FreeSurfaceVy_3D, **FreeSurfaceVy_2D, ***FreeSurfaceVz_3D, **FreeSurfaceVz_2D;
	Vec 			FreeSurfaceAtRankZero, InternalFreeSurface_3D_coords, InternalFreeSurface_3D_coords_zero, InternalFreeSurface_2D_coords;
	Vec 			FreeSurfaceVxAtRankZero, FreeSurfaceVyAtRankZero, FreeSurfaceVzAtRankZero, natural_ordering, temp;
	VecScatter		ctx, ctx_coords;
	DM				DMDA_InternalFreeSurface_3D, coordDA_3D, coordDA_2D, cda;
	DMDACoor3d		***coors_3D;
	DMDACoor2d		**coors_2D;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	/* Info on the 3D DMDA that contains the internal free surface */
	ierr = DMDAGetInfo(user->DA_SurfaceTopography,PETSC_NULL,&M,&N,&P,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);


	/* Create a 3D DMDA on rank zero to which we will copy the values */
	if (rank==0){
		ierr = DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,M, N, P, 1, 1,1,1,1,PETSC_NULL,PETSC_NULL, PETSC_NULL, &DMDA_InternalFreeSurface_3D);	CHKERRQ(ierr);
		ierr = DMDASetUniformCoordinates(DMDA_InternalFreeSurface_3D,0,1,0,1,0,1);		CHKERRQ(ierr); // initialize coordinates [will be reset later in this routine!]
	}
	ierr = DMGetCoordinates(user->DA_SurfaceTopography,&InternalFreeSurface_3D_coords);	// 3D coordinate vector


	/* We have to scatter the global values from all processors to rank 0
	 * For this to work, however, we have to first make sure that the values are in natural ordering, hence the GlobalToNatural commands
	 *
	 */

	/* Scatter global distributed values to local vector on proc 0 */
	DMDACreateNaturalVector(user->DA_SurfaceTopography, &natural_ordering);
	DMDAGlobalToNaturalBegin(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, natural_ordering);
	DMDAGlobalToNaturalEnd(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, natural_ordering);

	VecScatterCreateToZero(natural_ordering,&ctx,&FreeSurfaceAtRankZero);						// create a SEQ vector on rank 0 to which we'll scatter the values. On all other procs it's zero
	VecScatterBegin(ctx,natural_ordering,FreeSurfaceAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,natural_ordering,FreeSurfaceAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterDestroy(&ctx);
	VecDestroy(&natural_ordering);


	// Scatter surface velocity to rank 0
	DMDACreateNaturalVector(user->DA_SurfaceTopography, &natural_ordering);
	DMDAGlobalToNaturalBegin(user->DA_SurfaceTopography, 	user->SurfaceTopography_Vx, INSERT_VALUES, natural_ordering);
	DMDAGlobalToNaturalEnd(user->DA_SurfaceTopography, 		user->SurfaceTopography_Vx, INSERT_VALUES, natural_ordering);
	VecScatterCreateToZero(natural_ordering,&ctx,&FreeSurfaceVxAtRankZero);
	VecScatterBegin(ctx,natural_ordering,		FreeSurfaceVxAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,	natural_ordering,		FreeSurfaceVxAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterDestroy(&ctx);
	VecDestroy(&natural_ordering);

	DMDACreateNaturalVector(user->DA_SurfaceTopography, &natural_ordering);
	DMDAGlobalToNaturalBegin(user->DA_SurfaceTopography, user->SurfaceTopography_Vy, INSERT_VALUES, natural_ordering);
	DMDAGlobalToNaturalEnd(user->DA_SurfaceTopography, user->SurfaceTopography_Vy, INSERT_VALUES, natural_ordering);
	VecScatterCreateToZero(natural_ordering,&ctx,&FreeSurfaceVyAtRankZero);
	VecScatterBegin(ctx,natural_ordering,FreeSurfaceVyAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,natural_ordering,FreeSurfaceVyAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterDestroy(&ctx);
	VecDestroy(&natural_ordering);


	DMDACreateNaturalVector(user->DA_SurfaceTopography, &natural_ordering);
	DMDAGlobalToNaturalBegin(user->DA_SurfaceTopography, user->SurfaceTopography_Vz, INSERT_VALUES, natural_ordering);
	DMDAGlobalToNaturalEnd(user->DA_SurfaceTopography, user->SurfaceTopography_Vz, INSERT_VALUES, natural_ordering);
	VecScatterCreateToZero(natural_ordering,&ctx,&FreeSurfaceVzAtRankZero);
	VecScatterBegin(ctx,natural_ordering,FreeSurfaceVzAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,natural_ordering,FreeSurfaceVzAtRankZero,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterDestroy(&ctx);
	VecDestroy(&natural_ordering);


	/* Scatter coordinates to processor 0 */

	// change ordering into natural ordering
	DMGetCoordinateDM(user->DA_SurfaceTopography, &cda);
	DMDACreateNaturalVector(cda, &natural_ordering);
	DMDAGlobalToNaturalBegin(cda, InternalFreeSurface_3D_coords, INSERT_VALUES, natural_ordering);
	DMDAGlobalToNaturalEnd(cda, InternalFreeSurface_3D_coords, INSERT_VALUES, natural_ordering);

	// copy into temporary local vector
	VecScatterCreateToZero(natural_ordering,&ctx_coords,&temp);						// create a SEQ vector on rank 0 to which we'll scatter the values. On all other procs it's zero
	VecScatterBegin(ctx_coords,natural_ordering,temp,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx_coords,natural_ordering,temp,INSERT_VALUES,SCATTER_FORWARD);
	VecDestroy(&natural_ordering);



	if (rank==0){

		/* Copy values from 3D DMDA to 2D DA that contains the internal free surface and that will be passed over to the 2D erosion code
		 * The 3D DMDA contains the same values, but in the z-direction it has copies of the internal free surface, which is used in various parts
		 * of LaMEM, for example to track whether a tracer is above or below the free surface.
		 *
		 * The erosion code, on the other hand, should not have to be aware of these LaMEM-internal 'tricks',
		 * and should only receive a 2D DMDA that contains the topography at the LaMEM x and y resolution, and the change in topography during the last timestep.
		 * That is done here by copying values
		 *
		 */

		// copy local vector over to
		ierr = 	DMGetCoordinates(DMDA_InternalFreeSurface_3D,&InternalFreeSurface_3D_coords_zero);		CHKERRQ(ierr);
		ierr =	VecCopy(temp, InternalFreeSurface_3D_coords_zero);										CHKERRQ(ierr);	// copy temporary vector to vector with correct blocksize (3)

		// arrays that contain free surface in the 3D and 2D arrays
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceAtRankZero, 					&FreeSurface3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		InternalFreeSurfaceTopography, 			&FreeSurface2D ); 	CHKERRQ(ierr);

		ierr = DMDAVecGetArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVxAtRankZero, 				&FreeSurfaceVx_3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vx, 						&FreeSurfaceVx_2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVyAtRankZero, 				&FreeSurfaceVy_3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vy, 						&FreeSurfaceVy_2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVzAtRankZero, 				&FreeSurfaceVz_3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vz, 						&FreeSurfaceVz_2D ); 	CHKERRQ(ierr);



		ierr = DMGetCoordinateDM(DMDA_InternalFreeSurface_3D, 								&coordDA_3D);						CHKERRQ(ierr);

		ierr = DMDAVecGetArray(coordDA_3D, 				InternalFreeSurface_3D_coords_zero, 	&coors_3D ); 						CHKERRQ(ierr);

		ierr = DMGetCoordinateDM(DMDA_InternalFreeSurfaceOnRankZero, 							&coordDA_2D);						CHKERRQ(ierr);
		ierr = DMGetCoordinates(DMDA_InternalFreeSurfaceOnRankZero,							&InternalFreeSurface_2D_coords);	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(coordDA_2D, 		InternalFreeSurface_2D_coords, 					&coors_2D ); 						CHKERRQ(ierr);

		for (i=0; i<M; i++){
			for (j=0; j<N; j++){

				// Update topography in 2D DA
				FreeSurface2D[j][i] = 	FreeSurface3D[0][j][i];

				// update advection velocity at free surface to 2D DA
				FreeSurfaceVx_2D[j][i] = 	FreeSurfaceVx_3D[0][j][i];
				FreeSurfaceVy_2D[j][i] = 	FreeSurfaceVy_3D[0][j][i];
				FreeSurfaceVz_2D[j][i] = 	FreeSurfaceVz_3D[0][j][i];


				// Update coordinates in 2D DA
				coors_2D[j][i].x 	=	coors_3D[0][j][i].x;
				coors_2D[j][i].y 	=	coors_3D[0][j][i].y;

			}
		}

		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 	InternalFreeSurfaceTopography, 	&FreeSurface2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurface_3D, 		FreeSurfaceAtRankZero, 			&FreeSurface3D ); 	CHKERRQ(ierr);


		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vx, 						&FreeSurfaceVx_2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVxAtRankZero, 				&FreeSurfaceVx_3D ); 	CHKERRQ(ierr);

		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vy, 						&FreeSurfaceVy_2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVyAtRankZero, 				&FreeSurfaceVy_3D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurfaceOnRankZero, 		FreeSurface_Vz, 						&FreeSurfaceVz_2D ); 	CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(DMDA_InternalFreeSurface_3D, 			FreeSurfaceVzAtRankZero, 				&FreeSurfaceVz_3D ); 	CHKERRQ(ierr);

		ierr = DMDAVecRestoreArray(coordDA_2D, 			InternalFreeSurface_2D_coords, 			&coors_2D ); 				CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(coordDA_3D, 			InternalFreeSurface_3D_coords_zero, 	&coors_3D ); 				CHKERRQ(ierr);

		// Set coordinates of 2D DA
		ierr = DMSetCoordinates(DMDA_InternalFreeSurfaceOnRankZero,InternalFreeSurface_2D_coords);	CHKERRQ(ierr);

	}

	/* Cleaning up */
	if (rank==0){
		ierr = DMDestroy(&DMDA_InternalFreeSurface_3D);			CHKERRQ(ierr);
	}
	ierr = VecScatterDestroy(&ctx_coords);				CHKERRQ(ierr);
	ierr = VecDestroy(&FreeSurfaceAtRankZero);			CHKERRQ(ierr);
	ierr = VecDestroy(&FreeSurfaceVxAtRankZero);		CHKERRQ(ierr);
	ierr = VecDestroy(&FreeSurfaceVyAtRankZero);		CHKERRQ(ierr);
	ierr = VecDestroy(&FreeSurfaceVzAtRankZero);		CHKERRQ(ierr);
	ierr = VecDestroy(&temp);							CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Apply erosion to the top surface of the model, using either a simple erosion law, or a more sophisticated
 * erosion code such as cascade */
#undef __FUNCT__
#define __FUNCT__ "ApplyErosion"
PetscErrorCode ApplyErosion( DM da_coord, UserContext *user )
{
	PetscMPIInt     rank, size;
	PetscErrorCode  ierr;
	DM			    cda;
	DMDACoor3d		***coors;
	Vec				gc, global;
	PetscInt		xs,ys,zs,xm,ym, zm, ix,iy,iz, nnode_x, nnode_y, nnode_z, iproc;
	PetscInt        num,ErosionMethod;
	double 			**SurfaceTopography_x,**SurfaceTopography_y,**SurfaceTopography_z;
	double 			*Data_RecvX, *Data_RecvY, *Data_RecvZ;
	double 			*SurfaceVec_x, *SurfaceVec_y, *SurfaceVec_z;

	MPI_Status      status;

	PetscPrintf(PETSC_COMM_WORLD,"  Applying erosion to the top surface of the model \n");

	// Error verification -------------
	if (user->DimensionalUnits != 1){
		PetscPrintf(PETSC_COMM_WORLD,"  *** Dimensionless units are employed; the erosion algorithm requires dimensional units ***  \n");
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	// --------------------------------


	/* --------------------------------------------------------------*/
	// Allocate 2D arrays that holds topography on every processor.
	// Note: this part is not optimal for parallelisation!
	ierr = DMDAGetInfo(da_coord, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);

	ierr = PetscMalloc( sizeof(double *)*(size_t)nnode_y, &SurfaceTopography_x ); CHKERRQ(ierr);
	for (iy=0; iy<nnode_y; iy++){
		ierr = PetscMalloc( sizeof(double )*(size_t)nnode_x, &SurfaceTopography_x[iy] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(double *)*(size_t)nnode_y, &SurfaceTopography_y ); CHKERRQ(ierr);
	for (iy=0; iy<nnode_y; iy++){
		ierr = PetscMalloc( sizeof(double )*(size_t)nnode_x, &SurfaceTopography_y[iy] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(double *)*(size_t)nnode_y, &SurfaceTopography_z ); CHKERRQ(ierr);
	for (iy=0; iy<nnode_y; iy++){
		ierr = PetscMalloc( sizeof(double )*(size_t)nnode_x, &SurfaceTopography_z[iy] ); CHKERRQ(ierr);
	}

	for (iy=0; iy<nnode_y; iy++){
		for (ix=0; ix<nnode_x; ix++){
			// initialize to zero
			SurfaceTopography_x[iy][ix] = 0;
			SurfaceTopography_y[iy][ix] = 0;
			SurfaceTopography_z[iy][ix] = 0;
		}
	}
	/* --------------------------------------------------------------*/

	// We have to reconstruct the top surface on processor 0
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	/* Extract topography of the local portion of the grid ----------*/
	ierr = DMGetCoordinateDM(da_coord,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da_coord,			&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,				&coors); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_coord,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	num 	=	0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				if (iz== (nnode_z-1)){
					SurfaceTopography_x[iy][ix] = ((double) coors[iz][iy][ix].x*user->Characteristic.Length);
					SurfaceTopography_y[iy][ix] = ((double) coors[iz][iy][ix].y*user->Characteristic.Length);
					SurfaceTopography_z[iy][ix] = ((double) coors[iz][iy][ix].z*user->Characteristic.Length);

				}

			}
		}
	}

	ierr = DMDAVecRestoreArray(cda,gc,	&coors); CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_coord,	&global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);
	/* --------------------------------------------------------------*/

	// Allocate required arrays and put 2D data into 1D array
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &SurfaceVec_x); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &SurfaceVec_y); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &SurfaceVec_z); CHKERRQ(ierr);
	num = 0;
	for (iy=0; iy<nnode_y; iy++){
		for(ix=0; ix<nnode_x; ix++){
			SurfaceVec_x[num] = SurfaceTopography_x[iy][ix];
			SurfaceVec_y[num] = SurfaceTopography_y[iy][ix];
			SurfaceVec_z[num] = SurfaceTopography_z[iy][ix];
			num 			  = num+1;
		}
	}


	/* Collect topography data on processor 0 -------------------------------*/
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &Data_RecvX); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &Data_RecvY); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)( nnode_x*nnode_y)*sizeof(double), &Data_RecvZ); CHKERRQ(ierr);


	if (size>1){
		if (rank==0) {

			for (iproc=1; iproc<size; iproc++){

				// Receive Data from the various processors & add to Data on proc 0
				ierr = MPI_Recv(Data_RecvX,nnode_x*nnode_y, MPI_DOUBLE, iproc, 11, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);
				ierr = MPI_Recv(Data_RecvY,nnode_x*nnode_y, MPI_DOUBLE, iproc, 12, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);
				ierr = MPI_Recv(Data_RecvZ,nnode_x*nnode_y, MPI_DOUBLE, iproc, 13, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);
				num = 0;
				for (num=0; num<nnode_x*nnode_y; num++){
					SurfaceVec_x[num] = SurfaceVec_x[num] + Data_RecvX[num];
					SurfaceVec_y[num] = SurfaceVec_y[num] + Data_RecvY[num];
					SurfaceVec_z[num] = SurfaceVec_z[num] + Data_RecvZ[num];
				}

			}
		}
		else{

			// Send data to 0. Seems that MPI works better if this data is in the form of a 1D array
			ierr = MPI_Send(SurfaceVec_x,nnode_x*nnode_y, MPI_DOUBLE, 0, 11, PETSC_COMM_WORLD); CHKERRQ(ierr);
			ierr = MPI_Send(SurfaceVec_y,nnode_x*nnode_y, MPI_DOUBLE, 0, 12, PETSC_COMM_WORLD); CHKERRQ(ierr);
			ierr = MPI_Send(SurfaceVec_z,nnode_x*nnode_y, MPI_DOUBLE, 0, 13, PETSC_COMM_WORLD); CHKERRQ(ierr);

		}
	}


	/* at this stage, all coordinate information should be present on processor 0 */
	if ((rank==0) && (1==0)){
		num = 0;
		for (iy=0; iy<nnode_y; iy++){
			for(ix=0; ix<nnode_x; ix++){
				PetscPrintf(PETSC_COMM_WORLD,"   Surface topography and coordinates at  [%lld,%lld]:= [%g,%g,%g] \n",iy, ix,  (LLD)SurfaceVec_x[num],   (LLD)SurfaceVec_y[num],  SurfaceVec_z[num]);
				num = num+1;
			}
		}
	}
	/* --------------------------------------------------------------*/

	/* Perform surface erosion (call cascade) -----------------------*/
	if (rank==0){
		ErosionMethod = 2;
		if (ErosionMethod==1){
			/*  Apply a simple addition to the surface, mainly for testing purposes */
			for (num=0; num<nnode_y*nnode_x; num++){
				SurfaceVec_z[num] = SurfaceVec_z[num] + 1*(SurfaceVec_x[num])/100;
			}

		}
		else if (ErosionMethod==2){
#ifdef USE_CASCADE

			/* Call Cascade wrapper (which is a fortran routine), under the requirement that
			 * the code is compiled with the cascade option activated */
			SecYear 			= 3600*24*365.25;

			// Cascade parameters:
			fluvial_erosion   	= user->fluvial_erosion;
			diffusion_erosion 	= user->diffusion_erosion;
			baselevelx0 		= user->baselevelx0;
			baselevely0 		= user->baselevely0;
			baselevelx1 		= user->baselevelx1;
			baselevely1 		= user->baselevely1;

			dt_real 	 = dt*user->Characteristic.Time/SecYear;		// in years
			wrappercascadef_(&nnode_x,&nnode_y, &dt_real, SurfaceVec_x, SurfaceVec_y,
					SurfaceVec_z, &fluvial_erosion, &diffusion_erosion,
					&baselevelx0, &baselevelx1, &baselevely0, &baselevely1);
#endif
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	/* --------------------------------------------------------------*/


	/* Send back Surface from proc 0 to all other processors --------*/
	MPI_Bcast(SurfaceVec_x,nnode_x*nnode_y, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(SurfaceVec_y,nnode_x*nnode_y, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(SurfaceVec_z,nnode_x*nnode_y, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

	num = 0;
	for (iy=0; iy<nnode_y; iy++){
		for(ix=0; ix<nnode_x; ix++){
			SurfaceTopography_x[iy][ix] = SurfaceVec_x[num];
			SurfaceTopography_y[iy][ix] = SurfaceVec_y[num];
			SurfaceTopography_z[iy][ix] = SurfaceVec_z[num];
			num            			  	= num+1;
		}
	}

	//DMDestroy(cda);
	//VecDestroy(gc);
	//VecDestroy(global);

	/* --------------------------------------------------------------*/

	/* Put topography back to the local portion of the grid ----------*/
	ierr = DMGetCoordinateDM(da_coord,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da_coord,			&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,				&coors); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da_coord,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){
				if (iz== (nnode_z-1)){
					coors[iz][iy][ix].z = SurfaceTopography_z[iy][ix]/user->Characteristic.Length;
				}
			}
		}
	}

	ierr = DMDAVecRestoreArray(cda,gc,	&coors); CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_coord,	&global); CHKERRQ(ierr);

	ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);
	
	
	/* --------------------------------------------------------------*/

	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);
//	ierr = VecDestroy(global); CHKERRQ(ierr);

	// Free allocated arrays
	ierr = PetscFree(SurfaceVec_x); CHKERRQ(ierr);
	ierr = PetscFree(SurfaceVec_y); CHKERRQ(ierr);
	ierr = PetscFree(SurfaceVec_z); CHKERRQ(ierr);
	ierr = PetscFree(Data_RecvX); CHKERRQ(ierr);
	ierr = PetscFree(Data_RecvY); CHKERRQ(ierr);
	ierr = PetscFree(Data_RecvZ); CHKERRQ(ierr);


	for (iy=0; iy<nnode_y; iy++){
		ierr =  PetscFree( SurfaceTopography_x[iy] ); CHKERRQ(ierr);
	}
	ierr = PetscFree( SurfaceVec_x ); CHKERRQ(ierr);

	for (iy=0; iy<nnode_y; iy++){
		ierr = PetscFree( SurfaceTopography_y[iy] ); CHKERRQ(ierr);
	}
	ierr = PetscFree( SurfaceVec_y ); CHKERRQ(ierr);

	for (iy=0; iy<nnode_y; iy++){
		ierr = PetscFree( SurfaceTopography_z[iy] ); CHKERRQ(ierr);
	}
	ierr = PetscFree( SurfaceTopography_z ); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/* http://howto.wikia.com/wiki/Howto_compare_floating_point_numbers_in_the_C_programming_language */
PetscBool PetscCompareScalar(PetscScalar f1, PetscScalar f2)
{
	PetscScalar precision = 1.0e-6;

	if( ((f1 - precision) < f2) && ((f1 + precision) > f2) ) {
		return PETSC_TRUE;
	}
	else {
		return PETSC_FALSE;
	}
}


/* Remesh the grid taking into account the free and the bottom surface topography ===========================
 *
 *It uses the DA_SurfaceTopography and DA_BottomTopography to create a new grid at the current PROC. As such
 * this algorithm is more scalable than the previous algorithm used in LaMEM, as that one reconstructed the
 * full free surface on every PROC.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "RemeshGrid"
PetscErrorCode RemeshGrid( DM da, UserContext *user )
{
	PetscMPIInt     rank;
	PetscErrorCode 	ierr;
	DM				cda_SurfaceTopo, cda;
	DMDACoor3d		***coors_SurfaceTopo, ***coors;
	Vec				gc_SurfaceTopo, gc, Surface_FIELD, Bottom_FIELD;
	PetscInt		xs_Z, ys_Z, zs_Z, xm_Z, ym_Z, zm_Z, nnode_x, nnode_y, nnode_z, nsurface_levels;
	PetscInt		ix,iy,iz, xs,ys,zs,xm,ym,zm;
	PetscScalar		x_left, y_front, x_right, y_back, x_left_loc, y_front_loc, x_right_loc, y_back_loc, MeanHeight_loc, MeanHeight;
	PetscScalar		dx,dy, dz, ***BOTTOM_ZCOORD, ***SURFACE_ZCOORD;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	/* Update x and y-coordinates of the mesh, based on the coordinates of the top surface
		 & send the info to all other PROCs */
	ierr = DMDAGetInfo	(user->DA_SurfaceTopography, 0, 0,0, &nsurface_levels , 0,0,0,0,0,0,0,0,0); 	CHKERRQ(ierr);
	ierr = DMDAGetInfo	(da, 						0, &nnode_x,&nnode_y, &nnode_z , 0,0,0,0,0,0,0,0,0); 	CHKERRQ(ierr);
	ierr = DMDAGetCorners (user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z); 		CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,			&cda_SurfaceTopo	); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,		&gc_SurfaceTopo		); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,			&coors_SurfaceTopo	); 	CHKERRQ(ierr);


	x_left_loc  = 0;	y_front_loc = 0;
	if ((xs_Z==0) && (ys_Z==0)  && (zs_Z+zm_Z==nsurface_levels) ){
		x_left_loc 	= coors_SurfaceTopo[nsurface_levels-1][0][0].x;
		y_front_loc = coors_SurfaceTopo[nsurface_levels-1][0][0].y;
	}
	ierr = MPI_Allreduce(&x_left_loc, &x_left ,1, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&y_front_loc,&y_front,1, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);


	x_right_loc  = 0;	y_back_loc = 0;
	if ((xs_Z+xm_Z==nnode_x) && (ys_Z+ym_Z==nnode_y)  && (zs_Z+zm_Z==nsurface_levels) ){
		x_right_loc	= coors_SurfaceTopo[nsurface_levels-1][ys_Z+ym_Z-1][xs_Z+xm_Z-1].x;
		y_back_loc 	= coors_SurfaceTopo[nsurface_levels-1][ys_Z+ym_Z-1][xs_Z+xm_Z-1].y;
	}
	ierr = MPI_Allreduce(&x_right_loc, &x_right ,1, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&y_back_loc,	&y_back	 ,1, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);

	user->x_left 	= 	x_left;
	user->y_front 	= 	y_front;
	user->W  		=	x_right - x_left;
	user->L 		=	y_back  - y_front;

	/* Compute the mean elevation of the free surface (before remeshing the free surface)*/
	MeanHeight_loc = 0;
	for (iy=ys_Z; iy<ys_Z+ym_Z; iy++){
		for(ix=xs_Z; ix<xs_Z+xm_Z; ix++){

			if (zs_Z+zm_Z==nsurface_levels){
				MeanHeight_loc = MeanHeight_loc + coors_SurfaceTopo[nsurface_levels-1][iy][ix].z/( (PetscScalar) (nnode_x*nnode_y));
			}

		}
	}
	ierr = MPI_Allreduce(&MeanHeight_loc, &MeanHeight ,1, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);


	/* Create a regular x,y grid for the complete grid*/
	ierr = DMGetCoordinateDM(da,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,		&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// make sure to use Ghost corners!!!


	dx 	 = user->W/((double) (nnode_x-1));
	dy 	 = user->L/((double) (nnode_y-1));
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){

				coors[iz][iy][ix].x   = ((double) (ix))*dx + user->x_left;
				coors[iz][iy][ix].y   = ((double) (iy))*dy + user->y_front;

			}
		}
	}

	/* Cleanup 'real' grid */
	ierr 	= 	DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);
	//ierr 	= 	DMDestroy(cda); CHKERRQ(ierr);
	ierr 	= 	DASetCoordinatesFromLocalVector( da, gc ); CHKERRQ(ierr);		// send to all ghost points
	//ierr 	= 	VecDestroy(gc); CHKERRQ(ierr);


	/* Interpolate the free surface from the previous grid to the new regular grid
	 *	There should typically not be a huge difference
	 */
	// IGNORED FOR NOW


	/* Update the free surface info on all other procs */
	// IGNORED FOR NOW



	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
	//ierr 	= 	VecDestroy(gc_SurfaceTopo); CHKERRQ(ierr);
	//ierr 	=	DMDestroy(cda_SurfaceTopo); CHKERRQ(ierr);



	/* Create a new mesh based on the surface and bottom topographies */
	ierr = DMGetLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, Surface_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography, Surface_FIELD, &SURFACE_ZCOORD );CHKERRQ(ierr);

	ierr = DMGetLocalVector( user->DA_BottomTopography, &Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( user->DA_BottomTopography, user->BottomTopography, INSERT_VALUES, Bottom_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_BottomTopography, Bottom_FIELD, &BOTTOM_ZCOORD );CHKERRQ(ierr);


	ierr 	= 	DMGetCoordinateDM(da,			&cda); CHKERRQ(ierr);
	ierr 	= 	DMGetCoordinatesLocal(da,		&gc); CHKERRQ(ierr);
	ierr 	= 	DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);
	ierr 	= 	DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// make sure to use Ghost corners!!!
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				/* We simply assume a regular grid in the z-direction.
				 * This can be changed in the future to refine the mesh near the the surface for example;
				 *
				 */
				dz 		=   (SURFACE_ZCOORD[zs_Z][iy][ix]-BOTTOM_ZCOORD[zs_Z][iy][ix])/((double) (nnode_z-1));

				coors[iz][iy][ix].z = ((double) (iz))*dz + BOTTOM_ZCOORD[zs_Z][iy][ix];


			}
		}
	}

	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography, Surface_FIELD,&SURFACE_ZCOORD ); 						CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalBegin(user->DA_SurfaceTopography,Surface_FIELD,INSERT_VALUES,user->SurfaceTopography);	CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalEnd  (user->DA_SurfaceTopography,Surface_FIELD,INSERT_VALUES,user->SurfaceTopography);	CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_SurfaceTopography, &Surface_FIELD );									CHKERRQ(ierr);

	ierr 	= 	DMDAVecRestoreArray(user->DA_BottomTopography, Bottom_FIELD,&BOTTOM_ZCOORD ); 						CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalBegin(user->DA_BottomTopography,Bottom_FIELD,INSERT_VALUES,user->BottomTopography);		CHKERRQ(ierr);
	ierr 	=	DMLocalToGlobalEnd  (user->DA_BottomTopography,Bottom_FIELD,INSERT_VALUES,user->BottomTopography);		CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector( user->DA_BottomTopography, &Bottom_FIELD );									CHKERRQ(ierr);


	ierr = DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);

	/*=========================================================================================
	 * Ensure that also ghost coordinates are send to all other procs.
	 * This is particularly an issue if random noise is set as initial perturbation.
	 */
	/*=========================================================================================*/
	ierr = DASetCoordinatesFromLocalVector( da, gc ); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);




	PetscFunctionReturn(0);

}
/* ==================================================================================================== */


/* Remesh the grid taking into account the free surface, fully parallelized ===========================
 *
 * NOTE:  disadvantages of this routine are:
 * 	1) The free surface is reconstructed on every PROC. This is not scalable for large grids.
 *  2) The new grid is based on the surface topography only.
 *
 *	In june 2011, we therefore created a new routine to take bottom and top topography into accont.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "RemeshGrid_old"
PetscErrorCode RemeshGrid_old( DM da, UserContext *user )
{
	PetscMPIInt     rank;
	PetscErrorCode  ierr;
	DM				cda;
	DMDACoor3d		***coors;
	Vec				gc;
	PetscInt		xs,ys,zs,xm,ym, zm, ix,iy,iz, i, j;
	PetscInt		nnode_x, nnode_y, nnode_z, num;
	PetscInt		iy_start, ix_start;
	PetscScalar		*surfx_vec, *surfy_vec, *surfz_vec, dx, dy;
	PetscScalar		*surfx_vec_loc, *surfy_vec_loc, *surfz_vec_loc;
	PetscScalar		**Surface_x, 		**Surface_y, 		**Surface_z;
	PetscScalar		**Surface_x_new, 	**Surface_y_new,	**Surface_z_new;
	PetscScalar		MeanHeight,			facx,				facy;
	PetscScalar		MeanHeightNew,		DiffHeight,			facz;
	PetscScalar		dz,					z_middle,			angle,		max_angle;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


	/* Collect the surface grid from all CPUs ------------------------*/
	ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);


	// Allocate stuff
	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfx_vec); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfy_vec); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfz_vec); CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfx_vec_loc); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfy_vec_loc); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(nnode_x*nnode_y)*sizeof(PetscScalar), &surfz_vec_loc); CHKERRQ(ierr);
	for (i=0; i<nnode_x*nnode_y; i++){
		surfx_vec_loc[i] = 0;	surfy_vec_loc[i] = 0; 	surfz_vec_loc[i] = 0;
	}

	// Put local portion of surface grid in vector (if available)
	ierr = DMGetCoordinateDM(da,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,		&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	if ((zs+zm)==nnode_z){
		/* CPU has a portion of the surface grid */
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){
				surfx_vec_loc[ix*nnode_y + iy] = coors[nnode_z-1][iy][ix].x;
				surfy_vec_loc[ix*nnode_y + iy] = coors[nnode_z-1][iy][ix].y;
				surfz_vec_loc[ix*nnode_y + iy] = coors[nnode_z-1][iy][ix].z;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);

	// Sum all local portions
	ierr =MPI_Allreduce(surfx_vec_loc,surfx_vec,nnode_x*nnode_y, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr =MPI_Allreduce(surfy_vec_loc,surfy_vec,nnode_x*nnode_y, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr =MPI_Allreduce(surfz_vec_loc,surfz_vec,nnode_x*nnode_y, MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);


	// Update x and y-coordinates (this is necessary in case we deform the grid and remesh, e.g. in a detachment folding setup) */
	user->x_left  	= 	surfx_vec[0];
	user->y_front 	=	surfy_vec[0];
	user->W  		= 	surfx_vec[nnode_x*nnode_y-1]- surfx_vec[0];		// model width
	user->L  		= 	surfy_vec[nnode_x*nnode_y-1]- surfy_vec[0];		// model length





	// compute mean height
	MeanHeight = 0;
	for (i=0; i<nnode_x*nnode_y; i++){
		MeanHeight	=	MeanHeight + surfz_vec[i]/((double) nnode_x*nnode_y);
	}

	// Allocate 2D arrays that contain the free surface coordinates
	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_x ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_x[j] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_y ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_y[j] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_z ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_z[j] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_x_new ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_x_new[j] ); CHKERRQ(ierr);
	}
	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_y_new ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_y_new[j] ); CHKERRQ(ierr);
	}

	ierr = PetscMalloc( sizeof(PetscScalar *)*(size_t)nnode_y, &Surface_z_new ); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){
		ierr = PetscMalloc( sizeof(PetscScalar )*(size_t)nnode_x, &Surface_z_new[j] ); CHKERRQ(ierr);
	}

	dx = user->W/((double) (nnode_x-1));
	dy = user->L/((double) (nnode_y-1));
	num = 0;
	for (ix=0; ix<nnode_x; ix++){
		for (iy=0; iy<nnode_y; iy++){
			Surface_x[iy][ix] 		= surfx_vec[num];
			Surface_y[iy][ix] 		= surfy_vec[num];
			Surface_z[iy][ix] 		= surfz_vec[num];

			Surface_x_new[iy][ix]   = ((double) (ix))*dx + user->x_left;
			Surface_y_new[iy][ix]   = ((double) (iy))*dy + user->y_front;
			Surface_z_new[iy][ix] 	= surfz_vec[num];
			num = num+1;
		}
	}


	/* Free arrays we don't need anymore */
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);
	ierr = PetscFree(surfx_vec_loc); CHKERRQ(ierr);	ierr = PetscFree(surfy_vec_loc); CHKERRQ(ierr);	ierr = PetscFree(surfz_vec_loc); CHKERRQ(ierr);
	ierr = PetscFree(surfx_vec); CHKERRQ(ierr);			ierr = PetscFree(surfy_vec); CHKERRQ(ierr);			ierr = PetscFree(surfz_vec); CHKERRQ(ierr);
	/*----------------------------------------------------------------*/

	/* Remesh the surface grid, using linear interpolation */
	/* front and back boundaries */
	for (ix=1; ix<nnode_x-1; ix++){
		if (Surface_x_new[0][ix]>=Surface_x[0][ix]){
			facx = 	(Surface_x_new[0][ix]-Surface_x[0][ix])/(Surface_x[0][ix+1]-Surface_x[0][ix]);
			Surface_z_new[0][ix] = (1-facx)*Surface_z[0][ix] + facx*Surface_z[0][ix+1];
		}
		else{
			facx = 	(Surface_x_new[0][ix]-Surface_x[0][ix-1])/(Surface_x[0][ix]-Surface_x[0][ix-1]);
			Surface_z_new[0][ix] = (1-facx)*Surface_z[0][ix-1] + (facx)*Surface_z[0][ix];
		}

		if (Surface_x_new[nnode_y-1][ix]>=Surface_x[nnode_y-1][ix]){
			facx = 	(Surface_x_new[nnode_y-1][ix]-Surface_x[nnode_y-1][ix])/(Surface_x[nnode_y-1][ix+1]-Surface_x[nnode_y-1][ix]);
			Surface_z_new[nnode_y-1][ix] = (1-facx)*Surface_z[nnode_y-1][ix] + facx*Surface_z[nnode_y-1][ix+1];
		}
		else{
			facx = 	(Surface_x_new[nnode_y-1][ix]-Surface_x[nnode_y-1][ix-1])/(Surface_x[nnode_y-1][ix]-Surface_x[nnode_y-1][ix-1]);
			Surface_z_new[nnode_y-1][ix] = (1-facx)*Surface_z[nnode_y-1][ix-1] + (facx)*Surface_z[nnode_y-1][ix];
		}
	}


	/* left and right boundaries */
	for (iy=1; iy<nnode_y-1; iy++){
		if (Surface_y_new[iy][0]>=Surface_y[iy][0]){
			facy = 	(Surface_y_new[iy][0]-Surface_y[iy][0])/(Surface_y[iy+1][0]-Surface_y[iy][0]);
			Surface_z_new[iy][0] = (1-facy)*Surface_z[iy][0] + facy*Surface_z[iy+1][0];
		}
		else{
			facy = 	(Surface_y_new[iy][0]-Surface_y[iy-1][0])/(Surface_y[iy][0]-Surface_y[iy-1][0]);
			Surface_z_new[iy][0] = (1-facy)*Surface_z[iy-1][0] + (facy)*Surface_z[iy][0];
		}
		if (Surface_y_new[iy][nnode_x-1]>=Surface_y[iy][nnode_x-1]){
			facy = 	(Surface_y_new[iy][nnode_x-1]-Surface_y[iy][nnode_x-1])/(Surface_y[iy+1][nnode_x-1]-Surface_y[iy][nnode_x-1]);
			Surface_z_new[iy][nnode_x-1] = (1-facy)*Surface_z[iy][nnode_x-1] + facy*Surface_z[iy+1][nnode_x-1];
		}
		else{
			facy = 	(Surface_y_new[iy][nnode_x-1]-Surface_y[iy-1][nnode_x-1])/(Surface_y[iy][nnode_x-1]-Surface_y[iy-1][nnode_x-1]);
			Surface_z_new[iy][nnode_x-1] = (1-facy)*Surface_z[iy-1][nnode_x-1] + (facy)*Surface_z[iy][nnode_x-1];
		}
	}

	/* Now remesh all other coordinates */
	for (ix=1; ix<nnode_x-1; ix++){
		for (iy=1; iy<nnode_y-1; iy++){
			if (Surface_y_new[iy][ix]>=Surface_y[iy][ix]){
				facy 		= 	(Surface_y_new[iy][ix]-Surface_y[iy  ][ix])/(Surface_y[iy+1][ix]-Surface_y[iy][ix]);
				iy_start 	= 	iy;
			}
			else{
				facy 		= 	(Surface_y_new[iy][ix]-Surface_y[iy-1][ix])/(Surface_y[iy][ix]-Surface_y[iy-1][ix]);
				iy_start 	= 	iy-1;
			}
			if (Surface_x_new[iy][ix]>=Surface_x[iy][ix]){
				facx 		= 	(Surface_x_new[iy][ix]-Surface_x[iy][ix]  )/(Surface_x[iy][ix+1]-Surface_x[iy][ix]);
				ix_start 	= 	ix;
			}
			else{
				facx 		= 	(Surface_x_new[iy][ix]-Surface_x[iy][ix-1])/(Surface_x[iy][ix]-Surface_x[iy][ix-1]);
				ix_start 	= 	ix-1;
			}
			Surface_z_new[iy][ix] = (1-facx)*(1-facy)*Surface_z[iy_start  ][ix_start  ] +
					(  facx)*(1-facy)*Surface_z[iy_start  ][ix_start+1] +
					(1-facx)*(  facy)*Surface_z[iy_start+1][ix_start  ] +
					(  facx)*(  facy)*Surface_z[iy_start+1][ix_start+1];

		}
	}

	/* Limit the maximum angle in x and y-direction that can occur on the surface (to avoid numerical problems) */
	max_angle = 0;
	for (ix=0; ix<nnode_x; ix++){
		for (iy=1; iy<nnode_y; iy++){
			dz        =     (Surface_z_new[iy][ix]-Surface_z_new[iy-1][ix]);
			dx        =     (Surface_x_new[iy][ix]-Surface_x_new[iy-1][ix]);
			dy 	      =     (Surface_y_new[iy][ix]-Surface_y_new[iy-1][ix]);
			z_middle  = 0.5*(Surface_z_new[iy][ix]+Surface_z_new[iy-1][ix]);

			angle 	  = atan(dz/dy);
			if (PetscAbsScalar(angle)>max_angle){ max_angle = PetscAbsScalar(angle); }	// store maximum surface angle

			if (PetscAbsScalar(angle/M_PI*180.0)>user->MaximumSurfaceAngle){
				// correct the surface since critical angle is exceeded
				dz = tan(user->MaximumSurfaceAngle/180.0*M_PI)*dy;
				Surface_z_new[iy-1][ix] = z_middle - 0.5*dz;
				Surface_z_new[iy  ][ix] = z_middle + 0.5*dz;
			}


		}

	}
	PetscPrintf(PETSC_COMM_WORLD,"# Maximum Surface Angle Measured in y-direction = %g  Maximum Allowed = %g\n", max_angle/M_PI*180.0, user->MaximumSurfaceAngle);

	max_angle = 0;
	for (iy=0; iy<nnode_y; iy++){
		for (ix=1; ix<nnode_x; ix++){
			dz        =     (Surface_z_new[iy][ix]-Surface_z_new[iy][ix-1]);
			dx 	      =     (Surface_x_new[iy][ix]-Surface_x_new[iy][ix-1]);
			dy 	      =     (Surface_y_new[iy][ix]-Surface_y_new[iy][ix-1]);

			z_middle  = 0.5*(Surface_z_new[iy][ix]+Surface_z_new[iy][ix-1]);

			angle 	  = atan(dz/dx);
			if (PetscAbsScalar(angle)>max_angle){ max_angle = PetscAbsScalar(angle); }	// store maximum surface angle

			if (PetscAbsScalar(angle/M_PI*180.0)>user->MaximumSurfaceAngle){
				// correct the surface since critical angle is exceeded
				dz = tan(user->MaximumSurfaceAngle/180.0*M_PI)*dx;
				Surface_z_new[iy ][ix-1] = z_middle - 0.5*dz;
				Surface_z_new[iy ][ix  ] = z_middle + 0.5*dz;
			}

		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"# Maximum Surface Angle Measured in x-direction = %g  Maximum Allowed = %g\n", max_angle/M_PI*180.0, user->MaximumSurfaceAngle);



	// compute mean height of new surface
	MeanHeightNew = 0;
	for (ix=0; ix<nnode_x; ix++){
		for (iy=0; iy<nnode_y; iy++){
			MeanHeightNew = MeanHeightNew + Surface_z_new[iy][ix]/((double) nnode_x*nnode_y);
		}
	}
	DiffHeight = MeanHeight-MeanHeightNew;
	for (ix=0; ix<nnode_x; ix++){
		for (iy=0; iy<nnode_y; iy++){
			Surface_z_new[iy][ix] = Surface_z_new[iy][ix]	+ DiffHeight;
		}
	}

	dz =     (MeanHeight-user->z_bot)/((double) (nnode_z-1));
	ierr = DMGetCoordinateDM(da,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,		&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// make sure to use Ghost corners!!!
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				// Detect a case
				facz = 0;
				if  (iz> (nnode_z - user->NumSurfaceNodes)){
					facz = ((double)(iz-(nnode_z - user->NumSurfaceNodes))) /((double) (user->NumSurfaceNodes - 1));
				}

				coors[iz][iy][ix].x = Surface_x_new[iy][ix];
				coors[iz][iy][ix].y = Surface_y_new[iy][ix];
				coors[iz][iy][ix].z = ((double) (iz))*dz + user->z_bot + facz*(Surface_z_new[iy][ix]-MeanHeightNew);



			}
		}
	}

	ierr = DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);

	//ierr = DMDestroy(cda); CHKERRQ(ierr);

	/*=========================================================================================
	 * Ensure that also ghost coordinates are send to all other procs.
	 * This is particularly an issue if random noise is set as initial perturbation.
	 */
	/*=========================================================================================*/
	ierr = DASetCoordinatesFromLocalVector( da, gc ); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);

	/* Free allocated arrays */
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_x[j]);		CHKERRQ(ierr);	}	ierr=PetscFree(Surface_x		); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_y[j]);		CHKERRQ(ierr);	}	ierr=PetscFree(Surface_y		); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_z[j]);		CHKERRQ(ierr);	}	ierr=PetscFree(Surface_z		); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_x_new[j]);	CHKERRQ(ierr);	}	ierr=PetscFree(Surface_x_new	); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_y_new[j]);	CHKERRQ(ierr);	}	ierr=PetscFree(Surface_y_new	); CHKERRQ(ierr);
	for (j=0; j<nnode_y; j++){ ierr = PetscFree(Surface_z_new[j]);	CHKERRQ(ierr);	}	ierr=PetscFree(Surface_z_new	); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/*===================================================================================================== */
/* Compute the number of nodes on each processor, given a DM that defines the number of
 * elements in each direction on each processor. */
#undef __FUNCT__
#define __FUNCT__ "ComputeNeighbors"
PetscErrorCode ComputeNeighbors( DM da, PetscInt NeighborCPU[3][3][3] )
{
	const PetscMPIInt		*ranks;
	PetscInt 				i,ix,iy,iz;

	// get the neighbors of the current CPU
	DMDAGetNeighbors(da,&ranks);
	i = 0;
	for (iz=0; iz<3; iz++){
		for (iy=0; iy<3; iy++){
			for (ix=0; ix<3; ix++){
				NeighborCPU[iy][ix][iz] = ranks[i];
				i = i + 1;
			}
		}
	}

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/*===================================================================================================== */
/* Sets uniform coordinates for FDSTAG meshes															*/
#undef __FUNCT__
#define __FUNCT__ "SetUniformCoordinates_FDSTAG"
PetscErrorCode SetUniformCoordinates_FDSTAG( UserContext *user )
{
	PetscErrorCode ierr;
	PetscScalar		dx,dy,dz,xmin,xmax,ymin,ymax,zmin,zmax;

	dx 		=	((PetscScalar) (user->W))/((PetscScalar) (user->finest_nnode_x-1) );
	dy 		=	((PetscScalar) (user->L))/((PetscScalar) (user->finest_nnode_y-1) );
	dz		=	((PetscScalar) (user->H))/((PetscScalar) (user->finest_nnode_z-1) );
	xmin 	=	((PetscScalar) user->x_left) 			;
	xmax 	=	((PetscScalar) user->x_left) 	+ ((PetscScalar) user->W);
	ymin 	=	((PetscScalar) user->y_front)			;
	ymax 	=	((PetscScalar) user->y_front) + ((PetscScalar) user->L);
	zmin 	=	((PetscScalar) user->z_bot)				;
	zmax 	=	((PetscScalar) user->z_bot)	+ ((PetscScalar) user->H);

	/* Set coordinates of DM (for plotting later) */
	ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_CORNER,  xmin       ,xmax       ,ymin,		ymax,		zmin				 ,	zmax							);	CHKERRQ(ierr);

	if ((ymin+dy/2)<(ymax-dy/2.0)){
		ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_CENTER,    xmin+dx/2.0,xmax-dx/2.0,ymin+dy/2.0,ymax-dy/2.0,zmin+dz/2.0	,	zmax-dz/2.0);	CHKERRQ(ierr);
		ierr = 	DMDASetUniformCoordinates(user->DA_Pres, 		xmin+dx/2.0,xmax-dx/2.0,ymin+dy/2.0,ymax-dy/2.0,zmin+dz/2.0	,	zmax-dz/2.0);	CHKERRQ(ierr);
		ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_XZ_POINTS, xmin,xmax,ymin+dy/2.0,ymax-dy/2.0,				zmin					,	zmax							); CHKERRQ(ierr);
	}
	else {
		// if we have 1 element in y-direction, we might have issues if we don't do this (at least in PETSc 3.1)
		//

		ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_CENTER,    xmin+dx/2.0,xmax-dx/2.0,ymin,ymax,zmin+dz/2.0	,	zmax-dz/2.0);	CHKERRQ(ierr);
		ierr = 	DMDASetUniformCoordinates(user->DA_Pres, 		xmin+dx/2.0,xmax-dx/2.0,ymin,ymax,zmin+dz/2.0	,	zmax-dz/2.0);	CHKERRQ(ierr);

		ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_XZ_POINTS, xmin,xmax,ymin+dy/2.0,ymax,			zmin					,	zmax							); CHKERRQ(ierr);

	}

	ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_XY_POINTS, xmin,xmax,ymin,ymax,							zmin+dz/2.0	,	zmax-dz/2.0); CHKERRQ(ierr);
	ierr = 	DMDASetUniformCoordinates(user->FDSTAG.DA_YZ_POINTS, xmin+dx/2.0,xmax-dx/2.0,ymin,ymax,				zmin					,	zmax							); CHKERRQ(ierr);



	ierr = 	DMDASetUniformCoordinates(user->DA_Vel,  	xmin       ,xmax       ,ymin,		ymax,		zmin				 ,	zmax	);	CHKERRQ(ierr);


	// check this ! It causes an error when VTKOutpufiles == 0 since DA_Quadrature is only created when VTK_Outputfile==1 (in LaMEMInitialze l.985)
	// I don't know whether this needs to be here. tobi

	if (user->VTKOutputFiles==1){
			ierr = 	DMDASetUniformCoordinates(user->DA_Quadrature,  			xmin       ,xmax       ,ymin,		ymax,		zmin				 ,	zmax		);	CHKERRQ(ierr);
	}




	PetscFunctionReturn(0);
}
/* ==================================================================================================== */



/*===================================================================================================== */
/* Recreates the FDSTAG mesh based on the BG strainrate that is employed.								*/
#undef __FUNCT__
#define __FUNCT__ "DeformFDSTAGMeshWithBackgroundStrainrate"
PetscErrorCode DeformFDSTAGMeshWithBackgroundStrainrate( UserContext *user )
{
	PetscErrorCode ierr;
	PetscScalar		xmin,xmax,ymin,ymax,zmin,zmax;

	// This function is correct but terribly implicit.
	// Coordinate origin is constrained (has zero velocity in this settings).
	// The coordinate origin MUST be within domain!
	// Otherwise the domain will start moving out of the origin.

	xmin 	=	user->x_left 			;
	xmax 	=	user->x_left 	+ user->W;
	ymin 	=	user->y_front			;
	ymax 	=	user->y_front 	+ user->L;
	zmin 	=	user->z_bot				;
	zmax 	=	user->z_bot	 	+ user->H;

	// advect borders of grid based on BG strain rate
	xmin	=	xmin - user->dt*xmin*user->BC.Exx;
	xmax	=	xmax - user->dt*xmax*user->BC.Exx;
	ymin	=	ymin - user->dt*ymin*user->BC.Eyy;
	ymax	=	ymax - user->dt*ymax*user->BC.Eyy;
	zmin	=	zmin + user->dt*zmin*(user->BC.Exx + user->BC.Eyy);
	zmax	=	zmax + user->dt*zmax*(user->BC.Exx + user->BC.Eyy);

	/* updated coords */
	user->x_left  	= xmin;
	user->W  		= xmax-xmin;
	user->y_front  	= ymin;
	user->L 		= ymax-ymin;
	user->z_bot  	= zmin;
	user->H 		= zmax-zmin;

	/* recreate coordinates of FDSTAG mesh */
	ierr = SetUniformCoordinates_FDSTAG( user ); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */



