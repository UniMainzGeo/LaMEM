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

Output.c, contains the following subroutines:

WriteOutputFileMatlab			-	Write MATLAB output file for each CPU

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Output.h"
#include "Material.h"



/*==========================================================================================================*/
/* Write a MATLAB output file for every CPU - Matlab than puts the pieces together */
PetscErrorCode WriteOutputFileMatlab( UserContext *user, DM da_nodes, DM da_temp, Vec Velocity, Vec Temp, PetscInt itime, const PetscInt ElementType, const PetscInt ngp_vel, const char DirectoryName[])
{
	PetscMPIInt         rank, size;
	PetscErrorCode      ierr;
	PetscInt			mx,my,mz,cpu_x,cpu_y,cpu_z, i, phase;
	PetscInt			nel_x,nel_y,nel_z, xs,ys,zs,xm,ym,zm, intp;
	PetscInt			ix,iy,iz,num, iphase, istress;
	char				SaveFileName[PETSC_MAX_PATH_LEN];
	Vec					information, characteristic, coord, coord_x, coord_y, coord_z, local_Vel, Vx, Vy, Vz;
	Vec					mu, rho, intpx, intpy, intpz, T, Temp_local, G, C, k, Cp, Q, alpha, FK, phi, n, Phases;
	Vec					P, DevStress, DevStrainrate, Txx,Tyy,Tzz,Txy,Tyz,Txz, Exx,Eyy,Ezz,Exy,Eyz,Exz, NumParticles;
	Vec 				Strain, PlasticStrain;
	Mat					TimeDependentDataArray;
	PetscViewer			view_out;
	DM					cda;
	DMDACoor3d		 	***coords;
	PetscScalar			*coords_x, *coords_y, *coords_z, *Vx_array, *Vy_array, *Vz_array, *mu_array, *rho_array;
	PetscScalar			*intpx_array, *intpy_array, *intpz_array, *Temp_array, *NumParticles_array;
	PetscScalar			*G_array, *C_array, *k_array, *Cp_array, *Q_array, *Alpha_array, *FK_array, *n_array, *phi_array;
	PetscScalar			*Phase_array, *Pressure_array, *DevStrainrate_array, *DevStress_array;
	PetscScalar			*Txx_array, *Tyy_array, *Tzz_array, *Txz_array, *Txy_array, *Tyz_array;
	PetscScalar			*Exx_array, *Eyy_array, *Ezz_array, *Exz_array, *Exy_array, *Eyz_array;
	PetscScalar			*Strain_array, *PlasticStrain_array;
	Field 				***velocity;
	//MaterialsElement	***materials;
	PetscScalar			***temperature;
	GlobalTimeDependentData	TimeDependentData;
	PetscScalar ***materials_array;
	MaterialsElementDynamic material_data;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);



	sprintf(SaveFileName,"%s/%s.%lld.%lld.out",DirectoryName, user->OutputFile,(LLD)rank,(LLD)itime+1000000LL);  // construct the filename
	ierr = DMDAGetInfo(da_nodes, 0, &mx,    &my, &mz, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);

	ierr = DMDAGetInfo(user->DA_Processors, 0, &nel_x, &nel_y, &nel_z, &cpu_x,&cpu_y,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_SELF,&information); CHKERRQ(ierr);
	ierr = VecSetSizes(information,PETSC_DECIDE,40); CHKERRQ(ierr);
	ierr = VecSetFromOptions(information); CHKERRQ(ierr);
	ierr = VecSetValue(information,0 ,(PetscScalar)(user->finest_nnode_x), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,1 ,(PetscScalar)(user->finest_nnode_y), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,2 ,(PetscScalar)(user->finest_nnode_z), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,3 ,(PetscScalar)(nel_x), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,4 ,(PetscScalar)(nel_y), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,5 ,(PetscScalar)(nel_z), 	INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,6 ,user->mumax*user->Characteristic.Viscosity, 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,7 ,user->ampl2D*user->Characteristic.Length, 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,8 ,user->ampl3D*user->Characteristic.Length, 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,9 ,(PetscScalar)(ElementType), 				INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,10,(PetscScalar)(ngp_vel),INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,11,user->BC.Exx*user->Characteristic.Strainrate, 		INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,12,user->BC.Eyy*user->Characteristic.Strainrate, 		INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,13,NumMaterialPropsElem ,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,14,(PetscScalar)(cpu_x)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,15,(PetscScalar)(cpu_y)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,16,(PetscScalar)(cpu_z)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,17,(PetscScalar)(size)	,INSERT_VALUES); CHKERRQ(ierr);


	/* Save information about the local nodes-based arrays */
	ierr = DMDAGetCorners(da_nodes,      &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
	if ( (user->BC.LeftBound==3) && ((xs+xm)==user->finest_nnode_x-1) ) {
		xm = 	xm+1;
	}
	if ( (user->BC.FrontBound==3) && ((ys+ym)==user->finest_nnode_y-1) ) {
		ym =	ym+1;
	}

	ierr = VecSetValue(information,18,(PetscScalar)(xs)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,19,(PetscScalar)(ys)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,20,(PetscScalar)(zs)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,21,(PetscScalar)(xm)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,22,(PetscScalar)(ym)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,23,(PetscScalar)(zm)	,INSERT_VALUES); CHKERRQ(ierr);

	/* Save information about the local element-based arrays */
	ierr = DMDAGetCorners(user->DA_Processors,      &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
	ierr = VecSetValue(information,24,(PetscScalar)(xs)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,25,(PetscScalar)(ys)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,26,(PetscScalar)(zs)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,27,(PetscScalar)(xm)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,28,(PetscScalar)(ym)	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,29,(PetscScalar)(zm)	,INSERT_VALUES); CHKERRQ(ierr);

	/* time, timestep and angle of gravity */
	ierr = VecSetValue(information,30,user->time*user->Characteristic.Time	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,31,user->dt  *user->Characteristic.Time	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,32,user->GravityAngle	,INSERT_VALUES); CHKERRQ(ierr);





	ierr = VecAssemblyBegin(information); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(information); CHKERRQ(ierr);

	/* Save characteristic values */
	ierr = VecCreate(PETSC_COMM_SELF,&characteristic); CHKERRQ(ierr);
	ierr = VecSetSizes(characteristic,PETSC_DECIDE,40); CHKERRQ(ierr);
	ierr = VecSetFromOptions(characteristic); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,0, user->Characteristic.cmYear				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,1, user->Characteristic.Myrs					,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,2, user->Characteristic.MPa					,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,3, user->Characteristic.SecYear				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,4, user->Characteristic.km					,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,5, user->Characteristic.ThermalExpansivity	,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,6, user->Characteristic.Strainrate			,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,7, user->Characteristic.kg					,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,8, user->Characteristic.Density				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,9, user->Characteristic.Viscosity			,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,10,user->Characteristic.Temperature			,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,11,user->Characteristic.Velocity				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,12,user->Characteristic.Stress				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,13,user->Characteristic.Time					,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(characteristic,14,user->Characteristic.Length				,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(characteristic); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(characteristic); CHKERRQ(ierr);

	/* Create local vectors that contains coordinates and velocity @ nodal points */
	ierr = DMDAGetCorners(da_nodes,      &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
	if ( (user->BC.LeftBound==3) && ((xs+xm)==user->finest_nnode_x-1) ) {
		xm = 	xm+1;
	}
	if ( (user->BC.FrontBound==3) && ((ys+ym)==user->finest_nnode_y-1) ) {
		ym =	ym+1;
	}


	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_x); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_y); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&coords_z); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&Vx_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&Vy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&Vz_array); CHKERRQ(ierr);
#ifdef TEMPERATURE
	ierr = PetscMalloc((size_t)(xm*ym*zm)*sizeof(PetscScalar),&Temp_array); CHKERRQ(ierr);
#endif


	/* Get the coordinates from the DMDA and put them into vectors */
	ierr = DMGetCoordinateDM(da_nodes,&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da_nodes,&coord); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&coords); CHKERRQ(ierr);
	num = 0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){
				coords_x[num]   = coords[iz][iy][ix].x;
				coords_y[num]   = coords[iz][iy][ix].y;
				coords_z[num]   = coords[iz][iy][ix].z;

				// check for periodicity

				if ( (user->BC.LeftBound==3) && (ix==user->finest_nnode_x-1) ) {
					coords_x[num] = coords_x[num] + user->W;
				}
				if ( (user->BC.FrontBound==3) && (iy==user->finest_nnode_y-1) ) {
					coords_y[num] = coords_y[num] + user->L;
				}


				coords_x[num]   = coords_x[num]*user->Characteristic.Length;
				coords_y[num]   = coords_y[num]*user->Characteristic.Length;
				coords_z[num]   = coords_z[num]*user->Characteristic.Length;



				num			  	= num+1;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&coords); CHKERRQ(ierr);
	//ierr = VecDestroy(&coord); CHKERRQ(ierr);
	//ierr = DMDestroy(&cda); CHKERRQ(ierr);


	ierr = DMGetLocalVector(da_nodes,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da_nodes, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_nodes,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_nodes, local_Vel,	&velocity); CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = DMGetLocalVector(da_temp,&Temp_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_temp, Temp_local,	&temperature); CHKERRQ(ierr);
#endif


	num = 0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){
				//coords_x[num]   = coords[iz][iy][ix].x*user->Characteristic.Length;
				//coords_y[num]   = coords[iz][iy][ix].y*user->Characteristic.Length;
				//coords_z[num]   = coords[iz][iy][ix].z*user->Characteristic.Length;
				Vx_array[num]   = velocity[iz][iy][ix].Vx*user->Characteristic.Velocity;
				Vy_array[num]   = velocity[iz][iy][ix].Vy*user->Characteristic.Velocity;
				Vz_array[num]   = velocity[iz][iy][ix].Vz*user->Characteristic.Velocity;
#ifdef TEMPERATURE
				Temp_array[num] = temperature[iz][iy][ix]*user->Characteristic.Temperature;
#endif

				num			  = num+1;

			}
		}
	}
	ierr = DMDAVecRestoreArray(da_nodes, local_Vel,	&velocity); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da_nodes,&local_Vel); CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = DMDAVecRestoreArray(da_temp, Temp_local,	&temperature); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da_temp,&Temp_local); CHKERRQ(ierr);
#endif

	//ierr = DMDAVecRestoreArray(cda,coord,&coords); CHKERRQ(ierr);

	/* Store the time-averaged data in arrays */
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,8 + max_num_phases*6 + 2*6 + 3,itime,itime,PETSC_NULL,&TimeDependentDataArray); CHKERRQ(ierr);
	//	MatCreateSeqAIJ(PETSC_COMM_SELF,8 + 2*6,user->time_end,user->time_end,PETSC_NULL,&TimeDependentDataArray);
	for (i=0; i<itime; i++){
		TimeDependentData 	= 	user->TimeDependentData[i];
		ierr = MatSetValue(TimeDependentDataArray,0,i,TimeDependentData.Time, 		INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,1,i,TimeDependentData.Vrms, 		INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,2,i,TimeDependentData.Vx_max, 	INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,3,i,TimeDependentData.Vx_min, 	INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,4,i,TimeDependentData.Vy_max, 	INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,5,i,TimeDependentData.Vy_min, 	INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,6,i,TimeDependentData.Vz_max, 	INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValue(TimeDependentDataArray,7,i,TimeDependentData.Vz_min, 	INSERT_VALUES); CHKERRQ(ierr);

		for (iphase=0; iphase<max_num_phases; iphase++){
			num = 8 + iphase*6;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MinXCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MaxXCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MinYCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MaxYCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MinZCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
			ierr = MatSetValue(TimeDependentDataArray,num, i,TimeDependentData.MaxZCoordPhase[iphase], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
		}

		// store average stress & average strain rate
		for (istress=0; istress<6; istress++){
			ierr = MatSetValue(TimeDependentDataArray,num,i,TimeDependentData.DevStress[istress], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
		}
		for (istress=0; istress<6; istress++){
			ierr = MatSetValue(TimeDependentDataArray,num,i,TimeDependentData.DevStrainrate[istress], 	INSERT_VALUES); CHKERRQ(ierr);
			num=num+1;
		}

		// min, max & average topography
		ierr = MatSetValue(TimeDependentDataArray,num,i,TimeDependentData.MinTopography, 	INSERT_VALUES); CHKERRQ(ierr);
		num=num+1;
		ierr = MatSetValue(TimeDependentDataArray,num,i,TimeDependentData.MaxTopography, 	INSERT_VALUES); CHKERRQ(ierr);
		num=num+1;
		ierr = MatSetValue(TimeDependentDataArray,num,i,TimeDependentData.MeanTopography, 	INSERT_VALUES); CHKERRQ(ierr);



	}
	ierr = MatAssemblyBegin(TimeDependentDataArray,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(TimeDependentDataArray,		MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,coords_x,&coord_x); CHKERRQ(ierr);		ierr = VecAssemblyBegin(coord_x); CHKERRQ(ierr);	ierr = VecAssemblyEnd(coord_x); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,coords_y,&coord_y); CHKERRQ(ierr);		ierr = VecAssemblyBegin(coord_y); CHKERRQ(ierr);	ierr = VecAssemblyEnd(coord_y); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,coords_z,&coord_z); CHKERRQ(ierr);		ierr = VecAssemblyBegin(coord_z); CHKERRQ(ierr);	ierr = VecAssemblyEnd(coord_z); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,Vx_array,  &Vx     ); CHKERRQ(ierr);	ierr = VecAssemblyBegin(Vx); CHKERRQ(ierr);				ierr = VecAssemblyEnd(Vx); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,Vy_array,  &Vy     ); CHKERRQ(ierr);	ierr = VecAssemblyBegin(Vy); CHKERRQ(ierr);				ierr = VecAssemblyEnd(Vy); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,Vz_array,  &Vz     ); CHKERRQ(ierr);	ierr = VecAssemblyBegin(Vz); CHKERRQ(ierr);				ierr = VecAssemblyEnd(Vz); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xm*ym*zm,Temp_array,&T      ); CHKERRQ(ierr);	ierr = VecAssemblyBegin(T); CHKERRQ(ierr);				ierr = VecAssemblyEnd(T); CHKERRQ(ierr);

	/* Create local vectors that contain element-wise properties @ integration points */
	ierr = DMDAGetCorners(user->DA_Processors,  &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&mu_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&rho_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&G_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&C_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&k_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Cp_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Q_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Alpha_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&FK_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&n_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&phi_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Phase_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Pressure_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&DevStrainrate_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&DevStress_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Txx_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Tyy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Tzz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Txz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Tyz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Txy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Exx_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Eyy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Ezz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Exz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Eyz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Exy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&intpx_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&intpy_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&intpz_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&NumParticles_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&Strain_array); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&PlasticStrain_array); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);

	num = 0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, ix-xs,iy-ys,iz-zs, ngp_vel, materials_array );

				for (intp=0; intp<ngp_vel; intp++){
					Phase_array[num] =  material_data.Phases[intp];

					mu_array[num]  	 =  material_data.Viscosity[intp]*user->Characteristic.Viscosity;
					rho_array[num] 	 =  material_data.Density[intp]*user->Characteristic.Density;
					intpx_array[num] =  material_data.Coord[0][intp]*user->Characteristic.Length;
					intpy_array[num] =  material_data.Coord[1][intp]*user->Characteristic.Length;
					intpz_array[num] =  material_data.Coord[2][intp]*user->Characteristic.Length;
					G_array[num]  	 =  material_data.ElasticShearModule[intp]*user->Characteristic.Stress;


					NumParticles_array[num] 	= material_data.NumParticles[intp];
					Pressure_array[num] 		= material_data.Pressure[intp];
					DevStrainrate_array[num] 	= material_data.SecondInvariantDevStrainrate[intp];
					DevStress_array[num] 		= material_data.SecondInvariantDevStress[intp];
					Txx_array[num] 				= material_data.DevStress[0][intp];
					Tyy_array[num] 				= material_data.DevStress[1][intp];
					Tzz_array[num] 				= material_data.DevStress[2][intp];
					Txy_array[num] 				= material_data.DevStress[3][intp];
					Txz_array[num] 				= material_data.DevStress[4][intp];
					Tyz_array[num] 				= material_data.DevStress[5][intp];
					Exx_array[num] 				= material_data.DevStrainrate[0][intp];
					Eyy_array[num] 				= material_data.DevStrainrate[1][intp];
					Ezz_array[num] 				= material_data.DevStrainrate[2][intp];
					Exy_array[num] 				= material_data.DevStrainrate[3][intp];
					Exz_array[num] 				= material_data.DevStrainrate[4][intp];
					Eyz_array[num] 				= material_data.DevStrainrate[5][intp];

					Strain_array[num] 			= material_data.Strain[intp];
					PlasticStrain_array[num] 	= material_data.PlasticStrain[intp];


					// Look-up quantities (not sure we should keep outputting them; I leave it for now for compatibility reasons)
					phase 						= 	( (PetscInt) Phase_array[num] );
					n_array[num]  	 			=  user->PhaseProperties.n_exponent[phase];
					C_array[num]  	 			=  user->PhaseProperties.Cohesion[phase]*user->Characteristic.Stress;
					phi_array[num]   			=  user->PhaseProperties.FrictionAngle[phase];
					k_array[num]	 			=  user->PhaseProperties.T_Conductivity[phase]*user->Characteristic.T_conductivity;
					Cp_array[num]	 			=  user->PhaseProperties.HeatCapacity[phase] *user->Characteristic.HeatCapacity;
					Q_array[num]				=  user->PhaseProperties.RadioactiveHeat[phase] *user->Characteristic.RadioactiveHeat;
					Alpha_array[num] 			=  user->PhaseProperties.ThermalExpansivity[phase]*(1.0/user->Characteristic.Temperature);
					FK_array[num]	 			=  user->PhaseProperties.FrankKamenetskii[phase]*(1.0/user->Characteristic.Temperature);

					num			  	 = num+1;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);

	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,mu_array, &mu);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(mu);  CHKERRQ(ierr);				ierr = VecAssemblyEnd(mu);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,rho_array,&rho);  CHKERRQ(ierr);			ierr = VecAssemblyBegin(rho);  CHKERRQ(ierr);				ierr = VecAssemblyEnd(rho);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,n_array,	 &n);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(n);  CHKERRQ(ierr);					ierr = VecAssemblyEnd(n);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,G_array,  &G);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(G);  CHKERRQ(ierr);					ierr = VecAssemblyEnd(G);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,C_array,  &C);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(C);  CHKERRQ(ierr);					ierr = VecAssemblyEnd(C);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,phi_array,&phi);  CHKERRQ(ierr);			ierr = VecAssemblyBegin(phi);  CHKERRQ(ierr);				ierr = VecAssemblyEnd(phi);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,k_array,  &k);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(k);  CHKERRQ(ierr);					ierr = VecAssemblyEnd(k);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Cp_array, &Cp);  CHKERRQ(ierr);				ierr = VecAssemblyBegin(Cp);  CHKERRQ(ierr);				ierr = VecAssemblyEnd(Cp);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Q_array, &Q);  CHKERRQ(ierr);					ierr = VecAssemblyBegin(Q);  CHKERRQ(ierr);					ierr = VecAssemblyEnd(Q);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Alpha_array, &alpha);  CHKERRQ(ierr); ierr = VecAssemblyBegin(alpha);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(alpha);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,FK_array, 	&FK);  CHKERRQ(ierr);			ierr = VecAssemblyBegin(FK);  CHKERRQ(ierr);				ierr = VecAssemblyEnd(FK);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,intpx_array,&intpx);  CHKERRQ(ierr);	ierr = VecAssemblyBegin(intpx);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(intpx);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,intpy_array,&intpy);  CHKERRQ(ierr);	ierr = VecAssemblyBegin(intpy);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(intpy);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,intpz_array,&intpz);  CHKERRQ(ierr);	ierr = VecAssemblyBegin(intpz);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(intpz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Phase_array,&Phases);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Phases);  CHKERRQ(ierr);		ierr = VecAssemblyEnd(Phases);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,NumParticles_array,&NumParticles);  CHKERRQ(ierr); ierr = VecAssemblyBegin(NumParticles);  CHKERRQ(ierr);		ierr = VecAssemblyEnd(NumParticles);  CHKERRQ(ierr);


	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Pressure_array,&P);  CHKERRQ(ierr); ierr = VecAssemblyBegin(P);  CHKERRQ(ierr);ierr = VecAssemblyEnd(P);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,DevStress_array,&DevStress);  CHKERRQ(ierr); ierr = VecAssemblyBegin(DevStress);  CHKERRQ(ierr);ierr = VecAssemblyEnd(DevStress);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,DevStrainrate_array,&DevStrainrate);  CHKERRQ(ierr); ierr = VecAssemblyBegin(DevStrainrate);  CHKERRQ(ierr);ierr = VecAssemblyEnd(DevStrainrate);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Txx_array,&Txx);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Txx);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Txx);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Tyy_array,&Tyy);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Tyy);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Tyy);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Tzz_array,&Tzz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Tzz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Tzz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Txz_array,&Txz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Txz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Txz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Tyz_array,&Tyz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Tyz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Tyz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Txy_array,&Txy);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Txy);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Txy);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Exx_array,&Exx);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Exx);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Exx);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Eyy_array,&Eyy);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Eyy);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Eyy);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Ezz_array,&Ezz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Ezz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Ezz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Exz_array,&Exz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Exz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Exz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Eyz_array,&Eyz);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Eyz);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Eyz);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Exy_array,&Exy);  CHKERRQ(ierr); ierr = VecAssemblyBegin(Exy);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Exy);  CHKERRQ(ierr);

	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,Strain_array,&Strain);  CHKERRQ(ierr); 				ierr = VecAssemblyBegin(Strain);  CHKERRQ(ierr);ierr = VecAssemblyEnd(Strain);  CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,PlasticStrain_array,&PlasticStrain);  CHKERRQ(ierr); ierr = VecAssemblyBegin(PlasticStrain);  CHKERRQ(ierr);ierr = VecAssemblyEnd(PlasticStrain);  CHKERRQ(ierr);



	/* Write the actual file */
//	ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_SELF,SaveFileName,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileName,FILE_MODE_WRITE,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(view_out, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

	ierr = VecView(information, 		view_out); CHKERRQ(ierr);

	ierr = VecView(coord_x, 				view_out); CHKERRQ(ierr); 		// save coordinate array
	ierr = VecView(coord_y, 				view_out); CHKERRQ(ierr); 		// save coordinate array
	ierr = VecView(coord_z, 				view_out); CHKERRQ(ierr); 		// save coordinate array
	ierr = VecView(Vx,    					view_out); CHKERRQ(ierr); 		// save velocity   array
	ierr = VecView(Vy,    					view_out); CHKERRQ(ierr); 		// save velocity   array
	ierr = VecView(Vz,    					view_out); CHKERRQ(ierr); 		// save velocity   array
	ierr = VecView(mu,    					view_out); CHKERRQ(ierr); 		// save viscosity  array
	ierr = VecView(rho,   					view_out); CHKERRQ(ierr); 		// save density    array
	ierr = VecView(n,   					view_out); CHKERRQ(ierr); 		// save powerlaw exponent
	ierr = VecView(G,   					view_out); CHKERRQ(ierr); 		// save elastic shear module    array
	ierr = VecView(C,   					view_out); CHKERRQ(ierr); 		// save cohesion
	ierr = VecView(phi,   					view_out); CHKERRQ(ierr); 		// save friction angle
	ierr = VecView(k,   					view_out); CHKERRQ(ierr); 		// save thermal conductivity array
	ierr = VecView(Cp,   					view_out); CHKERRQ(ierr); 		// save heat capacity
	ierr = VecView(Q,   					view_out); CHKERRQ(ierr); 		// save heat production
	ierr = VecView(alpha, 					view_out); CHKERRQ(ierr); 		// save thermal expansivity
	ierr = VecView(FK, 					view_out); CHKERRQ(ierr); 		// save Frank-Kamenetskii
	ierr = VecView(intpx, 					view_out); CHKERRQ(ierr); 		// save x-coordinates array
	ierr = VecView(intpy, 					view_out); CHKERRQ(ierr); 		// save y-coordinates array
	ierr = VecView(intpz, 					view_out); CHKERRQ(ierr); 		// save z-coordinates array
	ierr = VecView(characteristic,			view_out); CHKERRQ(ierr); 		// save characteristic values
#ifdef TEMPERATURE
	ierr = VecView(T,						view_out); CHKERRQ(ierr); 		// save temperature array
#endif
	ierr = VecView(Phases,					view_out); CHKERRQ(ierr); 		// save phase distribution array
	ierr = VecView(P,						view_out); CHKERRQ(ierr); 		// save pressure array
	ierr = VecView(DevStress,				view_out); CHKERRQ(ierr); 		// save 2nd invariant of stress tensor
	ierr = VecView(DevStrainrate,			view_out); CHKERRQ(ierr); 		// save 2nd invariant of strainrate tensor
	ierr = VecView(Txx,					view_out); CHKERRQ(ierr); 		// save Tau_xx array
	ierr = VecView(Tyy,					view_out); CHKERRQ(ierr); 		// save Tau_yy array
	ierr = VecView(Tzz,					view_out); CHKERRQ(ierr); 		// save Tau_zz array
	ierr = VecView(Txz,					view_out); CHKERRQ(ierr); 		// save Tau_xz array
	ierr = VecView(Tyz,					view_out); CHKERRQ(ierr); 		// save Tau_yz array
	ierr = VecView(Txz,					view_out); CHKERRQ(ierr); 		// save Tau_xz array
	ierr = VecView(Exx,					view_out); CHKERRQ(ierr); 		// save Eps_xx array
	ierr = VecView(Eyy,					view_out); CHKERRQ(ierr); 		// save Eps_yy array
	ierr = VecView(Ezz,					view_out); CHKERRQ(ierr); 		// save Eps_zz array
	ierr = VecView(Exz,					view_out); CHKERRQ(ierr); 		// save Eps_xz array
	ierr = VecView(Eyz,					view_out); CHKERRQ(ierr); 		// save Eps_yz array
	ierr = VecView(Exz,					view_out); CHKERRQ(ierr); 		// save Eps_xz array
	ierr = MatView(TimeDependentDataArray, view_out); CHKERRQ(ierr);		// Store time-dependent data
	ierr = VecView(NumParticles,			view_out); CHKERRQ(ierr); 		// # of particles per integration point

	ierr = VecView(Strain,					view_out); CHKERRQ(ierr); 		// # of particles per integration point
	ierr = VecView(PlasticStrain,			view_out); CHKERRQ(ierr); 		// # of particles per integration point

	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"# Saved MATLAB output in file %s \n",SaveFileName); CHKERRQ(ierr);



	//DMDAVecRestoreArray(cda,coord,&coords); CHKERRQ(ierr);



	ierr = VecDestroy(&information); CHKERRQ(ierr);		ierr = VecDestroy(&characteristic); CHKERRQ(ierr);
	ierr = VecDestroy(&coord_x); CHKERRQ(ierr);				ierr = VecDestroy(&coord_y); CHKERRQ(ierr);								ierr = VecDestroy(&coord_z); CHKERRQ(ierr);
	ierr = VecDestroy(&Vx); CHKERRQ(ierr);							ierr = VecDestroy(&Vy); CHKERRQ(ierr);											ierr = VecDestroy(&Vz); CHKERRQ(ierr);
	ierr = VecDestroy(&mu); CHKERRQ(ierr);							ierr = VecDestroy(&rho); CHKERRQ(ierr);										ierr = VecDestroy(&intpx); CHKERRQ(ierr);
	ierr = VecDestroy(&intpy); CHKERRQ(ierr);					ierr = VecDestroy(&intpz); CHKERRQ(ierr);									ierr = VecDestroy(&G); CHKERRQ(ierr);
	ierr = VecDestroy(&C); CHKERRQ(ierr);							ierr = VecDestroy(&k); CHKERRQ(ierr);											ierr = VecDestroy(&Cp); CHKERRQ(ierr);
	ierr = VecDestroy(&Q); CHKERRQ(ierr);							ierr = VecDestroy(&alpha); CHKERRQ(ierr);									ierr = VecDestroy(&FK); CHKERRQ(ierr);
	ierr = VecDestroy(&n); CHKERRQ(ierr);							ierr = VecDestroy(&phi); CHKERRQ(ierr);									/*ierr = VecDestroy(local_Vel); CHKERRQ(ierr);*/
	ierr = VecDestroy(&Phases); CHKERRQ(ierr);					ierr = VecDestroy(&P); CHKERRQ(ierr);											ierr = VecDestroy(&DevStress); CHKERRQ(ierr);
	ierr = VecDestroy(&DevStrainrate); CHKERRQ(ierr);	ierr = VecDestroy(&Txx); CHKERRQ(ierr);										ierr = VecDestroy(&Tyy); CHKERRQ(ierr);
	ierr = VecDestroy(&Tzz); CHKERRQ(ierr);						ierr = VecDestroy(&Txy); CHKERRQ(ierr);										ierr = VecDestroy(&Tyz); CHKERRQ(ierr);
	ierr = VecDestroy(&Txz); CHKERRQ(ierr);						ierr = VecDestroy(&Exx); CHKERRQ(ierr);										ierr = VecDestroy(&Eyy); CHKERRQ(ierr);
	ierr = VecDestroy(&Ezz); CHKERRQ(ierr);						ierr = VecDestroy(&Exy); CHKERRQ(ierr);										ierr = VecDestroy(&Eyz); CHKERRQ(ierr);
	ierr = VecDestroy(&Exz); CHKERRQ(ierr);						ierr = MatDestroy(&TimeDependentDataArray); CHKERRQ(ierr); ierr = VecDestroy(&NumParticles); CHKERRQ(ierr);
	ierr = VecDestroy(&Strain); CHKERRQ(ierr);					ierr = VecDestroy(&PlasticStrain); CHKERRQ(ierr);
#ifdef TEMPERATURE
	ierr = VecDestroy(&T); CHKERRQ(ierr);				/*ierr = VecDestroy(&Temp_local); CHKERRQ(ierr);*/
#endif

	ierr = PetscFree(coords_x); CHKERRQ(ierr);    		ierr = PetscFree(coords_y); CHKERRQ(ierr);    	ierr = PetscFree(coords_z); CHKERRQ(ierr);
	ierr = PetscFree(Vx_array); CHKERRQ(ierr);  			ierr = PetscFree(Vy_array); CHKERRQ(ierr);   		ierr = PetscFree(Vz_array); CHKERRQ(ierr);
	ierr = PetscFree(mu_array); CHKERRQ(ierr);				ierr = PetscFree(rho_array); CHKERRQ(ierr);			ierr = PetscFree(intpx_array); CHKERRQ(ierr);
	ierr = PetscFree(intpy_array); CHKERRQ(ierr); 		ierr = PetscFree(intpz_array); CHKERRQ(ierr);
	ierr = PetscFree(G_array); CHKERRQ(ierr);					ierr = PetscFree(C_array); CHKERRQ(ierr);						ierr = PetscFree(k_array); CHKERRQ(ierr);
	ierr = PetscFree(Cp_array); CHKERRQ(ierr);				ierr = PetscFree(Q_array); CHKERRQ(ierr);						ierr = PetscFree(Alpha_array); CHKERRQ(ierr);
	ierr = PetscFree(FK_array); CHKERRQ(ierr);				ierr = PetscFree(n_array); CHKERRQ(ierr);						ierr = PetscFree(phi_array); CHKERRQ(ierr);
	ierr = PetscFree(Phase_array); CHKERRQ(ierr); 		ierr = PetscFree(Pressure_array); CHKERRQ(ierr);		ierr = PetscFree(DevStress_array); CHKERRQ(ierr);
	ierr = PetscFree(DevStrainrate_array); CHKERRQ(ierr); ierr = PetscFree(Txx_array); CHKERRQ(ierr); 		ierr = PetscFree(Tyy_array); CHKERRQ(ierr);
	ierr = PetscFree(Tzz_array); CHKERRQ(ierr);						ierr = PetscFree(Txy_array); CHKERRQ(ierr);			ierr = PetscFree(Txz_array); CHKERRQ(ierr);
	ierr = PetscFree(Tyz_array); CHKERRQ(ierr);						ierr = PetscFree(Exx_array); CHKERRQ(ierr); 		ierr = PetscFree(Eyy_array); CHKERRQ(ierr);
	ierr = PetscFree(Ezz_array); CHKERRQ(ierr);						ierr = PetscFree(Exy_array); CHKERRQ(ierr);			ierr = PetscFree(Exz_array); CHKERRQ(ierr);
	ierr = PetscFree(Eyz_array); CHKERRQ(ierr);						ierr = PetscFree(NumParticles_array); CHKERRQ(ierr);
	ierr = PetscFree(Strain_array); CHKERRQ(ierr);				ierr = PetscFree(PlasticStrain_array); CHKERRQ(ierr);
#ifdef TEMPERATURE
	ierr = PetscFree(Temp_array); CHKERRQ(ierr);
#endif

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
