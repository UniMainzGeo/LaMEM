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

Quadrature.c, contains the following subroutines:

IntegrationPoints				-	Compute the integration points
IntPointProperties				-	Compute integration point properties

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "Quadrature.h"
#include "Material.h"
#include "Elements.h"
#include "StokesOperators.h"
#include "LaMEM_Temperature.h"


/*==========================================================================================================*/
/* Computes integration points and weights for Q1 and Q2 elements*/
#undef __FUNCT__
#define __FUNCT__ "IntegrationPoints"
PetscErrorCode IntegrationPoints( const PetscInt nnel, const PetscInt nintp_1D, double IntPoint[3][MAX_ngp_vel], double IntWeight[] )
{
	double points[3], weights[3];
	PetscInt   intx,inty,intz, ngp_vel_i;

	if ((nnel==8) && (nintp_1D>1)){
		/* A Q1P0 or a Q1Q1 element */
		points[0]  = -0.577350269189626;  points[1]  =  0.577350269189626;  points[2]  =  0.0;
		weights[0] = 1.0;  weights[1] = 1.0;   weights[2] = 0.0;
	}
	else if ((nnel==8) && (nintp_1D==1)){
		/* FDSTAG formulation; here we only compute props @ the center */
		points[0]  = 0.0;  points[1]  = 0.0;  points[2]  = 0.0;
		weights[0] = 2.0;  weights[1] = 0.0;  weights[2] = 0.0;

	}
	else if (nnel==27){
		points[0]  = -0.774596669241483;  points[1]  =  0.0;  points[2]  =  0.774596669241483;
		weights[0] = 0.555555555555556;   weights[1] =  0.888888888888889;   weights[2] = 0.555555555555556;
	}

	ngp_vel_i = 0;
	for (intx=0; intx<nintp_1D; intx++){
		for (inty=0; inty<nintp_1D; inty++){
			for (intz=0; intz<nintp_1D; intz++){
				IntPoint[0][ngp_vel_i]	=	points[intx];
				IntPoint[1][ngp_vel_i]	=	points[inty];
				IntPoint[2][ngp_vel_i]	=	points[intz];
				IntWeight[ngp_vel_i]	=	weights[intx]*weights[inty]*weights[intz];
				ngp_vel_i = ngp_vel_i+1;
			}
		}
	}
	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Computes properties such as pressure, strainrate, deviatoric stress @ integration points
 * strictly speaking, this routine is for FE calculations, as this is the only place where we really use
 * properties @ integration points. In FDSTAG calculations, the routine is useful for visualization.
 *
 * */
#undef __FUNCT__
#define __FUNCT__ "IntPointProperties"
PetscErrorCode IntPointProperties( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Pressure, DM da_temp,
		Vec Temp, UserContext *user, PetscScalar dt, PetscInt ComputeFull )
{
	PetscErrorCode 				ierr;
	PetscInt       				i,j,k,ipres, iel_x, iel_y, iel_z,ii, nel_x, nel_y, nel_z;
	PetscInt		 			xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt		 			intp, num_elem;
	DMDACoor3d		 			***coords, coord_elem[MAX_nnel];
	DM			 				cda;
	Vec			 				gc;
	Field			 			***velocity;
	PetscScalar	 				V_element[MAX_edof], P_element[MAX_npres], Point[3],  DevStress[6], sizeintp, TotalMemoryOfIntp, MinMemoryOfIntp, MaxMemoryOfIntp;
	PetscScalar	 				mu, mu_viscous, IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel];
	PressureElem	 			***PressureNew;
	Vec							local_Vel, Temp_local, Pressure_local;
	//MaterialsElement			***materials;
	PointWiseInformation   		PointInformation;
	PetscInt 					ngp_vel, nnel, nintp_1D, npres;
	DAVPElementType 			element_type;
	PetscScalar 				***materials_array;
	MaterialsElementDynamic 	material_data;
	PetscScalar 				***PressureNew_array;
	P_array 					***PressureContinuous_array;
	PressureElemDynamic 		PressureNew_data;
	PetscScalar 				***temperature, T_element[MAX_edof_temp];
	PetscScalar 				MaximumT2nd_loc, MinimumT2nd_loc, MaximumE2nd_loc, MinimumE2nd_loc, MaximumP_loc, MinimumP_loc, Maximum_mu_loc, Minimum_mu_loc;
	PetscScalar 				MaximumT2nd, MinimumT2nd, MaximumE2nd, MinimumE2nd, MaximumP, MinimumP, Maximum_mu, Minimum_mu;
	PetscBool					flg;


	PetscFunctionBegin;


	element_type = C->type;
	ngp_vel  = C->ngp_vel;
	nnel     = C->nnel;
	npres    = C->npres;
	nintp_1D = C->nintp_1D;

	//coordinates
	//if (da != PETSC_NULL){
		ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);
	//}

	/* Copy pressure to local array, including ghostpoints */
	if (da_pres != PETSC_NULL){
		ierr = DMCreateLocalVector(da_pres,&Pressure_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da_pres, Pressure, INSERT_VALUES, Pressure_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da_pres,   Pressure, INSERT_VALUES, Pressure_local); CHKERRQ(ierr);

		ierr = DMDAVecGetArray(da_pres, 			Pressure_local, 		&PressureNew); CHKERRQ(ierr);
		if (element_type==DAVP_Q1Q1){
			// This needs
			ierr = DMDAVecGetArray(da_pres, 		Pressure_local, 		&PressureContinuous_array); CHKERRQ(ierr);
		}
		//else if (element_type==DAVP_FDSTAG){
		//	SETERRQ( PETSC_ERR_SUP, "Routine does not yet work for FDSTAG!" );
		//	ierr = MPI_Abort(PETSC_COMM_WORLD,1);  CHKERRQ(ierr);
		//}
		else {
			ierr = DMDAVecGetArray(da_pres, 		Pressure_local, 		&PressureNew_array); CHKERRQ(ierr);
		}
	}

	/* Copy global velocity solution to local processor, including ghostpoints */
	if (Velocity != PETSC_NULL){
		ierr = DMGetLocalVector(da,&local_Vel); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, local_Vel,	&velocity);		CHKERRQ(ierr);
	}

	/* Copy temperature to local processor */
	if (da_temp != PETSC_NULL){
		ierr = DMCreateLocalVector(da_temp,&Temp_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_temp, Temp_local,	&temperature); CHKERRQ(ierr);
	}

	// get material properties
	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);


	num_elem = 0;
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetInfo(user->DA_Processors,	0,&nel_x	,&nel_y,	&nel_z,		0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr); 	// # of elements in all directions

	/* Initialize variables which*/
	MaximumT2nd_loc 	=  -1e100;
	MinimumT2nd_loc 	=	1e100;
	MaximumE2nd_loc 	=  -1e100;
	MinimumE2nd_loc 	=	1e100;
	MaximumP_loc 		=  -1e100;
	MinimumP_loc 		=	1e100;
	Maximum_mu_loc 		=  -1e100;
	Minimum_mu_loc 		=	1e100;

	// Loop over all elements
	i=0; j=0; k=0;
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){

				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xsp,iel_y-ysp,iel_z-zsp, C->ngp_vel, materials_array );


				if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) || (element_type==DAVP_FDSTAG)){
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}
				else {
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type" );
				}

				/*
				i = element_counter_scale*iel_x;
				j = element_counter_scale*iel_y;
				k = element_counter_scale*iel_z;
				*/
				/*------------------------------------------------------
				* Form the local stiffness matrix
				* -----------------------------------------------------*/
				/* Extract coordinates of the local element in correct order */
				GetElementCoords(coord_elem, coords, i,j,k,1);
				CorrectElementCoordsForPeriodicity( coord_elem, i, j, user, 1);


				/* Stuff that can be precomputed */
				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				/* Extract nodal velocities of the current element */
				if (Velocity != PETSC_NULL){
					GetVelocityElement(velocity, V_element, i,j,k);
				}
				else {
					for (ii=0; ii<MAX_edof; ii++){
						V_element[ii] = 0.0;
					}
				}

				/* Extract Pressure shape function values of the current element */
				if (da_pres != PETSC_NULL){
					if ( element_type==DAVP_Q1Q1 ){
						/* Continuous pressure */
						GetPressureElement(PressureContinuous_array, P_element, iel_x,iel_y,iel_z);
					}
					else if ( element_type==DAVP_FDSTAG){
						// for FDSTAG, P is defined already at the center of the control volume
						P_element[0] = PressureNew_array[k][j][i];
					}
					else {
						/* Discontinuous pressure */
						LaMEMSetPressueElemDataMemoryFromArray( &PressureNew_data, iel_x,iel_y,iel_z, npres, PressureNew_array );
						for (ipres=0; ipres<npres; ipres++){
							P_element[ipres]  = PressureNew_data.P[ipres];
						}
					}
				}
				else {
					for (ii=0; ii<MAX_npres; ii++){
						P_element[ii] = 0.0;
					}
				}

				/* Extract Temperature nodal values of the current element */
				if (da_temp != PETSC_NULL){
					GetTemperatureElement(temperature    , T_element, i,j,k);
				}
				else {
					for (ii=0; ii<MAX_edof_temp; ii++){
						T_element[ii] = 0.0;
					}
				}

				/* Loop over integration points */
				for (intp=0; intp<ngp_vel; intp++){
					mu 	  		= 	material_data.Viscosity[intp];				 // 	Effective viscosity
					mu_viscous 	= 	material_data.TrueViscosity[intp];			 // 	Viscosity
//					G 			=  	material_data.ElasticShearModule[intp];       //  Elastic shear module
					Point[0]	= 	IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];

					for (i=0; i<6; i++){
						DevStress[i] = material_data.DevStress[i][intp];
					}

					/* Compute properties */
					ierr = ComputePointWiseProperties( C, &PointInformation, Point, V_element, P_element, T_element, T_element,
							DevStress, coord_elem, mu, mu_viscous, ComputeFull ); CHKERRQ(ierr);


					// Store information @ integration point

					/* Data that should be updated during iterations */
					material_data.Pressure[intp]   						= PointInformation.Pressure;
					material_data.SecondInvariantDevStress[intp]   		= PointInformation.SecondInvariantDeviatoricStress;
					material_data.SecondInvariantDevStrainrate[intp]   	= PointInformation.SecondInvariantDeviatoricStrainrate;
					material_data.Temperature[intp]   				= 	PointInformation.Temperature;

					if (ComputeFull==1){
						/* Data that is required only once per timestep 	*/
						material_data.Coord[0][intp]   					= 	PointInformation.x;
						material_data.Coord[1][intp]   					= 	PointInformation.y;
						material_data.Coord[2][intp]   					= 	PointInformation.z;

						material_data.DevStress[0][intp]   				=	PointInformation.DeviatoricStress[0];
						material_data.DevStress[1][intp]   				= 	PointInformation.DeviatoricStress[1];
						material_data.DevStress[2][intp]   				= 	PointInformation.DeviatoricStress[2];
						material_data.DevStress[3][intp]   				= 	PointInformation.DeviatoricStress[3];
						material_data.DevStress[4][intp]   				= 	PointInformation.DeviatoricStress[4];
						material_data.DevStress[5][intp]   				= 	PointInformation.DeviatoricStress[5];

						material_data.DevStrainrate[0][intp]   			= 	PointInformation.DeviatoricStrainRate[0];
						material_data.DevStrainrate[1][intp]   			= 	PointInformation.DeviatoricStrainRate[1];
						material_data.DevStrainrate[2][intp]   			= 	PointInformation.DeviatoricStrainRate[2];
						material_data.DevStrainrate[3][intp]   			= 	PointInformation.DeviatoricStrainRate[3];
						material_data.DevStrainrate[4][intp]   			= 	PointInformation.DeviatoricStrainRate[4];
						material_data.DevStrainrate[5][intp]   			= 	PointInformation.DeviatoricStrainRate[5];

						material_data.ShearHeat[intp] 					= 	PointInformation.ShearHeat;

						material_data.Strain[intp]						=	material_data.Strain[intp] + dt*PointInformation.SecondInvariantDeviatoricStrainrate;
						material_data.PlasticStrain[intp]				=	material_data.PlasticStrain[intp] + dt*PointInformation.PlasticStrainrate2ndInvariant;

						material_data.ElementVolume[intp]				=	PointInformation.ElementVolume;

						//material_data.Strain[intp]						=	PointInformation.PlasticStrainrate2ndInvariant;

					}

					// Keep track of max. and min. values
					if (PointInformation.SecondInvariantDeviatoricStrainrate>MaximumE2nd_loc){
						MaximumE2nd_loc = PointInformation.SecondInvariantDeviatoricStrainrate;
					}
					if (PointInformation.SecondInvariantDeviatoricStrainrate<MinimumE2nd_loc){
						MinimumE2nd_loc = PointInformation.SecondInvariantDeviatoricStrainrate;
					}
					if (PointInformation.SecondInvariantDeviatoricStress>MaximumT2nd_loc){
						MaximumT2nd_loc = PointInformation.SecondInvariantDeviatoricStress;
					}
					if (PointInformation.SecondInvariantDeviatoricStress<MinimumT2nd_loc){
						MinimumT2nd_loc = PointInformation.SecondInvariantDeviatoricStress;
					}
					if (PointInformation.Pressure>MaximumP_loc){
						MaximumP_loc 	= PointInformation.Pressure;
					}
					if (PointInformation.Pressure<MinimumP_loc){
						MinimumP_loc 	= PointInformation.Pressure;
					}
					if (mu>Maximum_mu_loc){
						Maximum_mu_loc 	= mu;
					}
					if (mu<Minimum_mu_loc){
						Minimum_mu_loc 	= mu;
					}


				}
				/*------------------------------------------------------
				* End of forming the local stiffness matrix
				* ----------------------------------------------------- */


				num_elem = num_elem+1;
			}
		}
	}

	/* restore pressure */
	if (da_pres != PETSC_NULL){
		ierr = DMDAVecRestoreArray(da_pres,			Pressure_local,			&PressureNew); CHKERRQ(ierr);
		if (element_type==DAVP_Q1Q1){
			ierr = DMDAVecRestoreArray(da_pres, 		Pressure_local, 			&PressureContinuous_array); CHKERRQ(ierr);
		}
		else {
			ierr = DMDAVecRestoreArray(da_pres, 		Pressure_local, 			&PressureNew_array); CHKERRQ(ierr);
		}
		ierr = VecDestroy(&Pressure_local);
	}


	ierr = DMDAVecRestoreArray(user->DA_Materials,   user->Materials, &materials_array); CHKERRQ(ierr);

	if (Velocity != PETSC_NULL){
		ierr = DMDAVecRestoreArray(da,local_Vel,&velocity);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(da,&local_Vel);  // YES !! This should match calls to DMGetLocalVector() //
	}

	ierr = DMDAVecRestoreArray(cda,gc,&coords);	CHKERRQ(ierr);
	//ierr = DMDestroy(&cda); CHKERRQ(ierr);
	//ierr = VecDestroy(&gc); CHKERRQ(ierr);

	/* Restore temperature */
	if (da_temp != PETSC_NULL){
		ierr = DMDAVecRestoreArray(da_temp, Temp_local,		&temperature	); CHKERRQ(ierr);
		ierr = VecDestroy(&Temp_local); CHKERRQ(ierr);
	}


	/* Collect global max. && min. values */
	ierr = MPI_Allreduce(&MaximumE2nd_loc,	&MaximumE2nd,	1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MinimumE2nd_loc,	&MinimumE2nd,	1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MaximumT2nd_loc,	&MaximumT2nd,	1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MinimumT2nd_loc,	&MinimumT2nd,	1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MaximumP_loc,	&MaximumP,		1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);     CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MinimumP_loc,	&MinimumP,		1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);     CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Maximum_mu_loc,	&Maximum_mu,	1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Minimum_mu_loc,	&Minimum_mu,	1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);

	/* Display result if required */
	PetscOptionsHasName(PETSC_NULL,"-DisplayStresses",&flg);
	if (flg){
		// yes, print some stuff
		PetscPrintf(PETSC_COMM_WORLD," T2nd:  Max=%8.8e    Min= %8.8e \n",MaximumT2nd,MinimumT2nd);
		PetscPrintf(PETSC_COMM_WORLD," E2nd:  Max=%8.8e    Min= %8.8e \n",MaximumE2nd,MinimumE2nd);
		PetscPrintf(PETSC_COMM_WORLD," P   :  Max=%8.8e    Min= %8.8e \n",MaximumP,   MinimumP   );
		PetscPrintf(PETSC_COMM_WORLD," Mu  :  Max=%8.8e    Min= %8.8e \n",Maximum_mu, Minimum_mu );

	}
	sizeintp = (double)((PetscInt)sizeof(MaterialsElement)/MAX_ngp_vel*ngp_vel*xmp*ymp*zmp )*1e-6;
	ierr = MPI_Allreduce(&sizeintp, &TotalMemoryOfIntp, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The total memory used by integration points is %e MB\n",TotalMemoryOfIntp);
	ierr = MPI_Allreduce(&sizeintp, &MaxMemoryOfIntp, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The max memory used by integration points is %e MB\n",MaxMemoryOfIntp);
	ierr = MPI_Allreduce(&sizeintp, &MinMemoryOfIntp, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The min memory used by integration points is %e MB\n",MinMemoryOfIntp);

	/* Clean up */
	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
