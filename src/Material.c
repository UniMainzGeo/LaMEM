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

#START_DOC#

contains the following subroutines:

\lamemfunction{\verb-SetMaterialProperties-}
Hardcoded setting of material properties

\lamemfunction{\verb-ComputeEffectiveViscosity-}
Computes Frank Kamenetskii effective viscosity

\lamemfunction{\verb-ComputePointWiseProperties-}
Compute properties such as velocity, strainrates, pressures etc.  at a local point, with given element coordinates




#END_DOC#

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
 */


#include "LaMEM.h"
#include "Material.h"
#include "Utils.h"
#include "StokesOperators.h"
#include "Elements.h"
#include "Quadrature.h"


# if 0
#undef __FUNCT__
#define __FUNCT__ "LaMEMSetMaterialDataMemoryFromArray_New"
PetscErrorCode LaMEMSetMaterialDataMemoryFromArray_New(
		MaterialsElementDynamic *data,
		const PetscInt i, const PetscInt j, const PetscInt k,
		const PetscInt ngp,
		PetscScalar ***list[NumMaterialPropsElem] )
{
	PetscScalar *block_ijk;
	PetscInt c;

	block_ijk = &list[k][j][i*ngp][0];


	data->Phases    = &block_ijk;		//c = c + ngp;

#if 0
	data->Phases    = &block_ijk[c];		c = c + ngp;
	data->Coord[0]  = &block_ijk[c];		c = c + ngp;
	data->Coord[1]  = &block_ijk[c];		c = c + ngp;
	data->Coord[2]  = &block_ijk[c];		c = c + ngp;
	data->Pressure  = &block_ijk[c];		c = c + ngp;
	data->Viscosity = &block_ijk[c];		c = c + ngp;
	
	data->NumParticles       = &block_ijk[c];		c = c + ngp;

	data->ElasticShearModule = &block_ijk[c];		c = c + ngp;


	data->Density        = &block_ijk[c];		c = c + ngp;


	data->DevStress[0] = &block_ijk[c];		c = c + ngp;
	data->DevStress[1] = &block_ijk[c];		c = c + ngp;
	data->DevStress[2] = &block_ijk[c];		c = c + ngp;
	data->DevStress[3] = &block_ijk[c];		c = c + ngp;
	data->DevStress[4] = &block_ijk[c];		c = c + ngp;
	data->DevStress[5] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[0] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[1] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[2] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[3] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[4] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[5] = &block_ijk[c];		c = c + ngp;

	data->SecondInvariantDevStress     = &block_ijk[c];		c = c + ngp;
	data->SecondInvariantDevStrainrate = &block_ijk[c];		c = c + ngp;

	data->Temperature 				    = &block_ijk[c];		c = c + ngp;
	data->Strain 				   		= &block_ijk[c];		c = c + ngp;
	data->PlasticStrain 				= &block_ijk[c];		c = c + ngp;
	data->ShearHeat 				    = &block_ijk[c];		c = c + ngp;
	data->PreExpFactor					= &block_ijk[c];		c = c + ngp;
	data->ActivationEnergy				= &block_ijk[c];		c = c + ngp;

	data->PlasticViscosity 				= &block_ijk[c];		c = c + ngp;
	data->Plastic 						= &block_ijk[c];		c = c + ngp;

	data->TrueViscosity 				= &block_ijk[c];		c = c + ngp;

	data->ElementVolume 				= &block_ijk[c];		c = c + ngp;

	if( c != NumMaterialPropsElem*ngp ) {
		PetscPrintf( PETSC_COMM_WORLD, "*** It appears you have not accessed all element material properties correctly *** \n");
	}
#endif

	PetscFunctionReturn(0);
}
#endif

#undef __FUNCT__
#define __FUNCT__ "LaMEMSetMaterialDataMemoryFromArray"
PetscErrorCode LaMEMSetMaterialDataMemoryFromArray(
		MaterialsElementDynamic *data,
		const PetscInt i, const PetscInt j, const PetscInt k,
		const PetscInt ngp,
		PetscScalar ***list )
{
	PetscScalar *block_ijk;
	PetscInt c;

	block_ijk = &list[k][j][i*NumMaterialPropsElem*ngp];


	c = 0;
	data->Phases    = &block_ijk[c];		c = c + ngp;
	data->Coord[0]  = &block_ijk[c];		c = c + ngp;
	data->Coord[1]  = &block_ijk[c];		c = c + ngp;
	data->Coord[2]  = &block_ijk[c];		c = c + ngp;
	data->Pressure  = &block_ijk[c];		c = c + ngp;
	data->Viscosity = &block_ijk[c];		c = c + ngp;

	data->NumParticles       = &block_ijk[c];		c = c + ngp;

	data->ElasticShearModule = &block_ijk[c];		c = c + ngp;


	data->Density        = &block_ijk[c];		c = c + ngp;


	data->DevStress[0] = &block_ijk[c];		c = c + ngp;
	data->DevStress[1] = &block_ijk[c];		c = c + ngp;
	data->DevStress[2] = &block_ijk[c];		c = c + ngp;
	data->DevStress[3] = &block_ijk[c];		c = c + ngp;
	data->DevStress[4] = &block_ijk[c];		c = c + ngp;
	data->DevStress[5] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[0] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[1] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[2] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[3] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[4] = &block_ijk[c];		c = c + ngp;
	data->DevStrainrate[5] = &block_ijk[c];		c = c + ngp;

	data->SecondInvariantDevStress     = &block_ijk[c];		c = c + ngp;
	data->SecondInvariantDevStrainrate = &block_ijk[c];		c = c + ngp;

	data->Temperature 				    = &block_ijk[c];		c = c + ngp;
	data->Strain 				   		= &block_ijk[c];		c = c + ngp;
	data->PlasticStrain 				= &block_ijk[c];		c = c + ngp;
	data->ShearHeat 				    = &block_ijk[c];		c = c + ngp;
	data->PreExpFactor					= &block_ijk[c];		c = c + ngp;
	data->ActivationEnergy				= &block_ijk[c];		c = c + ngp;

	data->PlasticViscosity 				= &block_ijk[c];		c = c + ngp;
	data->Plastic 						= &block_ijk[c];		c = c + ngp;

	data->TrueViscosity 				= &block_ijk[c];		c = c + ngp;

	data->ElementVolume 				= &block_ijk[c];		c = c + ngp;


	if( c != NumMaterialPropsElem*ngp ) {
		PetscPrintf( PETSC_COMM_WORLD, "*** It appears you have not accessed all element material properties correctly *** \n");
	}

	PetscFunctionReturn(0);
}


/*=======================================================================================================*/
/* Set material properties for each of the elements - currently done by hard-coding;
 * in future computed from particles */
#undef __FUNCT__
#define __FUNCT__ "SetMaterialProperties"
PetscErrorCode SetMaterialProperties( LaMEMVelPressureDA C, UserContext *user, DM DA_MG)
{
	PetscErrorCode 			ierr;
	PetscInt 				xs,ys,zs,xm,ym,zm,iel_x,iel_y,iel_z, intp,i,j,k, phase, Even;
	DMDACoor3d				***coords, coord_elem[MAX_nnel], CoordIntp;
	Vec						gc;
	DM						cda;
	PetscScalar	 			IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel], Point[3];
	PetscScalar				ShapeVel[MAX_nnel], **dhdsVel;
	DAVPElementType 		element_type;
	PetscInt				nnel, nintp_1D;
	PetscScalar 			***materials_array;
	MaterialsElementDynamic material_data;

	element_type = C->type;
	nnel         = C->nnel;
	nintp_1D     = C->nintp_1D;

	LaMEMCreate2dArray( 3, C->nnel, &dhdsVel, PETSC_NULL );

	ierr = DMGetCoordinateDM(DA_MG, &cda); CHKERRQ(ierr); //coordinates
	ierr = DMGetCoordinatesLocal(DA_MG,	&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, gc, &coords); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	/* Loop over elements */
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (iel_z=zs; iel_z<zs+zm; iel_z++){
		for (iel_y=ys; iel_y<ys+ym; iel_y++){
			for (iel_x=xs; iel_x<xs+xm; iel_x++){

				/* load data */
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xs,iel_y-ys,iel_z-zs, C->ngp_vel, materials_array );

				if( (element_type==DAVP_Q1P0) | (element_type==DAVP_Q1Q1) | (element_type==DAVP_FDSTAG) ) {
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
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented" );
				}

				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				// Extract coordinates of the local element in correct order
				GetElementCoords(coord_elem, coords, i,j,k, 1);

				// Make a loop over each integration point
				for (intp=0; intp<C->ngp_vel; intp++){
					// retrieve the coordinates of the integration point
					Point[0]= IntPoint[0][intp];
					Point[1]= IntPoint[1][intp];
					Point[2]= IntPoint[2][intp];
					ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);
					ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordIntp);


					material_data.Phases[intp] 				= 0;
					material_data.Viscosity[intp]  			= user->PhaseProperties.mu[0];
					material_data.Density[intp]  			= user->PhaseProperties.rho[0];
					material_data.ElasticShearModule[intp] 	= user->PhaseProperties.ElasticShearModule[0];
					material_data.Strain[intp]  			= 0.0;
					material_data.PlasticStrain[intp]  		= 0.0;
					material_data.Pressure[intp]  			= 0.0;
					material_data.Temperature[intp]  		= 0.0;
					material_data.ShearHeat[intp]  			= 0.0;
					material_data.PreExpFactor[intp]		= user->PhaseProperties.A[0];
					material_data.ActivationEnergy[intp]	= user->PhaseProperties.E[0];

					material_data.Coord[0][intp]   = CoordIntp.x;
					material_data.Coord[1][intp]   = CoordIntp.y;
					material_data.Coord[2][intp]   = CoordIntp.z;

					/* Set Phases depending on the case of Model Setup */
					if 		(user->Setup.Model==0){				// diapir setup
						if ((k+C->nnode_el_1D)>user->Setup.ind_Hi_diapir){
							material_data.Phases[intp] 			= 1;
						}
					}
					else if ((user->Setup.Model==1) || (user->Setup.Model==5)){ 			// single layer folding setup, growth rate test
						if ( ((k+C->nnode_el_1D)<=user->Setup.ind_fold_top) && ((k)>=(user->Setup.ind_fold_bot)) ){
							material_data.Phases[intp] 				= 1;
						}
					}
					else if (user->Setup.Model==4){  // multilayer detachment folding setup
						if ((k+C->nnode_el_1D)>user->Setup.ind_Hi_diapir){
							phase = 1;

							LaMEMMod( iel_z, 2, &Even);
							if (Even==0){
								phase = 1;
							}
							else {
								phase = 2;
							}
							material_data.Phases[intp] 			= phase;
						}
					}
					else if ((user->Setup.Model==-1) | (user->Setup.Model==3) | (user->Setup.Model==6) | (user->Setup.Model==2) | (user->Setup.Model==8) | (user->Setup.Model==9) | (user->Setup.Model==10)){  // from particles

					}
					else if (user->Setup.Model==7){ 			// inhomogeneity in center of cube (by Sarah)

						if ( ((k+C->nnode_el_1D)<=((user->nnode_z-1)/2+1)) && ((k)>=((user->nnode_z-1)/2-1)) &&
								((j+C->nnode_el_1D)<=((user->nnode_y-1)/2+1)) && ((j)>=((user->nnode_y-1)/2-1)) &&
								((i+C->nnode_el_1D)<=((user->nnode_x-1)/2+1)) && ((i)>=((user->nnode_x-1)/2-1))){

							material_data.Phases[intp] 				= 1;
						}

					}
					else{
						PetscPrintf(PETSC_COMM_WORLD,"  Unknown Setup! user->Setup.Model=%lld \n", (LLD)(user->Setup.Model));
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_Materials,  user->Materials, &materials_array); CHKERRQ(ierr);

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Reset stresses */
#undef __FUNCT__
#define __FUNCT__ "ResetStressesBeforeNonlinearIterations"
PetscErrorCode ResetStressesBeforeNonlinearIterations( LaMEMVelPressureDA C, UserContext *user)
{
	PetscErrorCode 			ierr;
	PetscInt 				xs,ys,zs,xm,ym,zm,iel_x,iel_y,iel_z, intp;
	PetscScalar 			***materials_array;
	MaterialsElementDynamic material_data;

	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	// Loop over elements
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (iel_z=zs; iel_z<zs+zm; iel_z++){
		for (iel_y=ys; iel_y<ys+ym; iel_y++){
			for (iel_x=xs; iel_x<xs+xm; iel_x++){
				/* load data */
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xs,iel_y-ys,iel_z-zs, C->ngp_vel, materials_array );
				// Make a loop over each integration point
				for (intp=0; intp<C->ngp_vel; intp++){
					material_data.SecondInvariantDevStress[intp] = 0;		// set dev stresses to zero here.
					//						material_data.PlasticViscosity[intp] 				= 0;
					//						material_data.SecondInvariantDevStrainrate[intp] 	= 0;

				}


			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials,  user->Materials, &materials_array); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Computes average viscosity within an element and sets it's value to all integration points if required   */
#undef __FUNCT__
#define __FUNCT__ "AverageViscosityDensityPerElement"
PetscErrorCode AverageViscosityDensityPerElement( double mu[MAX_ngp_vel], double mu_v[MAX_ngp_vel], double mu_p[MAX_ngp_vel],
		double *mu_average,  double rho[MAX_ngp_vel], UserContext *user, PetscInt ngp_vel)
{
	PetscInt 	intp;
	PetscScalar	rho_average, mu_v_average, mu_p_average;

	*mu_average  = mu_v_average = mu_p_average = 0;
	rho_average  = 0;
	if (user->MuMeanMethod==2){*mu_average=1.0; rho_average=1.0;}
	for (intp=0; intp<ngp_vel; intp++){
		if 		((user->MuMeanMethod==1) | (user->MuMeanMethod==0)){		// Arithmetic average
			*mu_average  = *mu_average  +  mu[intp]/( (double) ngp_vel);
			rho_average  = rho_average  +  rho[intp]/( (double) ngp_vel);
			mu_v_average = mu_v_average  + mu_v[intp]/( (double) ngp_vel);
			mu_p_average = mu_p_average  + mu_p[intp]/( (double) ngp_vel);

		}
		else if (user->MuMeanMethod==2){									// Geometric average
			*mu_average  = *mu_average*mu[intp];
			rho_average  = rho_average*rho[intp];
			mu_v_average = mu_v_average*mu_v[intp];
			mu_p_average = mu_p_average*mu_p[intp];

		}
		else if (user->MuMeanMethod==3){									// Harmonic average
			*mu_average  = *mu_average + 1.0/mu[intp];
			rho_average  = rho_average + 1.0/rho[intp];
			mu_v_average = mu_v_average + 1.0/mu_v[intp];
			mu_p_average = mu_p_average + 1.0/mu_p[intp];

		}
		else {
			PetscPrintf(PETSC_COMM_WORLD," Averaging method unknown!");
			MPI_Abort(PETSC_COMM_WORLD,1);
		}
	}
	if (user->MuMeanMethod==3){
		*mu_average  = (( (double) ngp_vel)/ (*mu_average) );
		rho_average  = (( (double) ngp_vel)/ (rho_average) );
		mu_v_average = (( (double) ngp_vel)/ (mu_v_average) );
		mu_p_average = (( (double) ngp_vel)/ (mu_p_average) );


	} // harmonic mean
	if (user->MuMeanMethod==2){
		*mu_average  = PetscPowScalar( (*mu_average),  1.0/( (double) ngp_vel) );

		rho_average  = PetscPowScalar( (rho_average),  1.0/( (double) ngp_vel) );

		mu_v_average = PetscPowScalar( (mu_v_average), 1.0/( (double) ngp_vel) );

		mu_p_average = PetscPowScalar( (mu_p_average), 1.0/( (double) ngp_vel) );

	} // geometric mean



	/* If required: set the viscosity of the element to the average viscosity
	 *
	 * Note that if MuMethod==0, an arithmetic average is computed per element (used
	 * in some of the preconditioners), but this is not set to the integration points.
	 */
	if (user->MuMeanMethod>0){
		for (intp=0; intp<ngp_vel; intp++){
			mu[intp] = *mu_average;
			rho[intp] = rho_average;
			mu_v[intp] = mu_v_average;
			mu_p[intp] = mu_p_average;
		}
	}



	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Average effective density from phases in the cell, both WITH and WITHOUT sticky air				   */
#undef __FUNCT__
#define __FUNCT__ "Average_EffectiveDensity_FDSTAG_Cell"
PetscErrorCode Average_EffectiveDensity_FDSTAG_Cell( const double VolumeFraction_phases[], UserContext *user,
		double T, double *EffectiveDensity_withAir, double *EffectiveDensity_withoutAir)
{
	PetscErrorCode 	ierr;
	PetscInt 		iphase, AirPhase;
	PetscScalar		VolumeFraction, rho_eff;
	PetscScalar	 	VolumeFraction_withoutAir;


	AirPhase 									=	user->ErosionParameters.StickyAirPhase;

	*EffectiveDensity_withAir 				=	0;
	*EffectiveDensity_withoutAir 				=	0;
	VolumeFraction_withoutAir 					=	0;
	for (iphase=0; iphase<user->num_phases; iphase++){
		VolumeFraction = VolumeFraction_phases[iphase];	// volume fraction of PHASE iphase


		rho_eff 					= 	user->PhaseProperties.rho[iphase];
		ierr 						= 	ComputeEffectiveDensity(iphase, user, T, &rho_eff);	CHKERRQ(ierr);
		*EffectiveDensity_withAir	+=	VolumeFraction*rho_eff;


		if (iphase != AirPhase){
			// sum over non-air phases

			VolumeFraction_withoutAir 		+= 	VolumeFraction;
			*EffectiveDensity_withoutAir	+= 	VolumeFraction*rho_eff;		// arithmetic averaging
		}

	}

	// Correct phase w/out air
	if (VolumeFraction_withoutAir>0)
	{
		*EffectiveDensity_withoutAir = *EffectiveDensity_withoutAir/VolumeFraction_withoutAir;		// required as the VolumeFraction_withoutAir will not sum to one
	}
	else{
		// only air - EffectiveDensity_withoutAir will be zero

	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Average effective thermal properties from phases in the cell, both WITH and WITHOUT sticky air			*/
#undef __FUNCT__
#define __FUNCT__ "Average_EffectiveThermalProperties_FDSTAG_Cell"
PetscErrorCode Average_EffectiveThermalProperties_FDSTAG_Cell(const double VolumeFraction_phases[], UserContext *user,
		PetscScalar *Effective_HeatCapacity_withAir, 		PetscScalar *Effective_HeatCapacity_withoutAir,
		PetscScalar *Effective_ThermalConductivity_withAir, PetscScalar *Effective_ThermalConductivity_withoutAir,
		PetscScalar *Effective_RadioactiveHeat_withAir, 	PetscScalar *Effective_RadioactiveHeat_withoutAir)
{
	PetscInt 		iphase, AirPhase;
	PetscScalar		VolumeFraction;
	PetscScalar	 	VolumeFraction_withoutAir;
	PetscScalar 	Cp_eff, k_eff, Q_eff;

	AirPhase 									=	user->ErosionParameters.StickyAirPhase;

	*Effective_HeatCapacity_withAir 			=	0;
	*Effective_HeatCapacity_withoutAir 			=	0;
	*Effective_ThermalConductivity_withAir 		=	0;
	*Effective_ThermalConductivity_withoutAir 	=	0;
	*Effective_RadioactiveHeat_withAir 			=	0;
	*Effective_RadioactiveHeat_withoutAir 		=	0;
	VolumeFraction_withoutAir 					=	0;

	for (iphase=0; iphase<user->num_phases; iphase++){
		VolumeFraction = VolumeFraction_phases[iphase];	// volume fraction of PHASE iphase

		Cp_eff 		=	user->PhaseProperties.HeatCapacity[iphase];
		k_eff 		=	user->PhaseProperties.T_Conductivity[iphase];
		Q_eff 		=	user->PhaseProperties.RadioactiveHeat[iphase];

		*Effective_HeatCapacity_withAir			+=	VolumeFraction*Cp_eff;
		*Effective_ThermalConductivity_withAir	+= 	VolumeFraction*k_eff;
		*Effective_RadioactiveHeat_withAir		+= 	VolumeFraction*Q_eff;

		if (iphase != AirPhase){
			// sum over non-air phases
			VolumeFraction_withoutAir 					+= 	VolumeFraction;
			*Effective_HeatCapacity_withoutAir			+= 	VolumeFraction*Cp_eff;		// arithmetic averaging
			*Effective_ThermalConductivity_withoutAir	+= 	VolumeFraction*k_eff;		// arithmetic averaging
			*Effective_RadioactiveHeat_withoutAir		+= 	VolumeFraction*Q_eff;		// arithmetic averaging


		}

	}

	// Correct phase w/out air
	if (VolumeFraction_withoutAir>0)
	{
		*Effective_HeatCapacity_withoutAir 			= *Effective_HeatCapacity_withoutAir/VolumeFraction_withoutAir;				// required as the VolumeFraction_withoutAir will not sum to one
		*Effective_ThermalConductivity_withoutAir 	= *Effective_ThermalConductivity_withoutAir/VolumeFraction_withoutAir;		// required as the VolumeFraction_withoutAir will not sum to one
		*Effective_RadioactiveHeat_withoutAir 		= *Effective_RadioactiveHeat_withoutAir/VolumeFraction_withoutAir;			// required as the VolumeFraction_withoutAir will not sum to one

	}
	else{
		// only air - EffectiveDensity_withoutAir will be zero

	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/




/*==========================================================================================================*/
/* Average effective viscosity from phases in the cell, both WITH and WITHOUT sticky air				   */
#undef __FUNCT__
#define __FUNCT__ "Average_EffectiveViscosity_FDSTAG_Cell"
PetscErrorCode Average_EffectiveViscosity_FDSTAG_Cell( const double VolumeFraction_phases[], UserContext *user,
		double T, double E2nd, double T2nd, double P, double P_lithos,
		double PlasticStrain,
		double *EffectiveViscosity_withAir, double *EffectiveViscosity_withoutAir)
{
	PetscErrorCode 	ierr;
	PetscInt 		iphase, AirPhase;
	PetscScalar		VolumeFraction, mu_viscous, mu_plastic, mu_eff;
	PetscScalar	 	VolumeFraction_withoutAir;


	AirPhase 									=	user->ErosionParameters.StickyAirPhase;

	*EffectiveViscosity_withAir 				=	0;
	*EffectiveViscosity_withoutAir 				=	0;
	VolumeFraction_withoutAir 					=	0;
	for (iphase=0; iphase<user->num_phases; iphase++){
		VolumeFraction = VolumeFraction_phases[iphase];	// volume fraction of PHASE iphase

		PetscInt PlasticityCutoff=0;

		/* Compute effective viscosity for this particular phase */
		ierr = ComputeEffectiveViscosity(iphase, user, T, E2nd, T2nd,  P, P_lithos, PlasticStrain, PlasticityCutoff,
				&mu_viscous, &mu_plastic, &mu_eff); CHKERRQ(ierr);

		*EffectiveViscosity_withAir	+= 	VolumeFraction*mu_viscous;		// arithmetic averaging

		if (iphase != AirPhase){
			// sum over non-air phases

			VolumeFraction_withoutAir 		+= 	VolumeFraction;
			*EffectiveViscosity_withoutAir	+= 	VolumeFraction*mu_viscous;		// arithmetic averaging
		}

	}

	// Correct phase w/out air
	if (VolumeFraction_withoutAir>0)
	{
		*EffectiveViscosity_withoutAir = *EffectiveViscosity_withoutAir/VolumeFraction_withoutAir;		// required as the VolumeFraction_withoutAir will not sum to one
	}
	else{
		// only air - EffectiveViscosity_withoutAir will be zero

	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Computes effective viscosity */
#undef __FUNCT__
#define __FUNCT__ "ComputeEffectiveViscosity"
PetscErrorCode ComputeEffectiveViscosity( PetscInt phase, UserContext *user, double T, double E2nd, double T2nd, double P, double P_lithos,
		double PlasticStrain, PetscInt PlasticityCutoff,
		double *mu_viscous, double *mu_plastic, double *mu_eff)
{
	PetscScalar 	exponent1, exponent2, Sigma_yield, Phi, Cohesion, N_exp, e0, P_upper, P_lower, A, E, R, T_dim, E2nd_dim;
	PetscScalar		Strain_begin, Strain_end, Cohesion_end, Phi_end, factor;
	PetscInt		ViscosityLaw, PlasticityLaw;


	/* ------------------------------------------------------------------------------------------------ */
	ViscosityLaw	=	user->PhaseProperties.ViscosityLaw[phase];

	if		(ViscosityLaw==1){
		/* Constant viscosity */
		*mu_viscous		=	user->PhaseProperties.mu[phase];

	}
	else if (ViscosityLaw==2){
		/* Powerlaw w/out T-dependency */
		N_exp 			=	user->PhaseProperties.n_exponent[phase];
		e0 				=	user->PhaseProperties.Powerlaw_e0[phase];

		exponent1 		= 	(1.0/N_exp)-1.0;

		if ((E2nd>0.0) && (PetscBool) isnan(E2nd) == PETSC_FALSE) {
			*mu_viscous 	= 	user->PhaseProperties.mu[phase]*PetscPowScalar( (E2nd/e0),exponent1);

		}
		else {
			*mu_viscous 	= 	user->PhaseProperties.mu[phase];
		}

	}
	else if (ViscosityLaw==3){
		/* Frank-Kamenetskii approximation */
		*mu_viscous	=	user->PhaseProperties.mu[phase]*(PetscExpScalar(-user->PhaseProperties.FrankKamenetskii[phase]*T));

	}
	else if (ViscosityLaw==4){
		/* temperature-dependent powerlaw viscosity */
		N_exp 			=	user->PhaseProperties.n_exponent[phase];
		e0 				=	user->PhaseProperties.Powerlaw_e0[phase]*user->Characteristic.Strainrate;
		A				=	user->PhaseProperties.A[phase]*((PetscPowScalar(user->Characteristic.Stress,(-user->PhaseProperties.n_exponent[phase])))/user->Characteristic.Time);
		E				=	user->PhaseProperties.E[phase]*user->Characteristic.Jmol;
		R				=	user->GasConstant;
		T_dim			=	T*user->Characteristic.Temperature;
		E2nd_dim		=	E2nd*user->Characteristic.Strainrate;

		exponent1 		= 	(1.0/N_exp)-1.0;
		exponent2 		= 	-(1.0/N_exp);

		if ((E2nd>0.0) && (PetscBool) isnan(E2nd) == PETSC_FALSE) {
			*mu_viscous 	= 	PetscPowScalar(A,exponent2)*PetscPowScalar((E2nd_dim/e0),exponent1)*PetscExpScalar(E/(R*N_exp*T_dim));
		}
		else {
			*mu_viscous 	= 	PetscPowScalar(A,exponent2)*PetscExpScalar(E/(R*N_exp*T_dim));
		}

		*mu_viscous = *mu_viscous/user->Characteristic.Viscosity;

	}
	else if (ViscosityLaw==5){
		/* temperature-dependent viscosity - without power-law */
		N_exp 			=	user->PhaseProperties.n_exponent[phase];
		e0 				=	user->PhaseProperties.Powerlaw_e0[phase]*user->Characteristic.Strainrate;
		A				=	user->PhaseProperties.A[phase]*((PetscPowScalar(user->Characteristic.Stress,(-user->PhaseProperties.n_exponent[phase])))/user->Characteristic.Time);
		E				=	user->PhaseProperties.E[phase]*user->Characteristic.Jmol;
		R				=	user->GasConstant;
		T_dim			=	T*user->Characteristic.Temperature;
		E2nd_dim		=	1e-16;

		exponent1 		= 	(1.0/N_exp)-1.0;
		exponent2 		= 	-(1.0/N_exp);

		*mu_viscous 	= 	PetscPowScalar(A,exponent2)*PetscPowScalar((E2nd_dim/e0),exponent1)*PetscExpScalar(E/(R*N_exp*T_dim));

		*mu_viscous = *mu_viscous/user->Characteristic.Viscosity;

	}
	else {
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP," There is something fuzzy about your VISCOSITY law! I have law=%i for phase=%i \n",ViscosityLaw, phase);
	}
	/* ------------------------------------------------------------------------------------------------ */

	/*
	 * If the element is non-plastic, we update the effective viscosity.
	 * If it is plastic already, we check whether the yield stress is exceeded already, and
	 * 	(1) if yes, update mu_plastic and mu_eff to represent new yield stress
	 *  (2) if no, leave as is (viscosity is fine already and reduced stresses to be at or just below yield)
	 */
	//if (PetscAbsScalar(mu_plastic_initial)<1e-15){
	*mu_eff 		=	*mu_viscous;
	//}


	/* ------------------------------------------------------------------------------------------------ */
	PlasticityLaw	=	user->PhaseProperties.PlasticityLaw[phase];
	if		(PlasticityLaw==0){
		/* No Plasticity */
		*mu_plastic =	0.0;

	}
	else if		(PlasticityLaw==1){
		/* Drucker-Prager */
		Cohesion		=   user->PhaseProperties.Cohesion[phase];
		Phi 			=  	user->PhaseProperties.FrictionAngle[phase];

		Strain_end 		=	user->PhaseProperties.Weakening_PlasticStrain_End[phase];

		if (Strain_end>0.0){
			// Strain weakening should be taken into account.

			Strain_begin 	=	user->PhaseProperties.Weakening_PlasticStrain_Begin[phase];
			Cohesion_end	=	user->PhaseProperties.CohesionAfterWeakening[phase];
			Phi_end 		=	user->PhaseProperties.FrictionAngleAfterWeakening[phase];

			if (PlasticStrain>Strain_end){
				Cohesion 	=	Cohesion_end;
				Phi 		=	Phi_end;
			}
			else if ( (PlasticStrain>Strain_begin) && (PlasticStrain<=Strain_end)){
				factor 		=	(PlasticStrain-Strain_begin)/(Strain_end-Strain_begin);
				Cohesion 	=	(1.0-factor)*Cohesion 	+ factor*Cohesion_end;
				Phi 		=	(1.0-factor)*Phi 		+ factor*Phi_end;
			}

		}


		if (E2nd>0.0){

			// Create a check that pressure is within reasonable limits
			//
			P_upper         =   -( P_lithos + Cohesion*cos(Phi/180*M_PI))/(sin(Phi/180*M_PI)-1);
			P_lower         =   -(-P_lithos + Cohesion*cos(Phi/180*M_PI))/(sin(Phi/180*M_PI)+1);
			if (P_lower<0.0){	P_lower=0.0; }
			if (PlasticityCutoff==1){
				if (P>P_upper){
					P 			=	P_upper;
				}
				else if (P<P_lower){
					P			=	P_lower;
				}

			}

			// Compute yield stress:
			Sigma_yield     =   P*sin( Phi/360.0*2.0*M_PI ) + Cohesion*cos( Phi/360.0*2.0*M_PI );

			// cutoff for small extensional stresses
			if (Sigma_yield<0.01*Cohesion){
				Sigma_yield = 0.01*Cohesion;
			}

			if( (T2nd>Sigma_yield) && (T2nd>0.0) ){
				*mu_plastic = 	Sigma_yield/(2.0*E2nd);
				*mu_eff 	= 	*mu_plastic;
			}

		}
		else {
			*mu_plastic = 0.0;
		}


	}
	else {
		PetscPrintf(PETSC_COMM_WORLD," There is something fuzzy about your PLASTICITY law! %lld \n", (LLD)PlasticityLaw);
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	/* ------------------------------------------------------------------------------------------------ */


	/* Perform cutoffs */
	if (*mu_eff<user->LowerViscosityCutoff){
		*mu_eff			=	user->LowerViscosityCutoff;
	}
	else if (*mu_eff>user->UpperViscosityCutoff){
		*mu_eff			=	user->UpperViscosityCutoff;
	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Computes effective density */
#undef __FUNCT__
#define __FUNCT__ "ComputeEffectiveDensity"
PetscErrorCode ComputeEffectiveDensity( PetscInt phase, UserContext *user, double T, double *rho_eff)
{
	PetscScalar 	alpha, T0, Ra;
	PetscInt		DensityLaw;


	DensityLaw	=	user->PhaseProperties.DensityLaw[phase];

	if (DensityLaw==1){

		/* Temperature-dependent density
		 * ------------------------------------------------------------------------------------------------ */
		alpha 			=	user->PhaseProperties.ThermalExpansivity[phase];
		T0 				=	0;			// should typically be 273 K -> Add this to the code !!!


		/* T-dependent density */
		*rho_eff		=	user->PhaseProperties.rho[phase]*(1.0 - alpha*(T-T0) );
	}
	else if (DensityLaw==2){
		// We do non-dimensional Rayleigh-Bernard convection in which case the rhs should be  -Ra*T  (and not rho*(1-dT))
		Ra 			=	user->PhaseProperties.Ra[phase];


		/* T-dependent density */
		*rho_eff		=	-Ra*T;

	}
	else if (DensityLaw==3){
		// We use an artificial temperature to set the density
		*rho_eff		=	T * user->Characteristic.Temperature / user->Characteristic.Density;

	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown density law");
	}


	//	PetscPrintf(PETSC_COMM_WORLD,"ComputeEffectiveDensity: alpha=%f, T=%f rho=%f \n",alpha, T, *rho_eff);



	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* Update Material properties @ FDSTAG points
 *
 * This updates for example effective viscosity and density at points where we need this info in order to compute the FDSTAG
 * stencil (for residual calculation and for computing the Jacobian).
 * This routine can be called during every  nonlinear iteration.
 *
 * Before calling the routine, we
 *
 * (1) need to compute history-variables at all relevant points by calling ComputeHistoryPropertiesFromParticles_FDSTAG
 * (2) need to update the second invariant of strain rate and stresses & interpolate it to all nodes (so it is available everywhere)
 *
 */
#undef __FUNCT__
#define __FUNCT__ "UpdateMaterialProperties_FDSTAG"
PetscErrorCode UpdateMaterialProperties_FDSTAG(UserContext *user)
{
	// ... p o r n o ...
	// "Real Programmer" writes code that every "idiot" can understand.
	// Best strategy is to consider yourself an "idiot" when you will be looking
	// at what you have written a week/month/year later.

	PetscErrorCode 			ierr;
	DM 						da, cda, cda_SurfaceTopo;
	Vec 					*Temperature_Vec, *Pressure_Vec, *E2nd_Vec, *T2nd_Vec, *PlasticStrain_Vec, *Strain_Vec;
	Vec 					*HeatCapacity_Vec,*Conductivity_Vec, *RadioactiveHeat_Vec;
	Vec						*EffectiveViscosity_Vec, *EffectiveDensity_Vec, *PhaseProportions_Vec[max_num_phases];
	Vec 					gc, gc_SurfaceTopo, LocalSurfaceTopography_vec;
	PetscInt 				ix,iy,iz,xs,ys,zs,xm,ym,zm, iphase, i, da_loop;
	PetscInt				ielx,iely,ielz, AirPhase;
	PetscInt 				xs_free,ys_free,zs_free,xm_free,ym_free,zm_free, ix_free, iy_free;
	DMDACoor3d				***coors, ***coors_SurfaceTopo;
	PetscScalar				VolumeFraction, E2nd, T2nd, P, P_lithos, PlasticStrain, T;
	PetscScalar				***Local_Viscosity,		***PhaseProportions_local[max_num_phases];
	PetscScalar				***Local_Temperature, 	***Local_Pressure, 		***Local_E2nd;
	PetscScalar				***Local_T2nd,			***Local_PlasticStrain, ***Local_Strain;
	PetscScalar				***Local_Density,		***Local_RadioactiveHeat, ***Local_Conductivity;
	PetscScalar				***Local_HeatCapacity,	***materials_fine_array,	***Local_NumParticles;
	PetscScalar				***LocalSurfaceTopography;
	PetscScalar			 	*VolumeFraction_phases;
	MaterialsElementDynamic material_fine_data;
	PetscMPIInt 			rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	HeatCapacity_Vec    = NULL;
	Conductivity_Vec    = NULL;
	RadioactiveHeat_Vec = NULL;

	/* PART 1: UPDATE PROPERTIES AT FDSTAG POINTS ----------------------------------------------------------*/

	/* We need the same information at Center, XY, XZ & YZ points.
	 * At Nodal points we also need information used in the temperature points.
	 * Instead of writing several subroutines, that do essentially the same,
	 * we here loop over all DA's & perform the same actions.
	 */
	ierr = PetscMalloc( (size_t)(user->num_phases)*sizeof(PetscScalar), &VolumeFraction_phases); CHKERRQ(ierr);

	for (da_loop=0; da_loop<5; da_loop++){

		/* Retrieve the correct vectors and DA's for this DMDA */
		if (da_loop==0){ 			/* Center points */
			da 						= 	user->FDSTAG.DA_CENTER;
			for (i=0; i<user->num_phases; i++){
				PhaseProportions_Vec[i]	= 	&(user->FDSTAG.Center_PhaseProportions[i]);
			}
			Temperature_Vec 		= 	&(user->FDSTAG.Center_Temperature			);
			Pressure_Vec 			=	&(user->FDSTAG.Center_Pressure				);
			E2nd_Vec 				=	&(user->FDSTAG.Center_E2nd					);
			T2nd_Vec 				=	&(user->FDSTAG.Center_T2nd					);
			PlasticStrain_Vec 		=	&(user->FDSTAG.Center_PlasticStrain			);
			Strain_Vec 				=	&(user->FDSTAG.Center_Strain				);
			EffectiveViscosity_Vec 	=	&(user->FDSTAG.Center_EffectiveViscosity 	);
			EffectiveDensity_Vec 	=	&(user->FDSTAG.Center_Density 				);
		}
		else if (da_loop==1){		/* XY points */
			da 						= 	user->FDSTAG.DA_XY_POINTS;
			for (i=0; i<user->num_phases; i++){
				PhaseProportions_Vec[i]	= 	&(user->FDSTAG.XYPoints_PhaseProportions[i]);
			}
			Temperature_Vec 		= 	&(user->FDSTAG.XYPoints_Temperature				);
			Pressure_Vec 			=	&(user->FDSTAG.XYPoints_Pressure				);
			E2nd_Vec 				=	&(user->FDSTAG.XYPoints_E2nd					);
			T2nd_Vec 				=	&(user->FDSTAG.XYPoints_T2nd					);
			PlasticStrain_Vec 		=	&(user->FDSTAG.XYPoints_PlasticStrain			);
			Strain_Vec 				=	&(user->FDSTAG.XYPoints_Strain					);
			EffectiveViscosity_Vec 	=	&(user->FDSTAG.XYPoints_EffectiveViscosity 		);
			EffectiveDensity_Vec 	=	&(user->FDSTAG.XYPoints_Density 				);
		}
		else if (da_loop==2){		/* XZ points */
			da 						= 	user->FDSTAG.DA_XZ_POINTS;
			for (i=0; i<user->num_phases; i++){
				PhaseProportions_Vec[i]	= 	&(user->FDSTAG.XZPoints_PhaseProportions[i]);
			}
			Temperature_Vec 		= 	&(user->FDSTAG.XZPoints_Temperature				);
			Pressure_Vec 			=	&(user->FDSTAG.XZPoints_Pressure				);
			E2nd_Vec 				=	&(user->FDSTAG.XZPoints_E2nd					);
			T2nd_Vec 				=	&(user->FDSTAG.XZPoints_T2nd					);
			PlasticStrain_Vec 		=	&(user->FDSTAG.XZPoints_PlasticStrain			);
			Strain_Vec 				=	&(user->FDSTAG.XZPoints_Strain					);
			EffectiveViscosity_Vec 	=	&(user->FDSTAG.XZPoints_EffectiveViscosity 		);
			EffectiveDensity_Vec 	=	&(user->FDSTAG.XZPoints_Density 				);
		}
		else if (da_loop==3){		/* YZ points */
			da 						= 	user->FDSTAG.DA_YZ_POINTS;
			for (i=0; i<user->num_phases; i++){
				PhaseProportions_Vec[i]	= 	&(user->FDSTAG.YZPoints_PhaseProportions[i]);
			}
			Temperature_Vec 		= 	&(user->FDSTAG.YZPoints_Temperature				);
			Pressure_Vec 			=	&(user->FDSTAG.YZPoints_Pressure				);
			E2nd_Vec 				=	&(user->FDSTAG.YZPoints_E2nd					);
			T2nd_Vec 				=	&(user->FDSTAG.YZPoints_T2nd					);
			PlasticStrain_Vec 		=	&(user->FDSTAG.YZPoints_PlasticStrain			);
			Strain_Vec 				=	&(user->FDSTAG.YZPoints_Strain					);
			EffectiveViscosity_Vec 	=	&(user->FDSTAG.YZPoints_EffectiveViscosity 		);
			EffectiveDensity_Vec 	=	&(user->FDSTAG.YZPoints_Density 				);
		}
		else if (da_loop==4){		/* Corner points */
			da 						= 	user->FDSTAG.DA_CORNER;
			for (i=0; i<user->num_phases; i++){
				PhaseProportions_Vec[i]	= 	&(user->FDSTAG.Corner_PhaseProportions[i]);
			}
			Temperature_Vec 		= 	&(user->FDSTAG.Corner_Temperature				);
			Pressure_Vec 			= 	&(user->FDSTAG.Corner_Pressure					);
			EffectiveDensity_Vec 	=	&(user->FDSTAG.Corner_Density 					);
			HeatCapacity_Vec 		=	&(user->FDSTAG.Corner_HeatCapacity 				);
			Conductivity_Vec 		=	&(user->FDSTAG.Corner_Conductivity 				);
			RadioactiveHeat_Vec 	=	&(user->FDSTAG.Corner_RadioactiveHeat 			);
		}

		/* Retrieve local vectors and arrays that will later contain
		 *
		 * 	- effective viscosity
		 * 	- density
		 *
		 * */

		for (i=0; i<user->num_phases; i++){
			ierr = DMDAVecGetArray(da,   *PhaseProportions_Vec[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr); // Phase proportions
		}
		ierr = DMDAVecGetArray(da,   	*Temperature_Vec, 				&Local_Temperature	); 		CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da,   	*Pressure_Vec, 					&Local_Pressure	 	); 		CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da,   *EffectiveDensity_Vec, 			&Local_Density	 	); 		CHKERRQ(ierr);

		if (da_loop<4){
			ierr = DMDAVecGetArray(da,   	*E2nd_Vec, 						&Local_E2nd	); 				CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   	*T2nd_Vec, 						&Local_T2nd	); 				CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   	*PlasticStrain_Vec, 			&Local_PlasticStrain); 		CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   	*Strain_Vec, 					&Local_Strain		); 		CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   *EffectiveViscosity_Vec, 	&Local_Viscosity	); 		CHKERRQ(ierr);
		}
		else if (da_loop==4){
			ierr = DMDAVecGetArray(da,   *HeatCapacity_Vec, 			&Local_HeatCapacity		);	CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   *Conductivity_Vec, 			&Local_Conductivity		);	CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   *RadioactiveHeat_Vec, 		&Local_RadioactiveHeat	);	CHKERRQ(ierr);
		}

		/* Get coordinates of the DA, including ghost points */
		ierr = DMGetCoordinateDM(da,			&cda); 		CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,   	&gc); 		CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,gc,			&coors); 	CHKERRQ(ierr);

		/* Make a loop over every FD node */
		ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
		for (iz=zs; iz<zs+zm; iz++){
			for (iy=ys; iy<ys+ym; iy++){
				for (ix=xs; ix<xs+xm; ix++){

					PetscScalar EffectiveDensity_withAir, 	EffectiveDensity_withoutAir;

					/* Retrieve properties for this node that are required to compute effective viscosity & plasticity */
					E2nd = 0; T2nd=0;
					if (da_loop<4){
						E2nd 			= 	Local_E2nd[iz][iy][ix];
						T2nd 			= 	Local_T2nd[iz][iy][ix];
						PlasticStrain 	= 	0;
					}
					P 				= 	Local_Pressure[iz][iy][ix];
					T 				=	Local_Temperature[iz][iy][ix];
					P_lithos 		=	0;		// improve this

					// Compute volume fractions of every phase in this FDSTAG cell
					for (iphase=0; iphase<user->num_phases; iphase++){
						VolumeFraction_phases[iphase] = PhaseProportions_local[iphase][iz][iy][ix];
					}

					/* Loop over all phases */

					/* Set properties to the relevant points */
					if (da_loop<4){
						PetscScalar EffectiveViscosity_withAir, EffectiveViscosity_withoutAir;

						//  Compute effective viscosity with and without air
						ierr = Average_EffectiveViscosity_FDSTAG_Cell(VolumeFraction_phases, user, T, E2nd, T2nd,  P, P_lithos, PlasticStrain,
								&EffectiveViscosity_withAir, &EffectiveViscosity_withoutAir); CHKERRQ(ierr);

						Local_Viscosity[iz][iy][ix] 	=		EffectiveViscosity_withAir;
					}
					else if (da_loop==4){	/* Thermal quantities @ corner nodes */
						PetscScalar 	Effective_HeatCapacity_withAir, 		Effective_HeatCapacity_withoutAir;
						PetscScalar 	Effective_ThermalConductivity_withAir, 	Effective_ThermalConductivity_withoutAir;
						PetscScalar 	Effective_RadioactiveHeat_withAir, 		Effective_RadioactiveHeat_withoutAir;


						// Effective thermal properties with and without air
						ierr = Average_EffectiveThermalProperties_FDSTAG_Cell(VolumeFraction_phases, user,
								&Effective_HeatCapacity_withAir, 		&Effective_HeatCapacity_withoutAir,
								&Effective_ThermalConductivity_withAir, &Effective_ThermalConductivity_withoutAir,
								&Effective_RadioactiveHeat_withAir, 	&Effective_RadioactiveHeat_withoutAir); CHKERRQ(ierr);


						Local_HeatCapacity[iz][iy][ix]		=		Effective_HeatCapacity_withAir;
						Local_Conductivity[iz][iy][ix]		=		Effective_ThermalConductivity_withAir;
						Local_RadioactiveHeat[iz][iy][ix]	=		Effective_RadioactiveHeat_withAir;

					}

					ierr = Average_EffectiveDensity_FDSTAG_Cell(VolumeFraction_phases, user, T, &EffectiveDensity_withAir, &EffectiveDensity_withoutAir); CHKERRQ(ierr);

					Local_Density[iz][iy][ix] 			=		EffectiveDensity_withAir;

				}
			}
		}

		if (user->ErosionParameters.UseInternalFreeSurface != 0){
			/* PART 2: TAKE INTERNAL FREE SURFACE INTO ACCOUNT -----------------------------------------------------*/
			// check location of free surface and determine if a node is above or below it. If it is above, it will
			// get the 'sticky air' properties. If it is below, we will average density & viscosity to take the Topo into account

			// NOTE: we use a double loop here (one over the da, one over the surface) and check every surface element with an inpolygon algorithm.
			// we could make it faster by making a rudimentary check before calling the inpolygon algorithm.

			ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
			ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,					&gc_SurfaceTopo			); 	CHKERRQ(ierr);
			ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);

			/* Copy GHOST points!!!! */
			ierr = DMGetLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);
			ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);

			ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
			ierr = DMDAGetGhostCorners(user->DA_SurfaceTopography,&xs_free,&ys_free,&zs_free,&xm_free,&ym_free,&zm_free); CHKERRQ(ierr);

			AirPhase = user->ErosionParameters.StickyAirPhase;

			for (iy=ys; iy<(ys+ym); iy++){		// loop over grid nodes			[which can be center, nodal, xy, etc. points]
				for (ix=xs; ix<(xs+xm); ix++){	// loop over da grid nodes
					PetscScalar 	x,y,z_FreeSurface;
					PetscBool		inside;

					/* x & y coordinates of da mesh*/
					x = coors[zs][iy][ix].x;
					y = coors[zs][iy][ix].y;

					// slightly correct points at boundaries, to prevent roundoff errors
					if (ix==0){
						x=x	+ 1e-3;
					}
					//else if (ix==user->nnode_x-1){
					else if (ix==xs+xm-1){

						x=x	- 1e-3;
					}
					if (iy==0){
						y=y	+ 1e-3;
					}
					else if (iy==ys+ym-1){
						y=y	- 1e-3;
					}

				//	PetscPrintf(PETSC_COMM_SELF,"rank=%i ix=%i, iy=%i, x,y=[%f,%f] xs+xm=%i xs=%i, user->nnode_x-1=%i\n",rank, ix,iy,x,y,xs+xm, xs,user->nnode_x-1);

					/* Loop over free surface [including ghost points!] */
					for (iy_free=ys_free; iy_free<(ys_free + ym_free - 1); iy_free++){
						for(ix_free=xs_free; ix_free<(xs_free + xm_free - 1); ix_free++){

							PetscScalar		p_x[5],  p_y[5],  p_z[5], mu_Air, rho_Air;

							/* Extract air properties */
							mu_Air 				=	user->PhaseProperties.mu[AirPhase];
							rho_Air 			=	user->PhaseProperties.rho[AirPhase];


							/* Extract a contour of the linear element and use an isinside routine to find out if the point is inside the element */
							p_x[0] = coors_SurfaceTopo		[zs_free][iy_free  ][ix_free  ].x;
							p_y[0] = coors_SurfaceTopo		[zs_free][iy_free  ][ix_free  ].y;
							p_z[0] = LocalSurfaceTopography	[zs_free][iy_free  ][ix_free  ];

							p_x[1] = coors_SurfaceTopo		[zs_free][iy_free  ][ix_free+1].x;
							p_y[1] = coors_SurfaceTopo		[zs_free][iy_free  ][ix_free+1].y;
							p_z[1] = LocalSurfaceTopography	[zs_free][iy_free  ][ix_free+1];

							p_x[2] = coors_SurfaceTopo		[zs_free][iy_free+1][ix_free+1].x;
							p_y[2] = coors_SurfaceTopo		[zs_free][iy_free+1][ix_free+1].y;
							p_z[2] = LocalSurfaceTopography	[zs_free][iy_free+1][ix_free+1];

							p_x[3] = coors_SurfaceTopo		[zs_free][iy_free+1][ix_free  ].x;
							p_y[3] = coors_SurfaceTopo		[zs_free][iy_free+1][ix_free  ].y;
							p_z[3] = LocalSurfaceTopography	[zs_free][iy_free+1][ix_free  ];
							p_x[4] = p_x[0];
							p_y[4] = p_y[0];
							p_z[4] = p_z[0];



							/* Is point inside square polygon? */
							inside = pnpoly(4, p_x, p_y, x,y);
							if (inside){
								/* yup, it's inside or on the edge.
								 * Now interpolate the height of the free surface */
								z_FreeSurface =  InterpolateWithin2DLinearElement(p_x, p_y, p_z, x, y);

								if ((z_FreeSurface > coors[zs][iy][ix].z) && (z_FreeSurface < coors[zs+zm-1][iy][ix].z)){

									/* Yup, free surface is on this PROC
									 *
									 * We program this such, that it takes into account a refinement in z-direction, which is why
									 * we here use a loop.
									 */
									PetscInt	iz_FreeSurface;
									PetscScalar	fac, dz;

									iz_FreeSurface = -1;
									for (iz=zs; iz<zs+zm-1; iz++){
										if ((z_FreeSurface > coors[iz][iy][ix].z) && (z_FreeSurface <= coors[iz+1][iy][ix].z)){
											iz_FreeSurface = iz;
										}
									}
									if (iz_FreeSurface<0){
										PetscPrintf(PETSC_COMM_WORLD,"iz_FreeSurface=%i",iz_FreeSurface);
									}


// which loop am I in, by the way?

									/* Compute the fraction that is attributed to the lower node [below free surface] and the one that goes to the upper node [air]*/
									dz 	= coors[iz_FreeSurface+1][iy][ix].z-coors[iz_FreeSurface][iy][ix].z;
									fac = (z_FreeSurface-coors[iz_FreeSurface][iy][ix].z)/dz;		// fraction that goes to the lower node

									if ((fac<0) || (fac>1.0)){
										SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "fac should be [0-1] and was %f \n",fac );
									}

									// Compute volume fractions of every phase in this FDSTAG cell
									for (iphase=0; iphase<user->num_phases; iphase++){
										VolumeFraction_phases[iphase] = PhaseProportions_local[iphase][iz_FreeSurface][iy][ix];
									}


									// Compute Effective viscosity without air phase
									if (da_loop<4){
										E2nd 			= 	Local_E2nd[iz_FreeSurface][iy][ix];
										T2nd 			= 	Local_T2nd[iz_FreeSurface][iy][ix];
										PlasticStrain 	= 	0;
									}
									P 					= 	Local_Pressure[iz_FreeSurface][iy][ix];
									T 					=	Local_Temperature[iz_FreeSurface][iy][ix];


									if (da_loop<4){
										PetscScalar EffectiveViscosity_withAir, EffectiveViscosity_withoutAir;

										//  Compute effective viscosity with and without air
										ierr = Average_EffectiveViscosity_FDSTAG_Cell(VolumeFraction_phases, user, T, E2nd, T2nd,  P, P_lithos, PlasticStrain,
												&EffectiveViscosity_withAir, &EffectiveViscosity_withoutAir); CHKERRQ(ierr);

										Local_Viscosity[iz_FreeSurface][iy][ix] 	=		    fac*EffectiveViscosity_withoutAir + (1-fac)*mu_Air;


									}
									else if (da_loop==4){	/* Thermal quantities @ corner nodes */
										PetscScalar 	k_Air, Cp_Air, Q_Air;
										PetscScalar 	Effective_HeatCapacity_withAir, 		Effective_HeatCapacity_withoutAir;
										PetscScalar 	Effective_ThermalConductivity_withAir, 	Effective_ThermalConductivity_withoutAir;
										PetscScalar 	Effective_RadioactiveHeat_withAir, 		Effective_RadioactiveHeat_withoutAir;


										// Effective thermal properties with and without air
										ierr = Average_EffectiveThermalProperties_FDSTAG_Cell(VolumeFraction_phases, user,
												&Effective_HeatCapacity_withAir, 		&Effective_HeatCapacity_withoutAir,
												&Effective_ThermalConductivity_withAir, &Effective_ThermalConductivity_withoutAir,
												&Effective_RadioactiveHeat_withAir, 	&Effective_RadioactiveHeat_withoutAir); CHKERRQ(ierr);


										k_Air 				=	user->PhaseProperties.T_Conductivity[AirPhase];
										Cp_Air				=	user->PhaseProperties.HeatCapacity[AirPhase];
										Q_Air				=	user->PhaseProperties.RadioactiveHeat[AirPhase];

										// Correct thermal material properties based on a volume-of-fluid approach that takes the location of the free
										// surface into account
										Local_HeatCapacity[iz_FreeSurface][iy][ix]		=		    fac*Effective_HeatCapacity_withoutAir 			+ (1-fac)*Cp_Air;
										Local_Conductivity[iz_FreeSurface][iy][ix]		=		    fac*Effective_ThermalConductivity_withoutAir	+ (1-fac)*k_Air;
										Local_RadioactiveHeat[iz_FreeSurface][iy][ix]	=		    fac*Effective_RadioactiveHeat_withoutAir 		+ (1-fac)*Q_Air;

									}

									// Density
									PetscScalar EffectiveDensity_withAir, 	EffectiveDensity_withoutAir;

									ierr = Average_EffectiveDensity_FDSTAG_Cell(VolumeFraction_phases, user, T, &EffectiveDensity_withAir, &EffectiveDensity_withoutAir); CHKERRQ(ierr);
									Local_Density[iz_FreeSurface  ][iy][ix] 			=		    fac*EffectiveDensity_withoutAir + (1-fac)*rho_Air;


									/* Set everything above to air properties */
									for (iz=iz_FreeSurface+1; iz<(zs+zm); iz++){
										if (da_loop<4){
											Local_Viscosity[iz][iy][ix] 	=		mu_Air;
										}
										else if (da_loop==4){	/* Thermal quantities @ corner nodes */
											PetscScalar k_Air, Cp_Air, Q_Air;

											k_Air 								=		user->PhaseProperties.T_Conductivity [AirPhase];
											Cp_Air								=		user->PhaseProperties.HeatCapacity   [AirPhase];
											Q_Air								=		user->PhaseProperties.RadioactiveHeat[AirPhase];

											Local_HeatCapacity[iz][iy][ix]		=		Cp_Air;
											Local_Conductivity[iz][iy][ix]		=		k_Air;
											Local_RadioactiveHeat[iz][iy][ix]	=		Q_Air;

										}
										Local_Density[iz][iy][ix] 				=		rho_Air;
									}
								}



								/* In case the whole PROC is above the free surface */
								if (coors[zs][iy][ix].z >= z_FreeSurface){

									for (iz=zs; iz<(zs+zm); iz++){
										/* All material properties get the 'air' phase */
										if (da_loop<4){
											Local_Viscosity[iz][iy][ix] 	=		mu_Air;
										}
										else if (da_loop==4){	/* Thermal quantities @ corner nodes */
											PetscScalar 	k_Air, Cp_Air, Q_Air;

											k_Air 								=		user->PhaseProperties.T_Conductivity[AirPhase];
											Cp_Air								=		user->PhaseProperties.HeatCapacity[AirPhase];
											Q_Air								=		user->PhaseProperties.RadioactiveHeat[AirPhase];

											Local_HeatCapacity[iz][iy][ix]		=		Cp_Air;
											Local_Conductivity[iz][iy][ix]		=		k_Air;
											Local_RadioactiveHeat[iz][iy][ix]	=		Q_Air;

										}
										Local_Density[iz][iy][ix] 				=		rho_Air;
									}

								}

							}	// isinside? well, I don't know ...




						}
					}

				}
			}

			ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
			ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
			ierr	= 	DMRestoreLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);

		}
		/* END OF PART 2: TAKE INTERNAL FREE SURFACE INTO ACCOUNT ----------------------------------------------*/



		/* Restore arrays */
		for (i=0; i<user->num_phases; i++ ){
			ierr = DMDAVecRestoreArray(da,   	*PhaseProportions_Vec[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr);
		}
		ierr = DMDAVecRestoreArray(da,   	*Temperature_Vec, 			&Local_Temperature	); 		CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(da,   	*Pressure_Vec, 				&Local_Pressure	 	); 		CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(da,   	*EffectiveDensity_Vec, 		&Local_Density		); 		CHKERRQ(ierr);

		if (da_loop<4){
			ierr = DMDAVecRestoreArray(da,   	*E2nd_Vec, 					&Local_E2nd			); 		CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   	*EffectiveViscosity_Vec, 	&Local_Viscosity	); 		CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   	*T2nd_Vec, 					&Local_T2nd			); 		CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   	*PlasticStrain_Vec, 		&Local_PlasticStrain); 		CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   	*Strain_Vec, 				&Local_Strain		); 		CHKERRQ(ierr);
		}
		else if (da_loop==4){
			ierr = DMDAVecRestoreArray(da,   *HeatCapacity_Vec, 			&Local_HeatCapacity		);	CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   *Conductivity_Vec, 			&Local_Conductivity		);	CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(da,   *RadioactiveHeat_Vec, 		&Local_RadioactiveHeat	);	CHKERRQ(ierr);
		}

		ierr 	= 	DMDAVecRestoreArray(cda,gc,&coors); 	CHKERRQ(ierr);

	}
	/* END OF PART 1 & 2: UPDATE PROPERTIES AT FDSTAG POINTS ---------------------------------------------------*/



	/* PART 3: UPDATE PROPERTIES ON INTEGRATION POINTS -----------------------------------------------------
	 *  This is done for consistency with the (FE-based) rest of the code (and in order to not have to change
	 *  all I/O routines)
	 */

	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_EffectiveViscosity , 	&Local_Viscosity		); 		CHKERRQ(ierr); // Viscosity @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_Density , 			&Local_Density	 		); 		CHKERRQ(ierr); // Density @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_NumParticles , 		&Local_NumParticles	 	); 		CHKERRQ(ierr); // # of particles @ center
	for (i=0; i<user->num_phases; i++){
		ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,   user->FDSTAG.Center_PhaseProportions[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr); // Phase proportions
	}

	ierr = DMDAVecGetArray(user->DA_Materials,   user->Materials, &materials_fine_array); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	/* Loop over all elements */
	for (ielz=zs; ielz<zs+zm; ielz++){
		for (iely=ys; iely<ys+ym; iely++){
			for (ielx=xs; ielx<xs+xm; ielx++){
				PetscScalar eta_mean, rho_mean, max_phase_fraction, num_particles;
				PetscInt 	ngp_vel=1, intp, phase_dominant;

				LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs,iely-ys,ielz-zs, ngp_vel, materials_fine_array );

				/* Extract data of current element */
				eta_mean 		= 	Local_Viscosity[ielz  ][iely  ][ielx  ];
				rho_mean 		=  	Local_Density[ielz  ][iely  ][ielx  ];
				num_particles 	=	Local_NumParticles[ielz  ][iely  ][ielx  ];

				/* Compute dominant phase in the cell */
				max_phase_fraction=0;	phase_dominant=0;
				for (iphase=0; iphase<user->num_phases; iphase++){
					VolumeFraction = PhaseProportions_local[iphase][ielz][iely][ielx];
					if (VolumeFraction>max_phase_fraction){
						max_phase_fraction	=	VolumeFraction;
						phase_dominant 		=	iphase;
					}
				}

				for (intp=0; intp<ngp_vel; intp++){  // loop over integration points

					material_fine_data.Viscosity[intp]		=	eta_mean;
					material_fine_data.Density[intp]		=	rho_mean;
					material_fine_data.TrueViscosity[intp]	=	eta_mean;
					material_fine_data.Phases[intp]			=	(PetscScalar) phase_dominant;

					material_fine_data.NumParticles[intp]	=	num_particles;


				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials,   user->Materials, &materials_fine_array); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_EffectiveViscosity , 	&Local_Viscosity	 ); 		CHKERRQ(ierr); // Viscosity @ center
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_Density , 			&Local_Density	 ); 			CHKERRQ(ierr); // Density @ center
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,  user->FDSTAG.Center_NumParticles , 		&Local_NumParticles	 	); 		CHKERRQ(ierr); // # of particles @ center

	for (i=0; i<user->num_phases; i++ ){
		ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,   	user->FDSTAG.Center_PhaseProportions[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr);
	}
	/* END OF PART 3: UPDATE PROPERTIES ON INTEGRATION POINTS ----------------------------------------------*/


	/* PART 4: UPDATE LOCAL VECTOR THAT CONTAINS PHASE INFO FOR CORNER POINTS */
	// This is used later to determine the dominant phase @ corners [in FDSTAG_DominantPhaseAtCorners].
	for (i=0; i<user->num_phases; i++ ){
		ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_PhaseProportions[i],INSERT_VALUES,user->FDSTAG.Corner_PhaseProportions_local[i]); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_PhaseProportions[i],INSERT_VALUES,user->FDSTAG.Corner_PhaseProportions_local[i]); CHKERRQ(ierr);
	}


	ierr = PetscFree(VolumeFraction_phases); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* This routine computes the dominant phase from the surrounding 8 corner nodes.
 * It is used
 * 	(1) In cases where the element is empty and we need to inject new particles.
 * 	(2) In cases where a 'sticky-air' element is subducted into the mantle and needs to be transformed into
 * 		a rock-phase.  In this case we determine the dominant non-air phase.
 *
 * 		The routine returns -1 as phase if it cannot detect a valid phase in the corner nodes either.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_DominantPhaseAtCorners"
PetscErrorCode FDSTAG_DominantPhaseAtCorners(UserContext *user, PetscInt ielx, PetscInt iely, PetscInt ielz,
		PetscInt *DominantPhase_WithAir, PetscInt *DominantPhase_WithoutAir, PetscInt *DominantPhase_WithoutAir_TopCell)
{
	PetscErrorCode 			ierr;
	PetscInt				i, xs,ys,zs,xm,ym,zm, iphase, StickyAirPhase, DomPhase, PhasesToBeExcluded[100], numPhasesToBeExcuded, num_exclude;
	PetscScalar 			***PhaseProportions_local[max_num_phases], PhaseFractionsWithAir[max_num_phases], PhaseFractionsWithoutAir[max_num_phases], PhaseFractionsWithoutAirTopCell[max_num_phases];
	PetscScalar				DominantFraction;
	PetscBool				found;

	/* Make a loop over every FD node */
	ierr = DMDAGetGhostCorners(user->FDSTAG.DA_CORNER,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	/* Get arrays of properties at the corners */
	for (i=0; i<user->num_phases; i++){

		/* We need to use a local vector WITH ghost points for this, that is created outside this subroutine [in UpdateMaterialProperties_FDSTAG]. */
		ierr = DMDAVecGetArray(user->FDSTAG.DA_CORNER,  user->FDSTAG.Corner_PhaseProportions_local[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr); // Phase proportions
	}


	/* Catch errors */
	if ((ielx<xs) || (ielx>=(xs+xm)) || (iely<ys) || (iely>=(ys+ym)) || (ielz<zs) || (ielz>=(zs+zm))){
		PetscPrintf( PETSC_COMM_SELF, "FDSTAG_DominantPhaseAtCorners: ielx=%i iely=%i ielz=%i xs=%i ys=%i zs=%i xs+xm=%i, ys+ym=%i, zs+zm=%i \n",ielx,iely, ielz,xs,ys,zs,xs+xm,ys+ym,zs+zm);


		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_SUP, "FDSTAG_DominantPhaseAtCorners: trying to find phase for the element that is not on the current processor, ielx=%i iely=%i ielz=%i",ielx,iely, ielz);
	}


	/* Sum the phase fractions of all phases of the surrounding 8 nodes */


	// Define which phases you want to exclude from computing the dominant phase in a cell (e.g., Air, Salt etc.)
	StickyAirPhase			= user->ErosionParameters.StickyAirPhase;

	// define additional phases on the command line that will be excluded while computing dominant phases
	numPhasesToBeExcuded = 100;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-PhasesToBeExcludedInDominantPhaseCalculation", PhasesToBeExcluded, &numPhasesToBeExcuded, &found); CHKERRQ(ierr);
	if (found){
	//	PetscPrintf(PETSC_COMM_WORLD," Excluding %i phases from Dominant Phase Calculation; StickyAirPhase=%i \n",numPhasesToBeExcuded,StickyAirPhase);
	}
	else{
		numPhasesToBeExcuded = 0;
	}
	PhasesToBeExcluded[numPhasesToBeExcuded] 	= StickyAirPhase;				// always add sticky air phase
	numPhasesToBeExcuded 						= numPhasesToBeExcuded+1;


	/* Initialize to zero */
	for (iphase=0; iphase<user->num_phases; iphase++){
		PhaseFractionsWithAir[iphase] 			= 0.0;
		PhaseFractionsWithoutAir[iphase] 		= 0.0;
		PhaseFractionsWithoutAirTopCell[iphase] = 0.0;
	}

	/* Compute fractions @ every corner */
	for (iphase=0; iphase<user->num_phases; iphase++){
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz  ][iely  ][ielx  ];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz  ][iely+1][ielx  ];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz  ][iely+1][ielx+1];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz  ][iely  ][ielx+1];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz+1][iely  ][ielx  ];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz+1][iely+1][ielx  ];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz+1][iely+1][ielx+1];
		PhaseFractionsWithAir[iphase] = PhaseFractionsWithAir[iphase] +	PhaseProportions_local[iphase][ielz+1][iely  ][ielx+1];

		// compute phase fraction @ top of cell
		PhaseFractionsWithoutAirTopCell[iphase] = PhaseFractionsWithoutAirTopCell[iphase] +	PhaseProportions_local[iphase][ielz+1][iely  ][ielx  ];
		PhaseFractionsWithoutAirTopCell[iphase] = PhaseFractionsWithoutAirTopCell[iphase] +	PhaseProportions_local[iphase][ielz+1][iely+1][ielx  ];
		PhaseFractionsWithoutAirTopCell[iphase] = PhaseFractionsWithoutAirTopCell[iphase] +	PhaseProportions_local[iphase][ielz+1][iely+1][ielx+1];
		PhaseFractionsWithoutAirTopCell[iphase] = PhaseFractionsWithoutAirTopCell[iphase] +	PhaseProportions_local[iphase][ielz+1][iely  ][ielx+1];

		//PetscPrintf(PETSC_COMM_SELF,"PhaseProportions_local[iphase=%i][ielz  ][iely  ][ielx  ]=%f \n",iphase,PhaseProportions_local[iphase][ielz  ][iely  ][ielx  ]);

		PhaseFractionsWithoutAir[iphase] = PhaseFractionsWithAir[iphase];

		// exclude some phase from computing the dominant phase
		for (num_exclude=0; num_exclude<numPhasesToBeExcuded; num_exclude++){
			if (iphase==PhasesToBeExcluded[num_exclude]){
				PhaseFractionsWithoutAir[iphase] 		=  	0.0;
				PhaseFractionsWithoutAirTopCell[iphase] = 	0.0;
			}
		}


	}

	/* Use this to find the most dominant fraction with air*/
	DominantFraction = 0; DomPhase 	 = -1;
	for (iphase=0; iphase<user->num_phases; iphase++){
		if (PhaseFractionsWithAir[iphase] > DominantFraction){
			DominantFraction 	= 	PhaseFractionsWithAir[iphase];
			DomPhase 			=	iphase;
		}
	}
	*DominantPhase_WithAir = DomPhase;

	/* Use this to find the most dominant fraction without air*/
	DominantFraction = 0; DomPhase 	 = -1;
	for (iphase=0; iphase<user->num_phases; iphase++){
		if (PhaseFractionsWithoutAir[iphase] > DominantFraction){
			DominantFraction 	= 	PhaseFractionsWithoutAir[iphase];
			DomPhase 			=	iphase;
		}
	}
	*DominantPhase_WithoutAir = DomPhase;


	/* Use this to find the most dominant fraction without air at the top of the cell*/
	DominantFraction = 0; DomPhase 	 = -1;
	for (iphase=0; iphase<user->num_phases; iphase++){
		if (PhaseFractionsWithoutAirTopCell[iphase] > DominantFraction){
			DominantFraction 	= 	PhaseFractionsWithoutAirTopCell[iphase];
			DomPhase 			=	iphase;
		}
	}
	*DominantPhase_WithoutAir_TopCell = DomPhase;


	/* Cleaning up */
	for (i=0; i<user->num_phases; i++ ){
		ierr = DMDAVecRestoreArray	(user->FDSTAG.DA_CORNER,  user->FDSTAG.Corner_PhaseProportions_local[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr); // Phase proportions
	}



	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
/* Compute properties such as velocity, strainrates, pressures etc.  at a local point, with given element coordinates. */
#undef __FUNCT__
#define __FUNCT__ "ComputePointWiseProperties"
PetscErrorCode ComputePointWiseProperties(
		LaMEMVelPressureDA C,
		PointWiseInformation *PointInformation, double Point[], double V_element[],
		double P_element[], double T_element[], double T_old_element[],
		double DevStress[], DMDACoor3d coord_elem[], double mu, double mu_viscous,
		PetscInt ComputeFull )
{
	PetscErrorCode 			ierr;
	DMDACoor3d		 		CoordPoint;
	PetscInt				i,j, npres;
	PetscScalar	 			ShapeVel[MAX_nnel],  **dhdsVel, ShapeP[MAX_npres];
	PetscScalar	 			**Jacob, **InvJacob, DetJacob, KinMtx[6][MAX_edof], **dhdPhys;
	PetscScalar	 			MatMtx[6][6], P, Temperature, Temperature_diff;
	PetscScalar	 			DevStrainRate[6], IncDevStress[6], FacInv[6], E2nd, T2nd;
	PetscInt edof, nnel;
	DAVPElementType element_type;


	nnel = C->nnel;
	npres = C->npres;
	edof = C->edof;
	element_type = C->type;

	ierr = LaMEMCreate2dArray( 3,nnel, &dhdPhys, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMCreate2dArray( 3,nnel, &dhdsVel, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMCreate2dArray( 3,3, &Jacob, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMCreate2dArray( 3,3, &InvJacob, PETSC_NULL ); CHKERRQ(ierr);


	/* Compute jacobian, shape functions and derivatives at the local particle coordinate */
	ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 		// Velocity shape function
	ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordPoint);				 	// Real coordinate of point

	if (element_type!=DAVP_FDSTAG){
		ComputeShapeFunctionPressure(ShapeP, 	CoordPoint,  Point);		 		// Pressure shape function
	}
	ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  	// Jacobian etc. of element
	ComputeKinmtx( nnel, edof, KinMtx, dhdPhys, dhdsVel, InvJacob);				// Kinematic matrix

	ComputeMaterialMatrix(mu, MatMtx);									 		// Material matrix


	if (DetJacob<0){

		PetscPrintf(PETSC_COMM_WORLD,"In routine ComputePointWiseProperties: Negative Jacobian on DetJacob=%g  \n",DetJacob);
		if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) ){
			PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
			PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
			PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
			PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x,
					coord_elem[8].x,coord_elem[9].x,coord_elem[10].x,coord_elem[11].x,coord_elem[12].x,coord_elem[13].x,coord_elem[14].x,coord_elem[15].x);

		}
		else {
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented" );
		}

		MPI_Abort(PETSC_COMM_WORLD,1);
	}



	// Compute strainrates @ integration point
	for( i = 0; i < 6; i++){ DevStrainRate[i] = 0.0; }
	for( i = 0; i < 6; i++){
		for( j = 0; j < edof; j++){
			DevStrainRate[i] = DevStrainRate[i] +  KinMtx[i][j]*V_element[j];
		}
	}
	for( i = 0; i < 6; i++){PointInformation->DeviatoricStrainRate[i] = DevStrainRate[i]; }


	// Compute increment in deviatoric stresses @ integration point
	for( i = 0; i < 6; i++){ IncDevStress[i] = 0.0; }
	for( i = 0; i < 6; i++){ for( j = 0; j < 6; j++){
		IncDevStress[i] = IncDevStress[i] +  MatMtx[i][j]*DevStrainRate[j];
	}}

	// Add increment in deviatoric stress to deviatoric stress
	for( i = 0; i < 6; i++){ DevStress[i] = 0.0; }
	for( i = 0; i < 6; i++){
		DevStress[i] = 0*DevStress[i] + IncDevStress[i];
	}
	for( i = 0; i < 6; i++){PointInformation->DeviatoricStress[i] = DevStress[i]; }


	// Compute invariants
	FacInv[0]=0.5; 	FacInv[1]=0.5;	FacInv[2]=0.5;
	FacInv[3]=1.0; 	FacInv[4]=1.0;	FacInv[5]=1.0;
	T2nd = 0.0; E2nd = 0.0;
	for (i=0; i<6; i++){
		T2nd = T2nd + FacInv[i]*DevStress[i]*DevStress[i];
		E2nd = E2nd + FacInv[i]*DevStrainRate[i]*DevStrainRate[i];
	}



	T2nd = PetscSqrtScalar(T2nd);
	E2nd = PetscSqrtScalar(E2nd);
	PointInformation->SecondInvariantDeviatoricStrainrate = E2nd;
	PointInformation->SecondInvariantDeviatoricStress     = T2nd;

	// Compute pressure
	if( element_type == DAVP_Q1P0 ) {
		P =  P_element[0];
	}
	else if( element_type == DAVP_Q2PM1L ) {
		P =  P_element[0] + P_element[1]*Point[0] + P_element[2]*Point[1] + P_element[3]*Point[2];
	}
	else if( element_type == DAVP_Q2PM1G ) {
		P =  P_element[0] + P_element[1]*CoordPoint.x + P_element[2]*CoordPoint.y + P_element[3]*CoordPoint.z;
	}
	else if (element_type==DAVP_Q1Q1){
		P = 0.0;
		for (i=0; i<npres; i++){
			P = P +  P_element[i]*ShapeP[i];
		}
	}
	else if ( element_type == DAVP_FDSTAG ) {
		P =  P_element[0];
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented" );
	}
	PointInformation->Pressure  = P;

	// Compute temperature
	Temperature = Temperature_diff = 0;
	for (i=0; i<nnel; i++){
		// this assumes, ofcourse, that the temperature shape fct is the same as velocity
		Temperature 	= Temperature 		+ 	ShapeVel[i]*T_element[i];							// temperature
		if (T_old_element != PETSC_NULL){
			Temperature_diff= Temperature_diff 	+  	ShapeVel[i]*(T_element[i]-T_old_element[i]);	// temperature difference
		}
		else{
			Temperature_diff=0;
		}

	}
	PointInformation->Temperature 		=	Temperature;
	PointInformation->Temperature_diff 	=	Temperature_diff;


	/* ==========================================================================================================
	 * The remaining part of the routine should only be done once every time-step (data is not required during iterations
	 */
	if (ComputeFull==1){
		/*  Compute shear heating  */
		PointInformation->ShearHeat = 0.0;

		// plastic strain rate
		PointInformation->PlasticStrainrate2ndInvariant = E2nd - T2nd/(2.0*mu_viscous);		// E=E_vis + E_pl  -> E_pl = E-E_vis = E-T2nd/2/mu_vis

		PointInformation->x         		= 	CoordPoint.x;
		PointInformation->y         		= 	CoordPoint.y;
		PointInformation->z         		= 	CoordPoint.z;

		PointInformation->ElementVolume		=	DetJacob;


	}
	/*
	 * ==========================================================================================================
	 */



	ierr = LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMDestroy2dArray(&Jacob, PETSC_NULL ); CHKERRQ(ierr);
	ierr = LaMEMDestroy2dArray(&InvJacob, PETSC_NULL ); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


