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

NonDimensionalisation.c, contains the following subroutines:

ComputeCharacteristicValues		-	Computes characteristic values used for nondimensionalisation
PerformNonDimensionalization	-	Scale variables by their characteristic values

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "NonDimensionalisation.h"


/*==========================================================================================================*/
void ComputeCharacteristicValues( UserContext *user )
{
	NonDimUnits Characteristic_;

	memset( (void*)&Characteristic_, 0, sizeof(NonDimUnits) );

	if (user->DimensionalUnits==0){
		/* No ND performed */
		Characteristic_.Length 		= 1.0;
		Characteristic_.Stress 		= 1.0;
		Characteristic_.Temperature	= 1.0;
		Characteristic_.Viscosity	= 1.0;

		PetscPrintf(PETSC_COMM_WORLD," Input units                    : Non-dimensional\n");

		/* Some useful values */
		Characteristic_.km          = 1.0;
		Characteristic_.SecYear     = 1.0;
		Characteristic_.cmYear      = 1.0;
		Characteristic_.Myrs        = 1.0;
		Characteristic_.MPa         = 1.0;
		Characteristic_.Jmol		= 1.0;
	}

	else if (user->DimensionalUnits==1){
		/* Dimensionalization is performed
		* in code: input parameters are in SI units
		*
		*/

		/* in most cases, these values are specified in the input file */
		if (user->InputParamFile){
			Characteristic_.Length 		= user->Characteristic.Length;			// m
			Characteristic_.Stress 		= user->Characteristic.Stress;			// Pa
			Characteristic_.Temperature	= user->Characteristic.Temperature;     // Celcius !!
			Characteristic_.Viscosity	= user->Characteristic.Viscosity;		// Pa.s
		}
		else{
			/* if we did not read an input file, set the values here */
			Characteristic_.Length 		= 10e3;		// m
			Characteristic_.Stress 		= 100e6;    // Pa
			Characteristic_.Temperature	= 1000;     // Celcius !!
			Characteristic_.Viscosity	= 1e21;		// Pa.s
		}

        PetscPrintf(PETSC_COMM_WORLD," Input units                    : Dimensional \n");
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Length 	= %g [m] \n",Characteristic_.Length);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Stress 	= %g [Pa] \n",Characteristic_.Stress);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Temperature 	= %g  [Celcius/K] \n",Characteristic_.Temperature);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Viscosity 	= %g [Pa.s] \n",Characteristic_.Viscosity);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Time 	= %g [s] \n",Characteristic_.Viscosity/Characteristic_.Stress);

		/* Some useful values */
		Characteristic_.km          = 1000;
		Characteristic_.SecYear     = 365.25*24*3600;
		Characteristic_.cmYear      = Characteristic_.SecYear*100;  // from m/s to cm/year
		Characteristic_.Myrs        = 1e6;
		Characteristic_.MPa         = 1e6;
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD," Something wrong in defining ND units! \n");
	}

	/* Derived values from the Characteristic_ values above */
	Characteristic_.Time   				= 	Characteristic_.Viscosity/Characteristic_.Stress;
	Characteristic_.Velocity				= 	Characteristic_.Length/Characteristic_.Time;
	Characteristic_.kg					= 	Characteristic_.Stress*Characteristic_.Length*(Characteristic_.Time*Characteristic_.Time);
	Characteristic_.Density				= 	Characteristic_.kg/(Characteristic_.Length*Characteristic_.Length*Characteristic_.Length);
	Characteristic_.Strainrate			=	1.0/Characteristic_.Time;
	Characteristic_.ThermalExpansivity 	=  	1.0/Characteristic_.Temperature;
	Characteristic_.Force			  	=  	Characteristic_.Stress*Characteristic_.Length*Characteristic_.Length;
	Characteristic_.Watt				  	=  	Characteristic_.Force*Characteristic_.Length/Characteristic_.Time;
	Characteristic_.Joule			  	=  	Characteristic_.Force*Characteristic_.Length;
	Characteristic_.HeatCapacity			=   Characteristic_.Joule/Characteristic_.kg/Characteristic_.Temperature;
	Characteristic_.T_conductivity	  	=  	Characteristic_.Watt/Characteristic_.Length/Characteristic_.Temperature;
	Characteristic_.RadioactiveHeat    	=  	Characteristic_.Watt/Characteristic_.Length/Characteristic_.Length/Characteristic_.Length;
	Characteristic_.Jmol					=	user->GasConstant*Characteristic_.Temperature;

	user->Characteristic = 	Characteristic_;
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Transforms input parameters into ND parameters */
void PerformNonDimensionalization( UserContext *user )
{
	PetscInt	i, iphase;
	PetscScalar SecYear;


	SecYear = 3600.0*24.0*365.25;

	user->W 	 = (user->W)/(user->Characteristic.Length);user->x_left = (user->x_left)/(user->Characteristic.Length);
	user->L 	 = user->L/user->Characteristic.Length;			user->y_front= user->y_front/user->Characteristic.Length;
	user->H 	 = user->H/user->Characteristic.Length;			user->z_bot  = user->z_bot/user->Characteristic.Length;
	user->ampl2D = user->ampl2D/user->Characteristic.Length;	user->ampl3D = user->ampl3D/user->Characteristic.Length;
	user->amplNoise = user->amplNoise/user->Characteristic.Length;
	user->Setup.Diapir_Hi = user->Setup.Diapir_Hi/user->Characteristic.Length;
	user->SurfaceNoiseAmplitude = user->SurfaceNoiseAmplitude/user->Characteristic.Length;

	user->BC.Vy_front 	= 	user->BC.Vy_front/user->Characteristic.Velocity;
	user->BC.Vy_back 	= 	user->BC.Vy_back/user->Characteristic.Velocity;
	user->BC.Vz_top 	= 	user->BC.Vz_top/user->Characteristic.Velocity;
	user->BC.Vz_bot 	= 	user->BC.Vz_bot/user->Characteristic.Velocity;
	user->BC.Vx_left 	= 	user->BC.Vx_left/user->Characteristic.Velocity;
	user->BC.Vx_right 	= 	user->BC.Vx_right/user->Characteristic.Velocity;
	user->BC.Exx 		= 	user->BC.Exx/user->Characteristic.Strainrate;
	user->BC.Eyy 		= 	user->BC.Eyy/user->Characteristic.Strainrate;
	user->time          = 	user->time/user->Characteristic.Time;
	user->dt            =   user->dt*user->Characteristic.SecYear/user->Characteristic.Time;
	user->dt_max        =   user->dt_max*user->Characteristic.SecYear/user->Characteristic.Time;
	user->dt_temp       =   user->dt_temp*user->Characteristic.SecYear/user->Characteristic.Time;
	user->Temp_top		=	user->Temp_top/user->Characteristic.Temperature;
	user->Temp_bottom	=	user->Temp_bottom/user->Characteristic.Temperature;

	// internal free surface & erosion on it
	user->ErosionParameters.InitialFreeSurfaceHeight 	= user->ErosionParameters.InitialFreeSurfaceHeight/user->Characteristic.Length;
	user->ErosionParameters.SedimentationRate 			= (user->ErosionParameters.SedimentationRate_cmYr*0.01/SecYear)/user->Characteristic.Velocity;
	user->ErosionParameters.SedimentLayerThicknessYears = (user->ErosionParameters.SedimentLayerThicknessYears)*SecYear/user->Characteristic.Time;

	user->Gravity 		=	user->Gravity/( user->Characteristic.Length/(user->Characteristic.Time*user->Characteristic.Time));

	/* Transform BCs in m/s and into ND units */
	user->Vx_Front = user->Vx_Front*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;
	user->Vx_Back  = user->Vx_Back*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;
	user->Vy_Front = user->Vy_Front*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;
	user->Vy_Back  = user->Vy_Back*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;


	/* transform PT-related parameters into ND units */
	for (iphase=0; iphase<user->num_phase_transitions; iphase++){
		user->PhaseTransitions[iphase].TransitionDepth 	= user->PhaseTransitions[iphase].TransitionDepth/user->Characteristic.Length;
		user->PhaseTransitions[iphase].TransitionP0 	= user->PhaseTransitions[iphase].TransitionP0/user->Characteristic.Stress;
		user->PhaseTransitions[iphase].TransitionAlpha 	= user->PhaseTransitions[iphase].TransitionAlpha/(user->Characteristic.Stress/user->Characteristic.Temperature);
	}

	user->LowerViscosityCutoff = user->LowerViscosityCutoff/user->Characteristic.Viscosity;
	user->UpperViscosityCutoff = user->UpperViscosityCutoff/user->Characteristic.Viscosity;

	/* Nondimensionalize material properties */
	for (i=0;i<user->num_phases; i++){

			user->PhaseProperties.mu[i]					=	user->PhaseProperties.mu[i] /user->Characteristic.Viscosity;
			user->PhaseProperties.FrankKamenetskii[i]   =	user->PhaseProperties.FrankKamenetskii[i]/(1.0/user->Characteristic.Temperature);
			user->PhaseProperties.Powerlaw_e0[i] 		=	user->PhaseProperties.Powerlaw_e0[i]/(user->Characteristic.Strainrate);

// I don't understand why "A" is interpreted as real pre-exponential factor of the stress-dependent dislocation creep.
// Why then all rheologies that use "A" variable, also define the reference strain rate "Powerlaw_e0"?
// What if I want a standard temperature-dependent power-low rheology?
// I reinterpret this as a reference viscosity similar to "mu", since no one else is using it (a.p.)
//			user->PhaseProperties.A[i]					=	user->PhaseProperties.A[i]/((pow(user->Characteristic.Stress,(-user->PhaseProperties.n_exponent[i])))/user->Characteristic.Time);
			user->PhaseProperties.A[i]					=	user->PhaseProperties.A[i]/user->Characteristic.Viscosity;

			user->PhaseProperties.E[i]					=	user->PhaseProperties.E[i]/user->Characteristic.Jmol;

			user->PhaseProperties.ElasticShearModule[i]		=	user->PhaseProperties.ElasticShearModule[i]/user->Characteristic.Stress;
			user->PhaseProperties.ElasticBulkModule[i]		=	user->PhaseProperties.ElasticBulkModule[i]/user->Characteristic.Stress;

			user->PhaseProperties.Cohesion[i]				=	user->PhaseProperties.Cohesion[i]/user->Characteristic.Stress;
			user->PhaseProperties.CohesionAfterWeakening[i]	=	user->PhaseProperties.CohesionAfterWeakening[i]/user->Characteristic.Stress;

			user->PhaseProperties.T_Conductivity[i]		=	user->PhaseProperties.T_Conductivity[i] /user->Characteristic.T_conductivity;
			user->PhaseProperties.HeatCapacity[i]		=	user->PhaseProperties.HeatCapacity[i]	/user->Characteristic.HeatCapacity;
			user->PhaseProperties.RadioactiveHeat[i]	=	user->PhaseProperties.RadioactiveHeat[i]/user->Characteristic.RadioactiveHeat;

			user->PhaseProperties.rho[i]				=	user->PhaseProperties.rho[i]/user->Characteristic.Density;
			user->PhaseProperties.Density_T0[i] 		=	user->PhaseProperties.Density_T0[i]/user->Characteristic.Temperature;
			user->PhaseProperties.ThermalExpansivity[i]	=	user->PhaseProperties.ThermalExpansivity[i]/(1.0/user->Characteristic.Temperature);

	}
	/* Nondimensionalize pushing block properties */
	user->Pushing.L_block 	 = user->Pushing.L_block/user->Characteristic.Length;
	user->Pushing.W_block 	 = user->Pushing.W_block/user->Characteristic.Length;
	user->Pushing.H_block 	 = user->Pushing.H_block/user->Characteristic.Length;
	user->Pushing.x_center_block 	 = user->Pushing.x_center_block/user->Characteristic.Length;
	user->Pushing.y_center_block 	 = user->Pushing.y_center_block/user->Characteristic.Length;
	user->Pushing.z_center_block 	 = user->Pushing.z_center_block/user->Characteristic.Length;
	for (i=0;i<user->Pushing.num_changes; i++){
		user->Pushing.V_push[i] = user->Pushing.V_push[i]*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;
		user->Pushing.omega[i]  = user->Pushing.omega[i]*(M_PI/180.0)/(user->Characteristic.SecYear/user->Characteristic.Time); //[radians]
	}
	for (i=0;i<user->Pushing.num_changes+1; i++){
		user->Pushing.time[i]	= user->Pushing.time[i]*1e6*user->Characteristic.SecYear/user->Characteristic.Time; //from Myr
	}

}
/*==========================================================================================================*/
