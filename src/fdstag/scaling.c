//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "scaling.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ScalingCreate"
PetscErrorCode ScalingCreate(
	Scaling     *scal,
	PetscInt     DimensionalUnits,
	PetscScalar  mass,
	PetscScalar  time,
	PetscScalar  length,
	PetscScalar  temperature,
	PetscScalar  force)
{

	PetscScalar yr = 1.0, Myr = 1.0, km = 1.0, cm = 1.0, cm_yr = 1.0, MPa = 1.0;

	if(!DimensionalUnits)
	{
		// clear primary units
		mass        =  1.0;
		time        =  1.0;
		length      =  1.0;
		temperature =  1.0;
		force       =  1.0;
	}
	else
	{
		yr    = 3600.0*24.0*365.0;
		Myr   = 1e6*yr;
		km    = 1e3;
		cm    = 1e-2;
		cm_yr = cm/yr;
		MPa   = 1e6;
	}

	// set primary characteristic units
	scal->mass        = mass;
	scal->time        = time;
	scal->length      = length;
	scal->temperature = temperature;
	scal->force       = force;

	// additional units
	scal->volume = length*length*length;
	scal->area   = length*length;

	// secondary units
	scal->velocity         = length/time;
	scal->stress           = force/scal->area;
	scal->strain_rate      = 1.0/time;
	scal->gravity_strength = force/mass;
	scal->energy           = force*length;
	scal->power            = scal->energy/time;
	scal->heat_flux        = scal->power/scal->area;
	scal->dissipation_rate = scal->power/scal->volume;

	// material parameters
	scal->density             = mass/scal->volume;
	scal->viscosity           = scal->stress*time;
	scal->cpecific_heat       = scal->energy/mass/temperature;
	scal->conductivity        = scal->power/length/temperature;
	scal->heat_production     = scal->power/mass;
	scal->expansivity         = 1.0/temperature;
	scal->pressure_sensivity  = temperature/scal->stress;

	scal->out_time     = time/Myr;
	scal->out_length   = length/km;
	scal->out_velocity = scal->velocity/cm_yr;
	scal->out_stress   = scal->stress/MPa;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n)
{
	// power law constant scaling: 1 / stress^n / time

	return 1.0/pow(scal->stress, n)/scal->time;
}
//---------------------------------------------------------------------------


