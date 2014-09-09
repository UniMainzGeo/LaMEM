//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "scaling.h"
//---------------------------------------------------------------------------
void ComputeScaling(
	Scaling    * scal,
	PetscScalar  mass,
	PetscScalar  time,
	PetscScalar  length,
	PetscScalar  temperature,
	PetscScalar  force)
{
	// primary characteristic units
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

	scal->phase               = 1.0;
}
//---------------------------------------------------------------------------
PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n)
{
	// power law constant scaling: 1 / stress^n / time

	return 1.0/pow(scal->stress, n)/scal->time;
}
//---------------------------------------------------------------------------


