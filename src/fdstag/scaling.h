//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __scaling_h__
#define __scaling_h__
//---------------------------------------------------------------------------
typedef struct
{
	// primary characteristic units
	PetscScalar mass;
	PetscScalar time;
	PetscScalar length;
	PetscScalar temperature; // Kelvin
	PetscScalar force;       // Newton's 2nd law can be violated for quasi-static problems !!!

	// additional units
	PetscScalar volume; // length^3
	PetscScalar area;   // length^2

	// secondary units
	PetscScalar velocity;         // length / time
	PetscScalar stress;           // force / area
	PetscScalar strain_rate;      // 1 / time
	PetscScalar gravity_strength; // force / mass
	PetscScalar energy;           // force * length
	PetscScalar power;            // energy / time
	PetscScalar heat_flux;        // power / area
	PetscScalar dissipation_rate; // power / volume

	// material parameters
	PetscScalar density;             // mass / volume
	PetscScalar viscosity;           // stress * time
	PetscScalar cpecific_heat;       // energy / mass / temperature
	PetscScalar conductivity;        // power / length / temperature
	PetscScalar heat_production;     // power / mass
	PetscScalar expansivity;         // 1 / temperature
	PetscScalar pressure_sensivity;  // temperature / stress

	// output scaling
	PetscScalar out_time;     // * -> s   -> Myr
	PetscScalar out_length;   // * -> m   -> km
	PetscScalar out_velocity; // * -> m/s -> cm/yr
	PetscScalar out_stress;   // * -> Pa  -> MPa
	PetscScalar out_phase;    // unit

} Scaling;
//---------------------------------------------------------------------------

PetscErrorCode ScalingCreate(
	Scaling     *scal,
	PetscInt     DimensionalUnits,
	PetscScalar  mass,
	PetscScalar  time,
	PetscScalar  length,
	PetscScalar  temperature,
	PetscScalar  force);

PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n);

//---------------------------------------------------------------------------


#endif
