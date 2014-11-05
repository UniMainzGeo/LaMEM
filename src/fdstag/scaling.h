//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __scaling_h__
#define __scaling_h__
//---------------------------------------------------------------------------

#define _lbl_sz_ 13

//---------------------------------------------------------------------------

typedef enum
{
	_NONE_, // non-dimensional
	_SI_,   // International System of units
	_GEO_   // geological-scale units

} UnitsType;

//---------------------------------------------------------------------------

typedef struct
{
	//=======================================================================
	// multiply with scale to get scaled output (normally SI units)
	// divide by scale to convert input into internal units
	//
	// units = none - input & output is non-dimensional (unit scaling is done)
	// units = si   - input & output is in SI units
	// units = geo  - input & output is in SI units, except:
	//
	//    time      - Myr
	//    length    - km
	//    velocity  - cm/yr
	//    stress    - MPa
	//    heat_flux - mW/m^2
	//
	// WARNING!
	//
	// * characteristic values must ALWAYS be provided in SI units
	//
	// * number of primary units is one more that usual
	//   Newton's 2nd law can be violated for quasi-static problems
	//   Gravity strength must be provided in the units [force/mass]
	//=======================================================================

	UnitsType   utype;

	// primary characteristic units
	PetscScalar mass;
	PetscScalar time;
	PetscScalar length;
	PetscScalar temperature; // Kelvin
	PetscScalar force;       // additional variable for quasi-static case

	// secondary units
	PetscScalar velocity;          // length / time
	PetscScalar stress;            // force / area
	PetscScalar strain_rate;       // 1 / time
	PetscScalar gravity_strength;  // force / mass
	PetscScalar energy;            // force * length
	PetscScalar power;             // energy / time
	PetscScalar heat_flux;         // power / area
	PetscScalar dissipation_rate;  // power / volume

	// material parameters
	PetscScalar density;            // mass / volume
	PetscScalar viscosity;          // stress * time
	PetscScalar cpecific_heat;      // energy / mass / temperature
	PetscScalar conductivity;       // power / length / temperature
	PetscScalar heat_production;    // power / mass
	PetscScalar expansivity;        // 1 / temperature
	PetscScalar pressure_sensivity; // temperature / stress

	PetscScalar phase;              // unit

	// output labels
	char lbl_time            [_lbl_sz_];
	char lbl_length          [_lbl_sz_];
	char lbl_temperature     [_lbl_sz_];
	char lbl_force           [_lbl_sz_];
	char lbl_velocity        [_lbl_sz_];
	char lbl_stress          [_lbl_sz_];
	char lbl_strain_rate     [_lbl_sz_];
	char lbl_heat_flux       [_lbl_sz_];
	char lbl_dissipation_rate[_lbl_sz_];
	char lbl_density         [_lbl_sz_];
	char lbl_viscosity       [_lbl_sz_];
	char lbl_phase           [_lbl_sz_];

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

PetscErrorCode ScalingReadFromFile(Scaling *scal, FILE *fp);

PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n);

//---------------------------------------------------------------------------


#endif
