//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __scaling_h__
#define __scaling_h__
//---------------------------------------------------------------------------

#define _lbl_sz_ 23

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
	// units = none - input & output is non-dimensional
	// units = si   - input & output is in SI units
	// units = geo  - input & output is in SI units, except:
	//
	//    time        - Myr
	//    length      - km
	//    velocity    - cm/yr
	//    stress      - MPa
	//    heat_flux   - mW/m^2
	//    Temperature - C
	//
	// WARNING!
	//
	// * characteristic values must ALWAYS be provided in SI units
	//
	// * in all dimensional cases (si & geo) angles are measured in degrees
	//   angular velocities are measured in degrees per unit time
	//
	// * number of primary units is one more than usual
	//   Newton's 2nd law can be violated for quasi-static problems
	//   Gravity strength must be provided in the units [force/mass]
	//=======================================================================

	UnitsType   utype;  // scaling type
	PetscScalar unit;   // always unit
	PetscScalar Tshift; // temperature shift (added on input, subtracted on output)

	// primary characteristic units
	PetscScalar mass;
	PetscScalar time;
	PetscScalar length;
	PetscScalar temperature;       // Kelvin (if dimensional)
	PetscScalar force;             // additional variable for quasi-static case
	PetscScalar angle;             // radian expressed in degrees (if dimensional)

	// secondary units
	PetscScalar velocity;          // length / time
	PetscScalar acceleration;      // length / time / time
	PetscScalar stress;            // force / area
	PetscScalar strain_rate;       // 1 / time
	PetscScalar gravity_strength;  // force / mass
	PetscScalar energy;            // force * length
	PetscScalar power;             // energy / time
	PetscScalar heat_flux;         // power / area
	PetscScalar dissipation_rate;  // power / volume
	PetscScalar angular_velocity;  // angle / time
	PetscScalar volumetric_force;  // force / volume

	// material parameters
	PetscScalar density;            // mass / volume
	PetscScalar viscosity;          // stress * time
	PetscScalar cpecific_heat;      // energy / mass / temperature
	PetscScalar conductivity;       // power / length / temperature
	PetscScalar heat_production;    // power / mass
	PetscScalar expansivity;        // 1 / temperature
	PetscScalar pressure_sensivity; // temperature / stress

	// output labels
	char lbl_unit            [_lbl_sz_];
	char lbl_angle           [_lbl_sz_];
	char lbl_time            [_lbl_sz_];
	char lbl_length          [_lbl_sz_];
	char lbl_temperature     [_lbl_sz_];
	char lbl_force           [_lbl_sz_];
	char lbl_velocity        [_lbl_sz_];
	char lbl_stress          [_lbl_sz_];
	char lbl_strain_rate     [_lbl_sz_];
	char lbl_heat_flux       [_lbl_sz_];
	char lbl_dissipation_rate[_lbl_sz_];
	char lbl_angular_velocity[_lbl_sz_];
	char lbl_volumetric_force[_lbl_sz_];
	char lbl_density         [_lbl_sz_];
	char lbl_viscosity       [_lbl_sz_];

} Scaling;
//---------------------------------------------------------------------------
// scaling routines
PetscErrorCode ScalingMain(Scaling *scal, MatParLim *matLim, Material_t  *phases, PetscInt numPhases, UserCtx *usr);

PetscErrorCode ScalingCreate(
	Scaling     *scal,
	PetscInt     DimensionalUnits,
	PetscScalar  mass,
	PetscScalar  time,
	PetscScalar  length,
	PetscScalar  temperature,
	PetscScalar  force);

PetscErrorCode ScalingReadFromFile(Scaling *scal, FILE *fp);

PetscScalar ScalingComputePowerLaw(Scaling * scal, PetscScalar n);

//---------------------------------------------------------------------------
// compute characteristic values
void ScalingCharValues(UserCtx *user);

// scaling of input parameters (UserCtx)
void ScalingInput(Scaling *scal, UserCtx *user);

// scaling material parameters
void ScalingMatProp(Scaling *scal, Material_t *phases, PetscInt numPhases);

// scaling material parameter limits
void ScalingMatParLim(Scaling *scal, MatParLim *matLim);

// OLD - scaling of input parameters (UserCtx)
void ScalingInputOLD(UserCtx *user);

// OLD - scaling material parameters
void ScalingMatPropOLD(UserCtx *user);

#endif
