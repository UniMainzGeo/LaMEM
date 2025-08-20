/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __scaling_h__
#define __scaling_h__

//---------------------------------------------------------------------------

struct FB;

//---------------------------------------------------------------------------

enum UnitsType
{
	_NONE_, // non-dimensional
	_SI_,   // International System of units
	_GEO_   // geological-scale units

};

//---------------------------------------------------------------------------

struct Scaling
{
	//=======================================================================
	// divide by scale to convert input into internal units
	// multiply with scale to get scaled output (normally SI units)
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
	// * material parameters must ALWAYS be provided in SI units
	//
	// * in all dimensional cases (si & geo) angles are measured in degrees
	//   angular velocities are measured in degrees per unit time
	//
	// * number of primary units is one more than usual
	//   Newton's 2nd law can be violated for quasi-static problems
	//
	//=======================================================================

	UnitsType   utype;  // scaling type
	PetscScalar unit;   // always unit
	PetscScalar Tshift; // temperature shift (added on input, subtracted on output)

	// primary characteristic units
	PetscScalar time;
	PetscScalar time_si;     // time in SI units for material parameter scaling
	PetscScalar length;
	PetscScalar length_si;   // length in SI units for material parameter scaling
	PetscScalar area_si;     // area   in SI units for material parameter scaling
	PetscScalar volume_si;   // volume in SI units for material parameter scaling
	PetscScalar temperature; // Kelvin (if dimensional)
	PetscScalar force;       // additional variable for quasi-static case
	PetscScalar angle;       // radian expressed in degrees (if dimensional)

	// secondary units
	PetscScalar velocity;          // length / time
	PetscScalar stress;            // force / area
	PetscScalar stress_si;         // stress in SI units for material parameter scaling
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

	// output labels
	char lbl_unit             [_lbl_sz_];
	char lbl_angle            [_lbl_sz_];
	char lbl_time             [_lbl_sz_];
	char lbl_length           [_lbl_sz_];
	char lbl_area_si          [_lbl_sz_];
	char lbl_temperature      [_lbl_sz_];
	char lbl_force            [_lbl_sz_];
	char lbl_velocity         [_lbl_sz_];
	char lbl_stress           [_lbl_sz_];
	char lbl_stress_si        [_lbl_sz_];
	char lbl_strain_rate      [_lbl_sz_];
	char lbl_gravity_strength [_lbl_sz_];
	char lbl_heat_flux        [_lbl_sz_];
	char lbl_dissipation_rate [_lbl_sz_];
	char lbl_angular_velocity [_lbl_sz_];
	char lbl_volumetric_force [_lbl_sz_];

	// material parameters labels
	char lbl_density          [_lbl_sz_];
	char lbl_viscosity        [_lbl_sz_];
	char lbl_cpecific_heat    [_lbl_sz_];
	char lbl_conductivity     [_lbl_sz_];
	char lbl_heat_production  [_lbl_sz_];
	char lbl_expansivity      [_lbl_sz_];

	// additional material parameters labels
	char lbl_diffusion_creep  [_lbl_sz_];
	char lbl_dislocation_creep[_lbl_sz_];
	char lbl_activation_energy[_lbl_sz_];
	char lbl_activation_volume[_lbl_sz_];
	char lbl_inverse_length   [_lbl_sz_];
	char lbl_inverse_stress   [_lbl_sz_];
	char lbl_gas_constant     [_lbl_sz_];

};
//---------------------------------------------------------------------------
// scaling routines

PetscErrorCode ScalingCreate(Scaling *scal, FB *fb, PetscBool PrintOutput = PETSC_FALSE);

//---------------------------------------------------------------------------
#endif
