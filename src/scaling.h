/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   scaling.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
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
	// * material parameters must ALWAYS be provided in SI units
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

	// input parameters
	PetscScalar inp_mass;
	PetscScalar inp_time;
	PetscScalar inp_length;
	PetscScalar inp_temperature;
	PetscScalar inp_force;

	// primary characteristic units
	PetscScalar mass;
	PetscScalar time;
	PetscScalar time_si;           // time in SI units for material parameter scaling
	PetscScalar length;
	PetscScalar length_si;         // length in SI units for material parameter scaling
	PetscScalar temperature;       // Kelvin (if dimensional)
	PetscScalar force;             // additional variable for quasi-static case
	PetscScalar angle;             // radian expressed in degrees (if dimensional)

	// secondary units
	PetscScalar velocity;          // length / time
	PetscScalar acceleration;      // length / time / time
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

PetscErrorCode ScalingCreate(Scaling *scal); //,PetscBool ExplicitSolver);

PetscErrorCode ScalingReadFromFile(Scaling *scal, FILE *fp);

//---------------------------------------------------------------------------

// scaling of input parameters (UserCtx)
void ScalingInput(Scaling *scal, UserCtx *user);

// scaling material parameters
void ScalingMatProp(Scaling *scal, Material_t *phases, PetscInt numPhases);

// scaling material parameter limits
void ScalingMatParLim(Scaling *scal, MatParLim *matLim);

void ScalingMeshSegDir(Scaling *scal, MeshSegInp *msi);

//---------------------------------------------------------------------------

#endif
