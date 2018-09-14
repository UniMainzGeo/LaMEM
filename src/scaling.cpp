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
 **    filename:   scaling.c
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
#include "LaMEM.h"
#include "scaling.h"
#include "parsing.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ScalingCreate"
PetscErrorCode ScalingCreate(Scaling *scal, FB *fb)
{
	// characteristic values must ALWAYS be given in SI units

	char        utype [_STR_LEN_];
	PetscScalar viscosity, stress, density;
	PetscScalar mass, time, length, temperature, force;
	PetscScalar angle, area, volume, energy, power;
	PetscScalar yr, Myr, km, cm, cm_yr, MPa, mW;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set unit scaling
	scal->unit = 1.0; sprintf(scal->lbl_unit, "[ ]");

	// read units type
	ierr = getStringParam(fb, _OPTIONAL_, "units", utype, "none"); CHKERRQ(ierr);

	// set units type
	if     (!strcmp(utype, "none")) scal->utype = _NONE_;
	else if(!strcmp(utype, "si"))   scal->utype = _SI_;
	else if(!strcmp(utype, "geo"))  scal->utype = _GEO_;
	else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect unit type: %s", utype);

	if(scal->utype == _NONE_)
	{
		//================
		// NON-DIMENSIONAL
		//================

		// temperature shift
		scal->Tshift = 0.0;

		// primary units
		scal->time                = 1.0;   sprintf(scal->lbl_time,             "[ ]");
		scal->time_si             = 1.0;
		scal->length              = 1.0;   sprintf(scal->lbl_length,           "[ ]");
		scal->length_si           = 1.0;
		scal->area_si             = 1.0;   sprintf(scal->lbl_area_si,          "[ ]");
		scal->temperature         = 1.0;   sprintf(scal->lbl_temperature,      "[ ]");
		scal->force               = 1.0;   sprintf(scal->lbl_force,            "[ ]");
		scal->angle               = 1.0;   sprintf(scal->lbl_angle,            "[ ]");

		// secondary units
		scal->velocity            = 1.0;   sprintf(scal->lbl_velocity,         "[ ]");
		scal->stress              = 1.0;   sprintf(scal->lbl_stress,           "[ ]");
		scal->stress_si           = 1.0;   sprintf(scal->lbl_stress_si,        "[ ]");
		scal->strain_rate         = 1.0;   sprintf(scal->lbl_strain_rate,      "[ ]");
		scal->gravity_strength    = 1.0;   sprintf(scal->lbl_gravity_strength, "[ ]");
		scal->energy              = 1.0;
		scal->power               = 1.0;
		scal->heat_flux           = 1.0;   sprintf(scal->lbl_heat_flux,        "[ ]");
		scal->dissipation_rate    = 1.0;   sprintf(scal->lbl_dissipation_rate, "[ ]");
		scal->angular_velocity    = 1.0;   sprintf(scal->lbl_angular_velocity, "[ ]");
		scal->volumetric_force    = 1.0;   sprintf(scal->lbl_volumetric_force, "[ ]");

		// material parameters
		scal->density             = 1.0;   sprintf(scal->lbl_density,          "[ ]");
		scal->viscosity           = 1.0;   sprintf(scal->lbl_viscosity,        "[ ]");
		scal->cpecific_heat       = 1.0;   sprintf(scal->lbl_cpecific_heat,    "[ ]");
		scal->conductivity        = 1.0;   sprintf(scal->lbl_conductivity,     "[ ]");
		scal->heat_production     = 1.0;   sprintf(scal->lbl_heat_production,  "[ ]");
		scal->expansivity         = 1.0;   sprintf(scal->lbl_expansivity,      "[ ]");

		sprintf(scal->lbl_diffusion_creep,   "[ ]");
		sprintf(scal->lbl_dislocation_creep, "[ ]");
		sprintf(scal->lbl_activation_energy, "[ ]");
		sprintf(scal->lbl_activation_volume, "[ ]");
		sprintf(scal->lbl_inverse_length,    "[ ]");
		sprintf(scal->lbl_inverse_stress,    "[ ]");
		sprintf(scal->lbl_gas_constant,      "[ ]");

		PetscFunctionReturn(0);
	}

	// read unit values
	temperature = 1.0;
	length      = 1.0;
	viscosity   = 1.0;
	stress      = 1.0;
	density     = 0.0;

	ierr = getScalarParam(fb, _REQUIRED_, "unit_temperature", &temperature, 1, 1.0);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "unit_length",      &length ,     1, 1.0);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "unit_viscosity",   &viscosity,   1, 1.0);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "unit_stress",      &stress,      1, 1.0);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "unit_density",     &density,     1, 1.0);  CHKERRQ(ierr);

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Scaling parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD,"   Temperature : %g [C/K] \n",    temperature);
	PetscPrintf(PETSC_COMM_WORLD,"   Length      : %g [m] \n",      length);
	PetscPrintf(PETSC_COMM_WORLD,"   Viscosity   : %g [Pa*s] \n",   viscosity);
	PetscPrintf(PETSC_COMM_WORLD,"   Stress      : %g [Pa] \n",     stress);

	if(density)
	{	PetscPrintf(PETSC_COMM_WORLD,"   Density     : %g [kg/m^3] \n", density);
		PetscPrintf(PETSC_COMM_WORLD,"   WRNING! Unconventional scaling is employed");
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// compute additional units
	angle  = 180.0/M_PI;
	area   = length*length;
	volume = area*length;
	time   = viscosity/stress;
	force  = stress*area;
	energy = force*length;
	power  = energy/time;

	// if density is not specified, enforce Second Newton's Law
	// WARINIG! this is not really necessary for quasi-static problems
	if(!density)
	{
		mass    = force/length*time*time;
		density = mass/volume;
	}
	else
	{
		mass = density*volume;
	}

	if(scal->utype == _SI_)
	{
		//=========================
		// SI UNITS (except angles)
		//=========================

		// temperature shift
		scal->Tshift = 0.0;

		// primary units
		scal->time                = time;                     sprintf(scal->lbl_time,             "[s]");
		scal->time_si             = time;
		scal->length              = length;                   sprintf(scal->lbl_length,           "[m]");
		scal->length_si           = length;
		scal->area_si             = area;                     sprintf(scal->lbl_area_si,          "[m^2]");
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[K]");
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @
		scal->volume              = volume;                   sprintf(scal->lbl_volume,        "[m^3]");

		// secondary units
		scal->velocity            = length/time;              sprintf(scal->lbl_velocity,         "[m/s]");
		scal->stress              = stress;                   sprintf(scal->lbl_stress,           "[Pa]");
		scal->stress_si           = stress;                   sprintf(scal->lbl_stress_si,        "[Pa]");
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->gravity_strength    = force/mass;               sprintf(scal->lbl_gravity_strength, "[m/s^2]");
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area;               sprintf(scal->lbl_heat_flux,        "[W/m^2]");
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/time;               sprintf(scal->lbl_angular_velocity, "[deg/s]"); // @
		scal->volumetric_force    = force/volume;             sprintf(scal->lbl_volumetric_force, "[N/m^3]");

		// material parameters
		scal->density             = density;                  sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = viscosity;                sprintf(scal->lbl_viscosity,        "[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;  sprintf(scal->lbl_cpecific_heat,    "[J/kg/K]");
		scal->conductivity        = power/length/temperature; sprintf(scal->lbl_conductivity,     "[W/m/K]");
		scal->heat_production     = power/mass;               sprintf(scal->lbl_heat_production,  "[W/m3]");
		scal->expansivity         = 1.0/temperature;          sprintf(scal->lbl_expansivity,      "[1/K]");

		sprintf(scal->lbl_diffusion_creep,   "[1/Pa/s]");
		sprintf(scal->lbl_dislocation_creep, "[1/Pa^n/s]");
		sprintf(scal->lbl_activation_energy, "[J/mol]");
		sprintf(scal->lbl_activation_volume, "[m^3/mol]");
		sprintf(scal->lbl_inverse_length,    "[1/m]");
		sprintf(scal->lbl_inverse_stress,    "[1/Pa]");
		sprintf(scal->lbl_gas_constant,      "[J/mol/K]");
	}
	else if(scal->utype == _GEO_)
	{
		//=================
		// GEOLOGICAL UNITS
		//=================

		// compute scaling factors
		yr     = 3600.0*24.0*365.25;
		Myr    = 1e6*yr;
		km     = 1e3;
		cm     = 1e-2;
		cm_yr  = cm/yr;
		MPa    = 1e6;
		mW     = 1e-3;

		// temperature shift
		scal->Tshift = 273.15;

		// primary units
		scal->time                = time/Myr;                 sprintf(scal->lbl_time,             "[Myr]");   // @
		scal->time_si             = time;
		scal->length              = length/km;                sprintf(scal->lbl_length,           "[km]");    // @
		scal->length_si           = length;
		scal->area_si             = area;                     sprintf(scal->lbl_area_si,          "[m^2]");
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[C]");     // @
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @
		scal->volume              = volume/(km*km*km);        sprintf(scal->lbl_volume,        "[km^3]");
		scal->volume_si           = volume;
		// secondary units
		scal->velocity            = length/time/cm_yr;        sprintf(scal->lbl_velocity,         "[cm/yr]"); // @
		scal->stress              = stress/MPa;               sprintf(scal->lbl_stress,           "[MPa]");   // @
		scal->stress_si           = stress;                   sprintf(scal->lbl_stress_si,        "[Pa]");
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->gravity_strength    = force/mass;               sprintf(scal->lbl_gravity_strength, "[m/s^2]");
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area/mW;            sprintf(scal->lbl_heat_flux,        "[mW/m^2]");  // @
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/(time/Myr);         sprintf(scal->lbl_angular_velocity, "[deg/Myr]"); // @
		scal->volumetric_force    = force/volume;             sprintf(scal->lbl_volumetric_force, "[N/m^3]");

		// material parameters
		scal->density             = density;                  sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = viscosity;                sprintf(scal->lbl_viscosity,        "[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;  sprintf(scal->lbl_cpecific_heat,    "[J/kg/K]");
		scal->conductivity        = power/length/temperature; sprintf(scal->lbl_conductivity,     "[W/m/K]");
		scal->heat_production     = power/mass;               sprintf(scal->lbl_heat_production,  "[W/m3]");
		scal->expansivity         = 1.0/temperature;          sprintf(scal->lbl_expansivity,      "[1/K]");

		sprintf(scal->lbl_diffusion_creep,   "[1/Pa/s]");
		sprintf(scal->lbl_dislocation_creep, "[1/Pa^n/s]");
		sprintf(scal->lbl_activation_energy, "[J/mol]");
		sprintf(scal->lbl_activation_volume, "[m^3/mol]");
		sprintf(scal->lbl_inverse_length,    "[1/m]");
		sprintf(scal->lbl_inverse_stress,    "[1/Pa]");
		sprintf(scal->lbl_gas_constant,      "[J/mol/K]");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
