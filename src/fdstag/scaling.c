//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "scaling.h"
#include "Parsing.h"
//---------------------------------------------------------------------------
// * replace scaling consistently (at the input level)
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ScalingReadFromFile"
PetscErrorCode ScalingReadFromFile(Scaling *scal, FILE *fp)
{
	PetscInt flg;
	char     utype [MAX_NAME_LEN];

	PetscFunctionBegin;

	// read model setup
	parse_GetString(fp, "units", utype, MAX_NAME_LEN-1, &flg);

	if(flg)
	{
		if     (!strcmp(utype, "none"))   scal->utype = _NONE_;
		else if(!strcmp(utype, "si"))     scal->utype = _SI_;
		else if(!strcmp(utype, "geo"))    scal->utype = _GEO_;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect type of units: %s", utype);
	}

	PetscFunctionReturn(0);
}
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
	// characteristic values must ALWAYS be given in SI units

	PetscScalar angle, area, volume, stress, energy, power;
	PetscScalar yr, Myr, km, cm, cm_yr, MPa, mW;

	if(DimensionalUnits && scal->utype ==_NONE_)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "set 'DimensionalUnits' & 'units' options coherently in the input file");
	}

	// unit
	scal->unit = 1.0; sprintf(scal->lbl_unit, "[ ]");

	if(scal->utype == _NONE_)
	{
		//================
		// NON-DIMENSIONAL
		//================

		// temperature shift
		scal->Tshift = 0.0;

		// primary units
		scal->mass                = 1.0;
		scal->time                = 1.0;   sprintf(scal->lbl_time,             "[ ]");
		scal->length              = 1.0;   sprintf(scal->lbl_length,           "[ ]");
		scal->temperature         = 1.0;   sprintf(scal->lbl_temperature,      "[ ]");
		scal->force               = 1.0;   sprintf(scal->lbl_force,            "[ ]");
		scal->angle               = 1.0;   sprintf(scal->lbl_angle,            "[ ]");

		// secondary units
		scal->velocity            = 1.0;   sprintf(scal->lbl_velocity,         "[ ]");
		scal->stress              = 1.0;   sprintf(scal->lbl_stress,           "[ ]");
		scal->strain_rate         = 1.0;   sprintf(scal->lbl_strain_rate,      "[ ]");
		scal->gravity_strength    = 1.0;
		scal->energy              = 1.0;
		scal->power               = 1.0;
		scal->heat_flux           = 1.0;   sprintf(scal->lbl_heat_flux,        "[ ]");
		scal->dissipation_rate    = 1.0;   sprintf(scal->lbl_dissipation_rate, "[ ]");
		scal->angular_velocity    = 1.0;   sprintf(scal->lbl_angular_velocity, "[ ]");

		// material parameters
		scal->density             = 1.0;   sprintf(scal->lbl_density,          "[ ]");
		scal->viscosity           = 1.0;   sprintf(scal->lbl_viscosity,        "[ ]");
		scal->cpecific_heat       = 1.0;
		scal->conductivity        = 1.0;
		scal->heat_production     = 1.0;
		scal->expansivity         = 1.0;
		scal->pressure_sensivity  = 1.0;

	}
	else if(scal->utype == _SI_)
	{
		//=========================
		// SI UNITS (except angles)
		//=========================

		// additional units
		angle  = 180.0/M_PI;
		area   = length*length;
		volume = area*length;
		stress = force/area;
		energy = force*length;
		power  = energy/time;

		// temperature shift
		scal->Tshift = 0.0;

		// primary units
		scal->mass                = mass;
		scal->time                = time;                     sprintf(scal->lbl_time,             "[s]");
		scal->length              = length;                   sprintf(scal->lbl_length  ,         "[m]");
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[K]");
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @

		// secondary units
		scal->velocity            = length/time;              sprintf(scal->lbl_velocity,         "[m/s]");
		scal->stress              = stress;                   sprintf(scal->lbl_stress ,          "[Pa]");
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->gravity_strength    = force/mass;
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area;               sprintf(scal->lbl_heat_flux,        "[W/m^2]");
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/time;               printf(scal->lbl_angular_velocity,  "[deg/s]"); // @

		// material parameters
		scal->density             = mass/volume;              sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = stress*time;              sprintf(scal->lbl_viscosity,        "[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;
		scal->conductivity        = power/length/temperature;
		scal->heat_production     = power/mass;
		scal->expansivity         = 1.0/temperature;
		scal->pressure_sensivity  = temperature/stress;

	}
	else if(scal->utype == _GEO_)
	{
		//=================
		// GEOLOGICAL UNITS
		//=================

		// additional units
		angle  = 180.0/M_PI;
		area   = length*length;
		volume = area*length;
		stress = force/area;
		energy = force*length;
		power  = energy/time;

		// additional scaling factors
		yr     = 3600.0*24.0*365.0;
		Myr    = 1e6*yr;
		km     = 1e3;
		cm     = 1e-2;
		cm_yr  = cm/yr;
		MPa    = 1e6;
		mW     = 1e-3;

		// temperature shift
		scal->Tshift = 273.15;

		// primary units
		scal->mass                = mass;
		scal->time                = time/Myr;                 sprintf(scal->lbl_time,             "[Myr]");   // @
		scal->length              = length/km;                sprintf(scal->lbl_length  ,         "[km]");    // @
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[C]");     // @
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @

		// secondary units
		scal->velocity            = length/time/cm_yr;        sprintf(scal->lbl_velocity,         "[cm/yr]"); // @
		scal->stress              = stress/MPa;               sprintf(scal->lbl_stress ,          "[MPa]");   // @
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->gravity_strength    = force/mass;
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area/mW;            sprintf(scal->lbl_heat_flux,        "[mW/m^2]");  // @
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/(time/Myr);         sprintf(scal->lbl_angular_velocity, "[deg/Myr]"); // @

		// material parameters
		scal->density             = mass/volume;              sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = stress*time;              sprintf(scal->lbl_viscosity,        "[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;
		scal->conductivity        = power/length/temperature;
		scal->heat_production     = power/mass;
		scal->expansivity         = 1.0/temperature;
		scal->pressure_sensivity  = temperature/stress;

	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n)
{
	// power law constant scaling: 1 / stress^n / time

	return 1.0/pow(scal->stress, n)/scal->time;
}
//---------------------------------------------------------------------------


