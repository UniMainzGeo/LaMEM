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
#define __FUNCT__ "ScalingCreate"
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
	PetscScalar volume, area, yr, km, cm, cm_yr, MPa, mW_msq;

	if(DimensionalUnits && scal->utype  ==_NONE_)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "set 'DimensionalUnits' & 'units' options coherently in the input file");
	}

	if(scal->utype == _NONE_)
	{
		//================
		// NON-DIMENSIONAL
		//================
		// primary units
		scal->mass                = 1.0;
		scal->time                = 1.0;
		scal->length              = 1.0;
		scal->temperature         = 1.0;
		scal->force               = 1.0;
		// secondary units
		scal->velocity            = 1.0;
		scal->stress              = 1.0;
		scal->strain_rate         = 1.0;
		scal->gravity_strength    = 1.0;
		scal->energy              = 1.0;
		scal->power               = 1.0;
		scal->heat_flux           = 1.0;
		scal->dissipation_rate    = 1.0;
		// material parameters
		scal->density             = 1.0;
		scal->viscosity           = 1.0;
		scal->cpecific_heat       = 1.0;
		scal->conductivity        = 1.0;
		scal->heat_production     = 1.0;
		scal->expansivity         = 1.0;
		scal->pressure_sensivity  = 1.0;
		// output labels
		sprintf(scal->lbl_time,             "[ ]");
		sprintf(scal->lbl_length,           "[ ]");
		sprintf(scal->lbl_temperature,      "[ ]");
		sprintf(scal->lbl_force,            "[ ]");
		sprintf(scal->lbl_velocity,         "[ ]");
		sprintf(scal->lbl_stress,           "[ ]");
		sprintf(scal->lbl_strain_rate,      "[ ]");
		sprintf(scal->lbl_heat_flux,        "[ ]");
		sprintf(scal->lbl_dissipation_rate, "[ ]");
		sprintf(scal->lbl_density,          "[ ]");
		sprintf(scal->lbl_viscosity,        "[ ]");
		sprintf(scal->lbl_phase,            "[ ]");
	}
	else
	{
		//============
		// DIMENSIONAL
		//============
		// primary units
		scal->mass                = mass;
		scal->time                = time;
		scal->length              = length;
		scal->temperature         = temperature;
		scal->force               = force;
		// additional units
		volume                    = length*length*length;
		area                      = length*length;
		// secondary units
		scal->velocity            = length/time;
		scal->stress              = force/area;
		scal->strain_rate         = 1.0/time;
		scal->gravity_strength    = force/mass;
		scal->energy              = force*length;
		scal->power               = scal->energy/time;
		scal->heat_flux           = scal->power/area;
		scal->dissipation_rate    = scal->power/volume;
		// material parameters
		scal->density             = mass/volume;
		scal->viscosity           = scal->stress*time;
		scal->cpecific_heat       = scal->energy/mass/temperature;
		scal->conductivity        = scal->power/length/temperature;
		scal->heat_production     = scal->power/mass;
		scal->expansivity         = 1.0/temperature;
		scal->pressure_sensivity  = temperature/scal->stress;

		// output labels
		sprintf(scal->lbl_temperature,      "[K]");
		sprintf(scal->lbl_force,            "[N]");
		sprintf(scal->lbl_strain_rate,      "[1/s]");
		sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		sprintf(scal->lbl_density,          "[kg/m^3]");
		sprintf(scal->lbl_viscosity,        "[Pa*s]");

		if(scal->utype == _SI_)
		{
			// output labels
			sprintf(scal->lbl_time,      "[s]");
			sprintf(scal->lbl_length  ,  "[m]");
			sprintf(scal->lbl_velocity,  "[m/s]");
			sprintf(scal->lbl_stress ,   "[Pa]");
			sprintf(scal->lbl_heat_flux, "[W/m^2]");
		}
		else if(scal->utype == _GEO_)
		{
			yr     = 3600.0*24.0*365.0;
			km     = 1e3;
			cm     = 1e-2;
			cm_yr  = cm/yr;
			MPa    = 1e6;
			mW_msq = 1e-3;
			// modify scales
			scal->time       = time/yr;                // internal -> yr
			scal->length     = length/km;              // internal -> km
			scal->velocity   = scal->velocity/cm_yr;   // internal -> cm/yr
			scal->stress     = scal->stress/MPa;       // internal -> MPa
			scal->heat_flux  = scal->heat_flux/mW_msq; // internal -> mW/m^2
			// output labels
			sprintf(scal->lbl_time,      "[yr]");
			sprintf(scal->lbl_length  ,  "[km]");
			sprintf(scal->lbl_velocity,  "[cm/yr]");
			sprintf(scal->lbl_stress  ,  "[MPa]");
			sprintf(scal->lbl_heat_flux, "[mW/m^2]");
		}
	}

	// phase
	scal->phase = 1.0;
	sprintf(scal->lbl_phase, "[ ]");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar ComputePowerLawScaling(Scaling * scal, PetscScalar n)
{
	// power law constant scaling: 1 / stress^n / time

	return 1.0/pow(scal->stress, n)/scal->time;
}
//---------------------------------------------------------------------------


