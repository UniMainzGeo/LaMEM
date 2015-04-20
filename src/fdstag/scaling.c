//---------------------------------------------------------------------------
//..........................   SCALING ROUTINES   ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Parsing.h"
#include "solVar.h"
#include "scaling.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ScalingReadFromFile"
PetscErrorCode ScalingReadFromFile(Scaling *scal, FILE *fp)
{
	PetscInt    flg;
	PetscInt    found;
	char        utype [MAX_NAME_LEN];
	PetscScalar length, viscosity, temperature, stress;
	PetscScalar acceleration, time, mass, force;

	PetscFunctionBegin;

	// read model setup
	parse_GetString(fp, "units", utype, MAX_NAME_LEN, &flg);

	if(flg)
	{
		if     (!strcmp(utype, "none"))   scal->utype = _NONE_;
		else if(!strcmp(utype, "si"))     scal->utype = _SI_;
		else if(!strcmp(utype, "geo"))    scal->utype = _GEO_;
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect type of units: %s", utype);
	}

	if(scal->utype != _NONE_)
	{
		// set
		length      = 0.0;
		viscosity   = 0.0;
		temperature = 0.0;
		stress      = 0.0;

		// read
		parse_GetDouble(fp, "Characteristic.Length",     &length ,     &found);
		parse_GetDouble(fp, "Characteristic.Viscosity",  &viscosity,   &found);
		parse_GetDouble(fp, "Characteristic.Temperature",&temperature, &found);
		parse_GetDouble(fp, "Characteristic.Stress",     &stress,      &found);

		// check
		if(length      == 0.0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Define Characteristic.Length parameter\n");
		if(viscosity   == 0.0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Define Characteristic.Viscosity parameter\n");
		if(temperature == 0.0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Define Characteristic.Temperature parameter\n");
		if(stress      == 0.0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Define Characteristic.Stress parameter\n");

		// print
		PetscPrintf(PETSC_COMM_WORLD," Characteristic Length      = %g [m] \n",    length);
		PetscPrintf(PETSC_COMM_WORLD," Characteristic Viscosity   = %g [Pa*s] \n", viscosity);
		PetscPrintf(PETSC_COMM_WORLD," Characteristic Temperature = %g [C/K] \n",  temperature);
		PetscPrintf(PETSC_COMM_WORLD," Characteristic Stress      = %g [Pa] \n",   stress);

		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

		// compute additional parameters
		// scaling is currently consistent with the Second Newton's Law
		// for quasi-static problems we can have additional free parameter!

		time         = viscosity/stress;
		force        = stress*length*length;
		acceleration = length/time/time;
		mass         = force/acceleration;

		// store extended input parameters
		scal->inp_mass        = mass;
		scal->inp_time        = time;
		scal->inp_length      = length;
		scal->inp_temperature = temperature;
		scal->inp_force       = force;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ScalingCreate"
PetscErrorCode ScalingCreate(Scaling *scal)
{
	// characteristic values must ALWAYS be given in SI units
	PetscScalar mass, time, length, temperature, force;
	PetscScalar angle, area, volume, stress, energy, power;
	PetscScalar yr, Myr, km, cm, cm_yr, MPa, mW;

	// access extended input parameters
	if(scal->utype != _NONE_)
	{
		mass        = scal->inp_mass;
		time        = scal->inp_time;
		length      = scal->inp_length;
		temperature = scal->inp_temperature;
		force       = scal->inp_force;
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
		scal->time_si             = 1.0;
		scal->length              = 1.0;   sprintf(scal->lbl_length,           "[ ]");
		scal->temperature         = 1.0;   sprintf(scal->lbl_temperature,      "[ ]");
		scal->force               = 1.0;   sprintf(scal->lbl_force,            "[ ]");
		scal->angle               = 1.0;   sprintf(scal->lbl_angle,            "[ ]");

		// secondary units
		scal->velocity            = 1.0;   sprintf(scal->lbl_velocity,         "[ ]");
		scal->stress              = 1.0;   sprintf(scal->lbl_stress,           "[ ]");
		scal->stress_si           = 1.0;
		scal->strain_rate         = 1.0;   sprintf(scal->lbl_strain_rate,      "[ ]");
		scal->acceleration        = 1.0;
		scal->gravity_strength    = 1.0;
		scal->energy              = 1.0;
		scal->power               = 1.0;
		scal->heat_flux           = 1.0;   sprintf(scal->lbl_heat_flux,        "[ ]");
		scal->dissipation_rate    = 1.0;   sprintf(scal->lbl_dissipation_rate, "[ ]");
		scal->angular_velocity    = 1.0;   sprintf(scal->lbl_angular_velocity, "[ ]");
		scal->volumetric_force    = 1.0;   sprintf(scal->lbl_volumetric_force, "[ ]");

		// material parameters
		scal->density             = 1.0;   sprintf(scal->lbl_density,          "[ ]");
		scal->viscosity           = 1.0;   sprintf(scal->lbl_viscosity,        "[ ]");
		scal->cpecific_heat       = 1.0;
		scal->conductivity        = 1.0;
		scal->heat_production     = 1.0;
		scal->expansivity         = 1.0;

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
		scal->time_si             = time;
		scal->length              = length;                   sprintf(scal->lbl_length  ,         "[m]");
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[K]");
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @

		// secondary units
		scal->velocity            = length/time;              sprintf(scal->lbl_velocity,         "[m/s]");
		scal->stress              = stress;                   sprintf(scal->lbl_stress ,          "[Pa]");
		scal->stress_si           = stress;
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->acceleration        = length/time/time; // m/s2
		scal->gravity_strength    = force/mass;
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area;               sprintf(scal->lbl_heat_flux,        "[W/m^2]");
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/time;               sprintf(scal->lbl_angular_velocity, "[deg/s]"); // @
		scal->volumetric_force    = force/volume;             sprintf(scal->lbl_volumetric_force, "[N/m^3]");

		// material parameters
		scal->density             = mass/volume;              sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = stress*time;              sprintf(scal->lbl_viscosity,        "[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;
		scal->conductivity        = power/length/temperature;
		scal->heat_production     = power/mass;
		scal->expansivity         = 1.0/temperature;

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
		scal->mass                = mass;
		scal->time                = time/Myr;                 sprintf(scal->lbl_time,             "[Myr]");   // @
		scal->time_si             = time;
		scal->length              = length/km;                sprintf(scal->lbl_length  ,         "[km]");    // @
		scal->temperature         = temperature;              sprintf(scal->lbl_temperature,      "[C]");     // @
		scal->force               = force;                    sprintf(scal->lbl_force,            "[N]");
		scal->angle               = angle;                    sprintf(scal->lbl_angle,            "[deg]");   // @

		// secondary units
		scal->velocity            = length/time/cm_yr;        sprintf(scal->lbl_velocity,         "[cm/yr]"); // @
		scal->stress              = stress/MPa;               sprintf(scal->lbl_stress ,          "[MPa]");   // @
		scal->stress_si           = stress;
		scal->strain_rate         = 1.0/time;                 sprintf(scal->lbl_strain_rate,      "[1/s]");
		scal->acceleration        = length/time/time;
		scal->gravity_strength    = force/mass;
		scal->energy              = energy;
		scal->power               = power;
		scal->heat_flux           = power/area/mW;            sprintf(scal->lbl_heat_flux,        "[mW/m^2]");  // @
		scal->dissipation_rate    = power/volume;             sprintf(scal->lbl_dissipation_rate, "[W/m^3]");
		scal->angular_velocity    = angle/(time/Myr);         sprintf(scal->lbl_angular_velocity, "[deg/Myr]"); // @
		scal->volumetric_force    = force/volume;             sprintf(scal->lbl_volumetric_force, "[N/m^3]");

		// material parameters
		scal->density             = mass/volume;              sprintf(scal->lbl_density,          "[kg/m^3]");
		scal->viscosity           = stress*time;              sprintf(scal->lbl_viscosity,        "log[Pa*s]");
		scal->cpecific_heat       = energy/mass/temperature;
		scal->conductivity        = power/length/temperature;
		scal->heat_production     = power/mass;
		scal->expansivity         = 1.0/temperature;

	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// scaling of input parameters (UserCtx)
void ScalingInput(Scaling *scal, UserCtx *user)
{
	PetscInt i;

	// domain
	user->W               /= scal->length;
	user->L               /= scal->length;
	user->H               /= scal->length;
	user->x_left          /= scal->length;
	user->y_front         /= scal->length;
	user->z_bot           /= scal->length;
	user->Setup_Diapir_Hi /= scal->length;

	// boundary conditions
	user->BC.Vy_front     /= scal->velocity;
	user->BC.Vy_back      /= scal->velocity;
	user->BC.Vz_top       /= scal->velocity;
	user->BC.Vz_bot       /= scal->velocity;
	user->BC.Vx_left      /= scal->velocity;
	user->BC.Vx_right     /= scal->velocity;
	user->BC.Exx          /= scal->strain_rate;
	user->BC.Eyy          /= scal->strain_rate;

	// time-stepping
	user->dt              /= scal->time;
	user->dt_max          /= scal->time;

	// temperature
	user->Temp_top        /= scal->temperature;
	user->Temp_bottom     /= scal->temperature;

	// gravity
	user->Gravity         /= scal->acceleration;

	user->LowerViscosityCutoff /= scal->viscosity;
	user->UpperViscosityCutoff /= scal->viscosity;
	user->InitViscosity        /= scal->viscosity;
//	user->PlastViscosity       /= scal->viscosity;

	user->DII_ref              /= scal->strain_rate;

	// pushing block parameters
	user->Pushing.L_block        /= scal->length;
	user->Pushing.W_block        /= scal->length;
	user->Pushing.H_block        /= scal->length;
	user->Pushing.x_center_block /= scal->length;
	user->Pushing.y_center_block /= scal->length;
	user->Pushing.z_center_block /= scal->length;

	for(i = 0; i < user->Pushing.num_changes; i++)
	{
		user->Pushing.V_push[i] /= scal->velocity;
		user->Pushing.omega[i]  /= scal->angular_velocity;
	}

	for (i=0; i < user->Pushing.num_changes+1; i++)
	{
		user->Pushing.time[i]  /= scal->time;
	}

	// scale mesh segment delimiters
	ScalingMeshSegDir(scal, &user->mseg_x);
	ScalingMeshSegDir(scal, &user->mseg_y);
	ScalingMeshSegDir(scal, &user->mseg_z);

}
//---------------------------------------------------------------------------
// scaling material parameters
// NOTE: [1] activation energy is not scaled
//       [2] activation volume is multiplied with characteristic stress in SI units
void ScalingMatProp(Scaling *scal, Material_t *phases, PetscInt numPhases)
{
	PetscInt     i;

	// scale dimensional parameters
	for(i = 0; i < numPhases; i++)
	{
		phases[i].rho     /= scal->density;

		// diffusion creep
		phases[i].Bd      *= scal->viscosity;
		phases[i].Vd      *= scal->stress_si;

		// dislocation creep (power-law)
		phases[i].Bn      *= pow(scal->stress_si, phases[i].n)*scal->time_si;
		phases[i].Vn      *= scal->stress_si;

		// Peierls creep
		phases[i].Bp      /=  scal->strain_rate;
		phases[i].Vp      *=  scal->stress_si;
		phases[i].taup    /=  scal->stress_si;

		// elasticity
		phases[i].G       /= scal->stress_si;
		phases[i].K       /= scal->stress_si;

		// plasticity
		phases[i].ch      /= scal->stress_si;
		phases[i].fr      /= scal->angle;

		// temperature
		phases[i].alpha   /= scal->expansivity;
		phases[i].Cp      /= scal->cpecific_heat;
		phases[i].k       /= scal->conductivity;
		phases[i].A       /= scal->heat_production;
	}
}
//---------------------------------------------------------------------------
// scaling material parameter limits
// NOTE: [1] gas constant is multiplied with characteristic temperature
void ScalingMatParLim(Scaling *scal, MatParLim *matLim)
{
	// scale gas constant with characteristic temperature
	matLim->Rugc *= scal->temperature;
}
//---------------------------------------------------------------------------
void ScalingMeshSegDir(Scaling *scal, MeshSegInp *msi)
{
	PetscInt i;

	// scale segment delimiters
	for(i = 0; i < msi->nsegs-1; i++)
	{
		msi->delims[i] /= scal->length;
	}
}
//---------------------------------------------------------------------------
