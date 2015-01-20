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
		scal->angular_velocity    = angle/time;               sprintf(scal->lbl_angular_velocity, "[deg/s]"); // @

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
// compute characteristic values - migrated from NonDimensionalisation.c
void ComputeCharValues( UserCtx *user )
{
	nonDimUnits Characteristic_;

	memset( (void*)&Characteristic_, 0, sizeof(nonDimUnits) );

	if (user->DimensionalUnits==0){
		// no non-dimensionalization performed
		Characteristic_.Length      = 1.0;
		Characteristic_.Stress      = 1.0;
		Characteristic_.Temperature = 1.0;
		Characteristic_.Viscosity   = 1.0;

		PetscPrintf(PETSC_COMM_WORLD," Input units                    : Non-dimensional\n");

		// some useful values
		Characteristic_.km          = 1.0;
		Characteristic_.SecYear     = 1.0;
		Characteristic_.cmYear      = 1.0;
		Characteristic_.Myrs        = 1.0;
		Characteristic_.MPa         = 1.0;
		Characteristic_.Jmol        = 1.0;
	}

	else if (user->DimensionalUnits==1){
		// dimensionalization is performed - input parameters are in SI units
		// in most cases, these values are specified in the input file
		if (user->InputParamFile){
			Characteristic_.Length      = user->Characteristic.Length;     // m
			Characteristic_.Stress      = user->Characteristic.Stress;     // Pa
			Characteristic_.Temperature = user->Characteristic.Temperature;// Celcius !!
			Characteristic_.Viscosity   = user->Characteristic.Viscosity;  // Pa.s
		}
		else{
			// if we did not read an input file, set the values here
			Characteristic_.Length      = 10e3;  // m
			Characteristic_.Stress      = 100e6; // Pa
			Characteristic_.Temperature = 1000;  // Celcius !!
			Characteristic_.Viscosity   = 1e21;  // Pa.s
		}

		PetscPrintf(PETSC_COMM_WORLD," Input units                    : Dimensional \n");
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Length 	= %g [m] \n",Characteristic_.Length);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Stress 	= %g [Pa] \n",Characteristic_.Stress);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Temperature 	= %g  [Celcius/K] \n",Characteristic_.Temperature);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Viscosity 	= %g [Pa.s] \n",Characteristic_.Viscosity);
		PetscPrintf(PETSC_COMM_WORLD,"                                  Characteristic Time 	= %g [s] \n",Characteristic_.Viscosity/Characteristic_.Stress);

		// Some useful values
		Characteristic_.km          = 1000;
		Characteristic_.SecYear     = 365.25*24*3600;
		Characteristic_.cmYear      = Characteristic_.SecYear*100;  // from m/s to cm/year
		Characteristic_.Myrs        = 1e6;
		Characteristic_.MPa         = 1e6;
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD," Something wrong in defining ND units! \n");
	}

	// Derived values from the Characteristic_ values above
	Characteristic_.Time               = Characteristic_.Viscosity/Characteristic_.Stress;
	Characteristic_.Velocity           = Characteristic_.Length/Characteristic_.Time;
	Characteristic_.kg                 = Characteristic_.Stress*Characteristic_.Length*(Characteristic_.Time*Characteristic_.Time);
	Characteristic_.Density            = Characteristic_.kg/(Characteristic_.Length*Characteristic_.Length*Characteristic_.Length);
	Characteristic_.Strainrate         = 1.0/Characteristic_.Time;
	Characteristic_.ThermalExpansivity = 1.0/Characteristic_.Temperature;
	Characteristic_.Force              = Characteristic_.Stress*Characteristic_.Length*Characteristic_.Length;
	Characteristic_.Watt               = Characteristic_.Force*Characteristic_.Length/Characteristic_.Time;
	Characteristic_.Joule              = Characteristic_.Force*Characteristic_.Length;
	Characteristic_.HeatCapacity       = Characteristic_.Joule/Characteristic_.kg/Characteristic_.Temperature;
	Characteristic_.T_conductivity     = Characteristic_.Watt/Characteristic_.Length/Characteristic_.Temperature;
	Characteristic_.RadioactiveHeat    = Characteristic_.Watt/Characteristic_.Length/Characteristic_.Length/Characteristic_.Length;
	Characteristic_.Jmol               = user->GasConstant*Characteristic_.Temperature;

	user->Characteristic = Characteristic_;
}
//---------------------------------------------------------------------------
// transforms input parameters into ND parameters - migrated from NonDimensionalisation.c
void PerformNonDimension(UserCtx *user)
{
	PetscInt i;

	// domain
	user->W               = user->W/user->Characteristic.Length;
	user->L               = user->L/user->Characteristic.Length;
	user->H               = user->H/user->Characteristic.Length;
	user->x_left          = user->x_left /user->Characteristic.Length;
	user->y_front         = user->y_front/user->Characteristic.Length;
	user->z_bot           = user->z_bot  /user->Characteristic.Length;
	user->Setup_Diapir_Hi = user->Setup_Diapir_Hi/user->Characteristic.Length;

	// boundary conditions
	user->BC.Vy_front   = user->BC.Vy_front/user->Characteristic.Velocity;
	user->BC.Vy_back    = user->BC.Vy_back/user->Characteristic.Velocity;
	user->BC.Vz_top     = user->BC.Vz_top/user->Characteristic.Velocity;
	user->BC.Vz_bot     = user->BC.Vz_bot/user->Characteristic.Velocity;
	user->BC.Vx_left    = user->BC.Vx_left/user->Characteristic.Velocity;
	user->BC.Vx_right   = user->BC.Vx_right/user->Characteristic.Velocity;
	user->BC.Exx        = user->BC.Exx/user->Characteristic.Strainrate;
	user->BC.Eyy        = user->BC.Eyy/user->Characteristic.Strainrate;

	// time-stepping
	user->dt            = user->dt*user->Characteristic.SecYear/user->Characteristic.Time;
	user->dt_max        = user->dt_max*user->Characteristic.SecYear/user->Characteristic.Time;

	// temperature
	user->Temp_top      = user->Temp_top/user->Characteristic.Temperature;
	user->Temp_bottom   = user->Temp_bottom/user->Characteristic.Temperature;

	// gravity
	user->Gravity       = user->Gravity/( user->Characteristic.Length/(user->Characteristic.Time*user->Characteristic.Time));

	// optimization
	user->LowerViscosityCutoff = user->LowerViscosityCutoff/user->Characteristic.Viscosity;
	user->UpperViscosityCutoff = user->UpperViscosityCutoff/user->Characteristic.Viscosity;
	user->InitViscosity        = user->InitViscosity/user->Characteristic.Viscosity;

	// material properties
	for (i = 0; i < user->num_phases; i++)
	{
		user->PhaseProperties.mu[i]                 = user->PhaseProperties.mu[i] /user->Characteristic.Viscosity;
		user->PhaseProperties.FrankKamenetskii[i]   = user->PhaseProperties.FrankKamenetskii[i]/(1.0/user->Characteristic.Temperature);
		user->PhaseProperties.Powerlaw_e0[i]        = user->PhaseProperties.Powerlaw_e0[i]/(user->Characteristic.Strainrate);

		// I don't understand why "A" is interpreted as real pre-exponential factor of the stress-dependent dislocation creep.
		// Why then all rheologies that use "A" variable, also define the reference strain rate "Powerlaw_e0"?
		// What if I want a standard temperature-dependent power-low rheology?
		// I reinterpret this as a reference viscosity similar to "mu", since no one else is using it (a.p.)
		// user->PhaseProperties.A[i]               = user->PhaseProperties.A[i]/((pow(user->Characteristic.Stress,(-user->PhaseProperties.n_exponent[i])))/user->Characteristic.Time);
		user->PhaseProperties.A[i]                  = user->PhaseProperties.A[i]/user->Characteristic.Viscosity;

		// Don't scale activation energy, only scale gas constant with characteristic temperature (a.p)
		// user->PhaseProperties.E[i]                   = user->PhaseProperties.E[i]/user->Characteristic.Jmol;

		user->PhaseProperties.ElasticShearModule[i]     = user->PhaseProperties.ElasticShearModule[i]/user->Characteristic.Stress;
		user->PhaseProperties.ElasticBulkModule[i]      = user->PhaseProperties.ElasticBulkModule[i]/user->Characteristic.Stress;

		user->PhaseProperties.Cohesion[i]               = user->PhaseProperties.Cohesion[i]/user->Characteristic.Stress;
		user->PhaseProperties.CohesionAfterWeakening[i] = user->PhaseProperties.CohesionAfterWeakening[i]/user->Characteristic.Stress;

		user->PhaseProperties.T_Conductivity[i]         = user->PhaseProperties.T_Conductivity[i] /user->Characteristic.T_conductivity;
		user->PhaseProperties.HeatCapacity[i]           = user->PhaseProperties.HeatCapacity[i]	/user->Characteristic.HeatCapacity;
		user->PhaseProperties.RadioactiveHeat[i]        = user->PhaseProperties.RadioactiveHeat[i]/user->Characteristic.RadioactiveHeat;

		user->PhaseProperties.rho[i]                    = user->PhaseProperties.rho[i]/user->Characteristic.Density;
		user->PhaseProperties.Density_T0[i]             = user->PhaseProperties.Density_T0[i]/user->Characteristic.Temperature;
		user->PhaseProperties.ThermalExpansivity[i]     = user->PhaseProperties.ThermalExpansivity[i]/(1.0/user->Characteristic.Temperature);
	}

	// pushing block parameters
	user->Pushing.L_block        = user->Pushing.L_block/user->Characteristic.Length;
	user->Pushing.W_block        = user->Pushing.W_block/user->Characteristic.Length;
	user->Pushing.H_block        = user->Pushing.H_block/user->Characteristic.Length;
	user->Pushing.x_center_block = user->Pushing.x_center_block/user->Characteristic.Length;
	user->Pushing.y_center_block = user->Pushing.y_center_block/user->Characteristic.Length;
	user->Pushing.z_center_block = user->Pushing.z_center_block/user->Characteristic.Length;
	for (i = 0; i < user->Pushing.num_changes; i++){
		user->Pushing.V_push[i] = user->Pushing.V_push[i]*1e-2/user->Characteristic.SecYear/user->Characteristic.Velocity;
		user->Pushing.omega[i]  = user->Pushing.omega[i]*(M_PI/180.0)/(user->Characteristic.SecYear/user->Characteristic.Time); //[radians]
	}
	for (i=0;i<user->Pushing.num_changes+1; i++){
		user->Pushing.time[i]   = user->Pushing.time[i]*1e6*user->Characteristic.SecYear/user->Characteristic.Time; //from Myr
	}
}
/*==========================================================================================================*/

