//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Parsing.h"
#include "Utils.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "input.h"
//---------------------------------------------------------------------------
// set default code parameters and read input file, if required
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitCode"
PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user)
{
	PetscMPIInt    size;
	char          *all_options;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED
	PetscOptionsGetAll( &all_options ); // copy all command line args

	// set default values for parameters
	ierr = FDSTAGSetDefaultValues(user); CHKERRQ(ierr);

	// initialize material properties
	ierr = FDSTAGInitMaterialProp(user); CHKERRQ(ierr);

	// read an input file if required - using the parser
	ierr = FDSTAGReadInputFile(jr, user); CHKERRQ(ierr);

	// read command line options
	ierr = FDSTAGReadCommLine(user); CHKERRQ(ierr);

	// set this option to monitor actual option usage
	PetscOptionsInsertString("-options_left");

	// error check and info print:
	// we need at least 2 cells in each direction (not 1) and dt>0! - IS THIS NECESSARY?
	if (user->nel_x==1){ user->nel_x=2; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the x-direction! I have increased this. \n");}
	if (user->nel_y==1){ user->nel_y=2; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the y-direction! I have increased this. \n");}
	if (user->dt==0.0) { user->dt=1e-3; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG dt cannot be zero. Value has been changed to dt = 1e-3. \n");}

	// print information about simulation
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Total # of cpu's               : %lld \n",(LLD)size);

	// info about particles if used
	if (user->ParticleInput==1)
	{
		// check whether a sufficient amount of particles are specified; if not increase it.
		if(user->NumPartX < 2) { user->NumPartX = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in x-direction; Increasing\n"); }
		if(user->NumPartY < 2) { user->NumPartY = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in y-direction; Increasing\n"); }
		if(user->NumPartZ < 2) { user->NumPartZ = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in z-direction; Increasing\n"); }
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of tracers/cell         : [%lld,%lld,%lld] \n",(LLD)(user->NumPartX),(LLD)(user->NumPartY),(LLD)(user->NumPartZ));

	// show initial timestep (WARNING! which years if nondimensional?)
	if (user->DimensionalUnits==1){
		PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g years \n",user->dt);
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g  \n",user->dt);
	}
	// boundary conditions
	PetscPrintf(PETSC_COMM_WORLD," BC employed                    : BC.[Left=%lld Right=%lld; Front=%lld Back=%lld; Lower=%lld Upper=%lld] \n",
			(LLD)(user->BC.LeftBound), (LLD)(user->BC.RightBound), (LLD)(user->BC.FrontBound), (LLD)(user->BC.BackBound), (LLD)(user->BC.LowerBound), (LLD)(user->BC.UpperBound) );

	// compute characteristic values - should be done in scaling.c
	ComputeCharValues(user);

	// compute material prop for default setup
	if (!(user->InputParamFile))
	{
		// define phase material properties for default setup (Falling Block Test)*/
		user->PhaseProperties.mu[0] = 1.0*user->Characteristic.Viscosity; user->PhaseProperties.rho[0] = 1.0*user->Characteristic.Density;
		user->PhaseProperties.mu[1] = 1.0*user->Characteristic.Viscosity; user->PhaseProperties.rho[1] = 2.0*user->Characteristic.Density;
	}

	// perform non-dimensionalisation of all input parameters - should be done in scaling.c
	PerformNonDimension(user);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// set default values for parameters
#undef __FUNCT__
#define __FUNCT__ "FDSTAGSetDefaultValues"
PetscErrorCode FDSTAGSetDefaultValues(UserCtx *user)
{
	// domain
	user->W                = 1.0;
	user->L                = 1.0;
	user->H                = 1.0;
	user->x_left           = -0.5*(user->W)*0.0;
	user->y_front          = -0.5*(user->L)*0.0;
	user->z_bot            = -0.5*(user->H)*0.0;

	// set default grid sizes
	user->nel_x            = 8;
	user->nel_y            = 8;
	user->nel_z            = 8;

	user->NumPartX         = 2; // Number of particles per cell in x-direction
	user->NumPartY         = 2; // Number of particles per cell in y-direction
	user->NumPartZ         = 2; // Number of particles per cell in z-direction

	// scaling
	user->DimensionalUnits = 0;

	// time-stepping
	user->time_end         = 1.0;
	user->save_timesteps   = 1;
	user->CFL              = 0.5;
	user->dt_max           = 1e6; // maximum timestep
	user->dt               = 1e-3;// initial timestep

	// temperature
	user->Temp_top         = 0.0;
	user->Temp_bottom      = 1.0;
	user->GasConstant      = 1.0;            // REMOVE

	user->Gravity          = 1.0;
	user->GravityAngle     = 0.0; // angle of gravity with z-axis (can be changed in the x-z plane)

	// input/output
	sprintf(user->OutputFile,                    "FallingBlock3D");
	sprintf(user->ParticleFilename,              "ParticlesInput.dat");
	sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");

	// FDSTAG Canonical Default Model Setup
	user->msetup            = BLOCK;

	// write partitioning
	user->SavePartitioning  = PETSC_FALSE;

	// linear solver settings
	user->SkipStokesSolver     = PETSC_FALSE;
	user->use_fdstag_canonical = PETSC_TRUE;  // request native staggered grid discretization

	// boundary conditions
	user->BC.Vy_front      = 0;
	user->BC.Vy_back       = 0;
	user->BC.Vz_bot        = 0;
	user->BC.Vz_top        = 0;
	user->BC.Vx_left       = 0;
	user->BC.Vx_right      = 0;
	user->BC.Exx           = 0;
	user->BC.Eyy           = 0;
	user->BC.InternalBound = 0; // 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.UpperBound    = 1; // 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.LowerBound    = 1; // 0-free surface, 1-free slip					, 2-no-slip
	user->BC.LeftBound     = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no slip
	user->BC.RightBound    = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.FrontBound    = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.BackBound     = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip

	user->ParticleInput      = 1;  // 0-do not use particles to track phases; 1-do use particles to track phases ???
	user->SaveParticles      = 0;  // Save particles or not?
	user->restart            = 1;  // Are we restarting a simulation?
	user->save_breakpoints   = 10; // After how many steps do we make a breakpoint?
	user->break_point_number = 0;  // The number of the breakpoint file

	// optimization
	user->LowerViscosityCutoff  = 1.0;  // lowermost viscosity cutoff in model
	user->UpperViscosityCutoff  = 1e30; // uppermost viscosity cutoff in model
	user->InitViscosity         = 0.0;  // initial viscosity

	user->mpi_group_id = 0;

	// phases
	user->num_phases   = 2;  // the default # of phases. In case we set props from the command line, and we use a folding setup combined with FDSTAG we might have to do something smarter (check the setup and increase this)

	// Default Pushing BC parameters
	user->AddPushing                    = 0;
	user->Pushing.reset_pushing_coord   = 0;
	user->Pushing.theta                 = 0.0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGReadInputFile"
PetscErrorCode FDSTAGReadInputFile(JacRes *jr, UserCtx *user)
{
	 // parse the input file

	FILE *fp;
	PetscInt found;
	double d_values[1000], data;
	PetscInt i_values[1000];
	PetscInt nv, i;
	char setup_name[PETSC_MAX_PATH_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// return if no input file
	if(!user->InputParamFile) PetscFunctionReturn(0);

	fp = fopen( user->ParamFile, "r" );
	if(!fp)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s", user->ParamFile);
	}

	// read number of cells
	parse_GetInt( fp, "nel_x",   &user->nel_x, &found );
	parse_GetInt( fp, "nel_y",   &user->nel_y, &found );
	parse_GetInt( fp, "nel_z",   &user->nel_z, &found );

	// only number of cells is a relevant parameter for FDSTAG
	user->nnode_x = user->nel_x+1;
	user->nnode_y = user->nel_y+1;
	user->nnode_z = user->nel_z+1;

	// Non-dimensionalization parameters
	ierr = ScalingReadFromFile(&jr->scal, fp); CHKERRQ(ierr);

	parse_GetInt   ( fp, "DimensionalUnits", &user->DimensionalUnits, &found );
	parse_GetDouble( fp, "Characteristic.Length",     &data, &found ); if (user->DimensionalUnits==1){ user->Characteristic.Length      = data;}// read data
	parse_GetDouble( fp, "Characteristic.Viscosity",  &data, &found ); if (user->DimensionalUnits==1){ user->Characteristic.Viscosity   = data;}// read data
	parse_GetDouble( fp, "Characteristic.Temperature",&data, &found ); if (user->DimensionalUnits==1){ user->Characteristic.Temperature = data;}// read data
	parse_GetDouble( fp, "Characteristic.Stress",     &data, &found ); if (user->DimensionalUnits==1){ user->Characteristic.Stress      = data;}// read data

	// domain parameters
	parse_GetDouble( fp, "L", &user->L, &found );
	parse_GetDouble( fp, "W", &user->W, &found );
	parse_GetDouble( fp, "H", &user->H, &found );
	parse_GetDouble( fp, "x_left",  &user->x_left,  &found );
	parse_GetDouble( fp, "y_front", &user->y_front, &found );
	parse_GetDouble( fp, "z_bot",   &user->z_bot,   &found );

	//==============
	// MESH SEGMENTS
	//==============
	ierr = ReadMeshSegDir(fp, "seg_x", user->x_left,  user->x_left  + user->W, &user->nel_x, &user->mseg_x, user->DimensionalUnits, user->Characteristic.Length); CHKERRQ(ierr);
	ierr = ReadMeshSegDir(fp, "seg_y", user->y_front, user->y_front + user->L, &user->nel_y, &user->mseg_y, user->DimensionalUnits, user->Characteristic.Length); CHKERRQ(ierr);
	ierr = ReadMeshSegDir(fp, "seg_z", user->z_bot,   user->z_bot   + user->H, &user->nel_z, &user->mseg_z, user->DimensionalUnits, user->Characteristic.Length); CHKERRQ(ierr);

	// read model setup
	parse_GetString(fp, "msetup", setup_name, PETSC_MAX_PATH_LEN-1, &found);
	if(found)
	{
		if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#ERROR! Incorrect model setup: %s", setup_name);
	}

	parse_GetString( fp, "OutputFile", user->OutputFile, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetInt( fp,    "save_breakpoints", &user->save_breakpoints, &found );
	parse_GetInt( fp,    "save_timesteps", &user->save_timesteps, &found );
	parse_GetInt( fp,    "time_end", &user->time_end, &found );
	parse_GetDouble( fp, "CFL", &user->CFL, &found );
	parse_GetDouble( fp, "dt_max", &user->dt_max, &found );

	// Particle related variables
	parse_GetInt( fp,    "ParticleInput", &user->ParticleInput, &found );
	parse_GetString( fp, "ParticleFilename", user->ParticleFilename, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetString( fp, "LoadInitialParticlesDirectory", user->LoadInitialParticlesDirectory, PETSC_MAX_PATH_LEN-1, &found );
	if (!found){
		sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetString( fp, "SaveInitialParticlesDirectory", user->SaveInitialParticlesDirectory, PETSC_MAX_PATH_LEN-1, &found );
	if (!found){
		sprintf(user->SaveInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetInt( fp,    "SaveParticles", &user->SaveParticles, &found );

	parse_GetInt( fp,    "NumPartX", &user->NumPartX, &found );
	parse_GetInt( fp,    "NumPartY", &user->NumPartY, &found );
	parse_GetInt( fp,    "NumPartZ", &user->NumPartZ, &found );

	// boundary conditions
	parse_GetDouble( fp, "BC.Eyy", &user->BC.Eyy, &found );
	parse_GetDouble( fp, "BC.Exx", &user->BC.Exx, &found );

	parse_GetInt( fp, "BC.InternalBound", &user->BC.InternalBound, &found );
	parse_GetInt( fp, "BC.LeftBound", &user->BC.LeftBound, &found );
	parse_GetInt( fp, "BC.RightBound", &user->BC.RightBound, &found );
	parse_GetInt( fp, "BC.FrontBound", &user->BC.FrontBound, &found );
	parse_GetInt( fp, "BC.BackBound", &user->BC.BackBound, &found );
	parse_GetInt( fp, "BC.LowerBound", &user->BC.LowerBound, &found );
	parse_GetInt( fp, "BC.UpperBound", &user->BC.UpperBound, &found );
	parse_GetDouble( fp, "Temp_top", &user->Temp_top, &found );
	parse_GetDouble( fp, "Temp_bottom", &user->Temp_bottom, &found );

	// optimization
	parse_GetDouble( fp, "LowerViscosityCutoff", &user->LowerViscosityCutoff, &found );
	parse_GetDouble( fp, "UpperViscosityCutoff", &user->UpperViscosityCutoff, &found );
	parse_GetDouble( fp, "InitViscosity",        &user->InitViscosity,        &found );

	parse_GetDouble( fp, "Gravity", &user->Gravity, &found );
	parse_GetDouble( fp, "GasConstant", &user->GasConstant, &found );

	// Pushing Parameters
	parse_GetInt( fp,    "AddPushing", &user->AddPushing, &found );
	parse_GetInt( fp,    "Pushing.num_changes", &user->Pushing.num_changes, &found );
	parse_GetInt( fp,    "Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, &found );
	parse_GetDouble( fp, "Pushing.L_block", &user->Pushing.L_block, &found );
	parse_GetDouble( fp, "Pushing.W_block", &user->Pushing.W_block, &found );
	parse_GetDouble( fp, "Pushing.H_block", &user->Pushing.H_block, &found );
	parse_GetDouble( fp, "Pushing.x_center_block", &user->Pushing.x_center_block, &found );
	parse_GetDouble( fp, "Pushing.y_center_block", &user->Pushing.y_center_block, &found );
	parse_GetDouble( fp, "Pushing.z_center_block", &user->Pushing.z_center_block, &found );
	parse_GetDoubleArray( fp, "Pushing.time",   &nv, d_values, &found );  for( i=0; i<user->Pushing.num_changes+1; i++ ) { user->Pushing.time[i]   = d_values[i];}
	parse_GetDoubleArray( fp, "Pushing.V_push", &nv, d_values, &found );  for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.V_push[i] = d_values[i];}
	parse_GetDoubleArray( fp, "Pushing.omega",  &nv, d_values, &found  ); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.omega[i]  = d_values[i];}
	parse_GetIntArray( fp, "Pushing.coord_advect",&nv, i_values, &found); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.coord_advect[i] = i_values[i];}
	parse_GetIntArray( fp, "Pushing.dir",         &nv, i_values, &found); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.dir[i]          = i_values[i];}

	fclose(fp);

	ierr = ReadMaterialProperties(user); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ReadMaterialProperties"
PetscErrorCode ReadMaterialProperties(UserCtx *user)
{
	PetscInt      found_data;
	PetscInt      n_phases,n_attrs, n_types;
	Phase         *phases, p;
	Attribute     *attrs, a;
	AttributeType *types, t;
	PetscInt      PP,AA,TT;
	PetscInt      attr_match, type_match;
	PetscInt      p_idx, a_idx, t_idx;

	PetscInt      *check_phase;
	PetscBool      empty_phase;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscMemzero(&user->PhaseProperties, sizeof(PhaseProps)); CHKERRQ(ierr);

	//-----------------------------------------------------------------------
	// Read material parameters from the input file (new, more general, input format for material properties)
	//
	// At this stage, the new data is read from the input file and added to the phase properties array.
	//-----------------------------------------------------------------------

	MaterialCreate( &user->PhaseMaterialProperties );
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ConstantViscosityReadFromFile,                     &found_data); // Constant viscosity params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &PowerLawViscosityReadFromFile,                     &found_data); // Power-law viscosity params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &TempDepViscosityReadFromFile,                      &found_data); // Temp-dep viscosity params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &TempDepNoPowerLawViscosityReadFromFile,            &found_data); // Temp-dep viscosity params without power law
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ConstantElasticityReadFromFile,                    &found_data); // Elasticity params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &TemperatureDependentDensityReadFromFile,           &found_data); // Density params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ConvectionDensityReadFromFile,                     &found_data); // Density params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ArtificialTemperatureDependentDensityReadFromFile, &found_data); // Density params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ConstantEnergyReadFromFile, 				         &found_data); // Energy equation params
	MaterialReadFromFile( user->PhaseMaterialProperties, user->ParamFile, &ConstantPlasticityReadFromFile,                    &found_data); // Plasticity params

	// Print an overview of Material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Phase material parameters read from %s: \n",user->ParamFile);

	MaterialGetAllPhases( user->PhaseMaterialProperties, &n_phases, &phases );

	// check total number of phases found
	if(n_phases > max_num_phases)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many phases in the input file! Actual: %lld, Maximum: %lld", (LLD)n_phases, (LLD)max_num_phases);
	}

	// initialize phase checking array
	ierr = makeIntArray(&check_phase, NULL, n_phases); CHKERRQ(ierr);

	for(PP = 0; PP < n_phases; PP++) check_phase[PP] = 0;

	// store total number of phases
	user->num_phases = n_phases;

	for(PP = 0; PP < n_phases; PP++)
	{
		p = phases[PP];

		// check phase numbering
		if(p->phase_number > n_phases-1)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Phase numbering out of bound! Phase ID: %lld, Maximum ID: %lld", (LLD)p->phase_number, (LLD)(n_phases-1));
		}
		check_phase[p->phase_number] = 1;

		p_idx = p->P_id;

		PhaseGetAllAttributes( p, &n_attrs, &attrs );

		for( AA=0; AA<n_attrs; AA++ ) {
			a = attrs[AA];
			a_idx = a->A_id;
			AttributeGetAllTypes( a, &n_types, &types );

			for( TT=0; TT<n_types; TT++ ) {
				t = types[TT];
				t_idx = t->T_id;

				AttributeCompare( a, "VISCOSITY", &attr_match );
				if( attr_match == _TRUE ) {
					// constant
					AttributeTypeCompare( t, "constant", &type_match );
					if( type_match == _TRUE ) {
						double eta0;

						MaterialGetConstantViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &eta0 );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), VISCOSITY, constant [%lld,%lld,%lld] = %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, eta0 );

						// Add to PhaseProperties
						user->PhaseProperties.ViscosityLaw[p->phase_number] = 1; 		// constant viscosity
						user->PhaseProperties.mu[p->phase_number] 			= eta0;

					}
					AttributeTypeCompare( t, "power_law", &type_match );
					if( type_match == _TRUE ) {
						double eta0;
						double N_exp;
						double e0;

						MaterialGetPowerLawViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &eta0, &N_exp, &e0 );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), VISCOSITY, power-law [%lld,%lld,%lld] = %e, %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, eta0, N_exp, e0 );


						// Add to PhaseProperties
						user->PhaseProperties.ViscosityLaw[p->phase_number] = 2; 		// powerlaw viscosity, given by eta=eta0*(e2nd/e0)^(1/n-1)
						user->PhaseProperties.mu[p->phase_number] 			= eta0;
						user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
						user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)

					}
					AttributeTypeCompare( t, "temp_dep", &type_match );
					if( type_match == _TRUE ) {
						double PreExpFactor;
						double N_exp;
						double e0;
						double ActivationEnergy;

						MaterialGetTempDepViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &PreExpFactor, &N_exp, &e0, &ActivationEnergy );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), VISCOSITY, temp-dep [%lld,%lld,%lld] = %e, %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, PreExpFactor, N_exp, e0, ActivationEnergy );

						// Add to PhaseProperties
						user->PhaseProperties.ViscosityLaw[p->phase_number] = 4;
						user->PhaseProperties.A[p->phase_number] 			= PreExpFactor;
						user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
						user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)
						user->PhaseProperties.E[p->phase_number]			= ActivationEnergy;
					}
					AttributeTypeCompare( t, "tempdep_nopowerlaw", &type_match );
					if( type_match == _TRUE ) {
						double PreExpFactor;
						double N_exp;
						double e0;
						double ActivationEnergy;

						MaterialGetTempDepNoPowerLawViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &PreExpFactor, &N_exp, &e0, &ActivationEnergy );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), VISCOSITY, tempdep_nopowerlaw [%lld,%lld,%lld] = %e, %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, PreExpFactor, N_exp, e0, ActivationEnergy );


						// Add to PhaseProperties
						user->PhaseProperties.ViscosityLaw[p->phase_number] = 5;
						user->PhaseProperties.A[p->phase_number] 			= PreExpFactor;
						user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
						user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)
						user->PhaseProperties.E[p->phase_number]			= ActivationEnergy;
					}
					// others
				}
				AttributeCompare( a, "ELASTICITY", &attr_match );
				if( attr_match == _TRUE ) {
					// constant
					AttributeTypeCompare( t, "constant", &type_match );
					if( type_match == _TRUE ) {
						double shear, bulk;

						MaterialGetConstantElasticityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &shear, &bulk );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), ELASTICITY, constant [%lld,%lld,%lld] = %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, shear,bulk );

						user->PhaseProperties.ElasticShearModule[p->phase_number] 		= shear; 	// add to PhaseProperties
						user->PhaseProperties.ElasticBulkModule[p->phase_number] 		= bulk; 	// add to PhaseProperties

					}
					// others
				}
				AttributeCompare( a, "DENSITY", &attr_match );
				if( attr_match == _TRUE ) {
					// temp dep.
					AttributeTypeCompare( t, "temperature_dependent", &type_match );
					if( type_match == _TRUE ) {
						double rho0, alpha, T0;

						MaterialGetTemperatureDependentDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0, &alpha, &T0 );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), DENSITY, temperature-dependent [%lld,%lld,%lld] = %e, %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0,alpha, T0 );


						user->PhaseProperties.DensityLaw[p->phase_number] 	= 1; 		// T-dependent density
						user->PhaseProperties.rho[p->phase_number] 			= rho0; 	// add to PhaseProperties
						user->PhaseProperties.ThermalExpansivity[p->phase_number] 		= alpha; 	// add to PhaseProperties
						user->PhaseProperties.Density_T0[p->phase_number] 	= T0; 		// add to PhaseProperties

					}
					// convection setup
					AttributeTypeCompare( t, "convection", &type_match );
					if( type_match == _TRUE ) {
						double rho0, Ra, T0;

						MaterialGetConvectionDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0, &Ra, &T0 );
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), DENSITY, convection [%lld,%lld,%lld] = %e, %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0,Ra, T0 );


						user->PhaseProperties.DensityLaw[p->phase_number] 	= 2; 				// Density used in Ra-convection simulations where the rhs of the force balance 2 is Ra*T*g instead of rho*g
						user->PhaseProperties.rho[p->phase_number] 			= rho0; 			// add to PhaseProperties
						user->PhaseProperties.Ra[p->phase_number] 			= Ra; 				// add to PhaseProperties
						user->PhaseProperties.Density_T0[p->phase_number] 	= T0; 				// add to PhaseProperties
					}
					// artificial temp dep.
					AttributeTypeCompare( t, "artificial_temperature_dependent", &type_match );
					if( type_match == _TRUE ) {
						double rho0;

						MaterialGetArtificialTemperatureDependentDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0);
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), DENSITY, artificial-temperature-dependent [%lld,%lld,%lld] = %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0);

						user->PhaseProperties.DensityLaw[p->phase_number] 	= 3; 		// artificial T-dependent density
						user->PhaseProperties.rho[p->phase_number] 			= rho0; 	// add to PhaseProperties
					}
					// others
				}
				AttributeCompare( a, "ENERGY", &attr_match );
				if( attr_match == _TRUE ) {
					// constant
					AttributeTypeCompare( t, "constant", &type_match );
					if( type_match == _TRUE ) {
						double  ThermalConductivity, HeatCapacity, RadioactiveHeat, ShearHeating;

						MaterialGetConstantEnergyParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &ThermalConductivity, &HeatCapacity, &RadioactiveHeat, &ShearHeating);
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), ENERGY, constant [%lld,%lld,%lld] = %e %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, ThermalConductivity, HeatCapacity, RadioactiveHeat, ShearHeating);


						user->PhaseProperties.T_Conductivity[p->phase_number] 		= ThermalConductivity; 	// add to PhaseProperties
						user->PhaseProperties.HeatCapacity[p->phase_number] 		= HeatCapacity; 		// add to PhaseProperties
						user->PhaseProperties.RadioactiveHeat[p->phase_number] 		= RadioactiveHeat; 		// add to PhaseProperties

					}
					// others
				}
				AttributeCompare( a, "PLASTICITY", &attr_match );
				if( attr_match == _TRUE ) {
					// constant
					AttributeTypeCompare( t, "DruckerPrager", &type_match );
					if( type_match == _TRUE ) {
						double  Cohesion, FrictionAngle;
						double Weakening_PlasticStrain_Begin,Weakening_PlasticStrain_End;
						double  CohesionAfterWeakening, FrictionAngleAfterWeakening;

						MaterialGetConstantPlasticityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &Cohesion, &FrictionAngle,
								&Weakening_PlasticStrain_Begin, &Weakening_PlasticStrain_End,
								&CohesionAfterWeakening, &FrictionAngleAfterWeakening);
						PetscPrintf(PETSC_COMM_WORLD,"    %s (id=%lld), PLASTICITY, constant [%lld,%lld,%lld] = %e %e %e %e %e %e\n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx,
								Cohesion, FrictionAngle, Weakening_PlasticStrain_Begin, Weakening_PlasticStrain_End, CohesionAfterWeakening, FrictionAngleAfterWeakening);

						user->PhaseProperties.PlasticityLaw[p->phase_number] 	= 1; 					// Drucker-Prager
						user->PhaseProperties.Cohesion[p->phase_number] 		= Cohesion; 			// add to PhaseProperties
						user->PhaseProperties.FrictionAngle[p->phase_number] 	= FrictionAngle; 		// add to PhaseProperties

						user->PhaseProperties.Weakening_PlasticStrain_Begin[p->phase_number] = Weakening_PlasticStrain_Begin; 			// add to PhaseProperties
						user->PhaseProperties.Weakening_PlasticStrain_End[p->phase_number] 	 = Weakening_PlasticStrain_End; 		  	// add to PhaseProperties
						user->PhaseProperties.CohesionAfterWeakening[p->phase_number] 	 	 = CohesionAfterWeakening; 		  						// add to PhaseProperties
						user->PhaseProperties.FrictionAngleAfterWeakening[p->phase_number] 	 = FrictionAngleAfterWeakening; 		  	// add to PhaseProperties

					}
					// others
				}
			}
		}
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// destroy attribute database right after use
	MaterialDestroy(user->PhaseMaterialProperties);

	// check empty phases
	empty_phase = PETSC_FALSE;

	for(PP = 0; PP < n_phases; PP++)
	{
		if(!check_phase[PP])
		{
			PetscPrintf(PETSC_COMM_WORLD, "Phase %lld is not initialized\n", (LLD)PP);

			empty_phase = PETSC_TRUE;
		}
	}

	ierr = PetscFree(check_phase); CHKERRQ(ierr);

	if(empty_phase == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect phase numbering");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ReadMeshSegDir"
PetscErrorCode ReadMeshSegDir(
	FILE        *fp,
	const char  *name,
	PetscScalar  beg,
	PetscScalar  end,
	PetscInt    *tncels,
	MeshSegInp  *msi,
	PetscInt     dim,
	PetscScalar  charLength)
{
	// read mesh refinement data for a direction from the input file
	// NOTE: parameter "tncels" passes negative number of segments & returns total number of cells

	PetscInt    i, jj, arsz, found;
	PetscScalar buff[3*MaxNumMeshSegs-1];

	PetscFunctionBegin;

	// check whether segments are not specified
	if((*tncels) > 0)
	{
		msi->nsegs = 0;
		PetscFunctionReturn(0);
	}

	// set number of segments
	msi->nsegs = -(*tncels);

	// read segments
	parse_GetDoubleArray(fp, name, &arsz, buff, &found);

	if(!found) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Mesh refinement segments are not specified\n");

	if(arsz != 3*msi->nsegs-1) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Incorrect number entries in the mesh refinement array\n");

	// load the data ... delimiters
	for(i = 0, jj = 0; i < msi->nsegs-1; i++, jj++) msi->delims[i] = buff[jj];

	// ... number of cells
	for(i = 0; i < msi->nsegs; i++, jj++) msi->ncells[i] = (PetscInt)buff[jj];

	// ... biases
	for(i = 0; i < msi->nsegs; i++, jj++) msi->biases[i] = buff[jj];

	// check the data ... delimiter sequence
	for(i = 1; i < msi->nsegs-1; i++)
	{
		if(msi->delims[i] <= msi->delims[i-1])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! refinement segments are unordered/overlapping\n");
		}
	}

	// ... delimiter bounds
	for(i = 0; i < msi->nsegs-1; i++)
	{
		if(msi->delims[i] <= beg || msi->delims[i] >= end)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Refinement segments out of bound\n");
		}
	}

	// ... number of cells
	for(i = 0; i < msi->nsegs; i++)
	{
		if(msi->ncells[i] <= 0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Number of cells must be non-negative\n");
		}
	}

	// ... biases
	for(i = 0; i < msi->nsegs; i++)
	{
		if(!msi->biases[i])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Bias factors must be non-zero\n");
		}
	}

	// ... compute total number of cells
	for(i = 0, (*tncels) = 0; i < msi->nsegs; i++) (*tncels) += msi->ncells[i];

	// nondimensionalize segment delimiters - after error checking
	if (dim) for(i = 0; i < msi->nsegs-1; i++) msi->delims[i] = msi->delims[i]/charLength;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// after PetscErrorCode  LaMEMInitializeMaterialProperties( UserContext *user ) in Utils.c line 430
// Set initial material properties
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitMaterialProp"
PetscErrorCode FDSTAGInitMaterialProp(UserCtx *user )
{
	PetscInt i;

	for( i=0; i<user->num_phases; i++ )
	{
		if (user->DimensionalUnits==1)
		{
			user->PhaseProperties.ViscosityLaw[i]                  = 1;     // constant viscosity
			user->PhaseProperties.mu[i]                            = 1e20;
			user->PhaseProperties.n_exponent[i]                    = 1.0;
			user->PhaseProperties.FrankKamenetskii[i]              = 0.0;
			user->PhaseProperties.Powerlaw_e0[i]                   = 1.0;
			user->PhaseProperties.A[i]                             = 1e-20;
			user->PhaseProperties.E[i]                             = 1.0;
			user->PhaseProperties.DensityLaw[i]                    = 1;     // T-dependent density
			user->PhaseProperties.rho[i]                           = 2800;
			user->PhaseProperties.Density_T0[i]                    = 273;   // Kelvin (temperature at which rho=rho0, in eq. rho=rho0*(1-alpha*(T-T0)))
			user->PhaseProperties.ThermalExpansivity[i]            = 0;     // no coupling mechanics - thermics
			user->PhaseProperties.PlasticityLaw[i]                 = 0;     // none
			user->PhaseProperties.Cohesion[i]                      = 1e100; // effectively switches off plasticity
			user->PhaseProperties.CohesionAfterWeakening[i]        = 1e100; // effectively switches off plasticity
			user->PhaseProperties.Weakening_PlasticStrain_Begin[i] = 0;
			user->PhaseProperties.Weakening_PlasticStrain_End[i]   = 0;
			user->PhaseProperties.FrictionAngle[i]                 = 0;     // effectively switches off plasticity
			user->PhaseProperties.FrictionAngleAfterWeakening[i]   = 0;     // effectively switches off plasticity
			user->PhaseProperties.ElasticShearModule[i]            = 1e100; // will make the model effectively viscous
			user->PhaseProperties.T_Conductivity[i]                = 3;     // reasonable value
			user->PhaseProperties.HeatCapacity[i]                  = 1050;  // reasonable value
			user->PhaseProperties.RadioactiveHeat[i]               = 0;
		}
		else
		{
			user->PhaseProperties.ViscosityLaw[i]                  = 1;     // constant viscosity
			user->PhaseProperties.mu[i]                            = 1;
			user->PhaseProperties.n_exponent[i]                    = 1.0;
			user->PhaseProperties.FrankKamenetskii[i]              = 0.0;
			user->PhaseProperties.Powerlaw_e0[i]                   = 1.0;
			user->PhaseProperties.A[i]                             = 1.0;
			user->PhaseProperties.E[i]                             = 1.0;
			user->PhaseProperties.DensityLaw[i]                    = 1;     // T-dependent density
			user->PhaseProperties.rho[i]                           = 1;
			user->PhaseProperties.Density_T0[i]                    = 0;     // Kelvin (temperature at which rho=rho0, in eq. rho=rho0*(1-alpha*(T-T0)))
			user->PhaseProperties.ThermalExpansivity[i]            = 0;     // no coupling mechanics - thermics
			user->PhaseProperties.PlasticityLaw[i]                 = 0;     // none
			user->PhaseProperties.Cohesion[i]                      = 1e100; // effectively switches off plasticity
			user->PhaseProperties.CohesionAfterWeakening[i]        = 1e100; // effectively switches off plasticity
			user->PhaseProperties.Weakening_PlasticStrain_Begin[i] = 0;
			user->PhaseProperties.Weakening_PlasticStrain_End[i]   = 0;
			user->PhaseProperties.FrictionAngle[i]                 = 0;     // effectively switches off plasticity
			user->PhaseProperties.FrictionAngleAfterWeakening[i]   = 0;     // effectively switches off plasticity
			user->PhaseProperties.ElasticShearModule[i]            = 1e100; // will make the model effectively viscous
			user->PhaseProperties.T_Conductivity[i]                = 1;     // reasonable value
			user->PhaseProperties.HeatCapacity[i]                  = 1;     // reasonable value
			user->PhaseProperties.RadioactiveHeat[i]               = 0;
		}
	}

	if (user->DimensionalUnits==1){
		// Add a higher density and viscosity to phase 1
		user->PhaseProperties.rho[1] = 3000;
		user->PhaseProperties.mu[1]  = 1e25;
	}
	else {
		// Add a higher density and viscosity to phase 1
		user->PhaseProperties.rho[1] = 2;
		user->PhaseProperties.mu[1]  = 1e3;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// change values from the command line
#undef __FUNCT__
#define __FUNCT__ "FDSTAGReadCommLine"
PetscErrorCode FDSTAGReadCommLine(UserCtx *user )
{
	PetscInt       i, nel_array[3], nel_input_max;
	PetscBool      found,flg;
	char           setup_name[PETSC_MAX_PATH_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscOptionsGetReal(PETSC_NULL,"-W"      , &user->W      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-L"      , &user->L      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-H"      , &user->H      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-y_front", &user->y_front, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-nel_x",   &user->nel_x,   PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_y",   &user->nel_y,   PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_z",   &user->nel_z,   PETSC_NULL);

	// alternative: specify the # of elements as -nel 8,16,32  which gives nel_x=8, nel_y=16, nel_z=32
	nel_input_max=3;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-nel", nel_array, &nel_input_max, &found); CHKERRQ(ierr);

	if (found==PETSC_TRUE) {
		user->nel_x = nel_array[0];
		user->nel_y = nel_array[1];
		user->nel_z = nel_array[2];
	}
	user->nnode_x = user->nel_x + 1;
	user->nnode_y = user->nel_y + 1;
	user->nnode_z = user->nel_z + 1;

	PetscOptionsGetInt(PETSC_NULL ,"-time_end",       &user->time_end,       PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-save_timesteps", &user->save_timesteps, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-CFL",            &user->CFL,            PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt",             &user->dt,             PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt_max",         &user->dt_max,         PETSC_NULL);

	// FDSTAG Canonical Model Setup
	PetscOptionsGetString(PETSC_NULL,"-msetup", setup_name, PETSC_MAX_PATH_LEN, &found);
	if(found == PETSC_TRUE)
	{	if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"ERROR! Incorrect model setup: %s", setup_name);
	}

	// boundary conditions
	PetscOptionsGetInt(PETSC_NULL, "-BC.InternalBound", &user->BC.InternalBound, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.UpperBound",    &user->BC.UpperBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.LowerBound",    &user->BC.LowerBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.LeftBound",     &user->BC.LeftBound,     PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.RightBound",    &user->BC.RightBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.FrontBound",    &user->BC.FrontBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.BackBound",     &user->BC.BackBound,     PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vy_front",      &user->BC.Vy_front,      PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vy_back",       &user->BC.Vy_back,       PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vz_top",        &user->BC.Vz_top,        PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vz_bot",        &user->BC.Vz_bot,        PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vx_left",       &user->BC.Vx_left,       PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vx_right",      &user->BC.Vx_right,      PETSC_NULL); // y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL,"-BC.Exx",           &user->BC.Exx,           PETSC_NULL); // Exx background strain-rate
	PetscOptionsGetReal(PETSC_NULL,"-BC.Eyy",           &user->BC.Eyy,           PETSC_NULL); // Eyy background strain-rate

	// number of markers
	PetscOptionsGetInt(PETSC_NULL,"-NumPartX",          &user->NumPartX,         PETSC_NULL); // # of tracers per cell in x-direction
	PetscOptionsGetInt(PETSC_NULL,"-NumPartY",          &user->NumPartY,         PETSC_NULL); // # of tracers per cell in y-direction
	PetscOptionsGetInt(PETSC_NULL,"-NumPartZ",          &user->NumPartZ,         PETSC_NULL); // # of tracers per cell in z-direction

	// flags
	PetscOptionsGetInt(PETSC_NULL,"-restart",          &user->restart,          PETSC_NULL); // # restart a simulation if possible?
	PetscOptionsGetInt(PETSC_NULL,"-save_breakpoints", &user->save_breakpoints, PETSC_NULL); // after how many steps do we create a breakpoint file?
	PetscOptionsGetInt(PETSC_NULL,"-SaveParticles",    &user->SaveParticles,    PETSC_NULL); // save particles to disk?

	PetscOptionsGetBool(PETSC_NULL,"-SavePartitioning",&user->SavePartitioning, PETSC_NULL);
	PetscOptionsGetBool(PETSC_NULL,"-SkipStokesSolver",&user->SkipStokesSolver, PETSC_NULL);

	// optimization
	PetscOptionsGetReal(PETSC_NULL,"-LowerViscosityCutoff", &user->LowerViscosityCutoff, PETSC_NULL); // lower viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL,"-UpperViscosityCutoff", &user->UpperViscosityCutoff, PETSC_NULL); // upper viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL,"-InitViscosity",        &user->InitViscosity,        PETSC_NULL); // upper viscosity cutoff

	// gravity
	PetscOptionsGetReal(PETSC_NULL ,"-GravityAngle",   &user->GravityAngle,     PETSC_NULL); // Gravity angle in x-z plane

	// pushing boundary conditions related parameters
	PetscOptionsGetInt(PETSC_NULL,"-AddPushing"                 , &user->AddPushing                 , PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.num_changes"        , &user->Pushing.num_changes        , PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.L_block"           , &user->Pushing.L_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.W_block"           , &user->Pushing.W_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.H_block"           , &user->Pushing.H_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.x_center_block"    , &user->Pushing.x_center_block     , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.y_center_block"    , &user->Pushing.y_center_block     , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.z_center_block"    , &user->Pushing.z_center_block     , PETSC_NULL);

	// pushing - array variables
	char matprop_opt[PETSC_MAX_PATH_LEN];
	flg = PETSC_FALSE;

	for(i=0;i<user->Pushing.num_changes;i++){
		// v_push in cm/year
		sprintf(matprop_opt,"-Pushing.V_push_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.V_push[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    V_push[%lld] = %g \n",(LLD)i,user->Pushing.V_push[i]);

		// rate of rotation deg/yr
		sprintf(matprop_opt,"-Pushing.omega_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.omega[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Omega[%lld] = %g \n",(LLD)i,user->Pushing.omega[i]);

		// advect block or not for each time segment
		sprintf(matprop_opt,"-Pushing.coord_advect_%lld",(LLD)i);
		ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.coord_advect[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Coord_advect[%lld] = %g \n",(LLD)i,user->Pushing.coord_advect[i]);

		// direction of pushing
				sprintf(matprop_opt,"-Pushing.dir_%lld",(LLD)i);
				ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.dir[i], &flg); CHKERRQ(ierr);
				if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Direction[%lld] = %g \n",(LLD)i,user->Pushing.dir[i]);
	}

	for(i=0;i<user->Pushing.num_changes+1;i++){
		// time intervals is a num_changes+1 array
		sprintf(matprop_opt,"-Pushing.time_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.time[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Time[%lld] = %g \n",(LLD)i,user->Pushing.time[i]);
	}

	// use LaMEM canonical
	PetscOptionsGetBool(PETSC_NULL,"-use_fdstag_canonical", &user->use_fdstag_canonical, PETSC_NULL );

	// get material properties from command line - de-activated
	//ierr = GetMaterialPropertiesFromCommandLine(user);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

