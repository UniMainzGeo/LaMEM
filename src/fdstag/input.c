//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "NonDimensionalisation.h"
#include "Parsing.h"
#include "Utils.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "input.h"
//---------------------------------------------------------------------------
// * add default solver options
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitCode"
PetscErrorCode FDSTAGInitCode(JacRes *jr, UserContext *user)
{
	// Set default code parameters and read input file, if required

	PetscMPIInt    size;
	PetscErrorCode ierr;
	PetscInt       i, nel_array[3], nel_input_max;
    PetscInt       n_int;
	PetscScalar    SecYear;
	PetscBool      found,flg;
	char          *all_options;
	char           setup_name[PETSC_MAX_PATH_LEN];

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED
	PetscOptionsGetAll( &all_options ); /* copy all command line args */

	SecYear 			=	3600*24*365.25;		// seconds per year

	/* Set default values */
	user->W 		  	= 		1.0;
	user->L 			= 		1.0;
	user->H 		 	= 		1.0;
	user->y_front 		= 		0.0;
	user->ampl2D 	  	= 		0.0;
	user->ampl3D  		= 		1e-2;
	user->amplNoise 	= 		0.0;
	user->Hinterface 	=		0.5;			// average interface height for diapir setup [0-1]
	user->mumax			= 		1.0;            // REMOVE?

	// set default grid sizes
	user->nel_x  		=       8;
	user->nel_y  		=       8;
	user->nel_z			=       8;
	user->refinex     	=     	2;              // REMOVE
	user->refiney   	=     	2;              // REMOVE
	user->refinez  		=  		2;              // REMOVE

	user->time_start 	=  		0.0;
	user->time_end 		= 		1.0;
	user->time_end_temp	= 		1;
	user->save_timesteps=		1;
	user->time		  	=   	0.0;
	user->CFL 		 	= 		0.5;
	user->Temp_top  		= 	0.0;
	user->Temp_bottom 		=	1.0;
	user->temp_initialize	=	0.0;            // REMOVE
	user->DimensionalUnits 	= 	0;
	user->Setup.Model 		=	2;              // // REMOVE 0-Diapir, 1-Single Layer Fold, 2-Falling block, 3-Particle file (defined below), 4-multilayer detachment fold, 6-subduction w/sticky air
	user->Gravity   		= 	1.0;
	user->GasConstant  		= 	1.0;            // REMOVE
	user->BC.Vy_front		=	0;
	user->BC.Vy_back		=	0;
	user->BC.Vz_bot 		=	0;
	user->BC.Vz_top   		=	0;
	user->BC.Vx_left		=	0;
	user->BC.Vx_right		=	0;
	user->BC.Exx			=	0;
	user->BC.Eyy		 	=	0;
	sprintf(user->OutputFile, 		"FallingBlock3D");
	sprintf(user->ParticleFilename, "ParticlesInput.dat");
	sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");

	// FDSTAG Canonical Model Setup
	user->msetup            = BLOCK;

	user->MatlabOutputFiles	=	1;		// write MATLAB output
	user->VTKOutputFiles	=	1;		// write VTK output
	user->AVDPhaseViewer  	= 	0;
	user->SavePartitioning  = 	PETSC_FALSE;		// write partitioning
	user->x_left 	  		= 	-0.5*(user->W)*0.0;
	user->y_front	  		= 	-0.5*(user->L)*0.0;
	user->z_bot 			= 	-0.5*(user->H)*0.0;

	// linear solver settings
	user->SkipStokesSolver     = PETSC_FALSE;
	user->StokesSolver         = 1; // 1 - Powell-Hesteness iterations; 2 - Schur Complement Reduction; 3 - Fully Coupled Solver; 4 - MatVec Test;
	user->VelocitySolver       = 1; // 0 - User defined; 1 - Direct (MUMPS); 2 - Galerkin geometric multigrid; 3 - Fieldsplit + Algebraic Multigrid (ML), 4 - GCR & HYPRE preconditioner on full velocity block
	user->VelocityTest         = PETSC_FALSE; // Request to perform single velocity solve for test purposes & quit
	user->ScaleSystem          = PETSC_FALSE; // Request to scale linear system before solution
	user->use_fdstag_canonical = PETSC_TRUE;  // request native staggered grid discretization

	user->internalBC_coord			= 0.0;
	user->internalBC_frontel		= 0.0;
	user->internalBC_backel			= 0.0;
	user->internalBC_node		    = 0.0;
	user->zdepth_BC_el			    = 0.0;  // REMOVE
	user->zdepth_BC_node		    = 0.0;  // REMOVE
	user->BC.InternalBound			=   0;	// 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.UpperBound				=   1;	// 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.LowerBound				=	1;	// 0-free surface, 1-free slip					, 2-no-slip
	user->BC.LeftBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no slip
	user->BC.RightBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.FrontBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.BackBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->NumParticlesToStartInjection = 0;  	// how many particles should an element have before we start injecting particles?
	user->ParticleInjectionPhase       = 0;  	// which is the phase of the injected particle?
	user->ParticleInput = 1;					// 0-do not use particles to track phases; 1-do use particles to track phases
	user->LoadInitialParticlesFromDisc = 0;		// Read the initial particles from disc
	user->remesh = 0;                           // REMOVE
	user->CriticalDiagonalRatio = 0.55;			// save range is [0.4-1.0] // REMOVE

	// Check if we are performing a benchmark with a build-in analytical benchmark
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_SolCx",   		  		&user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_FallingBlock",   		&user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_ArcTanFallingBlock",   &user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_VerticalDensityCollumn",   &user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);

	user->GridAdvectionMethod  = 0; 	 // 0-Fully Eulerian, 1-Fully Lagrangian, 2-ALE with remeshing @ surface layer      REMOVE
	user->EulerianAfterTimestep=-1;		//	if >0, you can specify after which timestep to switch to eulerian mode          REMOVE
	user->FactorSurfaceLayer   = 0.2; 	 // how thick is the surface layer?                                                 REMOVE
	user->num_subdt			   =  10;    // How many sub-timestep iterations if an ALE mode is selected?                    REMOVE

    user->SaveParticles 	   =  0;	 // Save particles or not?
	user->num_phase_transitions= 0;         // REMOVE
	user->InitialMeshFromFile  = 0;		 // In case you want to read an initial mesh from file
	user->InitialErosionSurfaceFromFile = 0; //
	user->InitialMantleLevel   = 10;     // for cases in which a lithosphere is modeled with hand-set crustal thickness     REMOVE
    
	user->dt_max               = 1e6;    // maximum timestep
	user->dt_temp              = 1e6;    // timestep to initialize temperature
	user->dt				   = 1e-3;	 // initial timestep
	user->Vx_Front             = 8.0;      // x-velocity at front boundary if BC==3 @ this boundary
	user->Vx_Back              = 0.0;      // x-velocity at back boundary if BC==3 @ this boundary
	user->Vy_Front             = 4.0;      // y-velocity at front boundary if BC==4 @ this boundary
	user->Vy_Back              = 0.0;      // y-velocity at back boundary if BC==4 @ this boundary
	user->Vy_Partx             = 0.5;    // Apply velocity condition at part of the x-domain, if BC=4
	user->Vy_Partz             = 0.0;      // Apply velocity condition at part of the x-domain, if BC=4

	user->MaximumSurfaceAngle  = 80.0;	 // maximum surface angle allowed                                               REMOVE
	user->MuMeanMethod         = 1;      // 0-compute mat. props @ integration points, 1-Arith. element-average, 2-Geom. element-average, 3-Harm. element-average REMOVE
    
	user->NumPartX			   = 2;      // Number of particles per cell in x-direction
	user->NumPartY			   = 2;      // Number of particles per cell in y-direction
	user->NumPartZ			   = 2;      // Number of particles per cell in z-direction

	user->restart			   = 1; 	 // Are we restarting a simulation?
	user->save_breakpoints	   = 10;		 // After how many steps do we make a breakpoint?
	if (user->AnalyticalBenchmark){
		user->restart 		   = 0;		// don't restart and don't create breakpoint files if we perform an analytical benchmark
		user->save_breakpoints = 0;
	}
	user->break_point_number   = 0;		 // The number of the breakpoint file
	user->incr_breakpoints 	   = 100;		 // After how many steps should an incremental breakpointfile be created?
	user->fileno   			   = 0;

	user->LowerViscosityCutoff  = 1.0;		// lowermost viscosity cutoff in model
	user->UpperViscosityCutoff  = 1e30;		// uppermost viscosity cutoff in model

	user->num_phases 		   = 2;  // the default # of phases. In case we set props from the command line, and we use a folding setup combined with FDSTAG we might have to do something smarter (check the setup and increase this)

	user->GravityAngle 		   	= 	0.0;		// angle of gravity with z-axis (can be changed in the x-z plane)
	user->FSSA					=	0.0;		// Free Surface Stabilization Algorithm parameter used in FDSTAG, to stabilize an internal free surface [should be between 0-1]

	user->ArtTemp               = PETSC_FALSE;

	// Default erosion parameters
	user->ApplyErosion 		  	=  0; 	 	// 0-no; 1-cascade
	user->SurfaceAngle 			=  10.0;      // Angle that the surface makes with horizontal    [degree]
	user->SurfaceNoiseAmplitude =  100;  	// Amplitude of noise on free surface 				[m     ]
	user->fluvial_erosion 		=  1.6e-2; 	// Fluvial erosion rate
	user->diffusion_erosion 	=  2.0e-2; 	// Diffusion erosion rate
	user->baselevelx0 			=  1;      	//	Outlet at front [if 1; no outlet=0]
	user->baselevelx1 			=  0;      	//
	user->baselevely0 			=  0;     	//
	user->baselevely1 			=  0;      	//
	user->ErosionParameters.UseInternalFreeSurface 	= 0;	// don't use by default
	user->ErosionParameters.StickyAirPhase 			= 1;	// sticky air phase
	user->ErosionParameters.ErosionModel 			= 0;	// none by default
	user->ErosionParameters.SedimentationModel 		= 0;	// no sedimentation by default; 1-constant rate sedimentation

	//Default FD erosion code parameters
	user->ErosionParameters.FE_ErosionCode.ResolutionFactorX 		=	2;				// how much larger is the FD erosion code resolution compared to the mechanical code?
	user->ErosionParameters.FE_ErosionCode.ResolutionFactorY 		=	2;				// how much larger is the FD erosion code resolution compared to the mechanical code?
	user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m 	=	1;				// what is the amplitude of the random noise on the erosion surface [in meters]?
	user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m 	=	100;			// By how much meters is the right side of the eroded surface uplifed [in meters]?
	user->ErosionParameters.FE_ErosionCode.dt 						=	100.0*SecYear;	// ideal erosion timestep
	user->ErosionParameters.FE_ErosionCode.rain_m_year 				=	0.3;			// rain in m/year
	user->ErosionParameters.FE_ErosionCode.k0 						=	3.2e-12;		// base erodability
	user->ErosionParameters.FE_ErosionCode.c 						=	1.0;			// prefactor for fluvial discharge
	user->ErosionParameters.FE_ErosionCode.n 						=	2.0;			// streampower exponent
    user->ErosionParameters.FE_ErosionCode.BC                       =   1;
    user->ErosionParameters.FE_ErosionCode.fill_lake                =   0;
    user->ErosionParameters.FE_ErosionCode.mode_river               =   0;
    user->ErosionParameters.FE_ErosionCode.nbre_river               =   1;
    user->ErosionParameters.FE_ErosionCode.rain_river_year 			=	0.0;			// rain in m/year


	// Default Output parameters
	user->Output.velocity				=   1;
	user->Output.temperature			=   1;
	user->Output.surface_topography		=   1;
	user->Output.bottom_topography		=   1;
	user->Output.quadrature				=   1;

	// Default Pushing BC parameters
	user->AddPushing					=	0;
	user->Pushing.reset_pushing_coord	= 	0;
	user->Pushing.theta					= 	0.0;

	// Default gravity field parameters
	user->GravityField.GetIt		=	0;
	user->GravityField.SaveDebug	=	0;
	user->GravityField.SaveRef		=	0;
	user->GravityField.SaveVTK		=	0;
	user->GravityField.UseNumerics	= 	PETSC_FALSE;
	user->GravityField.UseAnalytics	= 	PETSC_FALSE;
	user->GravityField.num_intp 	=	2;
	user->GravityField.survey_nx 	=	11;
	user->GravityField.survey_ny 	=	11;
	user->GravityField.survey_xs 	=	0.0;
	user->GravityField.survey_xm 	=	1.0;
	user->GravityField.survey_ys 	=	0.0;
	user->GravityField.survey_ym 	=	1.0;
	user->GravityField.survey_z		=	0.0;
	user->GravityField.ReferenceDensity = 2700.0;
	user->GravityField.StdDev	=	1.0;
	sprintf(user->GravityField.RefDatFile2load, "ReferenceData/GravityField_REF.bin");
	user->GravityField.LithColNum      =   1;
	user->GravityField.LithColDens[0]  =   2670.0;
	user->GravityField.LithColDepth[0] =   30.0e3;

	// Default isostasy parameters
	user->Isostasy.GetIt            = 0;
	user->Isostasy.SaveRef          = 0;
	user->Isostasy.corr_topo        = 0.0;
	user->Isostasy.ref_rho          = 3215.0; // kg/m3
	user->Isostasy.ref_xi          = 0; // []
	user->Isostasy.ref_yi          = 0; // []

	// Default surface velocity field parameters
	user->SurfVelField.GetIt		=	0;
	user->SurfVelField.SaveRef		=	0;
	user->SurfVelField.VxStdDev		=	1.0;
	user->SurfVelField.VyStdDev		=	1.0;
	user->SurfVelField.VzStdDev		=	1.0;
	sprintf(user->GravityField.RefDatFile2load, "ReferenceData/SurfVelField_REF.bin");


	// Default optimisation parameters
	user->Optimisation.GetIt		=	0;
	user->Optimisation.MisfitGravity=	0.0;
	user->Optimisation.MisfitSurfVel=	0.0;
	user->Optimisation.MisfitTiso   =	0.0;
	user->Optimisation.mpi_group_id	=	0;


	user->DA_Materials = PETSC_NULL;
	user->Materials    = PETSC_NULL;

	/* Initialize material properties */
	LaMEMInitializeMaterialProperties(user);

	/* Read an input file if required */
	if (user->InputParamFile){
		//	ReadInputFile(user);
		FDSTAGReadInputFile(jr, user);		// use the parser, to read the input file
	}

	/* Change values @ command prompt */

	PetscOptionsGetReal(PETSC_NULL,"-W"      ,	&user->W 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-L"      ,	&user->L 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-H"      ,	&user->H 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-He"      ,	&user->H 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-y_front" ,	&user->y_front 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-ampl2D" ,	&user->ampl2D 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-ampl3D" ,	&user->ampl3D 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-amplNoise",&user->amplNoise	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-mumax"  ,	&user->mumax 		, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-nel_x",	&user->nel_x 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_y",	&user->nel_y 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_z",	&user->nel_z 		, PETSC_NULL);

	nel_input_max=3;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-nel", nel_array, &nel_input_max, &found); CHKERRQ(ierr);
	// this allows us to also specify the # of elements as -nel 8,16,32  which gives nel_x=8, nel_y=16, nel_z=32
	if (found==PETSC_TRUE) {
		user->nel_x = nel_array[0];
		user->nel_y = nel_array[1];
		user->nel_z = nel_array[2];
	}

	PetscOptionsGetInt(PETSC_NULL ,"-refinex",		&user->refinex		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-refiney",		&user->refiney		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-refinez",		&user->refinez		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-time_end",	&user->time_end 	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-time_end_temp",	&user->time_end_temp 	, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-save_timesteps",	&user->save_timesteps 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-CFL"    ,		&user->CFL 			, PETSC_NULL);

	PetscOptionsGetReal(PETSC_NULL,"-dt"    ,		&user->dt 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt_max"    ,	&user->dt_max 			, PETSC_NULL);

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
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"ERROR! Incorrect model setup: %s", setup_name);
	}

	PetscOptionsGetInt(PETSC_NULL ,"-Setup.Model",	&user->Setup.Model 	, PETSC_NULL);	// 		0-diapir, 1-single layer folding
	PetscOptionsGetInt(PETSC_NULL ,"-GridAdvectionMethod",	&user->GridAdvectionMethod 	, PETSC_NULL);	// 		0 - eulerian, 2-ALE, 1-Lagrangian
	PetscOptionsGetBool( PETSC_NULL,"-SkipStokesSolver",&user->SkipStokesSolver,PETSC_NULL );

	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_coord",&user->internalBC_coord	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_frontel",&user->internalBC_frontel	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_backel",&user->internalBC_backel	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_node",&user->internalBC_node	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-zdepth_BC_el",&user->zdepth_BC_el	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-zdepth_BC_node",&user->zdepth_BC_node	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.InternalBound",&user->BC.InternalBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.UpperBound",&user->BC.UpperBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.LowerBound",&user->BC.LowerBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.LeftBound",	&user->BC.LeftBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.RightBound",&user->BC.RightBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.FrontBound",&user->BC.FrontBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.BackBound" ,&user->BC.BackBound		, PETSC_NULL);	//
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vy_front",	&user->BC.Vy_front 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vy_back",	&user->BC.Vy_back 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vz_top",	&user->BC.Vz_top 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vz_bot",	&user->BC.Vz_bot 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vx_left",	&user->BC.Vx_left 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vx_right",	&user->BC.Vx_right 	, PETSC_NULL);	//  	y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL ,"-BC.Exx",			&user->BC.Exx 		, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Eyy",			&user->BC.Eyy 		, PETSC_NULL);	//  	y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL ,"-Vx_Front",	&user->Vx_Front		, PETSC_NULL);	// Vx_Front in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vx_Back",		&user->Vx_Back		, PETSC_NULL);	// Vx_Back in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Front",	&user->Vx_Front		, PETSC_NULL);	// Vy_Front in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Back",		&user->Vx_Back		, PETSC_NULL);	// Vy_Back in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Partx",	&user->Vy_Partx		, PETSC_NULL);	// Vy_Partx for BC=4
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Partz",	&user->Vy_Partz		, PETSC_NULL);	// Vy_Partz for BC=4
	PetscOptionsGetInt(PETSC_NULL ,"-MuMeanMethod",	&user->MuMeanMethod	, PETSC_NULL);  // Specify element-averaging of effective material properties.
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartX",	&user->NumPartX	, PETSC_NULL);  		//	# of tracers per cell in x-direction
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartY",	&user->NumPartY	, PETSC_NULL);  		//	# of tracers per cell in y-direction
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartZ",	&user->NumPartZ	, PETSC_NULL);  		//	# of tracers per cell in z-direction
	PetscOptionsGetInt(PETSC_NULL ,"-restart",	&user->restart	, PETSC_NULL);  		//	# restart a simulation if possible?
	PetscOptionsGetInt(PETSC_NULL ,"-fileno",	&user->fileno	, PETSC_NULL);  		//	# restart a simulation with file no ...
	PetscOptionsGetInt(PETSC_NULL ,"-save_breakpoints",	&user->save_breakpoints	, PETSC_NULL);  	//	After how many steps do we create a breakpoint file?
	PetscOptionsGetInt(PETSC_NULL ,"-ApplyErosion",	      &user->ApplyErosion		, PETSC_NULL);  //	Apply erosion or not?
	PetscOptionsGetReal(PETSC_NULL ,"-fluvial_erosion",	  &user->fluvial_erosion		, PETSC_NULL);	//  fluvial erosion rate
	PetscOptionsGetReal(PETSC_NULL ,"-diffusion_erosion", &user->diffusion_erosion			, PETSC_NULL);	//  diffusion erosion rate
	PetscOptionsGetReal(PETSC_NULL ,"-SurfaceNoiseAmplitude", &user->SurfaceNoiseAmplitude	, PETSC_NULL);	// amplitude of noise @ surface [m]
	PetscOptionsGetReal(PETSC_NULL ,"-Hinterface", &user->Hinterface						, PETSC_NULL);	// Horizontal interface for diapir setup
	PetscOptionsGetInt(PETSC_NULL ,"-InitialErosionSurfaceFromFile", &user->InitialErosionSurfaceFromFile, PETSC_NULL);	// load initial erosion surface from file

	PetscOptionsGetReal(PETSC_NULL ,"-NonlinearIterationsAccuracy", &user->NonlinearIterationsAccuracy	, PETSC_NULL);		// accuracy of nonlinear iterations
	PetscOptionsGetInt(PETSC_NULL ,"-MaxNonlinearIterations",		&user->MaxNonlinearIterations		, PETSC_NULL);  	//	maximum number of nonlinear iterations

	PetscOptionsGetReal(PETSC_NULL ,"-LowerViscosityCutoff", &user->LowerViscosityCutoff	, PETSC_NULL);		// lower viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL ,"-UpperViscosityCutoff", &user->UpperViscosityCutoff	, PETSC_NULL);		// upper viscosity cutoff

	PetscOptionsGetReal(PETSC_NULL ,"-GravityAngle", &user->GravityAngle	, PETSC_NULL);		// Gravity angle in x-z plane


	PetscOptionsGetInt(PETSC_NULL ,"-incr_breakpoints",		&user->incr_breakpoints		, PETSC_NULL);  	//	save incr_breakpoints breakpoints
	PetscOptionsGetReal(PETSC_NULL ,"-MaximumSurfaceAngle", &user->MaximumSurfaceAngle	, PETSC_NULL);		// MaximumSurfaceAngle

	PetscOptionsGetInt(PETSC_NULL ,"-ParticleInjectionPhase",		&user->ParticleInjectionPhase	, PETSC_NULL);  		//	which phase do we inject?
	PetscOptionsGetInt(PETSC_NULL ,"-MatlabOutputFiles",			&user->MatlabOutputFiles		, PETSC_NULL);  		//	write matlab output or not?
	PetscOptionsGetInt(PETSC_NULL ,"-VTKOutputFiles",				&user->VTKOutputFiles			, PETSC_NULL);  		//	write VTK output or not?
	PetscOptionsGetInt(PETSC_NULL ,"-AVDPhaseViewer",				&user->AVDPhaseViewer			, PETSC_NULL);  		//	write AVDPhase output or not?
    PetscOptionsGetBool(PETSC_NULL ,"-SavePartitioning",			&user->SavePartitioning        	, PETSC_NULL);


	PetscOptionsGetInt(PETSC_NULL ,"-LoadInitialParticlesFromDisc",	&user->LoadInitialParticlesFromDisc			, PETSC_NULL);  		//	Load initial particles from file

	PetscOptionsGetReal(PETSC_NULL ,"-CriticalDiagonalRatio", 		&user->CriticalDiagonalRatio 	, PETSC_NULL);		// CriticalDiagonalRatio to induce remeshing if LaMEM is run in a lagrangian mode
	PetscOptionsGetInt(PETSC_NULL ,"-SaveParticles",				&user->SaveParticles			, PETSC_NULL);  		//	save particles to disk?
	PetscOptionsGetInt(PETSC_NULL ,"-EulerianAfterTimestep",				&user->EulerianAfterTimestep			, PETSC_NULL);  		//	switch to Eulerian mode after timestep ?? - decativated if <0 (default)

	PetscOptionsGetInt(PETSC_NULL  ,"-UseInternalFreeSurface",		&user->ErosionParameters.UseInternalFreeSurface	, PETSC_NULL);  	//	use the internal free surface
	PetscOptionsGetInt(PETSC_NULL  ,"-ErosionModel",				&user->ErosionParameters.ErosionModel	, PETSC_NULL);  			//	which erosion model do we employ? [0-none; 1-fast]


	/*set the FD erosion parameters from the command-line */
	PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.ResolutionFactorX",			&user->ErosionParameters.FE_ErosionCode.ResolutionFactorX		, PETSC_NULL);  			// how much larger is the FD erosion code resolution compared to the mechanical code?
	PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.ResolutionFactorY",			&user->ErosionParameters.FE_ErosionCode.ResolutionFactorY		, PETSC_NULL);  			// how much larger is the FD erosion code resolution compared to the mechanical code?
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.InitialRandomNoise_m", 	&user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.InitialUpliftedSide_m", 	&user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.dt", 						&user->ErosionParameters.FE_ErosionCode.dt 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.rain_m_year", 				&user->ErosionParameters.FE_ErosionCode.rain_m_year 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.k0", 						&user->ErosionParameters.FE_ErosionCode.k0 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.c", 						&user->ErosionParameters.FE_ErosionCode.c 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.n", 						&user->ErosionParameters.FE_ErosionCode.n 						, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.BC",                       &user->ErosionParameters.FE_ErosionCode.BC                      , PETSC_NULL);
    n_int = 100;
    found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.mode_river",               &user->ErosionParameters.FE_ErosionCode.mode_river              , &found);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.fill_lake",                &user->ErosionParameters.FE_ErosionCode.fill_lake               , PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.nbre_river",               &user->ErosionParameters.FE_ErosionCode.nbre_river              , PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.rain_river_year", 			&user->ErosionParameters.FE_ErosionCode.rain_river_year 		, PETSC_NULL);

    PetscOptionsGetRealArray(PETSC_NULL,"-FE_ErosionCode.location_river",user->ErosionParameters.FE_ErosionCode.location_river,&n_int, &flg);
    if (found){
    	/* print info if we select river parameters */

    	printf("==========================================================================================\n");
    	printf("in file utils.c \n");
    	printf("mode_river:%d, nbre_river:%d \n", user->ErosionParameters.FE_ErosionCode.mode_river, user->ErosionParameters.FE_ErosionCode.nbre_river);


    	for( i=0; i<user->ErosionParameters.FE_ErosionCode.nbre_river; i++ ) {
    		printf("location_river[%d]:%f\n",i,user->ErosionParameters.FE_ErosionCode.location_river[i]);
    	}

    	printf("==========================================================================================\n");
    }



	PetscOptionsGetInt(PETSC_NULL  ,"-StickyAirPhase",				&user->ErosionParameters.StickyAirPhase			, PETSC_NULL);  														 //	which phase is sticky air?
	PetscOptionsGetReal(PETSC_NULL ,"-FSSA", 						&user->FSSA 									, PETSC_NULL);		// FSSA parameter [should be between 0-1]
	PetscOptionsGetBool(PETSC_NULL, "-ArtificialTemperature", 	    &user->ArtTemp,    PETSC_NULL );

	PetscOptionsGetInt(PETSC_NULL  ,"-SedimentationModel",			&user->ErosionParameters.SedimentationModel	, PETSC_NULL);  			//	which sedimentation model do we employ? [0-none; 1-constant rate]
	PetscOptionsGetReal(PETSC_NULL ,"-InitialFreeSurfaceHeight", 	&user->ErosionParameters.InitialFreeSurfaceHeight 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-SedimentationRate_cmYr", 		&user->ErosionParameters.SedimentationRate_cmYr 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-SedimentLayerThicknessYears", &user->ErosionParameters.SedimentLayerThicknessYears 	, PETSC_NULL);  // we change the sediment phase every ?? years (mainly a visualization issue)
	PetscOptionsGetInt(PETSC_NULL  ,"-PhaseFirstSedimentedLayer",	&user->ErosionParameters.PhaseFirstSedimentedLayer		, PETSC_NULL);  // which is the first sediment layer?
	PetscOptionsGetInt(PETSC_NULL  ,"-PhaseLastSedimentedLayer",	&user->ErosionParameters.PhaseLastSedimentedLayer		, PETSC_NULL);  // which is the last sediment layer? After this, we sediment the first one again

	// --- SurfaceVelocity related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_SurfVelField", 			&user->SurfVelField.GetIt, 		PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-SurfVelField_SaveRef", 		&user->SurfVelField.SaveRef,	PETSC_NULL );
	PetscOptionsGetString(PETSC_NULL,"-SurfVelField_RefDatFile",  *(&user->SurfVelField.RefDatFile2load),PETSC_MAX_PATH_LEN,PETSC_NULL);
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VxStdDev",         &user->SurfVelField.VxStdDev,PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VyStdDev",         &user->SurfVelField.VyStdDev,PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VzStdDev",         &user->SurfVelField.VzStdDev,PETSC_NULL );
	 // --- Optimization related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_Misfit", 					&user->Optimisation.GetIt,		PETSC_NULL );

	// --- Gravity related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_GravityField",			&user->GravityField.GetIt, 	   	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveRef", 		&user->GravityField.SaveRef,	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveDebug",		&user->GravityField.SaveDebug,	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveVTK",		&user->GravityField.SaveVTK, 	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_survey_nx",		&user->GravityField.survey_nx, 	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_survey_ny", 		&user->GravityField.survey_ny, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_xs", 		&user->GravityField.survey_xs, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_xm", 		&user->GravityField.survey_xm, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_ys",	 	&user->GravityField.survey_ys, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_ym", 		&user->GravityField.survey_ym, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_z", 		&user->GravityField.survey_z,  	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_ReferenceDensity",&user->GravityField.ReferenceDensity,PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_num_intp", 		&user->GravityField.num_intp,  	PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL,"-GravityField_UseNumerics",	&user->GravityField.UseNumerics,PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL,"-GravityField_UseAnalytics",	&user->GravityField.UseAnalytics,PETSC_NULL );
	PetscOptionsGetString(PETSC_NULL,"-GravityField_RefDatFile",  *(&user->GravityField.RefDatFile2load),PETSC_MAX_PATH_LEN,PETSC_NULL);
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_StdDev",         &user->GravityField.StdDev,PETSC_NULL );
	ierr = GetLithColumnFromCommandLine(user);CHKERRQ(ierr);

	// --- Isostasy related parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-ComputeAiryIsostasy", 		&user->Isostasy.GetIt,   	PETSC_NULL );

	// --- Output related parameters ---
	PetscOptionsGetInt(PETSC_NULL,"-Output.velocity" 		    ,	&user->Output.velocity				, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.temperature" 		,	&user->Output.temperature			, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.surface_topography" 	,	&user->Output.surface_topography	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.bottom_topography" 	,	&user->Output.surface_topography	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.quadrature" 			,	&user->Output.quadrature			, PETSC_NULL);

	// --- Pushing BC related parameters ---
	PetscOptionsGetInt(PETSC_NULL,"-AddPushing" 			 ,	&user->AddPushing				, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.num_changes"	 ,	&user->Pushing.num_changes		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.reset_pushing_coord",	&user->Pushing.reset_pushing_coord		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.L_block"      	 ,	&user->Pushing.L_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.W_block"      	 ,	&user->Pushing.W_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.H_block"      	 ,	&user->Pushing.H_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.x_center_block" ,	&user->Pushing.x_center_block 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.y_center_block" ,	&user->Pushing.y_center_block 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.z_center_block" ,	&user->Pushing.z_center_block 	, PETSC_NULL);


	char 			matprop_opt[PETSC_MAX_PATH_LEN];
	flg = PETSC_FALSE;

	for(i=0;i<user->Pushing.num_changes;i++){
		// V_push in cm/year
		sprintf(matprop_opt,"-Pushing.V_push_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.V_push[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    V_push[%lld]	= %g \n",(LLD)i,user->Pushing.V_push[i]);

		// Rate of rotation deg/yr
		sprintf(matprop_opt,"-Pushing.omega_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.omega[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Omega[%lld]	= %g \n",(LLD)i,user->Pushing.omega[i]);

		// Options whether to advect block for each time segment
		sprintf(matprop_opt,"-Pushing.coord_advect_%lld",(LLD)i);
		ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.coord_advect[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Coord_advect[%lld]	= %g \n",(LLD)i,user->Pushing.coord_advect[i]);

		// Options to define the direction of pushing
				sprintf(matprop_opt,"-Pushing.dir_%lld",(LLD)i);
				ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.dir[i]	, &flg); 				CHKERRQ(ierr);
				if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Direction[%lld]	= %g \n",(LLD)i,user->Pushing.dir[i]);
	}

	for(i=0;i<user->Pushing.num_changes+1;i++){
		//Time is a num_changes+1 array
		sprintf(matprop_opt,"-Pushing.time_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.time[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    time[%lld]	= %g \n",(LLD)i,user->Pushing.time[i]);
	}

	// linear solver options
	PetscOptionsGetInt( PETSC_NULL, "-StokesSolver", 				&user->StokesSolver, 			PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL, "-VelocitySolver", 				&user->VelocitySolver, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-VelocityTest", 				&user->VelocityTest, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-ScaleSystem", 				&user->ScaleSystem, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-use_fdstag_canonical", 	    &user->use_fdstag_canonical,    PETSC_NULL );

	// --- Get material properties from command line
	ierr = GetMaterialPropertiesFromCommandLine(user);

	/* In case we use FDSTAG, we need at least 2 elements (not 1)! */
	// If FDSTAG is used, dt must be different to 0 to avoid zero division at "ComputeStiffnessMatrixRHSTemperature_FDSTAG"
	if( __ELEMENT_TYPE__ == ELEMENT_FDSTAG ){
		if (user->nel_x==1){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the x-direction! I have increased this. \n");
			user->nel_x=2;
		}
		if (user->nel_y==1){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the y-direction! I have increased this. \n");
			user->nel_y=2;
		}
		if (user->dt==0.0){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG dt cannot be zero. Value has been changed to dt = 1e-3. \n");
			user->dt=1e-3;
		}
	}

	// Ensure that sticky air phase is not larger than the max. number of phases
	if (user->ErosionParameters.StickyAirPhase > (user->num_phases-1)){
		PetscPrintf(PETSC_COMM_WORLD," The sticky air phase is %i but the maximum phase in the model setup is %i. I changed sticky air phase to %i. \n",user->ErosionParameters.StickyAirPhase, (user->num_phases-1),(user->num_phases-1));
		user->ErosionParameters.StickyAirPhase = (user->num_phases-1);
	}


	user->Setup.Diapir_Hi 	 = (user->H-user->z_bot)*user->Hinterface;	// for 0-diapir setup
	user->Setup.SingleFold_H = 1.0;						    			// for  1-Single-layer folding setup


	// set this option to monitor actual option usage
	PetscOptionsInsertString("-options_left");

	if (!(user->InputParamFile))
	{

		// defined in NonDimensionalisation.h

		ComputeCharacteristicValues(user);

		/* Define phase Material properties in case NO input file is specified (=falling-block test)*/
		/* viscosity                                     								 density 									*/
		user->PhaseProperties.mu[0]  = 1.0*user->Characteristic.Viscosity; 		user->PhaseProperties.rho[0] = 1.0*user->Characteristic.Density;
		user->PhaseProperties.mu[1] = user->mumax; 								user->PhaseProperties.rho[1] = 2.0*user->Characteristic.Density;

		PetscOptionsInsertString("-AddRandomNoiseParticles 0");		// will always give the same results
	}



	{
		// Define material parameters for folding benchmarks, which overrule values from the parameters file
		PetscScalar	R, nL, nM;
		R 	= 	100;
		nL 	=	1;
		nM 	=	1;
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_R",  &R	, &found);	// Viscosity contrast
		if (found== PETSC_TRUE){
			user->PhaseProperties.mu[1] 			=	R;
			user->PhaseProperties.mu[0] 			=	1;
		}
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_nL", &nL	, &found);	// powerlaw exponent of layer
		if (found== PETSC_TRUE){
			user->PhaseProperties.n_exponent[1] 	=	nL;
		}
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_nM", &nM	, &found);	// powerlaw exponent of matrix
		if (found== PETSC_TRUE){
			user->PhaseProperties.n_exponent[0] 	=	nM;
		}

	}

	/* print information about simulation */
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Total # of cpu's               : %lld \n",(LLD)size);

	if(user->GravityAngle != 90.0)
	{
		PetscPrintf(PETSC_COMM_WORLD," Gravity angle with z-axis : %g \n", user->GravityAngle);
	}

	/* Info about particles if used */

	if (user->ParticleInput==1)
	{
		// Check whether a sufficient amount of particles are specified; if not increase it.
		// We need at least 2x2x2 particles per cell

		if(user->NumPartX < 2) { user->NumPartX = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in x-direction; Increasing\n"); }
		if(user->NumPartY < 2) { user->NumPartY = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in y-direction; Increasing\n"); }
		if(user->NumPartZ < 2) { user->NumPartZ = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in z-direction; Increasing\n"); }
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of tracers/cell         : [%lld,%lld,%lld] \n",(LLD)(user->NumPartX),(LLD)(user->NumPartY),(LLD)(user->NumPartZ));


	/* Compute some useful stuff */
	if (user->GridAdvectionMethod != 2  || user->BC.UpperBound != 0){
		user->num_subdt			  = 1;
	}
	PetscOptionsGetInt(PETSC_NULL ,"-num_subdt",	&user->num_subdt	, PETSC_NULL);  		//	# modify number of sub-dt

	
	if  (user->ApplyErosion==0){
		// non zero surface angle is only relevant if erosion is applied to the model
		user->SurfaceAngle 		   = 0;
		user->SurfaceNoiseAmplitude= 0;
	}

	/* Check for periodic boundary conditions */
	if ( (user->BC.LeftBound==3) || (user->BC.RightBound==3) ){
		user->BC.LeftBound=3;  user->BC.RightBound=3;			// if one is periodic, they both should be
		user->BC.BCType_x = DM_BOUNDARY_PERIODIC;
	}
	else{
		user->BC.BCType_x = DM_BOUNDARY_NONE;
	}

	if ( (user->BC.FrontBound==3) || (user->BC.BackBound==3) ){
		user->BC.FrontBound=3;  user->BC.BackBound=3;			// if one is periodic, they both should be
		user->BC.BCType_y = DM_BOUNDARY_PERIODIC;
	}
	else{
		user->BC.BCType_y = DM_BOUNDARY_NONE;
	}

	if ( (user->BC.UpperBound==3) || (user->BC.LowerBound==3) ){
		PetscPrintf(PETSC_COMM_WORLD," Periodic BCs are not yet implemented in these directions! \n");
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	else{
		user->BC.BCType_z = DM_BOUNDARY_NONE;
	}



	/* Set the density of 'sticky-air' to zero if we eliminate the 'sticky-air' from the system */
	if (user->ErosionParameters.UseInternalFreeSurface==1){
		PetscInt 	AirPhase;
		PetscBool	eliminate_stickyair_from_system;

		eliminate_stickyair_from_system = PETSC_FALSE;
		ierr 		= 	PetscOptionsGetBool( PETSC_NULL, "-eliminate_stickyair_from_system", &eliminate_stickyair_from_system,PETSC_NULL ); CHKERRQ(ierr);
		AirPhase 	= 	user->ErosionParameters.StickyAirPhase; 		// sticky air phase

		if ((eliminate_stickyair_from_system) && (user->PhaseProperties.rho[AirPhase]>0)){

			PetscPrintf(PETSC_COMM_WORLD," You eliminate the sticky air phase from the system, but the sticky air density is %f kg/m3 which is inconsistent. I have set air density to 0 kg/m3 !! \n", user->PhaseProperties.rho[AirPhase]);
			user->PhaseProperties.rho[AirPhase]	=	0.0;		// set to zero


		}

	}

	// show initial timestep (WARNING! which years if nondimensional?)
    if (user->DimensionalUnits==1){
        PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g years \n",user->dt);
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g  \n",user->dt);
    }
    // Boundary conditions
    PetscPrintf(PETSC_COMM_WORLD," BC employed                    : BC.[Left=%lld Right=%lld; Front=%lld Back=%lld; Lower=%lld Upper=%lld] \n",
                (LLD)(user->BC.LeftBound), (LLD)(user->BC.RightBound), (LLD)(user->BC.FrontBound), (LLD)(user->BC.BackBound), (LLD)(user->BC.LowerBound), (LLD)(user->BC.UpperBound) );

    
	/* Compute characteristic values */
	ComputeCharacteristicValues(user);

	/* Perform non-dimensionalisation of all input parameters */
	PerformNonDimensionalization(user);

	/* Allocate arrays required for time-dependent data */
	ierr = PetscMalloc( (size_t)user->time_end*sizeof(GlobalTimeDependentData), 	&user->TimeDependentData); CHKERRQ(ierr);
	ierr = PetscMemzero(user->TimeDependentData, (size_t)user->time_end*sizeof(GlobalTimeDependentData)); CHKERRQ(ierr);

	user->ErosionParameters.HorizontalFreeSurfaceHeight = user->ErosionParameters.InitialFreeSurfaceHeight;

	ierr = PetscOptionsInsertString( all_options ); CHKERRQ(ierr); /* force command line args in to override defaults */

	ierr = PetscFree(all_options); CHKERRQ(ierr);



	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGReadInputFile"
PetscErrorCode FDSTAGReadInputFile(JacRes *jr, UserContext *user)
{
	 // Parse the input file

	FILE *fp;
	PetscInt found;
	double d_values[1000], data;
	PetscInt i_values[1000];
	PetscInt nv, iphase,i;
	const PetscInt max_vals = 1000;
	char setup_name[PETSC_MAX_PATH_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fp = fopen( user->ParamFile, "r" );
	if(!fp)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s", user->ParamFile);
	}

	/* read the # of elements @ coarsest level ; this override any specification of #nodes */
	parse_GetInt( fp, "nel_x",   &user->nel_x, &found );
	parse_GetInt( fp, "nel_y",   &user->nel_y, &found );
	parse_GetInt( fp, "nel_z",   &user->nel_z, &found );

	// only number of cells is a relevant parameter for FDSTAG
	user->nnode_x = user->nel_x+1;
	user->nnode_y = user->nel_y+1;
	user->nnode_z = user->nel_z+1;


	/* Characteristic values, used to non-dimensionalize parameters ------------------------------------------------------------- */

	ierr = ScalingReadFromFile(&jr->scal, fp); CHKERRQ(ierr);

	parse_GetInt( fp,    "DimensionalUnits", &user->DimensionalUnits, &found );
	parse_GetDouble( fp, "Characteristic.Length", &data, &found );			if (user->DimensionalUnits==1){		user->Characteristic.Length       = data;	}// read data
	parse_GetDouble( fp, "Characteristic.Viscosity", &data, &found );		if (user->DimensionalUnits==1){		user->Characteristic.Viscosity    = data;	}// read data
	parse_GetDouble( fp, "Characteristic.Temperature", &data, &found );		if (user->DimensionalUnits==1){		user->Characteristic.Temperature  = data; 	}// read data
	parse_GetDouble( fp, "Characteristic.Stress", &data, &found );			if (user->DimensionalUnits==1){		user->Characteristic.Stress  	  = data;	}// read data
	/* ------------------------------------------------------------------------------------------------------------------------- */

	parse_GetDouble( fp, "L", &user->L, &found );
	parse_GetDouble( fp, "W", &user->W, &found );
	parse_GetDouble( fp, "H", &user->H, &found );
	parse_GetDouble( fp, "x_left", &user->x_left, &found );
	parse_GetDouble( fp, "y_front", &user->y_front, &found );
	parse_GetDouble( fp, "z_bot", &user->z_bot, &found );

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
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#ERROR! Incorrect model setup: %s", setup_name);
	}

	parse_GetInt( fp,    "Setup.Model", &user->Setup.Model, &found );
	parse_GetDouble( fp, "ampl2D", &user->ampl2D, &found );
	parse_GetDouble( fp, "ampl3D", &user->ampl3D, &found );
	parse_GetDouble( fp, "amplNoise", &user->amplNoise, &found );
	parse_GetDouble( fp, "Hinterface", &user->Hinterface, &found );

	parse_GetString( fp, "OutputFile", user->OutputFile, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetInt( fp,    "save_timesteps", &user->save_timesteps, &found );
	parse_GetInt( fp,    "time_end", &user->time_end, &found );
	parse_GetInt( fp,    "time_end_temp", &user->time_end_temp, &found );
	parse_GetDouble( fp, "CFL", &user->CFL, &found );
	parse_GetDouble( fp, "dt_max", &user->dt_max, &found );
	parse_GetDouble( fp, "dt_temp", &user->dt_temp, &found );

	parse_GetDouble( fp, "BC.Eyy", &user->BC.Eyy, &found );
	parse_GetDouble( fp, "BC.Exx", &user->BC.Exx, &found );

	parse_GetInt( fp, "internalBC_coord", &user->internalBC_coord, &found );
	parse_GetInt( fp, "internalBC_frontel", &user->internalBC_frontel, &found );
	parse_GetInt( fp, "internalBC_backel", &user->internalBC_backel, &found );
	parse_GetInt( fp, "internalBC_node", &user->internalBC_node, &found );
	parse_GetInt( fp, "zdepth_BC_el", &user->zdepth_BC_el, &found );
	parse_GetInt( fp, "zdepth_BC_node", &user->zdepth_BC_node, &found );
	parse_GetInt( fp, "BC.InternalBound", &user->BC.InternalBound, &found );
	parse_GetInt( fp, "BC.LeftBound", &user->BC.LeftBound, &found );
	parse_GetInt( fp, "BC.RightBound", &user->BC.RightBound, &found );
	parse_GetInt( fp, "BC.FrontBound", &user->BC.FrontBound, &found );
	parse_GetInt( fp, "BC.BackBound", &user->BC.BackBound, &found );
	parse_GetInt( fp, "BC.LowerBound", &user->BC.LowerBound, &found );
	parse_GetInt( fp, "BC.UpperBound", &user->BC.UpperBound, &found );
	parse_GetDouble( fp, "Temp_top", &user->Temp_top, &found );
	parse_GetDouble( fp, "Temp_bottom", &user->Temp_bottom, &found );
	parse_GetInt( fp, "temp_initialize", &user->temp_initialize, &found );

	parse_GetInt( fp, "UseInternalFreeSurface", 		&user->ErosionParameters.UseInternalFreeSurface, 		&found );
	parse_GetInt( fp, "StickyAirPhase", 				&user->ErosionParameters.StickyAirPhase, 				&found );
	parse_GetDouble( fp, "FSSA", 	 					&user->FSSA, 											&found );	// FSSA parameter
	parse_GetDouble( fp, "InitialFreeSurfaceHeight", 	&user->ErosionParameters.InitialFreeSurfaceHeight, 		&found );

	// -- parameters related to the erosion model employed
	parse_GetInt( fp, "ErosionModel", 				&user->ErosionParameters.ErosionModel, 		&found );	// which erosion model do we employ?	[0=default=none]

	// parameters in case we use the FD_Erosion model
	parse_GetInt   ( fp, "FE_ErosionCode.ResolutionFactorX", 		&user->ErosionParameters.FE_ErosionCode.ResolutionFactorX, 		&found );
	parse_GetInt   ( fp, "FE_ErosionCode.ResolutionFactorY", 		&user->ErosionParameters.FE_ErosionCode.ResolutionFactorY, 		&found );
	parse_GetDouble( fp, "FE_ErosionCode.InitialRandomNoise_m",  	&user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m, 	&found );
	parse_GetDouble( fp, "FE_ErosionCode.InitialUpliftedSide_m",	&user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m, 	&found );
	parse_GetDouble( fp, "FE_ErosionCode.rain_m_year", 				&user->ErosionParameters.FE_ErosionCode.rain_m_year, 			&found );
	parse_GetDouble( fp, "FE_ErosionCode.k0", 						&user->ErosionParameters.FE_ErosionCode.k0, 					&found ); // k = k0 + c*q^n [q=fluvial discharge,all the others are parameters specified here]
	parse_GetDouble( fp, "FE_ErosionCode.c", 						&user->ErosionParameters.FE_ErosionCode.c, 						&found );
	parse_GetDouble( fp, "FE_ErosionCode.n", 						&user->ErosionParameters.FE_ErosionCode.n, 						&found );
	parse_GetDouble( fp, "FE_ErosionCode.dt", 						&user->ErosionParameters.FE_ErosionCode.dt, 					&found ); // in years
	user->ErosionParameters.FE_ErosionCode.dt = user->ErosionParameters.FE_ErosionCode.dt*3600*24*365.25;	// in seconds
	parse_GetInt   ( fp, "FE_ErosionCode.fill_lake",                &user->ErosionParameters.FE_ErosionCode.fill_lake,              &found );
	parse_GetInt   ( fp, "FE_ErosionCode.BC",                       &user->ErosionParameters.FE_ErosionCode.BC,                     &found );
	parse_GetInt   ( fp, "FE_ErosionCode.mode_river",               &user->ErosionParameters.FE_ErosionCode.mode_river,             &found );
	parse_GetInt   ( fp, "FE_ErosionCode.nbre_river",               &user->ErosionParameters.FE_ErosionCode.nbre_river,             &found );
	parse_GetDouble( fp, "FE_ErosionCode.rain_river_year", 			&user->ErosionParameters.FE_ErosionCode.rain_river_year, 		&found );


	parse_GetInt( fp, "SedimentationModel", 				&user->ErosionParameters.SedimentationModel, 		&found );	// which sedimentation model do we employ?
	parse_GetDouble( fp, "SedimentationRate_cmYr", 	 	&user->ErosionParameters.SedimentationRate_cmYr, 		&found );
	parse_GetInt( fp, 		"PhaseFirstSedimentedLayer", 		&user->ErosionParameters.PhaseFirstSedimentedLayer, 	&found );
	parse_GetInt( fp, 		"PhaseLastSedimentedLayer", 		&user->ErosionParameters.PhaseLastSedimentedLayer, 		&found );
	parse_GetDouble( fp, 	"SedimentLayerThicknessYears", 		&user->ErosionParameters.SedimentLayerThicknessYears, 	&found );

	parse_GetInt( fp,    "GridAdvectionMethod", &user->GridAdvectionMethod, &found );
	parse_GetDouble( fp, "FactorSurfaceLayer", &user->FactorSurfaceLayer, &found );
	parse_GetInt( fp,    "num_subdt", &user->num_subdt, &found );

	// linear solver options
	parse_GetInt( fp,    "StokesSolver", &user->StokesSolver, &found );
	parse_GetInt( fp,    "VelocitySolver", &user->VelocitySolver, &found );


	/* Particle related variables */
	parse_GetInt( fp,    "ParticleInput", &user->ParticleInput, &found );
	parse_GetInt( fp,    "LoadInitialParticlesFromDisc", &user->LoadInitialParticlesFromDisc, &found );
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

	parse_GetInt( fp,    "InitialMeshFromFile", &user->InitialMeshFromFile, &found );
	parse_GetString( fp, "InitialMeshFileName", user->InitialMeshFileName, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetInt( fp,    "InitialMantleLevel", &user->InitialMantleLevel, &found );

	parse_GetInt( fp,    "InitialErosionSurfaceFromFile", &user->InitialErosionSurfaceFromFile, &found );

	/* Read material properties - This manner of setting material properties will be disabled in the near future, and
	 * replaced with a more general routine
	 */


	parse_GetDouble( fp, "LowerViscosityCutoff", &user->LowerViscosityCutoff, &found );
	parse_GetDouble( fp, "UpperViscosityCutoff", &user->UpperViscosityCutoff, &found );

	parse_GetDouble( fp, "Gravity", &user->Gravity, &found );
	parse_GetInt( fp,    "PlasticityModel", &user->PlasticityModel, &found );
	parse_GetDouble( fp, "Xi", &user->Xi, &found );
	parse_GetDouble( fp, "GasConstant", &user->GasConstant, &found );

	// --- SurfaceVelocity related input parameters ---
	parse_GetInt( fp,    "get_SurfVelField", 	&user->SurfVelField.GetIt, 		&found );
	parse_GetInt( fp,    "SurfVelField_SaveRef", &user->SurfVelField.SaveRef,	&found );
	parse_GetString( fp, "SurfVelField_RefDatFile", user->SurfVelField.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetDouble( fp, "SurfVelField_VxStdDev", &user->SurfVelField.VxStdDev, &found );
	parse_GetDouble( fp, "SurfVelField_VyStdDev", &user->SurfVelField.VyStdDev, &found );
	parse_GetDouble( fp, "SurfVelField_VzStdDev", &user->SurfVelField.VzStdDev, &found );
	// --- Optimization related input parameters ---
	parse_GetInt( fp,    "get_Misfit", 				&user->Optimisation.GetIt, 		&found );

	// --- Gravity related input parameters ---
	parse_GetInt( fp,    "get_GravityField", 		&user->GravityField.GetIt, 		&found );
	parse_GetInt( fp,    "GravityField_SaveRef", 	&user->GravityField.SaveRef, 	&found );
	parse_GetInt( fp,    "GravityField_SaveDebug", 	&user->GravityField.SaveDebug, 	&found );
	parse_GetInt( fp,    "GravityField_SaveVTK", 	&user->GravityField.SaveVTK, 	&found );
	parse_GetInt( fp,    "GravityField_survey_nx", &user->GravityField.survey_nx	, &found );
	parse_GetInt( fp,    "GravityField_survey_ny", &user->GravityField.survey_ny	, &found );
	parse_GetDouble( fp, "GravityField_survey_xs", &user->GravityField.survey_xs	, &found );
	parse_GetDouble( fp, "GravityField_survey_xm", &user->GravityField.survey_xm	, &found );
	parse_GetDouble( fp, "GravityField_survey_ys", &user->GravityField.survey_ys	, &found );
	parse_GetDouble( fp, "GravityField_survey_ym", &user->GravityField.survey_ym	, &found );
	parse_GetDouble( fp, "GravityField_survey_z"	, &user->GravityField.survey_z 	, &found );
	parse_GetDouble( fp, "GravityField_ReferenceDensity", &user->GravityField.ReferenceDensity, &found );
	parse_GetInt( fp,    "GravityField_num_intp"	, &user->GravityField.num_intp		, &found );
	parse_GetString( fp, "GravityField_RefDatFile", user->GravityField.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetDouble( fp, "GravityField_StdDev", &user->GravityField.StdDev, &found );
	parse_GetInt( fp,    "GravityField_LithColNum", &user->GravityField.LithColNum, &found );
	parse_GetDoubleArray(fp,"GravityField_LithColDepth",&nv,d_values, &found );

	if (found!=0)
	{
		for( i=0; i<user->GravityField.LithColNum; i++ )
        {
            user->GravityField.LithColDepth[i] = d_values[i] ;
            PetscPrintf(PETSC_COMM_WORLD,"# LithColDepth[%lld] = %g \n",(LLD) i,user->GravityField.LithColDepth[i]);
		}
	}

	parse_GetDoubleArray( fp,    "GravityField_LithColDens",&nv,d_values, &found );
	if (found!=0){
        for( i=0; i<user->GravityField.LithColNum+1; i++ )
        {
            user->GravityField.LithColDens[i] = d_values[i] ;
            PetscPrintf(PETSC_COMM_WORLD,"# LithColDens[%lld] = %g \n",(LLD) i,user->GravityField.LithColDens[i]);
		}
	}

	// --- Isostasy related input parameters ---
	parse_GetInt( fp,    "get_AiryIsostasy", &user->Isostasy.GetIt, &found );
	parse_GetInt( fp,    "Isostasy_SaveRef", &user->Isostasy.SaveRef, &found );
	parse_GetInt( fp, "Isostasy_ref_xi", &user->Isostasy.ref_xi, &found );
	parse_GetInt( fp, "Isostasy_ref_yi", &user->Isostasy.ref_yi, &found );
	parse_GetDouble( fp, "Isostasy_corr_topo", &user->Isostasy.corr_topo,  &found );
	parse_GetDouble( fp, "Isostasy_TisoStdDev", &user->GravityField.StdDev, &found );
	parse_GetDouble( fp, "Isostasy_RefRho", &user->Isostasy.ref_rho, &found );
	parse_GetString( fp, "Isostasy_RefDatFile", user->Isostasy.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );


	parse_GetInt( fp,    "save_breakpoints", &user->save_breakpoints, &found );

	// --- Output related input parameters ---
	parse_GetInt( fp,    "Output.velocity", &user->Output.velocity, &found );
	parse_GetInt( fp,    "Output.temperature", &user->Output.temperature, &found );
	parse_GetInt( fp,    "Output.surface_topography", &user->Output.surface_topography, &found );
	parse_GetInt( fp,    "Output.bottom_topography", &user->Output.bottom_topography, &found );
	parse_GetInt( fp,    "Output.quadrature", &user->Output.quadrature, &found );

	// --- Pushing related input parameters ---
	parse_GetInt( fp,    "AddPushing", &user->AddPushing, &found );
	parse_GetInt( fp,    "Pushing.num_changes", &user->Pushing.num_changes, &found );
	parse_GetInt( fp,    "Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, &found );
	parse_GetDoubleArray( fp,    "Pushing.time",  &nv, d_values, &found );
	for( i=0; i<user->Pushing.num_changes+1; i++ ) {
		user->Pushing.time[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# time stage = %g \n",user->Pushing.time[i]);
	}
	parse_GetDoubleArray( fp, "Pushing.V_push", &nv, d_values, &found );
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.V_push[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# V_push stage = %g \n",user->Pushing.V_push[i]);
	}
	parse_GetDoubleArray( fp, "Pushing.omega",&nv, d_values, &found  );
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.omega[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# omega = %g \n",user->Pushing.omega[i]);
	}
	parse_GetIntArray( fp, "Pushing.coord_advect",&nv, i_values, &found);
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.coord_advect[i] = i_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# coord_advect = %d \n",user->Pushing.coord_advect[i]);
	}
	parse_GetIntArray( fp, "Pushing.dir",&nv, i_values, &found);
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.dir[i] = i_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# direction = %d \n",user->Pushing.dir[i]);
	}
	parse_GetDouble( fp, "Pushing.L_block", &user->Pushing.L_block, &found );
	parse_GetDouble( fp, "Pushing.W_block", &user->Pushing.W_block, &found );
	parse_GetDouble( fp, "Pushing.H_block", &user->Pushing.H_block, &found );
	parse_GetDouble( fp, "Pushing.x_center_block", &user->Pushing.x_center_block, &found );
	parse_GetDouble( fp, "Pushing.y_center_block", &user->Pushing.y_center_block, &found );
	parse_GetDouble( fp, "Pushing.z_center_block", &user->Pushing.z_center_block, &found );

	/* -------------------------------------------------------------------------------------------------------------------------
	 * Read phase transitions related information -
	 * The input structure of this is also likely to change but at a later stage
	 */
	parse_GetInt( fp,    "num_phase_transitions", &user->num_phase_transitions, &found );

	parse_GetIntAllInstances( fp,    "TransitionType", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionType = i_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "TransitionBelow", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionBelow = i_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionDepth", &nv, d_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionDepth = d_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionP0", &nv, d_values, max_vals, &found );		//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionP0 = d_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionAlpha", &nv, d_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionAlpha = d_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "InitialPhase", &nv, i_values, max_vals, &found );		//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].InitialPhase = i_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "TransformedPhase", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransformedPhase = i_values[iphase] ;
	}
	/* ------------------------------------------------------------------------------------------------------------------------- */

	fclose(fp);

	ierr = ReadMaterialProperties(user); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ReadMaterialProperties"
PetscErrorCode ReadMaterialProperties(UserContext *user)
{
	PetscInt      found_data;
	PetscInt      n_phases,n_attrs, n_types;
	Phase         *phases, p;
	Attribute     *attrs, a;
	AttributeType *types, t;
	PetscInt      PP,AA,TT;
	PetscInt      attr_match, type_match;
	PetscInt      p_idx, a_idx, t_idx;

	PetscFunctionBegin;

	//-----------------------------------------------------------------------
	// Read material parameters from the input file (new, more general, input format for material properties)
	//
	// At this stage, the new data is read from the input file and added to the phase properties array.
	//
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
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many phases in the input file! Actual: %lld, Maximum: %lld", n_phases, max_num_phases);
	}

	// store total number of phases
	user->num_phases = n_phases;

	for( PP=0; PP<n_phases; PP++ )
	{
		p = phases[PP];

		// check phase numbering
		if(p->phase_number != PP)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Phase numbering must be sequential in the input file! ID: %lld, Index: %lld", p->phase_number, PP);
		}

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

