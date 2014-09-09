/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

LaMEMLib.c, contains the following subroutine:
LaMEMLib - Main routine of LaMEM library

$Id$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"

#include "Version.h"

#include "Parsing.h"
#include "Quadrature.h"
#include "Solvers.h"
#include "Breakpoint.h"
#include "Output.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "Mesh.h"
#include "Assembly.h"
#include "Utils.h"
#include "LaMEM_Initialize.h"
#include "LaMEM_Temperature.h"
#include "LaMEM_Particles.h"
#include "LaMEM_AnalyticalSolutions.h"
#include "Assembly_FDSTAG.h"
#include "ApplyBoundaryConditions.h"
#include "GetGravityField.h"
#include "GetSurfaceVelocity.h"
#include "LaMEM_FE_ErosionCode.h"
#include "Addition_FDSTAG.h"

PetscInt __ELEMENT_TYPE__;
PetscInt __Q2_TYPE__;

//==========================================================================================================
// LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

#undef __FUNCT__
#define __FUNCT__ "LaMEMLib"
PetscErrorCode LaMEMLib(PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id)
{

	PetscInt found_data;
	PetscBool InputParamFile, ElemTypeSet, use_fdstag_canonical;
	char ParamFile[PETSC_MAX_PATH_LEN], ElemType[PETSC_MAX_PATH_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether input file is specified
	ierr = PetscOptionsGetString(PETSC_NULL, "-ParamFile", ParamFile, PETSC_MAX_PATH_LEN, &InputParamFile); CHKERRQ(ierr);

	// read additional PETSc options from input file
	if(InputParamFile == PETSC_TRUE)
	{
		ierr = PetscOptionsReadFromFile(ParamFile, &found_data, 1); CHKERRQ(ierr);
	}

	// check whether FDSTAG canonical must be used
	ierr = PetscOptionsHasName(PETSC_NULL, "-use_fdstag_canonical", &use_fdstag_canonical); CHKERRQ(ierr);

	if(use_fdstag_canonical == PETSC_TRUE)
	{
		// check whether basic element type is FDSTAG
		ierr = PetscOptionsGetString(PETSC_NULL, "-vpt_element", ElemType, PETSC_MAX_PATH_LEN, &ElemTypeSet);
		if(ElemTypeSet != PETSC_TRUE || strcmp(ElemType, "FDSTAG"))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Option -use_fdstag_canonical currently requires setting -vpt_element FDSTAG\n");
		}

		// using particles is required
#ifndef PARTICLES
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! recompile with -DPARTICLES to activate option -use_fdstag_canonical\n");
#endif

		// call FDSTAG solution routine
		ierr = LaMEMLib_FDSTAG(InputParamFile, ParamFile, LaMEM_OutputParameters, mpi_group_id); CHKERRQ(ierr);
	}
	else
	{
		// just call legacy solution routine for all other cases
		ierr = LaMEMLib_Legacy(InputParamFile, ParamFile, LaMEM_OutputParameters, mpi_group_id); CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "LaMEMLib_Legacy"
PetscErrorCode LaMEMLib_Legacy(PetscBool InputParamFile, const char *ParamFile, PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id)
{
	UserContext        user;
	PetscMPIInt        rank;
	PetscInt           mx, my, mz, SaveOrNot, i;
	PetscInt           itime, it_step;
	PetscBool          UseOld_InitializeParticles;
	PetscScalar        MaxVel, MinVel, dx, dy, dz, dmin, fac_time, velocity_error;
	PetscLogDouble     cputime_start,cputime_start0, cputime_end, cputime_start_tstep, cputime_start_nonlinear=0.0;
	Vec                rhs, rhs_local, sol, sol_local, rhs_add_BC, rhs_add_BC_loc, sol_old, sol_advect;
	Vec                Pressure, Pressure_local, Temp;
	Vec                rhs_p, rhs_p_local, rhs_add_pBC, rhs_add_pBC_loc,pv_rhs_push;
	Vec                ViscosityScaling;
#ifdef TEMPERATURE
	KSP                ksp_temp=PETSC_NULL;
	Vec                Temp_local, rhs_Temp, rhs_Temp_local, Temp_old;
#endif
	DAVPElementType    vpt_element_type;
	LaMEMVelPressureDA C;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscTime(&cputime_start0);

	// Start code
	PetscPrintf(PETSC_COMM_WORLD,"# -------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"#                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"#              Current Revision: %s - %s		     \n",__SVNREVISION__,__SVNDATE__);
	PetscPrintf(PETSC_COMM_WORLD,"#  Modified items: %s\n",__SVNMANCHANGES__);
	PetscPrintf(PETSC_COMM_WORLD,"#     Compiled: Date: %s - Time: %s - Mode:  %s	    \n",__DATE__,__TIME__,__OPTMODE__);
	PetscPrintf(PETSC_COMM_WORLD,"# -------------------------------------------------------------------------- \n");

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

	// Initialize context
	PetscMemzero( &user, sizeof(UserContext) );

	// set input file flag 
	user.InputParamFile = InputParamFile;

	if(InputParamFile == PETSC_TRUE)
	{
		strcpy(user.ParamFile, ParamFile);
	}

	ierr = LaMEMVelPressureDACreate(DAVP_FDSTAG, &C); CHKERRQ(ierr);

	// get the element type
	ierr = LaMEMVelPressureDAGetInfo( C, &vpt_element_type,0,0, 0,0,0,0, 0,0,0,0, 0,0 ); CHKERRQ(ierr);

	// define the yucky global variable
	if( vpt_element_type == DAVP_Q1P0 )
	{
		__ELEMENT_TYPE__ = ELEMENT_Q1P0;
	}
	if( vpt_element_type == DAVP_Q2PM1L )
	{
		__ELEMENT_TYPE__ = ELEMENT_Q2P1;
		__Q2_TYPE__      = ELEMENT_Q2P1_LOCAL;
	}
	if( vpt_element_type == DAVP_Q2PM1G )
	{
		__ELEMENT_TYPE__ = ELEMENT_Q2P1;
		__Q2_TYPE__      = ELEMENT_Q2P1_GLOBAL;
	}
	if( vpt_element_type == DAVP_Q1Q1 )
	{
		__ELEMENT_TYPE__ = ELEMENT_Q1Q1;
	}
	if( vpt_element_type == DAVP_FDSTAG )
	{
		__ELEMENT_TYPE__ = ELEMENT_FDSTAG;
	}

	// initialize variables
	ierr = 	InitializeCode(&user); CHKERRQ(ierr);

	// Give current LaMEM session a specific group ID
	user.Optimisation.mpi_group_id = *mpi_group_id;

	// Check LaMEM input variables, and stop code if impossible combinations are added
	ierr = LaMEM_Initialize_CheckInput(); CHKERRQ(ierr);

	//========================================================================================
	// Setting up the solver
	//========================================================================================

	PetscTime(&cputime_start);

	// Initialize DMDA and matrices, that:
	// (1) Store material properties,
	// (2) Distribute the structured grid on the various processors
	// (3) Are required for the velocity and pressure Stokes solution

	LaMEM_Initialize_StokesSolver(&user, vpt_element_type, C);

	// Initialize DMDA and matrices for temperature
#ifdef TEMPERATURE
	LaMEM_Initialize_TemperatureSolver(&user, vpt_element_type, C);
#endif

	// Initialize erosion solver and surface if required
	ierr = InitializeInternalErosionSurfaceOnRankZero(&user); CHKERRQ(ierr);

	// Call some initial functions related to the erosion code
	if(user.ErosionParameters.ErosionModel == 2)
	{
		if (user.InitialErosionSurfaceFromFile == 1)
		{
			ierr = LoadInitialErosionSurface(&user); CHKERRQ(ierr);
		}
		ierr = SaveInitialErosionSurface(&user,"InitialErosionSurface"); CHKERRQ(ierr);
	}

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"# Finished initializing Stokes matrices : %g s \n",cputime_end - cputime_start);
	PetscTime(&cputime_start);

	// Generate and initialize global Vectors for Stokes and Temperature
	ierr = DMCreateGlobalVector(user.DA_Vel, &rhs);              CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Vel, &rhs_local);        CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.DA_Vel, &rhs_add_BC);       CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Vel, &rhs_add_BC_loc);   CHKERRQ(ierr);
	ierr = VecSet(rhs,            0.0);                          CHKERRQ(ierr);
	ierr = VecSet(rhs_local,      0.0);                          CHKERRQ(ierr);
	ierr = VecSet(rhs_add_BC,     0.0);                          CHKERRQ(ierr);
	ierr = VecSet(rhs_add_BC_loc, 0.0);                          CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(user.DA_Vel, &sol);              CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Vel, &sol_local);        CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.DA_Vel, &sol_old);          CHKERRQ(ierr);
	ierr = VecSet(sol_old,   0.0);                               CHKERRQ(ierr);
	ierr = VecSet(sol_local, 0.0);                               CHKERRQ(ierr);
	ierr = VecSet(sol,       0.0);                               CHKERRQ(ierr);
	// later used to advect properties
	ierr = VecDuplicate(sol,  &sol_advect);                      CHKERRQ(ierr);
	ierr = VecSet(sol_advect, 0.0);                              CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(user.DA_Pres, &Pressure);        CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Pres, &Pressure_local);  CHKERRQ(ierr);
	ierr = VecSet(Pressure,       0.0);                          CHKERRQ(ierr);
	ierr = VecSet(Pressure_local, 0.0);                          CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(user.DA_Pres, &rhs_p);           CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Pres, &rhs_p_local);     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.DA_Pres, &rhs_add_pBC);     CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Pres, &rhs_add_pBC_loc); CHKERRQ(ierr);
	ierr = VecSet(rhs_p,           0.0);                         CHKERRQ(ierr);
	ierr = VecSet(rhs_p_local,     0.0);                         CHKERRQ(ierr);
	ierr = VecSet(rhs_add_pBC,     0.0);                         CHKERRQ(ierr);
	ierr = VecSet(rhs_add_pBC_loc, 0.0);                         CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = DMCreateGlobalVector(user.DA_Temp, &Temp);            CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.DA_Temp, &Temp_old);        CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Temp, &Temp_local);      CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user.DA_Temp, &rhs_Temp);        CHKERRQ(ierr);
	ierr = DMCreateLocalVector (user.DA_Temp, &rhs_Temp_local);  CHKERRQ(ierr);

	ierr = VecSet(Temp_local,     0.0);                          CHKERRQ(ierr);
	ierr = VecSet(Temp_old,       0.0);                          CHKERRQ(ierr);
	ierr = VecSet(rhs_Temp,       0.0);                          CHKERRQ(ierr);
	ierr = VecSet(Temp,           0.0);                          CHKERRQ(ierr);
	ierr = VecSet(rhs_Temp_local, 0.0);                          CHKERRQ(ierr);
#endif

	// Create a vector that holds the inverse of viscosity (for scaling)
	ierr = VecDuplicate(Pressure, &ViscosityScaling); CHKERRQ(ierr);

	// Generate fine-grid mesh
	ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);

	// Save Partitioning
	if(user.SavePartitioning)
	{
		SaveProcessorPartitioning(&user);
	}

	if (user.BC.InternalBound > 0)
	{
		DefineInternalBC(&user);
	}

	// Create a deformed mesh
	ierr = DMDASetUniformCoordinates(user.DA_Vel,user.x_left,user.x_left+user.W,user.y_front,user.y_front + user.L,user.z_bot,user.z_bot + user.H); CHKERRQ(ierr);
	ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = GenerateMeshFine(user.DA_Temp, &user); CHKERRQ(ierr);
#endif

	// Set the initial material properties
	ierr = SetMaterialProperties( C, &user, user.DA_Vel); CHKERRQ(ierr);

	// INITIALIZE PARTICLES-based routines
#ifdef PARTICLES

	if ((user.LoadInitialParticlesFromDisc!=0) && (vpt_element_type == DAVP_FDSTAG) )
	{
		// [0]. Load particles from InitialParticles directory for FDSTAG if LoadInitialParticlesFromDisc==1
		// [1]. Load particles from Matlab generated files if LoadInitialParticlesFromDisc==2

		ierr = LoadInitialParticlesFromDisc_FDSTAG(&user); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"# Successful loading of particles from disc for FDSTAG \n");

		if (user.LoadInitialParticlesFromDisc==2){
			// Find the element in which particles reside
			ierr = GetParticleNodes(C, user.DA_Vel, &user); CHKERRQ(ierr);		// Put the particle to the correct nodes
		}

	}
	else
	{
		// [2]. Particles are initialized internally - works with ALL OTHER cases for Model.Setup
		UseOld_InitializeParticles = PETSC_FALSE;
		ierr = PetscOptionsHasName(PETSC_NULL,"-UseOld_InitializeParticles",&UseOld_InitializeParticles); CHKERRQ(ierr);

		// For backwards compatibility, we here allow the use of the old InitializeParticles routine, by adding
		// -UseOld_InitializeParticles
		// to the command-line. This will be removed from LaMEM at some point in the future

		if (UseOld_InitializeParticles==PETSC_TRUE)
		{
			PetscPrintf(PETSC_COMM_WORLD,"# NOTE: using InitializeParticles instead of the (default) InitializeParticles_ElementWise to initialize Particles \n");
			ierr = InitializeParticles(user.DA_Vel, &user); 					CHKERRQ(ierr);		// Initialize tracers
		}
		else
		{
			ierr = InitializeParticles_ElementWise(C, user.DA_Vel, &user); 	CHKERRQ(ierr);			// Initialize tracers
		}


		// Find the element in which particles reside
		ierr = GetParticleNodes(C, user.DA_Vel, &user);	CHKERRQ(ierr); // Put the particle to the correct elements

		// Set tracer properties from
		ierr = SetInitialTracerProperties(user.DA_Vel, &user); CHKERRQ(ierr); // Set initial tracer phases

		// Read the initial mesh from file if asked for
		if (user.InitialMeshFromFile==1)
		{
			ierr = ReadMeshFromFile(user.DA_Vel, &user); CHKERRQ(ierr);
			ierr = InterpolateParticles(C,user.DA_Vel, &user); 				CHKERRQ(ierr);
		}

		if (user.Setup.Model==3 )
		{
			// This is the 'old' way of setting particle properties from a phase-grid that is regularly spaced
			PetscPrintf(PETSC_COMM_WORLD," You are using the OLD way of setting tracer properties. This piece of code will be REMOVED in a future version of LaMEM - Please change your input scripts!! \n");
			ierr = SetInitialTracerPhasesFromFile(user.DA_Vel, &user); 		CHKERRQ(ierr); 		// Read the tracer phase composition & temperature from file
		}

		if (user.LoadInitialParticlesFromDisc==1 )
		{
			// The 'new' way of setting particles
			ierr = LoadInitialParticlesFromDisc( &user); CHKERRQ(ierr);

			PetscPrintf(PETSC_COMM_WORLD,"# Successful loading of particles from disc \n");

		}

	}

#endif

#ifndef PARTICLES
	// Read the initial mesh from file if asked for
	if (user.InitialMeshFromFile==1)
	{
		ierr = ReadMeshFromFile(user.DA_Vel, &user); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"# Successful read initial mesh from disc \n");

	}
#endif

	// Read particles and properties from a breakpoint file if this was requested and if this is possible
	if (user.restart==1)
	{
		ierr = LoadBreakPoint(C, &user,user.DA_Vel, sol, Temp, Pressure, 0); CHKERRQ(ierr);
	}

	// Save initial mesh to file
	if (user.save_breakpoints > 0)
	{
		PetscPrintf(PETSC_COMM_WORLD," Saving initial mesh to file \n");
		ierr = SaveInitialMesh(&user,user.DA_Vel,"InitialMesh"); CHKERRQ(ierr);
	}


	// If we are using FDSTAG, the grid MUST be regular and undeformed.
	// Because of the Finite Element manner in which we create the mesh and particles,
	// we can actually initialize the mesh and particles in an irregular manner (or even read an arbitrary mesh from file) and set the initial particle
	// distribution quite accurately in this manner. Yet, at this stage we must ensure that the grid is really undeformed (otherwise we'll have
	// convergence issues with the code that are not easy to detect).
	if (vpt_element_type == DAVP_FDSTAG)
	{


		user.ampl2D    = 0;
		user.ampl3D    = 0;
		user.amplNoise = 0;
		ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);
#ifdef TEMPERATURE
		ierr = GenerateMeshFine(user.DA_Temp, &user); CHKERRQ(ierr);
#endif

	}

	if ((user.save_breakpoints >= 0) && (user.LoadInitialParticlesFromDisc==0))
	{
		// Save the initial particles to disk (every proc saves one file)
		// Deactivate this by adding -save_breakpoints -1  to the command line (e.g. while performing scaling tests)
		ierr = WriteInitialParticlesToDisc( &user ); CHKERRQ(ierr);	// Write the initial particles to the ./InitialParticles directory
	}


	// Interpolate temperature from particles 2 nodes (before starting the time step loop and in case we use the FEM) */
#ifdef TEMPERATURE
	if ((user.time_start==0) && (vpt_element_type != DAVP_FDSTAG))
	{
		ierr = ParticleTemperatureToNodes(C, &user, user.DA_Temp, Temp_local, Temp); CHKERRQ(ierr);
	}
	ierr = IntPointProperties( C, user.DA_Vel, PETSC_NULL, PETSC_NULL, PETSC_NULL, user.DA_Temp, Temp, &user, user.dt, 0 ); CHKERRQ(ierr);
#endif

	// ========================================================================================================================
	// Beginning of time step
	if (user.time_start == 0) { user.PlasticityCutoff = 1; }
	else                      { user.PlasticityCutoff = 0; }

	if(user.fileno > 0.0)
	{
		user.dt = user.dt_max;
	}

	//===============
	// TIME STEP LOOP
	//===============
	for (itime = user.time_start; itime < user.time_end; itime++)
	{
		PetscTime(&cputime_start_tstep);

		// save current iteration step
		user.itime = itime;

		//==========================================================================================
		if (user.EulerianAfterTimestep>0)
		{
			// only activated if  -EulerianAfterTimestep ??   is present on the command line
			if (itime>user.EulerianAfterTimestep)
			{
				user.GridAdvectionMethod  = 0; // Eulerian
				PetscPrintf(PETSC_COMM_WORLD," Activating EULERIAN mode \n");
			}
		}

		//==========================================================================================
		// Interpolate meshes material properties from tracers to grid
		PetscTime(&cputime_start);

		// Create new material parameters
#ifdef PARTICLES
		if (user.ParticleInput==1 )
		{
			if (vpt_element_type == DAVP_FDSTAG)
			{
				PetscBool  OutputVTK_ViscosityDensity;

				// Correct particles in case we employ an internal free surface
				ierr = CorrectPhasesForInternalFreeSurface( &user);			CHKERRQ(ierr);

				ierr = ComputeHistoryPropertiesFromParticles_FDSTAG(C, &user); CHKERRQ(ierr);	// compute history-dependent material properties from particles (using distance-based averaging)

				// Update material properties on FDSTAG grid
				ierr = UpdateMaterialProperties_FDSTAG(&user); CHKERRQ(ierr);

				OutputVTK_ViscosityDensity = PETSC_FALSE;
				ierr = PetscOptionsHasName(PETSC_NULL,"-OutputVTK_ViscosityDensity",&OutputVTK_ViscosityDensity); CHKERRQ(ierr);

				if (OutputVTK_ViscosityDensity)
				{
					OutputVTK_ViscosityDensityFields_FDSTAG( &user, Temp );		// write data to disk
				}
			}
			else
			{
				ierr = MaterialPropertiesFromTracers(C,&user); CHKERRQ(ierr);		// compute material properties from particles
			}
		}
#endif

		// Interpolate temperature from particles to nodes
#ifdef TEMPERATURE
		if (vpt_element_type != DAVP_FDSTAG)
		{	// only for FEM: FDSTAG does this in ComputePropertiesFromParticles_FDSTAG
			ierr = ParticleTemperatureToNodes(C, &user, user.DA_Temp, Temp_local, Temp_old); CHKERRQ(ierr);
		}
#endif
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD,"#  Interpolating meshes and temperatures took %g s \n",cputime_end-cputime_start);

		//==========================================================================================
		// Set material properties at integration points in case we are performing a benchmark
		if (user.AnalyticalBenchmark == PETSC_TRUE)
		{
			PetscBool  OutputVTK_ViscosityDensity;

			ierr = LaMEM_Initialize_AnalyticalSolution(&user, C); CHKERRQ(ierr);

			OutputVTK_ViscosityDensity = PETSC_FALSE;
			ierr = PetscOptionsHasName(PETSC_NULL,"-OutputVTK_ViscosityDensity",&OutputVTK_ViscosityDensity); CHKERRQ(ierr);

			if (OutputVTK_ViscosityDensity)
			{
				OutputVTK_ViscosityDensityFields_FDSTAG( &user, Temp );		// write data to disk
			}

		}

		//=========================================================================================
		//	NONLINEAR THERMO-MECHANICAL SOLVER
		//=========================================================================================
		if(!user.SkipStokesSolver)
		{
			ierr = ResetStressesBeforeNonlinearIterations(C, &user); CHKERRQ(ierr);	// reset stresses (required for plasticity)
			ierr = VecSet(sol_old, 0.0);                             CHKERRQ(ierr);

			if (user.temp_initialize==1 && itime > user.time_end_temp)
			{
				user.temp_initialize=0;
			}

			it_step 				= 	0;
			velocity_error 			= 	10000.0;
			user.PlasticityCutoff 	= 	1;

			//=====================
			// NONLINEAR ITERATIONS
			//=====================
			while (velocity_error>user.NonlinearIterationsAccuracy &&  user.temp_initialize==0)
			{
				PetscTime(&cputime_start_nonlinear);

				PetscTime(&cputime_start);

				// zero out solution and rhs vectors
				ierr = VecSet(rhs,             0.0); CHKERRQ(ierr);
				ierr = VecSet(rhs_add_BC,      0.0); CHKERRQ(ierr);
				ierr = VecSet(rhs_p,           0.0); CHKERRQ(ierr);
				ierr = VecSet(rhs_add_pBC,     0.0); CHKERRQ(ierr);
				ierr = VecSet(rhs_add_BC_loc,  0.0); CHKERRQ(ierr);
				ierr = VecSet(rhs_add_pBC_loc, 0.0); CHKERRQ(ierr);

				// Create vector for pushingBC - PV
				ierr = VecDuplicate(rhs_p,     &pv_rhs_push); CHKERRQ(ierr);
				ierr = VecSet(pv_rhs_push,     0.0);          CHKERRQ(ierr);

				//==========================================================================================
				// Compute stiffness matrix

				// Zero all matrices
				ierr = MatZeroEntries(user.VV_MAT);   CHKERRQ(ierr);
				ierr = MatZeroEntries(user.PP_MAT);   CHKERRQ(ierr);
				ierr = MatZeroEntries(user.PV_MAT);   CHKERRQ(ierr);
				ierr = MatZeroEntries(user.VP_MAT);   CHKERRQ(ierr);
				ierr = MatZeroEntries(user.approx_S); CHKERRQ(ierr);

				// Compute stiffness matrix & additions to RHS because of symmetrizing BC's
				ierr = ComputeStiffnessMatrixes( C, user.DA_Vel, user.DA_Pres,
						user.VV_MAT, user.PP_MAT, user.PV_MAT,
						user.VP_MAT, user, user.dt, rhs_add_BC, rhs_add_BC_loc, rhs_add_pBC,pv_rhs_push,
						ViscosityScaling, user.approx_S, &user.remesh ); CHKERRQ(ierr);

				ierr = VecAssemblyBegin(pv_rhs_push);  CHKERRQ(ierr);
				ierr = VecAssemblyEnd(pv_rhs_push); CHKERRQ(ierr);

				PetscTime(&cputime_end);
				PetscPrintf(PETSC_COMM_WORLD,"#  Forming stiffness matrixes took %g s \n",cputime_end-cputime_start);

				//=====================
				// LINEAR STOKES SOLVER
				//=====================

				PetscTime(&cputime_start);

				// Compute the rhs
				ierr = ComputeRHS( C, user.DA_Vel, user.DA_Pres, rhs, Pressure, user); CHKERRQ(ierr);

				// Add the additions because of symmetrizing to the RHS
				ierr = VecAXPY(rhs, 	1.0, rhs_add_BC); 	CHKERRQ(ierr);
				ierr = VecAXPY(rhs_p, 	1.0, rhs_add_pBC); 	CHKERRQ(ierr);

				// Add pushing BC - a component of pushing BC is also within ComputeStiffnessMatrixes
				ierr = AddPushingToModel(&user, user.VV_MAT, user.VP_MAT,rhs,rhs_p,pv_rhs_push); CHKERRQ(ierr);
				ierr = VecDestroy(&pv_rhs_push );    CHKERRQ(ierr);

				// Write stiffness matrix to disk, if requested by adding '-DumpStiffnessMatrixes' to the command line (for debugging)
				ierr = WriteStiffnessMatrixToDisk(user.DA_Vel, user.VV_MAT,
						user.VP_MAT, user.PV_MAT, user.PP_MAT, user.approx_S, rhs, rhs_p, sol, Pressure, &user); CHKERRQ(ierr);

				// Compute velocity solution
				// A range of solution techniques have been implemented, such as:

				// 1 - Powell-Hesteness iterations
				// 2 - Schur complement reduction (SCR) method
				// 3 - Fully coupled (FC) approach
				// 4 - No solver, but instead performs Matrix-Vector products like rhs = VV*sol [in order to test the theoretical scalability of LaMEM on parallel machines]

/*
				// DEBUGGING ==============
				{
					if (user.restart==0)
					{
						ierr = DebugSave2Bin_GlobalVec(rhs,										"rhs",							itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(rhs_p,									"rhs_p",						itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(sol,										"sol",							itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.SurfaceTopography,					"SurfaceTopography",			itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.FDSTAG.Center_Density,				"Center_Density",				itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.FDSTAG.Center_EffectiveViscosity,	"Center_EffectiveViscosity",	itime, 		"Initial");CHKERRQ(ierr);


						ierr = DebugSave2Bin_GlobalMat(user.VV_MAT,								"VV",							itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.PP_MAT,								"PP",							itime,	 	"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.PV_MAT,								"PV",							itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.VP_MAT,								"VP",							itime, 		"Initial");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.approx_S,							"approx_S",						itime, 		"Initial");CHKERRQ(ierr);

					}
					else
					{
						ierr = DebugSave2Bin_GlobalVec(rhs,										"rhs",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(rhs_p,									"rhs_p",						itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(sol,										"sol",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.SurfaceTopography,					"SurfaceTopography",			itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.FDSTAG.Center_Density,				"Center_Density",				itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalVec(user.FDSTAG.Center_EffectiveViscosity,	"Center_EffectiveViscosity",	itime, 		"Restart");CHKERRQ(ierr);

						ierr = DebugSave2Bin_GlobalMat(user.VV_MAT,								"VV",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.PP_MAT,								"PP",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.PV_MAT,								"PV",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.VP_MAT,								"VP",							itime, 		"Restart");CHKERRQ(ierr);
						ierr = DebugSave2Bin_GlobalMat(user.approx_S,							"approx_S",						itime, 		"Restart");CHKERRQ(ierr);

					}

				}
				// ========================
*/

				if (user.VelocityTest == PETSC_TRUE)
				{
					PetscPrintf(PETSC_COMM_WORLD,"#  ==================================================\n");
					PetscPrintf(PETSC_COMM_WORLD,"#  *** Testing performance of the Velocity Solver ***\n");
					PetscPrintf(PETSC_COMM_WORLD,"#  ==================================================\n");

					// test performance of a single velocity solve
					ierr = VelSolverTest(user.VV_MAT, sol, rhs, &user); CHKERRQ(ierr);

					// finalize computation
					goto cleanup;
				}

				if(user.StokesSolver==1)
				{
					// Powell-hesteness iterations
					ierr = VecSet(Pressure, 0.0); CHKERRQ(ierr);
					ierr = StokesSolve_PowellHestenes(
							user.VV_MAT, 	user.VP_MAT,
							user.PV_MAT, 	user.PP_MAT,
							sol,            Pressure,
							rhs, 			rhs_p); CHKERRQ(ierr);
				}
				else if(user.StokesSolver==2)
				{
					// Schur complement reduction (SCR) method

					PetscInt min_loc, max_loc;
					PetscScalar min_eta, max_eta;

					ierr = VecMax( ViscosityScaling, &max_loc, &max_eta ); CHKERRQ(ierr);
					ierr = VecMin( ViscosityScaling, &min_loc, &min_eta ); CHKERRQ(ierr);
					PetscPrintf( PETSC_COMM_WORLD, "SCR: Viscosity structure info: min(%lld)=%e : max(%lld)=%e : contrast=%e \n", (LLD)max_loc, 1.0/max_eta, (LLD)min_loc, 1.0/min_eta, max_eta/min_eta );

					// Schur complement reduction method with or without stabilization matrix (depends on element)
					if (vpt_element_type == DAVP_Q1Q1 )
					{	PetscPrintf( PETSC_COMM_WORLD, "Using stabilization matrix \n");
					// pass stabilization matrix
					ierr = StokesSolve_SCR(user.VV_MAT, user.VP_MAT, user.PV_MAT, user.PP_MAT, user.approx_S,
							sol, Pressure, rhs, rhs_p, &user); CHKERRQ(ierr);
					}
					else
					{	ierr = StokesSolve_SCR(user.VV_MAT, user.VP_MAT, user.PV_MAT, PETSC_NULL, user.approx_S,
							sol, Pressure, rhs, rhs_p, &user); CHKERRQ(ierr);
					}
				}
				else if (user.StokesSolver==3)
				{
					// Fully coupled Stokes solver

					PetscInt min_loc, max_loc;
					PetscScalar min_eta, max_eta;

					ierr = VecMax( ViscosityScaling, &max_loc, &max_eta ); CHKERRQ(ierr);
					ierr = VecMin( ViscosityScaling, &min_loc, &min_eta ); CHKERRQ(ierr);
					PetscPrintf( PETSC_COMM_WORLD, "FC: Viscosity structure info: min(%lld)=%e : max(%lld)=%e : contrast=%e \n", (LLD)max_loc, 1.0/max_eta, (LLD)min_loc, 1.0/min_eta, max_eta/min_eta );

					if ((vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG))
					{
						// WHY PP_MAT IS NONZERO FOR FDSTAG??

						ierr = StokesSolve_FC(user.VV_MAT, 	user.VP_MAT, user.PV_MAT, user.PP_MAT, user.approx_S,
								sol, Pressure, rhs, rhs_p, &user); CHKERRQ(ierr);
					}
					else
					{	ierr = StokesSolve_FC(user.VV_MAT, 	user.VP_MAT, user.PV_MAT, PETSC_NULL, user.approx_S,
							sol, Pressure, rhs, rhs_p, &user); CHKERRQ(ierr);
					}
				}
				else if (user.StokesSolver==4)
				{
					// We don't actually solve anything, but perform a number of MatVec products instead,
					//	mainly in order to test the scalability of LaMEM on parallel machines.

					PetscInt       numProd, iter;
					PetscRandom    rctx;
					PetscLogDouble time_1, time_2;

					ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
					ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

					// set random numbers to solution
					ierr = VecSetRandom(sol,rctx); CHKERRQ(ierr);

					numProd = 5000;

					PetscOptionsGetInt( PETSC_NULL ,"-numProd", &numProd, PETSC_NULL );
					PetscTime(&time_1);

					for (iter=0; iter<numProd; iter++){
						ierr = MatMult(user.VV_MAT,sol,rhs);		CHKERRQ(ierr); // rhs= sol*VV
					}
					PetscTime(&time_2);
					PetscPrintf(PETSC_COMM_WORLD,"#  Performed %i MatVec products of rhs = VV*sol in %f seconds [change number of products with -numProd] \n",numProd, time_2-time_1);

					// cleanup
					ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
					goto cleanup;

				}

				PetscTime(&cputime_end);
				PetscPrintf(PETSC_COMM_WORLD,"#  Finished Stokes solve after %g s \n",cputime_end-cputime_start);

				// Report info (primarily useful for debugging)
				PetscInt    loc;
				PetscScalar val;
				Vec 	    Div;
				Vec         mRes, prod;

				// compute divergence residual
				ierr = VecDuplicate(Pressure, &Div); CHKERRQ(ierr);
//				ierr = MatMultTranspose(user.VP_MAT,sol,Div); CHKERRQ(ierr);
				ierr = MatMult( user.PV_MAT,sol,Div); CHKERRQ(ierr);
				ierr = VecAXPY( Div, -1.0, rhs_p ); CHKERRQ(ierr);
				ierr = VecAbs(Div); CHKERRQ(ierr);

				// compute momentum residual
				ierr = VecDuplicate(rhs, &mRes);             CHKERRQ(ierr);
				ierr = VecCopy(rhs, mRes);                   CHKERRQ(ierr);
				ierr = VecDuplicate(rhs, &prod);             CHKERRQ(ierr);
				ierr = MatMult(user.VV_MAT, sol, prod);      CHKERRQ(ierr);
				ierr = VecAXPY(mRes, -1.0, prod);            CHKERRQ(ierr);
				ierr = MatMult(user.VP_MAT, Pressure, prod); CHKERRQ(ierr);
				ierr = VecAXPY(mRes, -1.0, prod);            CHKERRQ(ierr);
				ierr = VecAbs(mRes);                         CHKERRQ(ierr);

				PetscPrintf( PETSC_COMM_WORLD, "------------------------------------------\n" );
				PetscPrintf( PETSC_COMM_WORLD, "Stokes solver summary: \n" );

				PetscPrintf( PETSC_COMM_WORLD, "  divergence: \n" );
				ierr = VecMin( Div, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    Div_min = %12.12e \n", val );
				ierr = VecMax( Div, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    Div_max = %12.12e \n", val );
				if ( val>1.0e-3 ) 				// Error detection
				{
					PetscPrintf(PETSC_COMM_WORLD,"  *** Emergency stop! Maximum divergence is way too large; solver did not converge! \n");
					MPI_Abort(PETSC_COMM_WORLD,1);
				}
				ierr = VecNorm( Div, NORM_2, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |Div|_2 = %12.12e \n", val );
				ierr = VecNorm( Div, NORM_1, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |Div|_1 = %12.12e \n", val );

				PetscPrintf( PETSC_COMM_WORLD, "  momentum: \n" );
				ierr = VecMin( mRes, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    mRes_min = %12.12e \n", val );
				ierr = VecMax( mRes, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    mRes_max = %12.12e \n", val );
				ierr = VecNorm( mRes, NORM_2, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |mRes|_2 = %12.12e \n", val );
				ierr = VecNorm( mRes, NORM_1, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |mRes|_1 = %12.12e \n", val );

				PetscPrintf( PETSC_COMM_WORLD, "  velocity: \n" );
				ierr = VecMin( sol, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    u_min   = %12.12e \n", val );
				ierr = VecMax( sol, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    u_max   = %12.12e \n", val );
				ierr = VecNorm( sol, NORM_2, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |u|_2   = %12.12e \n", val );
				ierr = VecNorm( sol, NORM_1, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |u|_1   = %12.12e \n", val );

				PetscPrintf( PETSC_COMM_WORLD, "  pressure: \n" );
				ierr = VecMin( Pressure, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    p_min   = %12.12e \n", val );
				ierr = VecMax( Pressure, &loc, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    p_max   = %12.12e \n", val );
				ierr = VecNorm( Pressure, NORM_2, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |p|_2   = %12.12e \n", val );
				ierr = VecNorm( Pressure, NORM_1, &val ); CHKERRQ(ierr);
				PetscPrintf( PETSC_COMM_WORLD, "    |p|_1   = %12.12e \n", val );
				PetscPrintf( PETSC_COMM_WORLD, "------------------------------------------\n" );

				ierr = VecDestroy(&Div);  CHKERRQ(ierr);
				ierr = VecDestroy(&mRes); CHKERRQ(ierr);
				ierr = VecDestroy(&prod); CHKERRQ(ierr);

				//============================
				// END OF LINEAR STOKES SOLVER
				//============================

				// In FDSTAG 'sol' gets modified in order to make use of FEM-routines. Since we are going to use sol as initial guess
				// in the next timestep or when restarting a simulation we need a unmodified version of sol.
				// This is why we make a copy of sol that will be modified (sol_advect). Later on, we save sol to Breakpointfile NOT sol_advect.
				// (tobib 20.09.12)

				// NOTE:sol_advect is initialized earlier in the code, as Sarah added an option to do a few thermal steps BEFORE the first stokes step. BK, 28.9.2012
				ierr = VecCopy(sol,sol_advect);	CHKERRQ(ierr);

				if( vpt_element_type == DAVP_FDSTAG )
				{
					// If FDSTAG: interpolate velocities from cell-center to corner points, such that we can use
					// the Q1P0 routines to compute particle properties and for advection.
					FDSTAG_InterpolateVelocityPressureToCornerpoints(&user, sol_advect);
				}

				//==========================================================================================
				// Compute properties @ integration points, such as strain-rate, pressure, stress etc.
				// Here, we only compute arrays that are necessary for the nonlinear iterations (for speed)
				// After the iterations are finished, we update all variables (such as finite strain etc.)
				//==========================================================================================

				ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, user.DA_Temp, Temp, &user, user.dt, 0 ); CHKERRQ(ierr);

				if (user.NonlinearIterations==1)
				{
					// Compute error of nonlinear iterations
					{
						PetscScalar err_val;

						// Compute velocity error
						ierr = VecAXPBY(sol_old, 1.0, -1.0 ,sol_advect);   CHKERRQ(ierr); // sol_old = sol_advect-sol_old            (difference between two vectors)
						ierr = VecNorm( sol_old, NORM_2, &velocity_error); CHKERRQ(ierr); // velocity_error = |sol_advect-sol_old|_2 (norm2 of difference)
						ierr = VecNorm( sol_advect, NORM_2, &err_val);     CHKERRQ(ierr); // |sol_advect|_2							 (normalize over average velocity)
						velocity_error = velocity_error/err_val;           CHKERRQ(ierr);
					}

					PetscPrintf(PETSC_COMM_WORLD,"Timestep %lld, NONLINEAR ITERATION %lld [from %lld], accuracy level=%g, current velocity error: %g  \n \n",
							(LLD)itime, (LLD)it_step, (LLD)(user.MaxNonlinearIterations), user.NonlinearIterationsAccuracy, velocity_error);


					ierr = VecCopy(sol_advect, sol_old); CHKERRQ(ierr);

					it_step += 1;

					if (it_step>user.MaxNonlinearIterations)
					{
						velocity_error = 0.0;
						PetscPrintf(PETSC_COMM_WORLD," Aborting nonlinear iterations since maximum number of iterations is exceeded. \n \n \n");
					}

				}
				else
				{
					// Not performing nonlinear iterations
					velocity_error = 0.0;
				}

				if(it_step > 2)
				{
					user.PlasticityCutoff =	0;  // only during the first iterations of the first time step, this should be 1
				}


			}
			//============================
			// END OF NONLINEAR ITERATIONS
			//============================

			PetscTime(&cputime_end);
			PetscPrintf(PETSC_COMM_WORLD,"#  Solving for velocity took %g s \n",cputime_end - cputime_start_nonlinear);

		}
		//==========================================
		// END OF NONLINEAR THERMO-MECHANICAL SOLVER
		//==========================================

		//=========================================================================================
		// In case we perform an analytical benchmark, compare numerical and analytical results
		if (user.AnalyticalBenchmark)
		{

			// Update properties @ integration points (since for speed reasons we do not compute
			// all strainrate and stress components in the nonlinear iterations above).
			ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, PETSC_NULL, PETSC_NULL, &user, user.dt, 1 ); CHKERRQ(ierr);

			ierr = LaMEM_CompareNumerics_vs_AnalyticalSolution(&user, C, sol_advect, Pressure); CHKERRQ(ierr);

			// update INTP props once more (to correctly set P values)
			ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, user.DA_Temp, Temp, &user, user.dt, 1 ); CHKERRQ(ierr);
		}
		//==========================================================================================


		//=============================================================================================
		// debugging & testing the particle interpolation routine
/*
		// props @ intp
		IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, DA_TEMP, Temp, &user, user.dt, 1 );

		// write 2 disc
		WriteOutputFileMatlab(&user, user.DA_Vel, DA_TEMP, sol_advect, Temp, 111, C->ElementType, C->ngp_vel);

		// compute @ particles
		ComputePropertiesAtParticles(C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, DA_TEMP, Temp, Temp_old, &user, user.dt);

		// average to grid
		MaterialPropertiesFromTracers(C,&user, user.DA_Vel);

		// write 2 disc
		WriteOutputFileMatlab(&user, user.DA_Vel, DA_TEMP, sol_advect, Temp, 112, C->ElementType, C->ngp_vel);

		MPI_Abort(PETSC_COMM_WORLD,1);
*/
		//=============================================================================================


		//==========================================================================================
		// TEMPERATURE SOLVER WARNING! UNCOUPLED TEMPERATURE SOLUTION!
		//==========================================================================================
#ifdef TEMPERATURE


		if (user.temp_initialize==1)
		{
			user.dt = user.dt_temp;
		}

		// Compute temperature stiffness matrix
		PetscTime(&cputime_start);

		ierr = MatZeroEntries(user.TEMP_MAT); CHKERRQ(ierr); // Zero the stiffness matrix

		// Compute Temperature stiffness matrix
		ierr = ComputeStiffnessMatrixRHSTemperature( C, user.DA_Temp, user.DA_Vel, user.TEMP_MAT, user,
				Temp_local, Temp, rhs_Temp_local, rhs_Temp, user.dt ); 	CHKERRQ(ierr);

		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD,"#  Forming temperature stiffness matrixes took %g s \n",cputime_end-cputime_start);

		// Set stiffness matrix for temperature
		if (itime==user.time_start)
		{
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_temp); CHKERRQ(ierr);
		}
		ierr = KSPSetOptionsPrefix( ksp_temp, "temp_"); CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp_temp); CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp_temp,user.TEMP_MAT,user.TEMP_MAT); CHKERRQ(ierr);
		ierr = KSPSetInitialGuessNonzero(ksp_temp,PETSC_TRUE); CHKERRQ(ierr);

		// Compute solution; as diffusion is an 'easy' (and fast to compute) problem, we don't do anything fancy here
		PetscTime(&cputime_start);
		ierr = VecCopy(Temp, Temp_old);	         CHKERRQ(ierr);
		ierr = KSPSolve(ksp_temp,rhs_Temp,Temp); CHKERRQ(ierr);
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD,"#  Solving for temperature took %g s \n", cputime_end - cputime_start);

		ierr = WriteTemperatureStiffnessMatrixToDisk( user.DA_Temp, user.TEMP_MAT, rhs_Temp, Temp, &user); CHKERRQ(ierr);  // write to disk if -DeumpStiffnessMatrixes is added to the command-line
#endif
		//==========================================================================================

		//==========================================================================================
		// Update all properties @ integration points and compute properties @ particles
		PetscTime(&cputime_start);
#ifdef TEMPERATURE
		ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, user.DA_Temp, Temp, &user, user.dt, 1 ); CHKERRQ(ierr);
#ifdef PARTICLES
		if (user.ParticleInput==1)
		{
			ierr = ComputePropertiesAtParticles(C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, user.DA_Temp, Temp, Temp_old, &user, user.dt); CHKERRQ(ierr);
		}
#endif
#else
		ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, PETSC_NULL, PETSC_NULL, &user, user.dt, 1 ); CHKERRQ(ierr);
#ifdef PARTICLES
		if (user.ParticleInput==1)
		{
			ierr = ComputePropertiesAtParticles(C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, PETSC_NULL, PETSC_NULL, PETSC_NULL, &user, user.dt); CHKERRQ(ierr);
		}
#endif
#endif
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD,"#  Updating properties @ integration points and particles took %g s \n",cputime_end-cputime_start);

		//==========================================================================================


		//==========================================================================================
		// Advect the grid WARNING! EXPLICIT ADVECTION!
		//==========================================================================================

		PetscPrintf(PETSC_COMM_WORLD,"AdvectionMethod = %lld \n",(LLD)(user.GridAdvectionMethod));

		if (user.GridAdvectionMethod > 0){
			ierr = AdvectGrid(user.DA_Vel, sol_advect, user.DA_Vel, &user); CHKERRQ(ierr);
#ifdef TEMPERATURE
			ierr = AdvectGrid(user.DA_Vel, sol_advect, user.DA_Temp, &user); CHKERRQ(ierr);
#endif

		}

		// Print some info
		ierr = VecMax(sol_advect, PETSC_NULL, &MaxVel);	CHKERRQ(ierr);
		ierr = VecMin(sol_advect, PETSC_NULL, &MinVel); CHKERRQ(ierr);
		MaxVel = PetscMax(MaxVel, -MinVel);

		// Error detection
		if ( isnan(MaxVel) )
		{
			PetscPrintf(PETSC_COMM_WORLD,"  *** Emergency stop! Maximum velocity is NaN ***  \n");
			MPI_Abort(PETSC_COMM_WORLD,1);
		}


		//==========================================================================================
		// Advect tracers
#ifdef PARTICLES
		if (user.ParticleInput==1){
			PetscPrintf(PETSC_COMM_WORLD,"  Starting particle advection \n");
			MPI_Barrier(PETSC_COMM_WORLD);

			ierr = AdvectParticles(C,user.DA_Vel, &user, sol_advect, user.dt); CHKERRQ(ierr);

			PetscPrintf(PETSC_COMM_WORLD,"  Finished particle advection \n");
			MPI_Barrier(PETSC_COMM_WORLD);
		}
#endif
		//==========================================================================================

		//==========================================================================================
		// Update free surface information
		if ( vpt_element_type != DAVP_FDSTAG )
		{
			ierr = UpdateSurfaceAndBottomTopography_FEM( &user ); CHKERRQ(ierr); // Update surface & bottom topographies in case we use a FEM approach
		}


		if ( user.ErosionParameters.UseInternalFreeSurface == 1)
		{
			// In case we have an internal free surface, advect it
			ierr = AdvectAndUpdate_InternalFreeSurface( &user, sol_advect ); CHKERRQ(ierr);	 // Update internal free surface

			// Apply sedimentation to the internal free surface
			ierr = ApplySedimentationToInternalFreeSurface(&user); // Apply sedimentation to internal free surface if requested

			if(user.ErosionParameters.ErosionModel == 1)
			{
				ierr = ApplyInfinitelyFastErosionToInternalFreeSurface(&user); // Apply fast erosion to internal free surface if requested
			}
			else if (user.ErosionParameters.ErosionModel == 2)
			{
				// use FD based erosion code [serial only!]
				ierr = ApplySerialErosionCodeToInternalFreeSurface(&user); CHKERRQ(ierr);
			}

			ierr = CorrectPhasesForInternalFreeSurface( &user); CHKERRQ(ierr);
		}
		//==========================================================================================

		//==========================================================================================
		// Advect grid with BG deformation rate in case we use an FDSTAG approach with Exx=Eyy!=0
		if (vpt_element_type == DAVP_FDSTAG)
		{
			ierr = DeformFDSTAGMeshWithBackgroundStrainrate(&user);	CHKERRQ(ierr);
		}
		//==========================================================================================


		// Update time state
		user.time = user.time + user.dt;

		// Compute and store interesting properties such as Vrms
		ierr = ComputeGlobalProperties(user.DA_Vel, &user, itime, sol_advect, C); CHKERRQ(ierr);


		//==========================================================================================
		// Calculate gravity field + misfit

		if (user.GravityField.GetIt==1)
		{
			PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- START evaluating gravity field --------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- Check your results for reliability ----\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- This is a development version ---------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");

			ierr = GetGravityField(&user,C->ngp_vel,itime); CHKERRQ(ierr);

			PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- END evaluating gravity field ----------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
		}

		//==========================================================================================
		// Extract surface velocity + calculate misfit

		if (user.SurfVelField.GetIt==1)
		{
			PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- START extracting surface velocity -----\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- This is a development version ---------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- Check your results for reliability ----\n");
			PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");

			if (user.ErosionParameters.UseInternalFreeSurface == 1){
				//ierr = GetSurfaceVelocity_FreeSurface(&user,itime);CHKERRQ(ierr);
				ierr = GetSurfaceVelocity_InternalFreeSurface(&user,itime);CHKERRQ(ierr);
			}
			else{
				ierr = GetSurfaceVelocity(&user,user.DA_Vel,sol_advect);CHKERRQ(ierr);
			}

			PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#-- END extracting surface velocity -------\n");
			PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
		}


		if(user.Optimisation.GetIt==1)
		{
			if (user.SurfVelField.GetIt==1)
			{
				PetscPrintf(PETSC_COMM_WORLD,"==> Misfit SurfVelField (L1)/N = %g \n",user.Optimisation.SumAbsVel/user.Optimisation.NSurfVel);
				PetscPrintf(PETSC_COMM_WORLD,"==> Misfit SurfVelField (L2)/sqrt(N) = %g \n",sqrt(user.Optimisation.SumSqrsVel/user.Optimisation.NSurfVel));
			}
			else
			{
				user.Optimisation.SumAbsVel  = 0.0;
				user.Optimisation.SumSqrsVel = 0.0;
				user.Optimisation.NSurfVel   = 0.0;
			}
			if(user.GravityField.GetIt==1)
			{
				PetscPrintf(PETSC_COMM_WORLD,"==> Misfit GravityField (L1)/N = %g \n",user.Optimisation.SumAbsGrav/user.Optimisation.NGrav);
				PetscPrintf(PETSC_COMM_WORLD,"==> Misfit GravityField (L2)/sqrt(N) = %g \n",sqrt(user.Optimisation.SumSqrsGrav/user.Optimisation.NGrav));
			}
			else
			{
				user.Optimisation.SumAbsGrav  = 0.0;
				user.Optimisation.SumSqrsGrav = 0.0;
				user.Optimisation.NGrav       = 0.0;
			}

			PetscScalar OptMfit_L1,OptMfit_L2,OptMfit_ChiSquare;

			// Misfit based on L1 norm (arithmetic mean)
			OptMfit_L1 = (user.Optimisation.SumAbsGrav + user.Optimisation.SumAbsVel)/(user.Optimisation.NGrav+user.Optimisation.NSurfVel);

			// Misfit based on L2 norm (rms mean)
			OptMfit_L2 = sqrt((user.Optimisation.SumSqrsGrav + user.Optimisation.SumSqrsVel)/(user.Optimisation.NGrav+user.Optimisation.NSurfVel));

			// Chi-square misfit based
			OptMfit_ChiSquare = (user.Optimisation.SumSqrsGrav + user.Optimisation.SumSqrsVel)/(user.Optimisation.NGrav+user.Optimisation.NSurfVel);

			PetscPrintf(PETSC_COMM_WORLD,"==> Misfit (L1)/N = %g \n",OptMfit_L1);
			PetscPrintf(PETSC_COMM_WORLD,"==> Misfit (L2)/sqrt(N) = %g \n",OptMfit_L2);
			PetscPrintf(PETSC_COMM_WORLD,"==> Misfit chi-square = %g \n",OptMfit_ChiSquare);

			*LaMEM_OutputParameters = OptMfit_ChiSquare;

			PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
		}

		//==========================================================================================
		// Save data to disk

		if(user.save_timesteps != 0)
		{
			LaMEMMod( itime, user.save_timesteps, &SaveOrNot);
		}
		else
		{
			SaveOrNot = 2;
		}

		if(SaveOrNot == 0)
		{
			LaMEMView_QuadratureFields 	view;
			char 						*DirectoryName = NULL;

			// create directory
			if ((user.MatlabOutputFiles==1) || (user.VTKOutputFiles==1) || (user.AVDPhaseViewer))
			{
				asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
				ierr = LaMEM_CreateOutputDirectory(DirectoryName); CHKERRQ(ierr);
			}

			// --- Matlab output (a): regular files ---
			if (user.MatlabOutputFiles==1)
			{
				ierr = WriteOutputFileMatlab(&user, user.DA_Vel, user.DA_Temp, sol_advect, Temp, itime, C->ElementType, C->ngp_vel, DirectoryName); CHKERRQ(ierr);
			}

			// --- Paraview output (a): velocity temperature, topography ---
			if (user.VTKOutputFiles==1)
			{
				if (user.Output.velocity==1)
				{
				// Velocity mesh
				ierr = WriteVelocityOutputFile_VTS(&user,itime, sol_advect, DirectoryName); CHKERRQ(ierr);
				}
#ifdef TEMPERATURE
				if (user.Output.temperature==1)
				{
				// Temperature mesh
				ierr = WriteTemperatureOutputFile_VTS(&user,itime, Temp, DirectoryName); CHKERRQ(ierr);
				}
#endif
				// Surface, bottom topography and internal erosion surface if required
				ierr = WriteTopographyOutputFile_VTS(&user, itime, DirectoryName); CHKERRQ(ierr);

				if (user.Output.quadrature==1)
				{
				// Save VTS output @ quadrature points if required
				ierr = LaMEMView_QuadratureFieldsInit(&view); CHKERRQ(ierr);

				ierr = LaMEMViewQuadraturePoints_3DPVTS(&view,&user,user.DA_Processors,user.DA_Materials,user.Materials,C->nintp_1D,itime, DirectoryName); CHKERRQ(ierr);
				}
			}

			// --- Paraview output (b): phases ---
			if (user.AVDPhaseViewer)
			{	// Creates a Voronoi diagram and writes the phases as VTS files
				ierr = WritePhasesOutputFile_VTS(C, &user,itime, DirectoryName); CHKERRQ(ierr);
			}

#ifdef PARTICLES
			// --- Matlab output (b): particles (mainly debugging purposes) ---
			if (user.SaveParticles==1)
			{
				ierr = WriteParticlesToDisc(&user, itime); CHKERRQ(ierr);
			}
#endif

			// clean up
			if(DirectoryName) free(DirectoryName);
		}

		//==========================================================================================
		// compute length of the next timestep
		ierr     = VecMax(sol_advect, PETSC_NULL, &MaxVel); CHKERRQ(ierr);
		ierr     = VecMin(sol_advect, PETSC_NULL, &MinVel); CHKERRQ(ierr);
		MaxVel   = PetscMax(MaxVel, -MinVel);
		ierr     = DMDAGetInfo(user.DA_Vel, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
		dx       = user.W/((double) mx);
		dy       = user.L/((double) my);
		dz       = user.H/((double) mz);
		dmin     = PetscMin(dx,dy);
		dmin     = PetscMin(dmin,dz);
		fac_time = 1.0;
		if(itime < 10)
		{
			fac_time = ((double) itime+1)/((double) 10);
		}

		// compute Courant timestep
		user.dt = user.CFL*fac_time*dmin/MaxVel;

		// limit time step
		if (user.dt > user.dt_max) user.dt=user.dt_max;
		//==========================================================================================

		//==========================================================================================
		// Remesh grid if required
		if((user.GridAdvectionMethod == 2) || (user.remesh==1))
		{
			PetscPrintf(PETSC_COMM_WORLD," Start remeshing grid \n");

			ierr = RemeshGrid(user.DA_Vel, &user); CHKERRQ(ierr);

			ierr = UpdateSurfaceAndBottomTopography_FEM( &user ); CHKERRQ(ierr);

			ierr = SaveInitialMesh(&user,user.DA_Vel,"MeshAfterRemeshing"); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD," Finished remeshing grid \n");
		}

		//==========================================================================================
		// Map tracers to new grid in case we do NOT use a lagrangian mesh
#ifdef PARTICLES
		if ( ((user.GridAdvectionMethod != 1) && (user.ParticleInput==1)) || ((user.remesh==1) && (user.ParticleInput==1)) )
		{
			PetscPrintf(PETSC_COMM_WORLD," Starting GetParticleNodes \n");
			ierr = GetParticleNodes(C, user.DA_Vel, &user); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD," Finished GetParticleNodes \n");
		}
#endif
		//==========================================================================================

		//==========================================================================================
		// Perform phase transitions on particles (i.e. change the phase number of particles)
#ifdef PARTICLES
		//	ParticlePhaseTransitions(&user);
#endif
		//==========================================================================================

		//==========================================================================================
		// Print some information on screen
		PetscTime(&cputime_end);
		if (user.Characteristic.Length>1)
		{
			PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user.Characteristic.Velocity*user.Characteristic.cmYear);
			if (user.DimensionalUnits==1) { PetscPrintf(PETSC_COMM_WORLD," [cm/year]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user.Characteristic.Velocity);
			if (user.DimensionalUnits==1) { PetscPrintf(PETSC_COMM_WORLD," [m/s]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
		}

		if (user.Characteristic.Length>1)
		{
			// Most likely a setup that runs in natural lengthscales
			PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ",user.time*user.Characteristic.Time/user.Characteristic.SecYear, user.dt*user.Characteristic.Time/user.Characteristic.SecYear);
			if (user.DimensionalUnits==1) { PetscPrintf(PETSC_COMM_WORLD," [Years]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
		}
		else
		{
			// Lab timescale or non-dimensional units
			PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ",user.time*user.Characteristic.Time, user.dt*user.Characteristic.Time);
			if (user.DimensionalUnits==1) { PetscPrintf(PETSC_COMM_WORLD," [s]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
		}

		PetscPrintf(PETSC_COMM_WORLD,"# Finished timestep %lld out of %lld in %g s\n",(LLD)itime, (LLD)(user.time_end), cputime_end - cputime_start_tstep);
		PetscPrintf(PETSC_COMM_WORLD,"  \n");

		//==========================================================================================

		//==========================================================================================
		// Create breakpoint files, for restarting the code
		if (user.save_breakpoints > 0)
		{
			// Standard breakpoint file
			LaMEMMod( itime, user.save_breakpoints, &SaveOrNot);
			if ((SaveOrNot==0))
			{
				ierr = SaveBreakPoint(&user,user.DA_Vel, user.DA_Pres, sol, Temp, Pressure, itime, 0); CHKERRQ(ierr);
			}

			// Incremental breakpoint file
			LaMEMMod(itime, user.incr_breakpoints, &SaveOrNot);
			if ((itime != 0) && (SaveOrNot==0))
			{
				ierr = SaveBreakPoint(&user,user.DA_Vel, user.DA_Pres, sol, Temp, Pressure, itime, user.break_point_number); CHKERRQ(ierr);
				user.break_point_number = user.break_point_number+1;
			}
		}
		//==========================================================================================

	}
	//======================
	// END OF TIME STEP LOOP
	//======================

	// Write info to disk
	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"# Total time required: %g s \n",cputime_end - cputime_start0);

	// tidy up
cleanup:

#ifdef TEMPERATURE
	if( ksp_temp != PETSC_NULL ) { KSPDestroy( &ksp_temp ); }
#endif

	// DM's
	ierr = DMDestroy(&user.DA_Vel);        CHKERRQ(ierr);
	ierr = DMDestroy(&user.DA_Pres);       CHKERRQ(ierr);
	ierr = DMDestroy(&user.DA_Materials);  CHKERRQ(ierr);
	ierr = DMDestroy(&user.DA_Processors); CHKERRQ(ierr);

	if (user.VTKOutputFiles==1)
	{
		ierr = DMDestroy( &user.DA_Quadrature ); CHKERRQ(ierr);
	}
	ierr = DMDestroy(&user.DA_SurfaceTopography); CHKERRQ(ierr);
	ierr = DMDestroy(&user.DA_BottomTopography);  CHKERRQ(ierr);


#ifdef TEMPERATURE
	ierr = DMDestroy( &user.DA_Temp ); CHKERRQ(ierr);
#endif

	// Vectors
	ierr = VecDestroy(&user.Materials);            CHKERRQ(ierr);
	ierr = VecDestroy(&rhs );                      CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_local );                CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_add_BC );               CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_add_BC_loc );           CHKERRQ(ierr);
	ierr = VecDestroy(&user.SurfaceTopography);    CHKERRQ(ierr);
	ierr = VecDestroy(&user.SurfaceTopography_Vx); CHKERRQ(ierr);
	ierr = VecDestroy(&user.SurfaceTopography_Vy); CHKERRQ(ierr);
	ierr = VecDestroy(&user.SurfaceTopography_Vz); CHKERRQ(ierr);
	ierr = VecDestroy(&user.BottomTopography);     CHKERRQ(ierr);

	// Cleanup FD erosion code
	if ((user.ErosionParameters.ErosionModel==2) && (rank==0))
	{
		// not good! implement erosion-code-data-structure create and destroy functions
		ierr = DMDestroy (&user.ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode);CHKERRQ(ierr);
		ierr = VecDestroy(&user.ErosionParameters.FE_ErosionCode.ErosionSurface);   CHKERRQ(ierr);
	}

	if (vpt_element_type == DAVP_FDSTAG)
	{
		for (i=0; i<user.num_phases; i++)
		{
			ierr = VecDestroy(&user.FDSTAG.Center_PhaseProportions[i]);		 CHKERRQ(ierr);
			ierr = VecDestroy(&user.FDSTAG.Corner_PhaseProportions[i]);		 CHKERRQ(ierr);
			ierr = VecDestroy(&user.FDSTAG.Corner_PhaseProportions_local[i]);CHKERRQ(ierr);
			ierr = VecDestroy(&user.FDSTAG.XYPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
			ierr = VecDestroy(&user.FDSTAG.XZPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
			ierr = VecDestroy(&user.FDSTAG.YZPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
		}
		ierr = VecDestroy(&user.FDSTAG.Center_Temperature);                  CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_Pressure);                     CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_Strain);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_PlasticStrain);                CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_T2nd);                         CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_E2nd);                         CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_EffectiveViscosity);           CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_Density);                      CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Center_NumParticles);                 CHKERRQ(ierr);

		ierr = VecDestroy(&user.FDSTAG.Corner_Temperature);                  CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Corner_Pressure);                     CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Corner_Density);                      CHKERRQ(ierr);

		ierr = VecDestroy(&user.FDSTAG.XYPoints_Temperature);                CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_Pressure);                   CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_Strain);                     CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_PlasticStrain);              CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_T2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_E2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_EffectiveViscosity);         CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XYPoints_Density);                    CHKERRQ(ierr);

		ierr = VecDestroy(&user.FDSTAG.XZPoints_Temperature);                CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_Pressure);                   CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_Strain);                     CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_PlasticStrain);              CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_T2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_E2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_EffectiveViscosity);         CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.XZPoints_Density);                    CHKERRQ(ierr);

		ierr = VecDestroy(&user.FDSTAG.YZPoints_Temperature);                CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_Pressure);                   CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_Strain);                     CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_PlasticStrain);              CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_T2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_E2nd);                       CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_EffectiveViscosity);         CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.YZPoints_Density);                    CHKERRQ(ierr);

		ierr = DMDestroy(&user.FDSTAG.DA_CENTER);                            CHKERRQ(ierr);
		ierr = DMDestroy(&user.FDSTAG.DA_CORNER);                            CHKERRQ(ierr);
		ierr = DMDestroy(&user.FDSTAG.DA_XY_POINTS);                         CHKERRQ(ierr);
		ierr = DMDestroy(&user.FDSTAG.DA_XZ_POINTS);                         CHKERRQ(ierr);
		ierr = DMDestroy(&user.FDSTAG.DA_YZ_POINTS);                         CHKERRQ(ierr);
	}

	// particles
	ierr = PetscFree(user.ParticlesLocal); CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = VecDestroy(&Temp);           CHKERRQ(ierr);
	ierr = VecDestroy(&Temp_old);       CHKERRQ(ierr);
	ierr = VecDestroy(&Temp_local);     CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_Temp);       CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_Temp_local); CHKERRQ(ierr);

	if (vpt_element_type == DAVP_FDSTAG)
	{
		ierr = VecDestroy(&user.FDSTAG.Corner_HeatCapacity);	CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Corner_Conductivity);	CHKERRQ(ierr);
		ierr = VecDestroy(&user.FDSTAG.Corner_RadioactiveHeat);	CHKERRQ(ierr);
	}

#endif

	ierr = VecDestroy(&Pressure);        CHKERRQ(ierr);
	ierr = VecDestroy(&Pressure_local);  CHKERRQ(ierr);

	ierr = VecDestroy(&rhs_p);           CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_p_local);     CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_add_pBC);     CHKERRQ(ierr);
	ierr = VecDestroy(&rhs_add_pBC_loc); CHKERRQ(ierr);

	ierr = VecDestroy(&ViscosityScaling); CHKERRQ(ierr);
	ierr = VecDestroy(&sol);              CHKERRQ(ierr);
	ierr = VecDestroy(&sol_old);          CHKERRQ(ierr);
	ierr = VecDestroy(&sol_local);        CHKERRQ(ierr);
	ierr = VecDestroy(&sol_advect);       CHKERRQ(ierr);

	// Matrices
	ierr = MatDestroy(&user.VV_MAT); CHKERRQ(ierr);
	ierr = MatDestroy(&user.PP_MAT); CHKERRQ(ierr);

	ierr = MatDestroy(&user.VP_MAT); CHKERRQ(ierr);
	ierr = MatDestroy(&user.PV_MAT); CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = MatDestroy(&user.TEMP_MAT); CHKERRQ(ierr);
#endif
	ierr = MatDestroy(&user.approx_S); CHKERRQ(ierr);

	ierr = LaMEMVelPressureDADestroy(&C); CHKERRQ(ierr);


	ierr = PetscFree(user.ParticlesLocal);    CHKERRQ(ierr);
	ierr = PetscFree(user.TimeDependentData); CHKERRQ(ierr);

	if(user.InputParamFile)
	{
		MaterialDestroy(user.PhaseMaterialProperties);
	}

	PetscFunctionReturn(0);

}
//==========================================================================================================
// END OF LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

