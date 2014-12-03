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

$Id: LaMEMLib.c 4815 2013-10-06 23:00:33Z lapopov $
$Date:$

This version is compatible with PETSc version 3.5
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Version.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "Utils.h"
#include "LaMEM_Initialize.h"
#include "LaMEMLib_FDSTAG_private.h"

//#include "Breakpoint.h"
//#include "LaMEM_FE_ErosionCode.h"
//#include "LaMEM_AnalyticalSolutions.h"

#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "interface.h"
#include "multigrid.h"
#include "advect.h"
#include "marker.h"
#include "input.h"

//==========================================================================================================
// LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

#undef __FUNCT__
#define __FUNCT__ "LaMEMLib_FDSTAG"
PetscErrorCode LaMEMLib_FDSTAG(PetscBool InputParamFile, const char *ParamFile, PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id)
{

	PetscBool          done;
	UserContext        user;
	PetscInt           SaveOrNot;
//	PetscLogDouble     cputime_start, cputime_start0, cputime_end, cputime_start_tstep, cputime_start_nonlinear;
	PetscLogDouble     cputime_start, cputime_start_nonlinear, cputime_end;

	FDSTAG   fs;    // staggered-grid layout
	BCCtx    bc;    // boundary condition context
	JacRes   jr;    // Jacobian & residual context
	AdvCtx   actx;  // advection context
	PMat     pm;    // preconditioner matrix
	PCStokes pc;    // Stokes preconditioner
	SNES     snes;  // PETSc nonlinear solver
	NLSol    nl;    // nonlinear solver context
	PVOut    pvout; // paraview output driver

//	PetscViewer  viewer;
//	PetscBool    do_restart;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscTime(&cputime_start);

	// Start code
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"              Current Revision: %s - %s		     \n",__SVNREVISION__,__SVNDATE__);
	PetscPrintf(PETSC_COMM_WORLD,"  Modified items: %s\n",__SVNMANCHANGES__);
	PetscPrintf(PETSC_COMM_WORLD,"     Compiled: Date: %s - Time: %s - Mode:  %s	    \n",__DATE__,__TIME__,__OPTMODE__);
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"        STAGGERED-GRID FINITE DIFFERENCE CANONICAL IMPLEMENTATION           \n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");

	if(LaMEM_OutputParameters) LaMEM_OutputParameters = NULL;

	ierr = FDSTAGClear(&fs);    CHKERRQ(ierr);
	ierr = BCClear    (&bc);    CHKERRQ(ierr);
	ierr = JacResClear(&jr);    CHKERRQ(ierr);
	ierr = ADVClear   (&actx);  CHKERRQ(ierr);
	ierr = NLSolClear (&nl);    CHKERRQ(ierr);
	ierr = PVOutClear (&pvout); CHKERRQ(ierr);

	// Initialize context
	ierr = PetscMemzero(&user, sizeof(UserContext)); CHKERRQ(ierr);

	// set input file flag & name
	user.InputParamFile = InputParamFile;
	if(InputParamFile == PETSC_TRUE) strcpy(user.ParamFile, ParamFile);

	// initialize variables
	ierr = FDSTAGInitCode(&jr, &user); CHKERRQ(ierr);

	// Give current LaMEM session a specific group ID
	user.Optimisation.mpi_group_id = *mpi_group_id;

	//========================================================================================
	// few parameters that can be reused
	//========================================================================================

//	user.save_breakpoints = -1;
//	user.BC.InternalBound = -1;
//	user.restart = 0;
//	user.ErosionParameters.ErosionModel = 0;
//	user.ErosionParameters.UseInternalFreeSurface = 0;
//	user.SavePartitioning = 0;
//	user.LoadInitialParticlesFromDisc = 0;
//	user.remesh = 0;
//	user.InitialMeshFromFile == 1
//	user.Setup.Model == 3
//	user.EulerianAfterTimestep > 0
//	user.AnalyticalBenchmark == PETSC_TRUE
//	user.GridAdvectionMethod
//	user.NonlinearIterations==1
//	user.InitialErosionSurfaceFromFile == 1
//	user.ParticleInput == 1
//	user.fileno
//	user.time_start
//	user.time_end
//	user.time
//	user.dt
//	user.MatlabOutputFiles == 1
//	user.VTKOutputFiles == 1
//	user.AVDPhaseViewer
//  user.PlasticityCutoff

	//========================================================================================
	// Setting up the solver
	//========================================================================================

//	PetscTime(&cputime_start);

	// Initialize erosion solver and surface if required
//	ierr = InitializeInternalErosionSurfaceOnRankZero(&user); CHKERRQ(ierr);

	// Read/store erosion surface
//	if(user.ErosionParameters.ErosionModel == 2)
//	{
//		if (user.InitialErosionSurfaceFromFile == 1)
//		{
//			ierr = LoadInitialErosionSurface(&user); CHKERRQ(ierr);
//		}
//		ierr = SaveInitialErosionSurface(&user,"InitialErosionSurface"); CHKERRQ(ierr);
//	}


//	if (user.BC.InternalBound > 0)
//	{
//		DefineInternalBC(&user);
//	}


	// Read the initial mesh from file if asked for
//	if(user.InitialMeshFromFile == 1)
//	{
//		ierr = ReadMeshFromFile(user.DA_Vel, &user); CHKERRQ(ierr);
//		PetscPrintf(PETSC_COMM_WORLD,"# Successful read initial mesh from disc \n");
//	}

	// Read particles and properties from a breakpoint file if this was requested and if this is possible
//	if (user.restart == 1)
//	{
//		ierr = LoadBreakPoint(C, &user,user.DA_Vel, sol, Temp, Pressure, 0); CHKERRQ(ierr);
//	}

	//======================
	// SETUP DATA STRUCTURES
	//======================

	// initialize scaling object
	ierr = ScalingCreate(
		&jr.scal,
		user.DimensionalUnits,
		user.Characteristic.kg,
		user.Characteristic.Time,
		user.Characteristic.Length,
		user.Characteristic.Temperature,
		user.Characteristic.Force); CHKERRQ(ierr);

	// initialize material parameter limits
	ierr = SetMatParLim(&jr.matLim, &user); CHKERRQ(ierr);

	// initialize time stepping parameters
	ierr = TSSolSetUp(&jr.ts, &user); CHKERRQ(ierr);

	// create staggered grid object
	ierr = FDSTAGCreate(&fs, user.nnode_x, user.nnode_y, user.nnode_z,
		user.cpu_x, user.cpu_y, user.cpu_z); CHKERRQ(ierr);

	// generate coordinates of grid nodes/cells
	ierr = FDSTAGGenCoord(&fs, &user); CHKERRQ(ierr);

	// save processor partitioning
	if(user.SavePartitioning)
	{
		ierr = FDSTAGProcPartitioning(&fs, &user); CHKERRQ(ierr);

		// return immediately
		ierr = FDSTAGDestroy(&fs); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// print essential grid details
	ierr = FDSTAGView(&fs); CHKERRQ(ierr);

	// create boundary condition context
	ierr = BCCreate(&bc, &fs, &jr.ts, &jr.scal); CHKERRQ(ierr);

	// set background strain-rates
	ierr = BCSetStretch(&bc, &user); CHKERRQ(ierr);

	// set pushing block parameters
	ierr = BCSetPush(&bc, &user); CHKERRQ(ierr);

	// create Jacobian & residual evaluation context
	ierr = JacResCreate(&jr, &fs, &bc, user.num_phases, 0); CHKERRQ(ierr);

	// WARNING! NO TEMPERATURE! Set local temperature vector to unity (ad-hoc)
	ierr = VecSet(jr.lT, 1.0); CHKERRQ(ierr);

	// initialize material properties
	ierr = InitMaterialProps(&jr, &user); CHKERRQ(ierr);

	// initialize gravity acceleration
	jr.grav[0] = 0.0;
	jr.grav[1] = 0.0;
	jr.grav[2] = user.Gravity;

	// create advection context
	ierr = ADVCreate(&actx, &fs, &jr); CHKERRQ(ierr);

	// initialize markers
	ierr = ADVMarkInit(&actx, &user); CHKERRQ(ierr);

	// create Stokes preconditioner & matrix
	ierr = PMatCreate(&pm, &jr);    CHKERRQ(ierr);
	ierr = PCStokesCreate(&pc, pm); CHKERRQ(ierr);

	// create nonlinear solver
	ierr = NLSolCreate(&nl, pc, &snes); CHKERRQ(ierr);

	// create output object for all requested output variables
	ierr = PVOutCreate(&pvout, &jr, user.OutputFile); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," \n");

	//===============
	// TIME STEP LOOP
	//===============

//	PetscTime(&cputime_start_tstep);

	do
	{
		PetscPrintf(PETSC_COMM_WORLD,"Time step %lld -------------------------------------------------------- \n", (LLD)JacResGetStep(&jr));

		//==========================================================================================
		// Correct particles in case we employ an internal free surface
//		ierr = CorrectPhasesForInternalFreeSurface( &user);	CHKERRQ(ierr);

		//==========================================================================================
		// Set material properties at integration points in case we are performing a benchmark
//		if (user.AnalyticalBenchmark == PETSC_TRUE)
//		{
//			ierr = LaMEM_Initialize_AnalyticalSolution(&user, C); CHKERRQ(ierr);
//		}
		//==========================================================================================


		//=========================================================================================
		//	NONLINEAR THERMO-MECHANICAL SOLVER
		//=========================================================================================

		if(user.SkipStokesSolver != PETSC_TRUE)
		{
			PetscTime(&cputime_start_nonlinear);

			// initialize boundary constraint vectors
			ierr = BCApply(&bc, &fs); CHKERRQ(ierr);

			// compute inverse elastic viscosities
			ierr = JacResGetI2Gdt(&jr); CHKERRQ(ierr);

			// solve nonlinear system with SNES
			ierr = SNESSolve(snes, NULL, jr.gsol); CHKERRQ(ierr);

			// print analyze convergence/divergence reason & iteration count
			ierr = SNESPrintConvergedReason(snes); CHKERRQ(ierr);

			PetscTime(&cputime_end);

			PetscPrintf(PETSC_COMM_WORLD, " Nonlinear solve took %g s\n", cputime_end - cputime_start_nonlinear);
		}

		//==========================================
		// END OF NONLINEAR THERMO-MECHANICAL SOLVER
		//==========================================

		// view nonlinear residual
		ierr = JacResViewRes(&jr); CHKERRQ(ierr);


		//=========================================================================================
		// In case we perform an analytical benchmark, compare numerical and analytical results
//		if (user.AnalyticalBenchmark)
//		{

//			ierr = LaMEM_CompareNumerics_vs_AnalyticalSolution(&user, C, user.sol_advect, user.Pressure); CHKERRQ(ierr);

//		}
		//==========================================================================================

		// select new time step
		ierr = JacResGetCourantStep(&jr); CHKERRQ(ierr);

		//==========================================================================================
		// MARKER & FREE SURFACE ADVECTION + EROSION
		//==========================================================================================

		ierr = ADVAdvect(&actx); CHKERRQ(ierr);

		// advect pushing block
		ierr = BCAdvectPush(&bc, &jr.ts); CHKERRQ(ierr);


		//==========================================================================================
		// EROSION
		//==========================================================================================
//		if ( user.ErosionParameters.UseInternalFreeSurface == 1)
//		{
			// In case we have an internal free surface, advect it
//			ierr = AdvectAndUpdate_InternalFreeSurface( &user, user.sol_advect ); CHKERRQ(ierr);	 // Update internal free surface

			// Apply sedimentation to the internal free surface
//			ierr = ApplySedimentationToInternalFreeSurface(&user); // Apply sedimentation to internal free surface if requested

//			if(user.ErosionParameters.ErosionModel == 1)
//			{
//				ierr = ApplyInfinitelyFastErosionToInternalFreeSurface(&user); // Apply fast erosion to internal free surface if requested
//			}
//			else if (user.ErosionParameters.ErosionModel == 2)
//			{
				// use FD based erosion code [serial only!]
//				ierr = ApplySerialErosionCodeToInternalFreeSurface(&user); CHKERRQ(ierr);
//			}

//			ierr = CorrectPhasesForInternalFreeSurface( &user); CHKERRQ(ierr);
//		}
		//==========================================================================================


		// apply background strain-rate bc (dwindlar condition)
		ierr = FDSTAGStretch(&fs, bc.Exx, bc.Eyy, jr.ts.dt); CHKERRQ(ierr);


		//==========================================================================================


		// compute gravity misfits
//		ierr = CalculateMisfitValues(&user, C, itime, LaMEM_OutputParameters); CHKERRQ(ierr);

		//==========================================================================================
		// Save data to disk
		//==========================================================================================

// WARNING SORT OUT THIS MESS!

		if(user.save_timesteps != 0) LaMEMMod( JacResGetStep(&jr), user.save_timesteps, &SaveOrNot);
		else                         SaveOrNot = 2;

		if(SaveOrNot == 0)
		{
			char *DirectoryName = NULL;

			// create directory (encode current time & step number)
			asprintf(&DirectoryName, "Timestep_%1.6lld_%1.6e", (LLD)JacResGetStep(&jr), JacResGetTime(&jr));

			ierr = LaMEMCreateOutputDirectory(DirectoryName); CHKERRQ(ierr);

			// Matlab output (a) -> regular files
//			if (user.MatlabOutputFiles==1)
//			{
//				ierr = WriteOutputFileMatlab(&user, user.DA_Vel, user.DA_Temp, user.sol_advect, NULL, itime, C->ElementType, C->ngp_vel, DirectoryName); CHKERRQ(ierr);
//  		}

			// Matlab output (b) -> particles (mainly debugging purposes)
//			if (user.SaveParticles==1)
//			{
//				ierr = WriteParticlesToDisc(&user, itime); CHKERRQ(ierr);
//			}

			// Paraview output (a) -> topography only
//			if (user.VTKOutputFiles==1)
//			{
				// Surface, bottom topography and internal erosion surface if required
//				ierr = WriteTopographyOutputFile_VTS(&user, itime, DirectoryName); CHKERRQ(ierr);
//			}

			// Paraview output (b) -> phases only
//			if (user.AVDPhaseViewer)
//			{	// Creates a Voronoi diagram and writes the phases as VTS files
//				ierr = WritePhasesOutputFile_VTS(C, &user, itime, DirectoryName); CHKERRQ(ierr);
//			}

			// Paraview output
			ierr = PVOutWriteTimeStep(&pvout, &jr, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// clean up
			if(DirectoryName) free(DirectoryName);
		}


		//==========================================================================================
		// Perform phase transitions on particles (i.e. change the phase number of particles)
		// ParticlePhaseTransitions(&user);
		//==========================================================================================

		//==========================================================================================
		// Create breakpoint files, for restarting the code

//		if (user.save_breakpoints > 0)
//		{
		// Standard breakpoint file
//			LaMEMMod( itime, user.save_breakpoints, &SaveOrNot);
//			if ((SaveOrNot==0))
//			{
//				ierr = SaveBreakPoint(&user,user.DA_Vel, user.DA_Pres, sol, Temp, Pressure, itime, 0); CHKERRQ(ierr);
//			}

		// Incremental breakpoint file
//			LaMEMMod(itime, user.incr_breakpoints, &SaveOrNot);
//			if ((itime != 0) && (SaveOrNot==0))
//			{
//				ierr = SaveBreakPoint(&user,user.DA_Vel, user.DA_Pres, sol, Temp, Pressure, itime, user.break_point_number); CHKERRQ(ierr);
//				user.break_point_number = user.break_point_number+1;
//			}
//		}

		//==========================================================================================

		// store markers to disk
		ierr = ADVMarkSave(&actx, &user);  CHKERRQ(ierr);

		// update time state
		ierr = TSSolUpdate(&jr.ts, &jr.scal, &done); CHKERRQ(ierr);

	} while(done != PETSC_TRUE);

	//======================
	// END OF TIME STEP LOOP
	//======================

	// print total solution time
//	PetscTime(&cputime_end);
//	PetscPrintf(PETSC_COMM_WORLD,"# Total time required: %g s \n",cputime_end - cputime_start0);


//	ierr = PetscFree(user.TimeDependentData); CHKERRQ(ierr);


	// Cleanup FD erosion code
//	if ((user->ErosionParameters.ErosionModel==2) && (rank==0))
//	{
		// not good! implement erosion-code-data-structure create and destroy functions
//		ierr = DMDestroy (&user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode);CHKERRQ(ierr);
//		ierr = VecDestroy(&user->ErosionParameters.FE_ErosionCode.ErosionSurface);   CHKERRQ(ierr);
//	}

	// cleanup
	ierr = FDSTAGDestroy(&fs);   CHKERRQ(ierr);
	ierr = BCDestroy(&bc);       CHKERRQ(ierr);
	ierr = JacResDestroy(&jr);   CHKERRQ(ierr);
	ierr = ADVDestroy(&actx);    CHKERRQ(ierr);
	ierr = PCStokesDestroy(pc);  CHKERRQ(ierr);
	ierr = PMatDestroy(pm);      CHKERRQ(ierr);
	ierr = SNESDestroy(&snes);   CHKERRQ(ierr);
	ierr = NLSolDestroy(&nl);    CHKERRQ(ierr);
	ierr = PVOutDestroy(&pvout); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================
// END OF LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

