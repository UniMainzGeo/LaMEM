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

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Version.h"
#include "Output.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "Mesh.h"
#include "Utils.h"
#include "LaMEM_Initialize.h"
#include "LaMEM_Particles.h"
#include "Assembly_FDSTAG.h"
#include "LaMEMLib_FDSTAG_private.h"

//#include "Breakpoint.h"
//#include "LaMEM_FE_ErosionCode.h"
//#include "LaMEM_AnalyticalSolutions.h"

#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "interface.h"
#include "multigrid.h"
#include "check_fdstag.h"
#include "matrix.h"
#include "advect.h"
#include "marker.h"

//==========================================================================================================
// LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

#undef __FUNCT__
#define __FUNCT__ "LaMEMLib_FDSTAG"
PetscErrorCode LaMEMLib_FDSTAG(PetscBool InputParamFile, const char *ParamFile, PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id)
{

//	LaMEMVelPressureDA C;
	UserContext        user;
	PetscInt           SaveOrNot;
	PetscInt           itime;
//	DAVPElementType    vpt_element_type;
//	PetscLogDouble     cputime_start, cputime_start0, cputime_end, cputime_start_tstep, cputime_start_nonlinear;
	PetscLogDouble     cputime_end, cputime_start_nonlinear;

	AdvCtx         actx;   // advection context
	NLCtx          nlctx;  // nonlinear solver context
	FDSTAG         fs;     // staggered-grid layout
	JacResCtx      jrctx;  // fdstag Jacobian & residual context
	PVOut          pvout;  // fdstag paraview output driver
	BCCtx          cbc;    // boundary condition context (coupled)
	BCCtx          sbc;    // boundary condition context (split)

	SNES           snes;   // nonlinear solver
	SNESLineSearch snesls; // line search context
//	Vec            res;    // residual operator
//	Vec            sol;    // nonlinear solution vector
	BlockMat       bmat;   // block recinditioner matrix

	KSP            ksp;
	PC             pc;

//	PetscViewer  viewer;
//	PetscBool    do_restart;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(LaMEM_OutputParameters) LaMEM_OutputParameters = NULL;

//	PetscTime(&cputime_start0);

	// Start code
	PetscPrintf(PETSC_COMM_WORLD,"# -------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"#                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"#              Current Revision: %s - %s		     \n",__SVNREVISION__,__SVNDATE__);
	PetscPrintf(PETSC_COMM_WORLD,"#  Modified items: %s\n",__SVNMANCHANGES__);
	PetscPrintf(PETSC_COMM_WORLD,"#     Compiled: Date: %s - Time: %s - Mode:  %s	    \n",__DATE__,__TIME__,__OPTMODE__);
	PetscPrintf(PETSC_COMM_WORLD,"# -------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"#        STAGGERED-GRID FINITE DIFFERENCE CANONICAL IMPLEMENTATION           \n");
	PetscPrintf(PETSC_COMM_WORLD,"# -------------------------------------------------------------------------- \n");

	// Initialize context
	PetscMemzero(&user, sizeof(UserContext));

	// set input file flag & name
	user.InputParamFile = InputParamFile;
	if(InputParamFile == PETSC_TRUE) strcpy(user.ParamFile, ParamFile);

	// create element data structure (T_O___B_E___R_E_M_O_V_E_D)!
//	ierr = LaMEMVelPressureDACreate(DAVP_Q2PM1G, &C); CHKERRQ(ierr);
//	ierr = LaMEMVelPressureDAGetInfo( C, &vpt_element_type,0,0,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	__ELEMENT_TYPE__ = ELEMENT_FDSTAG;

	// initialize variables
	ierr = InitializeCode(&user); CHKERRQ(ierr);

	// Give current LaMEM session a specific group ID
	user.Optimisation.mpi_group_id = *mpi_group_id;

	// Check LaMEM input variables, and stop code if impossible combinations are added
	ierr = LaMEM_Initialize_CheckInput(); CHKERRQ(ierr);

	//========================================================================================
	// Test for unsupported options
	//========================================================================================

//	user.ParticleInput = 1;
//	user.save_breakpoints = -1;
//	user.BC.InternalBound = -1;
//	user.restart = 0;
//	user.ErosionParameters.ErosionModel = 0;
//	user.ErosionParameters.UseInternalFreeSurface = 0;
//	user.InitialMeshFromFile = 0;
//	user.SavePartitioning = 0;
//	user.InitialErosionSurfaceFromFile = 0;
//	user.LoadInitialParticlesFromDisc = 0;
//	user.remesh = 0;
	//========================================================================================
	// Setting up the solver
	//========================================================================================

//	PetscTime(&cputime_start);

	// Initialize DMDA and matrices, that:
	// (1) Store material properties,
	// (2) Distribute the structured grid on the various processors
	// (3) Are required for the velocity and pressure Stokes solution

//	ierr = LaMEM_Initialize_StokesSolver_FDSTAG(&user, vpt_element_type, C); CHKERRQ(ierr);


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

//	PetscTime(&cputime_end);
//	PetscPrintf(PETSC_COMM_WORLD,"# Finished initializing Stokes matrices : %g s \n",cputime_end - cputime_start);
//	PetscTime(&cputime_start);

	// create solution vectors
//	ierr = CreateSolutionVectors(&user); CHKERRQ(ierr);

	// Generate fine-grid mesh
//	ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);

	// Save Partitioning
//	if(user.SavePartitioning)
//	{
//		SaveProcessorPartitioning(&user);
//	}

//	if (user.BC.InternalBound > 0)
//	{
//		DefineInternalBC(&user);
//	}

	// Create a deformed mesh
//	ierr = DMDASetUniformCoordinates(user.DA_Vel,user.x_left,user.x_left+user.W,user.y_front,user.y_front + user.L,user.z_bot,user.z_bot + user.H); CHKERRQ(ierr);
//	ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);

	// Set the initial material properties
//	ierr = SetMaterialProperties(C, &user, user.DA_Vel); CHKERRQ(ierr);

	// INITIALIZE PARTICLES-based routines
//	if(user.LoadInitialParticlesFromDisc != 0)
//	{
//		// [0]. Load particles from InitialParticles directory for FDSTAG if LoadInitialParticlesFromDisc==1
//		// [1]. Load particles from Matlab generated files if LoadInitialParticlesFromDisc==2
//
//		ierr = LoadInitialParticlesFromDisc_FDSTAG(&user); CHKERRQ(ierr);
//
//		PetscPrintf(PETSC_COMM_WORLD,"# Successful loading of particles from disc for FDSTAG \n");
//
//		if(user.LoadInitialParticlesFromDisc == 2)
//		{
//			// Find the element in which particles reside
//			ierr = GetParticleNodes(C, user.DA_Vel, &user); CHKERRQ(ierr);		// Put the particle to the correct nodes
//		}
//
//	}
//	else
//	{
		// [2]. Particles are initialized internally - works with ALL OTHER cases for Model.Setup
//		ierr = InitializeParticles_ElementWise(C, user.DA_Vel, &user); 	CHKERRQ(ierr);			// Initialize tracers

		// Find the element in which particles reside
//		ierr = GetParticleNodes(C, user.DA_Vel, &user);	CHKERRQ(ierr); // Put the particle to the correct elements

		// Set tracer properties from
//		ierr = SetInitialTracerProperties(user.DA_Vel, &user); CHKERRQ(ierr); // Set initial tracer phases

		// Read the initial mesh from file if asked for
//		if (user.InitialMeshFromFile == 1)
//		{
//			ierr = ReadMeshFromFile(user.DA_Vel, &user); CHKERRQ(ierr);
//			ierr = InterpolateParticles(C,user.DA_Vel, &user); 				CHKERRQ(ierr);
//		}
//
//		if (user.LoadInitialParticlesFromDisc==1 )
//		{
//			// The 'new' way of setting particles
//			ierr = LoadInitialParticlesFromDisc( &user); CHKERRQ(ierr);
//
//			PetscPrintf(PETSC_COMM_WORLD,"# Successful loading of particles from disc \n");
//
//		}
//
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

	// Save initial mesh to file
//	if (user.save_breakpoints > 0)
//	{
//		PetscPrintf(PETSC_COMM_WORLD," Saving initial mesh to file \n");
//		ierr = SaveInitialMesh(&user,user.DA_Vel,"InitialMesh"); CHKERRQ(ierr);
//	}


	// If we are using FDSTAG, the grid MUST be regular and undeformed.
	// Because of the Finite Element manner in which we create the mesh and particles,
	// we can actually initialize the mesh and particles in an irregular manner (or even read an arbitrary mesh from file) and set the initial particle
	// distribution quite accurately in this manner. Yet, at this stage we must ensure that the grid is really undeformed (otherwise we'll have
	// convergence issues with the code that are not easy to detect).

//	user.ampl2D    = 0;
//	user.ampl3D    = 0;
//	user.amplNoise = 0;

//	ierr = GenerateMeshFine(user.DA_Vel, &user); CHKERRQ(ierr);

//	if ((user.save_breakpoints >= 0) && (user.LoadInitialParticlesFromDisc==0))
//	{
		// Save the initial particles to disk (every proc saves one file)
		// Deactivate this by adding -save_breakpoints -1  to the command line (e.g. while performing scaling tests)
//		ierr = WriteInitialParticlesToDisc( &user ); CHKERRQ(ierr);	// Write the initial particles to the ./InitialParticles directory
//	}

	// ========================================================================================================================

	// Beginning of time step
//	if (user.time_start == 0) { user.PlasticityCutoff = 1; }
//	else                      { user.PlasticityCutoff = 0; }

	if(user.fileno > 0.0)
	{
		user.dt = user.dt_max;
	}

	//===============
	// STAGGERED-GRID
	//===============

	// create staggered grid object
	ierr = FDSTAGCreate(&fs, user.nnode_x, user.nnode_y, user.nnode_z,
		user.cpu_x, user.cpu_y, user.cpu_z); CHKERRQ(ierr);

	// generate coordinates of grid nodes/cells
	ierr = FDSTAGGenCoord(&fs, &user); CHKERRQ(ierr);

	// create boundary condition context
	ierr = FDSTAGCreateBCCtx(&cbc, &fs); CHKERRQ(ierr);
	ierr = FDSTAGCreateBCCtx(&sbc, &fs); CHKERRQ(ierr);

	// create Jacobian & residual evaluation context
	ierr = FDSTAGCreateJacResCtx(&fs, &jrctx, user.num_phases, 0); CHKERRQ(ierr);

	// initialize boundary constraint vectors
	ierr = FDSTAGInitBC(&cbc, &fs, IDXCOUPLED);   CHKERRQ(ierr);
	ierr = FDSTAGInitBC(&sbc, &fs, IDXUNCOUPLED); CHKERRQ(ierr);

	// zero out global solution vector
	ierr = VecZeroEntries(jrctx.gsol); CHKERRQ(ierr);

	// WARNING! NON-DIMENSIONAL INPUT ASSUMED!
	ComputeScaling(&jrctx.scal, 1.0, 1.0, 1.0, 1.0, 1.0);

	// WARNING! NO TEMPERATURE! Set local temperature vector to unity (ad-hoc)
	ierr = VecSet(jrctx.lT, 1.0); CHKERRQ(ierr);

	// initialize material properties
	ierr = FDSTAGInitMaterialProps(&jrctx, &user); CHKERRQ(ierr);

	// initialize material parameter limits
	ierr = SetMatParLim(&jrctx.matLim, &user); CHKERRQ(ierr);

	// initialize gravity acceleration
	jrctx.grav[0] = 0.0;
	jrctx.grav[1] = 0.0;
	jrctx.grav[2] = user.Gravity;

	// create output object for all requested output variables
	ierr = PVOutCreate(&pvout, &fs, user.OutputFile); CHKERRQ(ierr);

	// create advection context
	ierr = ADVCreate(&actx); CHKERRQ(ierr);

	// initialize markers
	ierr = ADVMarkInit(&actx, &fs, &user); CHKERRQ(ierr);

	//================================
	// SETUP JACOBIAN & PRECONDITIONER
	//================================

	// create block preconditioner
	ierr = BlockMatCreate(&bmat, &fs, jrctx.gsol); CHKERRQ(ierr);

	// create residual & solution vectors
//	ierr = MatGetVecs(bmat.A, &sol, &res); CHKERRQ(ierr);

	// create nonlinear solver context
	ierr = NLCtxCreate(&nlctx, &bmat, &fs, &cbc, &sbc, &jrctx); CHKERRQ(ierr);

	//=======================
	// SETUP NONLINEAR SOLVER
	//=======================

	// create nonlinear solver
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);

	ierr = SNESSetType(snes, SNESNEWTONLS); CHKERRQ(ierr);

	ierr = SNESGetLineSearch(snes, &snesls); CHKERRQ(ierr);

	ierr = SNESLineSearchSetType(snesls, SNESLINESEARCHBASIC); CHKERRQ(ierr);

	// set initial guess (on subsequent steps the previous solution will be taken)
//	ierr = VecSet(jrctx.gsol, 0.0); CHKERRQ(ierr);

	// set residual evaluation function
	ierr = SNESSetFunction(snes, jrctx.gres, &FDSTAGFormResidual, &nlctx); CHKERRQ(ierr);

	// set Jacobian & preconditioner evaluation function
	ierr = SNESSetJacobian(snes, nlctx.Jac, NULL, &FDSTAGFormJacobian, &nlctx); CHKERRQ(ierr);
//	SNESSetPicard

	// set block stop test & residual monitor
//	ierr = SNESSetConvergenceTest(snes, SNESBlockStopTest, &nlctx, NULL); CHKERRQ(ierr);

	// setup snes options from command line
	ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

	//=====================================
	// SETUP LINEAR SOLVER & PRECONDITIONER
	//=====================================

	// retrieve linear solver
//	ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);

	// set initial guess flag
//	ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);

	// set block stop test & residual monitor
//	ierr = KSPSetConvergenceTest(ksp, &KSPBlockStopTest, &bmat, NULL);

	// set additional options
//	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

	// a bit more control on our side, we will be responsible for preconditioning
//	ierr = KSPGetPC(ksp, &pc);                    CHKERRQ(ierr);
//	ierr = PCSetType(pc, PCSHELL);                CHKERRQ(ierr);
//	ierr = PCShellSetContext(pc, &bmat);          CHKERRQ(ierr);
//	ierr = PCShellSetApply(pc, &ApplyFieldSplit); CHKERRQ(ierr);

	//===============
	// TIME STEP LOOP
	//===============

//	PetscTime(&cputime_start_tstep);

	for(itime = user.time_start; itime < user.time_end; itime++)
	{


		//==========================================================================================
		// Correct particles in case we employ an internal free surface
//		ierr = CorrectPhasesForInternalFreeSurface( &user);	CHKERRQ(ierr);
		// compute history-dependent material properties from particles (using distance-based averaging)
		// THIS IS WHERE PHASE RATIOS ARE INITIALIZED
//		ierr = ComputeHistoryPropertiesFromParticles_FDSTAG(C, &user); CHKERRQ(ierr);
		// Update material properties on FDSTAG grid
//		ierr = UpdateMaterialProperties_FDSTAG(&user); CHKERRQ(ierr);
		//==========================================================================================
		// Set material properties at integration points in case we are performing a benchmark
//		if (user.AnalyticalBenchmark == PETSC_TRUE)
//		{
//			ierr = LaMEM_Initialize_AnalyticalSolution(&user, C); CHKERRQ(ierr);
//		}
		//==========================================================================================

		// set time step
		jrctx.dt = user.dt;

		// copy phase ratios
//		ierr = FDSTAGInitPhaseRatios (&fs, &jrctx, &user); CHKERRQ(ierr);

		// project properties from markers to grid
		ierr = ADVProjHistMarkGrid(&actx, &fs, &jrctx); CHKERRQ(ierr);

		// compute inverse elastic viscosities
		ierr = FDSTAGetI2Gdt(&fs, &jrctx); CHKERRQ(ierr);

		//=========================================================================================
		//	NONLINEAR THERMO-MECHANICAL SOLVER
		//=========================================================================================

// ADHOC (uncomment test)
//		ierr = StrainRateInterpTest(&fs, &jrctx, &user, &pvout); CHKERRQ(ierr);
//		ierr = DoMGTests(&nlctx, &pvout); CHKERRQ(ierr);
//		ierr = DoDarcyTests(&nlctx, &user); CHKERRQ(ierr);

		//=========================================================================================

// ADHOC (delete residual evaluation here)
//		ierr = VecSet(jrctx.gsol, 0.0); CHKERRQ(ierr);
//		ierr = VecSet(jrctx.gres, 0.0); CHKERRQ(ierr);
//		ierr = FDSTAGFormResidual(NULL, jrctx.gsol, jrctx.gres, &nlctx); CHKERRQ(ierr);

// ADHOC (don't skip snes solve here)
		user.SkipStokesSolver = PETSC_TRUE;

		if(user.SkipStokesSolver != PETSC_TRUE)
		{
			PetscTime(&cputime_start_nonlinear);

			// solve nonlinear system with SNES
			ierr = SNESSolve(snes, NULL, jrctx.gsol); CHKERRQ(ierr);

			// print analyze convergence/divergence reason
//			ierr = SNESPrintConvergedReason(snes); CHKERRQ(ierr);

			// SNESGetIterationNumber(snes,&its);
			// SNESGetConvergedReason(snes,&reason);
			// PetscPrintf(PETSC_COMM_WORLD,"%s Number of nonlinear iterations = %D\n",SNESConvergedReasons[reason],its);

			PetscTime(&cputime_end);
			PetscPrintf(PETSC_COMM_WORLD,"#  Nonlinear solve took %g s \n", cputime_end - cputime_start_nonlinear);
		}
// ADHOC

		else
		{
			// assemble matrix & rhs
			ierr = VecSet(jrctx.gsol, 0.0); CHKERRQ(ierr);
			ierr = VecSet(jrctx.gres, 0.0); CHKERRQ(ierr);

			ierr = FDSTAGFormResidual(NULL, jrctx.gsol, jrctx.gres, &nlctx); CHKERRQ(ierr);

			ierr = VecScale(jrctx.gres, -1.0); CHKERRQ(ierr);

			ierr = BlockMatCompute(&bmat, &fs, &sbc, &jrctx); CHKERRQ(ierr);

			ierr = PowellHestenes(&bmat, jrctx.gres, jrctx.gsol); CHKERRQ(ierr);

			ierr = FDSTAGFormResidual(NULL, jrctx.gsol, jrctx.gres, &nlctx); CHKERRQ(ierr);

		}


		//==========================================
		// END OF NONLINEAR THERMO-MECHANICAL SOLVER
		//==========================================

		// access solution
//		ierr = BlockMatBlockToMonolithic(&bmat, sol); CHKERRQ(ierr);
//		ierr = VecCopy(bmat.wv, user.sol);            CHKERRQ(ierr);
//		ierr = VecCopy(bmat.wp, user.Pressure);       CHKERRQ(ierr);

		// In FDSTAG 'sol' gets modified in order to make use of FEM-routines. Since we are going to use sol as initial guess
		// in the next timestep or when restarting a simulation we need a unmodified version of sol.
		// This is why we make a copy of sol that will be modified (sol_advect). Later on, we save sol to Breakpointfile NOT sol_advect.
		// (tobib 20.09.12)

		// NOTE:sol_advect is initialized earlier in the code, as Sarah added an option to do a few thermal steps BEFORE the first stokes step. BK, 28.9.2012
//		ierr = VecCopy(user.sol, user.sol_advect);	CHKERRQ(ierr);

		// If FDSTAG: interpolate velocities from cell-center to corner points, such that we can use
		// the Q1P0 routines to compute particle properties and for advection.
//		FDSTAG_InterpolateVelocityPressureToCornerpoints(&user, user.sol_advect);

		//==========================================================================================
		// Compute properties @ integration points, such as strain-rate, pressure, stress etc.
		// Here, we only compute arrays that are necessary for the nonlinear iterations (for speed)
		// After the iterations are finished, we update all variables (such as finite strain etc.)
		//==========================================================================================

//		ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, Pressure, user.DA_Temp, Temp, &user, user.dt, 0 ); CHKERRQ(ierr);

		//=========================================================================================
		// In case we perform an analytical benchmark, compare numerical and analytical results
//		if (user.AnalyticalBenchmark)
//		{
			// Update properties @ integration points (since for speed reasons we do not compute
			// all strain-rate and stress components in the nonlinear iterations above).
//			ierr = IntPointProperties( C, user.DA_Vel, user.sol_advect, user.DA_Pres, user.Pressure, PETSC_NULL, PETSC_NULL, &user, user.dt, 1 ); CHKERRQ(ierr);

//			ierr = LaMEM_CompareNumerics_vs_AnalyticalSolution(&user, C, user.sol_advect, user.Pressure); CHKERRQ(ierr);

			// update INTP props once more (to correctly set P values)
//			ierr = IntPointProperties( C, user.DA_Vel, sol_advect, user.DA_Pres, user.Pressure, user.DA_Temp, Temp, &user, user.dt, 1 ); CHKERRQ(ierr);
//		}
		//==========================================================================================


		//==========================================================================================
		// Update all properties @ integration points and compute properties @ particles
		//==========================================================================================

//		PetscTime(&cputime_start);

//		ierr = ComputePropertiesAtParticles(C, user.DA_Vel, user.sol_advect, user.DA_Pres, user.Pressure, PETSC_NULL, PETSC_NULL, PETSC_NULL, &user, user.dt); CHKERRQ(ierr);

//		PetscTime(&cputime_end);

//		PetscPrintf(PETSC_COMM_WORLD,"#  Updating properties @ integration points and particles took %g s \n", cputime_end - cputime_start);

		//==========================================================================================

//		ierr = CheckVelocityError(&user);

		//==========================================================================================
		// Advect tracers
		//==========================================================================================

//		PetscPrintf(PETSC_COMM_WORLD,"  Starting particle advection \n");
//		MPI_Barrier(PETSC_COMM_WORLD);

//		ierr = AdvectParticles(C, user.DA_Vel, &user, user.sol_advect, user.dt); CHKERRQ(ierr);

//		PetscPrintf(PETSC_COMM_WORLD,"  Finished particle advection \n");
//		MPI_Barrier(PETSC_COMM_WORLD);

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

		//==========================================================================================

		// Advect grid with Background deformation rate in case we use an FDSTAG approach with Exx=Eyy!=0
//		ierr = DeformFDSTAGMeshWithBackgroundStrainrate(&user);	CHKERRQ(ierr);

		//==========================================================================================

		// Update time state
		user.time = user.time + user.dt;

		// Compute and store output properties such as velocity root mean square
//		ierr = ComputeGlobalProperties(user.DA_Vel, &user, itime, user.sol_advect, C); CHKERRQ(ierr);

		// compute gravity misfits
//		ierr = CalculateMisfitValues(&user, C, itime, LaMEM_OutputParameters); CHKERRQ(ierr);

		//==========================================================================================
		// Save data to disk
		//==========================================================================================

		if(user.save_timesteps != 0) LaMEMMod( itime, user.save_timesteps, &SaveOrNot);
		else                         SaveOrNot = 2;

// ADHOC (comment entire output section)

		if(SaveOrNot == 0)
		{
			char *DirectoryName = NULL;

			// create directory
			asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
			ierr = LaMEM_CreateOutputDirectory(DirectoryName); CHKERRQ(ierr);

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

			// Paraview output FDSTAG fields
			ierr = PVOutWriteTimeStep(&pvout, &jrctx, 0.0, itime); CHKERRQ(ierr);

			// clean up
			if(DirectoryName) free(DirectoryName);
		}

		//==========================================================================================
		// compute new time step length
		//==========================================================================================

//		ierr = CalculateTimeStep(&user, itime); CHKERRQ(ierr);

		//==========================================================================================
		// Map tracers to new grid
		//==========================================================================================

//		if ( ((user.GridAdvectionMethod != 1) ) || ((user.remesh==1) ) )
//		{
//			PetscPrintf(PETSC_COMM_WORLD," Starting GetParticleNodes \n");
//			ierr = GetParticleNodes(C, user.DA_Vel, &user); CHKERRQ(ierr);
//			PetscPrintf(PETSC_COMM_WORLD," Finished GetParticleNodes \n");
//		}
		//==========================================================================================

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

		// print time step duration
//		PetscTime(&cputime_end);
//		PetscPrintf(PETSC_COMM_WORLD,"# Finished timestep %lld out of %lld in %g s\n\n",(LLD)itime, (LLD)(user.time_end), cputime_end - cputime_start_tstep);


		// store marker to disk
		ierr = ADVMarkSave(&actx, &fs, &user);  CHKERRQ(ierr);



	}

	//======================
	// END OF TIME STEP LOOP
	//======================

	// print total solution time
//	PetscTime(&cputime_end);
//	PetscPrintf(PETSC_COMM_WORLD,"# Total time required: %g s \n",cputime_end - cputime_start0);

	// free memory
//	ierr = DestroySolutionObjects(&user, &C); CHKERRQ(ierr);

	ierr = PetscFree(user.TimeDependentData); CHKERRQ(ierr);

	if(user.InputParamFile)
	{
		MaterialDestroy(user.PhaseMaterialProperties);
	}


	ierr = FDSTAGDestroy(&fs);                CHKERRQ(ierr);
	ierr = FDSTAGDestroyBCCtx(&cbc);          CHKERRQ(ierr);
	ierr = FDSTAGDestroyBCCtx(&sbc);          CHKERRQ(ierr);
	ierr = FDSTAGDestroyJacResCtx(&jrctx);    CHKERRQ(ierr);
	ierr = PVOutDestroy(&pvout);              CHKERRQ(ierr);
	ierr = BlockMatDestroy(&bmat);            CHKERRQ(ierr);
	ierr = NLCtxDestroy(&nlctx);              CHKERRQ(ierr);
	ierr = SNESDestroy(&snes);                CHKERRQ(ierr);
	ierr = ADVDestroy(&actx);                 CHKERRQ(ierr);

//	ierr = VecDestroy(&sol);                  CHKERRQ(ierr);
//	ierr = VecDestroy(&res);                  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================
// END OF LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

/*
//	couple of abandoned tests
				if (user.VelocityTest == PETSC_TRUE)
				{
					PetscPrintf(PETSC_COMM_WORLD,"#  ==================================================\n");
					PetscPrintf(PETSC_COMM_WORLD,"#  *** Testing performance of the Velocity Solver ***\n");
					PetscPrintf(PETSC_COMM_WORLD,"#  ==================================================\n");

					// test performance of a single velocity solve
					ierr = VelSolverTest(user.VV_MAT, user.sol, user.rhs, &user); CHKERRQ(ierr);

					// finalize computation
					goto cleanup;
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
					cleanup

				}


*/

//==========================================================================================================

