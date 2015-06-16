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
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "paraViewOutBin.h"
#include "paraViewOutSurf.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "multigrid.h"
#include "advect.h"
#include "marker.h"
#include "input.h"
#include "matProps.h"
#include "break.h"

//==========================================================================================================
// LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

#undef __FUNCT__
#define __FUNCT__ "LaMEMLib_FDSTAG"
PetscErrorCode LaMEMLib_FDSTAG(void *echange_ctx)
{
	PetscBool          done;
	UserCtx            user;
	PetscInt           SaveOrNot;
//	PetscLogDouble     cputime_start, cputime_start0, cputime_end, cputime_start_tstep, cputime_start_nonlinear;
	PetscLogDouble     cputime_start, cputime_end, cputime_start_nonlinear, cputime_end_nonlinear;

	FDSTAG   fs;     // staggered-grid layout
	FreeSurf surf;   // free-surface grid
	BCCtx    bc;     // boundary condition context
	JacRes   jr;     // Jacobian & residual context
	AdvCtx   actx;   // advection context
	PMat     pm;     // preconditioner matrix
	PCStokes pc;     // Stokes preconditioner
	SNES     snes;   // PETSc nonlinear solver
	NLSol    nl;     // nonlinear solver context
	PVOut    pvout;  // paraview output driver
	PVSurf   pvsurf; // paraview output driver

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(echange_ctx) echange_ctx = NULL;

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

	// clear objects
	ierr = FDSTAGClear  (&fs);     CHKERRQ(ierr);
	ierr = FreeSurfClear(&surf);   CHKERRQ(ierr);
	ierr = BCClear      (&bc);     CHKERRQ(ierr);
	ierr = JacResClear  (&jr);     CHKERRQ(ierr);
	ierr = ADVClear     (&actx);   CHKERRQ(ierr);
	ierr = NLSolClear   (&nl);     CHKERRQ(ierr);
	ierr = PVOutClear   (&pvout);  CHKERRQ(ierr);
	ierr = PVSurfClear  (&pvsurf); CHKERRQ(ierr);
	ierr = PetscMemzero (&user, sizeof(UserCtx)); CHKERRQ(ierr);

	// initialize variables
	ierr = FDSTAGInitCode(&jr, &user); CHKERRQ(ierr);

	// check restart
	ierr = BreakCheck(&user); CHKERRQ(ierr);

	//======================
	// SETUP DATA STRUCTURES
	//======================

	// initialize and setup scaling object, perform scaling
	ierr = JacResInitScale(&jr, &user); CHKERRQ(ierr);

	// create staggered grid object
	ierr = FDSTAGCreate(&fs, user.nnode_x, user.nnode_y, user.nnode_z); CHKERRQ(ierr);

	// generate coordinates of grid nodes/cells
	ierr = FDSTAGGenCoord(&fs, &user); CHKERRQ(ierr);

	// save processor partitioning
	if(user.SavePartitioning)
	{
		ierr = FDSTAGProcPartitioning(&fs, jr.scal.length); CHKERRQ(ierr);

		// return immediately
		ierr = FDSTAGDestroy(&fs); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// print essential grid details
	ierr = FDSTAGView(&fs); CHKERRQ(ierr);

	// create boundary condition context
	ierr = BCCreate(&bc, &fs, &jr.ts, &jr.scal); CHKERRQ(ierr);

	// set boundary conditions parameters
	ierr = BCSetParam(&bc, &user); CHKERRQ(ierr);

	// overwrite grid info if restart and background strain-rates are applied - before marker init
	// get rid of this!
	if(user.restart == 1
	&&(bc.ExxAct    == PETSC_TRUE
	|| bc.EyyAct    == PETSC_TRUE))
	{
		ierr = BreakReadGrid(&user, &fs); CHKERRQ(ierr);
	}

	// set pushing block parameters
	ierr = BCSetPush(&bc, &user); CHKERRQ(ierr);

	// set parameters from PETSc options
	ierr = BCReadFromOptions(&bc); CHKERRQ(ierr);

	// create Jacobian & residual evaluation context
	ierr = JacResCreate(&jr, &fs, &bc); CHKERRQ(ierr);

	// create free surface grid
	ierr = FreeSurfCreate(&surf, &jr); CHKERRQ(ierr);

	// initialize free surface from breakpoints if restart
	if (user.restart == 1 && surf.UseFreeSurf == PETSC_TRUE) { ierr = BreakReadSurf(&fs, &surf); CHKERRQ(ierr); }

	// create advection context
	ierr = ADVCreate(&actx, &fs, &jr); CHKERRQ(ierr);

	// initialize markers
	ierr = ADVMarkInit(&actx, &user); CHKERRQ(ierr);


// ACHTUNG!
// change marker phase when crossing free surface
ierr = ADVMarkCrossFreeSurf(&actx, &surf); CHKERRQ(ierr);


	// update phase ratios taking into account actual free surface position
	ierr = FreeSurfGetAirPhaseRatio(&surf); CHKERRQ(ierr);

	// initialize temperature
	ierr = JacResInitTemp(&jr); CHKERRQ(ierr);

	// create Stokes preconditioner & matrix
	ierr = PMatCreate(&pm, &jr);    CHKERRQ(ierr);
	ierr = PCStokesCreate(&pc, pm); CHKERRQ(ierr);

	// create nonlinear solver
	ierr = NLSolCreate(&nl, pc, &snes); CHKERRQ(ierr);

	// create output object for all requested output variables
	ierr = PVOutCreate(&pvout, &jr, user.OutputFile); CHKERRQ(ierr);

	// create output object for the free surface
	ierr = PVSurfCreate(&pvsurf, &surf, user.OutputFile); CHKERRQ(ierr);

	// read breakpoint files if restart was requested and if is possible
	if (user.restart==1) { ierr = BreakRead(&user, &actx, &nl.jtype); CHKERRQ(ierr); }

	PetscPrintf(PETSC_COMM_WORLD," \n");

	//===============
	// TIME STEP LOOP
	//===============

//	PetscTime(&cputime_start_tstep);

	do
	{
		PetscPrintf(PETSC_COMM_WORLD,"Time step %lld -------------------------------------------------------- \n", (LLD)JacResGetStep(&jr));

		//====================================
		//	NONLINEAR THERMO-MECHANICAL SOLVER
		//====================================

		// initialize boundary constraint vectors
		ierr = BCApply(&bc); CHKERRQ(ierr);



// ACHTUNG!

// initialize temperature
ierr = JacResInitTemp(&jr); CHKERRQ(ierr);

// copy to local vector, apply bc constraints
ierr = JacResCopyTemp(&jr); CHKERRQ(ierr);



		// compute inverse elastic viscosities
		ierr = JacResGetI2Gdt(&jr); CHKERRQ(ierr);

		if(user.SkipStokesSolver != PETSC_TRUE)
		{
			PetscTime(&cputime_start_nonlinear);

			// solve nonlinear system with SNES
			ierr = SNESSolve(snes, NULL, jr.gsol); CHKERRQ(ierr);

			// print analyze convergence/divergence reason & iteration count
			ierr = SNESPrintConvergedReason(snes); CHKERRQ(ierr);

			PetscTime(&cputime_end_nonlinear);

			PetscPrintf(PETSC_COMM_WORLD, " Nonlinear solve took %g s\n", cputime_end_nonlinear - cputime_start_nonlinear);
		}
		else
		{
			// just evaluate initial residual
			ierr = FormResidual(snes, jr.gsol, jr.gres, &nl); CHKERRQ(ierr);
		}

		// switch off initial guess flag
		if(!JacResGetStep(&jr))
		{
			jr.matLim.initGuessFlg = PETSC_FALSE;
		}

		// view nonlinear residual
		ierr = JacResViewRes(&jr); CHKERRQ(ierr);

		// select new time step
		ierr = JacResGetCourantStep(&jr); CHKERRQ(ierr);

		//==========================================
		// MARKER & FREE SURFACE ADVECTION + EROSION
		//==========================================

		// advect free surface
		ierr = FreeSurfAdvect(&surf); CHKERRQ(ierr);

		// advect markers
		ierr = ADVAdvect(&actx); CHKERRQ(ierr);

		// apply background strain-rate "DWINDLAR" BC (Bob Shaw "Ship of Strangers")
		ierr = BCStretchGrid(&bc); CHKERRQ(ierr);

		// exchange markers between the processors (after mesh advection)
		ierr = ADVExchange(&actx); CHKERRQ(ierr);

		// apply erosion to the free surface
		ierr = FreeSurfAppErosion(&surf); CHKERRQ(ierr);

		// apply sedimentation to the free surface
		ierr = FreeSurfAppSedimentation(&surf); CHKERRQ(ierr);

		// change marker phase when crossing flat surface or free surface with fast sedimentation/erosion
		ierr = ADVMarkCrossFreeSurf(&actx, &surf); CHKERRQ(ierr);

		// remap markers onto (stretched) grid
		ierr = ADVRemap(&actx); CHKERRQ(ierr);

		// update phase ratios taking into account actual free surface position
		ierr = FreeSurfGetAirPhaseRatio(&surf); CHKERRQ(ierr);

		// advect pushing block
		ierr = BCAdvectPush(&bc); CHKERRQ(ierr);

		//==================
		// Save data to disk
		//==================

		// compute gravity misfits
//		ierr = CalculateMisfitValues(&user, C, itime, LaMEM_OutputParameters); CHKERRQ(ierr);

// WARNING! SORT OUT THIS MESS!

		if(user.save_timesteps != 0) LaMEMMod( JacResGetStep(&jr), user.save_timesteps, &SaveOrNot);
		else                         SaveOrNot = 2;

		if(SaveOrNot == 0)
		{
			char *DirectoryName = NULL;

			// create directory (encode current time & step number)
			asprintf(&DirectoryName, "Timestep_%1.6lld_%1.6e", (LLD)JacResGetStep(&jr), JacResGetTime(&jr));

			ierr = LaMEMCreateOutputDirectory(DirectoryName); CHKERRQ(ierr);

			// Paraview output (b) -> phases only
//			if (user.AVDPhaseViewer)
//			{	// Creates a Voronoi diagram and writes the phases as VTS files
//				ierr = WritePhasesOutputFile_VTS(C, &user, itime, DirectoryName); CHKERRQ(ierr);
//			}

			// grid ParaView output
			ierr = PVOutWriteTimeStep(&pvout, &jr, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// free surface ParaView output
			ierr = PVSurfWriteTimeStep(&pvsurf, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// clean up
			if(DirectoryName) free(DirectoryName);
		}

		// store markers to disk
		ierr = ADVMarkSave(&actx, &user);  CHKERRQ(ierr);

		// update time state
		ierr = TSSolUpdate(&jr.ts, &jr.scal, &done); CHKERRQ(ierr);

		// create BREAKPOINT files, for restarting the code
		if (user.save_breakpoints > 0) LaMEMMod(JacResGetStep(&jr)-1, user.save_breakpoints, &SaveOrNot);
		else                           SaveOrNot = 2;

		if (SaveOrNot == 0) { ierr = BreakWrite(&user, &actx, &surf, nl.jtype); CHKERRQ(ierr); }

	} while(done != PETSC_TRUE);

	//======================
	// END OF TIME STEP LOOP
	//======================

	// print total solution time
//	PetscTime(&cputime_end);
//	PetscPrintf(PETSC_COMM_WORLD,"# Total time required: %g s \n",cputime_end - cputime_start0);

	// cleanup
	ierr = FDSTAGDestroy(&fs);     CHKERRQ(ierr);
	ierr = FreeSurfDestroy(&surf); CHKERRQ(ierr);
	ierr = BCDestroy(&bc);         CHKERRQ(ierr);
	ierr = JacResDestroy(&jr);     CHKERRQ(ierr);
	ierr = ADVDestroy(&actx);      CHKERRQ(ierr);
	ierr = PCStokesDestroy(pc);    CHKERRQ(ierr);
	ierr = PMatDestroy(pm);        CHKERRQ(ierr);
	ierr = SNESDestroy(&snes);     CHKERRQ(ierr);
	ierr = NLSolDestroy(&nl);      CHKERRQ(ierr);
	ierr = PVOutDestroy(&pvout);   CHKERRQ(ierr);
	ierr = PVSurfDestroy(&pvsurf); CHKERRQ(ierr);

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD, " Simulation took %g s\n", cputime_end - cputime_start);

	PetscFunctionReturn(0);
}
//==========================================================================================================
// END OF LAMEM LIBRARY MODE ROUTINE
//==========================================================================================================

//========================================================================================
// few parameters that can be reused
//========================================================================================

//	user.save_breakpoints = -1;
//	user.restart = 0;
//	user.AnalyticalBenchmark == PETSC_TRUE
//	user.InitialErosionSurfaceFromFile == 1
//	user.AVDPhaseViewer

//========================================================================================

