/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   LaMEMLib.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
// LAMEM LIBRARY MODE ROUTINE
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tools.h"
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
#include "paraViewOutMark.h"
#include "input.h"
#include "matProps.h"
#include "objFunct.h"
#include "AVDView.h"
#include "break.h"
#include "parsing.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "LaMEMLib"
PetscErrorCode LaMEMLib(ModParam *IOparam)
{
	PetscBool          done;
	UserCtx            user;
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
	PVSurf   pvsurf; // paraview output driver for surface
	PVMark   pvmark; // paraview output driver for markers
	PVAVD    pvavd;  // paraview output driver for AVD
	ObjFunct objf;   // objective function

	PetscErrorCode ierr;
	PetscFunctionBegin;


	//=========================================================================

	PetscBool InputParamFile;

	PetscInt found_data;

	char ParamFile[MAX_PATH_LEN];

	// check whether input file is specified
	ierr = PetscOptionsGetString(PETSC_NULL, "-ParamFile", ParamFile, MAX_PATH_LEN, &InputParamFile); CHKERRQ(ierr);

	// read additional PETSc options from input file
	if(InputParamFile == PETSC_TRUE)
	{
		ierr = PetscOptionsReadFromFile(ParamFile, &found_data, 1); CHKERRQ(ierr);
	}

	//=========================================================================


	PetscTime(&cputime_start);

	// Start code
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"     Compiled: Date: %s - Time: %s 	    \n",__DATE__,__TIME__ );
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
	ierr = FDSTAGInitCode(&jr, &user, IOparam); CHKERRQ(ierr);

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

	// set pushing block parameters
	ierr = BCSetPush(&bc, &user); CHKERRQ(ierr);

	// set parameters from PETSc options
	ierr = BCReadFromOptions(&bc); CHKERRQ(ierr);

	// overwrite grid info if restart and background strain-rates are applied - before marker init
	// get rid of this!
	if(user.restart == 1
	&&(bc.ExxAct    == PETSC_TRUE
	|| bc.EyyAct    == PETSC_TRUE))
	{
		ierr = BreakReadGrid(&user, &fs); CHKERRQ(ierr);
	}

	// create Jacobian & residual evaluation context
	ierr = JacResCreate(&jr, &fs, &bc); CHKERRQ(ierr);

	// create free surface grid
	ierr = FreeSurfCreate(&surf, &jr, &user); CHKERRQ(ierr);

	// initialize free surface from breakpoints if restart
	if (user.restart == 1 && surf.UseFreeSurf == PETSC_TRUE) { ierr = BreakReadSurf(&fs, &surf); CHKERRQ(ierr); }

	ierr = BCSetupBoundVel(&bc, surf.InitLevel); CHKERRQ(ierr);

	// create advection context
	ierr = ADVCreate(&actx, &fs, &jr); CHKERRQ(ierr);

	// initialize markers
	ierr = ADVMarkInit(&actx, &user); CHKERRQ(ierr);

	// change marker phase when crossing free surface
	ierr = ADVMarkCrossFreeSurf(&actx, &surf, 0.05); CHKERRQ(ierr);

	// set air phase to properly treat marker advection & temperature diffusion
	if(surf.UseFreeSurf == PETSC_TRUE && jr.actTemp == PETSC_TRUE)
	{
		actx.AirPhase = surf.AirPhase;
		jr.AirPhase   = surf.AirPhase;
		actx.Ttop     = bc.Ttop;
	}

	// check thermal material parameters
	ierr = JacResCheckTempParam(&jr); CHKERRQ(ierr);

	// update phase ratios taking into account actual free surface position
	ierr = FreeSurfGetAirPhaseRatio(&surf); CHKERRQ(ierr);

	// create Stokes preconditioner & matrix
	ierr = PMatCreate(&pm, &jr);    CHKERRQ(ierr);
	ierr = PCStokesCreate(&pc, pm); CHKERRQ(ierr);

	// create nonlinear solver
	ierr = NLSolCreate(&nl, pc, &snes); CHKERRQ(ierr);

	// create output object for all requested output variables
	ierr = PVOutCreate(&pvout, &jr, user.OutputFile); CHKERRQ(ierr);

	// create output object for the free surface
	ierr = PVSurfCreate(&pvsurf, &surf, user.OutputFile); CHKERRQ(ierr);

	// create output object for the markers - for debugging
	ierr = PVMarkCreate(&pvmark, &actx, user.OutputFile); CHKERRQ(ierr);

	// AVD output driver
	ierr = PVAVDCreate(&pvavd, &actx, user.OutputFile); CHKERRQ(ierr);

	// read breakpoint files if restart was requested and if is possible
	if (user.restart==1) { ierr = BreakRead(&user, &actx, &pvout, &pvsurf, &pvmark, &pvavd, &nl.jtype); CHKERRQ(ierr); }

	// finish simulation if time_end was reached
	if (jr.ts.istep >= jr.ts.nstep)
	{
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
		ierr = PVMarkDestroy(&pvmark); CHKERRQ(ierr);
		ierr = PVAVDDestroy(&pvavd);   CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	//===================
	// OBJECTIVE FUNCTION
	//===================

	// create objective function object
	ierr = ObjFunctCreate(&objf, &surf); CHKERRQ(ierr);

	// transfer misfit value to IO structure
	IOparam->mfit = objf.errtot;

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

		// initialize temperature
		ierr = JacResInitTemp(&jr); CHKERRQ(ierr);

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

			PetscPrintf(PETSC_COMM_WORLD, " Nonlinear solve took %g (sec)\n", cputime_end_nonlinear - cputime_start_nonlinear);
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

		// prescribe velocity if rotation benchmark
		if (user.msetup == ROTATION) {ierr = JacResSetVelRotation(&jr); CHKERRQ(ierr);}

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

		// remap markers onto (stretched) grid
		ierr = ADVRemap(&actx, &surf); CHKERRQ(ierr);

		// update phase ratios taking into account actual free surface position
		// -- This routine requires a modification to also correct phase ratio's at edges and not just at corners --
		// it has been deactivated temporarily (affects convergence for salt-tectonics setups with brittle overburden)
		//ierr = FreeSurfGetAirPhaseRatio(&surf); CHKERRQ(ierr);

		// advect pushing block
		ierr = BCAdvectPush(&bc); CHKERRQ(ierr);

		//==================
		// Save data to disk
		//==================

		// compute gravity misfits
//		ierr = CalculateMisfitValues(&user, C, itime, LaMEM_OutputParameters); CHKERRQ(ierr);

		if(!(JacResGetStep(&jr) % user.save_timesteps))
		{
			char *DirectoryName = NULL;

			// redefine filename in case of inversion setup
			if (IOparam->use == 1)
			{
				asprintf(&DirectoryName, "Timestep_%1.6lld", (LLD)IOparam->mID);
			}
			else
			{
				// create directory (encode current time & step number)
				asprintf(&DirectoryName, "Timestep_%1.6lld_%1.6e", (LLD)JacResGetStep(&jr), JacResGetTime(&jr));
			}

			ierr = LaMEMCreateOutputDirectory(DirectoryName); CHKERRQ(ierr);

			// AVD phase output
			ierr = PVAVDWriteTimeStep(&pvavd, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// grid ParaView output
			ierr = PVOutWriteTimeStep(&pvout, &jr, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// free surface ParaView output
			ierr = PVSurfWriteTimeStep(&pvsurf, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// marker ParaView output
			ierr = PVMarkWriteTimeStep(&pvmark, DirectoryName, JacResGetTime(&jr), JacResGetStep(&jr)); CHKERRQ(ierr);

			// clean up
			if(DirectoryName) free(DirectoryName);
		}

		// store markers to disk
		ierr = ADVMarkSave(&actx, &user);  CHKERRQ(ierr);

		// update time state
		ierr = TSSolUpdate(&jr.ts, &jr.scal, &done); CHKERRQ(ierr);

		// create BREAKPOINT files, for restarting the code
		if(user.save_breakpoints>0 && !((JacResGetStep(&jr)-1) % user.save_breakpoints))
		{
			ierr = BreakWrite(&user, &actx, &surf, &pvout, &pvsurf, &pvmark, &pvavd, nl.jtype); CHKERRQ(ierr);
		}

		// check marker phases
		ierr = ADVCheckMarkPhases(&actx, jr.numPhases); CHKERRQ(ierr);

	} while(done != PETSC_TRUE);

	//======================
	// END OF TIME STEP LOOP
	//======================

	// print total solution time
//	PetscTime(&cputime_end);
//	PetscPrintf(PETSC_COMM_WORLD,"# Total time required: %g s \n",cputime_end - cputime_start0);

	// cleanup
	ierr = ObjFunctDestroy(&objf); CHKERRQ(ierr);
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
	ierr = PVMarkDestroy(&pvmark); CHKERRQ(ierr);
	ierr = PVAVDDestroy(&pvavd);   CHKERRQ(ierr);


	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD, " Simulation took %g (sec) \n", cputime_end - cputime_start);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// END OF LAMEM LIBRARY MODE ROUTINE
//---------------------------------------------------------------------------
