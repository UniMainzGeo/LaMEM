/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
// LAMEM LIBRARY MODE ROUTINE
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "phase.h"
#include "dike.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "tools.h"
#include "fdstag.h"
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
#include "Tensor.h"
#include "advect.h"
#include "marker.h"
#include "paraViewOutMark.h"
#include "paraViewOutAVD.h"
#include "objFunct.h"
#include "adjoint.h"
#include "paraViewOutPassiveTracers.h"
#include "LaMEMLib.h"
#include "phase_transition.h"
#include "passive_tracer.h"

// surface process
#include "fastscape.h"

//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibMain(void *param,PetscLogStage stages[4])
{
	LaMEMLib       lm;
	RunMode        mode;
	PetscBool      found;
	PetscInt       exists;
	char           str[_str_len_];
	PetscLogDouble cputime_start, cputime_end;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;       

	// start code
	ierr = PetscTime(&cputime_start); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"     Compiled: Date: %s - Time: %s 	    \n",__DATE__,__TIME__ );
	PetscPrintf(PETSC_COMM_WORLD,"     Version : 2.1.4 \n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"        STAGGERED-GRID FINITE DIFFERENCE CANONICAL IMPLEMENTATION           \n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");

	// read run mode
	mode = _NORMAL_;

	ierr = PetscOptionsGetCheckString("-mode", str, &found); CHKERRQ(ierr);

	if(found)
	{
		if     (!strcmp(str, "normal"))    mode = _NORMAL_;
		else if(!strcmp(str, "restart"))   mode = _RESTART_;
		else if(!strcmp(str, "dry_run"))   mode = _DRY_RUN_;
		else if(!strcmp(str, "save_grid")) mode = _SAVE_GRID_;
		else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect run mode type: %s", str);
	}

	// cancel restart if no database is available
	if(mode == _RESTART_)
	{
		ierr = DirCheck("./restart", &exists); CHKERRQ(ierr);

		if(!exists)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No restart database available (check -mode option)");
		}
	}

	//===========
	// INITIALIZE
	//===========

	// clear
	ierr = PetscMemzero(&lm, sizeof(LaMEMLib)); CHKERRQ(ierr);

	// setup cross-references between library objects
	ierr = LaMEMLibSetLinks(&lm); CHKERRQ(ierr);

	if(mode == _SAVE_GRID_)
	{
		// save grid & exit
		ierr = LaMEMLibSaveGrid(&lm); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}
	if(mode == _NORMAL_ || mode == _DRY_RUN_)
	{
		// create library objects
		ierr = LaMEMLibCreate(&lm, param); CHKERRQ(ierr);
	}
	else if(mode == _RESTART_)
	{
		// open restart database
		ierr = LaMEMLibLoadRestart(&lm); CHKERRQ(ierr);
	}

	//======
	// SOLVE
	//======

	if(mode == _DRY_RUN_)
	{
		// compute initial residual, output & stop
		ierr = LaMEMLibDryRun(&lm); CHKERRQ(ierr);
	}
	else if(mode == _NORMAL_ || mode == _RESTART_)
	{
		// solve coupled nonlinear equations
		ierr = LaMEMLibSolve(&lm, param,stages); CHKERRQ(ierr);
	}

	// destroy library objects
	ierr = LaMEMLibDestroy(&lm); CHKERRQ(ierr);

	PetscTime(&cputime_end);

	PetscPrintf(PETSC_COMM_WORLD, "Total solution time : %g (sec) \n", cputime_end - cputime_start);
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibCreate(LaMEMLib *lm, void *param )
{
	FB *fb;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(param) param = NULL;

	// load input file
	ierr = FBLoad(&fb, PETSC_TRUE); CHKERRQ(ierr);

	// create scaling object
	ierr = ScalingCreate(&lm->scal, fb, PETSC_TRUE);CHKERRQ(ierr);

	// create time stepping object
	ierr = TSSolCreate(&lm->ts, fb); 				CHKERRQ(ierr);

	// create parallel grid
	ierr = FDSTAGCreate(&lm->fs, fb); 				CHKERRQ(ierr);

	// create material database
	ierr = DBMatCreate(&lm->dbm, fb, PETSC_TRUE); 	CHKERRQ(ierr);

	// create free surface grid
	ierr = FreeSurfCreate(&lm->surf, fb); 			CHKERRQ(ierr);

	// create boundary condition context
	ierr = BCCreate(&lm->bc, fb); 					CHKERRQ(ierr);

	// create residual & Jacobian evaluation context
	ierr = JacResCreate(&lm->jr, fb); 				CHKERRQ(ierr);

	// create dike database
	ierr = DBDikeCreate(&lm->dbdike, &lm->dbm, fb, &lm->jr, PETSC_TRUE);   CHKERRQ(ierr);

	// initialize arrays for dynamic phase transition
	ierr = DynamicPhTr_Init(&lm->jr);			CHKERRQ(ierr);

	// create advection context
	ierr = ADVCreate(&lm->actx, fb); 				CHKERRQ(ierr);

	// create passive tracers
	ierr = ADVPtrPassive_Tracer_create(&lm->actx,fb);			CHKERRQ(ierr);

	// create output object for all requested output variables
	ierr = PVOutCreate(&lm->pvout, fb); 			CHKERRQ(ierr);

	// create output object for the free surface
	ierr = PVSurfCreate(&lm->pvsurf, fb); 			CHKERRQ(ierr);

	// create output object for the markers - for debugging
	ierr = PVMarkCreate(&lm->pvmark, fb); 			CHKERRQ(ierr);

	// create output object for the passive tracers
	ierr = PVPtrCreate(&lm->pvptr, fb);              CHKERRQ(ierr);

	// AVD output driver
	ierr = PVAVDCreate(&lm->pvavd, fb); 			CHKERRQ(ierr);

	// destroy file buffer
	ierr = FBDestroy(&fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSaveGrid(LaMEMLib *lm)
{
	FB *fb;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// load input file
	ierr = FBLoad(&fb, PETSC_TRUE); CHKERRQ(ierr);

	// create scaling object
	ierr = ScalingCreate(&lm->scal, fb, PETSC_TRUE); CHKERRQ(ierr);

	// create parallel grid
	ierr = FDSTAGCreate(&lm->fs, fb); CHKERRQ(ierr);

	// save processor partitioning
	ierr = FDSTAGSaveGrid(&lm->fs); CHKERRQ(ierr);

	// destroy parallel grid
	ierr = FDSTAGDestroy(&lm->fs); CHKERRQ(ierr);

	// destroy file buffer
	ierr = FBDestroy(&fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibLoadRestart(LaMEMLib *lm)
{
	FB              *fb;
	FILE            *fp;
	PetscLogDouble  t;
	PetscMPIInt     rank;
	PetscBool       found;
	char            restartFileName[_str_len_], *fileName;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PrintStart(&t, "Loading restart database", NULL);

	// get MPI processor rank
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// compile restart file name
	asprintf(&fileName, "./restart/rdb.%1.8lld.dat", (LLD)rank);

	// open restart file for reading in binary mode
	fp = fopen(fileName, "rb");

	if(fp == NULL)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open restart file %s\n", fileName);
	}

	// read LaMEM library database
	fread(lm, sizeof(LaMEMLib), 1, fp);

	// setup cross-references between library objects
	ierr = LaMEMLibSetLinks(lm); CHKERRQ(ierr);

	// staggered grid
	ierr = FDSTAGReadRestart(&lm->fs, fp); CHKERRQ(ierr);

	// free surface
	ierr = FreeSurfReadRestart(&lm->surf, fp); CHKERRQ(ierr);

	// boundary conditions context
	ierr = BCReadRestart(&lm->bc, fp); CHKERRQ(ierr);

	// solution variables
	ierr = JacResReadRestart(&lm->jr, fp); CHKERRQ(ierr);

	// markers
	ierr = ADVReadRestart(&lm->actx, fp); CHKERRQ(ierr);

	// passive tracers read restart
	ierr = ReadPassive_Tracers(&lm->actx,fp); CHKERRQ(ierr);

	// main output driver
	ierr = PVOutCreateData(&lm->pvout); CHKERRQ(ierr);

	// surface output driver
	ierr = PVSurfCreateData(&lm->pvsurf); CHKERRQ(ierr);

	// arrays for dynamic NotInAir phase_trans
	ierr = DynamicPhTr_ReadRestart(&lm->jr, fp); CHKERRQ(ierr);

	// read from input file, create arrays for dynamic diking, and read from restart file
	ierr = DynamicDike_ReadRestart(&lm->dbdike, &lm->dbm, &lm->jr, &lm->ts, fp);  CHKERRQ(ierr);
 
	// close temporary restart file
	fclose(fp);

	// free space
	free(fileName);

	// check whether restart input file is specified
	ierr = PetscOptionsGetCheckString("-RestartParamFile", restartFileName, &found); CHKERRQ(ierr);

	if(found == PETSC_TRUE)
	{
		// load restart input file
		ierr = FBLoad(&fb, PETSC_TRUE, restartFileName); CHKERRQ(ierr);

		// override material database
		ierr = DBMatCreate(&lm->dbm, fb, PETSC_TRUE); 	CHKERRQ(ierr);

		// destroy file buffer
		ierr = FBDestroy(&fb); CHKERRQ(ierr);
	}

	PrintDone(t);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSaveRestart(LaMEMLib *lm)
{
	// save new restart database, then delete the original

	FILE           *fp;
	PetscMPIInt    rank;
	char           *fileNameTmp;
	PetscLogDouble t;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(!TSSolIsRestart(&lm->ts)) PetscFunctionReturn(0);

	PrintStart(&t, "Saving restart database", NULL);

	// get MPI processor rank
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// compile actual & temporary restart file name
	asprintf(&fileNameTmp, "./restart-tmp/rdb.%1.8lld.dat", (LLD)rank);

	// create temporary restart directory
	ierr = DirMake("./restart-tmp"); CHKERRQ(ierr);

	// open temporary restart file for writing in binary mode
	fp = fopen(fileNameTmp, "wb");

	if(fp == NULL)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open restart file %s\n", fileNameTmp);
	}

	// write LaMEM library database
	fwrite(lm, sizeof(LaMEMLib), 1, fp);

	// staggered grid
	ierr = FDSTAGWriteRestart(&lm->fs, fp); CHKERRQ(ierr);

	// free surface
	ierr = FreeSurfWriteRestart(&lm->surf, fp); CHKERRQ(ierr);

	// boundary conditions context
	ierr = BCWriteRestart(&lm->bc, fp); CHKERRQ(ierr);

	// solution variables
	ierr = JacResWriteRestart(&lm->jr, fp); CHKERRQ(ierr);

	// markers
	ierr = ADVWriteRestart(&lm->actx, fp); CHKERRQ(ierr);

	// passive tracers
	ierr = Passive_Tracer_WriteRestart(&lm->actx, fp); CHKERRQ(ierr);

	// dynamic phase transition 
	ierr = DynamicPhTr_WriteRestart(&lm->jr, fp); CHKERRQ(ierr);

	// dynamic dike 
	ierr = DynamicDike_WriteRestart(&lm->jr, fp); CHKERRQ(ierr);

	// close temporary restart file
	fclose(fp);

	// delete existing restart database
	ierr = LaMEMLibDeleteRestart(); CHKERRQ(ierr);

	// push temporary database to actual
	ierr = DirRename("./restart-tmp", "./restart");

	// free space
	free(fileNameTmp);

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDeleteRestart()
{
	// delete existing restart database
	PetscMPIInt  rank;
	int          status;
	PetscInt     exists;
	char        *fileName;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get MPI processor rank
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	asprintf(&fileName, "./restart/rdb.%1.8lld.dat", (LLD)rank);

	// check for existing restart database
	ierr = DirCheck("./restart", &exists); CHKERRQ(ierr);

	if(exists)
	{
		// delete existing database
		status = remove(fileName);

		if(status && errno != ENOENT)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to delete file %s", fileName);
		}

		ierr = DirRemove("./restart"); CHKERRQ(ierr);
	}

	// free space
	free(fileName);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDestroy(LaMEMLib *lm)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = FDSTAGDestroy  (&lm->fs);     CHKERRQ(ierr);
	ierr = FreeSurfDestroy(&lm->surf);   CHKERRQ(ierr);
	ierr = BCDestroy      (&lm->bc);     CHKERRQ(ierr);
	ierr = JacResDestroy  (&lm->jr);     CHKERRQ(ierr);
	ierr = ADVPtrDestroy  (&lm->actx);   CHKERRQ(ierr);
	ierr = ADVDestroy     (&lm->actx);   CHKERRQ(ierr);
	ierr = PVOutDestroy   (&lm->pvout);  CHKERRQ(ierr);
	ierr = PVSurfDestroy  (&lm->pvsurf); CHKERRQ(ierr);

	ierr = DynamicPhTrDestroy (&lm->dbm); CHKERRQ(ierr);
	ierr = DynamicDike_Destroy(&lm->jr); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSetLinks(LaMEMLib *lm)
{
	//======================================================================
	// LaMEM library object initialization sequence
	//
	//                         Scaling
	//                            |
	//                          TSSol
	//                            |
	//                          DBMat/DBPropDike
	//                            |
	//                         FDSTAG
	//                            |
	//                         FreeSurf
	//                            |
	//                          BCCtx
	//                            |
	//                          JacRes
	//                            |
	//                          AdvCtx
	//                            |
	//              -----------------------------
	//              |       |          |        |
	//            PVOut   PVSurf     PVMark   PVAVD
	//======================================================================

	// setup cross-references between library objects

	// ... This is the house that Jack built ...

	PetscFunctionBeginUser;
	// TSSol
	lm->ts.scal     = &lm->scal;
	// DBMat
	lm->dbm.scal    = &lm->scal;
	// FDSTAG
	lm->fs.scal     = &lm->scal;
	// FreeSurf
	lm->surf.jr     = &lm->jr;
	// BCCtx
	lm->bc.scal     = &lm->scal;
	lm->bc.ts       = &lm->ts;
	lm->bc.fs       = &lm->fs;
	lm->bc.dbm      = &lm->dbm;
	lm->bc.jr       = &lm->jr;
	// JacRes
	lm->jr.scal     = &lm->scal;
	lm->jr.ts       = &lm->ts;
	lm->jr.fs       = &lm->fs;
	lm->jr.surf     = &lm->surf;
	lm->jr.bc       = &lm->bc;
	lm->jr.dbm      = &lm->dbm;
	lm->jr.dbdike   = &lm->dbdike;
	// AdvCtx
	lm->actx.fs     = &lm->fs;
	lm->actx.jr     = &lm->jr;
	lm->actx.surf   = &lm->surf;
	lm->actx.dbm    = &lm->dbm;
	lm->actx.Ptr    = &lm->Ptr;
	// PVOut
	lm->pvout.jr    = &lm->jr;
	// PVSurf
	lm->pvsurf.surf = &lm->surf;
	// PVMark
	lm->pvmark.actx = &lm->actx;
	// PVPTR
	lm->pvptr.actx  = &lm->actx;
	// PVAVD
	lm->pvavd.actx  = &lm->actx;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSaveOutput(LaMEMLib *lm)
{
	//==================
	// Save data to disk
	//==================

	Scaling        *scal;
	TSSol          *ts;
	PetscScalar    time;
	PetscInt       bgPhase, step;
	char           *dirName;
	PetscLogDouble t;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	scal = &lm->scal;
	ts   = &lm->ts;

	if(!TSSolIsOutput(ts)) PetscFunctionReturn(0);

	PrintStart(&t, "Saving output", NULL);

	time    = ts->time*scal->time;
	step    = ts->istep;
	bgPhase = lm->actx.bgPhase;



	// create directory (encode current time & step number)
	asprintf(&dirName, "Timestep_%1.8lld_%1.8e", (LLD)step, time);

	// create output directory
	ierr = DirMake(dirName); CHKERRQ(ierr);

	// AVD phase output
	ierr = PVAVDWriteTimeStep(&lm->pvavd, dirName, time); CHKERRQ(ierr);

	// grid ParaView output
	ierr = PVOutWriteTimeStep(&lm->pvout, dirName, time); CHKERRQ(ierr);

	// free surface ParaView output
	ierr = PVSurfWriteTimeStep(&lm->pvsurf, dirName, time); CHKERRQ(ierr);

	// marker ParaView output
	ierr = PVMarkWriteTimeStep(&lm->pvmark, dirName, time); CHKERRQ(ierr);

	// compute and output effective permeability
	ierr = JacResGetPermea(&lm->jr, bgPhase, step, lm->pvout.outfile); CHKERRQ(ierr);

	// passive tracers paraview output
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		// save .dat files// binary of passive tracers
		ierr = PVPtrWriteTimeStep(&lm->pvptr, dirName, time); CHKERRQ(ierr);

	}
	// clean up
	free(dirName);

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSolve(LaMEMLib *lm, void *param, PetscLogStage stages[4])
{
	PMat           pm;     // preconditioner matrix    (to be removed!)
	PCStokes       pc;     // Stokes preconditioner    (to be removed!)
	NLSol          nl;     // nonlinear solver context (to be removed!)
 	AdjGrad        aop;    // Adjoint options          (to be removed!)
	SNES           snes;   // PETSc nonlinear solver
	PetscInt       restart;
	PetscLogDouble t;


	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// create Stokes preconditioner, matrix and nonlinear solver
	ierr = PMatCreate(&pm, &lm->jr);    CHKERRQ(ierr);
	ierr = PCStokesCreate(&pc, pm);     CHKERRQ(ierr);
	ierr = NLSolCreate(&nl, pc, &snes); CHKERRQ(ierr);

	//==============
	// INITIAL GUESS
	//==============
	PetscCall(PetscLogStagePush(stages[0])); /* Start profiling stage*/

	ierr = LaMEMLibInitGuess(lm, snes); CHKERRQ(ierr);

	PetscCall(PetscLogStagePop()); /* Stop profiling stage*/

	if (param)
	{
		ierr = AdjointCreate(&aop, &lm->jr, (ModParam *)param); CHKERRQ(ierr);
	}

	//===============
	// TIME STEP LOOP
	//===============

	while(!TSSolIsDone(&lm->ts))
	{
		//====================================
		//	NONLINEAR THERMO-MECHANICAL SOLVER
		//====================================

		// apply phase transitions on particles
		ierr = Phase_Transition(&lm->actx); CHKERRQ(ierr);
		
		// initialize boundary constraint vectors
		ierr = BCApply(&lm->bc); CHKERRQ(ierr);

	
		// initialize temperature
		ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);

		// compute elastic parameters
		ierr = JacResGetI2Gdt(&lm->jr); CHKERRQ(ierr);

		// solve nonlinear equation system with SNES
		PetscTime(&t);

		PetscCall(PetscLogStagePush(stages[1])); /* Start profiling stage*/

		ierr = SNESSolve(snes, NULL, lm->jr.gsol); CHKERRQ(ierr);

		PetscCall(PetscLogStagePop()); /* Stop profiling stage*/
		// print analyze convergence/divergence reason & iteration count
		ierr = SNESPrintConvergedReason(snes, t); CHKERRQ(ierr);

		// view nonlinear residual
		ierr = JacResViewRes(&lm->jr); CHKERRQ(ierr);

		// Compute adjoint gradients every TS
		if (param)
		{
			
			ModParam      *IOparam;
			IOparam       = (ModParam *)param;	
			if (IOparam->use == _adjointgradients_ || IOparam->use == _gradientdescent_ || IOparam->use == _inversion_ )
			{	/* 	Compute the adjoint gradients 
				 	
					This is done here, as the adjoint should be cmputed with the current residual that does not take advection etc.
					into account. It does compute it every dt; one can perhaps only activate it for the last dt.
				*/
				ierr = AdjointObjectiveAndGradientFunction(&aop, &lm->jr, &nl, (ModParam *)param, snes, &lm->surf); CHKERRQ(ierr);
			}
		}

		//==========================================
		// MARKER & FREE SURFACE ADVECTION + EROSION
		//==========================================

		PetscCall(PetscLogStagePush(stages[2])); /* Start profiling stage*/

		// calculate current time step
		ierr = ADVSelectTimeStep(&lm->actx, &restart); CHKERRQ(ierr);
		
		// restart if fixed time step is larger than CFLMAX
		if(restart) continue;

		// advect free surface
		ierr = FreeSurfAdvect(&lm->surf); CHKERRQ(ierr);

		// advect markers
		ierr = ADVAdvect(&lm->actx); CHKERRQ(ierr);

		// apply background strain-rate "DWINDLAR" BC (Bob Shaw "Ship of Strangers")
		ierr = BCStretchGrid(&lm->bc); CHKERRQ(ierr);

		// exchange markers between the processors (after mesh advection)
		ierr = ADVExchange(&lm->actx); CHKERRQ(ierr);

		// Advect Passive tracers
		ierr = ADVAdvectPassiveTracer(&lm->actx); CHKERRQ(ierr);

		PetscCall(PetscLogStagePop()); /* Stop profiling stage*/

		int SurfaceMode;
		SurfaceMode = SURFACE;

		if( 1 == SurfaceMode )
		// using LaMEM original code to calculate topography
		{
			PetscPrintf(PETSC_COMM_WORLD, "\n Calculating surface process through LaMEM original code \n");
			// apply erosion to the free surface
			ierr = FreeSurfAppErosion(&lm->surf); CHKERRQ(ierr);

			// apply sedimentation to the free surface
			ierr = FreeSurfAppSedimentation(&lm->surf); CHKERRQ(ierr);
		}

		if( 2 == SurfaceMode)
		// Using FastScape to calculate topography
		{
			PetscPrintf(PETSC_COMM_WORLD, "\n Calculating surface process through FastScape \n");
			ierr = fastscape(&lm->surf, &lm->actx); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "\n FastScape Done \n");
		}


		// remap markers onto (stretched) grid
		ierr = ADVRemap(&lm->actx); CHKERRQ(ierr);

		// update phase ratios taking into account actual free surface position
		ierr = FreeSurfGetAirPhaseRatio(&lm->surf); CHKERRQ(ierr);

		//==================
		// Save data to disk
		//==================
	
		// update time stamp and counter
		ierr = TSSolStepForward(&lm->ts); CHKERRQ(ierr);

		PetscCall(PetscLogStagePush(stages[3])); /* Start profiling stage*/

		// grid & marker output
		ierr = LaMEMLibSaveOutput(lm); CHKERRQ(ierr);

		PetscCall(PetscLogStagePop()); /* Stop profiling stage*/

		// restart database
		ierr = LaMEMLibSaveRestart(lm); CHKERRQ(ierr);

	}

	//======================
	// END OF TIME STEP LOOP
	//======================

	if (param)
	{

		ModParam      *IOparam;
		IOparam       = (ModParam *)param;

		if(IOparam->use == _syntheticforwardrun_)
		{	// Assume this as a forward simulation and save the solution vector
	 		//VecDuplicate(lm->jr.gsol, &IOparam->xini);
			//VecCopy(lm->jr.gsol, IOparam->xini);
		}

		ierr = AdjointDestroy (&aop,  IOparam);  	CHKERRQ(ierr);

	}

	// destroy objects
	ierr = PCStokesDestroy(pc);    			CHKERRQ(ierr);
	ierr = PMatDestroy    (pm);    			CHKERRQ(ierr);
	ierr = SNESDestroy    (&snes); 			CHKERRQ(ierr);
	ierr = NLSolDestroy   (&nl);   			CHKERRQ(ierr);

	// save marker database
	ierr = ADVMarkSave(&lm->actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDryRun(LaMEMLib *lm)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize boundary constraint vectors
	ierr = BCApply(&lm->bc); CHKERRQ(ierr);

	// initialize temperature
	ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);

	// compute inverse elastic parameters (dependent on dt)
	ierr = JacResGetI2Gdt(&lm->jr); CHKERRQ(ierr);

	// evaluate initial residual
	ierr = JacResFormResidual(&lm->jr, lm->jr.gsol, lm->jr.gres); CHKERRQ(ierr);

	// save output for inspection
	ierr = LaMEMLibSaveOutput(lm); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibInitGuess(LaMEMLib *lm, SNES snes)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PetscLogDouble t;

	// initialize boundary constraint vectors
	ierr = BCApply(&lm->bc); CHKERRQ(ierr);

	// initialize temperature
	ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);

	// solve for steady-state temperature (if requested)
	ierr = LaMEMLibDiffuseTemp(lm); CHKERRQ(ierr);

	// initialize pressure
	ierr = JacResInitPres(&lm->jr,&lm->ts); CHKERRQ(ierr);

	// lithostatic pressure initializtaion
	ierr = JacResInitLithPres(&lm->jr, &lm->actx, &lm->ts); CHKERRQ(ierr);

	// compute inverse elastic parameters (dependent on dt)
	ierr = JacResGetI2Gdt(&lm->jr); CHKERRQ(ierr);

	if(lm->jr.ctrl.initGuess)
	{
		PetscPrintf(PETSC_COMM_WORLD, "============================== INITIAL GUESS =============================\n");
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		// solve nonlinear equation system with SNES
		PetscTime(&t);

		ierr = SNESSolve(snes, NULL, lm->jr.gsol); CHKERRQ(ierr);

		// print analyze convergence/divergence reason & iteration count
		ierr = SNESPrintConvergedReason(snes, t); CHKERRQ(ierr);

		// view nonlinear residual
		ierr = JacResViewRes(&lm->jr); CHKERRQ(ierr);

		// switch flag
		lm->jr.ctrl.initGuess = 0;
	}
	else
	{
		// evaluate initial residual
		ierr = JacResFormResidual(&lm->jr, lm->jr.gsol, lm->jr.gres); CHKERRQ(ierr);
	}

	// save output for inspection
	ierr = LaMEMLibSaveOutput(lm); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDiffuseTemp(LaMEMLib *lm)
{
	JacRes         *jr;
	TSSol          *ts;
	Controls       *ctrl;
	AdvCtx         *actx;
	PetscLogDouble t;
	PetscScalar    diff_step;
	PetscInt       i, num_steps;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ts       = &lm->ts; 
	jr      = &lm->jr;
	ctrl    = &jr->ctrl;
	actx    = &lm->actx;

	// check for infinite diffusion
	if (ctrl->actTemp && ctrl->actSteadyTemp && ts->istep==0)
	{
		PrintStart(&t,"Computing steady-state temperature distribution", NULL);

		// ignore existing temperature initialization
		ierr = VecZeroEntries(jr->lT); CHKERRQ(ierr);
		ierr = JacResApplyTempBC(jr); CHKERRQ(ierr);

		// compute steady-state temperature distribution
		ierr = LaMEMLibSolveTemp(lm, 0.0); CHKERRQ(ierr);

		// overwrite markers where T(phase) is set
		ierr = ADVMarkSetTempPhase(actx); CHKERRQ(ierr);

		// project temperature from markers to grid
		ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);

		// initialize temperature
		ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);
		
		PrintDone(t);
	}

	// check for additional limited diffusion
	if (ctrl->actTemp && ctrl->steadyTempStep && ts->istep==0)
	{
		PrintStart(&t,"Diffusing temperature", NULL);

		diff_step = ctrl->steadyTempStep;
		num_steps = 1;

		if (ctrl->steadyNumStep)
		{
			num_steps = ctrl->steadyNumStep;
			diff_step = diff_step/((PetscScalar) num_steps);
		}
		
		for(i=0;i<num_steps;i++)
		{
			// diffuse
			ierr = LaMEMLibSolveTemp(lm, diff_step); CHKERRQ(ierr);

			// reset temperature in anomalous phases every step
			if (ctrl->actHeatRech > 1)
			{
				// overwrite markers where T(phase) is set
				ierr = ADVMarkSetTempPhase(actx); CHKERRQ(ierr);

				// project temperature from markers to grid
				ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);
	
				// initialize temperature
				ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);
			}
		}

		// reset Temperature in anomalous phase
		if (ctrl->actHeatRech)
		{
			// overwrite markers where T(phase) is set
			ierr = ADVMarkSetTempPhase(actx); CHKERRQ(ierr);

			// project temperature from markers to grid
			ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);
	
			// initialize temperature
			ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);
		}
		
		PrintDone(t);		
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSolveTemp(LaMEMLib *lm, PetscScalar dt)
{
	JacRes         *jr;
	AdvCtx         *actx;
	KSP            tksp;
	
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	jr   = &lm->jr;
	actx = &lm->actx;
	
	// create temperature diffusion solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &tksp); CHKERRQ(ierr);

	// enable geometric multigrid
	PetscCall(KSPSetDM(tksp, jr->DA_T));
	PetscCall(KSPSetDMActive(tksp, PETSC_FALSE));

	// set options
	ierr = KSPSetOptionsPrefix(tksp,"its_");   CHKERRQ(ierr);
	ierr = KSPSetFromOptions(tksp);            CHKERRQ(ierr);

	// compute matrix and rhs
	// STEADY STATE solution is activated by setting time step to zero
	ierr = JacResGetTempRes(jr, dt); CHKERRQ(ierr);
	ierr = JacResGetTempMat(jr, dt); CHKERRQ(ierr);

	// solve linear system
	ierr = KSPSetOperators(tksp, jr->Att, jr->Att); CHKERRQ(ierr);
	ierr = KSPSetUp(tksp);                          CHKERRQ(ierr);
	ierr = KSPSolve(tksp, jr->ge, jr->dT);          CHKERRQ(ierr);

	// destroy initial temperature solver
	ierr = KSPDestroy(&tksp); CHKERRQ(ierr);

	// store computed temperature, enforce boundary constraints
	ierr = JacResUpdateTemp(jr); CHKERRQ(ierr);

	// copy temperature to markers
	ierr = ADVMarkSetTempVector(actx); CHKERRQ(ierr);

	// project temperature from markers to grid
	ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);
	
	// initialize temperature
	ierr = JacResInitTemp(&lm->jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//	ObjFunct objf;   // objective function
//	ierr = ObjFunctCreate(&objf, &IOparam, &lm->surf, fb); CHKERRQ(ierr);
//	ierr = ObjFunctDestroy(&objf); CHKERRQ(ierr);

//---------------------------------------------------------------------------
