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
#include "matData.h"
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
#include "phase_transition.h"
#include "passive_tracer.h"
#include "LaMEMLib.h"

//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibMain(void *param, FB *fb)
{
	LaMEMLib       lm;
	RunMode        mode;
	PetscBool      found;
	PetscInt       exists;
	char           str[_str_len_];
	PetscLogDouble cputime_start, cputime_end;

	
	PetscFunctionBeginUser;       

	// start code
	PetscCall(PetscTime(&cputime_start));

	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"                   Lithosphere and Mantle Evolution Model                   \n");
	PetscPrintf(PETSC_COMM_WORLD,"     Compiled: Date: %s - Time: %s 	    \n",__DATE__,__TIME__ );
	PetscPrintf(PETSC_COMM_WORLD,"     Version : 3.0.0 \n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD,"        STAGGERED-GRID FINITE DIFFERENCE CANONICAL IMPLEMENTATION           \n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");

	// read run mode
	mode = _NORMAL_;

	PetscCall(PetscOptionsGetCheckString("-mode", str, &found));

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
		PetscCall(DirCheck("./restart", &exists));

		if(!exists)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No restart database available (check -mode option)");
		}
	}

	//===========
	// INITIALIZE
	//===========

	// clear
	PetscCall(PetscMemzero(&lm, sizeof(LaMEMLib)));

	// setup cross-references between library objects
	PetscCall(LaMEMLibSetLinks(&lm));

	if(mode == _SAVE_GRID_)
	{
		// save grid & exit
		PetscCall(LaMEMLibSaveGrid(&lm, fb));

		PetscFunctionReturn(0);
	}
	if(mode == _NORMAL_ || mode == _DRY_RUN_)
	{
		// create library objects
		PetscCall(LaMEMLibCreate(&lm, param, fb));
	}
	else if(mode == _RESTART_)
	{
		// open restart database
		PetscCall(LaMEMLibLoadRestart(&lm, fb));
	}

	//======
	// SOLVE
	//======

	if(mode == _DRY_RUN_)
	{
		// compute initial residual, output & stop
		PetscCall(LaMEMLibDryRun(&lm));
	}
	else if(mode == _NORMAL_ || mode == _RESTART_)
	{
		// solve coupled nonlinear equations
		PetscCall(LaMEMLibSolve(&lm, param));
	}

	// destroy library objects
	PetscCall(LaMEMLibDestroy(&lm));

	PetscTime(&cputime_end);

	PetscPrintf(PETSC_COMM_WORLD, "Total solution time : %g (sec) \n", cputime_end - cputime_start);
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibCreate(LaMEMLib *lm, void *param, FB *fb)
{
	UNUSED(param);

	
	PetscFunctionBeginUser;

	// create scaling object
	PetscCall(ScalingCreate(&lm->scal, fb, PETSC_TRUE));

	// create time stepping object
	PetscCall(TSSolCreate(&lm->ts, fb));

	// create parallel grid
	PetscCall(FDSTAGCreate(&lm->fs, fb));

	// create material database
	PetscCall(DBMatCreate(&lm->dbm, fb, PETSC_TRUE));

	// create free surface grid
	PetscCall(FreeSurfCreate(&lm->surf, fb));

	// create boundary condition context
	PetscCall(BCCreate(&lm->bc, fb));

	// create residual & Jacobian evaluation context
	PetscCall(JacResCreate(&lm->jr, fb));

	// create dike database
	PetscCall(DBDikeCreate(&lm->dbdike, &lm->dbm, fb, &lm->jr, PETSC_TRUE));

	// initialize arrays for dynamic phase transition
	PetscCall(DynamicPhTr_Init(&lm->jr));

	// create advection context
	PetscCall(ADVCreate(&lm->actx, fb));

	// create passive tracers
	PetscCall(ADVPtrPassive_Tracer_create(&lm->actx,fb));

	// create output object for all requested output variables
	PetscCall(PVOutCreate(&lm->pvout, fb));

	// create output object for the free surface
	PetscCall(PVSurfCreate(&lm->pvsurf, fb));

	// create output object for the markers - for debugging
	PetscCall(PVMarkCreate(&lm->pvmark, fb));

	// create output object for the passive tracers
	PetscCall(PVPtrCreate(&lm->pvptr, fb));

	// AVD output driver
	PetscCall(PVAVDCreate(&lm->pvavd, fb));

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSaveGrid(LaMEMLib *lm, FB *fb)
{
	
	PetscFunctionBeginUser;

	// create scaling object
	PetscCall(ScalingCreate(&lm->scal, fb, PETSC_TRUE));

	// create parallel grid
	PetscCall(FDSTAGCreate(&lm->fs, fb));

	// save processor partitioning
	PetscCall(FDSTAGSaveGrid(&lm->fs));

	// destroy parallel grid
	PetscCall(FDSTAGDestroy(&lm->fs));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibLoadRestart(LaMEMLib *lm, FB *fb)
{
	FILE            *fp;
	PetscLogDouble  t;
	PetscMPIInt     rank;
	char            *fileName;

	
	PetscFunctionBeginUser;

	PrintStart(&t, "Loading restart database", NULL);

	// get MPI processor rank
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

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
	PetscCall(LaMEMLibSetLinks(lm));

	// staggered grid
	PetscCall(FDSTAGReadRestart(&lm->fs, fp));

	// free surface
	PetscCall(FreeSurfReadRestart(&lm->surf, fp));

	// boundary conditions context
	PetscCall(BCReadRestart(&lm->bc, fp));

	// solution variables
	PetscCall(JacResReadRestart(&lm->jr, fp));

	// markers
	PetscCall(ADVReadRestart(&lm->actx, fp));

	// passive tracers read restart
	PetscCall(ReadPassive_Tracers(&lm->actx,fp));

	// main output driver
	PetscCall(PVOutCreateData(&lm->pvout));

	// surface output driver
	PetscCall(PVSurfCreateData(&lm->pvsurf));
	
	// arrays for dynamic NotInAir phase_trans
	PetscCall(DynamicPhTr_ReadRestart(&lm->jr, fp));
	
	// read from input file, create arrays for dynamic diking, and read from restart file
	PetscCall(DynamicDike_ReadRestart(&lm->dbdike, &lm->dbm, &lm->jr, &lm->ts, fp, fb));
 
	// override material database
	PetscCall(DBMatCreate(&lm->dbm, fb, PETSC_TRUE));

	// close temporary restart file
	fclose(fp);

	// free space
	free(fileName);

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

	
	PetscFunctionBeginUser;

	if(!TSSolIsRestart(&lm->ts)) PetscFunctionReturn(0);

	PrintStart(&t, "Saving restart database", NULL);

	// get MPI processor rank
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

	// compile actual & temporary restart file name
	asprintf(&fileNameTmp, "./restart-tmp/rdb.%1.8lld.dat", (LLD)rank);

	// create temporary restart directory
	PetscCall(DirMake("./restart-tmp"));

	// open temporary restart file for writing in binary mode
	fp = fopen(fileNameTmp, "wb");

	if(fp == NULL)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open restart file %s\n", fileNameTmp);
	}

	// write LaMEM library database
	fwrite(lm, sizeof(LaMEMLib), 1, fp);

	// staggered grid
	PetscCall(FDSTAGWriteRestart(&lm->fs, fp));

	// free surface
	PetscCall(FreeSurfWriteRestart(&lm->surf, fp));

	// boundary conditions context
	PetscCall(BCWriteRestart(&lm->bc, fp));

	// solution variables
	PetscCall(JacResWriteRestart(&lm->jr, fp));

	// markers
	PetscCall(ADVWriteRestart(&lm->actx, fp));

	// passive tracers
	PetscCall(Passive_Tracer_WriteRestart(&lm->actx, fp));

	// dynamic phase transition 
	PetscCall(DynamicPhTr_WriteRestart(&lm->jr, fp));

	// dynamic dike 
	PetscCall(DynamicDike_WriteRestart(&lm->jr, fp));

	// close temporary restart file
	fclose(fp);

	// delete existing restart database
	PetscCall(LaMEMLibDeleteRestart());

	// push temporary database to actual
	PetscCall(DirRename("./restart-tmp", "./restart"));

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

	
	PetscFunctionBeginUser;

	// get MPI processor rank
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

	asprintf(&fileName, "./restart/rdb.%1.8lld.dat", (LLD)rank);

	// check for existing restart database
	PetscCall(DirCheck("./restart", &exists));

	if(exists)
	{
		// delete existing database
		status = remove(fileName);

		if(status && errno != ENOENT)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to delete file %s", fileName);
		}

		PetscCall(DirRemove("./restart"));
	}

	// free space
	free(fileName);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDestroy(LaMEMLib *lm)
{
	
	PetscFunctionBeginUser;

	PetscCall(FDSTAGDestroy      (&lm->fs));
	PetscCall(FreeSurfDestroy    (&lm->surf));
	PetscCall(BCDestroy          (&lm->bc));
	PetscCall(JacResDestroy      (&lm->jr));
	PetscCall(ADVPtrDestroy      (&lm->actx));
	PetscCall(ADVDestroy         (&lm->actx));
	PetscCall(PVOutDestroy       (&lm->pvout));
	PetscCall(PVSurfDestroy      (&lm->pvsurf));
	PetscCall(DynamicPhTrDestroy (&lm->dbm));
	PetscCall(DynamicDike_Destroy(&lm->jr));

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
	PetscCall(DirMake(dirName));

	// AVD phase output
	PetscCall(PVAVDWriteTimeStep(&lm->pvavd, dirName, time));

	// grid ParaView output
	PetscCall(PVOutWriteTimeStep(&lm->pvout, dirName, time));

	// free surface ParaView output
	PetscCall(PVSurfWriteTimeStep(&lm->pvsurf, dirName, time));

	// marker ParaView output
	PetscCall(PVMarkWriteTimeStep(&lm->pvmark, dirName, time));

	// compute and output effective permeability
	PetscCall(JacResGetPermea(&lm->jr, bgPhase, step, lm->pvout.outfile));

	// passive tracers paraview output
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		// save .dat files// binary of passive tracers
		PetscCall(PVPtrWriteTimeStep(&lm->pvptr, dirName, time));
	}
	// clean up
	free(dirName);

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibSolve(LaMEMLib *lm, void *param)
{
	SNES           snes;   // PETSc nonlinear solver
 	AdjGrad        aop;    // Adjoint options
	PetscInt       restart, track_stages;
	PetscLogDouble t;
	PetscLogStage  stages[4];

	
	PetscFunctionBeginUser;

	if(!param)
	{
		// normal mode only
		track_stages = 1;

		// name computational stages
		PetscCall(PetscLogStageRegister("Initial guess",  &stages[0]));
		PetscCall(PetscLogStageRegister("SNES solve",     &stages[1]));
		PetscCall(PetscLogStageRegister("Advect markers", &stages[2]));
		PetscCall(PetscLogStageRegister("I/O",            &stages[3]));
	}
	else
	{	// not for the inversion!
		track_stages = 0;
	}

	// create nonlinear solver
	PetscCall(NLSolCreate(&snes, &lm->jr));

	//==============
	// INITIAL GUESS
	//==============

	if(track_stages) { PetscCall(PetscLogStagePush(stages[0])); }

	PetscCall(LaMEMLibInitGuess(lm, snes));

	if(track_stages) { PetscCall(PetscLogStagePop()); }

	if(param)
	{
		PetscCall(AdjointCreate(&aop, &lm->jr, (ModParam *)param));
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
		PetscCall(Phase_Transition(&lm->actx));
		
		// initialize boundary constraint vectors
		PetscCall(BCApply(&lm->bc));

		// initialize temperature
		PetscCall(JacResInitTemp(&lm->jr));

		// compute elastic parameters
		PetscCall(JacResGetI2Gdt(&lm->jr));

		// solve nonlinear equation system with SNES
		PetscTime(&t);

		if(track_stages) { PetscCall(PetscLogStagePush(stages[1])); }

		PetscCall(SNESSolve(snes, NULL, lm->jr.gsol));

		if(track_stages) { PetscCall(PetscLogStagePop()); }

		// print analyze convergence/divergence reason & iteration count
		PetscCall(SNESPrintConvergedReason(snes, t));

		// view nonlinear residual
		PetscCall(JacResViewRes(&lm->jr));

		// Compute adjoint gradients every TS
		if(param)
		{
			ModParam *IOparam = (ModParam *)param;
			
			if(IOparam->use == _adjointgradients_ || IOparam->use == _gradientdescent_ || IOparam->use == _inversion_ )
			{
				//================================================================================================================
				// Compute the adjoint gradients
				//
				//  This is done here, as the adjoint should be computed with the current residual that does not take advection etc.
				//  into account. It does compute it every dt; one can perhaps only activate it for the last dt.
				//================================================================================================================

				PetscCall(AdjointObjectiveAndGradientFunction(&aop, (ModParam*)param, snes));
			}
		}

		//==========================================
		// MARKER & FREE SURFACE ADVECTION + EROSION
		//==========================================

		if(track_stages) { PetscCall(PetscLogStagePush(stages[2])); }

		// calculate current time step
		PetscCall(ADVSelectTimeStep(&lm->actx, &restart));
		
		// restart if fixed time step is larger than CFLMAX
		if(restart) continue;

		// advect free surface
		PetscCall(FreeSurfAdvect(&lm->surf));

		// advect markers
		PetscCall(ADVAdvect(&lm->actx));

		// apply background strain-rate "DWINDLAR" BC (Bob Shaw "Ship of Strangers")
		PetscCall(BCStretchGrid(&lm->bc));

		// exchange markers between the processors (after mesh advection)
		PetscCall(ADVExchange(&lm->actx));

		// Advect Passive tracers
		PetscCall(ADVAdvectPassiveTracer(&lm->actx));

		if(track_stages) { PetscCall(PetscLogStagePop()); }

		// apply erosion to the free surface
		PetscCall(FreeSurfAppErosion(&lm->surf));

		// apply sedimentation to the free surface
		PetscCall(FreeSurfAppSedimentation(&lm->surf));

		// apply topographic diffusion to the free surface
		PetscCall(FreeSurfAppTopoDiffusion(&lm->surf));

		// remap markers onto (stretched) grid
		PetscCall(ADVRemap(&lm->actx));

		// update phase ratios taking into account actual free surface position
		PetscCall(FreeSurfGetAirPhaseRatio(&lm->surf));

		//==================
		// Save data to disk
		//==================
	
		// update time stamp and counter
		PetscCall(TSSolStepForward(&lm->ts));

		if(track_stages) { PetscCall(PetscLogStagePush(stages[3])); }

		// grid & marker output
		PetscCall(LaMEMLibSaveOutput(lm));

		if(track_stages) { PetscCall(PetscLogStagePop()); }

		// restart database
		PetscCall(LaMEMLibSaveRestart(lm));

	}

	//======================
	// END OF TIME STEP LOOP
	//======================

	if(param)
	{
		ModParam *IOparam = (ModParam *)param;

		if(IOparam->use == _syntheticforwardrun_)
		{
			// Assume this as a forward simulation and save the solution vector
			// VecDuplicate(lm->jr.gsol, &IOparam->xini);
			// VecCopy(lm->jr.gsol, IOparam->xini);
		}

		PetscCall(AdjointDestroy(&aop, IOparam));
	}

	// destroy objects
	PetscCall(NLSolDestroy(&snes));

	// save marker database
	PetscCall(ADVMarkSave(&lm->actx));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibDryRun(LaMEMLib *lm)
{
	
	PetscFunctionBeginUser;

	// initialize boundary constraint vectors
	PetscCall(BCApply(&lm->bc));

	// initialize temperature
	PetscCall(JacResInitTemp(&lm->jr));

	// compute inverse elastic parameters (dependent on dt)
	PetscCall(JacResGetI2Gdt(&lm->jr));

	// evaluate initial residual
	PetscCall(JacResFormResidual(&lm->jr, lm->jr.gsol, lm->jr.gres));

	// save output for inspection
	PetscCall(LaMEMLibSaveOutput(lm));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode LaMEMLibInitGuess(LaMEMLib *lm, SNES snes)
{
	
	PetscFunctionBeginUser;

	PetscLogDouble t;

	// initialize boundary constraint vectors
	PetscCall(BCApply(&lm->bc));

	// initialize temperature
	PetscCall(JacResInitTemp(&lm->jr));

	// solve for steady-state temperature (if requested)
	PetscCall(LaMEMLibDiffuseTemp(lm));

	// initialize pressure
	PetscCall(JacResInitPres(&lm->jr,&lm->ts));

	// lithostatic pressure initializtaion
	PetscCall(JacResInitLithPres(&lm->jr, &lm->actx, &lm->ts));

	// compute inverse elastic parameters (dependent on dt)
	PetscCall(JacResGetI2Gdt(&lm->jr));

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	if(lm->jr.ctrl.initGuess)
	{
		PetscPrintf(PETSC_COMM_WORLD, "============================== INITIAL GUESS =============================\n");
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		// solve nonlinear equation system with SNES
		PetscTime(&t);

		PetscCall(SNESSolve(snes, NULL, lm->jr.gsol));

		// print analyze convergence/divergence reason & iteration count
		PetscCall(SNESPrintConvergedReason(snes, t));

		// view nonlinear residual
		PetscCall(JacResViewRes(&lm->jr));

		// switch flag
		lm->jr.ctrl.initGuess = 0;
	}
	else
	{
		// evaluate initial residual
		PetscCall(JacResFormResidual(&lm->jr, lm->jr.gsol, lm->jr.gres));
	}

	// save output for inspection
	PetscCall(LaMEMLibSaveOutput(lm));

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

	
	PetscFunctionBeginUser;

	// access context
	ts       = &lm->ts; 
	jr      = &lm->jr;
	ctrl    = &jr->ctrl;
	actx    = &lm->actx;

	// check for infinite diffusion
	if(ctrl->actTemp && ctrl->actSteadyTemp && ts->istep==0)
	{
		PrintStart(&t,"Computing steady-state temperature distribution", NULL);

		// ignore existing temperature initialization
		PetscCall(VecZeroEntries(jr->lT));
		PetscCall(JacResApplyTempBC(jr));

		// compute steady-state temperature distribution
		PetscCall(LaMEMLibSolveTemp(lm, 0.0));

		// overwrite markers where T(phase) is set
		PetscCall(ADVMarkSetTempPhase(actx));

		// project temperature from markers to grid
		PetscCall(ADVProjHistMarkToGrid(actx));

		// initialize temperature
		PetscCall(JacResInitTemp(&lm->jr));
		
		PrintDone(t);
	}

	// check for additional limited diffusion
	if(ctrl->actTemp && ctrl->steadyTempStep && ts->istep==0)
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
			PetscCall(LaMEMLibSolveTemp(lm, diff_step));

			// reset temperature in anomalous phases every step
			if (ctrl->actHeatRech > 1)
			{
				// overwrite markers where T(phase) is set
				PetscCall(ADVMarkSetTempPhase(actx));

				// project temperature from markers to grid
				PetscCall(ADVProjHistMarkToGrid(actx));
	
				// initialize temperature
				PetscCall(JacResInitTemp(&lm->jr));
			}
		}

		// reset Temperature in anomalous phase
		if(ctrl->actHeatRech)
		{
			// overwrite markers where T(phase) is set
			PetscCall(ADVMarkSetTempPhase(actx));

			// project temperature from markers to grid
			PetscCall(ADVProjHistMarkToGrid(actx));
	
			// initialize temperature
			PetscCall(JacResInitTemp(&lm->jr));
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
	Controls       *ctrl;
	KSP            tksp;
	PetscScalar    norm;
	PetscBool      set;
	PetscInt       ts_ksp_atol_auto;
	
	
	PetscFunctionBeginUser;

	// access context
	jr   = &lm->jr;
	actx = &lm->actx;
	ctrl = &jr->ctrl;

	// get automatic absolute tolerance initialization flag
	PetscCall(PetscOptionsHasName(NULL, NULL, "-ts_ksp_atol_auto", &set));
	
	if(set && ctrl->actTemp) { ts_ksp_atol_auto = 1; }
	else                     { ts_ksp_atol_auto = 0; }

	// create temperature diffusion solver
	PetscCall(KSPCreate(PETSC_COMM_WORLD, &tksp));

	// enable geometric multigrid
	PetscCall(KSPSetDM(tksp, jr->DA_T));

#if PETSC_VERSION_LT(3, 25, 0)
	PetscCall(KSPSetDMActive(tksp,                   PETSC_FALSE));
#else
	PetscCall(KSPSetDMActive(tksp, KSP_DMACTIVE_ALL, PETSC_FALSE));
#endif

	// set options
	PetscCall(KSPSetOptionsPrefix(tksp,"its_"));
	PetscCall(KSPSetFromOptions(tksp));

	// compute matrix and rhs
	// STEADY STATE solution is activated by setting time step to zero
	PetscCall(JacResGetTempRes(jr, dt));
	PetscCall(JacResGetTempMat(jr, dt));

	// update reference norm for automatic tolerance selection
	if(ts_ksp_atol_auto)
	{
		PetscCall(VecNorm(jr->ge, NORM_2, &norm));

		if(norm > jr->ts_ksp_ref_norm) { jr->ts_ksp_ref_norm = norm; }
	}

	// solve linear system
	PetscCall(KSPSetOperators(tksp, jr->Att, jr->Att));
	PetscCall(KSPSetUp(tksp));
	PetscCall(KSPSolve(tksp, jr->ge, jr->dT));

	// view solver
	if(!dt)
	{
		PetscCall(ViewSolver(tksp));
	}

	// destroy initial temperature solver
	PetscCall(KSPDestroy(&tksp));

	// store computed temperature, enforce boundary constraints
	PetscCall(JacResUpdateTemp(jr));

	// copy temperature to markers
	PetscCall(ADVMarkSetTempVector(actx));

	// project temperature from markers to grid
	PetscCall(ADVProjHistMarkToGrid(actx));
	
	// initialize temperature
	PetscCall(JacResInitTemp(&lm->jr));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
