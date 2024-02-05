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
//........................ LaMEM Library major context .......................
//---------------------------------------------------------------------------
#ifndef __LaMEMLib_h__
#define __LaMEMLib_h__
//---------------------------------------------------------------------------

enum RunMode
{
	//==================
	// simulation modes
	//==================

	_NORMAL_,    // start new simulation
	_RESTART_,   // start from restart database (if available)
	_DRY_RUN_,   // initialize model, output & stop
	_SAVE_GRID_, // write parallel grid to a file & stop

};

//---------------------------------------------------------------------------

struct LaMEMLib
{
	Scaling  scal;               // scaling
	TSSol    ts;                 // time-stepping controls
	DBMat    dbm;                // material database
    DBPropDike dbdike;           // dike database
	DBPropHeatZone dbheatzone;   // heat zone database
	FDSTAG   fs;                 // staggered-grid layout
	FreeSurf surf;               // free-surface grid
	BCCtx    bc;                 // boundary condition context
	AdvCtx   actx;               // advection context
	JacRes   jr;                 // Jacobian & residual context
	P_Tr     Ptr  ;              // passive tracers
	PVOut    pvout;              // paraview output driver
	PVSurf   pvsurf;             // paraview output driver for surface
	PVMark   pvmark;             // paraview output driver for markers
	PVAVD    pvavd;              // paraview output driver for AVD
	PVPtr    pvptr;              // paraview out passive tracers
};

//---------------------------------------------------------------------------
// LAMEM LIBRARY FUNCTIONS
//---------------------------------------------------------------------------

PetscErrorCode LaMEMLibCreate(LaMEMLib *lm, void *param);

PetscErrorCode LaMEMLibSaveGrid(LaMEMLib *lm);

PetscErrorCode LaMEMLibLoadRestart(LaMEMLib *lm);

PetscErrorCode LaMEMLibSaveRestart(LaMEMLib *lm);

PetscErrorCode LaMEMLibDeleteRestart();

PetscErrorCode LaMEMLibDestroy(LaMEMLib *lm);

PetscErrorCode LaMEMLibSetLinks(LaMEMLib *lm);

PetscErrorCode LaMEMLibSaveOutput(LaMEMLib *lm, PetscInt dirInd);

PetscErrorCode LaMEMLibSolve(LaMEMLib *lm, void *param);

PetscErrorCode LaMEMLibDryRun(LaMEMLib *lm);

PetscErrorCode LaMEMLibInitGuess(LaMEMLib *lm, SNES snes);

PetscErrorCode LaMEMLibSolveTemp(LaMEMLib *lm, PetscScalar dt);

PetscErrorCode LaMEMLibDiffuseTemp(LaMEMLib *lm);

//---------------------------------------------------------------------------
#endif
