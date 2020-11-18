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
 **    filename:   LaMEMLib.h
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
	Scaling  scal;   // scaling
	TSSol    ts;     // time-stepping controls
	DBMat    dbm;    // material database
	FDSTAG   fs;     // staggered-grid layout
	FreeSurf surf;   // free-surface grid
	BCCtx    bc;     // boundary condition context
	AdvCtx   actx;   // advection context
	JacRes   jr;     // Jacobian & residual context
	P_Tr     Ptr  ;  // passive tracers
	PVOut    pvout;  // paraview output driver
	PVSurf   pvsurf; // paraview output driver for surface
	PVMark   pvmark; // paraview output driver for markers
	PVAVD    pvavd;  // paraview output driver for AVD
	PVPtr    pvptr;  // paraview out passive tracers
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
