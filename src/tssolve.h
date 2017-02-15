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
 **    filename:   tssolve.h
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
//......................   TIME STEPPING PARAMETERS   .......................
//---------------------------------------------------------------------------
#ifndef __tssolve_h__
#define __tssolve_h__
//---------------------------------------------------------------------------

struct TSSol
{
	Scaling *scal;

	PetscScalar time_end;  // simulation end time
	PetscScalar dt;        // time step
	PetscScalar dt_min;    // minimum time step (declare divergence if lower value is attempted)
	PetscScalar dt_max;    // maximum time step
	PetscScalar dt_out;    // output step (output at least at fixed time intervals)
	PetscScalar inc_dt;    // time step increment per time step (fraction of unit)
	PetscScalar CFL;       // CFL (Courant-Friedrichs-Lewy) criterion
	PetscScalar CFLMAX;    // CFL criterion for elasticity
	PetscInt    nstep_max; // maximum allowed number of steps (lower bound: time_end/dt_max)
	PetscInt    nstep_out; // save output every n steps
	PetscInt    nstep_ini; // save output for n initial steps
	PetscInt    nstep_rdb; // save restart database every n steps
	PetscScalar dt_elast;  // elastic time step
	PetscScalar time;      // current time
	PetscInt    istep;     // current time step index
	PetscScalar time_out;  // output time stamp
	PetscInt    istep_out; // output time step index

};

//---------------------------------------------------------------------------

PetscErrorCode TSSolCreate(TSSol *ts, FB *fb);

PetscErrorCode TSSolSetupElasticity(TSSol *ts);

PetscScalar TSSolGetTime(TSSol *ts);

PetscInt TSSolIsOutput(TSSol *ts);

PetscInt TSSolIsRestart(TSSol *ts);

PetscInt TSSolIsDone(TSSol *ts);

PetscErrorCode TSSolStepForward(TSSol *ts);

PetscErrorCode TSSolCheckCFL(
	TSSol       *ts,
	PetscScalar  gidtmax,  // maximum global inverse time step
	PetscBool   *restart); // time step restart flag

PetscScalar TSSolGetTol(TSSol *ts);

PetscScalar TSSolGetCFLStep(TSSol *ts, PetscScalar CFL, PetscScalar gidtmax);

//---------------------------------------------------------------------------
#endif
