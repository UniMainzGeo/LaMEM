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

struct FB;
struct Scaling;

//---------------------------------------------------------------------------

struct TSSol
{
	Scaling *scal;

	PetscScalar dt;        // time step
	PetscScalar dt_next;   // next time step (CFL or dt_max)
	PetscScalar dt_min;    // minimum time step (declare divergence if lower value is attempted)
	PetscScalar dt_max;    // maximum time step
	PetscScalar dt_out;    // output step (output at least at fixed time intervals)
	PetscScalar inc_dt;    // time step increment per time step (fraction of unit)
	PetscScalar CFL;       // CFL (Courant-Friedrichs-Lewy) criterion
	PetscScalar CFLMAX;    // CFL tolerance for accepting fixed time steps
	PetscScalar time;      // current time
	PetscScalar time_out;  // output time stamp
	PetscScalar time_end;  // simulation end time
	PetscScalar tol;       // tolerance for time comparisons
	PetscInt    nstep_max; // maximum allowed number of steps (lower bound: time_end/dt_max)
	PetscInt    nstep_out; // save output every n steps
	PetscInt    nstep_ini; // save output for n initial steps
	PetscInt    nstep_rdb; // save restart database every n steps
	PetscInt    fix_dt;    // flag to keep time steps fixed for advection (elasticity, kinematic block BC)
	PetscInt    istep;     // time step counter

};

//---------------------------------------------------------------------------

PetscErrorCode TSSolCreate(TSSol *ts, FB *fb);

PetscInt TSSolIsDone(TSSol *ts);

PetscErrorCode TSSolStepForward(TSSol *ts);

PetscInt TSSolIsRestart(TSSol *ts);

PetscInt TSSolIsOutput(TSSol *ts);

PetscErrorCode TSSolGetCFLStep(
	TSSol       *ts,
	PetscScalar  gidtmax,  // maximum global inverse time step
	PetscInt    *restart); // time step restart flag

//---------------------------------------------------------------------------

// compute CFL time step with limit
#define GET_CFL_STEP(dt, dtmax, CFL, gidtmax) \
	{ (gidtmax) ? dt = CFL/gidtmax : dt = dtmax;  if(dt > dtmax) dt = dtmax; }

//---------------------------------------------------------------------------
#endif


