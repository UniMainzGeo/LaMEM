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
//........................   TIME STEPPING SOLVER   .........................
//---------------------------------------------------------------------------
#ifndef __tssolve_h__
#define __tssolve_h__
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt    nstep;   // maximum number of steps
	PetscScalar dtmax;   // maximum time step
//	PetscScalar dtmin;   // minimum time step
	PetscScalar Cmax;    // dimensionless Courant number (should be {significantly} less than unit)
//	PetscScalar timeEnd; // duration of simulation

	PetscInt    istep;   // current step index
	PetscScalar pdt;     // previous time step
	PetscScalar dt;      // current time step (to be defined)
	PetscScalar time;    // current time


/*
	// time-stepping
	PetscInt         save_timesteps;
	PetscScalar      CFL;
	PetscInt         time_end;
	PetscScalar      dt_max;
	PetscScalar      dt;

 */

} TSSol;

//---------------------------------------------------------------------------

PetscErrorCode TSSolSetUp(TSSol *ts);

PetscErrorCode TSSolUpdate(TSSol *ts, Scaling *scal, PetscBool *done);

//---------------------------------------------------------------------------
#endif
