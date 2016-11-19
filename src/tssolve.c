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
 **    filename:   tssolve.c
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
#include "LaMEM.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
//---------------------------------------------------------------------------
// * change time stepping for requested time interval, not integer number of steps
// * break-point files
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolSetUp"
PetscErrorCode TSSolSetUp(TSSol *ts)
{
	PetscFunctionBegin;
/*
	ts->nstep = usr->time_end; // maximum number of steps
	ts->dtmax = usr->dt_max;   // maximum time step
	ts->Cmax  = usr->CFL;      // Courant number
	ts->istep = 0;             // current step index
	ts->pdt   = 0.0;           // previous time step
	ts->dt    = usr->dt;       // current time step (to be defined)
	ts->time  = 0.0;
*/
	PetscPrintf(PETSC_COMM_WORLD, " CFL timestep factor            : %f \n", ts->Cmax);

	if(ts->Cmax > 1.0)
	{
		SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_USER, " Courant step length Cmax=%7.5f is larger than allowed (%7.5f).", ts->Cmax, 0.5);
	}

	if(ts->Cmax > 0.7)
	{
		PetscPrintf(PETSC_COMM_WORLD, " WARNING! Large Courant step length Cmax=%7.5f. Consider reducing.\n", ts->Cmax);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolUpdate"
PetscErrorCode TSSolUpdate(TSSol *ts, Scaling *scal, PetscBool *done)
{
	PetscFunctionBegin;

	//----------------------------------------
	// update time state, print info to screen
	//----------------------------------------

	// update time
	ts->time += ts->dt;

	// update time index
	ts->istep++;

	// print time info
	PetscPrintf(PETSC_COMM_WORLD," Time = %g%s, dt = %g%s \n",
		ts->time*scal->time, scal->lbl_time,
		ts->dt  *scal->time, scal->lbl_time);

	PetscPrintf(PETSC_COMM_WORLD," Finished timestep %lld out of %lld \n",(LLD)ts->istep-1, (LLD)ts->nstep-1);
	PetscPrintf(PETSC_COMM_WORLD," \n");

	// WARNING! CHECK TIME, NOT INDEX, IN THE FUTURE!

	// check to stop time-step loop
	if(ts->istep == ts->nstep) (*done) = PETSC_TRUE;
	else                       (*done) = PETSC_FALSE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
