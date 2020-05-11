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
//......................   TIME STEPPING PARAMETERS   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tssolve.h"
#include "parsing.h"
#include "scaling.h"
#include "tools.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolCreate"
PetscErrorCode TSSolCreate(TSSol *ts, FB *fb)
{
	Scaling     *scal;
	PetscScalar  time;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = ts->scal;
	time = scal->time;

	// set defaults
	ts->inc_dt    = 0.1;
	ts->CFL       = 0.5;
	ts->CFLMAX    = 0.8;
	ts->nstep_out = 1;
	ts->nstep_ini = 1;
	ts->tol       = 1e-8;

	// read parameters
	ierr = getScalarParam(fb, _OPTIONAL_, "time_end",  &ts->time_end,   1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "dt_max",    &ts->dt_max,     1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt",        &ts->dt,         1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_min",    &ts->dt_min,     1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_out",    &ts->dt_out,     1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "inc_dt",    &ts->inc_dt,     1, 1.0 );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFL",       &ts->CFL,        1, 1.0 );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFLMAX",    &ts->CFLMAX,     1, 1.0 );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_max", &ts->nstep_max,  1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_out", &ts->nstep_out,  1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_ini", &ts->nstep_ini,  1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_rdb", &ts->nstep_rdb,  1, -1  );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "time_tol",  &ts->tol,        1, 1.0 );  CHKERRQ(ierr);

	if(ts->CFL < 0.0 && ts->CFL > 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "CFL parameter must be between 0 and 1");
	}

	if(ts->CFLMAX < 0.0 && ts->CFLMAX > 1.0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "CFLMAX parameter must be between 0 and 1");
	}

	if(ts->CFL > ts->CFLMAX)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "CFL parameter should be smaller than CFLMAX");
	}

	if(!ts->time_end && !ts->nstep_max)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define at least one of the parameters: time_end, nstep_max");
	}

	// set defaults
	if(!ts->dt)        ts->dt        = ts->dt_max/5.0;
	if(!ts->dt_min)    ts->dt_min    = ts->dt_max/50.0;
	if(!ts->nstep_max) ts->nstep_max = 50*(PetscInt)PetscCeilReal(ts->time_end/ts->dt_max);
	if(!ts->time_end)  ts->time_end  = ((PetscScalar)ts->nstep_max)*ts->dt_max;

	if(ts->dt_min > ts->dt_max)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dt_max should be larger than dt_min");
	}

	if(!(ts->dt >= ts->dt_min && ts->dt <= ts->dt_max))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dt should lay between dt_min and dt_max");
	}

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Time stepping parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Simulation end time          : %g %s \n", ts->time_end*time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Maximum number of steps      : %lld \n", (LLD)ts->nstep_max);
	PetscPrintf(PETSC_COMM_WORLD, "   Time step                    : %g %s \n", ts->dt      *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Minimum time step            : %g %s \n", ts->dt_min  *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Maximum time step            : %g %s \n", ts->dt_max  *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Time step increase factor    : %g \n",    ts->inc_dt);
	PetscPrintf(PETSC_COMM_WORLD, "   CFL criterion                : %g \n",    ts->CFL);
    PetscPrintf(PETSC_COMM_WORLD, "   CFLMAX (fixed time steps)    : %g \n",    ts->CFLMAX);

	if(ts->dt_out)    PetscPrintf(PETSC_COMM_WORLD, "   Output time step             : %g %s \n", ts->dt_out  *time, scal->lbl_time);
	if(ts->nstep_out) PetscPrintf(PETSC_COMM_WORLD, "   Output every [n] steps       : %lld \n", (LLD)ts->nstep_out);
	if(ts->nstep_ini) PetscPrintf(PETSC_COMM_WORLD, "   Output [n] initial steps     : %lld \n", (LLD)ts->nstep_ini);
	if(ts->nstep_rdb) PetscPrintf(PETSC_COMM_WORLD, "   Save restart every [n] steps : %lld \n", (LLD)ts->nstep_rdb);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscInt TSSolIsDone(TSSol *ts)
{
	//=================================================
	// stop the simulation:
	//
	//  * end time is reached
	//  * maximum number of time steps is done
	//
	// plot time stamp otherwise
	//
	// function is called in the beginning of time step
	//=================================================

	Scaling     *scal;
	PetscScalar  time_end;
	PetscInt     done;

	// access context
	scal = ts->scal;

	// get end time (with tolerance)
	time_end = ts->time_end - ts->tol*ts->dt_max;

	if(ts->time  >= time_end
	|| ts->istep == ts->nstep_max)
	{
		PetscPrintf(PETSC_COMM_WORLD, "=========================== SOLUTION IS DONE! ============================\n");
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		done = 1;
	}
	else
	{
		// output time step information
		PrintStep(ts->istep + 1);
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD, "Current time        : %7.5f %s \n", ts->time*scal->time, scal->lbl_time);
		PetscPrintf(PETSC_COMM_WORLD, "Tentative time step : %7.5f %s \n", ts->dt  *scal->time, scal->lbl_time);
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		done = 0;
	}

	return done;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolStepForward"
PetscErrorCode TSSolStepForward(TSSol *ts)
{
	// function is called in the end of time step before output and restart

	PetscFunctionBegin;

	// update time
	ts->time += ts->dt;

	// update time step counter
	ts->istep++;

	// apply tentative time step (at latest here)
	ts->dt = ts->dt_next;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscInt TSSolIsRestart(TSSol *ts)
{
	// save restart database after fixed number of steps
	if(ts->nstep_rdb && !(ts->istep % ts->nstep_rdb)) return 1;

	return 0;
}
//---------------------------------------------------------------------------
PetscInt TSSolIsOutput(TSSol *ts)
{
	//==========================================
	// save output:
	//
	//  * for initial guess
	//  * for the fixed number of initial steps
	//  * after fixed number of steps
	//  * after fixed time interval
	//==========================================

	PetscScalar time_out;

	// get next output time (with tolerance)
	time_out = ts->time_out + ts->dt_out - ts->tol*ts->dt_max;

	// check output conditions
	if(!ts->istep
	|| (ts->nstep_ini &&   ts->istep <= ts->nstep_ini)
	|| (ts->nstep_out && !(ts->istep %  ts->nstep_out))
	|| (ts->dt_out    &&   ts->time  >= time_out))
	{
		// update output time stamp
		ts->time_out = ts->time;

		return 1;
	}

	return 0;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolGetCFLStep"
PetscErrorCode TSSolGetCFLStep(
	TSSol       *ts,
	PetscScalar  gidtmax, // maximum global inverse time step
	PetscInt    *restart) // time step restart flag
{
	Scaling     *scal;
	PetscScalar  dt_cfl, dt_cfl_max;

	PetscFunctionBegin;

	// get context
	scal = ts->scal;

	// set restart flag
	(*restart) = 0;

	// get CFL time step
	GET_CFL_STEP(dt_cfl, ts->dt_max, ts->CFL, gidtmax)

	// declare divergence if too small time step is required
	if(dt_cfl < ts->dt_min)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Time step is smaller than dt_min: %7.5f %s\n", ts->dt_min*scal->time, scal->lbl_time);
	}

	// check fixed time step restrictions
	if(ts->fix_dt)
	{
		// get CFLMAX time step
		GET_CFL_STEP(dt_cfl_max, ts->dt_max, ts->CFLMAX, gidtmax)

		// restart if fixed time step is too large (elasticity, kinematic block BC)
		if(ts->dt > dt_cfl_max)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Time step exceeds CFLMAX level: %7.5f %s\n", dt_cfl_max*scal->time, scal->lbl_time);
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
			PetscPrintf(PETSC_COMM_WORLD, "***********************   RESTARTING TIME STEP!   ************************\n");
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

			ts->dt = dt_cfl;

			(*restart) = 1;

			PetscFunctionReturn(0);
		}
		else if(ts->dt > dt_cfl)
		{
			// print warning
			PetscPrintf(PETSC_COMM_WORLD, "Time step exceeds CFL level: %7.5f %s\n", dt_cfl*scal->time, scal->lbl_time);
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
		}
	}

	// compute tentative time step
	ts->dt_next = ts->dt*(1.0 + ts->inc_dt);

	// check CFL limit
	if(ts->dt_next > dt_cfl) ts->dt_next = dt_cfl;

	// apply immediately if time step is not fixed (otherwise apply in the end of time step)
	if(!ts->fix_dt) ts->dt = ts->dt_next;

	// print time step information
	PetscPrintf(PETSC_COMM_WORLD, "Actual time step : %7.5f %s \n", ts->dt*scal->time, scal->lbl_time);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
