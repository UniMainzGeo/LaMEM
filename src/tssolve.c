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
	ierr = getScalarParam(fb, _OPTIONAL_, "time_end",  &ts->time_end,  1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "dt_max",    &ts->dt_max,    1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt",        &ts->dt,        1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_min",    &ts->dt_min,    1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_out",    &ts->dt_out,    1, time);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "inc_dt",    &ts->inc_dt,    1, 1.0 );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFL",       &ts->CFL,       1, 1.0 );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFLMAX",    &ts->CFLMAX,    1, 1.0 );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_max", &ts->nstep_max, 1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_out", &ts->nstep_out, 1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_ini", &ts->nstep_ini, 1, -1  );  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_rdb", &ts->nstep_rdb, 1, -1  );  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "time_tol",  &ts->tol,       1, 1.0 );  CHKERRQ(ierr);

	if(ts->CFL < 0.0 && ts->CFL > 1.0)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "CFL parameter must be between 0 and 1");
	}

	if(ts->CFLMAX < 0.0 && ts->CFLMAX > 1.0)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "CFLMAX parameter must be between 0 and 1");
	}

	if(!ts->time_end && !ts->nstep_max)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Define at least one of the parameters: time_end, nstep_max");
	}

	// set defaults
	if(!ts->dt)        ts->dt        = ts->dt_max/5.0;
	if(!ts->dt_min)    ts->dt_min    = ts->dt_max/50.0;
	if(!ts->nstep_max) ts->nstep_max = 50*(PetscInt)PetscCeilReal(ts->time_end/ts->dt_max);
	if(!ts->time_end)  ts->time_end  = ((PetscScalar)ts->nstep_max)*ts->dt_max;

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Time stepping parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Simulation end time     : %7.5f %s \n",  ts->time_end*time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Time step               : %7.5f %s \n",  ts->dt      *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Minimum time step       : %7.5f %s \n",  ts->dt_min  *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Maximum time step       : %7.5f %s \n",  ts->dt_max  *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Output time step        : %7.5f %s \n",  ts->dt_out  *time, scal->lbl_time);
	PetscPrintf(PETSC_COMM_WORLD, "   Time step increment     : %7.5f %s \n",  ts->inc_dt,        scal->lbl_unit);
	PetscPrintf(PETSC_COMM_WORLD, "   CFL criterion           : %7.5f %s \n",  ts->CFL,           scal->lbl_unit);
	PetscPrintf(PETSC_COMM_WORLD, "   CFLMAX for elasticity   : %7.5f %s \n",  ts->CFLMAX,        scal->lbl_unit);
	PetscPrintf(PETSC_COMM_WORLD, "   Maximum number of steps : %lld \n", (LLD)ts->nstep_max);
	PetscPrintf(PETSC_COMM_WORLD, "   Output interval         : %lld \n", (LLD)ts->nstep_out);
	PetscPrintf(PETSC_COMM_WORLD, "   Output initial steps    : %lld \n", (LLD)ts->nstep_ini);
	PetscPrintf(PETSC_COMM_WORLD, "   Restart interval        : %lld \n", (LLD)ts->nstep_rdb);
	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolSetupElasticity"
PetscErrorCode TSSolSetupElasticity(TSSol *ts)
{
	// activate elasticity

	PetscFunctionBegin;

	// set initial elastic time step
	ts->dt_elast = ts->dt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar TSSolGetScaledTime(TSSol *ts)
{
	return ts->time*ts->scal->time;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolStepForward"
PetscErrorCode TSSolStepForward(TSSol *ts)
{
	// This function is called before output

	PetscFunctionBegin;

	// update time
	ts->time += ts->dt;

	// update time step index
	ts->istep++;

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
	//  * for the fixed number of initial steps
	//  * after fixed number of steps
	//  * after fixed time interval
	//==========================================

	PetscScalar tol;

	// get tolerance
	tol = ts->tol*ts->dt_max;

	// check output conditions
	if((ts->nstep_ini &&   ts->istep < ts->nstep_ini)
	|| (ts->nstep_out && !(ts->istep % ts->nstep_out))
	|| (ts->dt_out    &&   ts->time  >= ts->time_out + ts->dt_out - tol))
	{
		// update output time stamp
		ts->time_out = ts->time;

		return 1;
	}

	return 0;
}
//---------------------------------------------------------------------------
PetscInt TSSolIsDone(TSSol *ts)
{
	//==========================================
	// stop the simulation:
	//
	//  * end time is reached
	//  * maximum number of time steps is done
	//
	// plot time stamp otherwise
	//==========================================

	Scaling     *scal;
	PetscScalar  tol;
	PetscInt     done;

	scal = ts->scal;

	// get tolerance
	tol = ts->tol*ts->dt_max;

	if(ts->time  >= ts->time_end - tol
	|| ts->istep == ts->nstep_max)
	{
		PetscPrintf(PETSC_COMM_WORLD, "===========================================================\n");
		PetscPrintf(PETSC_COMM_WORLD, "................... SOLUTION IS DONE!......................\n");
		PetscPrintf(PETSC_COMM_WORLD, "===========================================================\n");

		done = 1;
	}
	else
	{
		// output time step information
		PetscPrintf(PETSC_COMM_WORLD, "===========================================================\n");
		PetscPrintf(PETSC_COMM_WORLD, "........................ STEP: %lld .......................\n", (LLD)ts->istep);
		PetscPrintf(PETSC_COMM_WORLD, "===========================================================\n");
		PetscPrintf(PETSC_COMM_WORLD, " Current time: %7.5f %s \n", ts->time*scal->time, scal->lbl_time);

		done = 0;
	}

	return done;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolSelectStep"
PetscErrorCode TSSolSelectStep(
	TSSol       *ts,
	PetscScalar gidtmax,  // maximum global inverse time step
	PetscInt    *restart) // time step restart flag
{
	Scaling     *scal;
	PetscScalar  dt_cfl, dt_cfl_max;

	PetscFunctionBegin;

	scal = ts->scal;

	// set restart flag
	(*restart) = 0;

	// compute CFL time step
	dt_cfl = TSSolGetCFLStep(ts, ts->CFL, gidtmax);

	// declare divergence if too small time step is attempted
	if(dt_cfl < ts->dt_min)
	{
		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too small CFL time step: dt_cfl=%7.5f dt_min=%7.5f %s",
			dt_cfl*scal->time, ts->dt_min*scal->time, scal->lbl_time);
	}

	if(ts->dt_elast)
	{
		//============================
		// visco-elasto-plastic case
		//============================

		// check whether current elastic time step is larger than CFL
		if(ts->dt_elast > dt_cfl)
		{
			dt_cfl_max = TSSolGetCFLStep(ts, ts->CFLMAX, gidtmax);

			if(ts->dt_elast > dt_cfl_max)
			{
				ts->dt_elast = dt_cfl;

				(*restart) = 1;

				PetscPrintf(PETSC_COMM_WORLD, " Elastic time step exceeds CFLMAX level : %7.5f %s \n", dt_cfl_max, scal->lbl_time);
				PetscPrintf(PETSC_COMM_WORLD, " RESTARTING TIME STEP! \n");

				PetscFunctionReturn(0);
			}
			else
			{
				PetscPrintf(PETSC_COMM_WORLD, " Elastic time step exceeds CFL level : %7.5f %s \n", dt_cfl*scal->time, scal->lbl_time);
			}
		}

		// use current elastic step for advection
		ts->dt = ts->dt_elast;

		// compute elastic time step
		ts->dt_elast *= (1.0 + ts->inc_dt);
		if(ts->dt_elast > dt_cfl) ts->dt_elast = dt_cfl;
	}
	else
	{
		//===================
		// visco-plastic case
		//===================

		// compute advection step
		ts->dt *= (1.0 + ts->inc_dt);
		if(ts->dt > dt_cfl) ts->dt = dt_cfl;
	}

	// print time step information
	PetscPrintf(PETSC_COMM_WORLD, " Time increment : %7.5f %s \n", ts->dt*scal->time, scal->lbl_time);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar TSSolGetCFLStep(TSSol *ts, PetscScalar CFL, PetscScalar gidtmax)
{
	PetscScalar dt;

	if(gidtmax)         dt = CFL/gidtmax;
	else                dt = ts->dt_max;
	if(dt > ts->dt_max) dt = ts->dt_max;

	return dt;
}
//---------------------------------------------------------------------------
