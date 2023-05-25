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
//......................   TIME STEPPING PARAMETERS   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tssolve.h"
#include "parsing.h"
#include "scaling.h"
#include "tools.h"

//---------------------------------------------------------------------------
PetscErrorCode TSSolCreate(TSSol *ts, FB *fb)
{
	Scaling     *scal;
	PetscScalar  time;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
	ierr = getScalarParam(fb, _OPTIONAL_, "time_end",        &ts->time_end,   1,               time);          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "dt_max",          &ts->dt_max,     1,               time);          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt",              &ts->dt,         1,               time);          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_min",          &ts->dt_min,     1,               time);          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "dt_out",          &ts->dt_out,     1,               time);          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "inc_dt",          &ts->inc_dt,     1,               1.0 );          CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "num_dt_periods",  &ts->num_dtper,  1,               _max_periods_); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "time_dt_periods",  ts->t_dtper,    ts->num_dtper+1, time);           CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "step_dt_periods",  ts->dt_dtper,   ts->num_dtper+1, time);           CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFL",             &ts->CFL,        1,               1.0 );          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "CFLMAX",          &ts->CFLMAX,     1,               1.0 );          CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_max",       &ts->nstep_max,  1,               -1  );          CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_out",       &ts->nstep_out,  1,               -1  );          CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_ini",       &ts->nstep_ini,  1,               -1  );          CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "nstep_rdb",       &ts->nstep_rdb,  1,               -1  );          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "time_tol",        &ts->tol,        1,               1.0 );          CHKERRQ(ierr);

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
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dt should be between dt_min and dt_max");
	}

	if(ts->num_dtper)
	{
		ierr = TSSolMakeSchedule(ts);
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
		PetscPrintf(PETSC_COMM_WORLD, "Current time        : %7.8f %s \n", ts->time*scal->time, scal->lbl_time);
		PetscPrintf(PETSC_COMM_WORLD, "Tentative time step : %7.8f %s \n", ts->dt  *scal->time, scal->lbl_time);
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		done = 0;
	}

	return done;
}
//---------------------------------------------------------------------------
PetscErrorCode TSSolStepForward(TSSol *ts)
{
	// function is called in the end of time step before output and restart

	PetscFunctionBeginUser;

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
PetscErrorCode TSSolGetCFLStep(
	TSSol       *ts,
	PetscScalar  gidtmax, // maximum global inverse time step
	PetscInt    *restart) // time step restart flag
{
	Scaling     *scal;
	PetscScalar  dt_cfl, dt_cfl_max;
	PetscScalar *schedule;
	PetscInt     istep;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get context
	scal     = ts->scal;
	schedule = ts->schedule;
	istep    = ts->istep;

	// set restart flag
	(*restart) = 0;

	// get CFL time step
	GET_CFL_STEP(dt_cfl, ts->dt_max, ts->CFL, gidtmax)

	// declare divergence if too small time step is required
	if(dt_cfl < ts->dt_min)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Time step is smaller than dt_min: %7.5f %s\n", ts->dt_min*scal->time, scal->lbl_time);
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
	if(ts->num_dtper)
	{
		ts->dt_next = schedule[istep];

		// check CFL limit
		if(ts->dt_next > dt_cfl)
		{
			// adjust timestep
			ts->dt_next = dt_cfl;
			
			// adjust schedule
			ierr = TSSolAdjustSchedule(ts, dt_cfl, istep, schedule); CHKERRQ(ierr);
		} 
	}
	else
	{
		ts->dt_next = ts->dt*(1.0 + ts->inc_dt);

		// check CFL limit
		if(ts->dt_next > dt_cfl) ts->dt_next = dt_cfl;
	}

	// apply immediately if time step is not fixed (otherwise apply in the end of time step)
	if(!ts->fix_dt) ts->dt = ts->dt_next;

	// print time step information
	PetscPrintf(PETSC_COMM_WORLD, "Actual time step : %7.5f %s \n", ts->dt*scal->time, scal->lbl_time);

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode TSSolGetPeriodSteps(
	PetscScalar  dt_start, // timestep at the start of the period
	PetscScalar  dt_end,   // timestep at the end of the period
	PetscScalar  span,     // time span of period
	PetscScalar *dt,       // time steps in period
	PetscInt    &n)        // number of time steps
{
	PetscScalar  dt_avg, n_try, sum, err, corr;
	PetscInt     i;
	
	PetscFunctionBeginUser;

	// average timestep
	dt_avg = (dt_start + dt_end) / 2.0;

	// approximate number of steps
	n_try  = span / dt_avg;

	// actual number of steps
	n     = (PetscInt)max(1, (int)round(n_try));

	// make proposal for steps
	linSpace(dt_start,dt_end,n+1,dt);

	// how far are we off?
	sum    = 0;
	for(i = 0; i < n; i++) 
	{
		sum += dt[i];
	}
	err    = span - sum;

	// correction per step
	corr   = err / n;

	// add correction
	for(i = 0; i < n; i++)
	{
		dt[i] += corr;
	}

	// warning
	if(n < 2)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Warning: Only one transition step in time step schedule.\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode TSSolMakeSchedule(TSSol *ts)
{
	PetscScalar *schedule, *steps, *t, *dt_fix;
	PetscScalar  dt_start, dt_end, span;
	PetscInt     num_seg, iSeg, iter, n, i, maxSteps;
	
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access content
	num_seg   = ts->num_dtper;
	t         = ts->t_dtper;
	dt_fix    = ts->dt_dtper;
	maxSteps  = ts->nstep_max;

	// allocate
	ierr = PetscMalloc1((size_t)_max_num_steps_*sizeof(PetscScalar), &schedule); CHKERRQ(ierr);
	ierr = PetscMalloc1((size_t)_max_num_steps_*sizeof(PetscScalar), &steps);    CHKERRQ(ierr);
	ierr = PetscMemzero(schedule, (size_t)_max_num_steps_*sizeof(PetscScalar));  CHKERRQ(ierr);

	// loop through segments and build schedule
	iter = 0; n = 0;
	for(iSeg = 0; iSeg < num_seg; iSeg++)
	{
		// read input
		dt_start = dt_fix[iSeg];
		dt_end   = dt_fix[iSeg+1];
		span     = t[iSeg+1] - t[iSeg];

		// check input
		if(!(span > 0.0))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "time_dt_periods must be strinctly increasing.");
		}
		if(!(dt_start > 0.0) || !(dt_end > 0.0))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "step_dt_periods must be larger than 0.");
		}

		// get timesteps
		ierr = PetscMemzero(steps, (size_t)_max_num_steps_*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = TSSolGetPeriodSteps(dt_start, dt_end, span, steps, n);

		// add to schedule
		for(i = 0; i < n; i++)
		{
			schedule[iter] = steps[i];
			iter++;
		}
	}	
	schedule[iter] = dt_fix[iSeg];

	// use schedule
	maxSteps      = min(maxSteps, iter+1);
	ts->nstep_max = maxSteps;
	for(i = 0; i < maxSteps; i++)
	{
		ts->schedule[i] = schedule[i];
	}	

	// free memory
	ierr = PetscFree(steps);    CHKERRQ(ierr);
	ierr = PetscFree(schedule); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode TSSolAdjustSchedule(TSSol *ts, PetscScalar dt_cfl, PetscInt istep, PetscScalar *schedule)
{
	PetscScalar diff;
	PetscInt    maxSteps, i;

	PetscFunctionBeginUser;

	// access content
	maxSteps = ts->nstep_max;

	// difference between target and limit
	diff = schedule[istep] - dt_cfl;

	// adjust current step
	schedule[istep] -= diff;

	// adapt schedule
	if(diff < 0.25*schedule[istep+1])
	{
		// make next target bigger
		schedule[istep+1] += diff;
	}
	else
	{
		// squeeze in new time step to close the gap
		for(i = min(maxSteps, (PetscInt)(_max_num_steps_)-1); i > istep; i--)
		{
			schedule[i+1] = schedule[i];
		}
		schedule[istep+1] = diff;
		ts->nstep_max = maxSteps + 1;
	}	

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
