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
#ifndef __tssolve_h__
#define __tssolve_h__
//---------------------------------------------------------------------------

struct FB;
struct Scaling;

//---------------------------------------------------------------------------

struct TSSol
{
	Scaling *scal;

	PetscScalar dt;                        // time step
	PetscScalar dt_next;                   // next time step (CFL or dt_max)
	PetscScalar dt_min;                    // minimum time step (declare divergence if lower value is attempted)
	PetscScalar dt_max;                    // maximum time step (if CFL is larger, truncate)
	PetscScalar dt_out;                    // output step (output at least at fixed time intervals)
	PetscScalar inc_dt;                    // time step increment per time step (fraction of unit)
	PetscInt    num_dtper;                 // number of time stepping periods
	PetscScalar t_dtper[_max_periods_+1];  // timestamps where timestep should be fixed
	PetscScalar dt_dtper[_max_periods_+1]; // target timesteps ar timestamps above
	PetscScalar schedule[_max_num_steps_]; // time stepping schedule 
	PetscScalar CFL;                       // CFL (Courant-Friedrichs-Lewy) criterion
	PetscScalar CFLMAX;                    // CFL tolerance for accepting fixed time steps
	PetscScalar time;                      // current time
	PetscScalar time_out;                  // previous output time stamp
	PetscScalar time_end;                  // simulation end time
	PetscScalar tol;                       // tolerance for time comparisons
	PetscInt    nstep_max;                 // maximum allowed number of steps (lower bound: time_end/dt_max)
	PetscInt    nstep_out;                 // save output every n steps
	PetscInt    nstep_ini;                 // save output for n initial steps
	PetscInt    nstep_rdb;                 // save restart database every n steps
	PetscInt    fix_dt;                    // flag to keep time steps fixed for advection (elasticity, kinematic block BC)
	PetscInt    istep;                     // time step counter
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

PetscErrorCode TSSolGetPeriodSteps(
	PetscScalar  dt_start, // timestep at the start of the period
	PetscScalar  dt_end,   // timestep at the end of the period
	PetscScalar  span,     // time span of period
	PetscScalar *dt,       // time steps in period
	PetscInt    &n);       // number of time steps

PetscErrorCode TSSolMakeSchedule(TSSol *ts);

PetscErrorCode TSSolAdjustSchedule(TSSol *ts, PetscScalar dt_cfl, PetscInt istep, PetscScalar *schedule);

//---------------------------------------------------------------------------

// compute CFL time step with limit
#define GET_CFL_STEP(dt, dtmax, CFL, gidtmax) \
	{ (gidtmax) ? dt = CFL/gidtmax : dt = dtmax;  if(dt > dtmax) dt = dtmax; }

//---------------------------------------------------------------------------
#endif


