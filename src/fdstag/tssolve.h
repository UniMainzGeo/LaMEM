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

} TSSol;

//---------------------------------------------------------------------------

PetscErrorCode TSSolSetUp(TSSol *ts, UserCtx *usr);

PetscErrorCode TSSolUpdate(TSSol *ts, Scaling *scal, PetscBool *done);

//---------------------------------------------------------------------------
#endif
