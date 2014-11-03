//---------------------------------------------------------------------------
//........................   TIME STEPPING SOLVER   .........................
//---------------------------------------------------------------------------
#ifndef __tssolve_h__
#define __tssolve_h__
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar pdt;    // previous time step
	PetscScalar dt;     // current time step (to be defined)
	PetscScalar dtmin;  // minimum time step
	PetscScalar dtmax;  // maximum time step
	PetscScalar Cmax;   // dimensionless Courant number (should be {significantly} less than unit)

} TSSol;

//---------------------------------------------------------------------------

PetscErrorCode TSSolStepup(TSSol *ts, JacRes *jr);

PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax);

PetscErrorCode TSSolGetCourantStep(TSSol *ts, JacRes *jr);

//---------------------------------------------------------------------------
#endif
