//---------------------------------------------------------------------------
//........................   TIME STEPPING SOLVER   .........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * change time stepping for requested time interval, not integer number of steps
// * break-point files
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolSetUp"
PetscErrorCode TSSolSetUp(TSSol *ts, UserContext *usr)
{
	PetscFunctionBegin;

	ts->nstep = usr->time_end; // maximum number of steps
	ts->dtmax = usr->dt_max;   // maximum time step
	ts->Cmax  = usr->CFL;      // Courant number
	ts->istep = 0;             // current step index
	ts->pdt   = 0.0;           // previous time step
	ts->dt    = usr->dt;       // current time step (to be defined)
	ts->time  = 0.0;
    
    
    PetscPrintf(PETSC_COMM_WORLD, " CFL timestep factor            : %f \n", ts->Cmax);
    
	if(ts->Cmax > 0.5)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, " Courant step length Cmax=%7.5f is larger than allowed (%7.5f).", ts->Cmax, 0.5);
	}

	if(ts->Cmax > 0.3)
	{
		PetscPrintf(PETSC_COMM_WORLD, " WARNING! Large Courant step length Cmax=%7.5f. Consider reducing.\n", ts->Cmax);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getMaxInvStep1DLocal"
PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax)
{
	PetscScalar v, h, vmax, idt, idtmax;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, idx, ijk[3], jj, ln;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	idtmax = (*_idtmax);

	if(ds->h_uni != PETSC_TRUE)
	{
		// compute time step on variable spacing grid
		PetscScalar ***va;

		ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, gv, &va);                     CHKERRQ(ierr);

		START_STD_LOOP
		{
			// get velocity
			v = va[k][j][i];

			// prepare node index buffer
			ijk[0] = i-sx;
			ijk[1] = j-sy;
			ijk[2] = k-sz;

			// anisotropic direction-dependent criterion
			if(v >= 0.0)  idx = ijk[dir];
			else          idx = ijk[dir]-1;

			// get mesh step
			h = ds->ncoor[idx+1] - ds->ncoor[idx];

			// get inverse time step (safe to compute)
			idt = v/h;

			// update maximum inverse time step
			if(idt > idtmax) idtmax = idt;
		}
		END_STD_LOOP

		ierr = DMDAVecRestoreArray(da, gv, &va); CHKERRQ(ierr);
	}
	else
	{
		// compute time step on uniform spacing grid
		PetscScalar *va;

		// get maximum local velocity
		ierr = VecGetLocalSize(gv, &ln); CHKERRQ(ierr);
		ierr = VecGetArray(gv, &va);     CHKERRQ(ierr);

		vmax = 0.0;
		for(jj = 0; jj < ln; jj++) { v = PetscAbsScalar(va[jj]); if(v > vmax) vmax = v;	}

		ierr = VecRestoreArray(gv, &va); CHKERRQ(ierr);

		// get inverse time step
		idt = vmax/ds->h;

        // update maximum inverse time step
		if(idt > idtmax) idtmax = idt;
	}

	// return result
	(*_idtmax) = idtmax;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TSSolGetCourantStep"
PetscErrorCode TSSolGetCourantStep(TSSol *ts, JacRes *jr)
{
	//-------------------------------------
	// compute length of the next time step
	//-------------------------------------

	FDSTAG      *fs;
	PetscScalar dt, lidtmax, gidtmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs      = jr->fs;
	lidtmax = 0.0;

	// determine maximum local inverse time step
	ierr = getMaxInvStep1DLocal(&fs->dsx, fs->DA_X, jr->gvx, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsy, fs->DA_Y, jr->gvy, 0, &lidtmax); CHKERRQ(ierr);
	ierr = getMaxInvStep1DLocal(&fs->dsz, fs->DA_Z, jr->gvz, 0, &lidtmax); CHKERRQ(ierr);

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&lidtmax, &gidtmax, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gidtmax = lidtmax;
	}
    
	// compute time step
	gidtmax /= ts->Cmax;

	if(gidtmax < 1.0/ts->dtmax) dt = ts->dtmax;
	else                        dt = 1.0/gidtmax;

	// store new time step
	ts->pdt = ts->dt;
	ts->dt  = dt;

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

	//lbl_time

	// update time
	ts->time += ts->dt;

	// update time index
	ts->istep++;

    // print time info
	PetscPrintf(PETSC_COMM_WORLD," Time = %g%s, dt = %g%s \n",
		ts->time*scal->time, scal->lbl_time,
		ts->dt  *scal->time, scal->lbl_time);

    PetscPrintf(PETSC_COMM_WORLD," Finished timestep %lld out of %lld \n",(LLD)ts->istep, (LLD)ts->nstep);

    // WARNING! CHECK TIME, NOT INDEX, IN THE FUTURE!

    // check to stop time-step loop
    if(ts->istep == ts->nstep) (*done) = PETSC_TRUE;
    else                       (*done) = PETSC_FALSE;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
