
#include "LaMEM.h"
#include "AVD.h"
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "fdstag.h"
#include "bc.h"
#include "tools.h"
#include "phase_transition.h"
#include "phase.h"
#include "constEq.h"
#include "parsing.h"
#include "objFunct.h"
#include "surf.h"
#include "subgrid.h"
#include "tssolve.h"

// allocate storage for the passive tracers
// create initial distribution of markers (every processors knows everything)
// Assign the initial phase, find the closest marker
// Advection & interpolation
// 1. Interpolate the data from the grid (vx,vy,vz)
	// a. Advect them accordingly to the local velocity field
	// b. Select the marker that are effectively within the processor and put to zero all the quantity for the others and interpolate them
	// c. Sum all the quantity in all the processors
// 2. Communicate to the master processors all the data, and print the output



//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "ADVPtrPassive_Tracer_create"
PetscErrorCode ADVPtrPassive_Tracer_create(AdvCtx *actx, FB *fb)
{
/*
 *  This function creates all the vector required for tracing pressure, temperature, phase and x,y,z position.
 *  RecvBuf is a vector used only for the synching operation and it has any meaning.
 */

	P_Tr            *passive_tr;
	char             Condition_adv[_str_len_];
	PetscInt        nummark;
	PetscErrorCode  ierr;
	PetscFunctionBegin;

	if(!actx->jr->ctrl.Passive_Tracer)	PetscFunctionReturn(0);
	passive_tr = actx->Ptr;


    // Read input parameters
	ierr = getScalarParam(fb, _REQUIRED_, "PassiveTracer_Box",          passive_tr->box_passive_tracer,         6,  1.0);       CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, "PassiveTracer_Resolution",   passive_tr->passive_tracer_resolution,  3,  0);         CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "PassiveTracer_ActiveType",   Condition_adv,      "Always");                          CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "PassiveTracer_ActiveValue",  &passive_tr->value_condition,           1, 1.0);        CHKERRQ(ierr);

	if(strcmp(Condition_adv,"Always") && !passive_tr->value_condition)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "The value of the Passive_tracers advection must be specified \n");
	}

	if(!strcmp(Condition_adv,"Always"))
	{
		passive_tr->Condition_pr = _Always_;
	}
	else if(!strcmp(Condition_adv,"Melt_Fraction"))
	{
		passive_tr->Condition_pr = _Melt_Fr_;
	}
	else if(!strcmp(Condition_adv,"Temperature"))
	{
		passive_tr->Condition_pr = _Temp_ptr_;
		passive_tr->value_condition = (passive_tr->value_condition+actx->jr->scal->Tshift)/(actx->jr->scal->temperature);
	}
	else if(!strcmp(Condition_adv,"Time"))
	{
		passive_tr->Condition_pr = _Time_ptr_;
		passive_tr->value_condition = (passive_tr->value_condition)/(actx->jr->scal->time);

	}
	else if(!strcmp(Condition_adv,"Pressure"))
	{
		passive_tr->Condition_pr = _Pres_ptr_;
		passive_tr->value_condition = (passive_tr->value_condition)/(actx->jr->scal->stress);
	}

	nummark = passive_tr->passive_tracer_resolution[0]*passive_tr->passive_tracer_resolution[1]*passive_tr->passive_tracer_resolution[2];
	passive_tr->nummark = nummark;
	if (passive_tr->nummark>_max_passive_tracer)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "The total number of passive tracers must be lower than %d",_max_passive_tracer);
	}


     PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
     PetscPrintf(PETSC_COMM_WORLD,"Passive Tracers: \n");
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box x = [Left,Right] : %6f, %6f \n",passive_tr->box_passive_tracer[0],passive_tr->box_passive_tracer[1]);
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box y = [Front,Back] : %6f, %6f \n",passive_tr->box_passive_tracer[2],passive_tr->box_passive_tracer[3]);
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box z = [Bot, Top]   : %6f, %6f \n",passive_tr->box_passive_tracer[4],passive_tr->box_passive_tracer[5]);
     PetscPrintf(PETSC_COMM_WORLD,"   # of tracers in [x,y,z] direction       : [%d, %d, %d] \n",passive_tr->passive_tracer_resolution[0],passive_tr->passive_tracer_resolution[1],passive_tr->passive_tracer_resolution[2]);
     PetscPrintf(PETSC_COMM_WORLD,"   Total # of tracers                      : %d \n",nummark);
     PetscPrintf(PETSC_COMM_WORLD,"   Tracer advection activation type        : ");
    
	 if(passive_tr->Condition_pr==_Always_)
	 {
		 PetscPrintf(PETSC_COMM_WORLD,"Always active\n");
	 }
	 else
	 {
		 if(passive_tr->Condition_pr == _Melt_Fr_)     PetscPrintf(PETSC_COMM_WORLD,"Melt_Fraction > %g     \n",    passive_tr->value_condition );
		 if(passive_tr->Condition_pr == _Temp_ptr_)    PetscPrintf(PETSC_COMM_WORLD,"Temperature > %1.0f %s    \n",    passive_tr->value_condition*actx->jr->scal->temperature-actx->jr->scal->Tshift  ,actx->jr->scal->lbl_temperature    );
		 if(passive_tr->Condition_pr == _Time_ptr_)    PetscPrintf(PETSC_COMM_WORLD,"Time > %1.1f %s           \n",    passive_tr->value_condition*actx->jr->scal->time                                ,actx->jr->scal->lbl_time           );
		 if(passive_tr->Condition_pr == _Pres_ptr_)    PetscPrintf(PETSC_COMM_WORLD,"Pressure > %1.0f %s       \n",    passive_tr->value_condition*actx->jr->scal->stress                              ,actx->jr->scal->lbl_stress         );

	 }
	 PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");


	 // Allocate memory
	 ierr = ADVPtrReCreateStorage(actx); CHKERRQ(ierr);

	 // Initialize the initial coordinate distribution and phase
	 ierr =  ADVPassiveTracerInit(actx); CHKERRQ(ierr);

	 PetscFunctionReturn(0);
	}
// ---------------------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "ADVPtrReCreateStorage"
PetscErrorCode ADVPtrReCreateStorage(AdvCtx *actx)
{

	PetscErrorCode  ierr;
	PetscFunctionBegin;

		if(!actx->jr->ctrl.Passive_Tracer)	PetscFunctionReturn(0);
	// check whether current storage is insufficient

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->ID);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->ID); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->x);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->x); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark, &actx->Ptr->y);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->y); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark, &actx->Ptr->z);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->z); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->T);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->T); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->p);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->p); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->phase);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->phase); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Melt_fr);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->Melt_fr); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->C_advection);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->C_advection); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Recv);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->Recv); CHKERRQ(ierr);

		ierr = VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Melt_Grid);      CHKERRQ(ierr);
		ierr = VecZeroEntries(actx->Ptr->Melt_Grid); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------//

#undef __FUNCT__
#define __FUNCT__ "ADVPassiveTracerInit"
PetscErrorCode ADVPassiveTracerInit(AdvCtx *actx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->jr->ctrl.Passive_Tracer == 0 ) PetscFunctionReturn(0);

	ierr = ADVPtrInitCoord(actx); CHKERRQ(ierr);

	ierr = ADV_Assign_Phase(actx); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVPtrInitCoord"
PetscErrorCode ADVPtrInitCoord(AdvCtx *actx)
{
	//Initialize the passive tracer lagrangian grid. The initial passive tracer distribution is a rectangular grid, with a
	// a variable resolution. After initializing the coordinates, phase, temperature and pressure are interpolated from
	// the nearest marker (s.s.)

	PetscScalar  x, y, z, dx, dy, dz,nx,ny,nz;
	PetscInt     i, j, k;
	PetscInt     imark;
	PetscScalar  *Xp,*Yp,*Zp,*ID,*active;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	nx = actx->Ptr->passive_tracer_resolution[0];
	ny = actx->Ptr->passive_tracer_resolution[1];
	nz = actx->Ptr->passive_tracer_resolution[2];
	dx = (actx->Ptr->box_passive_tracer[1]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length))/nx;
	dy = (actx->Ptr->box_passive_tracer[3]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length))/ny;
	dz = (actx->Ptr->box_passive_tracer[5]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[4]/(actx->dbm->scal->length))/nz;


	// marker counter
	imark = 0;
	ierr = VecGetArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->ID, &ID)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->C_advection, &active)           ; CHKERRQ(ierr);



	// create uniform distribution of markers/cell for variable grid
	for(k = 0; k < actx->Ptr->passive_tracer_resolution[2]; k++)
	{
		// spacing of particles
		for(j = 0; j < actx->Ptr->passive_tracer_resolution[1]; j++)
		{
			for(i = 0; i < actx->Ptr->passive_tracer_resolution[0]; i++)
			{
				// spacing of particles
				// loop over markers in cells
				if(k==0)
				{
					z = actx->Ptr->box_passive_tracer[4]/(actx->dbm->scal->length) + dz/2;
				}
				else
				{
					z = actx->Ptr->box_passive_tracer[4]/(actx->dbm->scal->length) + dz/2 + k*dz;
				}
				if(j==0)
				{
					y = actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length) + dy/2;
				}
				else
				{
					y = actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length) + dy/2 + j*dy;
				}
				if(i==0)
				{
					x = actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length) + dx/2;
				}
				else
				{
					x = actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length) + dx/2+i*dx;
				}

				// set marker coordinates
				Xp[imark] = x;
				Yp[imark] = y;
				Zp[imark] = z;
				ID[imark] = i+ny*j+ny*nx*k;
				if(actx->Ptr->Condition_pr == _Always_)
				{
					active[imark] = 1.0;
				}
				else
				{
					active[imark] = 0.0;
				}

				// increment local counter
				imark++;
			}
		}

	}


	ierr = VecRestoreArray(actx->Ptr->x, &Xp)                  ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->y, &Yp)                  ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->z, &Zp)                  ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->ID, &ID)                 ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->C_advection, &active)    ; CHKERRQ(ierr);



	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADV_Assign_Phase"
PetscErrorCode ADV_Assign_Phase(AdvCtx *actx)
{
	// Initially the marker are phase-less. This routine assign both phase and
	// initial temperature to the passive tracers.

	FDSTAG      *fs;
	vector <spair>    dist;
	spair d;
	PetscScalar  X[3],Xm[3],*Xp,*Yp,*Zp,*Pr,*T,*phase;
	PetscInt     I, J, K,ii,numpassive,imark,ID,M,N,n;
	PetscScalar ex,bx,ey,by,ez,bz;


	PetscErrorCode ierr;
	PetscFunctionBegin;


	fs = actx->fs;

	numpassive = actx->Ptr->nummark;

	// marker counter
	imark = 0;
	// get context
	fs = actx->fs;
	M  = fs->dsx.ncels;
	N  = fs->dsy.ncels;

	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);


	dist.reserve(_mark_buff_sz_);
	ierr = VecGetArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->p, &Pr)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->T, &T)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->phase, &phase)           ; CHKERRQ(ierr);

	for(imark=0;imark<numpassive;imark++)
	{
		// get marker coordinates
			X[0] = Xp[imark];
			X[1] = Yp[imark];
			X[2] = Zp[imark];



			if(X[0] >= bx && X[0]< ex && X[1] >= by && X[1]< ey && X[2] >= bz && X[2]< ez)
			{
			// get host cell IDs in all directions
			ierr = Discret1DFindPoint(&fs->dsx, X[0], I); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsy, X[1], J); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsz, X[2], K); CHKERRQ(ierr);

			// compute and store consecutive index
			GET_CELL_ID(ID, I, J, K, M, N);

			dist.clear();


			n = actx->markstart[ID+1] - actx->markstart[ID];

			for (ii = actx->markstart[ID]; ii <actx->markstart[ID]+n; ii++)
			{
				Xm[0] =actx->markers[ii].X[0];
				Xm[1] =actx->markers[ii].X[1];
				Xm[2] =actx->markers[ii].X[2];

				d.first  = EDIST(Xm, X);
				d.second = ii;
				dist.push_back(d);
			}

			// sort markers by distance
			sort(dist.begin(), dist.end());

			// clone closest marker
			phase[imark]= actx->markers[dist.begin()->second].phase;
			T[imark]= actx->markers[dist.begin()->second].T;
			Pr[imark]= actx->markers[dist.begin()->second].p;
			}
			else
			{

				phase[imark]= -DBL_MAX;
				T[imark]= -DBL_MAX;
				Pr[imark]= -DBL_MAX;
			}

	}

	ierr = VecRestoreArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->p, &Pr)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->T, &T)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->phase, &phase)           ; CHKERRQ(ierr);

	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = Sync_Vector(actx->Ptr->p,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->p); CHKERRQ(ierr);

		ierr = Sync_Vector(actx->Ptr->T,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->T); CHKERRQ(ierr);

		ierr = Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->phase); CHKERRQ(ierr);
	}


	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVAdvectPassiveTracer"
PetscErrorCode ADVAdvectPassiveTracer(AdvCtx *actx)
{
/*
 * Warning 1: this routine was copied from the routine of marker advection.
 * 1st : if something change in the advection routine, it may creates some discrepancies
 * in this routine (this may cause the failing of t19_passive tracers)
 * 2nd : In order to mantain a certain degree of consistency between the routine it is necessary
 * to create a general function for the advection. On the other hand, a potential solution
 * Function: 1st part: Each timestep the function advect the passive tracers whose coordinate are
 * belonging to the current processor. The passive tracers that does not belong to the processor
 * are set to be equal to -DBL_MAX. After that, the x,y,z, Temperature, Pressure are syncronized in
 * all processor with an All_reduce operation.
 * 2nd part: The routine check if the passive tracer is below the free surface, changing eventually its phase
 * and following the same approach for the coordinate and P,T.
 */

	FDSTAG          *fs;
	JacRes          *jr;
	SolVarCell      *svCell;
	Material_t      *mat;
	PData           *Pd;
	PetscInt        sx, sy, sz, nx, ny,nz,rank;
	PetscInt        jj, I, J, K, II, JJ, KK, AirPhase, num_part,ID, n, ii, numActTracers ;
	PetscScalar     ex,bx,ey,by,ez,bz;
	PetscScalar     *ncx, *ncy, *ncz;
	PetscScalar     *ccx, *ccy, *ccz;
	PetscScalar     ***lvx, ***lvy, ***lvz, ***lp, ***lT;
	PetscScalar     vx, vy, vz, xc, yc, zc, xp, yp, zp, dt, Ttop, endx,endy,endz,begx,begy,begz,npx,npy,npz;
	PetscScalar     *Xp, *Yp,*Zp,*T,*Pr,*phase,*mf_ptr,*Active,dx,dy,dz,*melt_grid;
	PetscScalar     pShift;
	PetscScalar     Xm[3],X[3];
	PetscLogDouble t;
	vector <spair>    dist;
	spair d;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	AirPhase = -1;
	Ttop     =  0.0;


	

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// access context
	fs = actx->fs;
	jr = actx->jr;
	mat = jr->dbm->phases;
	Pd  = jr->Pd;

	if(jr->ctrl.Passive_Tracer == 0)  PetscFunctionReturn(0);
	
	PrintStart(&t, "Advection Passive tracers", NULL);

	if(actx->surf->UseFreeSurf)
		{
			AirPhase = actx->surf->AirPhase;
			Ttop     = actx->jr->bc->Ttop;
		}
	// current time step
	dt = jr->ts->dt;

//	if(jr->ts->istep == 1 )ierr = ADV_Assign_Phase(actx); CHKERRQ(ierr);
	if(jr->ctrl.pShift)
	{
		pShift = jr->ctrl.pShift;
	}
	else
	{
		pShift = 0.0;
	}
	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// node & cell coordinates
	ncx = fs->dsx.ncoor; ccx = fs->dsx.ccoor;
	ncy = fs->dsy.ncoor; ccy = fs->dsy.ccoor;
	ncz = fs->dsz.ncoor; ccz = fs->dsz.ccoor;

	begx = fs->dsx.gcrdbeg;
	endx = fs->dsx.gcrdend;

	begy = fs->dsy.gcrdbeg;
	endy = fs->dsy.gcrdend;

	begz = fs->dsz.gcrdbeg;
	endz = fs->dsz.gcrdend;

	// access velocity, pressure & temperature vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,  &lp) ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT) ; CHKERRQ(ierr);

	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);


	ierr = VecGetArray(actx->Ptr->x, &Xp)              ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->y, &Yp)              ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->z, &Zp)              ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->p, &Pr)              ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->T, &T)               ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->phase, &phase)       ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->Melt_fr, &mf_ptr)    ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->Melt_Grid, &melt_grid); CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->C_advection, &Active); CHKERRQ(ierr);

	// scan all markers
    numActTracers   = 0;
	num_part        = actx->Ptr->nummark;

	for(jj = 0; jj < num_part; jj++)
	{

		// get consecutive index of the host cell
		xp = Xp[jj];
		yp = Yp[jj];
		zp = Zp[jj];

		// get marker coordinates
		if(xp >= bx && xp< ex && yp >= by && yp< ey && zp >= bz && zp< ez)
		{
			ierr = Discret1DFindPoint(&fs->dsx, xp, I); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsy, yp, J); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsz, zp, K); CHKERRQ(ierr);

			// get coordinates of cell center
			xc = ccx[I];
			yc = ccy[J];
			zc = ccz[K];

			// map marker on the cells of X, Y, Z & center grids
			if(xp > xc) { II = I; } else { II = I-1; }
			if(yp > yc) { JJ = J; } else { JJ = J-1; }
			if(zp > zc) { KK = K; } else { KK = K-1; }

			// interpolate velocity, pressure & temperature
			vx = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, xp, yp, zp, ncx, ccy, ccz);
			vy = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, xp, yp, zp, ccx, ncy, ccz);
			vz = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, xp, yp, zp, ccx, ccy, ncz);

			// update pressure & temperature variables
			Pr[jj] = InterpLin3D(lp, II, JJ, K,  sx, sy, sz, xp, yp, zp, ccx, ccy, ncz) + pShift;
			T[jj]  = InterpLin3D(lT, II, JJ, K,  sx, sy, sz, xp, yp, zp, ccx, ccy, ncz);

			GET_CELL_ID(ID, I, J, K, nx, ny)

			svCell = &jr->svCell[ID];
			dx = SIZE_CELL(I,sx,fs->dsx);
			dy = SIZE_CELL(J,sy,fs->dsy);
			dz = SIZE_CELL(K,sy,fs->dsz);

			melt_grid[jj] = svCell->svBulk.mf;

			if(svCell->svBulk.mf>0.0)
			{
			  //check if the original phase saved is one that has a phase/melt law associated


				if(mat[PetscInt(phase[jj])].pdn)
				{
					ierr = setDataPhaseDiagram(Pd, Pr[jj], T[jj], mat[PetscInt(phase[jj])].pdn); CHKERRQ(ierr);
					mf_ptr[jj]= Pd->mf;
				}
				else
				{
					// Passive tracers are initialize during the initial stage of the simulation.
					// They can have a different phase as soon as the melting start.

					// sort markers by distance
					dist.clear();
					n = actx->markstart[ID+1] - actx->markstart[ID];

					for (ii = actx->markstart[ID]; ii <actx->markstart[ID]+n; ii++)
						{
							Xm[0] = actx->markers[ii].X[0];
							Xm[1] = actx->markers[ii].X[1];
							Xm[2] = actx->markers[ii].X[2];
							X[0]  = xp;
							X[1]  = yp;
							X[2]  = zp;

							if (mat[actx->markers[ii].phase].pdn)
							{
								d.first  = EDIST(Xm, X);
								d.second = ii;
								dist.push_back(d);
							}
						}
					sort(dist.begin(), dist.end());
					phase[jj] = actx->markers[dist.begin()->second].phase;

					ierr = setDataPhaseDiagram(Pd, Pr[jj], T[jj], mat[PetscInt(phase[jj])].pdn); CHKERRQ(ierr);

					mf_ptr[jj]=Pd->mf;

				}
			}
			else
			{
				mf_ptr[jj]=0.0;
			}


			if((Active[jj] == 0.0) && actx->Ptr->Condition_pr != _Always_)
			{
				ierr = Check_advection_condition(actx, jj, ID,xp,yp,zp,Pr[jj],T[jj],melt_grid[jj]); CHKERRQ(ierr);
			}


			// override temperature of air phase
			if(AirPhase != -1 && phase[jj] == AirPhase) T[jj] = Ttop;

			// advect marker

			if( Active[jj]==1.0)
			{
                numActTracers += 1; // keep track of the # of active tracers on this processor
				npx = xp + vx*dt;
				npy = yp + vy*dt;
				npz = zp + vz*dt;
			}
			else
			{
				npx = xp;
				npy = yp;
				npz = zp;
			}

			if(npz > endz)
				{
					npz = zp;
					Active[jj]=0.0;
				}
			else if(npz < begz)
				{
					npz = zp;
					Active[jj]=0.0;
				}

			if(npy > endy)
				{
					npy = yp;
					Active[jj]=0.0;
				}
			else if(npy < begy)
				{
				  	npy = yp;
					Active[jj]=0.0;
				}


			if(npx > endx)
				{
					npx = xp;
					Active[jj]=0.0;

				}
			else if(npx < begx)
				{
					npx = xp;
					Active[jj]=0.0;
				}



			Xp[jj]=npx;
			Yp[jj]=npy;
			Zp[jj]=npz;
		}
		else
		{
			Xp[jj]      =   -DBL_MAX;
			Yp[jj]      =   -DBL_MAX;
			Zp[jj]      =   -DBL_MAX;
			Pr[jj]      =   -DBL_MAX;
			T[jj]       =   -DBL_MAX;
			phase[jj]   =   -DBL_MAX;
			mf_ptr[jj]  =   -DBL_MAX;
			Active[jj]  =   -DBL_MAX;
			melt_grid[jj] = -DBL_MAX;

		}

	}


	ierr = VecRestoreArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->p, &Pr)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->T, &T)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->phase, &phase)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->Melt_fr, &mf_ptr)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->Melt_Grid, &melt_grid)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->C_advection, &Active)           ; CHKERRQ(ierr);


	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	// get local grid sizes
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);


	if(ISParallel(PETSC_COMM_WORLD))
	{

		// sync pressure
		ierr = Sync_Vector(actx->Ptr->p,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->p); CHKERRQ(ierr);

		// sync temperature
		ierr = Sync_Vector(actx->Ptr->T,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->T); CHKERRQ(ierr);

		// sync coordinate
		//x
		ierr = Sync_Vector(actx->Ptr->x,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->x); CHKERRQ(ierr);

		//y
		ierr = Sync_Vector(actx->Ptr->y,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->y); CHKERRQ(ierr);

		//z
		ierr = Sync_Vector(actx->Ptr->z,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->z); CHKERRQ(ierr);

		// sync melt fraction
		// melt fraction of the particle
		ierr = Sync_Vector(actx->Ptr->Melt_fr,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->Melt_fr); CHKERRQ(ierr);

		// melt fraction of the grid
		ierr = Sync_Vector(actx->Ptr->Melt_Grid,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->Melt_Grid); CHKERRQ(ierr);

		//sync advection condition
		ierr = Sync_Vector(actx->Ptr->C_advection,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->C_advection); CHKERRQ(ierr);

		//sync phase
		ierr = Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->phase); CHKERRQ(ierr);

		// number of active tracer in the whole domain
        PetscInt numActTracers_0;
        ierr = MPI_Reduce(&numActTracers, &numActTracers_0, 1, MPIU_INT, MPI_SUM, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        numActTracers   = numActTracers_0;       // sum of # of active tracers on root

	}

    // print output
    PetscPrintf(PETSC_COMM_WORLD,"\n Currently active tracers    :  %i \n", numActTracers);


	// Check whatever the marker are belonging to rocks phase or not

	ierr = ADVMarkCrossFreeSurfPassive_Tracers(actx); CHKERRQ(ierr);



	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark); CHKERRQ(ierr);

		ierr = VecCopy(actx->Ptr->Recv,actx->Ptr->phase); CHKERRQ(ierr);
	}

	PrintDone(t);
	

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------//

#undef __FUNCT__
#define __FUNCT__ "ADVMarkCrossFreeSurfPassive_Tracers"
PetscErrorCode ADVMarkCrossFreeSurfPassive_Tracers(AdvCtx *actx)
{
	// change marker passive tracers when crossing free surface

	Marker           *IP;
	FDSTAG          *fs;
	FreeSurf        *surf;
	Vec             vphase;
	PetscInt        sx, sy, sz;
	PetscInt        ii, jj, ID, I, J, K, L, AirPhase, phaseID, nmark, *markind, markid;
	PetscScalar     ***ltopo, ***phase, *ncx, *ncy, topo, xp, yp, zp, *IX,bz,ez,by,ey,bx,ex,Xm[3];
	PetscScalar *Xp, *Yp,*Zp,*phaseptr;
	spair           d;
	vector <spair>  dist;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free-surface cases only
	if(!actx->surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	surf      = actx->surf;
	fs        = actx->fs;
	L         = fs->dsz.rank;
	AirPhase  = surf->AirPhase;

	// starting indices & number of cells
	sx = fs->dsx.pstart;
	sy = fs->dsy.pstart;
	sz = fs->dsz.pstart;

	// grid coordinates
	ncx = fs->dsx.ncoor;
	ncy = fs->dsy.ncoor;

	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);


	// reserve marker distance buffer
	dist.reserve(_mark_buff_sz_);

	// request local vector for reference sedimentation phases
	ierr = DMGetLocalVector(fs->DA_CEN, &vphase);

	// compute reference sedimentation phases
	ierr = ADVGetSedPhase(actx, vphase); CHKERRQ(ierr);

	// access topography & phases
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &ltopo);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN,    vphase,      &phase);  CHKERRQ(ierr);

	ierr = VecGetArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->phase, &phaseptr)           ; CHKERRQ(ierr);

	// scan all markers
	for(jj = 0; jj < actx->Ptr->nummark; jj++)
	{
		// access next marker
		xp = Xp[jj];
		yp = Yp[jj];
		zp = Zp[jj];;
		// get consecutive index of the host cell



		if(xp >= bx && xp<ex && yp >= by && yp< ey && zp >= bz && zp< ez)
		{

			ierr = Discret1DFindPoint(&fs->dsx, xp, I); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsy, yp, J); CHKERRQ(ierr);
			ierr = Discret1DFindPoint(&fs->dsz, zp, K); CHKERRQ(ierr);

			GET_CELL_ID(ID, I, J, K, fs->dsx.ncels, fs->dsy.ncels)


			// compute surface topography at marker position
			topo = InterpLin2D(ltopo, I, J, L, sx, sy, xp, yp, ncx, ncy);

			// check whether rock marker is above the free surface
			if(phaseptr[jj] != AirPhase && zp > topo)
			{
				// erosion (physical or numerical) -> rock turns into air
				phaseptr[jj]= AirPhase;
			}

			// check whether air marker is below the free surface
			if(phaseptr[jj] == AirPhase && zp < topo)
			{
				if(surf->SedimentModel > 0)
				{
				// sedimentation (physical) -> air turns into a prescribed rock
					phaseptr[jj]= surf->phase;
				}
				else
				{
				// sedimentation (numerical) -> air turns into closest (reference) rock
					Xm[0]=xp;
					Xm[1]=yp;
					Xm[2]=zp;

				// get marker list in containing cell
					nmark   = actx->markstart[ID+1] - actx->markstart[ID];
					markind = actx->markind + actx->markstart[ID];

				// clear distance storage
					dist.clear();

					for(ii = 0; ii < nmark; ii++)
					{
					// get current marker
						markid = markind[ii];
						IP     = &actx->markers[markid];

						// sort out air markers
						if(IP->phase == AirPhase) continue;

						// get marker coordinates
						IX = IP->X;

						// store marker index and distance
						d.first  = EDIST(Xm, IX);
						d.second = markid;

						dist.push_back(d);
					}

					// find closest rock marker (if any)
					if(dist.size())
					{
					// sort rock markers by distance
						sort(dist.begin(), dist.end());

					// copy phase from closest marker
						IP = &actx->markers[dist.begin()->second];

						phaseptr[jj] = IP->phase;
					}
					else
					{
					// no local rock marker found, set phase to reference
						phaseID = (PetscInt)phase[sz+K][sy+J][sx+I];

						if(phaseID < 0)
						{
						SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect sedimentation phase");
						}

						phaseptr[jj] = phaseID;
					}
				}
				//=======================================================================
				// WARNING! At best clone history from nearest rock marker
				//=======================================================================

			}
		}
		else
		{
			phaseptr[jj] = -DBL_MAX;
		}
	}


	ierr = VecRestoreArray(actx->Ptr->x, &Xp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->y, &Yp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->z, &Zp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->phase, &phaseptr) ; CHKERRQ(ierr);
	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &ltopo);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN,    vphase,      &phase);  CHKERRQ(ierr);

	// restore phase vector
	ierr = DMRestoreLocalVector(fs->DA_CEN, &vphase); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Check_advection_condition"
PetscErrorCode Check_advection_condition(AdvCtx *actx, PetscInt jj, PetscInt ID, PetscScalar xp, PetscScalar yp, PetscScalar zp, PetscScalar P,PetscScalar T,PetscScalar mf)
{

	PetscScalar 		*phase,*Active;
	PetscScalar 		Xm[3],X[3];
	vector <spair>    	dist;
	spair 				d;
	PetscErrorCode 		ierr;

	PetscFunctionBegin;

	ierr = VecGetArray(actx->Ptr->C_advection, &Active); CHKERRQ(ierr);

	if(actx->Ptr->Condition_pr == _Time_ptr_)
	{
		if((actx->jr->ts->time>=actx->Ptr->value_condition )&& Active[jj] == 0.0)
		{
			Active[jj] = 1.0;
		}
	}
	else if(actx->Ptr->Condition_pr ==_Melt_Fr_)
	{
		if((mf>=actx->Ptr->value_condition) && Active[jj] == 0.0)
		{
			Active[jj] = 1.0;
		}
	}
	else if(actx->Ptr->Condition_pr ==_Temp_ptr_)
	{
		if((T>=actx->Ptr->value_condition) && Active[jj] == 0.0)
		{
			Active[jj] = 1.0;
		}
	}
	else if(actx->Ptr->Condition_pr ==_Pres_ptr_)
	{
		if((P>=actx->Ptr->value_condition) && Active[jj] == 0.0)
		{
			Active[jj] = 1.0;
		}
	}

	// overwrite the phase in case of delayed activation or if some condition are met

	if(((actx->Ptr->Condition_pr ==_Pres_ptr_)||(actx->Ptr->Condition_pr ==_Temp_ptr_)||(actx->Ptr->Condition_pr ==_Time_ptr_)) && Active[jj] == 1.0)
	{

		PetscInt n, ii;

		ierr = VecGetArray(actx->Ptr->phase, &phase); CHKERRQ(ierr);


		X[0]  = xp;
		X[1]  = yp;
		X[2]  = zp;

		dist.clear();
		n = actx->markstart[ID+1] - actx->markstart[ID];

		for (ii = actx->markstart[ID]; ii <actx->markstart[ID]+n; ii++)
			{
				Xm[0] =actx->markers[ii].X[0];
				Xm[1] =actx->markers[ii].X[1];
				Xm[2] =actx->markers[ii].X[2];


				d.first  = EDIST(Xm, X);
				d.second = ii;
				dist.push_back(d);

			}
		sort(dist.begin(), dist.end());
		phase[jj]= actx->markers[dist.begin()->second].phase;
		ierr = VecRestoreArray(actx->Ptr->phase, &phase); CHKERRQ(ierr);



	}

	ierr = VecRestoreArray(actx->Ptr->C_advection, &Active); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}



//----------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "ADVPtrDestroy"
PetscErrorCode ADVPtrDestroy(AdvCtx *actx)
{

	PetscFunctionBegin;


	// check whether current storage is insufficient

	VecDestroy(&actx->Ptr->ID);

	VecDestroy(&actx->Ptr->x);

	VecDestroy(&actx->Ptr->y);

	VecDestroy(&actx->Ptr->z);

	VecDestroy(&actx->Ptr->T);

	VecDestroy(&actx->Ptr->p);

	VecDestroy(&actx->Ptr->phase);

	VecDestroy(&actx->Ptr->Melt_fr);

	VecDestroy(&actx->Ptr->Melt_Grid);

	VecDestroy(&actx->Ptr->C_advection);

	VecDestroy(&actx->Ptr->Recv);



	PetscFunctionReturn(0);
}
// --------------------------------------------------------------------------------------- //

//-------------------------------------------------------------------------//

#undef __FUNCT__
#define __FUNCT__ "Passive_Tracer_WriteRestart"
PetscErrorCode Passive_Tracer_WriteRestart(AdvCtx *actx, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->jr->ctrl.Passive_Tracer)
	{

	// write solution vectors
	ierr = VecWriteRestart(actx->Ptr->x, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->y, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->z, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->p, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->T, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->phase, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->Melt_fr, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->Melt_Grid, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->C_advection, fp); CHKERRQ(ierr);
	ierr = VecWriteRestart(actx->Ptr->ID, fp); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

// --------------------------------------------------------------------------------------- //

#undef __FUNCT__
#define __FUNCT__ "ReadPassive_Tracers"
PetscErrorCode ReadPassive_Tracers(AdvCtx *actx, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read solution vectors
	if(actx->jr->ctrl.Passive_Tracer)
	{
		ierr = ADVPtrReCreateStorage(actx); CHKERRQ(ierr);

		ierr = VecReadRestart(actx->Ptr->x, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->y, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->z, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->p, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->T, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->phase, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->Melt_fr, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->Melt_Grid, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->C_advection, fp); CHKERRQ(ierr);
		ierr = VecReadRestart(actx->Ptr->ID, fp); CHKERRQ(ierr);
	}



	PetscFunctionReturn(0);
}
//---------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Sync_Vector"
PetscErrorCode Sync_Vector(Vec x,AdvCtx *actx ,PetscInt nummark)
{
	PetscScalar *recv,*send;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecZeroEntries(actx->Ptr->Recv); CHKERRQ(ierr);

	ierr = VecGetArray(x, &send)           ; CHKERRQ(ierr);
	ierr = VecGetArray(actx->Ptr->Recv, &recv)           ; CHKERRQ(ierr);

	ierr = MPI_Allreduce(send, recv, (PetscMPIInt)nummark, MPIU_SCALAR, MPI_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);

	ierr = VecRestoreArray(x, &send)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(actx->Ptr->Recv, &recv)           ; CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


//=========================================================
/*
#undef __FUNCT__
#define __FUNCT__ "Passive_tracers_save"
PetscErrorCode Passive_tracers_save(AdvCtx *actx)
{	// save new restart database, then delete the original

	Scaling        *scal;
	FILE           *fp;
	char           *fileName;
	PetscScalar    time;
	PetscInt        step;
	PetscInt       ii;
	PetscScalar    *xp,*yp,*zp,*P,*T,*phase,*ID,*mf_ptr,*Active;
	PetscErrorCode ierr;
	PetscFunctionBegin;


	if(actx->jr->ctrl.Passive_Tracer == 0) PetscFunctionReturn(0);

	scal = actx->jr->scal;
	step    = actx->jr->ts->istep;
	time    = actx->jr->ts->time*actx->jr->scal->time;

	if(step==0)
	{
		ierr = DirMake("./Passive_Tracers"); CHKERRQ(ierr);
	}

	if(ISRankZero(PETSC_COMM_WORLD))
	{
	// compile actual & temporary restart file name
		asprintf(&fileName, "./Passive_Tracers/PT_%1.8lld.dat",(LLD)step);

	// open temporary restart file for writing in binary mode
		fp = fopen(fileName, "w");
		fprintf(fp,"number_marker = %d \n ", actx->Ptr->nummark);

		fprintf(fp,"\n");


		fprintf(fp,"Time = %6f Timestep = %d \n",time,step);

		fprintf(fp,"nx = %d ny = %d nz = %d \n", actx->Ptr->passive_tracer_resolution[0],actx->Ptr->passive_tracer_resolution[1],actx->Ptr->passive_tracer_resolution[2]);

		fprintf(fp,"\n");

		fprintf(fp," # ID  X  Y  Z  P  T  PH MeltFr\r\n");

		ierr = VecGetArray(actx->Ptr->ID, &ID)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->x, &xp)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->y, &yp)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->z, &zp)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->p, &P)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->T, &T)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->phase, &phase)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->Melt_fr, &mf_ptr)           ; CHKERRQ(ierr);
		ierr = VecGetArray(actx->Ptr->C_advection, &Active)           ; CHKERRQ(ierr);


		for(ii=0;ii<actx->Ptr->nummark;ii++)
		{


			fprintf(fp," %d %3f   %3f  %3f  %2f  %2f  %d %6f %d \r\n",PetscInt(ID[ii]),xp[ii]*scal->length,yp[ii]*scal->length,zp[ii]*scal->length,P[ii]*scal->stress,T[ii]*scal->temperature - scal->Tshift, PetscInt(phase[ii]),mf_ptr[ii],PetscInt(Active[ii]));

		}
		ierr = VecRestoreArray(actx->Ptr->ID, &ID)         ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->x, &xp)          ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->y, &yp)          ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->z, &zp)          ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->p, &P)           ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->T, &T)           ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->phase, &phase)   ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->Melt_fr, &mf_ptr)           ; CHKERRQ(ierr);
		ierr = VecRestoreArray(actx->Ptr->C_advection, &Active)           ; CHKERRQ(ierr);


		fclose(fp);

		free(fileName);
	}

	PetscFunctionReturn(0);

}
*/
