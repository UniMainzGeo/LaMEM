/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

/*
 *  Created on: Jul 28, 2020
 *      Author: piccolo
 */

//---------------------------------------------------------------------------

#include "LaMEM.h"
#include "AVD.h"
#include "passive_tracer.h"
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
#include "interpolate.h"

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

PetscErrorCode ADVPtrPassive_Tracer_create(AdvCtx *actx, FB *fb)
{
/*
 *  This function creates all the vector required for tracing pressure, temperature, phase and x,y,z position.
 *  RecvBuf is a vector used only for the synching operation and it has any meaning.
 */

	P_Tr            *passive_tr;
	char             Condition_adv[_str_len_];
	PetscInt        nummark;
	PetscFunctionBeginUser;

	if(!actx->jr->ctrl.Passive_Tracer)	PetscFunctionReturn(0);
	passive_tr = actx->Ptr;


    // Read input parameters
	PetscCall(getScalarParam(fb, _REQUIRED_, "PassiveTracer_Box",          passive_tr->box_passive_tracer,         6,  1.0));
	PetscCall(getIntParam   (fb, _REQUIRED_, "PassiveTracer_Resolution",   passive_tr->passive_tracer_resolution,  3,  0));
	PetscCall(getStringParam(fb, _OPTIONAL_, "PassiveTracer_ActiveType",   Condition_adv,      "Always"));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "PassiveTracer_ActiveValue",  &passive_tr->value_condition,           1, 1.0));

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
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "The total number of passive tracers must be lower than %" PetscInt_FMT "",_max_passive_tracer);
	}


     PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
     PetscPrintf(PETSC_COMM_WORLD,"Passive Tracers: \n");
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box x = [Left,Right] : %6f, %6f \n",passive_tr->box_passive_tracer[0],passive_tr->box_passive_tracer[1]);
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box y = [Front,Back] : %6f, %6f \n",passive_tr->box_passive_tracer[2],passive_tr->box_passive_tracer[3]);
     PetscPrintf(PETSC_COMM_WORLD,"   Initial coordinate Box z = [Bot, Top]   : %6f, %6f \n",passive_tr->box_passive_tracer[4],passive_tr->box_passive_tracer[5]);
     PetscPrintf(PETSC_COMM_WORLD,"   # of tracers in [x,y,z] direction       : [%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT "] \n", passive_tr->passive_tracer_resolution[0],  passive_tr->passive_tracer_resolution[1],  passive_tr->passive_tracer_resolution[2]);
     PetscPrintf(PETSC_COMM_WORLD,"   Total # of tracers                      : %" PetscInt_FMT " \n", nummark);
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
	 PetscCall(ADVPtrReCreateStorage(actx));

	 // Initialize the initial coordinate distribution and phase
	 PetscCall(ADVPassiveTracerInit(actx));

	 PetscFunctionReturn(0);
	}
// ---------------------------------------------------------------------------------------------------------------------------//
PetscErrorCode ADVPtrReCreateStorage(AdvCtx *actx)
{

	PetscFunctionBeginUser;

		if(!actx->jr->ctrl.Passive_Tracer)	PetscFunctionReturn(0);
	// check whether current storage is insufficient

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->ID));
		PetscCall(VecZeroEntries(actx->Ptr->ID));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->x));
		PetscCall(VecZeroEntries(actx->Ptr->x));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark, &actx->Ptr->y));
		PetscCall(VecZeroEntries(actx->Ptr->y));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark, &actx->Ptr->z));
		PetscCall(VecZeroEntries(actx->Ptr->z));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark  ,&actx->Ptr->T));
		PetscCall(VecZeroEntries(actx->Ptr->T));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->p));
		PetscCall(VecZeroEntries(actx->Ptr->p));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->phase));
		PetscCall(VecZeroEntries(actx->Ptr->phase));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Melt_fr));
		PetscCall(VecZeroEntries(actx->Ptr->Melt_fr));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->C_advection));
		PetscCall(VecZeroEntries(actx->Ptr->C_advection));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Recv));
		PetscCall(VecZeroEntries(actx->Ptr->Recv));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->Melt_Grid));
		PetscCall(VecZeroEntries(actx->Ptr->Melt_Grid));

		PetscCall(VecCreateSeq(PETSC_COMM_SELF,actx->Ptr->nummark ,&actx->Ptr->APS));
		PetscCall(VecZeroEntries(actx->Ptr->APS));

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------//

PetscErrorCode ADVPassiveTracerInit(AdvCtx *actx)
{
	
	PetscFunctionBeginUser;

	if(actx->jr->ctrl.Passive_Tracer == 0 ) PetscFunctionReturn(0);

	PetscCall(ADVPtrInitCoord(actx));

	PetscCall(ADV_Assign_Phase(actx));


	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode ADVPtrInitCoord(AdvCtx *actx)
{
	//Initialize the passive tracer lagrangian grid. The initial passive tracer distribution is a rectangular grid, with a
	// a variable resolution. After initializing the coordinates, phase, temperature and pressure are interpolated from
	// the nearest marker (s.s.)

	PetscScalar  x, y, z, dx, dy, dz,nx,ny,nz;
	PetscInt     i, j, k;
	PetscInt     imark;
	PetscScalar  *Xp,*Yp,*Zp,*ID,*active;

	
	PetscFunctionBeginUser;

	nx = (PetscScalar) actx->Ptr->passive_tracer_resolution[0];
	ny = (PetscScalar) actx->Ptr->passive_tracer_resolution[1];
	nz = (PetscScalar) actx->Ptr->passive_tracer_resolution[2];
	dx = (actx->Ptr->box_passive_tracer[1]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length))/nx;
	dy = (actx->Ptr->box_passive_tracer[3]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length))/ny;
	dz = (actx->Ptr->box_passive_tracer[5]/(actx->dbm->scal->length)-actx->Ptr->box_passive_tracer[4]/(actx->dbm->scal->length))/nz;


	// marker counter
	imark = 0;
	PetscCall(VecGetArray(actx->Ptr->x, &Xp));
	PetscCall(VecGetArray(actx->Ptr->y, &Yp));
	PetscCall(VecGetArray(actx->Ptr->z, &Zp));
	PetscCall(VecGetArray(actx->Ptr->ID, &ID));
	PetscCall(VecGetArray(actx->Ptr->C_advection, &active));



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
					z = actx->Ptr->box_passive_tracer[4]/(actx->dbm->scal->length) + dz/2 + ((PetscScalar) k)*dz;
				}
				if(j==0)
				{
					y = actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length) + dy/2;
				}
				else
				{
					y = actx->Ptr->box_passive_tracer[2]/(actx->dbm->scal->length) + dy/2 + ((PetscScalar) j)*dy;
				}
				if(i==0)
				{
					x = actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length) + dx/2;
				}
				else
				{
					x = actx->Ptr->box_passive_tracer[0]/(actx->dbm->scal->length) + dx/2+ ((PetscScalar) i)*dx;
				}


				// set marker coordinates
				Xp[imark] = x;
				Yp[imark] = y;
				Zp[imark] = z;
				ID[imark] = ((PetscScalar) i) + ny*((PetscScalar) j) +ny*nx*((PetscScalar) k);

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


	PetscCall(VecRestoreArray(actx->Ptr->x, &Xp));
	PetscCall(VecRestoreArray(actx->Ptr->y, &Yp));
	PetscCall(VecRestoreArray(actx->Ptr->z, &Zp));
	PetscCall(VecRestoreArray(actx->Ptr->ID, &ID));
	PetscCall(VecRestoreArray(actx->Ptr->C_advection, &active));



	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode ADV_Assign_Phase(AdvCtx *actx)
{
	// Initially the marker are phase-less. This routine assign both phase and
	// initial temperature to the passive tracers.

	FDSTAG      *fs;
	vector <spair>    dist;
	spair d;
	Marker   *IP;
	PetscScalar  X[3],Xm[3],*Xp,*Yp,*Zp,*Pr,*T,*phase,*APS;
	PetscInt     I, J, K,ii,numpassive,imark,ID,nx,ny,n,*markind,id_m;
	PetscScalar ex,bx,ey,by,ez,bz;


	
	PetscFunctionBeginUser;


	fs = actx->fs;

	numpassive = actx->Ptr->nummark;

	PetscCall(ADVMapMarkToCells(actx));

	// get context
	fs = actx->fs;

	// starting indices & number of cells
	nx = fs->dsx.ncels;
	ny = fs->dsy.ncels;

	PetscCall(FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez));


	dist.reserve(_mark_buff_sz_);
	PetscCall(VecGetArray(actx->Ptr->x, &Xp));
	PetscCall(VecGetArray(actx->Ptr->y, &Yp));
	PetscCall(VecGetArray(actx->Ptr->z, &Zp));
	PetscCall(VecGetArray(actx->Ptr->p, &Pr));
	PetscCall(VecGetArray(actx->Ptr->T, &T));
	PetscCall(VecGetArray(actx->Ptr->APS, &APS));
	PetscCall(VecGetArray(actx->Ptr->phase, &phase));

	for(imark=0;imark<numpassive;imark++)
	{
		// get marker coordinates
			X[0] = Xp[imark];
			X[1] = Yp[imark];
			X[2] = Zp[imark];



			if(X[0] >= bx && X[0]< ex && X[1] >= by && X[1]< ey && X[2] >= bz && X[2]< ez)
			{
			// get host cell IDs in all directions
			PetscCall(Discret1DFindPoint(&fs->dsx, X[0], I));
			PetscCall(Discret1DFindPoint(&fs->dsy, X[1], J));
			PetscCall(Discret1DFindPoint(&fs->dsz, X[2], K));

			// compute and store consecutive index
			GET_CELL_ID(ID, I, J, K, nx, ny);

			dist.clear();


			n = actx->markstart[ID+1] - actx->markstart[ID];
			markind = actx->markind + actx->markstart[ID];

			for (ii = 0; ii < n; ii++)
			{
				id_m=markind[ii];
				Xm[0] = actx->markers[id_m].X[0];
				Xm[1] = actx->markers[id_m].X[1];
				Xm[2] = actx->markers[id_m].X[2];


				d.first  = EDIST(X, Xm);
				d.second = id_m;
				dist.push_back(d);
			}

			// sort markers by distance
			sort(dist.begin(), dist.end());
			IP = &actx->markers[dist.begin()->second];

			// clone closest marker
			phase[imark]= ((PetscScalar) IP->phase);
			T[imark]= IP->T;
			Pr[imark]= IP->p;
			APS[imark]= IP->APS;
			}
			else
			{

				phase[imark]= -DBL_MAX;
				T[imark]= -DBL_MAX;
				Pr[imark]= -DBL_MAX;
				APS[imark]= -DBL_MAX;
			}

	}

	PetscCall(VecRestoreArray(actx->Ptr->x, &Xp));
	PetscCall(VecRestoreArray(actx->Ptr->y, &Yp));
	PetscCall(VecRestoreArray(actx->Ptr->z, &Zp));
	PetscCall(VecRestoreArray(actx->Ptr->p, &Pr));
	PetscCall(VecRestoreArray(actx->Ptr->T, &T));
	PetscCall(VecRestoreArray(actx->Ptr->APS, &APS));
	PetscCall(VecRestoreArray(actx->Ptr->phase, &phase));

	if(ISParallel(PETSC_COMM_WORLD))
	{
		PetscCall(Sync_Vector(actx->Ptr->p,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->p));

		PetscCall(Sync_Vector(actx->Ptr->T,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->T));

		PetscCall(Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->phase));

		PetscCall(Sync_Vector(actx->Ptr->APS,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->APS));
	}


	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------
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
	PetscInt        sx, sy, sz, nx, ny,nz;
	PetscInt        jj, I, J, K, II, JJ, KK, AirPhase, num_part,ID, n, ii, numActTracers,*markind,id_m ;
	PetscScalar     ex,bx,ey,by,ez,bz;
	PetscScalar     *ncx, *ncy, *ncz;
	PetscScalar     *ccx, *ccy, *ccz;
	PetscScalar     ***lvx, ***lvy, ***lvz, ***lp, ***lT;
	PetscScalar     vx, vy, vz, xc, yc, zc, xp, yp, zp, dt, Ttop, endx,endy,endz,begx,begy,begz,npx,npy,npz;
	PetscScalar     *Xp, *Yp,*Zp,*T,*Pr,*phase,*mf_ptr,*Active,*melt_grid,*aps;
	PetscScalar     pShift;
	PetscScalar     Xm[3],X[3];
	PetscLogDouble t;
	PetscMPIInt 	rank;
	vector <spair>    dist;
	spair d;
	
	PetscFunctionBeginUser;

	AirPhase = -1;
	Ttop     =  0.0;

	PetscCallMPI(MPI_Comm_rank( MPI_COMM_WORLD, &rank));

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

	// initialize corners and edges for interpolation
	PetscCall(SetEdgeCornerXFace (fs, jr->lvx));
	PetscCall(SetEdgeCornerYFace (fs, jr->lvy));
	PetscCall(SetEdgeCornerZFace (fs, jr->lvz));
	PetscCall(SetEdgeCornerCenter(fs, jr->lp));
	PetscCall(SetEdgeCornerCenter(fs, jr->lT));

	// access velocity, pressure & temperature vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp,  &lp));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT));

	PetscCall(FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez));


	PetscCall(VecGetArray(actx->Ptr->x, &Xp));
	PetscCall(VecGetArray(actx->Ptr->y, &Yp));
	PetscCall(VecGetArray(actx->Ptr->z, &Zp));
	PetscCall(VecGetArray(actx->Ptr->p, &Pr));
	PetscCall(VecGetArray(actx->Ptr->T, &T));
	PetscCall(VecGetArray(actx->Ptr->phase, &phase));
	PetscCall(VecGetArray(actx->Ptr->Melt_fr, &mf_ptr));
	PetscCall(VecGetArray(actx->Ptr->Melt_Grid, &melt_grid));
	PetscCall(VecGetArray(actx->Ptr->APS, &aps));
	PetscCall(VecGetArray(actx->Ptr->C_advection, &Active));

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
			PetscCall(Discret1DFindPoint(&fs->dsx, xp, I));
			PetscCall(Discret1DFindPoint(&fs->dsy, yp, J));
			PetscCall(Discret1DFindPoint(&fs->dsz, zp, K));

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


			melt_grid[jj] = svCell->svBulk.mf;

			if(svCell->svBulk.mf>0.0)
			{
			  //check if the original phase saved is one that has a phase/melt law associated


				if(mat[PetscInt(phase[jj])].pdn[0] != '\0')
				{
					PetscCall(setDataPhaseDiagram(Pd, Pr[jj], T[jj], mat[PetscInt(phase[jj])].pdn));
					mf_ptr[jj]= Pd->mf;
				}
				else
				{
					// Passive tracers are initialize during the initial stage of the simulation.
					// They can have a different phase as soon as the melting start.

					// sort markers by distance
					dist.clear();
					n = actx->markstart[ID+1] - actx->markstart[ID];
					markind = actx->markind + actx->markstart[ID];


					for (ii = 0; ii < n; ii++)
					{
						id_m=markind[ii];
						Xm[0] = actx->markers[id_m].X[0];
						Xm[1] = actx->markers[id_m].X[1];
						Xm[2] = actx->markers[id_m].X[2];
						X[0]  = xp;
						X[1]  = yp;
						X[2]  = zp;

						if(mat[actx->markers[ii].phase].pdn[0] != '\0')
						{
							d.first  = EDIST(Xm, X);
							d.second = id_m;
							dist.push_back(d);
						}
					}
					sort(dist.begin(), dist.end());
					phase[jj] = (PetscScalar) actx->markers[dist.begin()->second].phase;

					PetscCall(setDataPhaseDiagram(Pd, Pr[jj], T[jj], mat[PetscInt(phase[jj])].pdn));

					mf_ptr[jj]=Pd->mf;

				}
			}
			else
			{
				mf_ptr[jj]=0.0;
			}


			if((Active[jj] == 0.0) && actx->Ptr->Condition_pr != _Always_)
			{
				PetscCall(Check_advection_condition(actx, jj, ID,xp,yp,zp,Pr[jj],T[jj],melt_grid[jj]));
			}

			// override temperature of air phase
			if(AirPhase != -1 && phase[jj] == ((PetscScalar) AirPhase)) T[jj] = Ttop;

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


	PetscCall(VecRestoreArray(actx->Ptr->x, &Xp));
	PetscCall(VecRestoreArray(actx->Ptr->y, &Yp));
	PetscCall(VecRestoreArray(actx->Ptr->z, &Zp));
	PetscCall(VecRestoreArray(actx->Ptr->p, &Pr));
	PetscCall(VecRestoreArray(actx->Ptr->T, &T));
	PetscCall(VecRestoreArray(actx->Ptr->phase, &phase));
	PetscCall(VecRestoreArray(actx->Ptr->Melt_fr, &mf_ptr));
	PetscCall(VecRestoreArray(actx->Ptr->Melt_Grid, &melt_grid));
	PetscCall(VecRestoreArray(actx->Ptr->APS, &aps));
	PetscCall(VecRestoreArray(actx->Ptr->C_advection, &Active));


	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT));

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));


	if(ISParallel(PETSC_COMM_WORLD))
	{

		// sync pressure
		PetscCall(Sync_Vector(actx->Ptr->p,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->p));

		// sync temperature
		PetscCall(Sync_Vector(actx->Ptr->T,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->T));

		// sync coordinate
		//x
		PetscCall(Sync_Vector(actx->Ptr->x,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->x));

		//y
		PetscCall(Sync_Vector(actx->Ptr->y,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->y));

		//z
		PetscCall(Sync_Vector(actx->Ptr->z,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->z));

		// sync melt fraction
		// melt fraction of the particle
		PetscCall(Sync_Vector(actx->Ptr->Melt_fr,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->Melt_fr));

		// melt fraction of the grid
		PetscCall(Sync_Vector(actx->Ptr->Melt_Grid,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->Melt_Grid));

		//sync advection condition
		PetscCall(Sync_Vector(actx->Ptr->C_advection,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->C_advection));

		//sync phase
		PetscCall(Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->phase));

		//sync APS
		PetscCall(Sync_Vector(actx->Ptr->APS,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->APS));

		// number of active tracer in the whole domain
        PetscInt numActTracers_0;
        PetscCallMPI(MPI_Reduce(&numActTracers, &numActTracers_0, 1, MPIU_INT, MPI_SUM, 0, PETSC_COMM_WORLD));
        numActTracers   = numActTracers_0;       // sum of # of active tracers on root

	}

    // print output
    PetscPrintf(PETSC_COMM_WORLD,"\n Currently active tracers    :  %" PetscInt_FMT " \n",  numActTracers);


	// Check whatever the marker are belonging to rocks phase or not

	PetscCall(ADVMarkCrossFreeSurfPassive_Tracers(actx));



	if(ISParallel(PETSC_COMM_WORLD))
	{
		PetscCall(Sync_Vector(actx->Ptr->phase,actx,actx->Ptr->nummark));

		PetscCall(VecCopy(actx->Ptr->Recv,actx->Ptr->phase));
	}

	PrintDone(t);
	

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------//

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


	
	PetscFunctionBeginUser;

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

	PetscCall(FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez));


	// reserve marker distance buffer
	dist.reserve(_mark_buff_sz_);

	// request local vector for reference sedimentation phases
	PetscCall(DMGetLocalVector(fs->DA_CEN, &vphase));

	// compute reference sedimentation phases
	PetscCall(ADVGetSedPhase(actx, vphase));

	// access topography & phases
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &ltopo));
	PetscCall(DMDAVecGetArray(fs->DA_CEN,    vphase,      &phase));

	PetscCall(VecGetArray(actx->Ptr->x, &Xp));
	PetscCall(VecGetArray(actx->Ptr->y, &Yp));
	PetscCall(VecGetArray(actx->Ptr->z, &Zp));
	PetscCall(VecGetArray(actx->Ptr->phase, &phaseptr));

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

			PetscCall(Discret1DFindPoint(&fs->dsx, xp, I));
			PetscCall(Discret1DFindPoint(&fs->dsy, yp, J));
			PetscCall(Discret1DFindPoint(&fs->dsz, zp, K));

			GET_CELL_ID(ID, I, J, K, fs->dsx.ncels, fs->dsy.ncels)


			// compute surface topography at marker position
			topo = InterpLin2D(ltopo, I, J, L, sx, sy, xp, yp, ncx, ncy);

			// check whether rock marker is above the free surface
            if(phaseptr[jj] != ((PetscScalar) AirPhase) && zp > topo)
			{
				// erosion (physical or numerical) -> rock turns into air
				phaseptr[jj]= ((PetscScalar) AirPhase);
			}

			// check whether air marker is below the free surface
			if(phaseptr[jj] == ((PetscScalar) AirPhase) && zp < topo)
			{
				if(surf->SedimentModel > 0)
				{
				// sedimentation (physical) -> air turns into a prescribed rock
					phaseptr[jj]= (PetscScalar) surf->phase;
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

						phaseptr[jj] = (PetscScalar) IP->phase;
					}
					else
					{
					// no local rock marker found, set phase to reference
						phaseID = (PetscInt)phase[sz+K][sy+J][sx+I];

						if(phaseID < 0)
						{
						SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect sedimentation phase");
						}

						phaseptr[jj] = (PetscScalar) phaseID;
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


	PetscCall(VecRestoreArray(actx->Ptr->x, &Xp));
	PetscCall(VecRestoreArray(actx->Ptr->y, &Yp));
	PetscCall(VecRestoreArray(actx->Ptr->z, &Zp));
	PetscCall(VecRestoreArray(actx->Ptr->phase, &phaseptr));
	// restore access
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &ltopo));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN,    vphase,      &phase));

	// restore phase vector
	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &vphase));

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------//
PetscErrorCode Check_advection_condition(AdvCtx *actx, PetscInt jj, PetscInt ID, PetscScalar xp, PetscScalar yp, PetscScalar zp, PetscScalar P,PetscScalar T,PetscScalar mf)
{

	PetscScalar 		*phase,*Active;
	PetscScalar 		Xm[3],X[3];
	vector <spair>    	dist;
	spair 				d;

	PetscFunctionBeginUser;

	PetscCall(VecGetArray(actx->Ptr->C_advection, &Active));

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

		PetscInt n, ii,id_m,*markind;

		PetscCall(VecGetArray(actx->Ptr->phase, &phase));


		X[0]  = xp;
		X[1]  = yp;
		X[2]  = zp;

		dist.clear();
		n = actx->markstart[ID+1] - actx->markstart[ID];
		markind = actx->markind + actx->markstart[ID];


		for (ii = 0; ii < n; ii++)
		{
			id_m=markind[ii];
			Xm[0] =actx->markers[id_m].X[0];
			Xm[1] =actx->markers[id_m].X[1];
			Xm[2] =actx->markers[id_m].X[2];


			d.first  = EDIST(Xm, X);
			d.second = id_m;
			dist.push_back(d);

		}
		sort(dist.begin(), dist.end());
		phase[jj]= (PetscScalar) actx->markers[dist.begin()->second].phase;
		PetscCall(VecRestoreArray(actx->Ptr->phase, &phase));



	}

	PetscCall(VecRestoreArray(actx->Ptr->C_advection, &Active));

	PetscFunctionReturn(0);
}



//----------------------------------------------------------------------------//
PetscErrorCode ADVPtrDestroy(AdvCtx *actx)
{

	PetscFunctionBeginUser;


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

	VecDestroy(&actx->Ptr->APS);

	VecDestroy(&actx->Ptr->C_advection);

	VecDestroy(&actx->Ptr->Recv);



	PetscFunctionReturn(0);
}
// --------------------------------------------------------------------------------------- //

//-------------------------------------------------------------------------//

PetscErrorCode Passive_Tracer_WriteRestart(AdvCtx *actx, FILE *fp)
{
	
	PetscFunctionBeginUser;

	if(actx->jr->ctrl.Passive_Tracer)
	{

	// write solution vectors
	PetscCall(VecWriteRestart(actx->Ptr->x, fp));
	PetscCall(VecWriteRestart(actx->Ptr->y, fp));
	PetscCall(VecWriteRestart(actx->Ptr->z, fp));
	PetscCall(VecWriteRestart(actx->Ptr->p, fp));
	PetscCall(VecWriteRestart(actx->Ptr->T, fp));
	PetscCall(VecWriteRestart(actx->Ptr->phase, fp));
	PetscCall(VecWriteRestart(actx->Ptr->Melt_fr, fp));
	PetscCall(VecWriteRestart(actx->Ptr->Melt_Grid, fp));
	PetscCall(VecWriteRestart(actx->Ptr->APS, fp));
	PetscCall(VecWriteRestart(actx->Ptr->C_advection, fp));
	PetscCall(VecWriteRestart(actx->Ptr->ID, fp));
	}

	PetscFunctionReturn(0);
}

// --------------------------------------------------------------------------------------- //

PetscErrorCode ReadPassive_Tracers(AdvCtx *actx, FILE *fp)
{
	
	PetscFunctionBeginUser;

	// read solution vectors
	if(actx->jr->ctrl.Passive_Tracer)
	{
		PetscCall(ADVPtrReCreateStorage(actx));

		PetscCall(VecReadRestart(actx->Ptr->x, fp));
		PetscCall(VecReadRestart(actx->Ptr->y, fp));
		PetscCall(VecReadRestart(actx->Ptr->z, fp));
		PetscCall(VecReadRestart(actx->Ptr->p, fp));
		PetscCall(VecReadRestart(actx->Ptr->T, fp));
		PetscCall(VecReadRestart(actx->Ptr->phase, fp));
		PetscCall(VecReadRestart(actx->Ptr->Melt_fr, fp));
		PetscCall(VecReadRestart(actx->Ptr->Melt_Grid, fp));
		PetscCall(VecReadRestart(actx->Ptr->APS, fp));
		PetscCall(VecReadRestart(actx->Ptr->C_advection, fp));
		PetscCall(VecReadRestart(actx->Ptr->ID, fp));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------
PetscErrorCode Sync_Vector(Vec x,AdvCtx *actx ,PetscInt nummark)
{
	PetscScalar *recv,*send;

	
	PetscFunctionBeginUser;

	PetscCall(VecZeroEntries(actx->Ptr->Recv));

	PetscCall(VecGetArray(x, &send));
	PetscCall(VecGetArray(actx->Ptr->Recv, &recv));

	PetscCallMPI(MPI_Allreduce(send, recv, (PetscMPIInt)nummark, MPIU_SCALAR, MPI_MAX,PETSC_COMM_WORLD));

	PetscCall(VecRestoreArray(x, &send));
	PetscCall(VecRestoreArray(actx->Ptr->Recv, &recv));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------
