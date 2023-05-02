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
 **    filename:   cvi.c
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
// Extended Advection and Conservative Velocity Interpolation Routines (CVI)
//---------------------------------------------------------------------------
// extended advection schemes for Runge-Kutta 2nd order and Runge-Kutta 4th order
// also testing velocity interpolation routines
// by A. Pusok, July 2015


#include "LaMEM.h"
#include "cvi.h"
#include "JacRes.h"
#include "fdstag.h"
#include "advect.h"
#include "surf.h"
#include "bc.h"
#include "tssolve.h"
#include "tools.h"

//---------------------------------------------------------------------------
PetscErrorCode ADVelAdvectMain(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION ROUTINE
	//=======================================================================
	AdvVelCtx vi;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// interpolate P,T - needs update
	ierr = ADVelInterpPT(actx); CHKERRQ(ierr);

	// velocity advection routine - with different velocity interpolations
	ierr = ADVelAdvectScheme(actx, &vi); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelInterpPT(AdvCtx *actx)
{
	// update p, T at marker positions

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt    nx, ny, sx, sy, sz;
	PetscInt    jj, ID, I, J, K, AirPhase;
	PetscScalar ***lp, ***lT, Ttop;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs = actx->fs;
	jr = actx->jr;

	AirPhase = -1;
	Ttop     =  0.0;

	if(actx->surf->UseFreeSurf)
	{
		AirPhase = actx->surf->AirPhase;
		Ttop     = actx->jr->bc->Ttop;
	}

	// starting indices & number of cells
	nx = fs->dsx.ncels;
	ny = fs->dsy.ncels;
	sx = fs->dsx.pstart;
	sy = fs->dsy.pstart;
	sz = fs->dsz.pstart;

	// access velocity, pressure & temperature vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	// scan all markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// access host cell solution variables
		svCell = &jr->svCell[ID];

		// update pressure & temperature variables
		P->p += lp[sz+K][sy+J][sx+I] - svCell->svBulk.pn;
		P->T += lT[sz+K][sy+J][sx+I] - svCell->svBulk.Tn;

		// override temperature of air phase
		if(AirPhase != -1 &&  P->phase == AirPhase) P->T = Ttop;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelAdvectScheme(AdvCtx *actx, AdvVelCtx *vi)
{
	PetscScalar  dt;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	//=======================================================================
	// VELOCITY ADVECTION ROUTINE
	//=======================================================================
	// create context
	ierr = ADVelCreate(actx, vi); CHKERRQ(ierr);

	// initialize marker positions
	ierr = ADVelInitCoord(actx, vi->interp, vi->nmark); CHKERRQ(ierr);

	// get current time step
	dt = actx->jr->ts->dt;

	//=======================================================================
	// START ADVECTION
	//=======================================================================
	// ---------------------------------
	// EULER (1st order)
	// ---------------------------------
	if(actx->advect == EULER)
	{
		// 1. Velocity interpolation
		ierr = ADVelInterpMain(vi); CHKERRQ(ierr);

		// 2. Update effective velocity
		ierr = ADVelCalcEffVel(vi->interp, vi->nmark, 1.0); CHKERRQ(ierr);

		// 3. New position
		ierr = ADVelAdvectCoord(vi->interp, vi->nmark, dt, 1); CHKERRQ(ierr);
	}

	// ---------------------------------
	// Runge-Kutta 2nd order in space
	// ---------------------------------
	else if(actx->advect == RUNGE_KUTTA_2)
	{
		// velocity interpolation A
		ierr = ADVelInterpMain(vi); CHKERRQ(ierr);

		// Runge-Kutta step to B
		ierr = ADVelRungeKuttaStep(vi, dt/2, 1.0, 0); CHKERRQ(ierr);

		// needed for mapping between vi and actx in parallel
		ierr = ADVelResetCoord(vi->interp, vi->nmark); CHKERRQ(ierr);
		ierr = ADVelExchange(vi); CHKERRQ(ierr);

		// final position
		ierr = ADVelAdvectCoord(vi->interp, vi->nmark, dt, 1); CHKERRQ(ierr);
	}

	//=======================================================================
	// END ADVECTION
	//=======================================================================

	// retrieve advected marker positions
	ierr = ADVelRetrieveCoord(actx, vi->interp, vi->nmark); CHKERRQ(ierr);

	// prepare indices for deletion in actx
	ierr = ADVelCollectIndices(actx, vi); CHKERRQ(ierr);

	// delete outside markers
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// clear memory
	ierr = ADVelDestroy(vi);      CHKERRQ(ierr);
	ierr = PetscFree(actx->idel); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelRungeKuttaStep(AdvVelCtx *vi, PetscScalar dt, PetscScalar a, PetscInt type)
{
	// routines to perform one runge-kutta step

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// 1. New position
	ierr = ADVelAdvectCoord(vi->interp, vi->nmark, dt, type); CHKERRQ(ierr);

	// 2. Delete marker outflow if it happens
	ierr = ADVelDeleteOutflow(vi); CHKERRQ(ierr);

	// 3. Exchange markers for interpolation
	ierr = ADVelExchange(vi); CHKERRQ(ierr);

	// 4. Velocity Interpolation
	ierr = ADVelInterpMain(vi); CHKERRQ(ierr);

	// 5. Update effective velocity
	ierr = ADVelCalcEffVel(vi->interp, vi->nmark, a); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelCreate(AdvCtx *actx, AdvVelCtx *vi)
{
	// create advection velocity context

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	vi->fs   = actx->fs;
	vi->jr   = actx->jr;
	vi->actx = actx;

	//=============
	// COMMUNICATOR
	//=============
	vi->icomm = actx->icomm;
	vi->nproc = actx->nproc;
	vi->iproc = actx->iproc;

	//========
	// STORAGE
	//========
	// allocate memory for interpolation routines
	vi->nmark = actx->nummark;
	vi->nbuff = actx->markcap;

	// allocate memory for interp markers
	ierr = PetscMalloc((size_t)vi->nbuff*sizeof(VelInterp), &vi->interp); CHKERRQ(ierr);
	ierr = PetscMemzero(vi->interp, (size_t)vi->nbuff*sizeof(VelInterp)); CHKERRQ(ierr);

	//========================
	// MARKER-CELL INTERACTION
	//========================
	ierr = makeIntArray(&vi->cellnum  , actx->cellnum, vi->nbuff       ); CHKERRQ(ierr);
	ierr = makeIntArray(&vi->markind  , NULL         , vi->nbuff       ); CHKERRQ(ierr);
	ierr = makeIntArray(&vi->markstart, NULL         , vi->fs->nCells+1); CHKERRQ(ierr);

	//=========
	// EXCHANGE
	//=========
	vi->sendbuf = NULL;
	vi->recvbuf = NULL;

	vi->nsend = 0;
	ierr = PetscMemzero(vi->nsendm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);
	ierr = PetscMemzero(vi->ptsend, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	vi->nrecv = 0;
	ierr = PetscMemzero(vi->nrecvm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);
	ierr = PetscMemzero(vi->ptrecv, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	vi->ndel = 0;
	vi->idel = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelDestroy(AdvVelCtx *vi)
{
	// destroy advection velocity context

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscFree(vi->interp);    CHKERRQ(ierr);
	ierr = PetscFree(vi->cellnum);   CHKERRQ(ierr);
	ierr = PetscFree(vi->markind);   CHKERRQ(ierr);
	ierr = PetscFree(vi->markstart); CHKERRQ(ierr);
	ierr = PetscFree(vi->sendbuf);   CHKERRQ(ierr);
	ierr = PetscFree(vi->recvbuf);   CHKERRQ(ierr);
	ierr = PetscFree(vi->idel);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelExchange(AdvVelCtx *vi)
{
	// exchange interpolated marker positions between the processors

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// count number of markers to be sent to each neighbor domain
	ierr = ADVelMapToDomains(vi); CHKERRQ(ierr);

	// communicate number of markers with neighbor processes
	ierr = ADVelExchangeNMark(vi); CHKERRQ(ierr);

	// create send and receive buffers for asynchronous MPI communication
	ierr = ADVelCreateMPIBuff(vi); CHKERRQ(ierr);

	// communicate markers with neighbor processes
	ierr = ADVelExchangeMark(vi); CHKERRQ(ierr);

	// store received markers, collect garbage
	ierr = ADVelCollectGarbage(vi); CHKERRQ(ierr);

	// free communication buffer
	ierr = ADVelDestroyMPIBuff(vi); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelInitCoord(AdvCtx *actx, VelInterp *interp, PetscInt n)
{
	PetscInt     jj;

	PetscFunctionBeginUser;

	// scan all markers
	for(jj = 0; jj < n; jj++)
	{
		// get marker coordinates - initial
		interp[jj].x0[0] = actx->markers[jj].X[0];
		interp[jj].x0[1] = actx->markers[jj].X[1];
		interp[jj].x0[2] = actx->markers[jj].X[2];

		// get marker coordinates - to later interpolate
		interp[jj].x[0] = actx->markers[jj].X[0];
		interp[jj].x[1] = actx->markers[jj].X[1];
		interp[jj].x[2] = actx->markers[jj].X[2];

		// initialize effective velocity
		interp[jj].v_eff[0] = 0.0;
		interp[jj].v_eff[1] = 0.0;
		interp[jj].v_eff[2] = 0.0;

		// save index
		interp[jj].ind = jj;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelCalcEffVel(VelInterp *interp, PetscInt n, PetscScalar a)
{
	PetscInt     jj;

	PetscFunctionBeginUser;

	// scan all markers
	for(jj = 0; jj < n; jj++)
	{
		interp[jj].v_eff[0] += a*interp[jj].v[0];
		interp[jj].v_eff[1] += a*interp[jj].v[1];
		interp[jj].v_eff[2] += a*interp[jj].v[2];
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelAdvectCoord(VelInterp *interp, PetscInt n, PetscScalar dt, PetscInt type)
{
	PetscInt     jj;

	PetscFunctionBeginUser;

	// scan all markers
	for(jj = 0; jj < n; jj++)
	{
		if (type==1)
		{
			// final advection of marker
			interp[jj].x[0] += interp[jj].v_eff[0]*dt;
			interp[jj].x[1] += interp[jj].v_eff[1]*dt;
			interp[jj].x[2] += interp[jj].v_eff[2]*dt;
		}
		else
		{
			// intermediate advection of marker
			interp[jj].x[0] = interp[jj].x0[0] + interp[jj].v[0]*dt;
			interp[jj].x[1] = interp[jj].x0[1] + interp[jj].v[1]*dt;
			interp[jj].x[2] = interp[jj].x0[2] + interp[jj].v[2]*dt;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelResetCoord(VelInterp *interp, PetscInt n)
{
	PetscInt     jj;

	PetscFunctionBeginUser;

	// scan all markers
	for(jj = 0; jj < n; jj++)
	{
		// reset coord
		interp[jj].x[0] = interp[jj].x0[0];
		interp[jj].x[1] = interp[jj].x0[1];
		interp[jj].x[2] = interp[jj].x0[2];
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelRetrieveCoord(AdvCtx *actx, VelInterp *interp, PetscInt n)
{
	PetscInt     jj, p;
	PetscFunctionBeginUser;

	// scan all markers
	for(jj = 0; jj < n; jj++)
	{
		// get marker coordinates
		p = interp[jj].ind;
		actx->markers[p].X[0] = interp[jj].x[0];
		actx->markers[p].X[1] = interp[jj].x[1];
		actx->markers[p].X[2] = interp[jj].x[2];

		// displacement
		actx->markers[p].U[0] += interp[jj].x[0] - interp[jj].x0[0];
		actx->markers[p].U[1] += interp[jj].x[1] - interp[jj].x0[1];
		actx->markers[p].U[2] += interp[jj].x[2] - interp[jj].x0[2];
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelCollectIndices(AdvCtx *actx, AdvVelCtx *vi)
{
	// collect indices to delete markers that went outside the domain during the runge-kutta steps

	PetscInt     jj, ind, ndel, *p;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// number to delete
	actx->ndel  = actx->nummark-vi->nmark;
	actx->nrecv = 0;

	// if no need for injection/deletion
	if (!actx->ndel) PetscFunctionReturn(0);

	// allocate storage
	if(actx->ndel) { ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel); CHKERRQ(ierr); }

	// allocate storage for mapping
	ierr = PetscMalloc((size_t)actx->nummark*sizeof(PetscInt), &p); CHKERRQ(ierr);
	ierr = PetscMemzero(p, (size_t)actx->nummark*sizeof(PetscInt)); CHKERRQ(ierr);

	// scan all advected markers
	for(jj = 0; jj < vi->nmark; jj++)
	{
		// get marker index
		ind = vi->interp[jj].ind;
		p[ind] = 1;
	}

	// scan mapping array and delete the markers
	for(jj = 0, ndel = 0; jj < actx->nummark; jj++)
	{
		if (p[jj]==0) actx->idel[ndel++] = jj;
	}

	// free memory
	ierr = PetscFree(p); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelDeleteOutflow(AdvVelCtx *vi)
{
	// check if advected positions are within the box bounds

	PetscInt     i, lrank, ndel;
	PetscMPIInt  grank;
	FDSTAG       *fs;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// scan all markers
	for(i = 0, ndel = 0; i < vi->nmark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, vi->interp[i].x, &lrank, &grank); CHKERRQ(ierr);

		// count markers outside
		if(grank == -1) ndel++;
	}

	// if no need for deletion return
	if (!ndel) PetscFunctionReturn(0);

	// save points to exclude
	vi->ndel    = ndel;

	// allocate storage
	ierr = PetscMalloc((size_t)ndel *sizeof(PetscInt), &vi->idel); CHKERRQ(ierr);

	// save markers indices to be deleted
	for(i = 0, ndel = 0; i < vi->nmark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, vi->interp[i].x, &lrank, &grank); CHKERRQ(ierr);

		// save markers outside
		if(grank == -1)
		{
			vi->idel[ndel++] = i;
		}
	}

	// delete outside markers
	ierr = ADVelCollectGarbage(vi); CHKERRQ(ierr);

	// clear
	ierr = PetscFree(vi->idel); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelMapToDomains(AdvVelCtx *vi)
{
	// count number of positions to be sent to each neighbor domain

	PetscInt     i, lrank, cnt;
	PetscMPIInt  grank;
	FDSTAG      *fs;

	PetscErrorCode  ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// clear send counters
	ierr = PetscMemzero(vi->nsendm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	// scan markers
	for(i = 0, cnt = 0; i < vi->nmark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, vi->interp[i].x, &lrank, &grank); CHKERRQ(ierr);

		if(grank != vi->iproc)
		{
			// count markers that should be sent to each neighbor
			vi->nsendm[lrank]++;
			cnt++;
		}
	}

	// store number of deleted markers
	vi->ndel = cnt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelExchangeNMark(AdvVelCtx *vi)
{
	// communicate number of markers with neighbor processes

	FDSTAG     *fs;
	PetscInt    k;
	PetscMPIInt scnt, rcnt;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// zero out message counters
	scnt = 0;
	rcnt = 0;

	// send number of markers to ALL neighbor processes (except self & non-existing)
	for(k = 0; k < _num_neighb_; k++)
	{
		if(fs->neighb[k] != vi->iproc && fs->neighb[k] != -1)
		{
			ierr = MPI_Isend(&vi->nsendm[k], 1, MPIU_INT,
					fs->neighb[k], 100, vi->icomm, &srequest[scnt++]); CHKERRQ(ierr);
		}
	}

	// receive number of markers from ALL neighbor processes (except self & non-existing)
	for(k = 0; k < _num_neighb_; k++)
	{
		if(fs->neighb[k] != vi->iproc && fs->neighb[k] != -1)
		{
			ierr = MPI_Irecv(&vi->nrecvm[k], 1, MPIU_INT,
					fs->neighb[k], 100, vi->icomm, &rrequest[rcnt++]); CHKERRQ(ierr);
		}
		else vi->nrecvm[k] = 0;
	}

	// wait until all communication processes have been terminated
	if(scnt) { ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }
	if(rcnt) { ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelCreateMPIBuff(AdvVelCtx *vi)
{
	// create send and receive buffers for asynchronous MPI communication

	FDSTAG     *fs;
	PetscInt    i, cnt, lrank;
	PetscMPIInt grank;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// compute buffer pointers
	vi->nsend = getPtrCnt(_num_neighb_, vi->nsendm, vi->ptsend);
	vi->nrecv = getPtrCnt(_num_neighb_, vi->nrecvm, vi->ptrecv);

	vi->sendbuf = NULL;
	vi->recvbuf = NULL;
	vi->idel    = NULL;

	// allocate exchange buffers & array of deleted (sent) marker indices
	if(vi->nsend) { ierr = PetscMalloc((size_t)vi->nsend*sizeof(VelInterp), &vi->sendbuf); CHKERRQ(ierr); }
	if(vi->nrecv) { ierr = PetscMalloc((size_t)vi->nrecv*sizeof(VelInterp), &vi->recvbuf); CHKERRQ(ierr); }
	if(vi->ndel)  { ierr = PetscMalloc((size_t)vi->ndel *sizeof(PetscInt ), &vi->idel);    CHKERRQ(ierr); }

	// copy markers to send buffer, store their indices
	for(i = 0, cnt = 0; i < vi->nmark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, vi->interp[i].x, &lrank, &grank); CHKERRQ(ierr);

		if(grank != vi->iproc)
		{
			// store marker in the send buffer
			vi->sendbuf[vi->ptsend[lrank]++] = vi->interp[i];

			// delete marker from the storage
			vi->idel[cnt++] = i;
		}
	}

	// rewind send buffer pointers
	rewindPtr(_num_neighb_, vi->ptsend);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelExchangeMark(AdvVelCtx *vi)
{
	// communicate markers with neighbor processes

	FDSTAG     *fs;
	PetscInt    k;
	PetscMPIInt scnt, rcnt, nbyte;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// zero out message counters
	scnt = 0;
	rcnt = 0;

	// send packages (if any) with markers to neighbor processes
	for(k = 0; k < _num_neighb_; k++)
	{
		if(vi->nsendm[k])
		{
			nbyte = (PetscMPIInt)(vi->nsendm[k]*(PetscInt)sizeof(VelInterp));

			ierr = MPI_Isend(&vi->sendbuf[vi->ptsend[k]], nbyte, MPI_BYTE,
					fs->neighb[k], 200, vi->icomm, &srequest[scnt++]); CHKERRQ(ierr);

		}
	}

	// receive packages (if any) with markers from neighbor processes
	for(k = 0; k < _num_neighb_; k++)
	{
		if(vi->nrecvm[k])
		{
			nbyte = (PetscMPIInt)(vi->nrecvm[k]*(PetscInt)sizeof(VelInterp));

			ierr = MPI_Irecv(&vi->recvbuf[vi->ptrecv[k]], nbyte, MPI_BYTE,
					fs->neighb[k], 200, vi->icomm, &rrequest[rcnt++]); CHKERRQ(ierr);
		}
	}

	// wait until all communication processes have been terminated
	if(scnt) { ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }
	if(rcnt) { ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelDestroyMPIBuff(AdvVelCtx *vi)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// destroy buffers
	ierr = PetscFree(vi->sendbuf); CHKERRQ(ierr);
	ierr = PetscFree(vi->recvbuf); CHKERRQ(ierr);
	ierr = PetscFree(vi->idel);  	 CHKERRQ(ierr);

	// reset values
	vi->nrecv = 0;
	vi->ndel  = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelCollectGarbage(AdvVelCtx *vi)
{
	// store received markers, collect garbage

	VelInterp   *interp, *recvbuf;
	PetscInt    *idel, nmark, nrecv, ndel;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access storage
	nmark   = vi->nmark;
	interp  = vi->interp;

	nrecv   = vi->nrecv;
	recvbuf = vi->recvbuf;

	ndel    = vi->ndel;
	idel    = vi->idel;

	// close holes in marker storage
	while(nrecv && ndel)
	{
		interp[idel[ndel-1]] = recvbuf[nrecv-1];
		nrecv--;
		ndel--;
	}

	if(nrecv)
	{
		// make sure space is enough
		ierr = ADVelReAllocStorage(vi, nmark + nrecv); CHKERRQ(ierr);

		// make sure we have a correct storage pointer
		interp = vi->interp;

		// put the rest in the end of marker storage
		while(nrecv)
		{
			interp[nmark++] = recvbuf[nrecv-1];
			nrecv--;
		}
	}

	if(ndel)
	{
		// collect garbage
		while(ndel)
		{
			if(idel[ndel-1] != nmark-1)
			{
				interp[idel[ndel-1]] = interp[nmark-1];
			}
			nmark--;
			ndel--;
		}
	}

	// store new number of markers
	vi->nmark = nmark;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelReAllocStorage(AdvVelCtx *vi, PetscInt nmark)
{
	// re-allocate memory to dynamic storage

	PetscInt     nbuff;
	VelInterp   *interp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// check whether current storage is insufficient
	if(nmark > vi->nbuff)
	{
		// delete host cell numbers
		ierr = PetscFree(vi->cellnum); CHKERRQ(ierr);

		// compute new capacity
		nbuff = (PetscInt)(_cap_overhead_*(PetscScalar)nmark);

		// allocate memory
		ierr = PetscMalloc((size_t)nbuff*sizeof(VelInterp), &interp); CHKERRQ(ierr);
		ierr = PetscMemzero(interp, (size_t)nbuff*sizeof(VelInterp)); CHKERRQ(ierr);

		// copy current data
		if(vi->nmark)
		{
			ierr = PetscMemcpy(interp, vi->interp, (size_t)vi->nmark*sizeof(VelInterp)); CHKERRQ(ierr);
		}

		// delete previous storage
		ierr = PetscFree(vi->interp); CHKERRQ(ierr);

		// save new capacity & storage
		vi->nbuff  = nbuff;
		vi->interp = interp;

		// allocate memory for host cell numbers
		ierr = PetscMalloc((size_t)nbuff*sizeof(PetscInt), &vi->cellnum); CHKERRQ(ierr);
		ierr = PetscMemzero(vi->cellnum, (size_t)nbuff*sizeof(PetscInt)); CHKERRQ(ierr);

		// allocate memory for id marker arranging per cell
		ierr = PetscMalloc((size_t)nbuff*sizeof(PetscInt), &vi->markind); CHKERRQ(ierr);
		ierr = PetscMemzero(vi->markind, (size_t)nbuff*sizeof(PetscInt)); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode ADVelMapMarkToCells(AdvVelCtx *vi)
{
	// maps markers to cells (local)

	FDSTAG      *fs;
	PetscScalar *X;
	PetscInt     i, ID, I, J, K, M, N;
	PetscInt    *numMarkCell, *m, p;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs = vi->fs;

	// get number of cells
	M = fs->dsx.ncels;
	N = fs->dsy.ncels;

	// loop over all local particles
	for(i = 0; i < vi->nmark; i++)
	{
		// get marker coordinates
		X = vi->interp[i].x;

		// get host cell IDs in all directions
		ierr = Discret1DFindPoint(&fs->dsx, X[0], I); CHKERRQ(ierr);
		ierr = Discret1DFindPoint(&fs->dsy, X[1], J); CHKERRQ(ierr);
		ierr = Discret1DFindPoint(&fs->dsz, X[2], K); CHKERRQ(ierr);

		// compute and store consecutive index
		GET_CELL_ID(ID, I, J, K, M, N);

		vi->cellnum[i] = ID;
	}

	// allocate marker counter array
	ierr = makeIntArray(&numMarkCell, NULL, fs->nCells); CHKERRQ(ierr);

	// count number of markers in the cells
	for(i = 0; i < vi->nmark; i++) numMarkCell[vi->cellnum[i]]++;

	// store starting indices of markers belonging to a cell
	vi->markstart[0] = 0;
	for(i = 1; i < fs->nCells+1; i++) vi->markstart[i] = vi->markstart[i-1]+numMarkCell[i-1];

	// allocate memory for id offset
	ierr = makeIntArray(&m, NULL, fs->nCells); CHKERRQ(ierr);

	// store marker indices belonging to a cell
	for(i = 0; i < vi->nmark; i++)
	{
		p = vi->markstart[vi->cellnum[i]];
		vi->markind[p + m[vi->cellnum[i]]] = i;
		m[vi->cellnum[i]]++;
	}

	// free memory
	ierr = PetscFree(m);           CHKERRQ(ierr);
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------
PetscErrorCode ADVelInterpMain(AdvVelCtx *vi)
{
	//=======================================================================
	// VELOCITY INTERPOLATION ROUTINE
	//=======================================================================
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if     (vi->actx->interp == STAG   )  { ierr = ADVelInterpSTAG   (vi); CHKERRQ(ierr); }
	else if(vi->actx->interp == MINMOD )  { ierr = ADVelInterpMINMOD (vi); CHKERRQ(ierr); }
	else if(vi->actx->interp == STAG_P )  { ierr = ADVelInterpSTAGP  (vi); CHKERRQ(ierr); }
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER," *** Unknown option for velocity interpolation scheme");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelInterpSTAG(AdvVelCtx *vi)
{
	// interpolate velocities from STAG points to markers

	FDSTAG      *fs;
	JacRes      *jr;
	PetscInt    sx, sy, sz, nx, ny, nmark;
	PetscInt    jj, ID, I, J, K, II, JJ, KK;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar *ccx, *ccy, *ccz;
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar xc, yc, zc, xp, yp, zp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// compute host cells for all markers
	ierr = ADVelMapMarkToCells(vi); CHKERRQ(ierr);

	// access context
	fs = vi->fs;
	jr = vi->jr;

	nmark = vi->nmark;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// node & cell coordinates
	ncx = fs->dsx.ncoor; ccx = fs->dsx.ccoor;
	ncy = fs->dsy.ncoor; ccy = fs->dsy.ccoor;
	ncz = fs->dsz.ncoor; ccz = fs->dsz.ccoor;

	// access velocity, pressure & temperature vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	// scan all markers
	for(jj = 0; jj < nmark; jj++)
	{
		// get marker coordinates
		xp = vi->interp[jj].x[0];
		yp = vi->interp[jj].x[1];
		zp = vi->interp[jj].x[2];

		// get consecutive index of the host cell
		ID = vi->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get coordinates of cell center
		xc = ccx[I];
		yc = ccy[J];
		zc = ccz[K];

		// map marker on the cells of X, Y, Z & center grids
		if(xp > xc) { II = I; } else { II = I-1; }
		if(yp > yc) { JJ = J; } else { JJ = J-1; }
		if(zp > zc) { KK = K; } else { KK = K-1; }

		// interpolate velocity, pressure & temperature
		vi->interp[jj].v[0] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, xp, yp, zp, ncx, ccy, ccz);
		vi->interp[jj].v[1] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, xp, yp, zp, ccx, ncy, ccz);
		vi->interp[jj].v[2] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, xp, yp, zp, ccx, ccy, ncz);
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelInterpMINMOD(AdvVelCtx *vi)
{
	// interpolation to nodes from fdstag points using a MINMOD limiter
	// then trilinear interpolation to markers
	// with velocity correction from Jenny et al (2001), Meyer and Jenny (2004) and Wang et al (2015)

	FDSTAG      *fs;
	JacRes      *jr;
	PetscInt    sx, sy, sz, nx, ny,nz;
	PetscInt    i, j, k, I, J, K;
	PetscInt    ID, pind, jj, n;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar *ccx, *ccy, *ccz;
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar vxn[8], vyn[8], vzn[8];
	PetscScalar C12, C23, C31, C10, C20, C30;
	PetscScalar xp, yp, zp;
	PetscScalar nxs, nys, nzs, nxe, nye, nze, dx, dy, dz;
	PetscScalar xpl, ypl, zpl;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// compute host cells for all markers
	ierr = ADVelMapMarkToCells(vi); CHKERRQ(ierr);

	// access context
	fs = vi->fs;
	jr = vi->jr;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart; nz = fs->dsz.ncels;

	// node & cell coordinates
	ncx = fs->dsx.ncoor; ccx = fs->dsx.ccoor;
	ncy = fs->dsy.ncoor; ccy = fs->dsy.ccoor;
	ncz = fs->dsz.ncoor; ccz = fs->dsz.ccoor;

	// access velocity, pressure & temperature vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	// scan all local cells
	START_STD_LOOP
	{
		// get local indices
		I = i-sx;
		J = j-sy;
		K = k-sz;

		// get start/end coordinates of nodes
		nxs = ncx[I]; nxe = ncx[I+1];
		nys = ncy[J]; nye = ncy[J+1];
		nzs = ncz[K]; nze = ncz[K+1];

		// get cell dimensions
		dx = nxe-nxs;
		dy = nye-nys;
		dz = nze-nzs;

		// Linear and minmod limiter
		PetscScalar U1, U2, U3;
		// ---------------------------------------------------------------------
		// Vx 0,2
		U1     = InterpLinMinmod(lvx[k-1][j-1][i], lvx[k][j-1][i], lvx[k+1][j-1][i], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U2     = InterpLinMinmod(lvx[k-1][j  ][i], lvx[k][j  ][i], lvx[k+1][j  ][i], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U3     = InterpLinMinmod(lvx[k-1][j+1][i], lvx[k][j+1][i], lvx[k+1][j+1][i], dz, ccz[K-1], ccz[K], ccz[K+1], 0);

		vxn[0] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		vxn[2] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 1);

		// Vx 4,6
		U1     = InterpLinMinmod(lvx[k-1][j-1][i], lvx[k][j-1][i], lvx[k+1][j-1][i], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U2     = InterpLinMinmod(lvx[k-1][j  ][i], lvx[k][j  ][i], lvx[k+1][j  ][i], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U3     = InterpLinMinmod(lvx[k-1][j+1][i], lvx[k][j+1][i], lvx[k+1][j+1][i], dz, ccz[K-1], ccz[K], ccz[K+1], 1);

		vxn[4] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		vxn[6] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 1);

		// Vx 1,3
		U1     = InterpLinMinmod(lvx[k-1][j-1][i+1], lvx[k][j-1][i+1], lvx[k+1][j-1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U2     = InterpLinMinmod(lvx[k-1][j  ][i+1], lvx[k][j  ][i+1], lvx[k+1][j  ][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U3     = InterpLinMinmod(lvx[k-1][j+1][i+1], lvx[k][j+1][i+1], lvx[k+1][j+1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);

		vxn[1] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		vxn[3] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 1);

		// Vx 5,7
		U1     = InterpLinMinmod(lvx[k-1][j-1][i+1], lvx[k][j-1][i+1], lvx[k+1][j-1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U2     = InterpLinMinmod(lvx[k-1][j  ][i+1], lvx[k][j  ][i+1], lvx[k+1][j  ][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U3     = InterpLinMinmod(lvx[k-1][j+1][i+1], lvx[k][j+1][i+1], lvx[k+1][j+1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);

		vxn[5] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		vxn[7] = InterpLinMinmod(U1, U2, U3, dy, ccy[J-1], ccy[J], ccy[J+1], 1);
		// ---------------------------------------------------------------------
		// Vy 0,1
		U1     = InterpLinMinmod(lvy[k-1][j][i-1], lvy[k][j][i-1], lvy[k+1][j][i-1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U2     = InterpLinMinmod(lvy[k-1][j][i  ], lvy[k][j][i  ], lvy[k+1][j][i  ], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U3     = InterpLinMinmod(lvy[k-1][j][i+1], lvy[k][j][i+1], lvy[k+1][j][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);

		vyn[0] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vyn[1] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vy 4,5
		U1     = InterpLinMinmod(lvy[k-1][j][i-1], lvy[k][j][i-1], lvy[k+1][j][i-1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U2     = InterpLinMinmod(lvy[k-1][j][i  ], lvy[k][j][i  ], lvy[k+1][j][i  ], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U3     = InterpLinMinmod(lvy[k-1][j][i+1], lvy[k][j][i+1], lvy[k+1][j][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);

		vyn[4] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vyn[5] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vy 2,3
		U1     = InterpLinMinmod(lvy[k-1][j+1][i-1], lvy[k][j+1][i-1], lvy[k+1][j+1][i-1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U2     = InterpLinMinmod(lvy[k-1][j+1][i  ], lvy[k][j+1][i  ], lvy[k+1][j+1][i  ], dz, ccz[K-1], ccz[K], ccz[K+1], 0);
		U3     = InterpLinMinmod(lvy[k-1][j+1][i+1], lvy[k][j+1][i+1], lvy[k+1][j+1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 0);

		vyn[2] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vyn[3] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vy 6,7
		U1     = InterpLinMinmod(lvy[k-1][j+1][i-1], lvy[k][j+1][i-1], lvy[k+1][j+1][i-1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U2     = InterpLinMinmod(lvy[k-1][j+1][i  ], lvy[k][j+1][i  ], lvy[k+1][j+1][i  ], dz, ccz[K-1], ccz[K], ccz[K+1], 1);
		U3     = InterpLinMinmod(lvy[k-1][j+1][i+1], lvy[k][j+1][i+1], lvy[k+1][j+1][i+1], dz, ccz[K-1], ccz[K], ccz[K+1], 1);

		vyn[6] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vyn[7] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);
		// ---------------------------------------------------------------------
		// Vz 0,1
		U1     = InterpLinMinmod(lvz[k][j-1][i-1], lvz[k][j][i-1], lvz[k][j+1][i-1], dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		U2     = InterpLinMinmod(lvz[k][j-1][i  ], lvz[k][j][i  ], lvz[k][j+1][i  ], dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		U3     = InterpLinMinmod(lvz[k][j-1][i+1], lvz[k][j][i+1], lvz[k][j+1][i+1], dy, ccy[J-1], ccy[J], ccy[J+1], 0);

		vzn[0] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vzn[1] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vz 2,3
		U1     = InterpLinMinmod(lvz[k][j-1][i-1], lvz[k][j][i-1], lvz[k][j+1][i-1], dy, ccy[J-1], ccy[J], ccy[J+1], 1);
		U2     = InterpLinMinmod(lvz[k][j-1][i  ], lvz[k][j][i  ], lvz[k][j+1][i  ], dy, ccy[J-1], ccy[J], ccy[J+1], 1);
		U3     = InterpLinMinmod(lvz[k][j-1][i+1], lvz[k][j][i+1], lvz[k][j+1][i+1], dy, ccy[J-1], ccy[J], ccy[J+1], 1);

		vzn[2] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vzn[3] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vz 4,5
		U1     = InterpLinMinmod(lvz[k+1][j-1][i-1], lvz[k+1][j][i-1], lvz[k+1][j+1][i-1], dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		U2     = InterpLinMinmod(lvz[k+1][j-1][i  ], lvz[k+1][j][i  ], lvz[k+1][j+1][i  ], dy, ccy[J-1], ccy[J], ccy[J+1], 0);
		U3     = InterpLinMinmod(lvz[k+1][j-1][i+1], lvz[k+1][j][i+1], lvz[k+1][j+1][i+1], dy, ccy[J-1], ccy[J], ccy[J+1], 0);

		vzn[4] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vzn[5] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);

		// Vz 6,7
		U1     = InterpLinMinmod(lvz[k+1][j-1][i-1], lvz[k+1][j][i-1], lvz[k+1][j+1][i-1], dy, ccy[J-1], ccy[J], ccy[J+1], 1);
		U2     = InterpLinMinmod(lvz[k+1][j-1][i  ], lvz[k+1][j][i  ], lvz[k+1][j+1][i  ], dy, ccy[J-1], ccy[J], ccy[J+1], 1);
		U3     = InterpLinMinmod(lvz[k+1][j-1][i+1], lvz[k+1][j][i+1], lvz[k+1][j+1][i+1], dy, ccy[J-1], ccy[J], ccy[J+1], 1);

		vzn[6] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 0);
		vzn[7] = InterpLinMinmod(U1, U2, U3, dx, ccx[I-1], ccx[I], ccx[I+1], 1);
		// ---------------------------------------------------------------------

		// calculate corrections
		C12 = dx/2/dy*(-vyn[0]+vyn[1]+vyn[2]-vyn[3]+vyn[4]-vyn[5]-vyn[6]+vyn[7]);
		C23 = dz/2/dx*(-vxn[0]+vxn[1]+vxn[2]-vxn[3]+vxn[4]-vxn[5]-vxn[6]+vxn[7]);
		C31 = dy/2/dz*(-vzn[0]+vzn[1]+vzn[2]-vzn[3]+vzn[4]-vzn[5]-vzn[6]+vzn[7]);

		C10 = dx/2/dz*(vzn[0]-vzn[4]+vzn[5]-vzn[1])       + dx/2/dy*(vyn[0]-vyn[2]+vyn[3]-vyn[1] + C31);
		C20 = dz/2/dx*(vxn[0]-vxn[1]+vxn[5]-vxn[4] + C12) + dz/2/dy*(vyn[0]-vyn[2]+vyn[6]-vyn[4]      );
		C30 = dy/2/dx*(vxn[0]-vxn[1]+vxn[3]-vxn[2])       + dy/2/dz*(vzn[0]-vzn[4]+vzn[6]-vzn[2] + C23);

		// get cell id
		GET_CELL_ID(ID, I, J, K, nx, ny);

		// get markers in cell
		n = vi->markstart[ID+1] - vi->markstart[ID];

		// scan cell markers
		for(jj = 0; jj < n; jj++)
		{
			// get marker index
			pind = vi->markind[vi->markstart[ID] + jj];

			// get marker coordinates
			xp = vi->interp[pind].x[0];
			yp = vi->interp[pind].x[1];
			zp = vi->interp[pind].x[2];

			// transform into local coordinates
			xpl = (xp-nxs)/(nxe-nxs);
			ypl = (yp-nys)/(nye-nys);
			zpl = (zp-nzs)/(nze-nzs);

			// interpolate velocities
			vi->interp[pind].v[0] = GenInterpLin3D(vxn,xpl,ypl,zpl);
			vi->interp[pind].v[1] = GenInterpLin3D(vyn,xpl,ypl,zpl);
			vi->interp[pind].v[2] = GenInterpLin3D(vzn,xpl,ypl,zpl);

			// add correction
			vi->interp[pind].v[0] += xpl*(1-xpl)*(C10 + zpl*C12);
			vi->interp[pind].v[1] += ypl*(1-ypl)*(C30 + xpl*C31);
			vi->interp[pind].v[2] += zpl*(1-zpl)*(C20 + ypl*C23);
		}
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ADVelInterpSTAGP(AdvVelCtx *vi)
{
	// interpolate velocities from STAG points to markers + pressure points

	FDSTAG      *fs;
	JacRes      *jr;
	PetscInt    sx, sy, sz, nx, ny, nz, nmark;
	PetscInt    jj, ID, I, J, K, II, JJ, KK;
	PetscInt    IN, JN, KN;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar *ccx, *ccy, *ccz;
	PetscScalar v[3], vp[3];
	PetscScalar ***lvx, ***lvy, ***lvz;
	PetscScalar xc, yc, zc, xp, yp, zp;
	PetscScalar vii[12], vxp[8], vyp[8], vzp[8];
	PetscScalar nxs, nxe, nys, nye, nzs, nze, xpl, ypl, zpl;
	PetscScalar A, B;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// compute host cells for all markers
	ierr = ADVelMapMarkToCells(vi); CHKERRQ(ierr);

	// set coefficients
	if(vi->actx->A) A = vi->actx->A;
	else            A = 2.0/3.0;

	B = 1.0 - A;

	// access context
	fs = vi->fs;
	jr = vi->jr;

	nmark = vi->nmark;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart; nz = fs->dsz.ncels;

	// node & cell coordinates
	ncx = fs->dsx.ncoor; ccx = fs->dsx.ccoor;
	ncy = fs->dsy.ncoor; ccy = fs->dsy.ccoor;
	ncz = fs->dsz.ncoor; ccz = fs->dsz.ccoor;

	// access velocity, pressure & temperature vectors
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	// scan all markers
	for(jj = 0; jj < nmark; jj++)
	{
		// get marker coordinates
		xp = vi->interp[jj].x[0];
		yp = vi->interp[jj].x[1];
		zp = vi->interp[jj].x[2];

		// get consecutive index of the host cell
		ID = vi->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get coordinates of cell center
		xc = ccx[I];
		yc = ccy[J];
		zc = ccz[K];

		// map marker on the cells of X, Y, Z & center grids
		if(xp > xc) { II = I; } else { II = I-1; }
		if(yp > yc) { JJ = J; } else { JJ = J-1; }
		if(zp > zc) { KK = K; } else { KK = K-1; }

		// interpolate velocity, pressure & temperature
		v[0] = InterpLin3D(lvx, I,  JJ, KK, sx, sy, sz, xp, yp, zp, ncx, ccy, ccz);
		v[1] = InterpLin3D(lvy, II, J,  KK, sx, sy, sz, xp, yp, zp, ccx, ncy, ccz);
		v[2] = InterpLin3D(lvz, II, JJ, K,  sx, sy, sz, xp, yp, zp, ccx, ccy, ncz);

		// check marker position relative to pressure cube
		// WARNING! ONLY FOR FREE SLIP BC

		// ----------------
		// VX
		// ----------------
		if(xp > xc) { IN = I+1; } else { IN = I; }

		if (IN == 0)
		{
			vii[0] = 2*lvx[sz+KK  ][sy+JJ  ][sx+IN]-lvx[sz+KK  ][sy+JJ  ][sx+IN+1]; vii[1]  = lvx[sz+KK  ][sy+JJ  ][sx+IN  ]; vii[2]  = lvx[sz+KK  ][sy+JJ  ][sx+IN+1];
			vii[3] = 2*lvx[sz+KK  ][sy+JJ+1][sx+IN]-lvx[sz+KK  ][sy+JJ+1][sx+IN+1]; vii[4]  = lvx[sz+KK  ][sy+JJ+1][sx+IN  ]; vii[5]  = lvx[sz+KK  ][sy+JJ+1][sx+IN+1];
			vii[6] = 2*lvx[sz+KK+1][sy+JJ  ][sx+IN]-lvx[sz+KK+1][sy+JJ  ][sx+IN+1]; vii[7]  = lvx[sz+KK+1][sy+JJ  ][sx+IN  ]; vii[8]  = lvx[sz+KK+1][sy+JJ  ][sx+IN+1];
			vii[9] = 2*lvx[sz+KK+1][sy+JJ+1][sx+IN]-lvx[sz+KK+1][sy+JJ+1][sx+IN+1]; vii[10] = lvx[sz+KK+1][sy+JJ+1][sx+IN  ]; vii[11] = lvx[sz+KK+1][sy+JJ+1][sx+IN+1];

			nxs = 2*ccx[IN]-ccx[IN+1]; nxe = ccx[IN];
		}
		else if (IN == nx)
		{
			vii[0] = lvx[sz+KK  ][sy+JJ  ][sx+IN-1]; vii[1]  = lvx[sz+KK  ][sy+JJ  ][sx+IN  ]; vii[2]  = 2*lvx[sz+KK  ][sy+JJ  ][sx+IN]-lvx[sz+KK  ][sy+JJ  ][sx+IN-1];
			vii[3] = lvx[sz+KK  ][sy+JJ+1][sx+IN-1]; vii[4]  = lvx[sz+KK  ][sy+JJ+1][sx+IN  ]; vii[5]  = 2*lvx[sz+KK  ][sy+JJ+1][sx+IN]-lvx[sz+KK  ][sy+JJ+1][sx+IN-1];
			vii[6] = lvx[sz+KK+1][sy+JJ  ][sx+IN-1]; vii[7]  = lvx[sz+KK+1][sy+JJ  ][sx+IN  ]; vii[8]  = 2*lvx[sz+KK+1][sy+JJ  ][sx+IN]-lvx[sz+KK+1][sy+JJ  ][sx+IN-1];
			vii[9] = lvx[sz+KK+1][sy+JJ+1][sx+IN-1]; vii[10] = lvx[sz+KK+1][sy+JJ+1][sx+IN  ]; vii[11] = 2*lvx[sz+KK+1][sy+JJ+1][sx+IN]-lvx[sz+KK+1][sy+JJ+1][sx+IN-1];

			nxs = ccx[IN-1]; nxe = 2*ccx[IN]-ccx[IN-1];
		}
		else
		{
			vii[0] = lvx[sz+KK  ][sy+JJ  ][sx+IN-1]; vii[1]  = lvx[sz+KK  ][sy+JJ  ][sx+IN  ]; vii[2]  = lvx[sz+KK  ][sy+JJ  ][sx+IN+1];
			vii[3] = lvx[sz+KK  ][sy+JJ+1][sx+IN-1]; vii[4]  = lvx[sz+KK  ][sy+JJ+1][sx+IN  ]; vii[5]  = lvx[sz+KK  ][sy+JJ+1][sx+IN+1];
			vii[6] = lvx[sz+KK+1][sy+JJ  ][sx+IN-1]; vii[7]  = lvx[sz+KK+1][sy+JJ  ][sx+IN  ]; vii[8]  = lvx[sz+KK+1][sy+JJ  ][sx+IN+1];
			vii[9] = lvx[sz+KK+1][sy+JJ+1][sx+IN-1]; vii[10] = lvx[sz+KK+1][sy+JJ+1][sx+IN  ]; vii[11] = lvx[sz+KK+1][sy+JJ+1][sx+IN+1];

			nxs = ccx[IN-1]; nxe = ccx[IN  ];
		}

		// get velocities and boundaries
		vxp[0] = (vii[0] + vii[1] )/2; vxp[1] = (vii[1]  + vii[2] )/2;
		vxp[2] = (vii[3] + vii[4] )/2; vxp[3] = (vii[4]  + vii[5] )/2;
		vxp[4] = (vii[6] + vii[7] )/2; vxp[5] = (vii[7]  + vii[8] )/2;
		vxp[6] = (vii[9] + vii[10])/2; vxp[7] = (vii[10] + vii[11])/2;

		// transform into local coordinates
		nys = ccy[JJ  ]; nye = ccy[JJ+1];
		nzs = ccz[KK  ]; nze = ccz[KK+1];

		xpl = (xp-nxs)/(nxe-nxs);
		ypl = (yp-nys)/(nye-nys);
		zpl = (zp-nzs)/(nze-nzs);

		// calculate velocity in pressure points
		vp[0] = GenInterpLin3D(vxp,xpl,ypl,zpl);

		// ----------------
		// VY
		// ----------------
		if(yp > yc) { JN = J+1; } else { JN = J; }

		if (JN == 0)
		{
			vii[0] = 2*lvy[sz+KK  ][sy+JN][sx+II  ]-lvy[sz+KK  ][sy+JN+1][sx+II  ]; vii[1]  = lvy[sz+KK  ][sy+JN][sx+II  ]; vii[2]  = lvy[sz+KK  ][sy+JN+1][sx+II  ];
			vii[3] = 2*lvy[sz+KK  ][sy+JN][sx+II+1]-lvy[sz+KK  ][sy+JN+1][sx+II+1]; vii[4]  = lvy[sz+KK  ][sy+JN][sx+II+1]; vii[5]  = lvy[sz+KK  ][sy+JN+1][sx+II+1];
			vii[6] = 2*lvy[sz+KK+1][sy+JN][sx+II  ]-lvy[sz+KK+1][sy+JN+1][sx+II  ]; vii[7]  = lvy[sz+KK+1][sy+JN][sx+II  ]; vii[8]  = lvy[sz+KK+1][sy+JN+1][sx+II  ];
			vii[9] = 2*lvy[sz+KK+1][sy+JN][sx+II+1]-lvy[sz+KK+1][sy+JN+1][sx+II+1]; vii[10] = lvy[sz+KK+1][sy+JN][sx+II+1]; vii[11] = lvy[sz+KK+1][sy+JN+1][sx+II+1];

			nys = 2*ccy[JN]-ccy[JN+1]; nye = ccy[JN];
		}
		else if (JN == ny)
		{
			vii[0] = lvy[sz+KK  ][sy+JN-1][sx+II  ]; vii[1]  = lvy[sz+KK  ][sy+JN][sx+II  ]; vii[2]  = 2*lvy[sz+KK  ][sy+JN][sx+II  ]-lvy[sz+KK  ][sy+JN-1][sx+II  ];
			vii[3] = lvy[sz+KK  ][sy+JN-1][sx+II+1]; vii[4]  = lvy[sz+KK  ][sy+JN][sx+II+1]; vii[5]  = 2*lvy[sz+KK  ][sy+JN][sx+II+1]-lvy[sz+KK  ][sy+JN-1][sx+II+1];
			vii[6] = lvy[sz+KK+1][sy+JN-1][sx+II  ]; vii[7]  = lvy[sz+KK+1][sy+JN][sx+II  ]; vii[8]  = 2*lvy[sz+KK+1][sy+JN][sx+II  ]-lvy[sz+KK+1][sy+JN-1][sx+II  ];
			vii[9] = lvy[sz+KK+1][sy+JN-1][sx+II+1]; vii[10] = lvy[sz+KK+1][sy+JN][sx+II+1]; vii[11] = 2*lvy[sz+KK+1][sy+JN][sx+II+1]-lvy[sz+KK+1][sy+JN-1][sx+II+1];

			nys = ccy[JN-1]; nye = 2*ccy[JN]-ccy[JN-1];
		}
		else
		{
			vii[0] = lvy[sz+KK  ][sy+JN-1][sx+II  ]; vii[1]  = lvy[sz+KK  ][sy+JN][sx+II  ]; vii[2]  = lvy[sz+KK  ][sy+JN+1][sx+II  ];
			vii[3] = lvy[sz+KK  ][sy+JN-1][sx+II+1]; vii[4]  = lvy[sz+KK  ][sy+JN][sx+II+1]; vii[5]  = lvy[sz+KK  ][sy+JN+1][sx+II+1];
			vii[6] = lvy[sz+KK+1][sy+JN-1][sx+II  ]; vii[7]  = lvy[sz+KK+1][sy+JN][sx+II  ]; vii[8]  = lvy[sz+KK+1][sy+JN+1][sx+II  ];
			vii[9] = lvy[sz+KK+1][sy+JN-1][sx+II+1]; vii[10] = lvy[sz+KK+1][sy+JN][sx+II+1]; vii[11] = lvy[sz+KK+1][sy+JN+1][sx+II+1];

			nys = ccy[JN-1]; nye = ccy[JN  ];
		}

		// get velocities and boundaries
		vyp[0] = (vii[0] + vii[1] )/2; vyp[1] = (vii[3]  + vii[4] )/2;
		vyp[2] = (vii[1] + vii[2] )/2; vyp[3] = (vii[4]  + vii[5] )/2;
		vyp[4] = (vii[6] + vii[7] )/2; vyp[5] = (vii[9]  + vii[10])/2;
		vyp[6] = (vii[7] + vii[8] )/2; vyp[7] = (vii[10] + vii[11])/2;

		// transform into local coordinates
		nxs = ccx[II  ]; nxe = ccx[II+1];
		nzs = ccz[KK  ]; nze = ccz[KK+1];

		xpl = (xp-nxs)/(nxe-nxs);
		ypl = (yp-nys)/(nye-nys);
		zpl = (zp-nzs)/(nze-nzs);

		// calculate velocity in pressure points
		vp[1] = GenInterpLin3D(vyp,xpl,ypl,zpl);

		// ----------------
		// VZ
		// ----------------
		if(zp > zc) { KN = K+1; } else { KN = K; }

		if (KN == 0)
		{
			vii[0] = 2*lvz[sz+KN][sy+JJ  ][sx+II  ]-lvz[sz+KN+1][sy+JJ  ][sx+II  ]; vii[1]  = lvz[sz+KN][sy+JJ  ][sx+II  ]; vii[2]  = lvz[sz+KN+1][sy+JJ  ][sx+II  ];
			vii[3] = 2*lvz[sz+KN][sy+JJ  ][sx+II+1]-lvz[sz+KN+1][sy+JJ  ][sx+II+1]; vii[4]  = lvz[sz+KN][sy+JJ  ][sx+II+1]; vii[5]  = lvz[sz+KN+1][sy+JJ  ][sx+II+1];
			vii[6] = 2*lvz[sz+KN][sy+JJ+1][sx+II  ]-lvz[sz+KN+1][sy+JJ+1][sx+II  ]; vii[7]  = lvz[sz+KN][sy+JJ+1][sx+II  ]; vii[8]  = lvz[sz+KN+1][sy+JJ+1][sx+II  ];
			vii[9] = 2*lvz[sz+KN][sy+JJ+1][sx+II+1]-lvz[sz+KN+1][sy+JJ+1][sx+II+1]; vii[10] = lvz[sz+KN][sy+JJ+1][sx+II+1]; vii[11] = lvz[sz+KN+1][sy+JJ+1][sx+II+1];

			nzs = 2*ccz[KN]-ccz[KN+1]; nze = ccz[KN];
		}
		else if (KN == nz)
		{
			vii[0] = lvz[sz+KN-1][sy+JJ  ][sx+II  ]; vii[1]  = lvz[sz+KN][sy+JJ  ][sx+II  ]; vii[2]  = 2*lvz[sz+KN][sy+JJ  ][sx+II  ]-lvz[sz+KN-1][sy+JJ  ][sx+II  ];
			vii[3] = lvz[sz+KN-1][sy+JJ  ][sx+II+1]; vii[4]  = lvz[sz+KN][sy+JJ  ][sx+II+1]; vii[5]  = 2*lvz[sz+KN][sy+JJ  ][sx+II+1]-lvz[sz+KN-1][sy+JJ  ][sx+II+1];
			vii[6] = lvz[sz+KN-1][sy+JJ+1][sx+II  ]; vii[7]  = lvz[sz+KN][sy+JJ+1][sx+II  ]; vii[8]  = 2*lvz[sz+KN][sy+JJ+1][sx+II  ]-lvz[sz+KN-1][sy+JJ+1][sx+II  ];
			vii[9] = lvz[sz+KN-1][sy+JJ+1][sx+II+1]; vii[10] = lvz[sz+KN][sy+JJ+1][sx+II+1]; vii[11] = 2*lvz[sz+KN][sy+JJ+1][sx+II+1]-lvz[sz+KN-1][sy+JJ+1][sx+II+1];

			nzs = ccz[KN-1]; nze = 2*ccz[KN]-ccz[KN-1];
		}
		else
		{
			vii[0] = lvz[sz+KN-1][sy+JJ  ][sx+II  ]; vii[1]  = lvz[sz+KN][sy+JJ  ][sx+II  ]; vii[2]  = lvz[sz+KN+1][sy+JJ  ][sx+II  ];
			vii[3] = lvz[sz+KN-1][sy+JJ  ][sx+II+1]; vii[4]  = lvz[sz+KN][sy+JJ  ][sx+II+1]; vii[5]  = lvz[sz+KN+1][sy+JJ  ][sx+II+1];
			vii[6] = lvz[sz+KN-1][sy+JJ+1][sx+II  ]; vii[7]  = lvz[sz+KN][sy+JJ+1][sx+II  ]; vii[8]  = lvz[sz+KN+1][sy+JJ+1][sx+II  ];
			vii[9] = lvz[sz+KN-1][sy+JJ+1][sx+II+1]; vii[10] = lvz[sz+KN][sy+JJ+1][sx+II+1]; vii[11] = lvz[sz+KN+1][sy+JJ+1][sx+II+1];

			nzs = ccz[KN-1]; nze = ccz[KN  ];
		}

		// get velocities and boundaries
		vzp[0] = (vii[0] + vii[1] )/2; vzp[1] = (vii[3]  + vii[4] )/2;
		vzp[2] = (vii[6] + vii[7] )/2; vzp[3] = (vii[9]  + vii[10])/2;
		vzp[4] = (vii[1] + vii[2] )/2; vzp[5] = (vii[4]  + vii[5] )/2;
		vzp[6] = (vii[7] + vii[8] )/2; vzp[7] = (vii[10] + vii[11])/2;

		// transform into local coordinates
		nxs = ccx[II  ]; nxe = ccx[II+1];
		nys = ccy[JJ  ]; nye = ccy[JJ+1];

		xpl = (xp-nxs)/(nxe-nxs);
		ypl = (yp-nys)/(nye-nys);
		zpl = (zp-nzs)/(nze-nzs);

		// calculate velocity in pressure points
		vp[2] = GenInterpLin3D(vzp,xpl,ypl,zpl);

		// interpolate velocity
		vi->interp[jj].v[0] = A*v[0] + B*vp[0];
		vi->interp[jj].v[1] = A*v[1] + B*vp[1];
		vi->interp[jj].v[2] = A*v[2] + B*vp[2];
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
