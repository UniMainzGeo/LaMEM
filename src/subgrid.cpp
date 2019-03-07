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
 **    filename:   subgrid.c
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
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "subgrid.h"
#include "advect.h"
#include "phase.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "bc.h"
#include "JacRes.h"
#include "surf.h"
#include "marker.h"
#include "AVD.h"
#include "cvi.h"
#include "tools.h"

/*
#START_DOC#
\lamemfunction{\verb- ADVCreate -}
Create advection context

\lamemfunction{\verb- ADVDestroy -}
Destroy advection context

\lamemfunction{\verb- ADVAdvect -}
Main advection routine

#END_DOC#
*/


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSubGrid"
PetscErrorCode ADVMarkSubGrid(AdvCtx *actx)
{

/*
	// check marker distribution and delete or inject markers if necessary
	FDSTAG         *fs;
	PetscScalar    xs[3], xe[3];
	PetscInt       ind, i, j, k, M, N;
	PetscInt       n, ninj, ndel;
	PetscLogDouble t0,t1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	ierr = PetscTime(&t0); CHKERRQ(ierr);

	// get number of cells
	M = fs->dsx.ncels;
	N = fs->dsy.ncels;

	// calculate storage
	ninj = 0;
	ndel = 0;
	for(i = 0; i < fs->nCells; i++)
	{
		// no of markers in cell
		n = actx->markstart[i+1] - actx->markstart[i];

		if (n < actx->nmin)
		{
			if ((actx->nmin - n) > n) ninj += n;
			else                      ninj += actx->nmin - n;
		}
		if (n > actx->nmax) ndel += n - actx->nmax;
	}

	// if no need for injection/deletion
	if ((!ninj) && (!ndel)) PetscFunctionReturn(0);

	actx->nrecv = ninj;
	actx->ndel  = ndel;

	// allocate memory
	if(ninj) { ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker),   &actx->recvbuf); CHKERRQ(ierr); }
	if(ndel) { ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel);    CHKERRQ(ierr); }

	actx->cinj = 0;
	actx->cdel = 0;
	ind        = 0;

	// inject/delete
	for(ind = 0; ind < fs->nCells; ind++)
	{
		// no of markers in cell
		n = actx->markstart[ind+1] - actx->markstart[ind];

		if ((n < actx->nmin) || (n > actx->nmax))
		{
			// expand i, j, k cell indices
			GET_CELL_IJK(ind, i, j, k, M, N);

			// get cell coordinates
			xs[0] = fs->dsx.ncoor[i]; xe[0] = fs->dsx.ncoor[i+1];
			xs[1] = fs->dsy.ncoor[j]; xe[1] = fs->dsy.ncoor[j+1];
			xs[2] = fs->dsz.ncoor[k]; xe[2] = fs->dsz.ncoor[k+1];

			// inject/delete markers
			ierr = AVDExecuteMarkerInjection(actx, n, xs, xe, ind); CHKERRQ(ierr);
		}
	}

	// store new markers
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// compute host cells for all the markers
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// update arrays for marker-cell interaction
	ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Marker control [%lld]: (AVD Cell) injected %lld markers and deleted %lld markers in %1.4e s\n",(LLD)actx->iproc, (LLD)ninj, (LLD)ndel, t1-t0);

	// clear
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->idel);    CHKERRQ(ierr);
*/
	PetscFunctionReturn(0);
}


/*
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVAdvect"
PetscErrorCode ADVAdvect(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION ROUTINE
	//=======================================================================

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->advect == ADV_NONE) PetscFunctionReturn(0);

	// project history INCREMENTS from grid to markers
	ierr = ADVProjHistGridToMark(actx); CHKERRQ(ierr);

	if(actx->advect != BASIC_EULER)
	{
		// advect markers (extended by Adina)
		ierr = ADVelAdvectMain(actx); CHKERRQ(ierr);
	}
	else
	{
		// advect markers (Forward Euler)
		ierr = ADVAdvectMark(actx); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVRemap"
PetscErrorCode ADVRemap(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION REMAPPING
	//=======================================================================

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->advect == ADV_NONE)
	{
		ierr = ADVUpdateHistADVNone(actx); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}
	if(actx->mctrl == CTRL_NONE)
	{
		// compute host cells for all the markers received
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);
	}
	else if(actx->mctrl == CTRL_BASIC)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Performing marker control (standard algorithm)\n");

		// compute host cells for all the markers received
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

		// check markers and inject/delete if necessary
		ierr = ADVMarkControl(actx); CHKERRQ(ierr);

		// check corners and inject 1 particle if empty
		ierr = ADVCheckCorners(actx); CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}
	else if(actx->mctrl == CTRL_AVD)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Performing marker control (updated algorithm)\n");

		// check markers and inject/delete if necessary in all control volumes
		ierr = AVDMarkerControl(actx); CHKERRQ(ierr);

		// compute host cells for all the markers received
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}
	else if(actx->mctrl == CTRL_SUB)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Performing marker control (subgrid algorithm)\n");

		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "... ooops subgrid marker control lacks implementation ...");

		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	// change marker phase when crossing flat surface or free surface with fast sedimentation/erosion
	ierr = ADVMarkCrossFreeSurf(actx); CHKERRQ(ierr);

	// project advected history from markers back to grid
	ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVExchange"
PetscErrorCode ADVExchange(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION EXCHANGE ROUTINE
	//
	// Exchange markers between the processors resulting from the position change
	//=======================================================================

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->advect == ADV_NONE) PetscFunctionReturn(0);

	// delete marker outflow if it happens
	ierr = ADVMarkDeleteOutflow(actx); CHKERRQ(ierr);

	// count number of markers to be sent to each neighbor domain
	ierr = ADVMapMarkToDomains(actx); CHKERRQ(ierr);

	// communicate number of markers with neighbor processes
	ierr = ADVExchangeNumMark(actx); CHKERRQ(ierr);

	// create send and receive buffers for asynchronous MPI communication
	ierr = ADVCreateMPIBuff(actx); CHKERRQ(ierr);

	// communicate markers with neighbor processes
	ierr = ADVExchangeMark(actx); CHKERRQ(ierr);

	// store received markers, collect garbage
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// free communication buffer
	ierr = ADVDestroyMPIBuff(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVAdvectMark"
PetscErrorCode ADVAdvectMark(AdvCtx *actx)
{
	// update marker positions from current velocities & time step
	// WARNING! Forward Euler Explicit algorithm
	// (need to implement more accurate schemes)

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt    sx, sy, sz, nx, ny;
	PetscInt    jj, ID, I, J, K, II, JJ, KK, AirPhase;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar *ccx, *ccy, *ccz;
	PetscScalar ***lvx, ***lvy, ***lvz, ***lp, ***lT;
	PetscScalar vx, vy, vz, xc, yc, zc, xp, yp, zp, dt, Ttop;

	AirPhase = -1;
	Ttop     =  0.0;

	if(actx->surf->UseFreeSurf)
	{
		AirPhase = actx->surf->AirPhase;
		Ttop     = actx->jr->bc->Ttop;
	}

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = actx->fs;
	jr = actx->jr;

	// current time step
	dt = jr->ts->dt;

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

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

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

		// access host cell solution variables
		svCell = &jr->svCell[ID];

		// update pressure & temperature variables
		P->p += lp[sz+K][sy+J][sx+I] - svCell->svBulk.pn;
		P->T += lT[sz+K][sy+J][sx+I] - svCell->svBulk.Tn;

		// override temperature of air phase
		if(AirPhase != -1 && P->phase == AirPhase) P->T = Ttop;

		// advect marker
		P->X[0] = xp + vx*dt;
		P->X[1] = yp + vy*dt;
		P->X[2] = zp + vz*dt;

		// update displacement
		P->U[0] += vx*dt;
		P->U[1] += vy*dt;
		P->U[2] += vz*dt;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx, &lvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy, &lvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz, &lvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,  &lp);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMapMarkToDomains"
PetscErrorCode ADVMapMarkToDomains(AdvCtx *actx)
{
	// count number of markers to be sent to each neighbor domain

	PetscInt     i, lrank, cnt;
	PetscMPIInt  grank;
	FDSTAG      *fs;

	PetscErrorCode  ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// clear send counters
	ierr = PetscMemzero(actx->nsendm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	// scan markers
	for(i = 0, cnt = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		if(grank == -1)
		{
			// currently all the markers must remain in the box
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "ERROR! Marker outflow is currently not implemented!");

			// otherwise, number of deleted markers should be updated here, i.e.:
			// cnt++;

			// we should also decide what to do with the outflow markers
			// periodic boundary conditions?
		}
		else if(grank != actx->iproc)
		{
			// count markers that should be sent to each neighbor
			actx->nsendm[lrank]++;
			cnt++;
		}
	}

	// store number of deleted markers
	actx->ndel = cnt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVExchangeNumMark"
PetscErrorCode ADVExchangeNumMark(AdvCtx *actx)
{
	// communicate number of markers with neighbor processes
	FDSTAG     *fs;
	PetscInt    k;
	PetscMPIInt scnt, rcnt;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// zero out message counters
	scnt = 0;
	rcnt = 0;

	// send number of markers to ALL neighbor processes (except self & non-existing)
	for(k = 0; k < _num_neighb_; k++)
	{
		if(fs->neighb[k] != actx->iproc && fs->neighb[k] != -1)
		{
			ierr = MPI_Isend(&actx->nsendm[k], 1, MPIU_INT,
				fs->neighb[k], 100, actx->icomm, &srequest[scnt++]); CHKERRQ(ierr);
		}
	}

	// receive number of markers from ALL neighbor processes (except self & non-existing)
	for(k = 0; k < _num_neighb_; k++)
	{
		if(fs->neighb[k] != actx->iproc && fs->neighb[k] != -1)
		{
			ierr = MPI_Irecv(&actx->nrecvm[k], 1, MPIU_INT,
				fs->neighb[k], 100, actx->icomm, &rrequest[rcnt++]); CHKERRQ(ierr);
		}
		else actx->nrecvm[k] = 0;
	}

	// wait until all communication processes have been terminated
	if(scnt) { ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }
	if(rcnt) { ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCreateMPIBuff"
PetscErrorCode ADVCreateMPIBuff(AdvCtx *actx)
{
	// create send and receive buffers for asynchronous MPI communication

	// NOTE! Currently the memory allocation model is fully dynamic.
	// Maybe it makes sense to introduce static model with reallocation.
	FDSTAG     *fs;
	PetscInt    i, cnt, lrank;
	PetscMPIInt grank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// compute buffer pointers
	actx->nsend = getPtrCnt(_num_neighb_, actx->nsendm, actx->ptsend);
	actx->nrecv = getPtrCnt(_num_neighb_, actx->nrecvm, actx->ptrecv);

	actx->sendbuf = NULL;
	actx->recvbuf = NULL;
	actx->idel    = NULL;

	// allocate exchange buffers & array of deleted (sent) marker indices
	if(actx->nsend) { ierr = PetscMalloc((size_t)actx->nsend*sizeof(Marker),   &actx->sendbuf); CHKERRQ(ierr); }
	if(actx->nrecv) { ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker),   &actx->recvbuf); CHKERRQ(ierr); }
	if(actx->ndel)  { ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel);    CHKERRQ(ierr); }

	// copy markers to send buffer, store their indices
	for(i = 0, cnt = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		if(grank == -1)
		{
			// currently this situation is prohibited in ADVMapMarkersDomains
			// but we already handle it here anyway

			// delete marker from the storage
			actx->idel[cnt++] = i;

		}
		else if(grank != actx->iproc)
		{
			// store marker in the send buffer
			actx->sendbuf[actx->ptsend[lrank]++] = actx->markers[i];

			// delete marker from the storage
			actx->idel[cnt++] = i;
		}
	}

	// rewind send buffer pointers
	rewindPtr(_num_neighb_, actx->ptsend);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVExchangeMark"
PetscErrorCode ADVExchangeMark(AdvCtx *actx)
{
	// communicate markers with neighbor processes
	FDSTAG     *fs;
	PetscInt    k;
	PetscMPIInt scnt, rcnt, nbyte;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// zero out message counters
	scnt = 0;
	rcnt = 0;

	// send packages (if any) with markers to neighbor processes
	for(k = 0; k < _num_neighb_; k++)
	{
		if(actx->nsendm[k])
		{
			nbyte = (PetscMPIInt)(actx->nsendm[k]*(PetscInt)sizeof(Marker));

			ierr = MPI_Isend(&actx->sendbuf[actx->ptsend[k]], nbyte, MPI_BYTE,
				fs->neighb[k], 200, actx->icomm, &srequest[scnt++]); CHKERRQ(ierr);

		}
	}

	// receive packages (if any) with markers from neighbor processes
	for(k = 0; k < _num_neighb_; k++)
	{
		if(actx->nrecvm[k])
		{
			nbyte = (PetscMPIInt)(actx->nrecvm[k]*(PetscInt)sizeof(Marker));

			ierr = MPI_Irecv(&actx->recvbuf[actx->ptrecv[k]], nbyte, MPI_BYTE,
				fs->neighb[k], 200, actx->icomm, &rrequest[rcnt++]); CHKERRQ(ierr);
		}
	}

	// wait until all communication processes have been terminated
	if(scnt) { ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }
	if(rcnt) { ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVDestroyMPIBuff"
PetscErrorCode ADVDestroyMPIBuff(AdvCtx *actx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// destroy buffers
	ierr = PetscFree(actx->sendbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscFree(actx->idel);  	 CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCollectGarbage"
PetscErrorCode ADVCollectGarbage(AdvCtx *actx)
{
	// store received markers, collect garbage

	Marker   *markers, *recvbuf;
	PetscInt *idel, nummark, nrecv, ndel;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access storage
	nummark = actx->nummark;
	markers = actx->markers;

	nrecv   = actx->nrecv;
	recvbuf = actx->recvbuf;

	ndel    = actx->ndel;
	idel    = actx->idel;

	// close holes in marker storage
	while(nrecv && ndel)
	{
		markers[idel[ndel-1]] = recvbuf[nrecv-1];
		nrecv--;
		ndel--;
	}

	if(nrecv)
	{
		// make sure space is enough
		ierr = ADVReAllocStorage(actx, nummark + nrecv); CHKERRQ(ierr);

		// make sure we have a correct storage pointer
		markers = actx->markers;

		// put the rest in the end of marker storage
		while(nrecv)
		{
			markers[nummark++] = recvbuf[nrecv-1];
			nrecv--;
		}
	}

	if(ndel)
	{
		// collect garbage
		while(ndel)
		{
			if(idel[ndel-1] != nummark-1)
			{
				markers[idel[ndel-1]] = markers[nummark-1];
			}
			nummark--;
			ndel--;
		}
	}

	// store new number of markers
	actx->nummark = nummark;

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMapMarkToCells"
PetscErrorCode ADVMapMarkToCells(AdvCtx *actx)
{
	// computes local numbers of the host cells containing markers
	// NOTE: this routine MUST be called for the local markers only

	FDSTAG      *fs;
	PetscScalar *X;
	PetscInt     i, ID, I, J, K, M, N, P;

	PetscFunctionBegin;

	fs = actx->fs;

	// get number of cells
	M = fs->dsx.ncels;
	N = fs->dsy.ncels;
	P = fs->dsz.ncels;

	// loop over all local particles
	for(i = 0; i < actx->nummark; i++)
	{
		// get marker coordinates
		X = actx->markers[i].X;

		// find I, J, K indices by bisection algorithm
		I = FindPointInCell(fs->dsx.ncoor, 0, M, X[0]);
		J = FindPointInCell(fs->dsy.ncoor, 0, N, X[1]);
		K = FindPointInCell(fs->dsz.ncoor, 0, P, X[2]);

		// compute and store consecutive index
		GET_CELL_ID(ID, I, J, K, M, N);

		actx->cellnum[i] = ID;
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVUpdateMarkCell"
PetscErrorCode ADVUpdateMarkCell(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction

	FDSTAG      *fs;
	PetscInt    *numMarkCell, *m, i, p;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// allocate marker counter array
	ierr = makeIntArray(&numMarkCell, NULL, fs->nCells); CHKERRQ(ierr);

	// count number of markers in the cells
	for(i = 0; i < actx->nummark; i++) numMarkCell[actx->cellnum[i]]++;

	// store starting indices of markers belonging to a cell
	actx->markstart[0] = 0;
	for(i = 1; i < fs->nCells+1; i++) actx->markstart[i] = actx->markstart[i-1]+numMarkCell[i-1];

	// allocate memory for id offset
	ierr = makeIntArray(&m, NULL, fs->nCells); CHKERRQ(ierr);

	// store marker indices belonging to a cell
	for(i = 0; i < actx->nummark; i++)
	{
		p = actx->markstart[actx->cellnum[i]];
		actx->markind[p + m[actx->cellnum[i]]] = i;
		m[actx->cellnum[i]]++;
	}

	// free memory
	ierr = PetscFree(m);           CHKERRQ(ierr);
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCheckCorners"
PetscErrorCode ADVCheckCorners(AdvCtx *actx)
{
	// check corner marker distribution
	// if empty insert one marker in center of corner
	FDSTAG         *fs;
	BCCtx          *bc;
	Marker         *P, *markers;
	NumCorner      *numcorner;
	PetscInt       i, j, ii, jj, kk, ind, nx, ny, nz, lx[3];
	PetscInt       n, p, I, J, K;
	PetscInt       ID, ineigh, Ii, Ji, Ki, indcell[27], nummark[27];
	PetscInt       ninj = 0, nind = 0, sind = 0;
	PetscScalar    dx[3], xp[3], xc[3], xs[3], xe[3], *X;
	PetscScalar    x, sumind = 0.0;
	PetscScalar    cf_rand;
	PetscRandom    rctx;
	PetscLogDouble t0,t1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bc = actx->jr->bc;

	ierr = PetscTime(&t0); CHKERRQ(ierr);

	fs     = actx->fs;
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nz     = fs->dsz.ncels;

	// allocate memory for corners
	ierr = PetscMalloc((size_t)fs->nCells*sizeof(NumCorner), &numcorner); CHKERRQ(ierr);
	ierr = PetscMemzero(numcorner, (size_t)fs->nCells*sizeof(NumCorner)); CHKERRQ(ierr);

	// count number of markers in corners
	for(i = 0; i < fs->nCells; i++)
	{
		// get I, J, K cell indices
		GET_CELL_IJK(i, I, J, K, nx, ny);

		// get coordinates of cell center
		xc[0] = fs->dsx.ccoor[I];
		xc[1] = fs->dsy.ccoor[J];
		xc[2] = fs->dsz.ccoor[K];

		// load markers in cell
		n = actx->markstart[i+1] - actx->markstart[i];
		p = actx->markstart[i];

		for(ii = 0; ii < n; ii++)
		{
			P = &actx->markers[actx->markind[p+ii]];

			// get marker coordinates
			xp[0] = P->X[0];
			xp[1] = P->X[1];
			xp[2] = P->X[2];

			// check position in cell
			if ((xp[0] < xc[0]) && (xp[1] < xc[1]) && (xp[2] < xc[2])) numcorner[i].s[0]++; // LSW
			if ((xp[0] > xc[0]) && (xp[1] < xc[1]) && (xp[2] < xc[2])) numcorner[i].s[1]++; // LSE
			if ((xp[0] < xc[0]) && (xp[1] > xc[1]) && (xp[2] < xc[2])) numcorner[i].s[2]++; // LNW
			if ((xp[0] > xc[0]) && (xp[1] > xc[1]) && (xp[2] < xc[2])) numcorner[i].s[3]++; // LNE

			if ((xp[0] < xc[0]) && (xp[1] < xc[1]) && (xp[2] > xc[2])) numcorner[i].s[4]++; // USW
			if ((xp[0] > xc[0]) && (xp[1] < xc[1]) && (xp[2] > xc[2])) numcorner[i].s[5]++; // USE
			if ((xp[0] < xc[0]) && (xp[1] > xc[1]) && (xp[2] > xc[2])) numcorner[i].s[6]++; // UNW
			if ((xp[0] > xc[0]) && (xp[1] > xc[1]) && (xp[2] > xc[2])) numcorner[i].s[7]++; // UNE
		}
	}

	// count how much to allocate memory for new markers
	for(i = 0; i < fs->nCells; i++)
	{
		for(j = 0; j < 8; j++)
		{
			if (numcorner[i].s[j] == 0) ninj++;
		}
	}

	// if no need for new markers
	if (!ninj)
	{
		// clear memory
		ierr = PetscFree(numcorner); CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// initialize
	lx[0] = -1;
	lx[1] = 0;
	lx[2] = 1;

	// allocate memory for new markers
	actx->nrecv = ninj;
	ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->recvbuf, (size_t)actx->nrecv*sizeof(Marker)); CHKERRQ(ierr);

	// initialize the random number generator
	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

	// inject markers in corners
	for(i = 0; i < fs->nCells; i++)
	{
		for(j = 0; j < 8; j++)
		{
			if (numcorner[i].s[j] == 0)
			{
				// expand I, J, K cell indices
				GET_CELL_IJK(i, I, J, K, nx, ny);

				// get coordinates of cell center
				xc[0] = fs->dsx.ccoor[I];
				xc[1] = fs->dsy.ccoor[J];
				xc[2] = fs->dsz.ccoor[K];

				// get coordinates of cell node start
				xs[0] = fs->dsx.ncoor[I];
				xs[1] = fs->dsy.ncoor[J];
				xs[2] = fs->dsz.ncoor[K];

				// get coordinates of cell node end
				xe[0] = fs->dsx.ncoor[I+1];
				xe[1] = fs->dsy.ncoor[J+1];
				xe[2] = fs->dsz.ncoor[K+1];

				// GET MARKER NUMBER FROM LOCAL NEIGHBOURS
				ineigh = 0;
				n      = 0;
				for (kk = 0; kk < 3; kk++)
				{
					// get Ki index
					Ki = K + lx[kk];

					for (jj = 0; jj < 3; jj++)
					{
						// get Ji index
						Ji = J + lx[jj];

						for (ii = 0; ii < 3; ii++)
						{
							// get Ii index
							Ii = I + lx[ii];

							// check if within local domain bounds
							if ((Ki > -1) && (Ki < nz) && (Ji > -1) && (Ji < ny) && (Ii > -1) && (Ii < nx))
							{
								// get single index
								GET_CELL_ID(ID, Ii, Ji, Ki, nx, ny);

								// get markers number
								n += actx->markstart[ID+1] - actx->markstart[ID];

								// save markers number and single indices of cells
								indcell[ineigh] = ID;
								nummark[ineigh] = actx->markstart[ID+1] - actx->markstart[ID];
							}
							else
							{
								// save markers number and single indices of cells
								indcell[ineigh] = -1;
								nummark[ineigh] =  0;
							}

							// increase counter
							ineigh++;
						}
					}
				}

				// allocate memory for markers in cell
				ierr = PetscMalloc((size_t)n*sizeof(Marker),&markers); CHKERRQ(ierr);
				ierr = PetscMemzero(markers,(size_t)n*sizeof(Marker)); CHKERRQ(ierr);

				// load markers from all local neighbours
				jj = 0;
				for (ineigh = 0; ineigh < 27; ineigh++)
				{
					if ((indcell[ineigh]>-1) && (nummark[ineigh]>0))
					{
						for (ii = 0; ii < nummark[ineigh]; ii++)
						{
							// get index
							ind = actx->markind[actx->markstart[indcell[ineigh]] + ii];

							// save marker
							markers[jj] = actx->markers[ind];
							jj++;
						}
					}
				}

				// create coordinate of corner
				if (j == 0) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 1) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 2) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}
				if (j == 3) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xs[2]+xc[2])*0.5;}

				if (j == 4) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 5) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xs[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 6) { xp[0] = (xs[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}
				if (j == 7) { xp[0] = (xe[0]+xc[0])*0.5; xp[1] = (xe[1]+xc[1])*0.5; xp[2] = (xe[2]+xc[2])*0.5;}

				// add some random noise
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[0] += (cf_rand-0.5)*((xc[0]-xs[0])*0.5)*0.5;
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[1] += (cf_rand-0.5)*((xc[1]-xs[1])*0.5)*0.5;
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				xp[2] += (cf_rand-0.5)*((xc[2]-xs[2])*0.5)*0.5;

				// calculate the closest (parent marker)
				for (ii = 0; ii < n; ii++)
				{
					X  = markers[ii].X;
					dx[0] = X[0] - xp[0];
					dx[1] = X[1] - xp[1];
					dx[2] = X[2] - xp[2];
					x  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

					if (ii == 0)    { sumind = x; sind = ii; }
					if (x < sumind) { sumind = x; sind = ii; }
				}

				// create new marker
				actx->recvbuf[nind]      = markers[sind];
				actx->recvbuf[nind].X[0] = xp[0];
				actx->recvbuf[nind].X[1] = xp[1];
				actx->recvbuf[nind].X[2] = xp[2];

				// override marker phase (if necessary)
				ierr = BCOverridePhase(bc, i, actx->recvbuf + nind); CHKERRQ(ierr);

				// increase counter
				nind++;

				// free memory
				ierr = PetscFree(markers); CHKERRQ(ierr);
			}
		}
	}

	// destroy random context
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

	// store new markers
	actx->ndel = 0;
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// compute host cells for all the markers
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// update arrays for marker-cell interaction
	ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

	// print info
	ierr = PetscTime(&t1); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Marker control [%lld]: (Corners ) injected %lld markers in %1.4e s \n",(LLD)actx->iproc, (LLD)ninj, t1-t0);

	// clear
	ierr = PetscFree(numcorner);     CHKERRQ(ierr);
	ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkDeleteOutflow"
PetscErrorCode ADVMarkDeleteOutflow(AdvCtx *actx)
{
	// checks if markers are within the box bounds
	PetscInt     i, lrank, ndel;
	PetscMPIInt  grank;
	FDSTAG      *fs;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// scan markers
	for(i = 0, ndel = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		// count markers outside
		if(grank == -1) ndel++;
	}

	// if no need for deletion return
	if (!ndel) PetscFunctionReturn(0);

	// allocate storage
	actx->ndel  = ndel;
	ierr = PetscMalloc((size_t)actx->ndel *sizeof(PetscInt), &actx->idel); CHKERRQ(ierr);

	// save markers indices to be deleted
	for(i = 0, ndel = 0; i < actx->nummark; i++)
	{
		// get global & local ranks of a marker
		ierr = FDSTAGGetPointRanks(fs, actx->markers[i].X, &lrank, &grank); CHKERRQ(ierr);

		// save markers outside
		if(grank == -1) actx->idel[ndel++] = i;
	}

	// delete outside markers
	actx->nrecv   = 0;
	actx->recvbuf = NULL;
	ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

	// clear
	ierr = PetscFree(actx->idel);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVProjHistMarkToGrid"
PetscErrorCode ADVProjHistMarkToGrid(AdvCtx *actx)
{
	// Project the following history fields from markers to grid:

	// - phase ratios (centers and edges)
	// - pressure     (centers)
	// - temperature  (centers)
	// - APS          (centers and edges)
	// - stress       (centers or edges)
	// - displacement (centers)

	FDSTAG   *fs;
	JacRes   *jr;
	PetscInt  ii, jj, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs        = actx->fs;
	jr        = actx->jr;
	numPhases = actx->dbm->numPhases;

	// check marker phases
	ierr = ADVCheckMarkPhases(actx); CHKERRQ(ierr);

	//======
	// CELLS
	//======

	ierr = ADVInterpMarkToCell(actx); CHKERRQ(ierr);

	//======
	// EDGES
	//======

	// NOTE: edge phase ratio computation algorithm is the worst possible.
	// The xy, xz, yz edge points phase ratios are first computed locally,
	// and then assembled separately for each phase. This step involves
	// excessive communication, which is proportional to the number of phases.

	// *****************************************
	// SO PLEASE KEEP NUMBER OF PHASES MINIMIZED
	// *****************************************

	// compute edge phase ratios (consecutively)
	for(ii = 0; ii < numPhases; ii++)
	{
		ierr = ADVInterpMarkToEdge(actx, ii, _PHASE_); CHKERRQ(ierr);
	}

	// normalize phase ratios
	for(jj = 0; jj < fs->nXYEdg; jj++)  { ierr = getPhaseRatio(numPhases, jr->svXYEdge[jj].phRat, &jr->svXYEdge[jj].ws); CHKERRQ(ierr); }
	for(jj = 0; jj < fs->nXZEdg; jj++)  { ierr = getPhaseRatio(numPhases, jr->svXZEdge[jj].phRat, &jr->svXZEdge[jj].ws); CHKERRQ(ierr); }
	for(jj = 0; jj < fs->nYZEdg; jj++)  { ierr = getPhaseRatio(numPhases, jr->svYZEdge[jj].phRat, &jr->svYZEdge[jj].ws); CHKERRQ(ierr); }

	// interpolate history stress to edges
	ierr = ADVInterpMarkToEdge(actx, 0, _STRESS_); CHKERRQ(ierr);

	// interpolate plastic strain to edges
	ierr = ADVInterpMarkToEdge(actx, 0, _APS_); CHKERRQ(ierr);

	// update phase ratios taking into account actual free surface position
	ierr = FreeSurfGetAirPhaseRatio(actx->surf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVInterpMarkToCell"
PetscErrorCode ADVInterpMarkToCell(AdvCtx *actx)
{
	// marker-to-grid projection (cell nodes)

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     ii, jj, ID, I, J, K;
	PetscInt     nx, ny, nCells, numPhases;
	PetscScalar  xp, yp, zp, wxc, wyc, wzc, w = 0.0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs        = actx->fs;
	jr        = actx->jr;
	numPhases = actx->dbm->numPhases;

	// number of cells
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nCells = fs->nCells;

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
		svCell->svBulk.pn = 0.0;
		svCell->svBulk.Tn = 0.0;
		svCell->svDev.APS = 0.0;
		svCell->svDev.ATS = 0.0;
		svCell->hxx       = 0.0;
		svCell->hyy       = 0.0;
		svCell->hzz       = 0.0;
		svCell->U[0]      = 0.0;
		svCell->U[1]      = 0.0;
		svCell->U[2]      = 0.0;
	}

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables
		svCell->svBulk.pn += w*P->p;
		svCell->svBulk.Tn += w*P->T;
		svCell->svDev.APS += w*P->APS;
		svCell->svDev.ATS += w*P->ATS;
		svCell->hxx       += w*P->S.xx;
		svCell->hyy       += w*P->S.yy;
		svCell->hzz       += w*P->S.zz;
		svCell->U[0]      += w*P->U[0];
		svCell->U[1]      += w*P->U[1];
		svCell->U[2]      += w*P->U[2];

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		// normalize history variables
		svCell->svBulk.pn /= w;
		svCell->svBulk.Tn /= w;
		svCell->svDev.APS /= w;
		svCell->svDev.ATS /= w;
		svCell->hxx       /= w;
		svCell->hyy       /= w;
		svCell->hzz       /= w;
		svCell->U[0]      /= w;
		svCell->U[1]      /= w;
		svCell->U[2]      /= w;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVInterpMarkToEdge"
PetscErrorCode ADVInterpMarkToEdge(AdvCtx *actx, PetscInt iphase, InterpCase icase)
{
	// marker-to-grid projection (edge nodes)

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	PetscScalar  UPXY, UPXZ, UPYZ;
	PetscInt     nx, ny, sx, sy, sz;
	PetscInt     jj, ID, I, J, K, II, JJ, KK;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;
	PetscScalar  xc, yc, zc, xp, yp, zp, wxc, wyc, wzc, wxn, wyn, wzn;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;
	jr = actx->jr;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// clear local vectors
	ierr = VecZeroEntries(jr->ldxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldyz); CHKERRQ(ierr);

	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// set interpolated fields to defaults
	UPXY = 1.0; UPXZ = 1.0; UPYZ = 1.0;

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// perform phase ID test
		if(icase == _PHASE_ && P->phase != iphase) continue;

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get coordinates of cell center
		xc = fs->dsx.ccoor[I];
		yc = fs->dsy.ccoor[J];
		zc = fs->dsz.ccoor[K];

		// map marker on the control volumes of edge nodes
		if(xp > xc) { II = I+1; } else { II = I; }
		if(yp > yc) { JJ = J+1; } else { JJ = J; }
		if(zp > zc) { KK = K+1; } else { KK = K; }

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get interpolation weights in node control volumes
		wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);
		wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);
		wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);

		if      (icase == _STRESS_) { UPXY = P->S.xy; UPXZ = P->S.xz; UPYZ = P->S.yz; }
		else if (icase == _APS_)    { UPXY = P->APS;  UPXZ = P->APS;  UPYZ = P->APS;  }

		// update required fields from marker to edge nodes
		lxy[sz+K ][sy+JJ][sx+II] += wxn*wyn*wzc*UPXY;
		lxz[sz+KK][sy+J ][sx+II] += wxn*wyc*wzn*UPXZ;
		lyz[sz+KK][sy+JJ][sx+I ] += wxc*wyn*wzn*UPYZ;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// assemble global vectors
	LOCAL_TO_GLOBAL(fs->DA_XY, jr->ldxy, jr->gdxy)
	LOCAL_TO_GLOBAL(fs->DA_XZ, jr->ldxz, jr->gdxz)
	LOCAL_TO_GLOBAL(fs->DA_YZ, jr->ldyz, jr->gdyz)

	// access 1D layouts of global vectors
	ierr = VecGetArray(jr->gdxy, &gxy);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdxz, &gxz);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdyz, &gyz);  CHKERRQ(ierr);

	// copy (normalized) data to the residual context
	if(icase == _PHASE_)
	{
		for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].phRat[iphase] = gxy[jj];
		for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].phRat[iphase] = gxz[jj];
		for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].phRat[iphase] = gyz[jj];
	}
	else if(icase == _STRESS_)
	{
		for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].h = gxy[jj]/jr->svXYEdge[jj].ws;
		for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].h = gxz[jj]/jr->svXZEdge[jj].ws;
		for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].h = gyz[jj]/jr->svYZEdge[jj].ws;
	}
	else if(icase == _APS_)
	{
		for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].svDev.APS = gxy[jj]/jr->svXYEdge[jj].ws;
		for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].svDev.APS = gxz[jj]/jr->svXZEdge[jj].ws;
		for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].svDev.APS = gyz[jj]/jr->svYZEdge[jj].ws;
	}

	// restore access
	ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCheckMarkPhases"
PetscErrorCode ADVCheckMarkPhases(AdvCtx *actx)
{
	// check phases of markers
	Marker    *P;
	PetscInt  jj;
	PetscInt  numPhases;

	PetscFunctionBegin;

	numPhases = actx->dbm->numPhases;

	// scan all markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access marker
		P = &actx->markers[jj];

		// check marker phase
		if ((P->phase < 0) || (P->phase > numPhases-1))
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, " Detected markers with wrong phase! \n");
		}
	}

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
