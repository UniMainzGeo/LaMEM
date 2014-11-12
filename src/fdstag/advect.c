//---------------------------------------------------------------------------
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Utils.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "advect.h"
#include "constEq.h"

/*
#START_DOC#
\lamemfunction{\verb- ADVCreate -}
Create advection context

\lamemfunction{\verb- ADVDestroy -}
Destroy advection context

\lamemfunction{\verb- ADVAdvect -}
Main advection routine


\lamemfunction{\verb- ADVTestPoint -}
Check the location of the marker. Does it belong to an other cpu now ?


\lamemfunction{\verb- ADVCheckOwner -}
ADVCheckOwner is a subroutine used by ADVTestPoint to determine the cpu rank
of a particular point at coordinate (x,y,z)


\lamemfunction{\verb- ADVExchangeNumMarkers -}
Exchange the number of markers the will be sent in a next step.
This function call non-blocking \verb-MPI_Isend- and \verb-MPI_Irecv-.


\lamemfunction{\verb- ADVExchangeMarkers -}
This is the actual routine to exchange markers between cpus.
Its communication pattern is also based on non-blocking \verb-MPI_Isend- and \verb-MPI_Irecv-.

\lamemfunction{\verb- ADVCreateMPIBuffer -}
Here we create a dynamic \verb-MPI_buffer- depending on the number of markers that
have to be sent to each cpu, and that are going to be received from other cpus

\lamemfunction{\verb- ADVDestroyMPIBuffer -}
This destroy the allocated dynamic buffer space


\lamemfunction{\verb- ADVRemapping -}
Adina, please add a few comments

\lamemfunction{\verb- FDSTAGetVorticity -}
Anton, please add a few comments

#END_DOC#
*/
//---------------------------------------------------------------------------
// * add different advection methods (echo to output)
// * add different types of GRID->MARKER interpolation (echo to output)
//   (currently piece-wise constant, alternative - linear)
// * check weights of distance-dependent MARKER->GRID interpolation
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVClear"
PetscErrorCode ADVClear(AdvCtx *actx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(actx, sizeof(AdvCtx)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCreate"
PetscErrorCode ADVCreate(AdvCtx *actx, FDSTAG *fs, JacRes *jr)
{
	// create advection context

	PetscMPIInt nproc, iproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	actx->fs = fs;
	actx->jr = jr;

	//=============
	// COMMUNICATOR
	//=============

	ierr = MPI_Comm_dup(PETSC_COMM_WORLD, &actx->icomm); CHKERRQ(ierr);

	ierr = MPI_Comm_size(actx->icomm, &nproc); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(actx->icomm, &iproc); CHKERRQ(ierr);

	actx->nproc = (PetscInt)nproc;
	actx->iproc = (PetscInt)iproc;

	//========
	// STORAGE
	//========

	actx->nummark = 0;
	actx->markcap = 0;
	actx->markers = NULL;
	actx->cellnum = NULL;

	//=========
	// EXCHANGE
	//=========

	actx->sendbuf = NULL;
	actx->recvbuf = NULL;

	actx->nsend = 0;
	ierr = PetscMemzero(actx->nsendm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->ptsend, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	actx->nrecv = 0;
	ierr = PetscMemzero(actx->nrecvm, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->ptrecv, _num_neighb_*sizeof(PetscInt)); CHKERRQ(ierr);

	actx->ndel = 0;
	actx->idel = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVDestroy"
PetscErrorCode ADVDestroy(AdvCtx *actx)
{
	// destroy advection context

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_free(&actx->icomm); CHKERRQ(ierr);
	ierr = PetscFree(actx->markers);    CHKERRQ(ierr);
	ierr = PetscFree(actx->cellnum);    CHKERRQ(ierr);
	ierr = PetscFree(actx->sendbuf);    CHKERRQ(ierr);
	ierr = PetscFree(actx->recvbuf);    CHKERRQ(ierr);
	ierr = PetscFree(actx->idel);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVReAllocStorage"
PetscErrorCode ADVReAllocStorage(AdvCtx *actx, PetscInt nummark)
{
	// WARNING! This is a very crappy approach. Make sure the overhead is
	// large enough to prevent memory reallocations. Do marker management
	// before reallocating, or implement different memory model (e.g. paging,
	// or fixed maximum number markers per cell + deleting excessive markers.
	// The latter has an advantage of maintaining memory locality).

	PetscInt  markcap;
	Marker   *markers;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether current storage is insufficient
	if(nummark > actx->markcap)
	{
		// delete host cell numbers
		ierr = PetscFree(actx->cellnum); CHKERRQ(ierr);

		// compute new capacity
		markcap = (PetscInt)(_cap_overhead_*(PetscScalar)nummark);

		// allocate memory for markers
		ierr = PetscMalloc((size_t)markcap*sizeof(Marker), &markers); CHKERRQ(ierr);
		ierr = PetscMemzero(markers, (size_t)markcap*sizeof(Marker)); CHKERRQ(ierr);

		// copy current data
		if(actx->nummark)
		{
			ierr = PetscMemcpy(markers, actx->markers, (size_t)actx->nummark*sizeof(Marker)); CHKERRQ(ierr);
		}

		// delete previous marker storage
		ierr = PetscFree(actx->markers); CHKERRQ(ierr);

		// save new capacity & storage
		actx->markcap = markcap;
		actx->markers = markers;

		// allocate memory for host cell numbers
		ierr = PetscMalloc((size_t)markcap*sizeof(PetscInt), &actx->cellnum); CHKERRQ(ierr);
		ierr = PetscMemzero(actx->cellnum, (size_t)markcap*sizeof(PetscInt)); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVAdvect"
PetscErrorCode ADVAdvect(AdvCtx *actx)
{
	//=======================================================================
	// MAJOR ADVECTION ROUTINE
	//
	// WARNING!
	// Currently only implements Forward Euler Explicit algorithm.
	//=======================================================================

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// project history INCREMENTS from grid to markers
	ierr = ADVProjHistGridToMark(actx); CHKERRQ(ierr);

	// advect markers (Forward Euler)
	ierr = ADVAdvectMark(actx); CHKERRQ(ierr);

	// exchange markers between the processors
	ierr = ADVExchange(actx); CHKERRQ(ierr);

	// project advected history from markers back to grid
	ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVExchange"
PetscErrorCode ADVExchange(AdvCtx *actx)
{
	// exchange markers between the processors resulting from the position change

	PetscErrorCode ierr;
	PetscFunctionBegin;

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

	// compute host cells for all the markers received
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVProjHistGridToMark"
PetscErrorCode ADVProjHistGridToMark(AdvCtx *actx)
{

	// project history INCREMENTS from grid to markers

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = ADVInterpFieldToMark(actx, _APS_);       CHKERRQ(ierr);

	ierr = ADVInterpFieldToMark(actx, _STRESS_);    CHKERRQ(ierr);

	ierr = ADVInterpFieldToMark(actx, _VORTICITY_); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVInterpFieldToMark"
PetscErrorCode ADVInterpFieldToMark(AdvCtx *actx, InterpCase icase)
{
	//=======================================================================
	// interpolate increments of a history field to markers
	//
	// WARNING! vorticity case MUST be called after stress case (rotation)
	//=======================================================================

	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	Tensor2RN    R;
	Tensor2RS    SR;
	SolVarCell  *svCell;
	PetscScalar  UPXY, UPXZ, UPYZ;
	PetscInt     nx, ny, sx, sy, sz;
	PetscInt     jj, ID, I, J, K, II, JJ, KK;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;
	PetscScalar  xc, yc, zc, xp, yp, zp, wx, wy, wz, dt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;
	jr = actx->jr;

	// current time step
	dt = jr->ts.dt;

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// copy history increments into edge buffers
	if(icase == _VORTICITY_)
	{
		// compute current vorticity field
		ierr = JacResGetVorticity(jr); CHKERRQ(ierr);
	}
	else
	{
		// access 1D layouts of global vectors
		ierr = VecGetArray(jr->gdxy, &gxy);  CHKERRQ(ierr);
		ierr = VecGetArray(jr->gdxz, &gxz);  CHKERRQ(ierr);
		ierr = VecGetArray(jr->gdyz, &gyz);  CHKERRQ(ierr);

		if(icase == _STRESS_)
		{
			for(jj = 0; jj < fs->nXYEdg; jj++) gxy[jj] = jr->svXYEdge[jj].s - jr->svXYEdge[jj].h;
			for(jj = 0; jj < fs->nXZEdg; jj++) gxz[jj] = jr->svXZEdge[jj].s - jr->svXZEdge[jj].h;
			for(jj = 0; jj < fs->nYZEdg; jj++) gyz[jj] = jr->svYZEdge[jj].s - jr->svYZEdge[jj].h;
		}
		else if(icase == _APS_)
		{
			for(jj = 0; jj < fs->nXYEdg; jj++) gxy[jj] = jr->svXYEdge[jj].svDev.PSR;
			for(jj = 0; jj < fs->nXZEdg; jj++) gxz[jj] = jr->svXZEdge[jj].svDev.PSR;
			for(jj = 0; jj < fs->nYZEdg; jj++) gyz[jj] = jr->svYZEdge[jj].svDev.PSR;
		}

		// restore access
		ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
		ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
		ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

		// communicate boundary values
		GLOBAL_TO_LOCAL(fs->DA_XY, jr->gdxy, jr->ldxy);
		GLOBAL_TO_LOCAL(fs->DA_XZ, jr->gdxz, jr->ldxz);
		GLOBAL_TO_LOCAL(fs->DA_YZ, jr->gdyz, jr->ldyz);
	}

	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

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

		// get coordinates of cell center
		xc = fs->dsx.ccoor[I];
		yc = fs->dsy.ccoor[J];
		zc = fs->dsz.ccoor[K];

		// map marker on the control volumes of edge nodes
		if(xp > xc) { II = I+1; } else { II = I; }
		if(yp > yc) { JJ = J+1; } else { JJ = J; }
		if(zp > zc) { KK = K+1; } else { KK = K; }

		// access buffer
		UPXY = lxy[sz+K ][sy+JJ][sx+II];
		UPXZ = lxz[sz+KK][sy+J ][sx+II];
		UPYZ = lyz[sz+KK][sy+JJ][sx+I ];

		// access host cell solution variables
		svCell = &jr->svCell[ID];

		// update history fields on markers
		if(icase == _STRESS_)
		{
			P->S.xx += svCell->sxx - svCell->hxx;
			P->S.yy += svCell->syy - svCell->hyy;
			P->S.zz += svCell->szz - svCell->hzz;
			P->S.xy += UPXY;
			P->S.xz += UPXZ;
			P->S.yz += UPYZ;
		}
		else if(icase == _APS_)
		{
			P->APS += dt*sqrt(svCell->svDev.PSR + UPXY + UPXZ + UPYZ);
		}
		else if(icase == _VORTICITY_)
		{
			// interpret vorticity components
			wx = UPYZ;
			wy = UPXZ;
			wz = UPXY;

			// compute rotation matrix from local vorticity field
			GetRotationMatrix(&R, dt, wx, wy, wz);

			// rotate history stress
			RotateStress(&R, &P->S, &SR);

			// store rotated stress on the marker
			Tensor2RSCopy(&SR, &P->S);
		}
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

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
	PetscInt    jj, ID, I, J, K, II, JJ, KK;
	PetscScalar *ncx, *ncy, *ncz;
	PetscScalar *ccx, *ccy, *ccz;
	PetscScalar ***lvx, ***lvy, ***lvz, ***lp, ***lT;
	PetscScalar vx, vy, vz, p, T, xc, yc, zc, xp, yp, zp, dt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = actx->fs;
	jr = actx->jr;

	// current time step
	dt = jr->ts.dt;

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
		p  = InterpLin3D(lp,  II, JJ, KK, sx, sy, sz, xp, yp, zp, ccx, ccy, ccz);
		T  = InterpLin3D(lT,  II, JJ, KK, sx, sy, sz, xp, yp, zp, ccx, ccy, ccz);

		// access host cell solution variables
		svCell = &jr->svCell[ID];

		// update pressure & temperature variables
		P->p += p - svCell->svBulk.pn;
		P->T += T - svCell->svBulk.Tn;

		// advect marker
		P->X[0] = xp + vx*dt;
		P->X[1] = yp + vy*dt;
		P->X[2] = zp + vz*dt;
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
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Marker outflow is currently not implemented!");

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
#define __FUNCT__ "ADVProjHistMarkToGrid"
PetscErrorCode ADVProjHistMarkToGrid(AdvCtx *actx)
{
	// Project the following history fields from markers to grid:

	// - phase ratios (centers and edges)
	// - pressure     (centers)
	// - temperature  (centers)
	// - APS          (centers and edges)
	// - stress       (centers or edges)

	FDSTAG   *fs;
	JacRes   *jr;
	PetscInt  ii, jj;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;
	jr = actx->jr;

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
	for(ii = 0; ii < jr->numPhases; ii++)
	{
		ierr = ADVInterpMarkToEdge(actx, ii, _PHASE_); CHKERRQ(ierr);
	}

	// normalize phase ratios
	for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].ws = getPhaseRatio(jr->numPhases, jr->svXYEdge[jj].phRat);
	for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].ws = getPhaseRatio(jr->numPhases, jr->svXZEdge[jj].phRat);
	for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].ws = getPhaseRatio(jr->numPhases, jr->svYZEdge[jj].phRat);

	// interpolate history stress to edges
	ierr = ADVInterpMarkToEdge(actx, 0, _STRESS_); CHKERRQ(ierr);

	// interpolate plastic strain to edges
	ierr = ADVInterpMarkToEdge(actx, 0, _APS_); CHKERRQ(ierr);

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
	PetscInt     nx, ny, nCells;
	PetscScalar  xp, yp, zp, wxc, wyc, wzc, w;

	PetscFunctionBegin;

	fs = actx->fs;
	jr = actx->jr;

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
		for(ii = 0; ii < jr->numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
		svCell->svBulk.pn = 0.0;
		svCell->svBulk.Tn = 0.0;
		svCell->svDev.APS = 0.0;
		svCell->hxx       = 0.0;
		svCell->hyy       = 0.0;
		svCell->hzz       = 0.0;
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
		svCell->hxx       += w*P->S.xx;
		svCell->hyy       += w*P->S.yy;
		svCell->hzz       += w*P->S.zz;

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		w = getPhaseRatio(jr->numPhases, svCell->phRat);

		// normalize history variables
		svCell->svBulk.pn /= w;
		svCell->svBulk.Tn /= w;
		svCell->svDev.APS /= w;
		svCell->hxx       /= w;
		svCell->hyy       /= w;
		svCell->hzz       /= w;
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
// service functions
//-----------------------------------------------------------------------------
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[])
{
	// compute pointers from counts, return total count

	PetscInt i, tcnt = 0;

	for(i = 0; i < n; i++)
	{
		ptr[i] = tcnt;
		tcnt  += counts[i];
	}
	return tcnt;
}
//---------------------------------------------------------------------------
void rewindPtr(PetscInt n, PetscInt ptr[])
{
	// rewind pointers after using them as access iterators

	PetscInt i, prev = 0, next;

	for(i = 0; i < n; i++)
	{
		next   = ptr[i];
		ptr[i] = prev;
		prev   = next;
	}
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getPhaseRatio"
PetscScalar getPhaseRatio(PetscInt n, PetscScalar *v)
{
	// compute phase ratio array

	PetscInt    i;
	PetscScalar sum = 0.0;

	for(i = 0; i < n; i++) sum  += v[i];

	if(sum == 0.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, " Empty control volume");

	for(i = 0; i < n; i++) v[i] /= sum;

	return sum;
}
//-----------------------------------------------------------------------------
