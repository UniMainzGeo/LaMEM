//---------------------------------------------------------------------------
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Utils.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "advect.h"

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
/* Macro for distance-dependent interpolation from markers to edge nodes.
 * Parameters:
 * _TEST_                 - optional loop test
 * _UPXY_, _UPXZ_, _UPYZ_ - fields to be updated from marker to local vectors
 * _CP_                   - field to be copied from global vectors into context
 */
#define INTERP_MARKER_TO_EDGES(_TEST_, _UPXY_, _UPXZ_, _UPYZ_, _CP_) \
	/* clear local vectors */ \
	ierr = VecZeroEntries(jrctx->ldxy); CHKERRQ(ierr); \
	ierr = VecZeroEntries(jrctx->ldxz); CHKERRQ(ierr); \
	ierr = VecZeroEntries(jrctx->ldyz); CHKERRQ(ierr); \
	/* access 3D layouts of local vectors */ \
	ierr = DMDAVecGetArray(fs->DA_XY, jrctx->ldxy, &lxy); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(fs->DA_XZ, jrctx->ldxz, &lxz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(fs->DA_YZ, jrctx->ldyz, &lyz); CHKERRQ(ierr); \
	/* scan ALL markers*/ \
	for(jj = 0; jj < actx->nummark; jj++) \
	{	/* access next marker */ \
		P = &actx->markers[jj]; \
		/* perform optional loop test*/ \
		_TEST_ \
		/* get consecutive index of the host cell */ \
		ID = actx->cellnum[jj]; \
		/* expand I, J, K cell indices */ \
		GET_CELL_IJK(ID, Ic, Jc, Kc, nx, ny) \
		/* get marker coordinates */ \
		xp = P->X[0]; \
		yp = P->X[1]; \
		zp = P->X[2]; \
		/* get coordinates of cell center */ \
		xc = fs->dsx.ccoor[Ic]; \
		yc = fs->dsy.ccoor[Jc]; \
		zc = fs->dsz.ccoor[Kc]; \
		/* map marker on the control volumes of edge nodes */ \
		if(xp > xc) In = Ic+1; else In = Ic; \
		if(yp > yc) Jn = Jc+1; else Jn = Jc; \
		if(zp > zc) Kn = Kc+1; else Kn = Kc; \
		/* get interpolation weights in cell control volumes */ \
		wxc = WEIGHT_POINT_CELL(Ic, xp, fs->dsx); \
		wyc = WEIGHT_POINT_CELL(Jc, yp, fs->dsy); \
		wzc = WEIGHT_POINT_CELL(Kc, zp, fs->dsz); \
		/* get interpolation weights in node control volumes */ \
		wxn = WEIGHT_POINT_CELL(In, xp, fs->dsx); \
		wyn = WEIGHT_POINT_CELL(Jn, yp, fs->dsy); \
		wzn = WEIGHT_POINT_CELL(Kn, zp, fs->dsz); \
		/* update required fields from marker to edge nodes */ \
		lxy[Kc+sz][Jn+sy][In+sx] += wxn*wyn*wzc*_UPXY_; \
		lxz[Kn+sz][Jc+sy][In+sx] += wxn*wyc*wzn*_UPXZ_; \
		lyz[Kn+sz][Jn+sy][Ic+sx] += wxc*wyn*wzn*_UPYZ_; \
	} \
	/* restore access */ \
	ierr = DMDAVecRestoreArray(fs->DA_XY, jrctx->ldxy, &lxy); CHKERRQ(ierr); \
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jrctx->ldxz, &lxz); CHKERRQ(ierr); \
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jrctx->ldyz, &lyz); CHKERRQ(ierr); \
	/* assemble global vectors */ \
	LOCAL_TO_GLOBAL(fs->DA_XY, jrctx->ldxy, jrctx->gdxy) \
	LOCAL_TO_GLOBAL(fs->DA_XZ, jrctx->ldxz, jrctx->gdxz) \
	LOCAL_TO_GLOBAL(fs->DA_YZ, jrctx->ldyz, jrctx->gdyz) \
	/* access 1D layouts of global vectors */ \
	ierr = VecGetArray(jrctx->gdxy, &gxy);  CHKERRQ(ierr); \
	ierr = VecGetArray(jrctx->gdxz, &gxz);  CHKERRQ(ierr); \
	ierr = VecGetArray(jrctx->gdyz, &gyz);  CHKERRQ(ierr); \
	/* copy normalized data to residual context */ \
	for(jj = 0; jj < fs->nXYEdg; jj++) jrctx->svXYEdge[jj]._CP_ = gxy[jj]/jrctx->svXYEdge[jj].ws; \
	for(jj = 0; jj < fs->nXZEdg; jj++) jrctx->svXZEdge[jj]._CP_ = gxz[jj]/jrctx->svXZEdge[jj].ws; \
	for(jj = 0; jj < fs->nYZEdg; jj++) jrctx->svYZEdge[jj]._CP_ = gyz[jj]/jrctx->svYZEdge[jj].ws; \
	/* restore access */ \
	ierr = VecRestoreArray(jrctx->gdxy, &gxy); CHKERRQ(ierr); \
	ierr = VecRestoreArray(jrctx->gdxz, &gxz); CHKERRQ(ierr); \
	ierr = VecRestoreArray(jrctx->gdyz, &gyz); CHKERRQ(ierr);
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCreate"
PetscErrorCode ADVCreate(AdvCtx *actx)
{
	// create advection context

	PetscMPIInt nproc, iproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
#define __FUNCT__ "ADVReAllocateStorage"
PetscErrorCode ADVReAllocateStorage(AdvCtx *actx, PetscInt nummark)
{
	// WARNING! This is a very crappy approach. Make sure the overhead is
	// large enough to prevent memory reallocations. Do marker management
	// before reallocating, or implement different memory model (e.g. paging).

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
PetscErrorCode ADVAdvect(AdvCtx *actx, FDSTAG *fs)
{

	PetscErrorCode ierr;

	PetscFunctionBegin;

	// count number of markers to be sent to each neighbor domain
	ierr = ADVMapMarkersDomains(actx, fs); CHKERRQ(ierr);

	// communicate number of markers with neighbor processes
	ierr = ADVExchangeNumMarkers(actx, fs); CHKERRQ(ierr);

	// create send and receive buffers for asynchronous MPI communication
	ierr = ADVCreateMPIBuffer(actx, fs); CHKERRQ(ierr);

	// communicate markers with neighbor processes
	ierr = ADVExchangeMarkers(actx, fs); CHKERRQ(ierr);

	// store received markers, collect garbage
	ierr = ADVCollectGarbage(actx);

	// free communication buffer
	ierr = ADVDestroyMPIBuffer(actx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "ADVAdvectMarkers"
PetscErrorCode ADVAdvectMarkers(AdvCtx *actx, FDSTAG *fs, JacResCtx *jrctx)
{
//	PetscErrorCode ierr;
	PetscFunctionBegin;


	// map markers on control volumes of the X, Y, Z -face nodes

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVProjHistMarkGrid"
PetscErrorCode ADVProjHistMarkGrid(AdvCtx *actx, FDSTAG *fs, JacResCtx *jrctx)
{
	// Project the following history fields from markers to grid:

	// - phase ratios (centers and edges)
	// - pressure     (centers)
	// - temperature  (centers)
	// - APS          (centers and edges)
	// - stress       (centers or edges)

	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     nx, ny, nz, sx, sy, sz, nCells;
	PetscInt     ii, jj, ID, Ic, Jc, Kc, In, Jn, Kn;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;
	PetscScalar  xc, yc, zc, xp, yp, zp, wxc, wyc, wzc, wxn, wyn, wzn, w;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get number of cells
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	nCells = nx*ny*nz;

	//======
	// CELLS
	//======

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jrctx->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < jrctx->numPhases; ii++) svCell->phRat[ii] = 0.0;

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
		GET_CELL_IJK(ID, Ic, Jc, Kc, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(Ic, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(Jc, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(Kc, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jrctx->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables
		svCell->svBulk.pn += w*P->p;
		svCell->svBulk.Tn += w*P->T;
		svCell->svDev.APS += w*P->APS;
		svCell->hxx       += w*P->s.xx;
		svCell->hyy       += w*P->s.yy;
		svCell->hzz       += w*P->s.zz;

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jrctx->svCell[jj];

		// normalize phase ratios
		w = normVect(jrctx->numPhases, svCell->phRat);

		// normalize history variables
		svCell->svBulk.pn /= w;
		svCell->svBulk.Tn /= w;
		svCell->svDev.APS /= w;
		svCell->hxx       /= w;
		svCell->hyy       /= w;
		svCell->hzz       /= w;
	}

	//======
	// EDGES
	//======

	// NOTE: edge phase ratio computation algorithm is the worst possible.
	// The xy, xz, yz edge points phase ratios are first computed locally,
	// and then assembled separately for each phase. This step involves
	// excessive communication, which is proportional to the number of phases.

	// initialize sum of interpolation weights
	for(jj = 0; jj < fs->nXYEdg; jj++) jrctx->svXYEdge[jj].ws = 1.0;
	for(jj = 0; jj < fs->nXZEdg; jj++) jrctx->svXZEdge[jj].ws = 1.0;
	for(jj = 0; jj < fs->nYZEdg; jj++) jrctx->svYZEdge[jj].ws = 1.0;

	// define loop test for phase ratio calculation
	#define _TEST_ if(P->phase != ii) continue;

	// compute edge phase ratios (consecutively)
	for(ii = 0; ii < jrctx->numPhases; ii++)
	{
		INTERP_MARKER_TO_EDGES(_TEST_, 1.0, 1.0, 1.0, phRat[ii])
	}

	// normalize phase ratios
	for(jj = 0; jj < fs->nXYEdg; jj++) jrctx->svXYEdge[jj].ws = normVect(jrctx->numPhases, jrctx->svXYEdge[jj].phRat);
	for(jj = 0; jj < fs->nXZEdg; jj++) jrctx->svXZEdge[jj].ws = normVect(jrctx->numPhases, jrctx->svXZEdge[jj].phRat);
	for(jj = 0; jj < fs->nYZEdg; jj++) jrctx->svYZEdge[jj].ws = normVect(jrctx->numPhases, jrctx->svYZEdge[jj].phRat);

	// clear loop test
	#undef  _TEST_
	#define _TEST_

	// interpolate history stress to edges
	INTERP_MARKER_TO_EDGES(_TEST_, P->s.xy, P->s.xz, P->s.yz, h)

	// interpolates plastic strain to edges
	INTERP_MARKER_TO_EDGES(_TEST_, P->APS, P->APS, P->APS, svDev.APS)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMapMarkersDomains"
PetscErrorCode ADVMapMarkersDomains(AdvCtx *actx, FDSTAG *fs)
{
	// count number of markers to be sent to each neighbor domain

	PetscInt    i, lrank, cnt;
	PetscMPIInt grank;

	PetscErrorCode  ierr;
	PetscFunctionBegin;

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
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Marker outflow is currently not implemented!\n");

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
#define __FUNCT__ "ADVExchangeNumMarkers"
PetscErrorCode ADVExchangeNumMarkers(AdvCtx *actx, FDSTAG *fs)
{
	// communicate number of markers with neighbor processes

	PetscInt    k;
	PetscMPIInt scnt, rcnt;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
	ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr);
	ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVCreateMPIBuffer"
PetscErrorCode ADVCreateMPIBuffer(AdvCtx *actx, FDSTAG *fs)
{
	// create send and receive buffers for asynchronous MPI communication

	// NOTE! Currently the memory allocation model is fully dynamic.
	// Maybe it makes sense to introduce static model with reallocation.

	PetscInt    i, cnt, lrank;
	PetscMPIInt grank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// compute buffer pointers
	actx->nsend = getPtrCnt(_num_neighb_, actx->nsendm, actx->ptsend);
	actx->nrecv = getPtrCnt(_num_neighb_, actx->nrecvm, actx->ptrecv);

	// allocate exchange buffers
	ierr = PetscMalloc((size_t)actx->nsend*sizeof(Marker), &actx->sendbuf); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);

	// allocate array of the deleted (sent) marker indices
	ierr = PetscMalloc((size_t)actx->ndel*sizeof(PetscInt), &actx->idel); CHKERRQ(ierr);

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
#define __FUNCT__ "ADVExchangeMarkers"
PetscErrorCode ADVExchangeMarkers(AdvCtx *actx, FDSTAG *fs)
{
	// communicate markers with neighbor processes

	PetscInt    k;
	PetscMPIInt scnt, rcnt, nbyte;
	MPI_Request srequest[_num_neighb_];
	MPI_Request rrequest[_num_neighb_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
	ierr = MPI_Waitall(scnt, srequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr);
	ierr = MPI_Waitall(rcnt, rrequest, MPI_STATUSES_IGNORE); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVDestroyMPIBuffer"
PetscErrorCode ADVDestroyMPIBuffer(AdvCtx *actx)
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
		ierr = ADVReAllocateStorage(actx, nummark + nrecv); CHKERRQ(ierr);

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
#define __FUNCT__ "ADVMapMarkersCells"
PetscErrorCode ADVMapMarkersCells(AdvCtx *actx, FDSTAG *fs)
{
	// computes local numbers of the host cells containing markers
	// NOTE: this routine MUST be called for the local markers only

	PetscScalar *X;
	PetscInt     i, ID, I, J, K, M, N, P;

	PetscFunctionBegin;

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
#define __FUNCT__ "FDSTAGetVorticity"
PetscErrorCode FDSTAGetVorticity(
	FDSTAG *fs,
	Vec lvx,  Vec lvy,  Vec lvz, // local (ghosted) velocities
	Vec gwx,  Vec gwy,  Vec gwz) // global vorticity components
{
	// Compute components of vorticity vector
	// (instantaneous rotation rate around three coordinate axis).
	// Take care of rotation direction and sign convention.
	// Throughout LaMEM, right-handed coordinate system is assumed!

	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar dvxdy, dvydx, dvxdz, dvzdx, dvydz, dvzdy;
	PetscScalar ***vx, ***vy, ***vz;
	PetscScalar ***wx, ***wy, ***wz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  gwz,  &wz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  gwy,  &wy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  gwx,  &wx);  CHKERRQ(ierr);

	//-------------------------------
	// xy edge points (wz)
	//-------------------------------

	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvxdy = (vx[k][j][i] - vx[k][j-1][i])/SIZE_NODE(j, sy, fs->dsy);
		dvydx = (vy[k][j][i] - vy[k][j][i-1])/SIZE_NODE(i, sx, fs->dsx);

		// positive (counter-clockwise) rotation around Z axis X -> Y

		wz[k][j][i] = dvydx - dvxdy;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (wy)
	//-------------------------------

	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvxdz = (vx[k][j][i] - vx[k-1][j][i])/SIZE_NODE(k, sz, fs->dsz);
		dvzdx = (vz[k][j][i] - vz[k][j][i-1])/SIZE_NODE(i, sx, fs->dsx);

		// positive (counter-clockwise) rotation around Y axis Z -> X

		wy[k][j][i] = dvxdz - dvzdx;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (wx)
	//-------------------------------

	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		dvydz = (vy[k][j][i] - vy[k-1][j][i])/SIZE_NODE(k, sz, fs->dsz);
		dvzdy = (vz[k][j][i] - vz[k][j-1][i])/SIZE_NODE(j, sy, fs->dsy);

		// positive (counter-clockwise) rotation around X axis Y -> Z

		wx[k][j][i] = dvzdy - dvydz;
	}
	END_STD_LOOP

	// restore velocity & strain rate component vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   lvx,  &vx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   lvy,  &vy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   lvz,  &vz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  gwz,  &wz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  gwy,  &wy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  gwx,  &wx);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
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
//---------------------------------------------------------------------------


/*
	for(ii = 0; ii < jrctx->numPhases; ii++)
	{
		// clear
		ierr = VecZeroEntries(jrctx->ldxy); CHKERRQ(ierr);
		ierr = VecZeroEntries(jrctx->ldxz); CHKERRQ(ierr);
		ierr = VecZeroEntries(jrctx->ldyz); CHKERRQ(ierr);

		// access 3D layouts of local vectors
		ierr = DMDAVecGetArray(fs->DA_XY, jrctx->ldxy, &lxy); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_XZ, jrctx->ldxz, &lxz); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_YZ, jrctx->ldyz, &lyz); CHKERRQ(ierr);

		// compute local contributions
		for(jj = 0; jj < actx->nummark; jj++)
		{
			// skip the marker which phase number is different than ii
			if(actx->markers[jj].phase != ii) continue;

			// get consecutive index of the host cell
			ID = actx->cellnum[jj];

			// expand I, J, K cell indices
			GET_CELL_IJK(ID, Ic, Jc, Kc, M, N)

			// get marker coordinates
			xm  = actx->markers[jj].X[0];
			ym  = actx->markers[jj].X[1];
			zm  = actx->markers[jj].X[2];

			// get coordinates of cell center
			xc = fs->dsx.ccoor[Ic];
			yc = fs->dsy.ccoor[Jc];
			zc = fs->dsz.ccoor[Kc];

			// map marker on the control volumes of edge nodes
			if(xm > xc) In = Ic+1; else In = Ic;
			if(ym > yc) Jn = Jc+1; else Jn = Jc;
			if(zm > zc) Kn = Kc+1; else Kn = Kc;

			// get interpolation weights in cell control volumes
			wxc = WEIGHT_POINT_CELL(Ic, xm, fs->dsx);
			wyc = WEIGHT_POINT_CELL(Jc, ym, fs->dsy);
			wzc = WEIGHT_POINT_CELL(Kc, zm, fs->dsz);

			// get interpolation weights in node control volumes
			wxn = WEIGHT_POINT_CELL(In, xm, fs->dsx);
			wyn = WEIGHT_POINT_CELL(Jn, ym, fs->dsy);
			wzn = WEIGHT_POINT_CELL(Kn, zm, fs->dsz);

			// update phase ratios
			lxy[In][Jn][Kc] += wxn*wyn*wzc;
			lxz[In][Jc][Kn] += wxn*wyc*wzn;
			lyz[Ic][Jn][Kn] += wxc*wyn*wzn;
		}

		// restore access
		ierr = DMDAVecRestoreArray(fs->DA_XY, jrctx->ldxy, &lxy); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_XZ, jrctx->ldxz, &lxz); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_YZ, jrctx->ldyz, &lyz); CHKERRQ(ierr);

		// assemble
		LOCAL_TO_GLOBAL(fs->DA_XY, jrctx->ldxy, jrctx->gdxy)
		LOCAL_TO_GLOBAL(fs->DA_XZ, jrctx->ldxz, jrctx->gdxz)
		LOCAL_TO_GLOBAL(fs->DA_YZ, jrctx->ldyz, jrctx->gdyz)

		// access 1D layouts of global vectors
		ierr = VecGetArray(jrctx->gdxy, &gxy);  CHKERRQ(ierr);
		ierr = VecGetArray(jrctx->gdxz, &gxz);  CHKERRQ(ierr);
		ierr = VecGetArray(jrctx->gdyz, &gyz);  CHKERRQ(ierr);

		// copy data to residual context
		for(jj = 0; jj < fs->nXYEdg; jj++) jrctx->svXYEdge[jj].phRat[ii] = gxy[jj];
		for(jj = 0; jj < fs->nXZEdg; jj++) jrctx->svXZEdge[jj].phRat[ii] = gxz[jj];
		for(jj = 0; jj < fs->nYZEdg; jj++) jrctx->svYZEdge[jj].phRat[ii] = gyz[jj];

		// restore access
		ierr = VecRestoreArray(jrctx->gdxy, &gxy); CHKERRQ(ierr);
		ierr = VecRestoreArray(jrctx->gdxz, &gxz); CHKERRQ(ierr);
		ierr = VecRestoreArray(jrctx->gdyz, &gyz); CHKERRQ(ierr);
	}
*/
