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
//........   PARALLEL STAGGERED GRID USING PETSC DISTRIBUTED ARRAYS  ........
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "parsing.h"
#include "scaling.h"
#include "tools.h"

//---------------------------------------------------------------------------
// MeshSeg1D functions
//---------------------------------------------------------------------------
PetscErrorCode MeshSeg1DReadParam(
	MeshSeg1D  *ms,
	PetscScalar leng,
	PetscScalar gtol,
	const char *dir,
	FB         *fb)
{
	PetscInt    i, tcels, uniform;
	PetscInt    ncells[_max_num_segs_];
	PetscScalar avgsz, sz;
	char        *nseg, *nel, *coord, *bias, *periodic;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize
	ierr = PetscMemzero(ms, sizeof(MeshSeg1D)); CHKERRQ(ierr);

	ms->nsegs = 1;

	for(i = 0; i < _max_num_segs_; i++)
	{
		ms->biases[i] = 1.0;
		ncells    [i] = 0.0;
	}

	// compose option keys
	asprintf(&nseg,     "nseg_%s",     dir);
	asprintf(&nel,      "nel_%s",      dir);
	asprintf(&coord,    "coord_%s",    dir);
	asprintf(&bias,     "bias_%s",     dir);
	asprintf(&periodic, "periodic_%s", dir);

	// read parameters
	ierr = getIntParam   (fb, _OPTIONAL_, nseg,     &ms->nsegs,    1,           _max_num_segs_);  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, nel,       ncells,       ms->nsegs,   _max_num_cells_); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, periodic, &ms->periodic, 1,           1);               CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, coord,     ms->xstart,   ms->nsegs+1, leng);            CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, bias,      ms->biases,   ms->nsegs,   1.0 );            CHKERRQ(ierr);

	// compute starting node indices
	for(i = 0, tcels = 0; i < ms->nsegs; i++)
	{
		ms->istart[i] = tcels;
		tcels        += ncells[i];
	}
	ms->istart[ms->nsegs] = tcels;

	// check total number of cells
	if(tcels < 2)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Less than two cells are specified in the %s - direction\n", dir);
	}

	// check ordering and bias factors
	for(i = 0; i < ms->nsegs; i++)
	{
		if(ms->xstart[i] >= ms->xstart[i+1])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unordered coordinates in parameter %s\n", coord);
		}
		if(ms->biases[i] != 1.0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Non-unit bias ratios are deprecated (%s)\n", bias);
		}
	}

	// check for uniform grid
	uniform = 1;
	avgsz   = (ms->xstart[ms->nsegs] - ms->xstart[0])/(PetscScalar)tcels;

	for(i = 0; i < ms->nsegs; i++)
	{
		sz = (ms->xstart[i+1] - ms->xstart[i])/(PetscScalar)ncells[i];

		if(ms->biases[i] != 1.0 || PetscAbsScalar(avgsz-sz) > gtol*avgsz)
		{
			uniform = 0; break;
		}
	}

	// set grid parameters
	ms->tcels   = tcels;
	ms->uniform = uniform;

	// free keys
	free(nseg);
	free(nel);
	free(coord);
	free(bias);
	free(periodic);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MeshSeg1DGenCoord(
	MeshSeg1D   *ms,     // segments description
	PetscInt     iseg,   // segment index
	PetscInt     nl,     // number of nodes to be generated
	PetscInt     istart, // index of the first node
	PetscScalar *crd)    // coordinates of the nodes
{
	// (partially) mesh a segment with (optionally) biased element size

	PetscInt    i, N, M, sum;
	PetscScalar xstart, xclose, bias, avgSz, begSz, endSz, dx;

	PetscFunctionBeginUser;

	// total number of nodes in segment (including both ends)
	N = ms->istart[iseg+1] - ms->istart[iseg] + 1;

	// total number of cells
	M = N-1;

	// starting & closing coordinates
	xstart = ms->xstart[iseg];
	xclose = ms->xstart[iseg+1];

	// bias (last to first cell size ratio > 1 -> growing)
	bias = ms->biases[iseg];

	// average cell size
	avgSz = (xclose - xstart)/(PetscScalar)M;

	// uniform case
	if(bias == 1.0)
	{
		// generate coordinates of local nodes
		for(i = 0; i < nl; i++)
		{
			crd[i] = xstart + (PetscScalar)(istart + i)*avgSz;
		}
	}
	// non-uniform case
	else
	{
		// cell size limits
		begSz = 2.0*avgSz/(1.0 + bias);
		endSz = bias*begSz;

		// cell size increment (negative for bias < 1)
		dx = (endSz - begSz)/(PetscScalar)(M-1);

		// get accumulated sum of increments
		for(i = 0, sum = 0; i < istart; i++) sum += i;

		// generate coordinates of local nodes
		for(i = 0; i < nl; i++)
		{
			crd[i] = xstart + (PetscScalar)(istart + i)*begSz + (PetscScalar)sum*dx;
			sum += istart + i;
		}
	}

	// override last node coordinate
	if(istart+nl == N) crd[nl-1] = xclose;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Discret1D functions
//---------------------------------------------------------------------------
PetscErrorCode Discret1DCreate(
		Discret1D  *ds,
		PetscInt    nproc,     // number of processors
		PetscInt    rank,      // processor rank
		PetscInt   *nnodProc,  // number of nodes per processor
		PetscInt    color,     // column color
		PetscMPIInt grprev,    // global rank of previous process
		PetscMPIInt grnext,    // global rank of next process
		PetscScalar gtol,      // geometric tolerance
		const char *dir)       // direction label
{
	PetscInt i, cnt;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize
	ierr = PetscMemzero(ds, sizeof(Discret1D)); CHKERRQ(ierr);

	// number of processors
	ds->nproc = nproc;

	// rank of current processor
	ds->rank = (PetscMPIInt)rank;

	// index of first node (cell) on all processors + last index
	ierr = makeIntArray(&ds->starts, 0, nproc+1); CHKERRQ(ierr);

	for(i = 0, cnt = 0; i < nproc; i++)
	{	ds->starts[i] = cnt;
		cnt          += nnodProc[i];
	}
	ds->starts[nproc] = cnt-1;

	// index of first node (cell) on this processors
	ds->pstart = ds->starts[ds->rank];

	// total number of nodes
	ds->tnods = cnt;

	// total number of cells
	ds->tcels = cnt-1;

	// NUMBER OF NODES / CELLS

	// number of local nodes
	ds->nnods = nnodProc[rank];

	// number of local cells
	if(grnext != -1) ds->ncels = nnodProc[rank];
	else             ds->ncels = nnodProc[rank] - 1;

	// check number of cells in local grid
	if(ds->ncels < 2)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Less than two local cells in the %s - direction\n", dir);
	}

	// coordinates of local nodes + 1 layer (left) & 2 layers (right) of ghost points
	// NOTE: on the last processor there is only one ghost point from the right

	ierr = makeScalArray(&ds->nbuff, 0, ds->ncels+3); CHKERRQ(ierr);
	ds->ncoor = ds->nbuff + 1;

	// coordinates of local cells + 1 layer (both sides) of ghost points
	ierr = makeScalArray(&ds->cbuff, 0, ds->ncels+2); CHKERRQ(ierr);
	ds->ccoor = ds->cbuff + 1;

	// global rank of previous process (-1 if none)
	ds->grprev = grprev;

	// global rank of next process (-1 if none)
	ds->grnext = grnext;

	// column color
	ds->color = (PetscMPIInt) color;

	// column communicator
	ds->comm = MPI_COMM_NULL;

	// geometric tolerance
	ds->gtol = gtol;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DDestroy(Discret1D *ds)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// free memory buffers
	ierr = PetscFree(ds->nbuff);        CHKERRQ(ierr);
	ierr = PetscFree(ds->cbuff);        CHKERRQ(ierr);
	ierr = PetscFree(ds->starts);       CHKERRQ(ierr);
	ierr = Discret1DFreeColumnComm(ds); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DReadRestart(Discret1D *ds, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = makeIntArray (&ds->starts, NULL, ds->nproc + 1); CHKERRQ(ierr);
	ierr = makeScalArray(&ds->nbuff,  NULL, ds->ncels + 3); CHKERRQ(ierr);
	ierr = makeScalArray(&ds->cbuff,  NULL, ds->ncels + 2); CHKERRQ(ierr);

   	fread(ds->starts, sizeof(PetscInt   )*(size_t)(ds->nproc + 1), 1, fp);
	fread(ds->nbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 3), 1, fp);
	fread(ds->cbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 2), 1, fp);

	ds->ncoor = ds->nbuff + 1;
	ds->ccoor = ds->cbuff + 1;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DWriteRestart(Discret1D *ds, FILE *fp)
{
	PetscFunctionBeginUser;

	fwrite(ds->starts, sizeof(PetscInt   )*(size_t)(ds->nproc + 1), 1, fp);
	fwrite(ds->nbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 3), 1, fp);
	fwrite(ds->cbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 2), 1, fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DGetNumCells(Discret1D *ds, PetscInt **ncelProc)
{
	// get number of cells per processor

	PetscInt i, *l;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = makeIntArray(&l, NULL, ds->nproc); CHKERRQ(ierr);

	for(i = 0; i < ds->nproc; i++)
	{
		l[i] = ds->starts[i+1] - ds->starts[i];
	}

	(*ncelProc) = l;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DGenCoord(Discret1D *ds, MeshSeg1D *ms)
{
	PetscInt     i, n, nl, pstart, istart;
	PetscScalar *crd;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// compute number of nodes to be generated locally
	pstart = ds->pstart;
	crd    = ds->ncoor;
	n      = ds->nnods;

	// correct numbers if we need to include internal ghost points
	if(ds->grprev != -1) { pstart--; crd--; n++; }
	if(ds->grnext != -1) { n += 2; }

	// apply local filter & expand segment data
	for(i = 0; n; i++)
	{
		// compute number of nodes within this segment
		nl = ms->istart[i+1] - pstart + 1;

		// skip the non-overlapping segments
		if(nl < 0) continue;

		// correct if the rest of the mesh completely fits into the segment
		if(nl > n) nl = n;

		// compute starting index within the segment
		istart = pstart - ms->istart[i];

		// generate nodal coordinates for the local part of the segment
		ierr = MeshSeg1DGenCoord(ms, i, nl, istart, crd); CHKERRQ(ierr);

		// update the rest of the local mesh to be generated
		pstart += nl;
		crd    += nl;
		n      -= nl;
	}

	// generate ghost points and cell center coordinates
	ierr = Discret1DCompleteCoord(ds); CHKERRQ(ierr);

	// set uniform grid flag
	ds->uniform = ms->uniform;

	// set periodic periodic topology flag
	ds->periodic = ms->periodic;

	// set global grid coordinate bounds
	ds->gcrdbeg = ms->xstart[0];
	ds->gcrdend = ms->xstart[ms->nsegs];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DCoarsenCoord(Discret1D *coarse, Discret1D *fine)
{
	PetscInt    i, nn;
	PetscMPIInt cnt;
	MPI_Request request[4];
	PetscScalar sprev, rprev, snext, rnext;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// copy data
	coarse->uniform  = fine->uniform;  // uniform grid flag
	coarse->periodic = fine->periodic; // periodic topology flag
	coarse->gcrdbeg  = fine->gcrdbeg;  // global grid coordinate bound (begin)
	coarse->gcrdend  = fine->gcrdend;  // global grid coordinate bound (end)

	// check whether mesh is coarsened
	if(coarse->ncels == fine->ncels)
	{
		// copy node coordinate buffer
		for(i = 0, nn = fine->ncels+3; i < nn; i++)
			coarse->nbuff[i] = fine->nbuff[i];

		// copy cell coordinate buffer
		for(i = 0, nn = fine->ncels+2; i < nn; i++)
			coarse->cbuff[i] = fine->cbuff[i];
	}
	else
	{
		// get number of coarse local nodes
		nn = coarse->ncels + 1;

		// store coarse local node coordinates
		for(i = 0; i < nn; i++)
			coarse->ncoor[i] = fine->ncoor[2*i];

		// exchange ghost point coordinates
		cnt = 0;

		if(fine->grprev != -1)
		{
			sprev = fine->ncoor[2];

			ierr = MPI_Isend(&sprev, 1, MPIU_SCALAR, fine->grprev, 700, PETSC_COMM_WORLD, &request[cnt++]); CHKERRQ(ierr);
			ierr = MPI_Irecv(&rprev, 1, MPIU_SCALAR, fine->grprev, 700, PETSC_COMM_WORLD, &request[cnt++]); CHKERRQ(ierr);
		}

		if(fine->grnext != -1)
		{
			snext = fine->ncoor[fine->ncels-2];

			ierr = MPI_Isend(&snext, 1, MPIU_SCALAR, fine->grnext, 700, PETSC_COMM_WORLD, &request[cnt++]); CHKERRQ(ierr);
			ierr = MPI_Irecv(&rnext, 1, MPIU_SCALAR, fine->grnext, 700, PETSC_COMM_WORLD, &request[cnt++]); CHKERRQ(ierr);
		}

		// wait until all communication processes have been terminated
		if(cnt) { ierr = MPI_Waitall(cnt, request, MPI_STATUSES_IGNORE); CHKERRQ(ierr); }

		if(fine->grprev != -1) { coarse->ncoor[-1] = rprev; }
		if(fine->grnext != -1) { coarse->ncoor[nn] = rnext; }

		// generate ghost points and cell center coordinates
		ierr = Discret1DCompleteCoord(coarse); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DCompleteCoord(Discret1D *ds)
{
	// generate ghost points and cell center coordinates

	PetscInt    i;
	PetscScalar A, B, C;

	PetscFunctionBeginUser;

	// set boundary ghost coordinates
	if(ds->grprev == -1)
	{
		A = ds->ncoor[0];
		B = ds->ncoor[1];
		C = A - (B - A);
		ds->ncoor[-1] = C;
	}
	if(ds->grnext == -1)
	{
		A = ds->ncoor[ds->nnods-2];
		B = ds->ncoor[ds->nnods-1];
		C = B + (B - A);
		ds->ncoor[ds->nnods] = C;
	}

	// compute coordinates of the cell centers including ghosts
	for(i = -1; i < ds->ncels+1; i++)
		ds->ccoor[i] = (ds->ncoor[i] + ds->ncoor[i+1])/2.0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DStretch(Discret1D *ds, PetscScalar eps, PetscScalar ref)
{
	// stretch grid with constant stretch factor about reference point
	// x_new = x_old + eps*(x_old - x_ref)

	PetscInt i;

	PetscFunctionBeginUser;

	// recompute (stretch) node coordinates in the buffer
	for(i = 0; i < ds->ncels + 3; i++) ds->nbuff[i] += eps*(ds->nbuff[i] - ref);

	// recompute cell coordinates
	for(i = -1; i < ds->ncels+1; i++)
		ds->ccoor[i] = (ds->ncoor[i] + ds->ncoor[i+1])/2.0;

	// recompute global coordinate bounds
	ds->gcrdbeg *= (1.0 + eps);
	ds->gcrdend *= (1.0 + eps);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DGetColumnComm(Discret1D *ds)
{
	// This function is called every time the column communicator is needed.
	// Nothing is done if communicator already exists or in sequential case.

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(ds->nproc != 1 && ds->comm == MPI_COMM_NULL)
	{
		ierr = MPI_Comm_split(PETSC_COMM_WORLD, ds->color, ds->rank, &ds->comm); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DFreeColumnComm(Discret1D *ds)
{
	// This function is called either in the destructor or when it's likely
	// that communicator is no longer necessary. Calling it is safe, because
	// the constructor will be called anyways when necessary.

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(ds->comm != MPI_COMM_NULL)
	{
		ierr = MPI_Comm_free(&ds->comm); CHKERRQ(ierr);

		ds->comm = MPI_COMM_NULL;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DGatherCoord(Discret1D *ds, PetscScalar **coord)
{
	// gather coordinate array on rank zero of PETSC_COMM_WORLD
	// WARNING! the array only exists on rank zero of PETSC_COMM_WORLD
	// WARNING! the array must be destroyed after use!

	PetscInt     i;
	PetscScalar *pcoord;
	PetscMPIInt *recvcnts;
	PetscMPIInt *recvdisp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	pcoord   = NULL;
	recvcnts = NULL;
	recvdisp = NULL;

	// create column communicator
	ierr = Discret1DGetColumnComm(ds); CHKERRQ(ierr);

	// check for sequential case
	if(ds->nproc == 1)
	{
		// copy coordinates on rank zero of PETSC_COMM_WORLD
		if(ISRankZero(PETSC_COMM_WORLD))
		{
			ierr = makeScalArray(&pcoord, ds->ncoor, ds->tnods); CHKERRQ(ierr);
		}
	}
	else
	{
		// gather coordinates on ranks zero of column communicator
		if(ISRankZero(ds->comm))
		{
			// allocate coordinates
			ierr = makeScalArray(&pcoord, NULL, ds->tnods); CHKERRQ(ierr);

			// allocate receive counts
			ierr = makeMPIIntArray(&recvcnts, NULL, ds->nproc) ; CHKERRQ(ierr);

			// allocate receive displacements
			ierr = makeMPIIntArray(&recvdisp, NULL, ds->nproc) ; CHKERRQ(ierr);

			// compute receive counts
			for(i = 0; i < ds->nproc; i++) recvcnts[i] = (PetscMPIInt)(ds->starts[i+1] - ds->starts[i]);

			// ds->starts[ds->nproc] stores index of last node (not total number of nodes)
			recvcnts[ds->nproc-1]++;

			// store receive displacements
			for(i = 0; i < ds->nproc; i++) recvdisp[i] = (PetscMPIInt)ds->starts[i];
		}

		// gather coordinates
		ierr = MPI_Gatherv(ds->ncoor, (PetscMPIInt)ds->nnods, MPIU_SCALAR,
			pcoord, recvcnts, recvdisp, MPIU_SCALAR, 0, ds->comm); CHKERRQ(ierr);

		// free memory
		if(!ISRankZero(PETSC_COMM_WORLD))
		{	ierr = PetscFree(pcoord);   CHKERRQ(ierr); }
			ierr = PetscFree(recvcnts); CHKERRQ(ierr);
			ierr = PetscFree(recvdisp); CHKERRQ(ierr);
	}

	// return coordinates
	(*coord) = pcoord;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DCheckMG(Discret1D *ds, const char *dir, PetscInt *_ncors)
{
	PetscInt sz, ncors;

	PetscFunctionBeginUser;

	// check whether local grid size is an even number
	if(ds->ncels % 2)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Local grid size is an odd number in %s-direction", dir);
	}

	// check uniform local grid size (constant on all processors)
	if(ds->tcels % ds->nproc)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Uniform local grid size doesn't exist in %s-direction", dir);
	}

	// compare actual grid size with uniform value
	if(ds->tcels/ds->nproc != ds->ncels)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Local grid size is not constant on all processors in %s-direction", dir);
	}

	// determine maximum number of coarsening steps (enforce at least two coarse grid cells per processor)
	sz    = ds->ncels;
	ncors = 0;
	while(!(sz % 2) && sz > 2) { sz /= 2; ncors++; }

	// return
	(*_ncors) = ncors;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DgetMaxInvStep(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax)
{
	// get maximum inverse time step on local domain

	PetscScalar v, h, vmax, idt, idtmax;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, idx, ijk[3], jj, ln;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize
	idtmax = (*_idtmax);

	if(!ds->uniform)
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

		// get uniform mesh step
		h = (ds->gcrdend - ds->gcrdbeg)/(PetscScalar)ds->tcels;

		// get inverse time step
		idt = vmax/h;

		// update maximum inverse time step
		if(idt > idtmax) idtmax = idt;
	}

	// return result
	(*_idtmax) = idtmax;

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
PetscErrorCode Discret1DFindPoint(Discret1D *ds, PetscScalar x, PetscInt &ID)
{
	// find index of a cell containing point (local points only)

	PetscScalar  *px, dx, tol;
	PetscInt      n, M, L, R;

	PetscFunctionBeginUser;

	n   =  ds->ncels;
	px  =  ds->ncoor;
	dx  = (px[n] - px[0])/((PetscScalar)n);
	tol =  ds->gtol*dx;

	// check bounds
	if(x < px[0] - tol || x > px[n] + tol)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Non-local point cannot be mapped to local cell");
	}

	if(ds->uniform)
	{
		// get cell index
		ID = (PetscInt)PetscFloorReal((x - px[0])/dx);

		// check bounds
		if(ID < 0)   ID = 0;
		if(ID > n-1) ID = n-1;
	}
	else
	{
		// binary search
		L = 0;
		R = n;

		while((R - L) > 1)
		{
			M = (L + R)/2;
			if(px[M] <= x) L = M;
			if(px[M] >= x) R = M;
		}

		ID = L;

		if(ID < 0 || ID > n-1)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Out-of-bound cell index occurred while mapping point to cell");
		}
	}

	PetscFunctionReturn(0);
/*
	Discret1D       ds;
	PetscInt        ID;
	PetscErrorCode 	ierr;
	PetscScalar     x[]     = { 0.0, 0.3, 0.8, 1.4, 2.4, 2.9, 3.2, 3.5 };
	PetscInt        uniform = 0;
	PetscScalar     x[]     = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 };
	PetscInt        uniform = 1;
	PetscScalar     p = 2.0 - 2.0*DBL_EPSILON;
	ds.ncels   = 7;
	ds.ncoor   = x;
	ds.uniform = uniform;
	ds.gtol    = 1e-9;
	ierr = Discret1DFindPoint(&ds, p, ID); CHKERRQ(ierr);

 */
}
//---------------------------------------------------------------------------
// DOFIndex functions
//---------------------------------------------------------------------------
PetscErrorCode DOFIndexCreate(DOFIndex *dof, DM DA_CEN, DM DA_X, DM DA_Y, DM DA_Z)
{
	// compute number of local dof and starting indices

	PetscInt nx, ny, nz, NUM[2], SUM[3];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get local number of dof
	ierr = DMDAGetCorners(DA_X,   NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr); dof->lnvx = nx*ny*nz;
	ierr = DMDAGetCorners(DA_Y,   NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr); dof->lnvy = nx*ny*nz;
	ierr = DMDAGetCorners(DA_Z,   NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr); dof->lnvz = nx*ny*nz;
	ierr = DMDAGetCorners(DA_CEN, NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr); dof->lnp  = nx*ny*nz;

	dof->lnv = dof->lnvx +  dof->lnvy +  dof->lnvz;

	NUM[0] = dof->lnv;
	NUM[1] = dof->lnp;

	// compute prefix sums
	ierr = MPI_Scan(NUM, SUM, 2, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

	// set starting indices
	dof->stv = SUM[0] - dof->lnv;
	dof->stp = SUM[1] - dof->lnp;

    dof->ln = dof->lnv + dof->lnp;
    dof->st = dof->stv + dof->stp;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// FDSTAG functions
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGCreate(FDSTAG *fs, FB *fb)
{
	// Create object with all necessary arrays to handle FDSTAG discretization.

	// NOTE: velocity components have one layer of boundary ghost points.
	// The idea is that velocity vectors should contain sufficient information
	// to compute strain/rates/stresses/residuals including boundary conditions.

	Scaling         *scal;
	PetscMPIInt      rank;
	const PetscInt  *plx, *ply, *plz;
	PetscInt        *lx,  *ly,  *lz;
	PetscInt         rx,   ry,   rz;
	PetscInt         cx,   cy,   cz;
	PetscInt         Nx,   Ny,   Nz;
	PetscInt         Px,   Py,   Pz;
	MeshSeg1D        msx,  msy,  msz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	scal = fs->scal;

	// set & read geometry tolerance
	fs->gtol = 1e-6;
	ierr = getScalarParam(fb, _OPTIONAL_, "gtol", &fs->gtol, 1, 1.0); CHKERRQ(ierr);

	// set number of processors
	Px = PETSC_DECIDE;
	Py = PETSC_DECIDE;
	Pz = PETSC_DECIDE;

	// fix number of processors in all directions
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_x", &Px, 1, _max_num_procs_); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_y", &Py, 1, _max_num_procs_); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_z", &Pz, 1, _max_num_procs_); CHKERRQ(ierr);

	// read mesh parameters
	ierr = MeshSeg1DReadParam(&msx, scal->length, fs->gtol, "x", fb); CHKERRQ(ierr);
	ierr = MeshSeg1DReadParam(&msy, scal->length, fs->gtol, "y", fb); CHKERRQ(ierr);
	ierr = MeshSeg1DReadParam(&msz, scal->length, fs->gtol, "z", fb); CHKERRQ(ierr);

	// get total number of nodes
	Nx = msx.tcels + 1;
	Ny = msy.tcels + 1;
	Nz = msz.tcels + 1;

	// partition central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz-1, Px, Py, Pz, 1, 1, 0, 0, 0, &fs->DA_CEN); CHKERRQ(ierr);

	// get actual number of processors (can be different compared to given)
	ierr = DMDAGetInfo(fs->DA_CEN, 0, 0, 0, 0, &Px, &Py, &Pz, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	// get number of cells per processor
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &plx, &ply, &plz); CHKERRQ(ierr);

	ierr = makeIntArray(&lx, plx, Px); CHKERRQ(ierr);
	ierr = makeIntArray(&ly, ply, Py); CHKERRQ(ierr);
	ierr = makeIntArray(&lz, plz, Pz); CHKERRQ(ierr);

	// get number of nodes per processor (only different on the last processor)
	lx[Px-1]++; ly[Py-1]++; lz[Pz-1]++;

	// create corner, face and edge DMDA objects
	ierr = FDSTAGCreateDMDA(fs, Nx, Ny, Nz, Px, Py, Pz, lx, ly, lz); CHKERRQ(ierr);

	// setup indexing data
	ierr = DOFIndexCreate(&fs->dof, fs->DA_CEN, fs->DA_X, fs->DA_Y, fs->DA_Z); CHKERRQ(ierr);

	// get MPI processor rank
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// determine i-j-k ranks of processor
	getLocalRank(&rx, &ry, &rz, rank, Px, Py);

	// compute column colors
	cx = ry + rz*Py; // global index in YZ-plane
	cy = rx + rz*Px; // global index in XZ-plane
	cz = rx + ry*Px; // global index in XY-plane

	// set discretization / domain decomposition data
	ierr = Discret1DCreate(&fs->dsx, Px, rx, lx, cx,
			getGlobalRank(rx-1, ry, rz, Px, Py, Pz),
			getGlobalRank(rx+1, ry, rz, Px, Py, Pz),
			fs->gtol, "x"); CHKERRQ(ierr);

	ierr = Discret1DCreate(&fs->dsy, Py, ry, ly, cy,
			getGlobalRank(rx, ry-1, rz, Px, Py, Pz),
			getGlobalRank(rx, ry+1, rz, Px, Py, Pz),
			fs->gtol, "y"); CHKERRQ(ierr);

	ierr = Discret1DCreate(&fs->dsz, Pz, rz, lz, cz,
			getGlobalRank(rx, ry, rz-1, Px, Py, Pz),
			getGlobalRank(rx, ry, rz+1, Px, Py, Pz),
			fs->gtol, "z"); CHKERRQ(ierr);

	// delete temporary arrays
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	// set number of local grid points
	ierr = FDSTAGSetNum(fs); CHKERRQ(ierr);

	// get ranks of neighbor processes
	ierr = FDSTAGGetNeighbProc(fs); CHKERRQ(ierr);

	// generate coordinates
	ierr = Discret1DGenCoord(&fs->dsx, &msx); CHKERRQ(ierr);
	ierr = Discret1DGenCoord(&fs->dsy, &msy); CHKERRQ(ierr);
	ierr = Discret1DGenCoord(&fs->dsz, &msz); CHKERRQ(ierr);

	// print essential grid details
	ierr = FDSTAGView(fs); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGReadRestart(FDSTAG *fs, FILE *fp)
{
	PetscInt *lx,  *ly,  *lz;
	PetscInt  Nx,   Ny,   Nz;
	PetscInt  Px,   Py,   Pz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = Discret1DReadRestart(&fs->dsx, fp); CHKERRQ(ierr);
	ierr = Discret1DReadRestart(&fs->dsy, fp); CHKERRQ(ierr);
	ierr = Discret1DReadRestart(&fs->dsz, fp); CHKERRQ(ierr);

	// get total number of nodes
	Nx = fs->dsx.tnods;
	Ny = fs->dsy.tnods;
	Nz = fs->dsz.tnods;

	// get number of processes
	Px = fs->dsx.nproc;
	Py = fs->dsy.nproc;
	Pz = fs->dsz.nproc;

	// get number cells per processor
	lx = NULL;
	ly = NULL;
	lz = NULL;
	
	ierr = Discret1DGetNumCells(&fs->dsx, &lx); CHKERRQ(ierr);
	ierr = Discret1DGetNumCells(&fs->dsy, &ly); CHKERRQ(ierr);
	ierr = Discret1DGetNumCells(&fs->dsz, &lz); CHKERRQ(ierr);

	// central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_CEN); CHKERRQ(ierr);

	// get number of nodes per processor (only different on the last processor)
	lx[Px-1]++; ly[Py-1]++; lz[Pz-1]++;

	// create corner, face and edge DMDA objects
	ierr = FDSTAGCreateDMDA(fs, Nx, Ny, Nz, Px, Py, Pz, lx, ly, lz); CHKERRQ(ierr);

	// setup indexing data
	ierr = DOFIndexCreate(&fs->dof, fs->DA_CEN, fs->DA_X, fs->DA_Y, fs->DA_Z); CHKERRQ(ierr);

	// delete temporary arrays
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	fs->dsx.comm = MPI_COMM_NULL;
	fs->dsy.comm = MPI_COMM_NULL;
	fs->dsz.comm = MPI_COMM_NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGWriteRestart(FDSTAG *fs, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = Discret1DWriteRestart(&fs->dsx, fp); CHKERRQ(ierr);
	ierr = Discret1DWriteRestart(&fs->dsy, fp); CHKERRQ(ierr);
	ierr = Discret1DWriteRestart(&fs->dsz, fp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGCoarsen(FDSTAG *coarse, FDSTAG *fine)
{
	PetscInt         i;
	PetscInt         Nx,    Ny,    Nz;
	PetscInt         Px,    Py,    Pz;
	const PetscInt  *plx,  *ply,  *plz;
	PetscInt        *lx,   *ly,   *lz;
	Discret1D       *fdsx, *fdsy, *fdsz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// clear memory
	ierr = PetscMemzero(coarse, sizeof(FDSTAG)); CHKERRQ(ierr);

	// copy data
	coarse->scal = fine->scal;
	coarse->gtol = fine->gtol;
	for(i = 0; i < _num_neighb_; i++) { coarse->neighb[i] = fine->neighb[i]; }

	// get number of cells & processors in the fine grid
	ierr = DMDAGetInfo(fine->DA_CEN, 0, &Nx, &Ny, &Nz, &Px, &Py, &Pz, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	// get number of cells per processor in fine grid
	ierr = DMDAGetOwnershipRanges(fine->DA_CEN, &plx, &ply, &plz); CHKERRQ(ierr);

	ierr = makeIntArray(&lx, plx, Px); CHKERRQ(ierr);
	ierr = makeIntArray(&ly, ply, Py); CHKERRQ(ierr);
	ierr = makeIntArray(&lz, plz, Pz); CHKERRQ(ierr);

	if(Nx > 2) { Nx /= 2;  for(i = 0; i < Px; i++) { lx[i] /= 2; } }
	if(Ny > 2) { Ny /= 2;  for(i = 0; i < Py; i++) { ly[i] /= 2; } }
	if(Nz > 2) { Nz /= 2;  for(i = 0; i < Pz; i++) { lz[i] /= 2; } }

	// central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &coarse->DA_CEN); CHKERRQ(ierr);

	// get total number of nodes
	Nx++; Ny++; Nz++;

	// get number of nodes per processor (only different on the last processor)
	lx[Px-1]++; ly[Py-1]++; lz[Pz-1]++;

	// create corner, face and edge DMDA objects
	ierr = FDSTAGCreateDMDA(coarse, Nx, Ny, Nz, Px, Py, Pz, lx, ly, lz); CHKERRQ(ierr);

	// create index arrays
	ierr = DOFIndexCreate(&coarse->dof, coarse->DA_CEN, coarse->DA_X, coarse->DA_Y, coarse->DA_Z); CHKERRQ(ierr);

	// set discretization / domain decomposition data
	fdsx = &fine->dsx;
	fdsy = &fine->dsy;
	fdsz = &fine->dsz;

	ierr = Discret1DCreate(&coarse->dsx, Px, fdsx->rank, lx, fdsx->color,
			fdsx->grprev, fdsx->grnext, coarse->gtol, "x"); CHKERRQ(ierr);

	ierr = Discret1DCreate(&coarse->dsy, Py, fdsy->rank, ly, fdsy->color,
			fdsy->grprev, fdsy->grnext, coarse->gtol, "y"); CHKERRQ(ierr);

	ierr = Discret1DCreate(&coarse->dsz, Pz, fdsz->rank, lz, fdsz->color,
			fdsz->grprev, fdsz->grnext, coarse->gtol, "z"); CHKERRQ(ierr);

	// clear temporary storage
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

	// set number of local grid points
	ierr = FDSTAGSetNum(coarse); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGCoarsenCoord(FDSTAG *coarse, FDSTAG *fine)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// coarsen coordinates
	ierr = Discret1DCoarsenCoord(&coarse->dsx, &fine->dsx); CHKERRQ(ierr);
	ierr = Discret1DCoarsenCoord(&coarse->dsy, &fine->dsy); CHKERRQ(ierr);
	ierr = Discret1DCoarsenCoord(&coarse->dsz, &fine->dsz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGDestroy(FDSTAG * fs)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// destroy distributed arrays
	ierr = DMDestroy(&fs->DA_CEN);     CHKERRQ(ierr);
	ierr = DMDestroy(&fs->DA_COR);     CHKERRQ(ierr);

	ierr = DMDestroy(&fs->DA_XY);      CHKERRQ(ierr);
	ierr = DMDestroy(&fs->DA_XZ);      CHKERRQ(ierr);
	ierr = DMDestroy(&fs->DA_YZ);      CHKERRQ(ierr);

	ierr = DMDestroy(&fs->DA_X);       CHKERRQ(ierr);
	ierr = DMDestroy(&fs->DA_Y);       CHKERRQ(ierr);
	ierr = DMDestroy(&fs->DA_Z);       CHKERRQ(ierr);

	// destroy discretization data
	ierr = Discret1DDestroy(&fs->dsx); CHKERRQ(ierr);
	ierr = Discret1DDestroy(&fs->dsy); CHKERRQ(ierr);
	ierr = Discret1DDestroy(&fs->dsz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGCreateDMDA(FDSTAG *fs,
	PetscInt  Nx, PetscInt  Ny, PetscInt  Nz,
	PetscInt  Px, PetscInt  Py, PetscInt  Pz,
	PetscInt *lx, PetscInt *ly, PetscInt *lz)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// corners (DA_COR) no boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_COR); CHKERRQ(ierr);

	// XY edges (DA_XY) no boundary ghost points (1-layer stencil box)
	lz[Pz-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_XY); CHKERRQ(ierr);
	lz[Pz-1]++;

	// XZ edges (DA_XZ) no boundary ghost points (1-layer stencil box)
	ly[Py-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_XZ); CHKERRQ(ierr);
	ly[Py-1]++;

	// YZ edges (DA_YZ) no boundary ghost points (1-layer stencil box)
	lx[Px-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_YZ); CHKERRQ(ierr);
	lx[Px-1]++;

	// X face (DA_X) with boundary ghost points (1-layer stencil box)
	ly[Py-1]--; lz[Pz-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_X); CHKERRQ(ierr);
	ly[Py-1]++; lz[Pz-1]++;

	// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; lz[Pz-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_Y); CHKERRQ(ierr);
	lx[Px-1]++; lz[Pz-1]++;

	// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; ly[Py-1]--;
	ierr = DMDACreate3DSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_Z); CHKERRQ(ierr);
	lx[Px-1]++; ly[Py-1]++;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGSetNum(FDSTAG *fs)
{
	// set number of local grid points

	PetscInt nnx, nny, nnz;
	PetscInt ncx, ncy, ncz;

	PetscFunctionBeginUser;

	// compute local number of grid points
	nnx = fs->dsx.nnods; ncx = fs->dsx.ncels;
	nny = fs->dsy.nnods; ncy = fs->dsy.ncels;
	nnz = fs->dsz.nnods; ncz = fs->dsz.ncels;

	fs->nCells = ncx*ncy*ncz;
	fs->nCorns = nnx*nny*nnz;
	fs->nXYEdg = nnx*nny*ncz;
	fs->nXZEdg = nnx*ncy*nnz;
	fs->nYZEdg = ncx*nny*nnz;
	fs->nXFace = nnx*ncy*ncz;
	fs->nYFace = ncx*nny*ncz;
	fs->nZFace = ncx*ncy*nnz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGGetNeighbProc(FDSTAG *fs)
{
	// return an array with the global ranks of adjacent processes (including itself)

	PetscInt i, j, k, rx, ry, rz, Px, Py, Pz, ptx, pty, ptz, cnt;
	PetscFunctionBeginUser;

	// get ranks
	rx = fs->dsx.rank;
	ry = fs->dsy.rank;
	rz = fs->dsz.rank;

	// get number processors
	Px = fs->dsx.nproc;
	Py = fs->dsy.nproc;
	Pz = fs->dsz.nproc;

	// get periodic topology flags
	ptx = fs->dsx.periodic;
	pty = fs->dsy.periodic;
	ptz = fs->dsz.periodic;

	// clear counter
	cnt = 0;

	for(k = -1; k < 2; k++)
	{	for(j = -1; j < 2; j++)
		{	for(i = -1; i < 2; i++)
			{
				fs->neighb[cnt++] = getGlobalRankPeriodic(rx+i, ry+j, rz+k, Px, Py, Pz, ptx, pty, ptz);
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGGetPointRanks(FDSTAG *fs, PetscScalar *X, PetscInt *lrank, PetscMPIInt *grank)
{
	// get local & global ranks of a domain containing a point (only neighbors are checked)

	PetscInt    rx, ry, rz;
	PetscScalar bx, by, bz;
	PetscScalar ex, ey, ez;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get local coordinate bounds
	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	// gel local relative ranks
	GET_POINT_REL_RANK(rx, X[0], bx, ex);
	GET_POINT_REL_RANK(ry, X[1], by, ey);
	GET_POINT_REL_RANK(rz, X[2], bz, ez);


	(*lrank) = rx + 3*ry + 9*rz;
	(*grank) = fs->neighb[(*lrank)];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGGetAspectRatio(FDSTAG *fs, PetscScalar *maxAspRat)
{
	// compute maximum aspect ratio in the grid

	PetscScalar dx, dy, dz, rt, lrt, grt;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	lrt = 0.0;

	START_STD_LOOP
	{
		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		if(dx > dy) rt = dx/dy; else rt = dy/dx; if(rt > lrt) lrt = rt;
		if(dx > dz) rt = dx/dz; else rt = dz/dx; if(rt > lrt) lrt = rt;
		if(dy > dz) rt = dy/dz; else rt = dz/dy; if(rt > lrt) lrt = rt;
	}
	END_STD_LOOP

	// get global aspect ratio
	if(ISParallel(PETSC_COMM_WORLD))
	{
		// exchange
		ierr = MPI_Allreduce(&lrt, &grt, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

	}
	else
	{
		// there is no difference between global & local values
		grt = lrt;
	}

	// store the result
	(*maxAspRat) = grt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGView(FDSTAG *fs)
{
	// print & check essential grid details

	PetscMPIInt nproc;
	PetscScalar bx, by, bz;
	PetscScalar ex, ey, ez;
	PetscScalar maxAspRat, chLen;
	PetscInt    px, py, pz, cx, cy, cz, nx, ny, nz, nVelDOF, nCells;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	chLen = fs->scal->length;

	px = fs->dsx.nproc;  cx = fs->dsx.tcels;  nx = fs->dsx.tnods;
	py = fs->dsy.nproc;  cy = fs->dsy.tcels;  ny = fs->dsy.tnods;
	pz = fs->dsz.nproc;  cz = fs->dsz.tcels;  nz = fs->dsz.tnods;

	nCells  = cx*cy*cz;
	nVelDOF = nx*cy*cz + cx*ny*cz + cx*cy*nz;

	ierr = FDSTAGGetAspectRatio(fs, &maxAspRat); CHKERRQ(ierr);

	ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Grid parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Total number of cpu                  : %lld \n", (LLD)nproc);
	PetscPrintf(PETSC_COMM_WORLD, "   Processor grid  [nx, ny, nz]         : [%lld, %lld, %lld]\n", (LLD)px, (LLD)py, (LLD)pz);
	PetscPrintf(PETSC_COMM_WORLD, "   Fine grid cells [nx, ny, nz]         : [%lld, %lld, %lld]\n", (LLD)cx, (LLD)cy, (LLD)cz);
	PetscPrintf(PETSC_COMM_WORLD, "   Number of cells                      :  %lld\n", (LLD)nCells);
	PetscPrintf(PETSC_COMM_WORLD, "   Number of faces                      :  %lld\n", (LLD)nVelDOF);
	PetscPrintf(PETSC_COMM_WORLD, "   Maximum cell aspect ratio            :  %7.5f\n", maxAspRat);
	PetscPrintf(PETSC_COMM_WORLD, "   Lower coordinate bounds [bx, by, bz] : [%g, %g, %g]\n", bx*chLen, by*chLen, bz*chLen);
	PetscPrintf(PETSC_COMM_WORLD, "   Upper coordinate bounds [ex, ey, ez] : [%g, %g, %g]\n", ex*chLen, ey*chLen, ez*chLen);


	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	if(maxAspRat > 10.0) PetscPrintf(PETSC_COMM_WORLD, " Don't expect any magic with this aspect ratio %g ...\n", maxAspRat);
	if(maxAspRat > 30.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, " Everything has a limit, reduce this aspect ratio: %g ...\n", maxAspRat);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGGetLocalBox(
	FDSTAG      *fs,
	PetscScalar *bx,
	PetscScalar *by,
	PetscScalar *bz,
	PetscScalar *ex,
	PetscScalar *ey,
	PetscScalar *ez)
{
	PetscFunctionBeginUser;

	if(bx) (*bx) = fs->dsx.ncoor[0];
	if(by) (*by) = fs->dsy.ncoor[0];
	if(bz) (*bz) = fs->dsz.ncoor[0];

	if(ex) (*ex) = fs->dsx.ncoor[fs->dsx.ncels];
	if(ey) (*ey) = fs->dsy.ncoor[fs->dsy.ncels];
	if(ez) (*ez) = fs->dsz.ncoor[fs->dsz.ncels];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGGetGlobalBox(
	FDSTAG      *fs,
	PetscScalar *bx,
	PetscScalar *by,
	PetscScalar *bz,
	PetscScalar *ex,
	PetscScalar *ey,
	PetscScalar *ez)
{
	PetscFunctionBeginUser;

	if(bx) (*bx) = fs->dsx.gcrdbeg;
	if(by) (*by) = fs->dsy.gcrdbeg;
	if(bz) (*bz) = fs->dsz.gcrdbeg;

	if(ex) (*ex) = fs->dsx.gcrdend;
	if(ey) (*ey) = fs->dsy.gcrdend;
	if(ez) (*ez) = fs->dsz.gcrdend;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FDSTAGSaveGrid(FDSTAG *fs)
{
	int            fid;
	char           *fname;
	PetscScalar    *xc, *yc, *zc, chLen;
	PetscMPIInt    rank;
	PetscLogDouble t;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PrintStart(&t, "Saving processor partitioning", NULL);

	// characteristic length
	chLen = fs->scal->length;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// gather global coord
	ierr = Discret1DGatherCoord(&fs->dsx, &xc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsy, &yc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsz, &zc); CHKERRQ(ierr);

	if(rank == 0)
	{
		// save file
		asprintf(&fname, "ProcessorPartitioning_%lldcpu_%lld.%lld.%lld.bin",
			(LLD)(fs->dsx.nproc*fs->dsy.nproc*fs->dsz.nproc),
			(LLD)fs->dsx.nproc, (LLD)fs->dsy.nproc, (LLD)fs->dsz.nproc);

		PetscBinaryOpen(fname, FILE_MODE_WRITE, &fid);

		PetscBinaryWrite(fid, &fs->dsx.nproc, 1,               PETSC_INT);
		PetscBinaryWrite(fid, &fs->dsy.nproc, 1,               PETSC_INT);
		PetscBinaryWrite(fid, &fs->dsz.nproc, 1,               PETSC_INT);
		PetscBinaryWrite(fid, &fs->dsx.tnods, 1,               PETSC_INT);
		PetscBinaryWrite(fid, &fs->dsy.tnods, 1,               PETSC_INT);
		PetscBinaryWrite(fid, &fs->dsz.tnods, 1,               PETSC_INT);
		PetscBinaryWrite(fid, fs->dsx.starts, fs->dsx.nproc+1, PETSC_INT);
		PetscBinaryWrite(fid, fs->dsy.starts, fs->dsy.nproc+1, PETSC_INT);
		PetscBinaryWrite(fid, fs->dsz.starts, fs->dsz.nproc+1, PETSC_INT);
		PetscBinaryWrite(fid, &chLen,         1,               PETSC_SCALAR);
		PetscBinaryWrite(fid, xc,             fs->dsx.tnods,   PETSC_SCALAR);
		PetscBinaryWrite(fid, yc,             fs->dsy.tnods,   PETSC_SCALAR);
		PetscBinaryWrite(fid, zc,             fs->dsz.tnods,   PETSC_SCALAR);

		PetscBinaryClose(fid);
		free(fname);

		ierr = PetscFree(xc); CHKERRQ(ierr);
		ierr = PetscFree(yc); CHKERRQ(ierr);
		ierr = PetscFree(zc); CHKERRQ(ierr);
	}

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DMDACreate3DSetUp(MPI_Comm comm,
	DMBoundaryType bx, DMBoundaryType by, DMBoundaryType bz, DMDAStencilType stencil_type,
	PetscInt M, PetscInt N, PetscInt P, PetscInt m, PetscInt n, PetscInt p,
	PetscInt dof, PetscInt s, const PetscInt lx[], const PetscInt ly[], const PetscInt lz[], DM *da)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = DMDACreate3d(comm, bx, by, bz, stencil_type, M, N, P, m, n, p, dof, s, lx, ly, lz, da); CHKERRQ(ierr);

	ierr = DMSetFromOptions((*da)); CHKERRQ(ierr);
	ierr = DMSetUp((*da));          CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

