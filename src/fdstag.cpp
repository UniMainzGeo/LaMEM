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
 **    filename:   fdstag.c
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
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DReadParam"
PetscErrorCode MeshSeg1DReadParam(
	MeshSeg1D  *ms,
	PetscScalar leng,
	PetscScalar gtol,
	const char *dir,
	FB         *fb)
{
	PetscInt    i, tcels, uniform;
	PetscInt    ncells[MaxNumSegs];
	PetscScalar avgsz, sz;
	char        *nseg, *nel, *coord, *bias;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	ierr = PetscMemzero(ms, sizeof(MeshSeg1D)); CHKERRQ(ierr);

	ms->nsegs = 1;

	for(i = 0; i < MaxNumSegs; i++)
	{
		ms->biases[i] = 1.0;
		ncells    [i] = 0.0;
	}

	// compose option keys
	asprintf(&nseg,  "nseg_%s",  dir);
	asprintf(&nel,   "nel_%s",   dir);
	asprintf(&coord, "coord_%s", dir);
	asprintf(&bias,  "bias_%s",  dir);

	// read parameters
	ierr = getIntParam   (fb, _OPTIONAL_, nseg,  &ms->nsegs,  1,           MaxNumSegs);  CHKERRQ(ierr);
	ierr = getIntParam   (fb, _REQUIRED_, nel,    ncells,     ms->nsegs,   MaxNumCells); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, coord,  ms->xstart, ms->nsegs+1, leng);        CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, bias,   ms->biases, ms->nsegs,   1.0 );        CHKERRQ(ierr);

	// compute starting node indices
	for(i = 0, tcels = 0; i < ms->nsegs; i++)
	{
		ms->istart[i] = tcels;
		tcels        += ncells[i];
	}
	ms->istart[ms->nsegs] = tcels;

	// check ordering
	for(i = 0; i < ms->nsegs; i++)
	{
		if(ms->xstart[i] >= ms->xstart[i+1])
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unordered coordinates in parameter %s\n", coord);
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DGenCoord"
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

	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "Discret1DCreate"
PetscErrorCode Discret1DCreate(
		Discret1D  *ds,
		PetscInt    nproc,     // number of processors
		PetscInt    rank,      // processor rank
		PetscInt   *nnodProc,  // number of nodes per processor
		PetscInt    color,     // column color
		PetscMPIInt grprev,    // global rank of previous process
		PetscMPIInt grnext)    // global rank of next process
{
	PetscInt       i, cnt;
	PetscErrorCode ierr;
	PetscFunctionBegin;

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

	// coordinates of local nodes + 1 layer (left) & 2 layers (right) of ghost points
	// NOTE: on the last processor there is only one ghost point from the right

	if(grnext != -1) ds->bufsz = ds->nnods+3;
	else             ds->bufsz = ds->nnods+2;
	ierr = makeScalArray(&ds->nbuff, 0, ds->bufsz); CHKERRQ(ierr);
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DDestroy"
PetscErrorCode Discret1DDestroy(Discret1D *ds)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;

	// free memory buffers
	ierr = PetscFree(ds->nbuff);        CHKERRQ(ierr);
	ierr = PetscFree(ds->cbuff);        CHKERRQ(ierr);
	ierr = PetscFree(ds->starts);       CHKERRQ(ierr);
	ierr = Discret1DFreeColumnComm(ds); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DReadRestart"
PetscErrorCode Discret1DReadRestart(Discret1D *ds, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = makeIntArray (&ds->starts, NULL, ds->nproc + 1); CHKERRQ(ierr);
	ierr = makeScalArray(&ds->nbuff,  NULL, ds->bufsz    ); CHKERRQ(ierr);
	ierr = makeScalArray(&ds->cbuff,  NULL, ds->ncels + 2); CHKERRQ(ierr);

   	fread(ds->starts, sizeof(PetscInt   )*(size_t)(ds->nproc + 1), 1, fp);
	fread(ds->nbuff,  sizeof(PetscScalar)*(size_t)(ds->bufsz    ), 1, fp);
	fread(ds->cbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 2), 1, fp);

	ds->ncoor = ds->nbuff + 1;
	ds->ccoor = ds->cbuff + 1;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DWriteRestart"
PetscErrorCode Discret1DWriteRestart(Discret1D *ds, FILE *fp)
{
	PetscFunctionBegin;

	fwrite(ds->starts, sizeof(PetscInt   )*(size_t)(ds->nproc + 1), 1, fp);
	fwrite(ds->nbuff,  sizeof(PetscScalar)*(size_t)(ds->bufsz    ), 1, fp);
	fwrite(ds->cbuff,  sizeof(PetscScalar)*(size_t)(ds->ncels + 2), 1, fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGetNumCells"
PetscErrorCode Discret1DGetNumCells(Discret1D *ds, PetscInt **ncelProc)
{
	// get number of cells per processor

	PetscInt i, *l;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = makeIntArray(&l, NULL, ds->nproc); CHKERRQ(ierr);

	for(i = 0; i < ds->nproc; i++)
	{
		l[i] = ds->starts[i+1] - ds->starts[i];
	}

	(*ncelProc) = l;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGenCoord"
PetscErrorCode Discret1DGenCoord(Discret1D *ds, MeshSeg1D *ms)
{
	PetscInt     i, n, nl, pstart, istart;
	PetscScalar *crd, A, B, C;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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

	// set uniform grid flag
	ds->uniform = ms->uniform;

	// set global grid coordinate bounds
	ds->gcrdbeg = ms->xstart[0];
	ds->gcrdend = ms->xstart[ms->nsegs];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DStretch"
PetscErrorCode Discret1DStretch(Discret1D *ds, PetscScalar eps)
{
	// stretch grid with constant stretch factor about coordinate origin.
	// x_new = x_old + eps*x_old

	PetscInt i;

	PetscFunctionBegin;

	// recompute (stretch) node coordinates in the buffer
	for(i = 0; i < ds->bufsz; i++) ds->nbuff[i] *= (1.0 + eps);

	// recompute cell coordinates
	for(i = -1; i < ds->ncels+1; i++)
		ds->ccoor[i] = (ds->ncoor[i] + ds->ncoor[i+1])/2.0;

	// recompute global coordinate bounds
	ds->gcrdbeg *= (1.0 + eps);
	ds->gcrdend *= (1.0 + eps);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGetColumnComm"
PetscErrorCode Discret1DGetColumnComm(Discret1D *ds)
{
	// This function is called every time the column communicator is needed.
	// Nothing is done if communicator already exists or in sequential case.

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(ds->nproc != 1 && ds->comm == MPI_COMM_NULL)
	{
		ierr = MPI_Comm_split(PETSC_COMM_WORLD, ds->color, ds->rank, &ds->comm); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DFreeColumnComm"
PetscErrorCode Discret1DFreeColumnComm(Discret1D *ds)
{
	// This function is called either in the destructor or when it's likely
	// that communicator is no longer necessary. Calling it is safe, because
	// the constructor will be called anyways when necessary.

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(ds->comm != MPI_COMM_NULL)
	{
		ierr = MPI_Comm_free(&ds->comm); CHKERRQ(ierr);

		ds->comm = MPI_COMM_NULL;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGatherCoord"
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
	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "Discret1DCheckMG"
PetscErrorCode Discret1DCheckMG(Discret1D *ds, const char *dir, PetscInt *_ncors)
{
	PetscInt sz, ncors;

	PetscFunctionBegin;

	// check whether local grid size is an even number
	if(ds->ncels % 2)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Local grid size is an odd number in %s-direction", dir);
	}

	// check uniform local grid size (constant on all processors)
	if(ds->tcels % ds->nproc)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Uniform local grid size doesn't exist in %s-direction", dir);
	}

	// compare actual grid size with uniform value
	if(ds->tcels/ds->nproc != ds->ncels)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Local grid size is not constant on all processors in %s-direction", dir);
	}

	// determine maximum number of coarsening steps
	sz    = ds->ncels;
	ncors = 0;
	while(!(sz % 2)) { sz /= 2; ncors++; }

	// return
	(*_ncors) = ncors;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DgetMaxInvStep"
PetscErrorCode Discret1DgetMaxInvStep(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax)
{
	// get maximum inverse time step on local domain

	PetscScalar v, h, vmax, idt, idtmax;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, idx, ijk[3], jj, ln;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
// DOFIndex functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DOFIndexCreate"
PetscErrorCode DOFIndexCreate(DOFIndex *dof, DM DA_CEN, DM DA_X, DM DA_Y, DM DA_Z)
{
	// compute & set global indices of local & ghost nodes

	// **********************************************************************
	// NOTE:
	// for the ghost points, store negative global index of the primary DOF
	// instead of -1
	// **********************************************************************

	PetscInt nx, ny, nz, NUM[2], SUM[3];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear index mode
	dof->idxmod = IDXNONE;

	// store distributed arrays
	dof->DA_X   = DA_X;
	dof->DA_Y   = DA_Y;
	dof->DA_Z   = DA_Z;
	dof->DA_CEN = DA_CEN;

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

	// create index vectors (ghosted)
	ierr = DMCreateLocalVector(DA_X,   &dof->ivx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(DA_Y,   &dof->ivy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(DA_Z,   &dof->ivz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(DA_CEN, &dof->ip);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DOFIndexDestroy"
PetscErrorCode DOFIndexDestroy(DOFIndex *dof)
{
	PetscErrorCode 	 ierr;
	PetscFunctionBegin;

	// destroy index vectors (ghosted)
	ierr = VecDestroy(&dof->ivx); CHKERRQ(ierr);
	ierr = VecDestroy(&dof->ivy); CHKERRQ(ierr);
	ierr = VecDestroy(&dof->ivz); CHKERRQ(ierr);
	ierr = VecDestroy(&dof->ip);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DOFIndexCompute"
PetscErrorCode DOFIndexCompute(DOFIndex *dof, idxtype idxmod)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, stv, stp;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set global indices of the local and ghost nodes (including boundary)
	ierr = VecSet(dof->ivx, -1.0); CHKERRQ(ierr);
	ierr = VecSet(dof->ivy, -1.0); CHKERRQ(ierr);
	ierr = VecSet(dof->ivz, -1.0); CHKERRQ(ierr);
	ierr = VecSet(dof->ip,  -1.0); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(dof->DA_X,   dof->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dof->DA_Y,   dof->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dof->DA_Z,   dof->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dof->DA_CEN, dof->ip,  &ip);   CHKERRQ(ierr);

	//=======================================================
	// compute interlaced global numbering of the local nodes
	//=======================================================

	if     (idxmod == IDXCOUPLED)   { stv = dof->st;  stp = dof->st + dof->lnv; }
	else if(idxmod == IDXUNCOUPLED) { stv = dof->stv; stp = dof->stp;           }

	//---------
	// X-points
	//---------

	ierr = DMDAGetCorners(dof->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivx[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------

	ierr = DMDAGetCorners(dof->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivy[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	ierr = DMDAGetCorners(dof->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivz[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	ierr = DMDAGetCorners(dof->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ip[k][j][i] = (PetscScalar)stp; stp++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(dof->DA_X,   dof->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_Y,   dof->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_Z,   dof->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_CEN, dof->ip,  &ip);   CHKERRQ(ierr);

	// get ghost point indices
	LOCAL_TO_LOCAL(dof->DA_X,   dof->ivx)
	LOCAL_TO_LOCAL(dof->DA_Y,   dof->ivy)
	LOCAL_TO_LOCAL(dof->DA_Z,   dof->ivz)
	LOCAL_TO_LOCAL(dof->DA_CEN, dof->ip)

	// store index mode
	dof->idxmod = idxmod;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// FDSTAG functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCreate"
PetscErrorCode FDSTAGCreate(FDSTAG *fs, FB *fb)
{
	// Create object with all necessary arrays to handle FDSTAG discretization.

	// NOTE: velocity components have one layer of boundary ghost points.
	// The idea is that velocity vectors should contain sufficient information
	// to compute strain/rates/stresses/residuals including boundary conditions.

	Scaling          *scal;
	PetscMPIInt      rank;
	PetscInt         nnx, nny, nnz;
	PetscInt         ncx, ncy, ncz;
	const PetscInt  *plx, *ply, *plz;
	PetscInt        *lx,  *ly,  *lz;
	PetscInt         rx,   ry,   rz;
	PetscInt         cx,   cy,   cz;
	PetscInt         Nx,   Ny,   Nz;
	PetscInt         Px,   Py,   Pz;
	MeshSeg1D        msx,  msy,  msz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = fs->scal;

	// set & read geometry tolerance
	fs->gtol = 1e-9;
	ierr = getScalarParam(fb, _OPTIONAL_, "gtol", &fs->gtol, 1, 1.0); CHKERRQ(ierr);

	// set number of processors
	Px = PETSC_DECIDE;
	Py = PETSC_DECIDE;
	Pz = PETSC_DECIDE;

	// fix number of processors in all directions
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_x", &Px, 1, MaxNumProcs); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_y", &Py, 1, MaxNumProcs); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "cpu_z", &Pz, 1, MaxNumProcs); CHKERRQ(ierr);

	// read mesh parameters
	ierr = MeshSeg1DReadParam(&msx, scal->length, fs->gtol, "x", fb); CHKERRQ(ierr);
	ierr = MeshSeg1DReadParam(&msy, scal->length, fs->gtol, "y", fb); CHKERRQ(ierr);
	ierr = MeshSeg1DReadParam(&msz, scal->length, fs->gtol, "z", fb); CHKERRQ(ierr);

	// get total number of nodes
	Nx = msx.tcels + 1;
	Ny = msy.tcels + 1;
	Nz = msz.tcels + 1;

	// partition central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
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
			getGlobalRank(rx+1, ry, rz, Px, Py, Pz)); CHKERRQ(ierr);

	ierr = Discret1DCreate(&fs->dsy, Py, ry, ly, cy,
			getGlobalRank(rx, ry-1, rz, Px, Py, Pz),
			getGlobalRank(rx, ry+1, rz, Px, Py, Pz)); CHKERRQ(ierr);

	ierr = Discret1DCreate(&fs->dsz, Pz, rz, lz, cz,
			getGlobalRank(rx, ry, rz-1, Px, Py, Pz),
			getGlobalRank(rx, ry, rz+1, Px, Py, Pz)); CHKERRQ(ierr);

	// delete temporary arrays
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);

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
#undef __FUNCT__
#define __FUNCT__ "FDSTAGReadRestart"
PetscErrorCode FDSTAGReadRestart(FDSTAG *fs, FILE *fp)
{
	PetscInt *lx,  *ly,  *lz;
	PetscInt  Nx,   Ny,   Nz;
	PetscInt  Px,   Py,   Pz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
	ierr = Discret1DGetNumCells(&fs->dsx, &lx); CHKERRQ(ierr);
	ierr = Discret1DGetNumCells(&fs->dsy, &ly); CHKERRQ(ierr);
	ierr = Discret1DGetNumCells(&fs->dsz, &lz); CHKERRQ(ierr);

	// central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
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
#undef __FUNCT__
#define __FUNCT__ "FDSTAGWriteRestart"
PetscErrorCode FDSTAGWriteRestart(FDSTAG *fs, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = Discret1DWriteRestart(&fs->dsx, fp); CHKERRQ(ierr);
	ierr = Discret1DWriteRestart(&fs->dsy, fp); CHKERRQ(ierr);
	ierr = Discret1DWriteRestart(&fs->dsz, fp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGDestroy"
PetscErrorCode FDSTAGDestroy(FDSTAG * fs)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

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

	// destroy indexing data
	ierr = DOFIndexDestroy(&fs->dof);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCreateDMDA"
PetscErrorCode FDSTAGCreateDMDA(FDSTAG *fs,
	PetscInt  Nx, PetscInt  Ny, PetscInt  Nz,
	PetscInt  Px, PetscInt  Py, PetscInt  Pz,
	PetscInt *lx, PetscInt *ly, PetscInt *lz)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// corners (DA_COR) no boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_COR); CHKERRQ(ierr);

	// XY edges (DA_XY) no boundary ghost points (1-layer stencil box)
	lz[Pz-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_XY); CHKERRQ(ierr);
	lz[Pz-1]++;

	// XZ edges (DA_XZ) no boundary ghost points (1-layer stencil box)
	ly[Py-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_XZ); CHKERRQ(ierr);
	ly[Py-1]++;

	// YZ edges (DA_YZ) no boundary ghost points (1-layer stencil box)
	lx[Px-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_YZ); CHKERRQ(ierr);
	lx[Px-1]++;

	// X face (DA_X) with boundary ghost points (1-layer stencil box)
	ly[Py-1]--; lz[Pz-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_X); CHKERRQ(ierr);
	ly[Py-1]++; lz[Pz-1]++;

	// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; lz[Pz-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz-1, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_Y); CHKERRQ(ierr);
	lx[Px-1]++; lz[Pz-1]++;

	// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; ly[Py-1]--;
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz, Px, Py, Pz, 1, 1, lx, ly, lz, &fs->DA_Z); CHKERRQ(ierr);
	lx[Px-1]++; ly[Py-1]++;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGetNeighbProc"
PetscErrorCode FDSTAGGetNeighbProc(FDSTAG *fs)
{
	// return an array with the global ranks of adjacent processes (including itself)

	PetscInt i, j, k, rx, ry, rz, Px, Py, Pz, cnt;
	PetscFunctionBegin;

	// get ranks
	rx = fs->dsx.rank;
	ry = fs->dsy.rank;
	rz = fs->dsz.rank;

	// get number processors
	Px = fs->dsx.nproc;
	Py = fs->dsy.nproc;
	Pz = fs->dsz.nproc;

	// clear counter
	cnt = 0;

	for(k = -1; k < 2; k++)
	{	for(j = -1; j < 2; j++)
		{	for(i = -1; i < 2; i++)
			{
				fs->neighb[cnt++] = getGlobalRank(rx+i, ry+j, rz+k, Px, Py, Pz);
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGetPointRanks"
PetscErrorCode FDSTAGGetPointRanks(FDSTAG *fs, PetscScalar *X, PetscInt *lrank, PetscMPIInt *grank)
{
	// get local & global ranks of a domain containing a point (only neighbors are checked)

	PetscInt    rx, ry, rz;
	PetscScalar bx, by, bz;
	PetscScalar ex, ey, ez;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGetAspectRatio"
PetscErrorCode FDSTAGGetAspectRatio(FDSTAG *fs, PetscScalar *maxAspRat)
{
	// compute maximum aspect ratio in the grid

	PetscMPIInt nproc;
	PetscScalar dx, dy, dz, rt, lrt, grt;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get number of processors
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);

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
	if(nproc != 1)
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
#undef __FUNCT__
#define __FUNCT__ "FDSTAGView"
PetscErrorCode FDSTAGView(FDSTAG *fs)
{
	// print & check essential grid details

	PetscMPIInt nproc;
	PetscScalar bx, by, bz;
	PetscScalar ex, ey, ez;
	PetscScalar maxAspRat, chLen;
	PetscInt    px, py, pz, cx, cy, cz, nx, ny, nz, nVelDOF, nCells;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
	if(maxAspRat > 30.0) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, " Everything has a limit, reduce this aspect ratio: %g ...\n", maxAspRat);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGetLocalBox"
PetscErrorCode FDSTAGGetLocalBox(
	FDSTAG      *fs,
	PetscScalar *bx,
	PetscScalar *by,
	PetscScalar *bz,
	PetscScalar *ex,
	PetscScalar *ey,
	PetscScalar *ez)
{
	PetscFunctionBegin;

	if(bx) (*bx) = fs->dsx.ncoor[0];
	if(by) (*by) = fs->dsy.ncoor[0];
	if(bz) (*bz) = fs->dsz.ncoor[0];

	if(ex) (*ex) = fs->dsx.ncoor[fs->dsx.ncels];
	if(ey) (*ey) = fs->dsy.ncoor[fs->dsy.ncels];
	if(ez) (*ez) = fs->dsz.ncoor[fs->dsz.ncels];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGetGlobalBox"
PetscErrorCode FDSTAGGetGlobalBox(
	FDSTAG      *fs,
	PetscScalar *bx,
	PetscScalar *by,
	PetscScalar *bz,
	PetscScalar *ex,
	PetscScalar *ey,
	PetscScalar *ez)
{
	PetscFunctionBegin;

	if(bx) (*bx) = fs->dsx.gcrdbeg;
	if(by) (*by) = fs->dsy.gcrdbeg;
	if(bz) (*bz) = fs->dsz.gcrdbeg;

	if(ex) (*ex) = fs->dsx.gcrdend;
	if(ey) (*ey) = fs->dsy.gcrdend;
	if(ez) (*ez) = fs->dsz.gcrdend;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGSaveGrid"
PetscErrorCode FDSTAGSaveGrid(FDSTAG *fs)
{
	int            fid;
	char           *fname;
	PetscScalar    *xc, *yc, *zc, chLen;
	PetscMPIInt    rank;
	PetscLogDouble t;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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

		PetscBinaryWrite(fid, &fs->dsx.nproc, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &fs->dsy.nproc, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &fs->dsz.nproc, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &fs->dsx.tnods, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &fs->dsy.tnods, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &fs->dsz.tnods, 1,               PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, fs->dsx.starts, fs->dsx.nproc+1, PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, fs->dsy.starts, fs->dsy.nproc+1, PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, fs->dsz.starts, fs->dsz.nproc+1, PETSC_INT,    PETSC_FALSE);
		PetscBinaryWrite(fid, &chLen,         1,               PETSC_SCALAR, PETSC_FALSE);
		PetscBinaryWrite(fid, xc,             fs->dsx.tnods,   PETSC_SCALAR, PETSC_FALSE);
		PetscBinaryWrite(fid, yc,             fs->dsy.tnods,   PETSC_SCALAR, PETSC_FALSE);
		PetscBinaryWrite(fid, zc,             fs->dsz.tnods,   PETSC_SCALAR, PETSC_FALSE);

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
#undef __FUNCT__
#define __FUNCT__ "DMDACreate3dSetUp"
PetscErrorCode DMDACreate3dSetUp(MPI_Comm comm,
	DMBoundaryType bx, DMBoundaryType by, DMBoundaryType bz, DMDAStencilType stencil_type,
	PetscInt M, PetscInt N, PetscInt P, PetscInt m, PetscInt n, PetscInt p,
	PetscInt dof, PetscInt s, const PetscInt lx[], const PetscInt ly[], const PetscInt lz[], DM *da)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = DMDACreate3d(comm, bx, by, bz, stencil_type, M, N, P, m, n, p, dof, s, lx, ly, lz, da); CHKERRQ(ierr);

	ierr = DMSetFromOptions((*da)); CHKERRQ(ierr);
	ierr = DMSetUp((*da));          CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


/*
// Split points into slots according to weights
// Each slot gets at least one point
void splitPointSlot(
	PetscInt     points,
	PetscInt     n,
	PetscScalar *weights,
	PetscInt    *slots)
{
	PetscScalar min, max, diff;
	PetscInt    i, res, jj;

	// compute initial distribution & residual
	res = points;
	for(i = 0; i < n; i++)
	{	slots[i] = 1 + (PetscInt)round((double)(points-n) * (double)weights[i]);
		res -= slots[i];
	}
	// distribute residual
	while(res)
	{	// add one point to a slot with smallest residual occupation
		if(res > 0)
		{	min = 1.0;
			jj  = 0;
			for(i = 0; i < n; i++)
			{	diff = (PetscScalar)slots[i]/(PetscScalar)(points-res) - weights[i];
				if(diff < min)
				{	min = diff;
					jj  = i;
				}
			}
			slots[jj]++;
			res--;
		}
		// subtract one point from a slot with largest residual occupation
		// except the slots that already have only one point
		if(res < 0)
		{	max = -1.0;
			jj  = 0;
			for(i = 0; i < n; i++)
			{	if(slots[i] == 1) continue;
				diff = (PetscScalar)slots[i]/(PetscScalar)(points-res) - weights[i];
				if(diff > max)
				{	max = diff;
					jj  = i;
				}
			}
			slots[jj]--;
			res++;
		}
	}
}
*/
//---------------------------------------------------------------------------
