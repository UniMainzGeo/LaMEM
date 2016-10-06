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
#include "tools.h"
//---------------------------------------------------------------------------
// * unify coupled & decoupled indexing objects
//---------------------------------------------------------------------------
// MeshSeg1D functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DCreate"
PetscErrorCode MeshSeg1DCreate(
	MeshSeg1D  *ms,
	PetscScalar beg,
	PetscScalar end,
	PetscInt    tncels,
	MeshSegInp *msi)
{
	PetscInt i, istart;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set number of segments
	if(msi->nsegs) ms->nsegs = msi->nsegs;
	else           ms->nsegs = 1;

	// allocate space
	ierr = makeIntArray (&ms->istart, NULL, ms->nsegs+1); CHKERRQ(ierr);
	ierr = makeScalArray(&ms->xstart, NULL, ms->nsegs+1); CHKERRQ(ierr);
	ierr = makeScalArray(&ms->biases, NULL, ms->nsegs);   CHKERRQ(ierr);

	// expand the input segments
	if(msi->nsegs)
	{
		// coordinate & index delimiters
		for(i = 0, istart=0; i < msi->nsegs-1; i++)
		{
			istart += msi->ncells[i];
			ms->istart[i+1] = istart;
			ms->xstart[i+1] = msi->delims[i];
		}

		// biases
		for(i = 0; i < msi->nsegs; i++) ms->biases[i] = msi->biases[i];

	}
	// create uniform mesh by default
	else
	{
		ms->biases[0] = 1.0;
	}

	// set mesh boundaries
	ms->istart[0]         = 0;
	ms->istart[ms->nsegs] = tncels;
	ms->xstart[0]         = beg;
	ms->xstart[ms->nsegs] = end;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DDestroy"
PetscErrorCode MeshSeg1DDestroy(MeshSeg1D *ms)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscFree(ms->istart); CHKERRQ(ierr);
	ierr = PetscFree(ms->xstart); CHKERRQ(ierr);
	ierr = PetscFree(ms->biases); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DStretch"
PetscErrorCode MeshSeg1DStretch(MeshSeg1D *ms, PetscScalar eps)
{
	// Stretch grid with constant stretch factor about coordinate origin.
	// x_new = x_old + eps*x_old

	PetscInt i;

	PetscFunctionBegin;

	// recompute (stretch) coordinates
	for(i = 0; i < ms->nsegs+1; i++) ms->xstart[i] *= (1.0 + eps);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#define SAVENODE(pstart, inode, cnt, n, crd, x) \
	{ if(inode >= pstart && inode < pstart + n) { crd[cnt++] = x; if(cnt == n) return 0; } inode++; }
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeshSeg1DGenCoord"
PetscErrorCode MeshSeg1DGenCoord(
	MeshSeg1D   *ms,     // segments description
	PetscInt     pstart, // starting node index
	PetscInt     n,      // number of nodes to be generated
	PetscScalar *crd)    // coordinates of the nodes
{
	PetscInt    i, ns, iseg, inode, cnt, sum;
	PetscScalar x, dx, xbeg, xend, bias, avgSz, begSz, endSz;

	inode = 0;
	cnt   = 0;

	for(iseg = 0; iseg < ms->nsegs; iseg++)
	{
		// coordinate bounds
		xbeg = ms->xstart[iseg];
		xend = ms->xstart[iseg+1];

		// number of nodes & cells in the segment
		ns = ms->istart[iseg+1] - ms->istart[iseg];

		// bias coefficient (last to first cell size ratio > 1 -> growing)
		bias = ms->biases[iseg];

		// average cell size
		avgSz = (xend - xbeg)/(PetscScalar)(ns);

		// first node
		SAVENODE(pstart, inode, cnt, n, crd, xbeg);

		// uniform case
		if(bias == 1.0)
		{
			// generate coordinates of local nodes
			for(i = 1; i < ns; i++)
			{
				x = xbeg + (PetscScalar)i*avgSz;

				SAVENODE(pstart, inode, cnt, n, crd, x);
			}
		}
		// non-uniform case
		else
		{
			// cell size limits
			begSz = 2.0*avgSz/(1.0 + bias);
			endSz = bias*begSz;

			// cell size increment (negative for bias < 1)
			dx = (endSz - begSz)/(PetscScalar)(ns-1);

			// generate coordinates of local nodes
			for(i = 1, sum = 0; i < ns; i++)
			{
				x = xbeg + (PetscScalar)i*begSz + (PetscScalar)sum*dx;
				sum += i;

				SAVENODE(pstart, inode, cnt, n, crd, x);
			}
		}
	}

	// last node
	SAVENODE(pstart, inode, cnt, n, crd, xend);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar MeshSeg1DGetUniStep(MeshSeg1D *ms)
{
	return (ms->xstart[ms->nsegs] - ms->xstart[0])/(PetscScalar)ms->istart[ms->nsegs];
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

	// column color
	ds->color = (PetscMPIInt) color;

	// column communicator
	ds->comm = MPI_COMM_NULL;

	// global rank of previous process (-1 if none)
	ds->grprev = grprev;

	// global rank of next process (-1 if none)
	ds->grnext = grnext;

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
#define __FUNCT__ "Discret1DGenCoord"
PetscErrorCode Discret1DGenCoord(Discret1D *ds, MeshSeg1D *ms)
{
	PetscInt     i, n, pstart;
	PetscScalar *crd, A, B, C;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// compute number of nodes to be generated locally
	crd    = ds->ncoor;
	n      = ds->nnods;
	pstart = ds->pstart;

	// correct numbers if we need to include internal ghost points
	if(ds->grprev != -1) { pstart--; crd--; n++; }
	if(ds->grnext != -1) { n += 2; }

	ierr = MeshSeg1DGenCoord(ms, pstart, n, crd); CHKERRQ(ierr);

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

	// compute extreme cell sizes
	ierr = Discret1DGetMinMaxCellSize(ds, ms); CHKERRQ(ierr);

	// exchange domain bounds
	ierr = Discret1DExcahngeBounds(ds); CHKERRQ(ierr);

	// free column communicator (optimization)
	ierr = Discret1DFreeColumnComm(ds); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGetMinMaxCellSize"
PetscErrorCode Discret1DGetMinMaxCellSize(Discret1D *ds, MeshSeg1D *ms)
{
	// globally compute and set extreme cell sizes

	PetscInt    i;
	PetscScalar h, sz, lminsz, lmaxsz, gminsz, gmaxsz, rtol;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create communicator if necessary
	ierr = Discret1DGetColumnComm(ds);  CHKERRQ(ierr);

	// set tolerance
	rtol = 1e-8;

	// get local values
	lminsz = lmaxsz = SIZE_CELL(0, 0, (*ds));

	for(i = 1; i < ds->ncels; i++)
	{
		sz = SIZE_CELL(i, 0, (*ds));

		if(sz < lminsz) lminsz = sz;
		if(sz > lmaxsz) lmaxsz = sz;
	}

	// sort out sequential case
	if(ds->nproc == 1)
	{
		gminsz = lminsz;
		gmaxsz = lmaxsz;
	}
	else
	{	// synchronize
		ierr = MPI_Allreduce(&lminsz, &gminsz, 1, MPIU_SCALAR, MPI_MIN, ds->comm); CHKERRQ(ierr);
		ierr = MPI_Allreduce(&lmaxsz, &gmaxsz, 1, MPIU_SCALAR, MPI_MAX, ds->comm); CHKERRQ(ierr);
	}

	// detect uniform mesh & store result
	h = MeshSeg1DGetUniStep(ms);

	if(PetscAbsScalar(gmaxsz-gminsz) < rtol*h)
	{
		ds->h_uni = h;
		ds->h_min = h;
		ds->h_max = h;
	}
	else
	{
		ds->h_uni = -1.0;
		ds->h_min =  gminsz;
		ds->h_max =  gmaxsz;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DExcahngeBounds"
PetscErrorCode Discret1DExcahngeBounds(Discret1D *ds)
{
	// exchange coordinate bounds to be exactly the same on neighboring processors

	MPI_Request srequest, rrequest;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create communicator if necessary
	ierr = Discret1DGetColumnComm(ds); CHKERRQ(ierr);

	// set coordinate bounds
	ds->crdbeg = ds->ncoor[0];
	ds->crdend = ds->ncoor[ds->ncels];

	// exchange coordinate bounds
	if(ds->grnext != -1)
	{
		ierr = MPI_Isend(&ds->crdend, 1, MPIU_SCALAR, ds->rank+1, 0, ds->comm, &srequest); CHKERRQ(ierr);
	}

	if(ds->grprev != -1)
	{
		ierr = MPI_Irecv(&ds->crdbeg, 1, MPIU_SCALAR, ds->rank-1, 0, ds->comm, &rrequest); CHKERRQ(ierr);
	}

	if(ds->grnext != -1)
	{
		ierr =  MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	}

	if(ds->grprev != -1)
	{
		ierr =  MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DStretch"
PetscErrorCode Discret1DStretch(Discret1D *ds, MeshSeg1D *ms, PetscScalar eps)
{
	// stretch grid with constant stretch factor about coordinate origin.
	// x_new = x_old + eps*x_old

	PetscInt    i;
	PetscScalar h;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// recompute segment data
	ierr = MeshSeg1DStretch(ms, eps); CHKERRQ(ierr);

	// recompute (stretch) node coordinates in the buffer
	for(i = 0; i < ds->bufsz; i++) ds->nbuff[i] *= (1.0 + eps);

	// recompute cell coordinates
	for(i = -1; i < ds->ncels+1; i++)
		ds->ccoor[i] = (ds->ncoor[i] + ds->ncoor[i+1])/2.0;

	// update mesh steps
	if(ds->h_uni < 0.0)
	{
		ds->h_min *= (1.0 + eps);
		ds->h_max *= (1.0 + eps);
	}
	else
	{
		h         = MeshSeg1DGetUniStep(ms);
		ds->h_uni = h;
		ds->h_min = h;
		ds->h_max = h;
	}

	// recompute coordinate bounds
	ds->crdbeg *= (1.0 + eps);
	ds->crdend *= (1.0 + eps);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DView"
PetscErrorCode Discret1DView(Discret1D *ds, const char *name)
{
	PetscInt    i;
	PetscMPIInt grank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &grank); CHKERRQ(ierr);
	ierr = PetscPrintf (PETSC_COMM_WORLD, "\n\n\n === %s ===\n", name); CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "   ***grank=%lld  rank=%lld  start=%lld  tnods=%lld  nnods=%lld  ncels=%lld\n",
		(LLD)grank, (LLD)ds->rank, (LLD)ds->starts[ds->rank], (LLD)ds->tnods, (LLD)ds->nnods, (LLD)ds->ncels); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
	ierr = PetscPrintf (PETSC_COMM_WORLD, "====================================================\n\n\n"); CHKERRQ(ierr);

	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "   ***grank=%lld  ", (LLD)grank); CHKERRQ(ierr);

	for(i = 0; i < ds->bufsz; i++)
	{	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%f ", ds->nbuff[i]); CHKERRQ(ierr);
	}
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
	ierr = PetscPrintf (PETSC_COMM_WORLD, "====================================================\n\n\n"); CHKERRQ(ierr);

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
#define __FUNCT__ "FDSTAGClear"
PetscErrorCode FDSTAGClear(FDSTAG  *fs)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(fs, sizeof(FDSTAG)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCreate"
PetscErrorCode FDSTAGCreate(
	FDSTAG  *fs,
	PetscInt Nx, PetscInt Ny, PetscInt Nz)
{
	// Create object with all necessary arrays to handle FDSTAG discretization.

	// NOTE: velocity components have one layer of boundary ghost points.
	// The idea is that velocity vectors should contain sufficient information
	// to compute strain/rates/stresses/residuals including boundary conditions.

	PetscInt         nnx, nny, nnz;
	PetscInt         ncx, ncy, ncz;
	PetscInt         ndof, nlayer;
	const PetscInt  *plx, *ply, *plz;
	PetscInt        *lx,  *ly,  *lz;
	PetscInt         rx,   ry,   rz;
	PetscInt         cx,   cy,   cz;
	PetscInt         Px,   Py,   Pz;

	PetscErrorCode 	 ierr;
	PetscFunctionBegin;

	ndof   = 1;
	nlayer = 1;

	// get number of processors
	Px = PETSC_DECIDE;
	Py = PETSC_DECIDE;
	Pz = PETSC_DECIDE;

	PetscOptionsGetInt(NULL, NULL, "-cpu_x", &Px, NULL); // fix # of processors in x-direction
	PetscOptionsGetInt(NULL, NULL, "-cpu_y", &Py, NULL); // fix # of processors in y-direction
	PetscOptionsGetInt(NULL, NULL, "-cpu_z", &Pz, NULL); // fix # of processors in z-direction

	// partition central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz-1, Px, Py, Pz, ndof, nlayer, 0, 0, 0, &fs->DA_CEN); CHKERRQ(ierr);

	// get actual number of processors in every direction (can be different compared to given)
	ierr = DMDAGetInfo(fs->DA_CEN, 0, 0, 0, 0, &Px, &Py, &Pz, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	// get number of points per processor
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &plx, &ply, &plz); CHKERRQ(ierr);
	ierr = makeIntArray(&lx, plx, Px); CHKERRQ(ierr);
	ierr = makeIntArray(&ly, ply, Py); CHKERRQ(ierr);
	ierr = makeIntArray(&lz, plz, Pz); CHKERRQ(ierr);

	// increment number of points on last processor
	lx[Px-1]++; ly[Py-1]++; lz[Pz-1]++;

	// corners (DA_COR) no boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_COR); CHKERRQ(ierr);

	// XY edges (DA_XY) no boundary ghost points (1-layer stencil box)
	lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz-1, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_XY); CHKERRQ(ierr);
	lz[Pz-1]++;

	// XZ edges (DA_XZ) no boundary ghost points (1-layer stencil box)
	ly[Py-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_XZ); CHKERRQ(ierr);
	ly[Py-1]++;

	// YZ edges (DA_YZ) no boundary ghost points (1-layer stencil box)
	lx[Px-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_YZ); CHKERRQ(ierr);
	lx[Px-1]++;


	// X face (DA_X) with boundary ghost points (1-layer stencil box)
	ly[Py-1]--; lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz-1, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_X); CHKERRQ(ierr);
	ly[Py-1]++; lz[Pz-1]++;

	// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz-1, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_Y); CHKERRQ(ierr);
	lx[Px-1]++; lz[Pz-1]++;

	// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; ly[Py-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz, Px, Py, Pz, ndof, nlayer, lx, ly, lz, &fs->DA_Z); CHKERRQ(ierr);
	lx[Px-1]++; ly[Py-1]++;

	// get processor ranks
	ierr = DMDAGetProcessorRank(fs->DA_CEN, &rx, &ry, &rz, 0); CHKERRQ(ierr);

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

	// clear temporary storage
	ierr = PetscFree(lx);  CHKERRQ(ierr);
	ierr = PetscFree(ly);  CHKERRQ(ierr);
	ierr = PetscFree(lz);  CHKERRQ(ierr);

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

	// setup indexing data
	ierr = DOFIndexCreate(&fs->dof, fs->DA_CEN, fs->DA_X, fs->DA_Y, fs->DA_Z); CHKERRQ(ierr);

	// compute number of local and ghost points
	nnx = fs->dsx.nnods+2; ncx = fs->dsx.ncels+2;
	nny = fs->dsy.nnods+2; ncy = fs->dsy.ncels+2;
	nnz = fs->dsz.nnods+2; ncz = fs->dsz.ncels+2;
/*
	fs->nCellsGh = ncx*ncy*ncz;
	fs->nXFaceGh = nnx*ncy*ncz;
	fs->nYFaceGh = ncx*nny*ncz;
	fs->nZFaceGh = ncx*ncy*nnz;

	// compute number of local & ghost DOF
	fs->numdofGh = fs->nXFaceGh + fs->nYFaceGh + fs->nZFaceGh + fs->nCellsGh;
*/

	// get ranks of neighbor processes
	ierr =  FDSTAGGetNeighbProc(fs); CHKERRQ(ierr);

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

	// destroy mesh segment data
	ierr = MeshSeg1DDestroy(&fs->msx); CHKERRQ(ierr);
	ierr = MeshSeg1DDestroy(&fs->msy); CHKERRQ(ierr);
	ierr = MeshSeg1DDestroy(&fs->msz); CHKERRQ(ierr);

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
#define __FUNCT__ "FDSTAGGenCoord"
PetscErrorCode FDSTAGGenCoord(FDSTAG *fs, UserCtx *usr)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// generate mesh segment data
	ierr = MeshSeg1DCreate(&fs->msx, usr->x_left,  usr->x_left  + usr->W, usr->nel_x, &usr->mseg_x); CHKERRQ(ierr);
	ierr = MeshSeg1DCreate(&fs->msy, usr->y_front, usr->y_front + usr->L, usr->nel_y, &usr->mseg_y); CHKERRQ(ierr);
	ierr = MeshSeg1DCreate(&fs->msz, usr->z_bot,   usr->z_bot   + usr->H, usr->nel_z, &usr->mseg_z); CHKERRQ(ierr);

	// generate coordinates
	ierr = Discret1DGenCoord(&fs->dsx, &fs->msx); CHKERRQ(ierr);
	ierr = Discret1DGenCoord(&fs->dsy, &fs->msy); CHKERRQ(ierr);
	ierr = Discret1DGenCoord(&fs->dsz, &fs->msz); CHKERRQ(ierr);

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

	PetscInt  rx, ry, rz;

	PetscFunctionBegin;

	GET_POINT_RANK(X[0], rx, fs->dsx);
	GET_POINT_RANK(X[1], ry, fs->dsy);
	GET_POINT_RANK(X[2], rz, fs->dsz);

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

	PetscScalar maxAspRat;
	PetscInt    px, py, pz, cx, cy, cz, nx, ny, nz, nVelDOF, nCells;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	px = fs->dsx.nproc;  cx = fs->dsx.tcels;  nx = fs->dsx.tnods;
	py = fs->dsy.nproc;  cy = fs->dsy.tcels;  ny = fs->dsy.tnods;
	pz = fs->dsz.nproc;  cz = fs->dsz.tcels;  nz = fs->dsz.tnods;

	nCells  = cx*cy*cz;
	nVelDOF = nx*cy*cz + cx*ny*cz + cx*cy*nz;

	ierr = FDSTAGGetAspectRatio(fs, &maxAspRat); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, " Processor grid  [nx, ny, nz]   : [%lld, %lld, %lld]\n", (LLD)px, (LLD)py, (LLD)pz);
	PetscPrintf(PETSC_COMM_WORLD, " Fine grid cells [nx, ny, nz]   : [%lld, %lld, %lld]\n", (LLD)cx, (LLD)cy, (LLD)cz);
	PetscPrintf(PETSC_COMM_WORLD, " Number of cells                :  %lld\n", (LLD)nCells);
	PetscPrintf(PETSC_COMM_WORLD, " Number of velocity DOF         :  %lld\n", (LLD)nVelDOF);
	PetscPrintf(PETSC_COMM_WORLD, " Maximum cell aspect cell ratio :  %7.5f\n", maxAspRat);

	if(maxAspRat > 1e100) PetscPrintf(PETSC_COMM_WORLD, " hmm... sorry if it bothers anyone, but IMHO the aspect ratio is a bit too large ...\n");

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

	if(bx) (*bx) = fs->msx.xstart[0];
	if(by) (*by) = fs->msy.xstart[0];
	if(bz) (*bz) = fs->msz.xstart[0];

	if(ex) (*ex) = fs->msx.xstart[fs->msx.nsegs];
	if(ey) (*ey) = fs->msy.xstart[fs->msy.nsegs];
	if(ez) (*ez) = fs->msz.xstart[fs->msz.nsegs];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGProcPartitioning"
PetscErrorCode FDSTAGProcPartitioning(FDSTAG *fs, PetscScalar chLen)
{
	int         fid;
	char        *fname;
	PetscScalar *xc, *yc, *zc;
    PetscMPIInt rank;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
	PetscPrintf(PETSC_COMM_WORLD,"# Save processor partitioning \n");

	// gather global coord
	ierr = Discret1DGatherCoord(&fs->dsx, &xc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsy, &yc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsz, &zc); CHKERRQ(ierr);

	if(rank == 0)
	{
        PetscPrintf(PETSC_COMM_SELF,"# Save processor partitioning file on rank 0 \n");
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
    MPI_Barrier(PETSC_COMM_WORLD);
    
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
