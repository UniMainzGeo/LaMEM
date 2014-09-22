//---------------------------------------------------------------------------
//........   PARALLEL STAGGERED GRID USING PETSC DISTRIBUTED ARRAYS  ........
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "Utils.h"
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
	PetscInt       i, n, cnt;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	memset(ds, 0, sizeof(ds));

	// number of processors
	ds->nproc = nproc;

	// rank of current processor
	ds->rank = rank;

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
	if(grnext != -1) n = ds->nnods+3;
	else             n = ds->nnods+2;
	ierr = makeScalArray(&ds->nbuff, 0, n); CHKERRQ(ierr);
	ds->ncoor = ds->nbuff + 1;

	// coordinates of local cells + 1 layer (both sides) of ghost points
	ierr = makeScalArray(&ds->cbuff, 0, ds->ncels+2); CHKERRQ(ierr);
	ds->ccoor = ds->cbuff + 1;

	// column color
	ds->color = (PetscMPIInt) color;

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
	ierr = PetscFree(ds->nbuff);  CHKERRQ(ierr);
	ierr = PetscFree(ds->cbuff);  CHKERRQ(ierr);
	ierr = PetscFree(ds->starts); CHKERRQ(ierr);

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DView"
PetscErrorCode Discret1DView(Discret1D *ds, const char *name)
{
	PetscInt    i, n;
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

	n = ds->nnods+1; if(ds->grnext != -1) n++;

	for(i = -1; i < n; i++)
	{	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%f ", ds->ncoor[i]); CHKERRQ(ierr);
	}
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
	ierr = PetscPrintf (PETSC_COMM_WORLD, "====================================================\n\n\n"); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGetMinCellSize"
PetscErrorCode Discret1DGetMinCellSize(Discret1D *ds, PetscScalar *sz)
{
	// compute and return the size of smallest local cell

	PetscInt    i;
	PetscScalar lsz, msz;

	PetscFunctionBegin;

	msz = SIZE_CELL(0, 0, (*ds));

	for(i = 1; i < ds->ncels; i++)
	{
		lsz = SIZE_CELL(i, 0, (*ds));
		if(lsz < msz) msz = lsz;
	}

	(*sz) = msz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Discret1DGetColumnComm"
PetscErrorCode Discret1DGetColumnComm(Discret1D *ds, MPI_Comm *comm)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_split(PETSC_COMM_WORLD, ds->color, ds->rank, comm); CHKERRQ(ierr);

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

	MPI_Comm     comm;
	PetscScalar *pcoord   = NULL;
	PetscInt    *recvcnts = NULL, i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check for sequential case
	if(ds->nproc == 1)
	{
		// return coordinates on rank zero of PETSC_COMM_WORLD
		if(ISRankZero(PETSC_COMM_WORLD))
		{
			ierr = makeScalArray(&pcoord, NULL, ds->tnods); CHKERRQ(ierr);

			for(i = 0; i < ds->tnods; i++) pcoord[i] = ds->ncoor[i];

			(*coord) = pcoord;
		}

		PetscFunctionReturn(0);
	}

	// create column communicator
	ierr = Discret1DGetColumnComm(ds, &comm); CHKERRQ(ierr);

	// gather coordinates on zero ranks of column communicator
	if(ISRankZero(comm))
	{
		// allocate coordinates
		ierr = makeScalArray(&pcoord, NULL, ds->tnods); CHKERRQ(ierr);

		// allocate receive counts
		ierr = makeIntArray(&recvcnts, NULL, ds->nproc); CHKERRQ(ierr);

		// compute receive counts
		for(i = 0; i < ds->nproc-1; i++) recvcnts[i] = ds->starts[i+1] - ds->starts[i];

		// last element stores index of the last node (not total number of nodes)
		recvcnts[ds->nproc-1] = ds->starts[ds->nproc] - ds->starts[ds->nproc-1] + 1;
	}

	// gather coordinates
	ierr = MPI_Gatherv(ds->ncoor, ds->nnods, MPIU_SCALAR, pcoord, recvcnts, ds->starts, MPIU_SCALAR, 0, comm); CHKERRQ(ierr);

	// return coordinates on rank zero of PETSC_COMM_WORLD
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		(*coord) = pcoord;
	}
	else
	{
		ierr = PetscFree(pcoord); CHKERRQ(ierr);
	}

	ierr = MPI_Comm_free(&comm); CHKERRQ(ierr);

	ierr = PetscFree(recvcnts); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// DOFIndex functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DOFIndexCreate"
PetscErrorCode DOFIndexCreate(DOFIndex *id, FDSTAG *fs, idxtype idxmod)
{

	// compute & set global indices of local & ghost nodes

	//**********************************************************************
	// NOTE:
	// for the ghost points, store negative global index of the primary DOF
	// instead of -1
	//**********************************************************************

	PetscInt    ln, lnp, start, startp, sum;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	Vec         givx, givy, givz, gip;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// store index mode
	id->idxmod = idxmod;

	// set local number DOF
	if(idxmod == IDXCOUPLED)
	{
		lnp    = -1;
		startp = -1;

		ln = fs->nXFace + fs->nYFace + fs->nZFace + fs->nCells;
	}

	if(idxmod == IDXUNCOUPLED)
	{
		lnp = fs->nCells;

		ierr = MPI_Scan(&lnp, &sum, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

		startp = sum - lnp;

		ln = fs->nXFace + fs->nYFace + fs->nZFace;

	}

	// compute prefix sum
	ierr = MPI_Scan(&ln, &sum, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

	// compute index of the first DOF on this processor
	start = sum - ln;

	// store number DOF & first indices
	id->numdof  = ln;
	id->istart  = start;
	id->numdofp = lnp;
	id->istartp = startp;

	// create index vectors
	ierr = DMCreateGlobalVector(fs->DA_X,   &givx); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Y,   &givy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_Z,   &givz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_CEN, &gip);  CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   givx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   givy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   givz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, gip,  &ip);   CHKERRQ(ierr);

	//=======================================================
	// compute interlaced global numbering of the local nodes
	//=======================================================

	//---------
	// X-points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ivx[k][j][i] = (PetscScalar)start; start++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ivy[k][j][i] = (PetscScalar)start; start++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		ivz[k][j][i] = (PetscScalar)start; start++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	if(idxmod == IDXCOUPLED)
	{
		startp = start;
	}

	START_STD_LOOP
	{
		ip[k][j][i] = (PetscScalar)startp; startp++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   givx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   givy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   givz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, gip,  &ip);   CHKERRQ(ierr);

	// create index vectors (ghosted)
	ierr = DMCreateLocalVector(fs->DA_X,   &id->ivx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &id->ivy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &id->ivz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &id->ip);  CHKERRQ(ierr);

	// set global indices of the local and ghost nodes (including boundary)
	ierr = VecSet(id->ivx, -1.0); CHKERRQ(ierr);
	ierr = VecSet(id->ivy, -1.0); CHKERRQ(ierr);
	ierr = VecSet(id->ivz, -1.0); CHKERRQ(ierr);
	ierr = VecSet(id->ip,  -1.0); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(fs->DA_X,   givx, id->ivx)
	GLOBAL_TO_LOCAL(fs->DA_Y,   givy, id->ivy)
	GLOBAL_TO_LOCAL(fs->DA_Z,   givz, id->ivz)
	GLOBAL_TO_LOCAL(fs->DA_CEN, gip,  id->ip)

	// destroy index vectors
	ierr = VecDestroy(&givx); CHKERRQ(ierr);
	ierr = VecDestroy(&givy); CHKERRQ(ierr);
	ierr = VecDestroy(&givz); CHKERRQ(ierr);
	ierr = VecDestroy(&gip);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DOFIndexDestroy"
PetscErrorCode DOFIndexDestroy(DOFIndex *id)
{
	PetscErrorCode 	 ierr;
	PetscFunctionBegin;

	// destroy index vectors (ghosted)
	ierr = VecDestroy(&id->ivx); CHKERRQ(ierr);
	ierr = VecDestroy(&id->ivy); CHKERRQ(ierr);
	ierr = VecDestroy(&id->ivz); CHKERRQ(ierr);
	ierr = VecDestroy(&id->ip);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// FDSTAG functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCreate"
PetscErrorCode FDSTAGCreate(
	FDSTAG  *fs,
	PetscInt Nx, PetscInt Ny, PetscInt Nz,
	PetscInt Px, PetscInt Py, PetscInt Pz)
{
	// Create object with all necessary arrays to handle FDSTAG discretization.

	// NOTE: velocity components have one layer of boundary ghost points.
	// The idea is that velocity vectors should contain sufficient information
	// to compute strain/rates/stresses/residuals including boundary conditions.

	PetscInt         nnx, nny, nnz;
	PetscInt         ncx, ncy, ncz;
	PetscInt         dof, nlayer;
	const PetscInt  *plx, *ply, *plz;
	PetscInt        *lx,  *ly,  *lz;
	PetscInt         rx,   ry,   rz;
	PetscInt         cx,   cy,   cz;

	PetscErrorCode 	 ierr;
	PetscFunctionBegin;

	dof    = 1;
	nlayer = 1;

	// get number of processors
	Px = PETSC_DECIDE;
	Py = PETSC_DECIDE;
	Pz = PETSC_DECIDE;

//	PetscOptionsGetInt(PETSC_NULL, "-da_proc_x", &Px, PETSC_NULL); // fix # of processors in x-direction
//	PetscOptionsGetInt(PETSC_NULL, "-da_proc_y", &Py, PETSC_NULL); // fix # of processors in y-direction
//	PetscOptionsGetInt(PETSC_NULL, "-da_proc_z", &Pz, PETSC_NULL); // fix # of processors in z-direction

	// partition central points (DA_CEN) with boundary ghost points (1-layer stencil box)
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz-1, Px, Py, Pz, dof, 1, 0, 0, 0, &fs->DA_CEN); CHKERRQ(ierr);

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
		Nx, Ny, Nz, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_COR); CHKERRQ(ierr);

	// XY edges (DA_XY) no boundary ghost points (1-layer stencil box)
	lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny, Nz-1, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_XY); CHKERRQ(ierr);
	lz[Pz-1]++;

	// XZ edges (DA_XZ) no boundary ghost points (1-layer stencil box)
	ly[Py-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_XZ); CHKERRQ(ierr);
	ly[Py-1]++;

	// YZ edges (DA_YZ) no boundary ghost points (1-layer stencil box)
	lx[Px-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_YZ); CHKERRQ(ierr);
	lx[Px-1]++;


	// X face (DA_X) with boundary ghost points (1-layer stencil box)
	ly[Py-1]--; lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx, Ny-1, Nz-1, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_X); CHKERRQ(ierr);
	ly[Py-1]++; lz[Pz-1]++;

	// Y face (DA_Y) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; lz[Pz-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny, Nz-1, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_Y); CHKERRQ(ierr);
	lx[Px-1]++; lz[Pz-1]++;

	// Z face (DA_Z) with boundary ghost points (1-layer stencil box)
	lx[Px-1]--; ly[Py-1]--;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
		Nx-1, Ny-1, Nz, Px, Py, Pz, dof, nlayer, lx, ly, lz, &fs->DA_Z); CHKERRQ(ierr);
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
	ierr = DOFIndexCreate(&fs->cdof, fs, IDXCOUPLED);   CHKERRQ(ierr);
	ierr = DOFIndexCreate(&fs->udof, fs, IDXUNCOUPLED); CHKERRQ(ierr);

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
	ierr = DOFIndexDestroy(&fs->cdof); CHKERRQ(ierr);
	ierr = DOFIndexDestroy(&fs->udof); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGGenCoord"
PetscErrorCode FDSTAGGenCoord(FDSTAG *fs, UserContext *usr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// generate mesh segment data
	ierr = MeshSeg1DCreate(&fs->msx, usr->x_left,  usr->x_left  + usr->W, usr->nel_x, &usr->mseg_x); CHKERRQ(ierr);
	ierr = MeshSeg1DCreate(&fs->msy, usr->y_front, usr->y_front + usr->L, usr->nel_y, &usr->mseg_y); CHKERRQ(ierr);
	ierr = MeshSeg1DCreate(&fs->msz, usr->z_bot,   usr->z_bot   + usr->H, usr->nel_z, &usr->mseg_z); CHKERRQ(ierr);

	// generate coordinates
	ierr = Discret1DGenCoord(&fs->dsx, &fs->msx);
	ierr = Discret1DGenCoord(&fs->dsy, &fs->msy);
	ierr = Discret1DGenCoord(&fs->dsz, &fs->msz);

	// compute sizes of the smallest cells
	ierr = FDSTAGGetMinCellSize(fs);

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
#define __FUNCT__ "FDSTAGGetMinCellSize"
PetscErrorCode FDSTAGGetMinCellSize(FDSTAG *fs)
{
	// compute minimum cell size in each direction
	PetscMPIInt nproc;
	PetscScalar lsz[3], gsz[3];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get number of processors
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);

	// compute local cell sizes
	ierr = Discret1DGetMinCellSize(&fs->dsx, &lsz[0]); CHKERRQ(ierr);
	ierr = Discret1DGetMinCellSize(&fs->dsy, &lsz[1]); CHKERRQ(ierr);
	ierr = Discret1DGetMinCellSize(&fs->dsz, &lsz[2]); CHKERRQ(ierr);

	// compute & set global cell sizes
	if(nproc != 1)
	{
		// exchange
		ierr = MPI_Allreduce(lsz, gsz, 3, MPIU_SCALAR, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);

		// store the result
		fs->mdx = gsz[0];
		fs->mdy = gsz[1];
		fs->mdz = gsz[2];
	}
	else
	{	// store the result
		fs->mdx = lsz[0];
		fs->mdy = lsz[1];
		fs->mdz = lsz[2];
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGProcPartitioning"
PetscErrorCode FDSTAGProcPartitioning(FDSTAG *fs, UserContext *user)
{
	PetscScalar *xc,*yc,*zc;
	int         fid;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"# Save processor partitioning \n");

	// gather global coord
	ierr = Discret1DGatherCoord(&fs->dsx,&xc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsy,&yc); CHKERRQ(ierr);
	ierr = Discret1DGatherCoord(&fs->dsz,&zc); CHKERRQ(ierr);

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		// save file
		asprintf(&fname,"ProcessorPartitioning_%lldcpu_%lld.%lld.%lld.bin",(LLD)(fs->dsx.nproc*fs->dsy.nproc*fs->dsz.nproc),(LLD)fs->dsx.nproc,(LLD)fs->dsy.nproc,(LLD)fs->dsz.nproc);
		PetscBinaryOpen(fname,FILE_MODE_WRITE,&fid);
		PetscBinaryWrite(fid,&fs->dsx.nproc,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&fs->dsy.nproc,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&fs->dsz.nproc,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&fs->dsx.tnods,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&fs->dsy.tnods,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&fs->dsz.tnods,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,fs->dsx.starts,fs->dsx.nproc+1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,fs->dsy.starts,fs->dsy.nproc+1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,fs->dsz.starts,fs->dsz.nproc+1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fid,&user->Characteristic.Length,1,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fid,xc,fs->dsx.tnods,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fid,yc,fs->dsy.tnods,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fid,zc,fs->dsz.tnods,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryClose(fid);
		free(fname);

		ierr = PetscFree(xc); CHKERRQ(ierr);
		ierr = PetscFree(yc); CHKERRQ(ierr);
		ierr = PetscFree(zc); CHKERRQ(ierr);
	}

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
