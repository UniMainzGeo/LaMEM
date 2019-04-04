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
 **    filename:   fdstag.h
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
#ifndef __fdstag_h__
#define __fdstag_h__
//---------------------------------------------------------------------------

struct FB;
struct Scaling;

//---------------------------------------------------------------------------
// mesh segments data (temporary structure)
//---------------------------------------------------------------------------

struct MeshSeg1D
{
	PetscInt    nsegs;                    // number of segments
	PetscInt    istart[_max_num_segs_+1]; // indices of the first nodes plus last index
	PetscScalar xstart[_max_num_segs_+1]; // coordinates of the first nodes plus total size
	PetscScalar biases[_max_num_segs_  ]; // biases for each segment
	PetscInt    tcels;                    // total number of cells
	PetscInt    uniform;                  // uniform grid flag

};

//---------------------------------------------------------------------------

PetscErrorCode MeshSeg1DReadParam(
	MeshSeg1D  *ms,
	PetscScalar leng,
	PetscScalar gtol,
	const char *dir,
	FB         *fb);

// (partially) mesh a segment with (optionally) biased element size
PetscErrorCode MeshSeg1DGenCoord(
	MeshSeg1D   *ms,     // segments description
	PetscInt     iseg,   // segment index
	PetscInt     nl,     // number of nodes to be generated
	PetscInt     istart, // index of the first node
	PetscScalar *crd);   // coordinates of the nodes


//---------------------------------------------------------------------------
// finite difference discretization / domain decomposition data for single direction
struct Discret1D
{
	PetscInt      nproc;   // number of processors
	PetscMPIInt   rank;    // rank of current processor

	PetscInt     *starts;  // index of first node (cell) on all processors + last index
	PetscInt      pstart;  // index of first node (cell) on this processors

	PetscInt      tnods;   // total number of nodes
	PetscInt      tcels;   // total number of cells (tnods-1)

	PetscInt      nnods;   // number of local nodes
	PetscInt      ncels;   // number of local cells

	PetscScalar  *ncoor;   // coordinates of local nodes (+ 1 layer of ghost points)
	PetscScalar  *ccoor;   // coordinates of local cells (+ 1 layer of ghost points)
	PetscScalar  *nbuff;   // memory buffer for node coordinates
	PetscScalar  *cbuff;   // memory buffer for cells coordinates
	PetscInt      bufsz;   // size of node buffer

	PetscMPIInt   grprev;  // global rank of previous process (-1 for first processor)
	PetscMPIInt   grnext;  // global rank of next process (-1 for last processor)

	PetscMPIInt   color;   // color of processor column in base direction
	MPI_Comm      comm;    // column communicator

	PetscInt      uniform; // uniform grid flag

	PetscScalar   gcrdbeg; // global grid coordinate bound (begin)
	PetscScalar   gcrdend; // global grid coordinate bound (end)

	PetscScalar   gtol;    // geometric tolerance
};

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
	PetscScalar gtol);     // geometric tolerance

PetscErrorCode Discret1DDestroy(Discret1D *ds);

PetscErrorCode Discret1DReadRestart(Discret1D *ds, FILE *fp);

PetscErrorCode Discret1DWriteRestart(Discret1D *ds, FILE *fp);

// get number of cells per processor
PetscErrorCode Discret1DGetNumCells(Discret1D  *ds, PetscInt **ncelProc);

// generate local coordinates
PetscErrorCode Discret1DGenCoord(Discret1D *ds, MeshSeg1D *ms);

// stretch grid with constant stretch factor about reference point
PetscErrorCode Discret1DStretch(Discret1D *ds,  PetscScalar eps, PetscScalar ref);

// create 1D communicator of the processor column in the base direction
PetscErrorCode Discret1DGetColumnComm(Discret1D *ds);

// destroy 1D communicator
PetscErrorCode Discret1DFreeColumnComm(Discret1D *ds);

// gather coordinate array on rank zero of PETSC_COMM_WORLD
// WARNING! the array only exists on rank zero of PETSC_COMM_WORLD
// WARNING! the array must be destroyed after use!
PetscErrorCode Discret1DGatherCoord(Discret1D *ds, PetscScalar **coord);

// check multigrid restrictions, get maximum number of coarsening steps
PetscErrorCode Discret1DCheckMG(Discret1D *ds, const char *dir, PetscInt *_ncors);

// get maximum inverse time step on local domain (CFL)
PetscErrorCode Discret1DgetMaxInvStep(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax);

// find index of a cell containing point (local points only)
PetscErrorCode Discret1DFindPoint(Discret1D *ds, PetscScalar x, PetscInt &ID);

//---------------------------------------------------------------------------

enum idxtype { IDXNONE, IDXCOUPLED, IDXUNCOUPLED };

// global indexing of the DOF
struct DOFIndex
{
	//=====================================================================
	//
	// index vectors contain global DOF numbers
	// boundary ghost points are marked by -1
	//
	//=====================================================================

	idxtype  idxmod;            // indexing mode
	DM       DA_CEN;            // central points
	DM       DA_X, DA_Y, DA_Z;  // face points
	PetscInt lnvx, lnvy, lnvz;  // local number of DOF
	PetscInt lnv, lnp, ln;      // ...
	PetscInt stv, stp, st;      // starting indices (stv & stp - decoupled layout)
	Vec      ivx, ivy, ivz, ip; // index vectors (ghosted)

};

//---------------------------------------------------------------------------
// DOFIndex functions
//---------------------------------------------------------------------------

PetscErrorCode DOFIndexCreate(DOFIndex *dof, DM DA_CEN, DM DA_X, DM DA_Y, DM DA_Z);

PetscErrorCode DOFIndexDestroy(DOFIndex *dof);

PetscErrorCode DOFIndexCompute(DOFIndex *dof, idxtype idxmod);

//---------------------------------------------------------------------------
// staggered grid data structure
struct FDSTAG
{
	Scaling  *scal;

	// local discretization data (coordinates, indexing & domain decomposition)
	Discret1D dsx;
	Discret1D dsy;
	Discret1D dsz;

	// distributed arrays (partitioning and communication layouts)
	DM DA_CEN;              // central points
	DM DA_COR;              // corner points
	DM DA_XY, DA_XZ, DA_YZ; // edges
	DM DA_X,  DA_Y,  DA_Z;  // face velocities & residuals

	DOFIndex dof; // global variable indexing

	// local number of local grid points
	PetscInt nCells;  // cells
	PetscInt nCorns;  // corners
	PetscInt nXYEdg;  // XY-edges
	PetscInt nXZEdg;  // XZ-edges
	PetscInt nYZEdg;  // YZ-edges
	PetscInt nXFace;  // X-faces
	PetscInt nYFace;  // Y-faces
	PetscInt nZFace;  // Z-faces

	PetscMPIInt neighb[_num_neighb_]; // global ranks of neighboring process

	PetscScalar gtol; // relative geometry tolerance

};

//---------------------------------------------------------------------------
// FDSTAG functions
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGCreate(FDSTAG *fs, FB *fb);

PetscErrorCode FDSTAGDestroy(FDSTAG *fs);

PetscErrorCode FDSTAGReadRestart(FDSTAG *fs, FILE *fp);

PetscErrorCode FDSTAGWriteRestart(FDSTAG *fs, FILE *fp);

PetscErrorCode FDSTAGCreateDMDA(FDSTAG *fs,
	PetscInt  Nx, PetscInt  Ny, PetscInt  Nz,
	PetscInt  Px, PetscInt  Py, PetscInt  Pz,
	PetscInt *lx, PetscInt *ly, PetscInt *lz);

// return an array with the global ranks of adjacent processes (including itself)
PetscErrorCode FDSTAGGetNeighbProc(FDSTAG *fs);

// get local & global ranks of a domain containing a point (only neighbors are checked)
PetscErrorCode FDSTAGGetPointRanks(FDSTAG *fs, PetscScalar *X, PetscInt *lrank, PetscMPIInt *grank);

// compute maximum aspect ratio in the grid
PetscErrorCode FDSTAGGetAspectRatio(FDSTAG *fs, PetscScalar *maxAspRat);

// print & check essential grid details
PetscErrorCode FDSTAGView(FDSTAG *fs);

PetscErrorCode FDSTAGGetLocalBox(FDSTAG *fs,
	PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
	PetscScalar *ex, PetscScalar *ey, PetscScalar *ez);

PetscErrorCode FDSTAGGetGlobalBox(FDSTAG *fs,
	PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
	PetscScalar *ex, PetscScalar *ey, PetscScalar *ez);

// save grid coordinates and processor partitioning to disk
PetscErrorCode FDSTAGSaveGrid(FDSTAG *fs);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// get sub-domain ranks, starting node IDs, and number of nodes for output
#define GET_OUTPUT_RANGE(r, n, s, ds) { r = ds.rank; s = ds.starts[r]; n = ds.starts[r+1] - s + 1; }

// get loop bounds for node discretization
#define GET_NODE_RANGE(n, s, ds) { n = ds.nnods; s = ds.pstart; }

// get loop bounds for cell discretization
#define GET_CELL_RANGE(n, s, ds) { n = ds.ncels; s = ds.pstart; }

// get loop bounds for node discretization (including BOUNDARY ghost points)
#define GET_NODE_RANGE_GHOST_BND(n, s, ds) { n = ds.nnods; s = ds.pstart; if(ds.grprev == -1) { s--; n++; } if(ds.grnext == -1) n++; }

// get loop bounds for cell discretization (including BOUNDARY ghost points)
#define GET_CELL_RANGE_GHOST_BND(n, s, ds) { n = ds.ncels; s = ds.pstart; if(ds.grprev == -1) { s--; n++; } if(ds.grnext == -1) n++; }

// get loop bounds for node discretization (including INTERNAL ghost points)
#define GET_NODE_RANGE_GHOST_INT(n, s, ds) { n = ds.nnods + 2; s = ds.pstart - 1; if(ds.grprev == -1) { s++; n--; } if(ds.grnext == -1) n--; }

// get loop bounds for cell discretization (including INTERNAL ghost points)
#define GET_CELL_RANGE_GHOST_INT(n, s, ds) { n = ds.ncels + 2; s = ds.pstart - 1; if(ds.grprev == -1) { s++; n--; } if(ds.grnext == -1) n--; }

// get loop bounds for node discretization (including ALL ghost points)
#define GET_NODE_RANGE_GHOST_ALL(n, s, ds) { n = ds.nnods + 2; s = ds.pstart - 1; }

// get loop bounds for cell discretization (including ALL ghost points)
#define GET_CELL_RANGE_GHOST_ALL(n, s, ds) { n = ds.ncels + 2; s = ds.pstart - 1; }

//---------------------------------------------------------------------------

// get coordinate of i-th CELL (center)
#define COORD_CELL(i, s, ds) (ds.ccoor[(i-s)])

// get coordinate of i-th NODE
#define COORD_NODE(i, s, ds) (ds.ncoor[(i-s)])

// get size of i-th CELL control volume (distance between two bounding nodes)
#define SIZE_CELL(i, s, ds) (ds.ncoor[(i-s)+1] - ds.ncoor[(i-s)])

// get size of i-th NODE control volume (distance between two neighboring cell centers)
#define SIZE_NODE(i, s, ds) (ds.ccoor[(i-s)] - ds.ccoor[(i-s)-1])

// get interpolation weight for the end of i-th CELL control volume (w_beg = 1 - w_end)
#define WEIGHT_CELL(i, s, ds) ((ds.ccoor[(i-s)] - ds.ncoor[(i-s)])/(ds.ncoor[(i-s)+1] - ds.ncoor[(i-s)]))

// get interpolation weight for the end of i-th NODE control volume (w_beg = 1 - w_end)
#define WEIGHT_NODE(i, s, ds) ((ds.ncoor[(i-s)] - ds.ccoor[(i-s)-1])/(ds.ccoor[(i-s)] - ds.ccoor[(i-s)-1]))

// get interpolation weight for a point in the i-th CELL control volume (local index)
#define WEIGHT_POINT_CELL(i, x, ds) (1.0 - PetscAbsScalar(x - ds.ccoor[i])/(ds.ncoor[i+1] - ds.ncoor[i]))

// get interpolation weight for a point in the i-th NODE control volume (local index)
#define WEIGHT_POINT_NODE(i, x, ds) (1.0 - PetscAbsScalar(x - ds.ncoor[i])/(ds.ccoor[i] - ds.ccoor[i-1]))

// return relative rank of a point
#define GET_POINT_REL_RANK(r, x, bx, ex) { r = 1; if(x < bx) { r--; } else if(x >= ex) { r++; } }

// get consecutive index from I, J, K indices
#define GET_CELL_ID(ID, i, j, k, m, n) { ID = i + (j)*(m) + (k)*(m)*(n); }

// get I, J, K indices from consecutive index
#define GET_CELL_IJK(ID, i, j, k, m, n) \
	(k) = (ID)/((m)*(n));               \
	(j) = (ID - (k)*(m)*(n))/m;         \
	(i) =  ID - (k)*(m)*(n) - (j)*(m);

//---------------------------------------------------------------------------

// initialize standard access loop
#define START_STD_LOOP \
	for(k = sz; k < sz+nz; k++) \
	{	for(j = sy; j < sy+ny; j++) \
		{	for(i = sx; i < sx+nx; i++) \
			{

// finalize standard access loop
#define END_STD_LOOP \
			} \
		} \
	}

//---------------------------------------------------------------------------

// initialize plane access loop
#define START_PLANE_LOOP \
	for(j = sy; j < sy+ny; j++) \
	{	for(i = sx; i < sx+nx; i++) \
		{

// finalize plane access loop
#define END_PLANE_LOOP \
		} \
	}

//---------------------------------------------------------------------------

// scatter operation (two-vectors)
#define GLOBAL_TO_LOCAL(dm, gvec, lvec) \
	ierr = DMGlobalToLocalBegin(dm, gvec, INSERT_VALUES, lvec); CHKERRQ(ierr); \
	ierr = DMGlobalToLocalEnd  (dm, gvec, INSERT_VALUES, lvec); CHKERRQ(ierr);

// scatter operation (one-vector)
#define LOCAL_TO_LOCAL(dm, lvec) \
	ierr = DMLocalToLocalBegin(dm, lvec, INSERT_VALUES, lvec); CHKERRQ(ierr); \
	ierr = DMLocalToLocalEnd  (dm, lvec, INSERT_VALUES, lvec); CHKERRQ(ierr);

// assembly operation
#define LOCAL_TO_GLOBAL(dm, lvec, gvec) \
	ierr = VecZeroEntries(gvec); CHKERRQ(ierr); \
	ierr = DMLocalToGlobalBegin(dm, lvec, ADD_VALUES, gvec); CHKERRQ(ierr); \
	ierr = DMLocalToGlobalEnd  (dm, lvec, ADD_VALUES, gvec); CHKERRQ(ierr);

// create and initialize local vector, scatter ghost values, access array
#define GET_INIT_LOCAL_VECTOR(dm, gvec, lvec, array) \
	ierr = DMGetLocalVector(dm, &lvec); CHKERRQ(ierr); \
	ierr = DMGlobalToLocalBegin(dm, gvec, INSERT_VALUES, lvec); CHKERRQ(ierr); \
	ierr = DMGlobalToLocalEnd  (dm, gvec, INSERT_VALUES, lvec); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(dm, lvec, &array); CHKERRQ(ierr);

// close access & return local vector
#define RESTORE_LOCAL_VECTOR(dm, lvec, array) \
	ierr = DMDAVecRestoreArray(dm, lvec, &array); CHKERRQ(ierr); \
	ierr = DMRestoreLocalVector(dm, &lvec);

//-----------------------------------------------------------------------------
// WRAPPERS
//-----------------------------------------------------------------------------

PetscErrorCode DMDACreate3dSetUp(MPI_Comm comm,
	DMBoundaryType bx, DMBoundaryType by, DMBoundaryType bz, DMDAStencilType stencil_type,
	PetscInt M, PetscInt N, PetscInt P, PetscInt m, PetscInt n, PetscInt p,
	PetscInt dof, PetscInt s, const PetscInt lx[], const PetscInt ly[], const PetscInt lz[], DM *da);


//---------------------------------------------------------------------------
#endif

