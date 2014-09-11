//---------------------------------------------------------------------------
//...................   MATERIAL ADVECTION ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __advect_h__
#define __advect_h__
//---------------------------------------------------------------------------

#define _cap_overhead_ 1.3

//---------------------------------------------------------------------------

// Set of variables that must be tracked during the advection steps.
// During the nonlinear iteration, history accumulates on the integration points.
// During the advection a cumulative update of the history variables is added to
// the individual histories of the markers locating in the integration point
// control volume. There are four types of control volumes for different stress
// components: XY-, XZ-, YZ- edges, and cell centers. Each particle should be
// mapped on all types of control volumes to update corresponding components.

//---------------------------------------------------------------------------
// Advection context
typedef struct
{

	//=============
	// COMMUNICATOR
	//=============
	MPI_Comm  icomm;   // distinct communicator for communicating markers
	PetscInt  nproc;   // total number of processors
	PetscInt  iproc;   // processor rank

	//========
	// STORAGE
	//========
	PetscInt  nummark; // local number of markers
	PetscInt  markcap; // capacity of marker storage
	Marker 	 *markers; // storage for local markers
	PetscInt *cellnum; // host cells local numbers

	//=========
	// EXCHANGE
	//=========
	Marker   *sendbuf; // send buffer
	Marker   *recvbuf; // receive buffer

	PetscInt  nsend;                // total number of markers to be sent (local)
	PetscInt  nsendm[_num_neighb_]; // number of markers to be sent to each process
	PetscInt  ptsend[_num_neighb_]; // send buffer pointers

	PetscInt  nrecv;                // total number of markers to be received (local)
	PetscInt  nrecvm[_num_neighb_]; // number of markers to be received from each process
	PetscInt  ptrecv[_num_neighb_]; // receive buffer pointers

	PetscInt  ndel; // number of markers to be deleted from storage
	PetscInt *idel;	// indices of markers to be deleted

	// Entscheidung treffen:

	// 1. Viscosities are computed in the centers & then averaged to edges (BY FAR THE SIMPLEST SOLUTION!!!)
	// 2. Synchronize SEPARATELY every phase for a given edge set xy, xz, or yz (overwhelming communication)
	// 3. Synchronize SIMULTANEOUSLY all the phases for a given edge set xy, xz, or yz (large memory requirements)
	// 4. Duplicate the markers in the overlapping control volumes near the inter-processor boundaries (a compromise, but still more memory)

} AdvCtx;

//---------------------------------------------------------------------------

PetscErrorCode ADVCreate(AdvCtx *actx);

PetscErrorCode ADVDestroy(AdvCtx *actx);

PetscErrorCode ADVReAllocateStorage(AdvCtx *actx, PetscInt capacity);

PetscErrorCode ADVAdvect(AdvCtx *actx, FDSTAG *fs);

PetscErrorCode ADVAdvectMarkers(AdvCtx *actx, FDSTAG *fs, JacResCtx *jrctx);

// project history variables from markers to grid
PetscErrorCode ADVProjHistMarkGrid(AdvCtx *actx, FDSTAG *fs, JacResCtx *jrctx);

// count number of markers to be sent to each neighbor domain
PetscErrorCode ADVMapMarkersDomains(AdvCtx *actx, FDSTAG *fs);

// communicate number of markers with neighbor processes
PetscErrorCode ADVExchangeNumMarkers(AdvCtx *actx, FDSTAG *fs);

// create send and receive buffers for asynchronous MPI communication
PetscErrorCode ADVCreateMPIBuffer(AdvCtx *actx, FDSTAG *fs);

// communicate markers with neighbor processes
PetscErrorCode ADVExchangeMarkers(AdvCtx *actx, FDSTAG *fs);

// free communication buffer
PetscErrorCode ADVDestroyMPIBuffer(AdvCtx *actx);

// store received markers, collect garbage
PetscErrorCode ADVCollectGarbage(AdvCtx *actx);

// find host cells for local markers
PetscErrorCode ADVMapMarkersCells(AdvCtx *actx, FDSTAG *fs);

PetscErrorCode FDSTAGetVorticity(
	FDSTAG *fs,
	Vec lvx,  Vec lvy,  Vec lvz,  // local (ghosted) velocities
	Vec gwx,  Vec gwy,  Vec gwz); // global vorticity components

//-----------------------------------------------------------------------------
// call this function for local markers only!
static inline PetscInt FindPointInCell(
	PetscScalar *px, // node coordinates
	PetscInt     L,  // index of the leftmost node
	PetscInt     R,  // index of the rightmost node
	PetscScalar  x)  // point coordinate
{
	// get initial guess assuming uniform grid
	PetscInt M = L + (PetscInt)((x-px[L])/((px[R]-px[L])/(PetscScalar)(R-L)));

	if(M <= L) return L;
	if(M >= R) return R-1;

	if(px[M]   <= x) L=M;
	if(px[M+1] >= x) R=M+1;

	while((R-L) > 1)
	{
		M = (L+R)/2;
		if(px[M] <= x) L=M;
		if(px[M] >= x) R=M;

	}
	return(L);
}
//-----------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------

// compute pointers from counts, return total count
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]);

// rewind pointers after using them as access iterators
void rewindPtr(PetscInt n, PetscInt ptr[]);

// normalize vector by the inverse sum of its elements
static inline PetscScalar normVect(PetscInt n, PetscScalar *v)
{
	// normalize vector by the inverse sum of its elements

	PetscInt    i;
	PetscScalar sum = 0.0;

	for(i = 0; i < n; i++) sum  += v[i];
	for(i = 0; i < n; i++) v[i] /= sum;

	return sum;
}
//-----------------------------------------------------------------------------

/*

// vorticity components
Vec gwx,  gwy,  gwz; // global vorticity components
*/


//---------------------------------------------------------------------------
#endif
