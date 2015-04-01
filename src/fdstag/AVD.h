//---------------------------------------------------------------------------
//..........   Routines based on Approximate Voronoi Diagram (AVD)  .........
//---------------------------------------------------------------------------
#ifndef __AVD_h__
#define __AVD_h__
//---------------------------------------------------------------------------

#define AVD_CELL_MASK      -2
#define AVD_CELL_UNCLAIMED -1

typedef struct
{
	PetscInt    ind;   // single cell index
	PetscInt    i,j,k; // i,j,k  cell index
	PetscScalar x[3];  // coordinates of center
	PetscInt    p;     // marker index
	PetscBool   done;  // flag
	PetscInt    col;   // colour for half-centroid

} AVDCell;

typedef struct
{
	PetscInt    p;                      // marker index
	PetscInt    ind;                    // index
	PetscInt    length;                 // current length of chain
	PetscInt    nclaimed;               // claimed cells in the current cycle
	PetscInt    tclaimed;               // total no of cells claimed
	PetscInt    ibound;                 // no. of new boundary cells
	PetscInt    iclaim;                 // no. of new claimed cells
	PetscInt    *bound;                 // new boundary cells
	PetscInt    *claim;                 // new claimed cells
	PetscBool   done;                   // flag
	PetscInt    gind;                   // marker index in actx->markers
	PetscScalar xc[3];                  // centroid coordinates
	PetscScalar xh, yh, zh;             // half-axis of the centroid

} AVDChain;

typedef struct
{
	PetscInt    mmin, mmax;       // limit number of markers
	PetscScalar xs[3],xe[3];      // coordinate limits of the Voronoi diagram
	PetscScalar dx,dy,dz;         // spacing
	PetscInt    nx,ny,nz;         // grid cells
	PetscInt    buffer;           // buffer
	AVDCell     *cell;            // voronoi grid (size of nx*ny*nz)
	AVDChain    *chain;           // voronoi chain for every point (size of npoints)
	Marker      *points;          // points that we want to compute voronoi diagram (size of npoints)
	PetscInt    npoints;          // no. markers

} AVD3D;

//---------------------------------------------------------------------------
// basic AVD routines
PetscErrorCode AVDCreate     (AVD3D *A);
PetscErrorCode AVDDestroy    (AVD3D *A);
PetscErrorCode AVDCellInit   (AVD3D *A);
PetscErrorCode AVDClaimCells (AVD3D *A, const PetscInt ip);
PetscErrorCode AVDUpdateChain(AVD3D *A, const PetscInt ip);
PetscErrorCode AVDReAlloc    (AVDChain *chain,PetscInt buffer);

// routines for marker control
PetscErrorCode AVDLoadPoints            (AdvCtx *actx, AVD3D *A, PetscInt ind);
PetscErrorCode AVDInjectDeletePoints    (AdvCtx *actx, AVD3D *A);
PetscErrorCode AVDExecuteMarkerInjection(AdvCtx *actx, PetscInt npoints, PetscScalar xs[3], PetscScalar xe[3], PetscInt ind);
//---------------------------------------------------------------------------
static inline PetscScalar AVDDistanceTest(PetscScalar x0[3],PetscScalar x1[3],PetscScalar x2[3])
{
	return (x1[0]+x2[0]-2*x0[0])*(x1[0]-x2[0]) + (x1[1]+x2[1]-2*x0[1])*(x1[1]-x2[1]) + (x1[2]+x2[2]-2*x0[2])*(x1[2]-x2[2]);
}
//---------------------------------------------------------------------------
#endif
