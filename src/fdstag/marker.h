//---------------------------------------------------------------------------
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// input polygon data
typedef struct
{
	PetscInt     dir; // normal vector of polygon plane
	PetscInt   ax[2]; // axis that span the polygon plane
	PetscInt   phase; // phase that the polygon defines
	PetscInt    type; // type can be of additive or assigning nature
	PetscInt     num; // number of polygon slices defining the volume
	PetscInt     len; // number of nodes of polygon
	PetscInt    idxs; // index of first polygon slice
	PetscInt    gidx; // global plane index (consistent with planes of markers)
	PetscInt    lidx; // local plane index (consistent with planes of markers)
	PetscInt       n; // number of polygon nodes
	PetscInt   nmark; // number of markers in current volume
	PetscScalar   *X; // coordinates of the polygon x1,y1,x2,y2,...xn,yn

} Polygon2D;

// internal polygon context to avoid allocations inside "polyin"
typedef struct
{
	PetscScalar   *Xtemp;
	PetscScalar   *nodetemp;
	PetscScalar   *cnect;
	PetscBool     *cn;
	PetscBool     *on;
	PetscScalar   *x;
	PetscScalar   *y;
	PetscInt      *idx;
} PolyCtx;



//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, UserCtx *user);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, UserCtx *user);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx, UserCtx *user);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFilePolygons (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBands        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDomes        (AdvCtx *actx, UserCtx *user);

//---------------------------------------------------------------------------

// service functions

void ADVMarkSecIdx(AdvCtx *actx, UserCtx *user, PetscInt dir, PetscInt Nslice, PetscInt *idx);
void inpoly(PolyCtx *polydat, PetscInt N, PetscScalar *X, PetscScalar *node, PetscInt Nnode, PetscBool *in, PetscBool *bnd);

PetscErrorCode CreatePolyCtx(PolyCtx *polydat, PetscInt N, PetscInt Nnode);
PetscErrorCode DestroyPolyCtx(PolyCtx polydat);
void qsindex (PetscScalar  *a, PetscInt *idx , PetscInt lo, PetscInt hi);

//---------------------------------------------------------------------------

// definitions

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

//---------------------------------------------------------------------------
#endif
