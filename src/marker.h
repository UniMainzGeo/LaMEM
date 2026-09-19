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
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__

//---------------------------------------------------------------------------

struct FB;
struct AdvCtx;
struct Marker;
struct Material_t;

//---------------------------------------------------------------------------

// input volume data
typedef struct
{
	PetscInt     dir; // normal vector of polygon plane
	PetscInt   ax[2]; // axis that span the polygon plane
	PetscInt   phase; // phase that the polygon defines
	PetscInt    type; // type can be of additive or assigning nature
	PetscInt     num; // number of polygon slices defining the volume
} Volume3D;

// input polygon data
typedef struct
{
	PetscInt     len; // number of nodes of polygon
	PetscInt    idxs; // index of first polygon slice
	PetscInt    gidx; // global plane index (consistent with planes of markers)
	PetscInt    lidx; // local plane index (consistent with planes of markers)
	PetscInt       n; // number of polygon nodes
} Polygon2D;

//---------------------------------------------------------------------------

// geometric primitives

// geometric primitive type.
// Stored in GeomPrim.type so that the setPhase function pointer can be restored
// after a binary restart (function addresses are not valid across runs).
enum GeomPrimType
{
	_GEOM_NONE_ = 0,
	_GEOM_SPHERE_,
	_GEOM_ELLIPSOID_,
	_GEOM_BOX_,
	_GEOM_RIDGE_,
	_GEOM_LAYER_,
	_GEOM_HEX_,
	_GEOM_CYLINDER_

};

typedef struct GeomPrim GeomPrim;

struct GeomPrim
{
	PetscInt    phase;
	// sphere & cylinder & ellipsoid
	PetscScalar center[3];
	// sphere & cylinder
	PetscScalar radius;
	// ellipsoid
	PetscScalar axes[3];
	// cylinder
	PetscScalar base[3], cap[3];
	// box & hex
	PetscScalar bounds[6], coord[24];
	// layer
	PetscScalar top;
	PetscScalar bot;
	PetscInt    cosine;
	PetscScalar amplitude;
	PetscScalar wavelength;
	PetscScalar rand_amplitude;
	// ridge
	PetscScalar v_spread;
	PetscScalar x_oblique;
	PetscScalar ridgeseg_x[2];
	PetscScalar ridgeseg_y[2];
	PetscScalar x_ridgeLeft;
	PetscScalar x_ridgeRight;
	PetscScalar y_ridgeFront;
	PetscScalar y_ridgeBack;
	PetscScalar thermalAgeRidge;
	PetscScalar age0;                   // thermal age @ ridge
	PetscScalar maxAge;                 // maximum thermal Age a plate can have [say 80 Myrs on Earth]
	// temperature
	PetscInt    setTemp;
	PetscScalar cstTemp;
	PetscScalar topTemp, botTemp;
	PetscScalar thermalAge;
	PetscScalar kappa;

	// mid-run injection (optional, controlled by the t_inject parameter)
	PetscInt    type;                      // GeomPrimType, used to restore setPhase after restart
	PetscInt    numInject;                 // number of entries in t_inject (default 1, t_inject[0] = 0)
	PetscScalar t_inject[_max_inj_times_]; // injection times (non-dimensional); 0 = apply at initialization
	PetscInt    done[_max_inj_times_];     // 1 once the corresponding t_inject entry has been applied

	void (*setPhase)(GeomPrim*, Marker*);
};

// set primitive type together with the matching setPhase function pointer
void GeomPrimSetType(GeomPrim *geom, PetscInt type);

// does the primitive have to be applied during initialization (t_inject == 0)?
PetscInt GeomPrimAtInit(GeomPrim *geom);

// does the primitive have to be injected during the simulation (t_inject > 0)?
PetscInt GeomPrimDeferred(GeomPrim *geom);

// name of the primitive type (for diagnostic output)
const char * GeomPrimGetName(PetscInt type);

void setPhaseSphere(GeomPrim *sphere, Marker *P);

void setPhaseEllipsoid(GeomPrim *ellipsoid, Marker *P);

void setPhaseBox(GeomPrim *box, Marker *P);

void setPhaseRidge(GeomPrim *ridge, Marker *P);

void setPhaseLayer(GeomPrim *layer, Marker *P);

void setPhaseHex(GeomPrim *hex, Marker *P);

void setPhaseCylinder(GeomPrim *cylinder, Marker *P);

void HexGetBoundingBox(
    PetscScalar *coord,   // hex coordinates
    PetscScalar *bounds); // bounding box

PetscInt TetPointTest(
    PetscScalar *coord, // tetrahedron coordinates
    PetscInt    *ii,    // corner indices
    PetscScalar *xp,    // point coordinate
    PetscScalar  tol);  // relative tolerance

void computeTemperature(GeomPrim *geom, Marker *P, PetscScalar *T );

//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, FB *fb);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx);

// perturb marker coordinates after initialization
PetscErrorCode ADVMarkPerturb(AdvCtx *actx);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx);

// initialize temperature on markers based on linear gradient
PetscErrorCode ADVMarkSetTempGrad(AdvCtx *actx);

// initialize temperature on markers based on phase temperature
PetscErrorCode ADVMarkSetTempPhase(AdvCtx *actx);

// initialize temperature on markers redundantly form file
PetscErrorCode ADVMarkSetTempFile(AdvCtx *actx, FB *fb);

// initialize temperature on markers from vector
PetscErrorCode ADVMarkSetTempVector(AdvCtx *actx);

// Load and set data from phase diagram
PetscErrorCode LoadPhaseDiagram(AdvCtx *actx, Material_t  *phases, PetscInt i);

// read control polygons
struct CtrlP
{
	PetscInt    ID[_max_ctrl_poly_];
	PetscInt    VolID[_max_ctrl_poly_];
	PetscInt    Pos[_max_ctrl_poly_];
	PetscScalar Sx[_max_ctrl_poly_];
	PetscScalar Sy[_max_ctrl_poly_];
};

PetscErrorCode ADVMarkReadCtrlPoly(FB *fb, CtrlP *CtrlPoly, PetscInt &VolID, PetscInt &nCP);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitGeom    (AdvCtx *actx, FB *fb);
PetscErrorCode ADVMarkInitFiles   (AdvCtx *actx, FB *fb);
PetscErrorCode ADVMarkInitPolygons(AdvCtx *actx, FB *fb);

//---------------------------------------------------------------------------

// Mid-run injection of geometric primitives (t_inject)

// read all geometric primitive blocks from the input file
PetscErrorCode ADVMarkReadGeom(AdvCtx *actx, FB *fb, GeomPrim *geom, GeomPrim **pgeom, PetscInt *ngeom_);

// store the primitives that carry a t_inject > 0 entry in the advection context
PetscErrorCode ADVMarkStoreInjectGeom(AdvCtx *actx, GeomPrim **pgeom, PetscInt ngeom);

// read deferred primitives for setups that do not use geometric primitives for the
// initial geometry (msetup = files, e.g. GeophysicalModelGenerator, or msetup = polygons)
PetscErrorCode ADVMarkInitInjectGeom(AdvCtx *actx, FB *fb);

// apply due injections; called every time step from LaMEMLibSolve
PetscErrorCode ADVMarkInjectGeom(AdvCtx *actx);

//---------------------------------------------------------------------------

// service functions
void ADVMarkSecIdx(AdvCtx *actx, PetscInt dir, PetscInt Nslice, PetscInt *idx);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#ifndef max
#define max(a,b) (a >= b ? a : b)
#define min(a,b) (a <= b ? a : b)
#endif

#define GET_GEOM(p, s, i, n) if(i < n) { p = &s[i++]; } \
    else { SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many geometric primitives! Max allowed: %" PetscInt_FMT "", n); }

//---------------------------------------------------------------------------
#endif
