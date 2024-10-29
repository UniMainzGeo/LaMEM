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
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#ifndef __surf_h__
#define __surf_h__

//---------------------------------------------------------------------------

struct FB;
struct InterpFlags;
struct FDSTAG;
struct JacRes;

//---------------------------------------------------------------------------

// free surface grid

struct FreeSurf
{
	JacRes *jr;             // global residual context
	DM      DA_SURF;        // free surface grid
	Vec     ltopo, gtopo;   // topography vectors                (local and global)
	Vec     vx, vy, vz;     // velocity vectors                  (local)
	Vec     vpatch, vmerge; // patch and merged velocity vectors (global)

	Vec vz_fs;
	Vec gtopo_fs;
	VecScatter ctx;

	// flags/parameters
	PetscInt    UseFreeSurf; // free surface activation flag
	PetscInt    phaseCorr;   // free surface phase correction flag
	PetscScalar InitLevel;   // initial level
	PetscInt    AirPhase;    // air phase number
	PetscScalar MaxAngle;    // maximum angle with horizon (smoothed if larger)

	// erosion/sedimentation parameters
	PetscInt    ErosionModel;               // [0-none, 1-infinitely fast, 2-prescribed rate...]
	PetscInt    SedimentModel;              // [0-none, 1-prescribed rate, 2-gaussian margin...]
	PetscInt    numLayers;                  // number of sediment layers
	PetscInt    numErPhs;                   // number of erosion phases
	PetscScalar timeDelims[_max_sed_layers_-1];  // sediment layers time delimiters
	PetscScalar timeDelimsEr[_max_er_phases_-1]; // sediment layers time delimiters
	PetscScalar erRates[_max_er_phases_];        // erosion rates
	PetscScalar erLevels[_max_er_phases_];       // erosion levels
	PetscScalar sedRates[_max_sed_layers_  ];    // sedimentation rates
	PetscScalar sedLevels[_max_sed_layers_];     // sedimentation levels
	PetscScalar sedRates2nd[_max_sed_layers_  ]; // sedimentation rates
	PetscInt    sedPhases[_max_sed_layers_  ];   // sediment layers phase numbers
	PetscScalar marginO[2];                 // lateral coordinates of continental margin - origin
	PetscScalar marginE[2];                 // lateral coordinates of continental margin - 2nd point
	PetscScalar hUp;                        // up dip thickness of sediment cover
	PetscScalar hDown;                      // down dip thickness of sediment cover
	PetscScalar dTrans;                     // half of transition zone

	// run-time parameters
	PetscScalar avg_topo; // average topography (updated by all functions changing topography)
	PetscInt    phase;    // current sediment phase

};

//---------------------------------------------------------------------------

PetscErrorCode FreeSurfCreate(FreeSurf *surf, FB *fb);

PetscErrorCode FreeSurfCreateData(FreeSurf *surf);

PetscErrorCode FreeSurfGetAvgTopo(FreeSurf *surf);

PetscErrorCode FreeSurfReadRestart(FreeSurf *surf, FILE *fp);

PetscErrorCode FreeSurfWriteRestart(FreeSurf *surf, FILE *fp);

PetscErrorCode FreeSurfDestroy(FreeSurf *surf);

// advect topography on the free surface mesh
PetscErrorCode FreeSurfAdvect(FreeSurf *surf);

// get single velocity component on the free surface
PetscErrorCode FreeSurfGetVelComp(
	FreeSurf *surf,
	PetscErrorCode (*interp) (FDSTAG *, Vec, Vec, InterpFlags),
	Vec vcomp_grid, Vec vcomp_surf);

// advect/interpolate topography of the free surface
PetscErrorCode FreeSurfAdvectTopo(FreeSurf *surf);

// smooth topography if maximum angle with horizon is exceeded
PetscErrorCode FreeSurfSmoothMaxAngle(FreeSurf *surf);

// correct phase ratios based on actual position of the free surface
PetscErrorCode FreeSurfGetAirPhaseRatio(FreeSurf *surf);

// apply erosion to the free surface
PetscErrorCode FreeSurfAppErosion(FreeSurf *surf);

// apply sedimentation to the free surface
PetscErrorCode FreeSurfAppSedimentation(FreeSurf *surf);

// Set topography from file
PetscErrorCode FreeSurfSetTopoFromFile(FreeSurf *surf, FB *fb);

// Set initial perturbation
PetscErrorCode FreeSurfSetInitialPerturbation(FreeSurf *surf);

//---------------------------------------------------------------------------
// SERVICE FUNCTIONS
//---------------------------------------------------------------------------

PetscInt InterpolateTriangle(
	PetscScalar *x,   // x-coordinates of triangle
	PetscScalar *y,   // y-coordinates of triangle
	PetscScalar *f,   // interpolated field
	PetscInt    *i,   // indices of triangle corners
	PetscScalar  xp,  // x-coordinate of point
	PetscScalar  yp,  // y-coordinate of point
	PetscScalar  tol, // relative tolerance
	PetscScalar *fp); // field value in the point

PetscScalar IntersectTriangularPrism(
	PetscScalar *x,     // x-coordinates of prism base
	PetscScalar *y,     // y-coordinates of prism base
	PetscScalar *z,     // z-coordinates of prism top surface
	PetscInt    *i,     // indices of base corners
	PetscScalar  vcell, // total volume of cell
	PetscScalar  bot,   // z-coordinate of bottom plane
	PetscScalar  top,   // z-coordinate of top plane
	PetscScalar  tol);  // relative tolerance

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// NOTE! this macro computes double of actual area
#define GET_AREA_TRIANG(x1, x2, x3, y1, y2, y3) PetscAbsScalar((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))

// NOTE! this macro computes double of actual volume
#define GET_VOLUME_PRISM(x1, x2, x3, y1, y2, y3, z1, z2, z3, level) \
	((z1+z2+z3)/3.0 > level ? ((z1+z2+z3)/3.0-level)*PetscAbsScalar((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3)) : 0)

#define INTERSECT_EDGE(x1, y1, z1, x2, y2, z2, xp, yp, zp, level, dh) \
	zp = level; \
	w  = z1; if(z2 < w) w = z2; if(zp < w) zp = w; \
	w  = z1; if(z2 > w) w = z2; if(zp > w) zp = w; \
	w  = 0.0; \
	if(PetscAbsScalar(z2-z1) > dh) w = (zp-z1)/(z2-z1); \
	xp = x1 + w*(x2-x1); \
	yp = y1 + w*(y2-y1);

//---------------------------------------------------------------------------
#endif
