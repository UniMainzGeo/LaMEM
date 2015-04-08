//---------------------------------------------------------------------------
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#ifndef __surf_h__
#define __surf_h__
//---------------------------------------------------------------------------

// free surface grid

typedef struct
{
	JacRes *jr;               // global residual context
	DM      DA_SURF;          // free surface grid
	Vec     topo, vx, vy, vz; // topography and velocity vectors (local)
	Vec     wa, wb;           // work vectors                    (global)

	// flags/parameters
	PetscBool   UseFreeSurf;
	PetscScalar InitLevel;
	PetscInt    AirPhase;
	PetscScalar MaxAngle;

} FreeSurf;

//---------------------------------------------------------------------------

PetscErrorCode FreeSurfClear(FreeSurf *surf);

PetscErrorCode FreeSurfCreate(FreeSurf *surf, JacRes *jr);

PetscErrorCode FreeSurfReadFromOptions(FreeSurf *surf);

PetscErrorCode FreeSurfDestroy(FreeSurf *surf);

// advect topography on the free surface mesh
PetscErrorCode FreeSurfAdvect(FreeSurf *surf);

// get single velocity component on the free surface
PetscErrorCode FreeSurfGetVelComp(
	FreeSurf *surf,
	PetscErrorCode (*interp)(FDSTAG *, Vec, Vec, InterpFlags),
	Vec vcomp_grid, Vec vcomp_surf);

PetscErrorCode FreeSurfGetTopo(FreeSurf *surf);

PetscInt InterpTriangle(
	PetscScalar *x,   // x-coordinates of triangle
	PetscScalar *y,   // y-coordinates of triangle
	PetscScalar *f,   // interpolated field
	PetscInt    *i,   // indices of triangle corners
	PetscScalar  xp,  // x-coordinate of point
	PetscScalar  yp,  // y-coordinate of point
	PetscScalar *fp); // field value in the point

//---------------------------------------------------------------------------

#define GET_AREA_TRIANG(x1, x2, x3, y1, y2, y3) PetscAbsScalar((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))

//---------------------------------------------------------------------------
#endif

/*
// map uniform free surface onto non-uniform computational grid
PetscErrorCode FreeSurfGetPartition(
	Discret1D    *ds,  // discretization
	PetscScalar   beg, // starting coordinate
	PetscScalar   len, // domain length
	PetscScalar   h,   // target mesh step
	PetscInt     *n,   // total number of nodes
	PetscInt    **l);  // free surface partitioning vector
*/
