//---------------------------------------------------------------------------
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#ifndef __surf_h__
#define __surf_h__
//---------------------------------------------------------------------------

PetscErrorCode FreeSurfCreate(FDSTAG *fs, UserContext *user);

PetscErrorCode FreeSurfDestroy(UserContext *user);

// map uniform free surface onto non-uniform computational grid
PetscErrorCode FreeSurfGetPartition(
	Discret1D    *ds,  // discretization
	PetscScalar   beg, // starting coordinate
	PetscScalar   len, // domain length
	PetscScalar   h,   // target mesh step
	PetscInt     *n,   // total number of nodes
	PetscInt    **l);  // free surface partitioning vector

// project velocities from the grid on the free surface
PetscErrorCode FreeSurfGetVel(UserContext *user);

// advect topography on the free surface mesh
PetscErrorCode FreeSurfAdvect(FDSTAG *fs, UserContext *user);

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
