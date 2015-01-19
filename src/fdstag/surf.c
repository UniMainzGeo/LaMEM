//---------------------------------------------------------------------------
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "surf.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * stair-case type of free surface
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfCreate"
PetscErrorCode FreeSurfCreate(FDSTAG *fs, UserCtx *user)
{
	// here we must make sure that local part of the free surface
	// COMPLETELY overlaps the local part of computational domain

	PetscScalar z_surface;
	PetscScalar xs, ys, zs;
	PetscScalar dx, dy, dz;
	PetscInt    nx, ny, *lx, *ly;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get mesh sizes
	xs = user->x_left;  dx = user->W;
	ys = user->y_front; dy = user->L;
	zs = user->z_bot;   dz = user->H;

	// get partitioning of the free surface in the XY plane
	ierr = FreeSurfGetPartition(&fs->dsx, xs, dx, fs->dsx.h_min, &nx, &lx); CHKERRQ(ierr);
	ierr = FreeSurfGetPartition(&fs->dsy, ys, dy, fs->dsy.h_min, &ny, &ly); CHKERRQ(ierr);

	// create a DMDA that holds the surface topography
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		user->BC.BCType_x, user->BC.BCType_y, user->BC.BCType_z,
		DMDA_STENCIL_BOX, nx, ny, fs->dsz.nproc,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, NULL, &user->DA_SurfaceTopography); CHKERRQ(ierr);

	// set field name for output
	ierr = DMDASetFieldName(user->DA_SurfaceTopography, 0, "SurfaceTopography"); CHKERRQ(ierr);

	// set uniform grid spacing (refinement is not supported)
	ierr = DMDASetUniformCoordinates(user->DA_SurfaceTopography,
		xs,  xs + dx,
		ys,  ys + dy,
		zs,  zs + dz); CHKERRQ(ierr);

	// create free surface topography vector
	ierr = DMCreateGlobalVector(user->DA_SurfaceTopography,	&user->SurfaceTopography); 	CHKERRQ(ierr);

	// by default the free surface is coincident with the top mesh boundary
	z_surface = user->z_bot	+ user->H;

/*	// otherwise use specified free surface height
	if(user->ErosionParameters.UseInternalFreeSurface == 1)
	{
		z_surface = user->ErosionParameters.InitialFreeSurfaceHeight;
	}*/

	// initialize (internal) surface topography
	ierr = VecSet(user->SurfaceTopography, z_surface); CHKERRQ(ierr);

	// create free surface velocity vectors (Vx, Vy, Vz)
	ierr = VecDuplicate(user->SurfaceTopography, &user->SurfaceTopography_Vx); CHKERRQ(ierr);
	ierr = VecDuplicate(user->SurfaceTopography, &user->SurfaceTopography_Vy); CHKERRQ(ierr);
	ierr = VecDuplicate(user->SurfaceTopography, &user->SurfaceTopography_Vz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfDestroy"
PetscErrorCode FreeSurfDestroy(UserCtx *user)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = DMDestroy(&user->DA_SurfaceTopography);  CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography);    CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vx); CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vy); CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetPartition"
PetscErrorCode FreeSurfGetPartition(
	Discret1D    *ds,  // discretization
	PetscScalar   beg, // starting coordinate
	PetscScalar   len, // domain length
	PetscScalar   h,   // target mesh step
	PetscInt     *n,   // total number of nodes
	PetscInt    **l)   // free surface partitioning vector
{
	MPI_Comm    comm;
	PetscScalar step, tol, rtol;
	PetscInt    i, first, last, lnum, nnod, ncel, sum, *part;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	rtol = 1e-8;

	// compute total number of cells & nodes
	ncel = (PetscInt)ceil(len/h);
	nnod = ncel + 1;

	if(ds->nproc == 1)
	{
		// single processor case
		part = NULL;
	}
	else
	{
		// compute actual mesh step & geometric tolerance
		step = len/(PetscScalar)ncel;
		tol  = rtol*step;

		// get index of first local node
		if(ds->grprev != -1) first = (PetscInt)floor((ds->ncoor[0] - beg + tol)/step);
		else                 first = 0;

		// get index of last local node
		if(ds->grnext != -1) last = (PetscInt)floor((ds->ncoor[ds->ncels] - beg - tol)/step);
		else                 last = ncel;

		// get local number of nodes
		lnum = last - first + 1;

		// create column communicator
		ierr = Discret1DGetColumnComm(ds, &comm); CHKERRQ(ierr);

		// create partitioning vector
		ierr = makeIntArray(&part, NULL, ds->nproc); CHKERRQ(ierr);

		// gather partitioning vector
		ierr = MPI_Allgather(&lnum, 1, MPIU_INT, part, 1, MPIU_INT, comm); CHKERRQ(ierr);

		// free
		ierr = MPI_Comm_free(&comm); CHKERRQ(ierr);

		// checksum
		for(i = 0, sum = 0; i < ds->nproc; i++) sum += part[i];

		if(sum != nnod) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Inconsistent free surface partitioning");
	}

	// return partitioning
	(*n) = nnod;
	(*l) = part;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetVel"
PetscErrorCode FreeSurfGetVel(UserCtx *user)
{
	// project velocities from the grid on the free surface

	//PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAdvect"
PetscErrorCode FreeSurfAdvect(FDSTAG *fs, UserCtx *user)
{
	// advect topography on the free surface mesh

	DM          cda;
	Vec         gc, ltopo, lvx, lvy, lvz;
	PetscInt    i, j, sx, sy, nx, ny, iz;
	DMDACoor3d  ***coors;
	PetscScalar ***topo, ***vx, ***vy, ***vz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get local (ghosted) coordinates of the free surface
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography, &cda);    CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography, &gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, gc, &coors);                       CHKERRQ(ierr);

	// create & access local (ghosted) versions of the topograhy & velocity vectors
	GET_INIT_LOCAL_VECTOR(user->DA_SurfaceTopography, user->SurfaceTopography,    ltopo, topo)
	GET_INIT_LOCAL_VECTOR(user->DA_SurfaceTopography, user->SurfaceTopography_Vx, lvx,   vx)
	GET_INIT_LOCAL_VECTOR(user->DA_SurfaceTopography, user->SurfaceTopography_Vy, lvy,   vy)
	GET_INIT_LOCAL_VECTOR(user->DA_SurfaceTopography, user->SurfaceTopography_Vz, lvz,   vz)









	// get index ranges
	ierr = DMDAGetCorners(user->DA_SurfaceTopography, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	// get vertical index (redundant storage)
	iz = fs->dsz.rank;

	// loop over all local nodes
	for(j = sy; j < sy + ny; j++)
	{	for(i = sx; i < sx + nx; i++)
		{

		}
	}

	// restore access
	ierr = DMDAVecRestoreArray(cda, gc, &coors); CHKERRQ(ierr);

	// return local vectors
	RESTORE_LOCAL_VECTOR(user->DA_SurfaceTopography, ltopo,topo)
	RESTORE_LOCAL_VECTOR(user->DA_SurfaceTopography, lvx,  vx)
	RESTORE_LOCAL_VECTOR(user->DA_SurfaceTopography, lvy,  vy)
	RESTORE_LOCAL_VECTOR(user->DA_SurfaceTopography, lvz,  vz)


	// update free surface coordinates
//	ierr = DMGetCoordinates(user->DA_SurfaceTopography, &global); CHKERRQ(ierr);
//	ierr = DMLocalToGlobalBegin(cda, gc, INSERT_VALUES, global);  CHKERRQ(ierr);
//	ierr = DMLocalToGlobalEnd  (cda, gc, INSERT_VALUES, global);  CHKERRQ(ierr);


/*
	DMDASetUniformCoordinates(da,0.0,1.0,0.0,1.0,0.0,1.0);
	DMGetCoordinateDM(da,&cda);
	DMGetCoordinatesLocal(da,&gc);
	DMDAVecGetArray(cda,gc,&coors);
	DMDAGetCorners(cda,&start,0,0,&m,0,0);
	for (i=start; i<start+m; i++)
	{
		if(i % 2)
		{
			coors[i] = coors[i-1] + .1*(coors[i+1] - coors[i-1]);
		}
	}
	DMDAVecRestoreArray(cda,gc,&coors);
	DMGetCoordinates(da,&global);
	DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global);
	DMLocalToGlobalEnd(cda,gc,INSERT_VALUES,global);
 */



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscInt InterpTriangle(
	PetscScalar *x,  // x-coordinates of triangle
	PetscScalar *y,  // y-coordinates of triangle
	PetscScalar *f,  // interpolated field
	PetscInt    *i,  // indices of triangle corners
	PetscScalar  xp, // x-coordinate of point
	PetscScalar  yp, // y-coordinate of point
	PetscScalar *fp) // field value in the point
{
	PetscScalar xa, xb, xc, ya, yb, yc, la, lb, lc, A;

	// function returns:
	// 1 - if point is inside, fp contains the interpolated value
	// 0 - otherwise

	// access coordinates
	xa = x[i[0]];
	xb = x[i[1]];
	xc = x[i[2]];
	ya = y[i[0]];
	yb = y[i[1]];
	yc = y[i[2]];

	// triangle p-b-c (a-vertex)
	la = GET_AREA_TRIANG(xp, xb, xc, yp, yb, yc);

	// triangle p-c-a (b-vertex)
	lb = GET_AREA_TRIANG(xp, xc, xa, yp, yc, ya);

	// triangle p-a-b (c-vertex)
	lc = GET_AREA_TRIANG(xp, xa, xb, yp, ya, yb);

	// triangle a-b-c (total area)
	A  = GET_AREA_TRIANG(xa, xb, xc, ya, yb, yc);

	// perform point test
	if(la + lb + lc > A*(1.0 + FLT_EPSILON)) return 0;

	// perform interpolation
	la /= A;
	lb /= A;
	lc /= A;

	(*fp) = la*f[i[0]] + lb*f[i[1]] + lc*f[i[2]];

	return 1;
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FreeSurfInterpTopoPoint"
PetscErrorCode FreeSurfInterpTopoPoint(PetscScalar *crd_stencil,  PetscScalar *topo_stencil, PetscScalar *crd_)
{
	// search a point within 2 x 2 cell stencil subdivided into 8 triangles
	// interpolate scalar field

	PetscErrorCode ierr;
	PetscFunctionBegin;


	PetscFunctionReturn(0);
}
*/



/*

//============================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "SetSinusoidalPerturbation"
// There are some cases in which we set an initial sinusoidal perturbation on the free surface
PetscErrorCode SetSinusoidalPerturbation(PetscScalar SinusoidalFreeSurfaceAmplitude, UserCtx *user)
{
	PetscErrorCode ierr;
	PetscScalar    ***LocalSurfaceTopography, x,z, maxVec;
	PetscInt	   xs_Z, ys_Z, zs_Z, xm_Z, ym_Z, zm_Z, ix, iy;
	DM			   cda_SurfaceTopo;
	Vec			   gc_SurfaceTopo;
	DMDACoor3d	   ***coors_SurfaceTopo;
	// nondimensionalize
	SinusoidalFreeSurfaceAmplitude = SinusoidalFreeSurfaceAmplitude/user->Characteristic.Length;
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,	&cda_SurfaceTopo); 	                         CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography, &gc_SurfaceTopo); 	                     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo, &coors_SurfaceTopo); 	                     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,user->SurfaceTopography, &LocalSurfaceTopography); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs_Z,&ys_Z,&zs_Z,&xm_Z,&ym_Z,&zm_Z);               CHKERRQ(ierr);
	for (iy=ys_Z; iy<ys_Z+ym_Z; iy++)
	{	for(ix=xs_Z; ix<xs_Z+xm_Z; ix++)
		{	// Extract x,y,z coordinates of surface topography
			x = coors_SurfaceTopo[zs_Z][iy][ix].x;
//			y = coors_SurfaceTopo[zs_Z][iy][ix].y;
			z = LocalSurfaceTopography[zs_Z][iy][ix];
			z = z + cos(x/user->W*2*PETSC_PI)*SinusoidalFreeSurfaceAmplitude;
			// set topography
			LocalSurfaceTopography[zs_Z][iy][ix] = z;
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);
	VecMax(user->SurfaceTopography,PETSC_NULL, &maxVec);
	PetscPrintf(PETSC_COMM_WORLD,"max topo = %f ", maxVec*user->Characteristic.Length/1000.0);
	PetscFunctionReturn(0);
}


*/
