//---------------------------------------------------------------------------
//.............................. FREE SURFACE ...............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "Utils.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "advect.h"
#include "interpolate.h"
#include "surf.h"
//---------------------------------------------------------------------------
// * stair-case type of free surface
// ...
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfClear"
PetscErrorCode FreeSurfClear(FreeSurf *surf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(surf, sizeof(FreeSurf)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfCreate"
PetscErrorCode FreeSurfCreate(FreeSurf *surf, JacRes *jr)
{
	FDSTAG         *fs;
	const PetscInt *lx, *ly;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = FreeSurfReadFromOptions(surf, &jr->scal); CHKERRQ(ierr);

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	// store context
	surf->jr = jr;

	// access context
	fs = jr->fs;

	// get grid partitioning in X & Y directions
	ierr = DMDAGetOwnershipRanges(fs->DA_COR, &lx, &ly, NULL); CHKERRQ(ierr);

	// create redundant free surface DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_BOX,
		fs->dsx.tnods, fs->dsy.tnods, fs->dsz.nproc,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, NULL, &surf->DA_SURF); CHKERRQ(ierr);

	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->topo); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vx);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vy);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (surf->DA_SURF, &surf->vz);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(surf->DA_SURF, &surf->wa);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(surf->DA_SURF, &surf->wb);   CHKERRQ(ierr);

	// set initial internal free surface level
	ierr = VecSet(surf->topo, surf->InitLevel); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfReadFromOptions"
PetscErrorCode FreeSurfReadFromOptions(FreeSurf *surf, Scaling *scal)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read output flags
	ierr = PetscOptionsGetBool  (NULL, "-surf_use",       &surf->UseFreeSurf,NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, "-surf_level",     &surf->InitLevel,  NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (NULL, "-surf_air_phase", &surf->AirPhase,   NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, "-surf_max_angle", &surf->MaxAngle,   NULL); CHKERRQ(ierr);

	// nondimensionalize
	surf->InitLevel /= scal->length;
	surf->MaxAngle  /= scal->angle;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfDestroy"
PetscErrorCode FreeSurfDestroy(FreeSurf *surf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	ierr = DMDestroy (&surf->DA_SURF); CHKERRQ(ierr);
	ierr = VecDestroy(&surf->topo);    CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vx);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vy);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->vz);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->wa);      CHKERRQ(ierr);
	ierr = VecDestroy(&surf->wb);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfAdvect"
PetscErrorCode FreeSurfAdvect(FreeSurf *surf)
{
	// advect topography on the free surface mesh

	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(surf->UseFreeSurf != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	jr = surf->jr;

	// get surface velocities
	ierr = FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy); CHKERRQ(ierr);
	ierr = FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz); CHKERRQ(ierr);

	// advect topography
	ierr = FreeSurfGetTopo(surf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetVelComp"
PetscErrorCode FreeSurfGetVelComp(
	FreeSurf *surf,
	PetscErrorCode (*interp)(FDSTAG *, Vec, Vec, InterpFlags),
	Vec vcomp_grid, Vec vcomp_surf)
{
	// project velocity component from grid faces on the free surface

	JacRes      *jr;
	FDSTAG      *fs;
	Discret1D   *dsz;
	InterpFlags iflag;
	PetscInt    i, j, nx, ny, sx, sy, sz, level, K;
	PetscScalar ***topo, ***vsurf, ***vgrid, *vpatch, *vmerge, z, w;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr    = surf->jr;
	fs    = jr->fs;
	dsz   = &fs->dsz;
	level = dsz->rank;

	// create column communicator
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

	// set interpolation flags
	iflag.update    = PETSC_FALSE;
	iflag.use_bound = PETSC_TRUE;

	// interpolate velocity component from grid faces to corners
	ierr = interp(fs, vcomp_grid, jr->lbcor, iflag); CHKERRQ(ierr);

	// load ghost values
	LOCAL_TO_LOCAL(fs->DA_COR, jr->lbcor)

	// clear buffer vector for surface velocity patch
	ierr = VecZeroEntries(surf->wa); CHKERRQ(ierr);

	// access topograpy, grid and surface velocity
	ierr = DMDAVecGetArray(fs->DA_COR,    jr->lbcor,  &vgrid); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->wa,   &vsurf); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->topo, &topo);  CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// get topography
		z = topo[level][j][i];

		// check whether point belongs to domain
		if(z >= dsz->crdbeg && z < dsz->crdend)
		{
			// find containing cell
			K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, z);

			// get interpolation weight
			w = (z - dsz->ncoor[K])/(dsz->ncoor[K+1] - dsz->ncoor[K]);

			// interpolate velocity
			vsurf[level][j][i] = (1.0 - w)*vgrid[sz+K][j][i] + w*vgrid[sz+K+1][j][i];
		}
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_COR,    jr->lbcor,  &vgrid); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->wa,   &vsurf); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->topo, &topo);  CHKERRQ(ierr);

	// merge velocity patches
	// compute ghosted version of the velocity component
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(surf->wa, &vpatch); CHKERRQ(ierr);
		ierr = VecGetArray(surf->wb, &vmerge); CHKERRQ(ierr);

		ierr = MPI_Allreduce(vpatch, vmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);

		ierr = VecRestoreArray(surf->wa, &vpatch); CHKERRQ(ierr);
		ierr = VecRestoreArray(surf->wb, &vmerge); CHKERRQ(ierr);

		// compute ghosted version of the velocity component
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->wb, vcomp_surf);
	}
	else
	{
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->wa, vcomp_surf);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FreeSurfGetTopo"
PetscErrorCode FreeSurfGetTopo(FreeSurf *surf)
{
	// advect topography on the free surface mesh

	JacRes      *jr;
	FDSTAG      *fs;
	PetscInt    I, I1, I2, J, J1, J2;
	PetscInt    i, j, jj, found, nx, ny, sx, sy, L, mx, my;
	PetscScalar cx[13], cy[13], cz[13];
	PetscScalar X, X1, X2, Y, Y1, Y2, Z, Exx, Eyy, step;
	PetscScalar ***advect, ***topo, ***vx, ***vy, ***vz;

	// local search grid triangulation
	PetscInt tria [] =
	{
		// first layer
		4, 5,  12, // 0
		4, 12, 7,  // 1
		4, 7,  11, // 2
		4, 11, 3,  // 3
		4, 3,  9,  // 4
		4, 9,  1,  // 5
		4, 1,  10, // 6
		4, 10, 5,  // 7

		// second layer
		5, 8, 12,  // 8
		8, 7, 12,  // 9
		7, 6, 11,  // 10
		6, 3, 11,  // 11
		3, 0, 9,   // 12
		0, 1, 9,   // 13
		1, 2, 10,  // 14
		2, 5, 10   // 15
	};

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	jr   = surf->jr;
	fs   = jr->fs;
	Exx  = jr->bc->Exx;
	Eyy  = jr->bc->Eyy;
	step = jr->ts.dt;
	mx   = fs->dsx.tnods;
	my   = fs->dsy.tnods;
	L    = fs->dsz.rank;

	// access surface topograpy and velocity
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->wa,    &advect); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->topo,  &topo);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx,    &vx);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy,    &vy);     CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz,    &vz);     CHKERRQ(ierr);

	// scan all free surface local points
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// get node coordinates
		X  = COORD_NODE(i,   sx, fs->dsx);
		X1 = COORD_NODE(i-1, sx, fs->dsx);
		X2 = COORD_NODE(i+1, sx, fs->dsx);
		Y  = COORD_NODE(j,   sy, fs->dsy);
		Y1 = COORD_NODE(j-1, sy, fs->dsy);
		Y2 = COORD_NODE(j+1, sy, fs->dsy);

		// get node indices
		I  = i;
		I1 = I-1; if(I1 == -1) I1 = I;
		I2 = I+1; if(I2 == mx) I2 = I;
		J  = j;
		J1 = J-1; if(J1 == -1) J1 = J;
		J2 = J+1; if(J2 == my) J2 = J;

		// compute deformed grid x-coordinates
		cx[0]  = step*vx[L][J1][I1] + X1;
		cx[1]  = step*vx[L][J1][I ] + X;
		cx[2]  = step*vx[L][J1][I2] + X2;
		cx[3]  = step*vx[L][J ][I1] + X1;
		cx[4]  = step*vx[L][J ][I ] + X;
		cx[5]  = step*vx[L][J ][I2] + X2;
		cx[6]  = step*vx[L][J2][I1] + X1;
		cx[7]  = step*vx[L][J2][I ] + X;
		cx[8]  = step*vx[L][J2][I2] + X2;
		cx[9]  = (cx[0] + cx[1] + cx[3] + cx[4])/4.0;
		cx[10] = (cx[1] + cx[2] + cx[4] + cx[5])/4.0;
		cx[11] = (cx[3] + cx[4] + cx[6] + cx[7])/4.0;
		cx[12] = (cx[4] + cx[5] + cx[7] + cx[8])/4.0;

		// compute deformed grid y-coordinates
		cy[0]  = step*vy[L][J1][I1] + Y1;
		cy[1]  = step*vy[L][J1][I ] + Y1;
		cy[2]  = step*vy[L][J1][I2] + Y1;
		cy[3]  = step*vy[L][J ][I1] + Y;
		cy[4]  = step*vy[L][J ][I ] + Y;
		cy[5]  = step*vy[L][J ][I2] + Y;
		cy[6]  = step*vy[L][J2][I1] + Y2;
		cy[7]  = step*vy[L][J2][I ] + Y2;
		cy[8]  = step*vy[L][J2][I2] + Y2;
		cy[9]  = (cy[0] + cy[1] + cy[3] + cy[4])/4.0;
		cy[10] = (cy[1] + cy[2] + cy[4] + cy[5])/4.0;
		cy[11] = (cy[3] + cy[4] + cy[6] + cy[7])/4.0;
		cy[12] = (cy[4] + cy[5] + cy[7] + cy[8])/4.0;

		// compute deformed grid z-coordinates
		cz[0]  = step*vz[L][J1][I1] + topo[L][J1][I1];
		cz[1]  = step*vz[L][J1][I ] + topo[L][J1][I ];
		cz[2]  = step*vz[L][J1][I2] + topo[L][J1][I2];
		cz[3]  = step*vz[L][J ][I1] + topo[L][J ][I1];
		cz[4]  = step*vz[L][J ][I ] + topo[L][J ][I ];
		cz[5]  = step*vz[L][J ][I2] + topo[L][J ][I2];
		cz[6]  = step*vz[L][J2][I1] + topo[L][J2][I1];
		cz[7]  = step*vz[L][J2][I ] + topo[L][J2][I ];
		cz[8]  = step*vz[L][J2][I2] + topo[L][J2][I2];
		cz[9]  = (cz[0] + cz[1] + cz[3] + cz[4])/4.0;
		cz[10] = (cz[1] + cz[2] + cz[4] + cz[5])/4.0;
		cz[11] = (cz[3] + cz[4] + cz[6] + cz[7])/4.0;
		cz[12] = (cz[4] + cz[5] + cz[7] + cz[8])/4.0;

		// compute updated node position if background strain rate is defined
		X *= (1.0 + step*Exx);
		Y *= (1.0 + step*Eyy);

		// find point in the deformed grid, interpolate topography
		found = 0;

		for(jj = 0; jj < 16; jj++)
		{
			found = InterpTriangle(cx, cy, cz, tria + 3*jj, X, Y, &Z);

			if(found) break;
		}

		if(!found)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Problems with topography advection");
		}

		// store advected topography
		advect[L][J][I] = Z;

	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->wa,    &advect); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->topo,  &topo);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx,    &vx);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy,    &vy);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz,    &vz);     CHKERRQ(ierr);

	// compute ghosted version of the surface topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->wa, surf->topo);

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
//---------------------------------------------------------------------------
// get partitioning of the free surface in the XY plane
// ierr = FreeSurfGetPartition(&fs->dsx, xs, dx, fs->dsx.h_min, &nx, &lx); CHKERRQ(ierr);
// ierr = FreeSurfGetPartition(&fs->dsy, ys, dy, fs->dsy.h_min, &ny, &ly); CHKERRQ(ierr);
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

		if(sum != nnod) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Inconsistent free surface partitioning");
	}

	// return partitioning
	(*n) = nnod;
	(*l) = part;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
