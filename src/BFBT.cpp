#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matrix.h"
#include "surf.h"

#include "BFBT.h"

//---------------------------------------------------------------------------

#define SCATTER_FIELD(da, vec, FIELD) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec)

#define GET_VISC_TOTAL buff[k][j][i] = jr->svCell[iter++].svDev.eta;


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetViscMat"
PetscErrorCode JacResGetViscMat(PMat pm)
{
	// assemble viscosity scaling matrices K and C

	JacRes 		*jr;
	PMatBlock   *P;

	FDSTAG     *fs;
	BCCtx      *bc;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bvx, fvx, bvy, fvy, bvz, fvz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;

 	PetscScalar v[7], cf[6];
	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcvx, ***bcvy, ***bcvz, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	jr 	 = pm->jr;
	fs   = jr->fs;
	bc   = jr->bc;
	num  = bc->pNumSPC;
	list = bc->pSPCList;
	P 	 = (PMatBlock*) pm->data;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_VISC_TOTAL);

	// clear matrix coefficients
	ierr = MatZeroEntries(P->K); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X, bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, bc->bcvz,  &bcvz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// check index bounds and TPC multipliers
		//
		//		-----y_i,j+1-----
		//		|				|
		//		x_i,j	i,j		x_i+1,j
		//		|				|
		//		-----y_i,j-------
		//
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcvx[k][j][i-1] != DBL_MAX) cf[0] = 0.0;}
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcvx[k][j][i+1] != DBL_MAX) cf[1] = 0.0;}
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcvy[k][j-1][i] != DBL_MAX) cf[2] = 0.0;}
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcvy[k][j+1][i] != DBL_MAX) cf[3] = 0.0;}
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcvz[k-1][j][i] != DBL_MAX) cf[4] = 0.0;}
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcvz[k+1][j][i] != DBL_MAX) cf[5] = 0.0;}

		// compute average viscosities
		bvx = (lk[k][j][i] + lk[k][j][Im1])/2.0;      fvx = (lk[k][j][i] + lk[k][j][Ip1])/2.0;
		bvy = (lk[k][j][i] + lk[k][Jm1][i])/2.0;      fvy = (lk[k][j][i] + lk[k][Jp1][i])/2.0;
		bvz = (lk[k][j][i] + lk[Km1][j][i])/2.0;      fvz = (lk[k][j][i] + lk[Kp1][j][i])/2.0;

		// inverse viscosities:
		bvx = 1/sqrt(bvx);  		fvx = 1/sqrt(fvx);
		bvy = 1/sqrt(bvy);  		fvy = 1/sqrt(fvy);
		bvz = 1/sqrt(bvz);  		fvz = 1/sqrt(fvz);

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC multipliers
		v[0] = -bvx/bdx/dx*cf[0];
		v[1] = -fvx/fdx/dx*cf[1];
		v[2] = -bvy/bdy/dy*cf[2];
		v[3] = -fvy/fdy/dy*cf[3];
		v[4] = -bvz/bdz/dz*cf[4];
		v[5] = -fvz/fdz/dz*cf[5];
		v[6] =  (bvx/bdx + fvx/fdx)/dx
			 +  (bvy/bdy + fvy/fdy)/dy
		     +  (bvz/bdz + fvz/fdz)/dz;

		// set matrix coefficients
		ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz);  CHKERRQ(ierr);

	// assemble K matrix
	ierr = MatAIJAssemble(P->K, num, list, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CopyViscosityToScalingVector"
PetscErrorCode CopyViscosityToScalingVector(Vec a, Vec b, Vec c, Vec ScalingVec)
{
	PetscInt     as,  bs,  cs;
	PetscScalar *ap, *bp, *cp, *sv;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecGetLocalSize(a, &as); CHKERRQ(ierr);
	ierr = VecGetLocalSize(b, &bs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(c, &cs); CHKERRQ(ierr);

	// access vectors
	ierr = VecGetArray(a, &ap); CHKERRQ(ierr);
	ierr = VecGetArray(b, &bp); CHKERRQ(ierr);
	ierr = VecGetArray(c, &cp); CHKERRQ(ierr);
	ierr = VecGetArray(ScalingVec, &sv); CHKERRQ(ierr);

	ierr = PetscMemcpy(sv,    	 ap, (size_t)as*sizeof(PetscScalar)); CHKERRQ(ierr);
	ierr = PetscMemcpy(sv+as, 	 bp, (size_t)bs*sizeof(PetscScalar)); CHKERRQ(ierr);
	ierr = PetscMemcpy(sv+as+bs, cp, (size_t)cs*sizeof(PetscScalar)); CHKERRQ(ierr);

	ierr = VecRestoreArray(a, &ap); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &bp); CHKERRQ(ierr);
	ierr = VecRestoreArray(c, &cp); CHKERRQ(ierr);
	ierr = VecRestoreArray(ScalingVec, &sv); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


