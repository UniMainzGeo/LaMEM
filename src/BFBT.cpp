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
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetViscRes"
PetscErrorCode JacResGetViscRes(JacRes *jr, PetscScalar dt)
{
	// compute viscosity residual vector
	// STEADY STATE solution is activated by setting time step to zero

	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
 	PetscScalar bvx, fvx, bvy, fvy, bvz, fvz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bdpdx, bdpdy, bdpdz, fdpdx, fdpdy, fdpdz;
 	PetscScalar dx, dy, dz;
	PetscScalar invdt, Vc, Pc, Hr, Ha;
	PetscScalar ***vr, ***lV, ***lk, ***hxy, ***hxz, ***hyz, ***buff, *e,***P;;
	PetscScalar ***vx,***vy,***vz;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs    = jr->fs;
	bc    = jr->bc;
	num   = bc->tNumSPC;
	list  = bc->tSPCList;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	// compute inverse time step
	if(dt) invdt = 1.0/dt;
	else   invdt = 0.0;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_VISC_TOTAL)

	// access work vectors
	ierr = DMDAVecGetArray(jr->DA_V,   jr->vr,   &vr);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lV,   &lV);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx) ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy) ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz) ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &P );  CHKERRQ(ierr);



	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// access
		Vc  = lV[k][j][i]; // current viscosity
		Pc  = P[k][j][i] ; // Current Pressure

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average viscosities
		bvx = (lk[k][j][i] + lk[k][j][Im1])/2.0;      fvx = (lk[k][j][i] + lk[k][j][Ip1])/2.0;
		bvy = (lk[k][j][i] + lk[k][Jm1][i])/2.0;      fvy = (lk[k][j][i] + lk[k][Jp1][i])/2.0;
		bvz = (lk[k][j][i] + lk[Km1][j][i])/2.0;      fvz = (lk[k][j][i] + lk[Kp1][j][i])/2.0;

		// inverse viscosity:
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

		// original balance equation:

		// rho*Cp*(Tc - Tn)/dt = (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz + Hr + rho*A

		// to get positive diagonal in the preconditioner matrix
		// put right hand side to the left, which gives the following:

		vr[k][j][i] = - (fvx - bvx)/dx - (fvy - bvy)/dy - (fvz - bvz)/dz;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_V,   jr->vr,   &vr);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lV,   &lV);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   jr->lvx,     &vx) ;  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,     &vy) ;  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,     &vz) ;  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &P)  ;  CHKERRQ(ierr);


	// impose primary temperature constraints
	ierr = VecGetArray(jr->ge, &e); CHKERRQ(ierr);

	for(i = 0; i < num; i++) e[list[i]] = 0.0;

	ierr = VecRestoreArray(jr->ge, &e); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetViscMat"
PetscErrorCode JacResGetViscMat(PMat pm)
{
	// assemble temperature preconditioner matrix
	// STEADY STATE solution is activated by setting time step to zero
	// COMPLETE SINGLE-POINT CONSTRIANT IMLEMENTATION !!!

	JacRes 		*jr;
	PMatBlock   *P;
	PetscInt	WI, WJ; // Indices of WMat


	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bvx, fvx, bvy, fvy, bvz, fvz;
	PetscScalar viscx, viscy, viscz, visc_center;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;

 	PetscScalar v[7], cf[6], invdt;

	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcv, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	jr 	 = pm->jr;
	fs   = jr->fs;
	bc   = jr->bc;
	num  = bc->tNumSPC;
	list = bc->tSPCList;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_VISC_TOTAL);

	// clear matrix coefficients
	ierr = MatZeroEntries(jr->Att); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcv,  &bcv); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcv[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcv[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcv[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcv[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcv[k-1][j][i] != DBL_MAX) cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcv[k+1][j][i] != DBL_MAX) cf[5] = -1.0; }

		// compute average viscosities
		bvx = (lk[k][j][i] + lk[k][j][Im1])/2.0;      fvx = (lk[k][j][i] + lk[k][j][Ip1])/2.0;
		bvy = (lk[k][j][i] + lk[k][Jm1][i])/2.0;      fvy = (lk[k][j][i] + lk[k][Jp1][i])/2.0;
		bvz = (lk[k][j][i] + lk[Km1][j][i])/2.0;      fvz = (lk[k][j][i] + lk[Kp1][j][i])/2.0;

		// inverse viscosities:
		bvx = 1/sqrt(bvx);  		fvx = 1/sqrt(fvx);
		bvy = 1/sqrt(bvy);  		fvy = 1/sqrt(fvy);
		bvz = 1/sqrt(bvz);  		fvz = 1/sqrt(fvz);

		// compute viscosities in cell-center
		viscx = (fvx+bvx)/2.0;
		viscy = (fvy+bvy)/2.0;
		viscz = (fvz+bvz)/2.0;
		visc_center = (viscx + viscy + viscz)/3.0;     																	 // /3.0 ??

		// store viscosities
		WI = i + (j-1)*nx + (k-1)*ny*nx; // compute indices of WMat
		WJ = j + (i-1)*ny + (k-1)*ny*nx;
		P->WMat[WI][WJ] = visc_center; 	 // store viscosity of the cell in weighting matrix

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
		ierr = MatSetValuesStencil(jr->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

		// NOTE! since only TPC are active, no SPC modification is necessary
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcv, &bcv);  CHKERRQ(ierr);

	// assemble K matrix
	ierr = MatAIJAssemble(jr->K, num, list, 1.0); CHKERRQ(ierr);

	// assemble C vector
	ierr = LumpMatrixToVector(pm); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LumpMatrixToVector"
PetscErrorCode LumpMatrixToVector(PMat pm)
{
	PMatBlock   *P;
	DOFIndex 	*dof;
	FDSTAG		*fs;
	JacRes		*jr;

	PetscInt	RowSum, lnv;
	Vec			C;
	Mat			WMat;

	C 		= P->C;
	WMat 	= P->WMat;

	fs  	= jr->fs;
	dof 	= &fs->dof;
	lnv 	= dof->lnv;

	// lumping
	for(i=1, i<lnv, i++)
	{
		RowSum = 0;
		for(j=1, j<lnv, j++)
		{
			RowSum = RowSum + WMat[i][j];
		}
		C[i] = RowSum;
	}

	P->C 	= C;

	PetscFunctionReturn(0);

}

