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
//..................... MATRIX EVALUATION CONTEXT  .........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "matData.h"
#include "tssolve.h"
#include "fdstag.h"
#include "bc.h"
#include "JacRes.h"
#include "tools.h"

//---------------------------------------------------------------------------
PetscErrorCode MatDataCreate(MatData *md, JacRes *jr, idxtype idxmod)
{
	
	PetscFunctionBeginUser;

	// set coarse grid flag
	md->coarsened = 0;

	// copy data
	md->fs       = jr->fs;
	md->bcvx     = jr->bc->bcvx;
	md->bcvy     = jr->bc->bcvy;
	md->bcvz     = jr->bc->bcvz;
	md->bcp      = jr->bc->bcp;
	md->fssa     = jr->ctrl.FSSA;
	md->rescal   = jr->ctrl.rescal;
	md->grav[0]  = jr->ctrl.grav[0];
	md->grav[1]  = jr->ctrl.grav[1];
	md->grav[2]  = jr->ctrl.grav[2];
	md->idxmod   = idxmod;

	// allocate storage
	PetscCall(MatDataCreateData(md));

	// compute index vectors
	PetscCall(MatDataComputeIndex(md));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataCreateData(MatData *md)
{
	FDSTAG   *fs;
	DOFIndex *dof;

	
	PetscFunctionBeginUser;

	fs  =  md->fs;
	dof = &fs->dof;

	// create index vectors
	PetscCall(DMCreateLocalVector(fs->DA_X,   &md->ivx));
	PetscCall(DMCreateLocalVector(fs->DA_Y,   &md->ivy));
	PetscCall(DMCreateLocalVector(fs->DA_Z,   &md->ivz));
	PetscCall(DMCreateLocalVector(fs->DA_CEN, &md->ip));

	// create parameter vectors
	PetscCall(DMCreateLocalVector (fs->DA_CEN, &md->Kb));
	PetscCall(DMCreateLocalVector (fs->DA_CEN, &md->rho));
	PetscCall(DMCreateLocalVector (fs->DA_CEN, &md->eta));
	PetscCall(DMCreateGlobalVector(fs->DA_XY,  &md->etaxy));
	PetscCall(DMCreateGlobalVector(fs->DA_XZ,  &md->etaxz));
	PetscCall(DMCreateGlobalVector(fs->DA_YZ,  &md->etayz));

	// single points constraints (SPC)
	PetscCall(makeIntArray(&md->SPCListMat, NULL, dof->ln));
	PetscCall(makeIntArray(&md->SPCListVec, NULL, dof->ln));

	if(md->coarsened)
	{
		// create boundary condition vectors
		PetscCall(DMCreateLocalVector(fs->DA_X,   &md->bcvx));
		PetscCall(DMCreateLocalVector(fs->DA_Y,   &md->bcvy));
		PetscCall(DMCreateLocalVector(fs->DA_Z,   &md->bcvz));
		PetscCall(DMCreateLocalVector(fs->DA_CEN, &md->bcp));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataDestroy(MatData *md)
{
	
	PetscFunctionBeginUser;

	PetscCall(VecDestroy(&md->ivx));
	PetscCall(VecDestroy(&md->ivy));
	PetscCall(VecDestroy(&md->ivz));
	PetscCall(VecDestroy(&md->ip));

	PetscCall(VecDestroy(&md->Kb));
	PetscCall(VecDestroy(&md->rho));
	PetscCall(VecDestroy(&md->eta));
	PetscCall(VecDestroy(&md->etaxy));
	PetscCall(VecDestroy(&md->etaxz));
	PetscCall(VecDestroy(&md->etayz));

	PetscCall(PetscFree(md->SPCListMat));
	PetscCall(PetscFree(md->SPCListVec));

	if(md->coarsened)
	{
		PetscCall(FDSTAGDestroy(md->fs));
		PetscCall(PetscFree    (md->fs));

		PetscCall(VecDestroy(&md->bcvx));
		PetscCall(VecDestroy(&md->bcvy));
		PetscCall(VecDestroy(&md->bcvz));
		PetscCall(VecDestroy(&md->bcp));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataCoarsen(MatData *coarse, MatData *fine)
{
	
	PetscFunctionBeginUser;

	// set coarse grid flag
	coarse->coarsened = 1;

	// copy data
	coarse->idxmod  = fine->idxmod;
	coarse->fssa    = fine->fssa;
	coarse->grav[0] = fine->grav[0];
	coarse->grav[1] = fine->grav[1];
	coarse->grav[2] = fine->grav[2];

	// allocate staggered grid context
	PetscCall(PetscMalloc(sizeof(FDSTAG), &coarse->fs));

	// coarsen staggered grid
	PetscCall(FDSTAGCoarsen(coarse->fs, fine->fs));

	// allocate storage
	PetscCall(MatDataCreateData(coarse));

	// compute index vectors
	PetscCall(MatDataComputeIndex(coarse));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataComputeIndex(MatData *md)
{
	// compute & set global indices of local & ghost nodes

	// **********************************************************************
	// NOTE:
	// for the ghost points, store negative global index of the primary DOF
	// instead of -1
	// **********************************************************************

	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, stv=0, stp=0;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	
	PetscFunctionBeginUser;

	fs  =  md->fs;
	dof = &fs->dof;

	// set global indices of the local and ghost nodes (including boundary)
	PetscCall(VecSet(md->ivx, -1.0));
	PetscCall(VecSet(md->ivy, -1.0));
	PetscCall(VecSet(md->ivz, -1.0));
	PetscCall(VecSet(md->ip,  -1.0));

	// access index vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   md->ivx, &ivx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   md->ivy, &ivy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   md->ivz, &ivz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->ip,  &ip));

	//=======================================================
	// compute interlaced global numbering of the local nodes
	//=======================================================

	if     (md->idxmod == _IDX_COUPLED_)  { stv = dof->st;  stp = dof->st + dof->lnv; }
	else if(md->idxmod == _IDX_BLOCK_)    { stv = dof->stv; stp = dof->stp;           }

	//---------
	// X-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		ivx[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		ivy[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		ivz[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		ip[k][j][i] = (PetscScalar)stp; stp++;
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   md->ivx, &ivx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   md->ivy, &ivy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   md->ivz, &ivz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->ip,  &ip));

	// get ghost point indices
	LOCAL_TO_LOCAL(fs->DA_X,   md->ivx)
	LOCAL_TO_LOCAL(fs->DA_Y,   md->ivy)
	LOCAL_TO_LOCAL(fs->DA_Z,   md->ivz)
	LOCAL_TO_LOCAL(fs->DA_CEN, md->ip)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataSetup(MatData *md, JacRes *jr)
{
	
	PetscFunctionBeginUser;

	// update time step
	md->dt = jr->ts->dt;

	// update material parameters
	PetscCall(MatDataInitParam(md, jr));

	// update SPC constraints
	PetscCall(MatDataListSPC(md));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrict(MatData *coarse, MatData *fine, PetscInt MG2D)
{
	
	PetscFunctionBeginUser;

	// update time step
	coarse->dt = fine->dt;

	// coarsen coordinates
	PetscCall(FDSTAGCoarsenCoord(coarse->fs, fine->fs));

	if(MG2D)
	{
		// coarsen material parameters
		PetscCall(MatDataRestrictParam2D(coarse, fine));

		// coarsen boundary conditions
		PetscCall(MatDataRestrictBC2D(coarse, fine));
	}
	else
	{
		// coarsen material parameters
		PetscCall(MatDataRestrictParam3D(coarse, fine));

		// coarsen boundary conditions
		PetscCall(MatDataRestrictBC3D(coarse, fine));
	}

	// update SPC constraints on coarse grid
	PetscCall(MatDataListSPC(coarse));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataInitParam(MatData *md, JacRes *jr)
{
	// initialize parameters on fine grid

	FDSTAG     *fs;
 	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar ***Kb, ***rho, ***eta, ***etaxy, ***etaxz, ***etayz;
	PetscScalar IKdt, dt;

	
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;
	dt = md->dt;

	// initialize ghost points
	PetscCall(VecZeroEntries(md->Kb));
	PetscCall(VecZeroEntries(md->rho));
	PetscCall(VecZeroEntries(md->eta));

	// access parameter vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->Kb,    &Kb));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->rho,   &rho));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, md->eta,   &eta));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  md->etaxy, &etaxy));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &etaxz));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  md->etayz, &etayz));

	//-------------
	// cell centers
	//-------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		IKdt         = jr->svCell[iter].svBulk.IKdt;
		rho[k][j][i] = jr->svCell[iter].svBulk.rho;
		eta[k][j][i] = jr->svCell[iter].svDev.eta;
		Kb [k][j][i] = 1.0/(IKdt*dt);

		iter++;
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		etaxy[k][j][i] = jr->svXYEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		etaxz[k][j][i] = jr->svXZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		etayz[k][j][i] = jr->svYZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &Kb));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &rho));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &eta));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &etaxy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &etaxz));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &etayz));

	// exchange ghost points
	LOCAL_TO_LOCAL(md->fs->DA_CEN, md->Kb)
	LOCAL_TO_LOCAL(md->fs->DA_CEN, md->rho)
	LOCAL_TO_LOCAL(md->fs->DA_CEN, md->eta)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrictParam3D(MatData *coarse, MatData *fine)
{
	// restrict parameters from fine to coarse grid

	PetscInt    mnx, mny, mnz;
	PetscInt    I, J, K, Im1, Jm1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***fKb, ***frho, ***feta;
	PetscScalar ***Kb,  ***rho,  ***eta, ***etaxy, ***etaxz, ***etayz;
	PetscScalar sum;

	
	PetscFunctionBeginUser;

	// initialize index bounds in fine grid
	mnx = fine->fs->dsx.tnods - 1;
	mny = fine->fs->dsy.tnods - 1;
	mnz = fine->fs->dsz.tnods - 1;

	// initialize ghost points in coarse grid
	PetscCall(VecZeroEntries(coarse->Kb));
	PetscCall(VecZeroEntries(coarse->rho));
	PetscCall(VecZeroEntries(coarse->eta));

	// access parameter vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->Kb,  &fKb));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->rho, &frho));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->eta, &feta));

	// access parameter vectors in coarse grid
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, coarse->Kb,    &Kb));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, coarse->rho,   &rho));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, coarse->eta,   &eta));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_XY,  coarse->etaxy, &etaxy));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_XZ,  coarse->etaxz, &etaxz));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_YZ,  coarse->etayz, &etayz));

	//---------------------------
	// cell centers (coarse grid)
	//---------------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// compute average bulk modulus
		sum = fKb[K  ][J  ][I  ]
		+     fKb[K  ][J  ][I+1]
		+     fKb[K  ][J+1][I  ]
		+     fKb[K  ][J+1][I+1]
		+     fKb[K+1][J  ][I  ]
		+     fKb[K+1][J  ][I+1]
		+     fKb[K+1][J+1][I  ]
		+     fKb[K+1][J+1][I+1];

		Kb[k][j][i] = sum/8.0;

		// compute average density
		sum = frho[K  ][J  ][I  ]
		+     frho[K  ][J  ][I+1]
		+     frho[K  ][J+1][I  ]
		+     frho[K  ][J+1][I+1]
		+     frho[K+1][J  ][I  ]
		+     frho[K+1][J  ][I+1]
		+     frho[K+1][J+1][I  ]
		+     frho[K+1][J+1][I+1];

		rho[k][j][i] = sum/8.0;

		// compute average viscosity
		sum = feta[K  ][J  ][I  ]
		+     feta[K  ][J  ][I+1]
		+     feta[K  ][J+1][I  ]
		+     feta[K  ][J+1][I+1]
		+     feta[K+1][J  ][I  ]
		+     feta[K+1][J  ][I+1]
		+     feta[K+1][J+1][I  ]
		+     feta[K+1][J+1][I+1];

		eta[k][j][i] = sum/8.0;
	}
	END_STD_LOOP

	//-----------------------------
	// xy edge points (coarse grid)
	//-----------------------------

	PetscCall(DMDAGetCorners(coarse->fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i; Im1 = I-1;
		J = 2*j; Jm1 = J-1;
		K = 2*k;

		// check index bounds
		if(I   > mnx) { I--;   }
		if(J   > mny) { J--;   }
		if(Im1 < 0)   { Im1++; }
		if(Jm1 < 0)   { Jm1++; }

		// compute average viscosity
		sum = feta[K  ][Jm1][Im1]
		+     feta[K  ][J  ][Im1]
		+     feta[K  ][Jm1][I  ]
		+     feta[K  ][J  ][I  ]
		+     feta[K+1][Jm1][Im1]
		+     feta[K+1][J  ][Im1]
		+     feta[K+1][Jm1][I  ]
		+     feta[K+1][J  ][I  ];

		etaxy[k][j][i] = sum/8.0;
	}
	END_STD_LOOP

	//-----------------------------
	// xz edge points (coarse grid)
	//-----------------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i; Im1 = I-1;
		J = 2*j;
		K = 2*k; Km1 = K-1;

		// check index bounds
		if(I   > mnx) { I--;   }
		if(K   > mnz) { K--;   }
		if(Im1 < 0)   { Im1++; }
		if(Km1 < 0)   { Km1++; }

		// compute average viscosity
		sum = feta[Km1][J  ][Im1]
		+     feta[K  ][J  ][Im1]
		+     feta[Km1][J  ][I  ]
		+     feta[K  ][J  ][I  ]
		+     feta[Km1][J+1][Im1]
		+     feta[K  ][J+1][Im1]
		+     feta[Km1][J+1][I  ]
		+     feta[K  ][J+1][I  ];

		etaxz[k][j][i] = sum/8.0;
	}
	END_STD_LOOP

	//-----------------------------
	// yz edge points (coarse grid)
	//-----------------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j; Jm1 = J-1;
		K = 2*k; Km1 = K-1;

		// check index bounds
		if(J   > mny) { J--;   }
		if(K   > mnz) { K--;   }
		if(Jm1 < 0)   { Jm1++; }
		if(Km1 < 0)   { Km1++; }

		// compute average viscosity
		sum = feta[Km1][Jm1][I  ]
		+     feta[K  ][Jm1][I  ]
		+     feta[Km1][J  ][I  ]
		+     feta[K  ][J  ][I  ]
		+     feta[Km1][Jm1][I+1]
		+     feta[K  ][Jm1][I+1]
		+     feta[Km1][J  ][I+1]
		+     feta[K  ][J  ][I+1];

		etayz[k][j][i] = sum/8.0;
	}
	END_STD_LOOP

	// restore access (fine grid)
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->Kb,  &fKb));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->rho, &frho));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->eta, &feta));

	// restore access (coarse grid)
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->Kb,    &Kb));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->rho,   &rho));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->eta,   &eta));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_XY,  coarse->etaxy, &etaxy));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_XZ,  coarse->etaxz, &etaxz));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_YZ,  coarse->etayz, &etayz));

	// exchange ghost points in coarse grid
	LOCAL_TO_LOCAL(coarse->fs->DA_CEN, coarse->Kb)
	LOCAL_TO_LOCAL(coarse->fs->DA_CEN, coarse->rho)
	LOCAL_TO_LOCAL(coarse->fs->DA_CEN, coarse->eta)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrictParam2D(MatData *coarse, MatData *fine)
{
	// restrict parameters from fine to coarse grid
	//
	// WARNING!
	//
	// material parameters are not accessed on coarse grids in 2D case
	// implement 2D material parameter coarsening if this changes

	UNUSED(fine);

	
	PetscFunctionBeginUser;

	PetscCall(VecZeroEntries(coarse->Kb));
	PetscCall(VecZeroEntries(coarse->rho));
	PetscCall(VecZeroEntries(coarse->eta));
	PetscCall(VecZeroEntries(coarse->etaxy));
	PetscCall(VecZeroEntries(coarse->etaxz));
	PetscCall(VecZeroEntries(coarse->etayz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrictBC3D(MatData *coarse, MatData *fine)
{
	// restrict boundary condition vectors from fine grid to coarse grid

	// Constrained DOF stores parent DOF index in the boundary condition vector.
	// Parent DOF index is the only nonzero that is set in the row of R-matrix
	// and column of P-matrix to impose the constraints in a coarse grid operator
	// automatically. The finest grid uses standard boundary condition vectors.

	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	
	PetscFunctionBeginUser;

	// mark all variables unconstrained
	PetscCall(VecSet(coarse->bcvx, DBL_MAX));
	PetscCall(VecSet(coarse->bcvy, DBL_MAX));
	PetscCall(VecSet(coarse->bcvz, DBL_MAX));
	PetscCall(VecSet(coarse->bcp,  DBL_MAX));

	// access index vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fine->ivx, &ivx));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fine->ivy, &ivy));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fine->ivz, &ivz));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->ip,  &ip));

	// access boundary condition vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fine->bcvx, &fbcvx));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->bcp,  &fbcp));

	// access boundary condition vectors in coarse grid
	PetscCall(DMDAVecGetArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp));

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict constraint
		if(fbcvx[K  ][J  ][I] != DBL_MAX
		&& fbcvx[K  ][J+1][I] != DBL_MAX
		&& fbcvx[K+1][J  ][I] != DBL_MAX
		&& fbcvx[K+1][J+1][I] != DBL_MAX)
		{
			// store parent DOF index
			cbcvx[k][j][i] = ivx[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict constraint
		if(fbcvy[K  ][J][I  ] != DBL_MAX
		&& fbcvy[K  ][J][I+1] != DBL_MAX
		&& fbcvy[K+1][J][I  ] != DBL_MAX
		&& fbcvy[K+1][J][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvy[k][j][i] = ivy[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict constraint
		if(fbcvz[K][J  ][I  ] != DBL_MAX
		&& fbcvz[K][J  ][I+1] != DBL_MAX
		&& fbcvz[K][J+1][I  ] != DBL_MAX
		&& fbcvz[K][J+1][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvz[k][j][i] = ivz[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// P-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		// restrict constraint
		if(fbcp[K  ][J  ][I  ] != DBL_MAX
		&& fbcp[K  ][J  ][I+1] != DBL_MAX
		&& fbcp[K  ][J+1][I  ] != DBL_MAX
		&& fbcp[K  ][J+1][I+1] != DBL_MAX
		&& fbcp[K+1][J  ][I  ] != DBL_MAX
		&& fbcp[K+1][J  ][I+1] != DBL_MAX
		&& fbcp[K+1][J+1][I  ] != DBL_MAX
		&& fbcp[K+1][J+1][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcp[k][j][i] = ip[K][J][I];
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fine->ivx, &ivx));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fine->ivy, &ivy));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fine->ivz, &ivz));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->ip,  &ip));

	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fine->bcvx, &fbcvx));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->bcp,  &fbcp));

	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp));

	// exchange ghost point constraints
	LOCAL_TO_LOCAL(coarse->fs->DA_X,   coarse->bcvx)
	LOCAL_TO_LOCAL(coarse->fs->DA_Y,   coarse->bcvy)
	LOCAL_TO_LOCAL(coarse->fs->DA_Z,   coarse->bcvz)
	LOCAL_TO_LOCAL(coarse->fs->DA_CEN, coarse->bcp)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrictBC2D(MatData *coarse, MatData *fine)
{
	// restrict boundary condition vectors from fine grid to coarse grid

	// Constrained DOF stores parent DOF index in the boundary condition vector.
	// Parent DOF index is the only nonzero that is set in the row of R-matrix
	// and column of P-matrix to impose the constraints in a coarse grid operator
	// automatically. The finest grid uses standard boundary condition vectors.

	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	
	PetscFunctionBeginUser;

	// mark all variables unconstrained
	PetscCall(VecSet(coarse->bcvx, DBL_MAX));
	PetscCall(VecSet(coarse->bcvy, DBL_MAX));
	PetscCall(VecSet(coarse->bcvz, DBL_MAX));
	PetscCall(VecSet(coarse->bcp,  DBL_MAX));

	// access index vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fine->ivx, &ivx));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fine->ivy, &ivy));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fine->ivz, &ivz));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->ip,  &ip));

	// access boundary condition vectors in fine grid
	PetscCall(DMDAVecGetArray(fine->fs->DA_X,   fine->bcvx, &fbcvx));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy));
	PetscCall(DMDAVecGetArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz));
	PetscCall(DMDAVecGetArray(fine->fs->DA_CEN, fine->bcp,  &fbcp));

	// access boundary condition vectors in coarse grid
	PetscCall(DMDAVecGetArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz));
	PetscCall(DMDAVecGetArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp));

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J =   j;
		K = 2*k;

		// restrict constraint
		if(fbcvx[K  ][J][I] != DBL_MAX
		&& fbcvx[K+1][J][I] != DBL_MAX)
		{
			// store parent DOF index
			cbcvx[k][j][i] = ivx[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J =   j;
		K = 2*k;

		// restrict constraint
		if(fbcvy[K  ][J][I  ] != DBL_MAX
		&& fbcvy[K  ][J][I+1] != DBL_MAX
		&& fbcvy[K+1][J][I  ] != DBL_MAX
		&& fbcvy[K+1][J][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvy[k][j][i] = ivy[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J =   j;
		K = 2*k;

		// restrict constraint
		if(fbcvz[K][J][I  ] != DBL_MAX
		&& fbcvz[K][J][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcvz[k][j][i] = ivz[K][J][I];
		}
	}
	END_STD_LOOP

	//-----------------------
	// P-points (coarse grid)
	//-----------------------
	PetscCall(DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J =   j;
		K = 2*k;

		// restrict constraint
		if(fbcp[K  ][J][I  ] != DBL_MAX
		&& fbcp[K  ][J][I+1] != DBL_MAX
		&& fbcp[K+1][J][I  ] != DBL_MAX
		&& fbcp[K+1][J][I+1] != DBL_MAX)
		{
			// store parent DOF index
			cbcp[k][j][i] = ip[K][J][I];
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fine->ivx, &ivx));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fine->ivy, &ivy));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fine->ivz, &ivz));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->ip,  &ip));

	PetscCall(DMDAVecRestoreArray(fine->fs->DA_X,   fine->bcvx, &fbcvx));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz));
	PetscCall(DMDAVecRestoreArray(fine->fs->DA_CEN, fine->bcp,  &fbcp));

	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz));
	PetscCall(DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp));

	// exchange ghost point constraints
	LOCAL_TO_LOCAL(coarse->fs->DA_X,   coarse->bcvx)
	LOCAL_TO_LOCAL(coarse->fs->DA_Y,   coarse->bcvy)
	LOCAL_TO_LOCAL(coarse->fs->DA_Z,   coarse->bcvz)
	LOCAL_TO_LOCAL(coarse->fs->DA_CEN, coarse->bcp)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataListSPC(MatData *md)
{
	// create SPC constraint lists

	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscInt    vShift=0, pShift=0;
	PetscInt    iter, numSPC, *SPCListMat, *SPCListVec;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz;

	
	PetscFunctionBeginUser;

	// access context
	fs         = md->fs;
	dof        = &fs->dof;
	SPCListMat = md->SPCListMat;
	SPCListVec = md->SPCListVec;

	// clear constraints
	PetscCall(PetscMemzero(SPCListMat, sizeof(PetscInt)*(size_t)dof->ln));
	PetscCall(PetscMemzero(SPCListVec, sizeof(PetscInt)*(size_t)dof->ln));

	// access vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, md->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, md->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, md->bcvz, &bcvz));

	iter   = 0;
	numSPC = 0;

	//---------
	// X points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvx, SPCListVec, numSPC, iter)

        iter++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvy, SPCListVec, numSPC, iter)

        iter++;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvz, SPCListVec, numSPC, iter)

		iter++;
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, md->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, md->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, md->bcvz, &bcvz));

	// store velocity list
	md->vNumSPC     = numSPC;
	md->vSPCListVec = SPCListVec;
	md->vSPCListMat = SPCListMat;

	// WARNING! primary pressure constraints are not implemented, otherwise compute here
	md->pNumSPC = 0;

	// store total number of SPC constraints
	md->numSPC = numSPC;

	// get index shift from local index space (vector) to global index space (matrix)
	if(md->idxmod == _IDX_COUPLED_) { vShift = dof->st;  pShift = dof->st;             }
	if(md->idxmod == _IDX_BLOCK_)   { vShift = dof->stv; pShift = dof->stp - dof->lnv; }

	for(i = 0; i < md->vNumSPC; i++) md->vSPCListMat[i] = md->vSPCListVec[i] + vShift;
	for(i = 0; i < md->pNumSPC; i++) md->pSPCListMat[i] = md->pSPCListVec[i] + pShift;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

