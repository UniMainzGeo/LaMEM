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
	PetscErrorCode ierr;
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
	ierr = MatDataCreateData(md);

	// compute index vectors
	ierr = MatDataComputeIndex(md); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataCreateData(MatData *md)
{
	FDSTAG   *fs;
	DOFIndex *dof;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  =  md->fs;
	dof = &fs->dof;

	// create index vectors
	ierr = DMCreateLocalVector(fs->DA_X,   &md->ivx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Y,   &md->ivy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_Z,   &md->ivz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(fs->DA_CEN, &md->ip);  CHKERRQ(ierr);

	// create parameter vectors
	ierr = DMCreateLocalVector (fs->DA_CEN, &md->Kb);    CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &md->rho);   CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &md->eta);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XY,  &md->etaxy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XZ,  &md->etaxz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_YZ,  &md->etayz); CHKERRQ(ierr);

	// single points constraints (SPC)
	ierr = makeIntArray(&md->SPCListMat, NULL, dof->ln); CHKERRQ(ierr);
	ierr = makeIntArray(&md->SPCListVec, NULL, dof->ln); CHKERRQ(ierr);

	if(md->coarsened)
	{
		// create boundary condition vectors
		ierr = DMCreateLocalVector(fs->DA_X,   &md->bcvx); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_Y,   &md->bcvy); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_Z,   &md->bcvz); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(fs->DA_CEN, &md->bcp);  CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataDestroy(MatData *md)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = VecDestroy(&md->ivx); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ivy); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ivz); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ip);  CHKERRQ(ierr);

	ierr = VecDestroy(&md->Kb);    CHKERRQ(ierr);
	ierr = VecDestroy(&md->rho);   CHKERRQ(ierr);
	ierr = VecDestroy(&md->eta);   CHKERRQ(ierr);
	ierr = VecDestroy(&md->etaxy); CHKERRQ(ierr);
	ierr = VecDestroy(&md->etaxz); CHKERRQ(ierr);
	ierr = VecDestroy(&md->etayz); CHKERRQ(ierr);

	ierr = PetscFree(md->SPCListMat); CHKERRQ(ierr);
	ierr = PetscFree(md->SPCListVec); CHKERRQ(ierr);

	if(md->coarsened)
	{
		ierr = FDSTAGDestroy(md->fs); CHKERRQ(ierr);
		ierr = PetscFree    (md->fs); CHKERRQ(ierr);

		ierr = VecDestroy(&md->bcvx); CHKERRQ(ierr);
		ierr = VecDestroy(&md->bcvy); CHKERRQ(ierr);
		ierr = VecDestroy(&md->bcvz); CHKERRQ(ierr);
		ierr = VecDestroy(&md->bcp);  CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataCoarsen(MatData *coarse, MatData *fine)
{
	PetscErrorCode ierr;
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
	ierr = PetscMalloc(sizeof(FDSTAG), &coarse->fs); CHKERRQ(ierr);

	// coarsen staggered grid
	ierr = FDSTAGCoarsen(coarse->fs, fine->fs); CHKERRQ(ierr);

	// allocate storage
	ierr = MatDataCreateData(coarse);

	// compute index vectors
	ierr = MatDataComputeIndex(coarse); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  =  md->fs;
	dof = &fs->dof;

	// set global indices of the local and ghost nodes (including boundary)
	ierr = VecSet(md->ivx, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ivy, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ivz, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ip,  -1.0); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   md->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   md->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   md->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->ip,  &ip);   CHKERRQ(ierr);

	//=======================================================
	// compute interlaced global numbering of the local nodes
	//=======================================================

	if     (md->idxmod == _IDX_COUPLED_)  { stv = dof->st;  stp = dof->st + dof->lnv; }
	else if(md->idxmod == _IDX_BLOCK_)    { stv = dof->stv; stp = dof->stp;           }

	//---------
	// X-points
	//---------
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivx[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivy[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivz[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ip[k][j][i] = (PetscScalar)stp; stp++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   md->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   md->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   md->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->ip,  &ip);   CHKERRQ(ierr);

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
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// update time step
	md->dt = jr->ts->dt;

	// update material parameters
	ierr = MatDataInitParam(md, jr); CHKERRQ(ierr);

	// update SPC constraints
	ierr = MatDataListSPC(md); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatDataRestrict(MatData *coarse, MatData *fine, PetscInt MG2D)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// update time step
	coarse->dt = fine->dt;

	// coarsen coordinates
	ierr = FDSTAGCoarsenCoord(coarse->fs, fine->fs); CHKERRQ(ierr);

	if(MG2D)
	{
		// coarsen material parameters
		ierr = MatDataRestrictParam2D(coarse, fine); CHKERRQ(ierr);

		// coarsen boundary conditions
		ierr = MatDataRestrictBC2D(coarse, fine); CHKERRQ(ierr);
	}
	else
	{
		// coarsen material parameters
		ierr = MatDataRestrictParam3D(coarse, fine); CHKERRQ(ierr);

		// coarsen boundary conditions
		ierr = MatDataRestrictBC3D(coarse, fine); CHKERRQ(ierr);
	}

	// update SPC constraints on coarse grid
	ierr = MatDataListSPC(coarse); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs = md->fs;
	dt = md->dt;

	// initialize ghost points
	ierr = VecZeroEntries(md->Kb);  CHKERRQ(ierr);
	ierr = VecZeroEntries(md->rho); CHKERRQ(ierr);
	ierr = VecZeroEntries(md->eta); CHKERRQ(ierr);

	// access parameter vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, md->Kb,    &Kb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->rho,   &rho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, md->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  md->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  md->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  md->etayz, &etayz); CHKERRQ(ierr);

	//-------------
	// cell centers
	//-------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		etaxy[k][j][i] = jr->svXYEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		etaxz[k][j][i] = jr->svXZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		etayz[k][j][i] = jr->svYZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->Kb,    &Kb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->rho,   &rho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, md->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  md->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  md->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  md->etayz, &etayz); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize index bounds in fine grid
	mnx = fine->fs->dsx.tnods - 1;
	mny = fine->fs->dsy.tnods - 1;
	mnz = fine->fs->dsz.tnods - 1;

	// initialize ghost points in coarse grid
	ierr = VecZeroEntries(coarse->Kb);  CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->rho); CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->eta); CHKERRQ(ierr);

	// access parameter vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->Kb,  &fKb);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->rho, &frho); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	// access parameter vectors in coarse grid
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, coarse->Kb,    &Kb);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, coarse->rho,   &rho);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, coarse->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_XY,  coarse->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_XZ,  coarse->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_YZ,  coarse->etayz, &etayz); CHKERRQ(ierr);

	//---------------------------
	// cell centers (coarse grid)
	//---------------------------
	ierr = DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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

	ierr = DMDAGetCorners(coarse->fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->Kb,  &fKb);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->rho, &frho); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	// restore access (coarse grid)
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->Kb,    &Kb);    CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->rho,   &rho);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_XY,  coarse->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_XZ,  coarse->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_YZ,  coarse->etayz, &etayz); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = VecZeroEntries(coarse->Kb);    CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->rho);   CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->eta);   CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->etaxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->etaxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(coarse->etayz); CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// mark all variables unconstrained
	ierr = VecSet(coarse->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcp,  DBL_MAX); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->bcp,  &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp);  CHKERRQ(ierr);

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->bcp,  &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp);  CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// mark all variables unconstrained
	ierr = VecSet(coarse->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(coarse->bcp,  DBL_MAX); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->bcp,  &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp);  CHKERRQ(ierr);

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(coarse->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(coarse->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->bcvx, &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->bcvy, &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->bcvz, &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->bcp,  &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(coarse->fs->DA_X,   coarse->bcvx, &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Y,   coarse->bcvy, &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_Z,   coarse->bcvz, &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coarse->fs->DA_CEN, coarse->bcp,  &cbcp);  CHKERRQ(ierr);

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

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	fs         = md->fs;
	dof        = &fs->dof;
	SPCListMat = md->SPCListMat;
	SPCListVec = md->SPCListVec;

	// clear constraints
	ierr = PetscMemzero(SPCListMat, sizeof(PetscInt)*(size_t)dof->ln); CHKERRQ(ierr);
	ierr = PetscMemzero(SPCListVec, sizeof(PetscInt)*(size_t)dof->ln); CHKERRQ(ierr);

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X, md->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y, md->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z, md->bcvz, &bcvz); CHKERRQ(ierr);

	iter   = 0;
	numSPC = 0;

	//---------
	// X points
	//---------

	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvx, SPCListVec, numSPC, iter)

        iter++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------

	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvy, SPCListVec, numSPC, iter)

        iter++;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------

	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		LIST_SPC_IND(bcvz, SPCListVec, numSPC, iter)

		iter++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X, md->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y, md->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z, md->bcvz, &bcvz); CHKERRQ(ierr);

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

