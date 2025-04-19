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
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "matData.h"
#include "tssolve.h"
#include "fdstag.h"
#include "bc.h"
#include "JacRes.h"
#include "tools.h"

//---------------------------------------------------------------------------
PetscErrorCode MatDataSetFromOptions(MatData *md)
{
	char pname[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set defaults
	sprintf(pname, "user");
	md->pgamma = 1.0;

	// read options
	ierr = PetscOptionsGetString(NULL, NULL, "-jp_type",    pname, _str_len_, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-jp_pgamma", &md->pgamma,       NULL); CHKERRQ(ierr);

	if     (!strcmp(pname, "mg") || !strcmp(pname, "user")) md->type = _MONOLITHIC_;
	else if(!strcmp(pname, "bf"))                           md->type = _BLOCK_;
	else    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes preconditioner type (jp_type): %s", pname);

	if(md->pgamma < 1.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Penalty parameter is less than unit (jp_pgamma)");

	PetscPrintf(PETSC_COMM_WORLD, "Preconditioner parameters: \n");
	if     (md->type == _MONOLITHIC_) PetscPrintf(PETSC_COMM_WORLD, "   Matrix type                   : monolithic\n");
	else if(md->type == _BLOCK_)      PetscPrintf(PETSC_COMM_WORLD, "   Matrix type                   : block\n");
	if     (md->pgamma > 1.0)         PetscPrintf(PETSC_COMM_WORLD, "   Penalty parameter (pgamma)    : %e\n", md->pgamma);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


PetscErrorCode MatDataCreate(MatData *md, JacRes *jr)
{
	FDSTAG      *fs;
	BCCtx       *bc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  = jr->fs;
	bc  = jr->bc;


	// read options
	ierr = MatDataSetFromOptions(md); CHKERRQ(ierr);

	md->fs = fs;

	ierr = MatDataCreateIndex(md); CHKERRQ(ierr);

	md->bcvx = bc->bcvx;
	md->bcvy = bc->bcvy;
	md->bcvz = bc->bcvz;
	md->bcp  = bc->bcp;

	md->numSPC  = bc->numSPC;  md->SPCList  = bc->SPCList;
	md->vNumSPC = bc->vNumSPC; md->vSPCList = bc->vSPCList;
	md->pNumSPC = bc->pNumSPC; md->pSPCList = bc->pSPCList;



	// create viscosity vectors
	ierr = DMCreateLocalVector (fs->DA_CEN, &md->eta);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XY,  &md->etaxy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_XZ,  &md->etaxz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(fs->DA_YZ,  &md->etayz); CHKERRQ(ierr);

	Vec         K, rho, eta, etaxy, etaxz, etayz; // parameter vectors



	md->dt     = jr->ts->dt;      // time step
	md->fssa   = jr->ctrl.FSSA;   // density gradient penalty parameter



	ierr = PetscMemcpy(md->grav, jr->ctrl.grav, 3*sizeof(PetscScalar)); CHKERRQ(ierr);


	md->coarsened = 0;

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

PetscErrorCode MatDataInitParam(MatData *md, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode MatDataDestroy(MatData *md)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;


	// destroy index vectors (ghosted)
	ierr = VecDestroy(&md->ivx); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ivy); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ivz); CHKERRQ(ierr);
	ierr = VecDestroy(&md->ip);  CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode MatDataCreateIndex(MatData *md)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, stv=0, stp=0;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// create index vectors (ghosted)
	ierr = DMCreateLocalVector(md->fs->DA_X,   &md->ivx); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(md->fs->DA_Y,   &md->ivy); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(md->fs->DA_Z,   &md->ivz); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(md->fs->DA_CEN, &md->ip);  CHKERRQ(ierr);

	// set global indices of the local and ghost nodes (including boundary)
	ierr = VecSet(md->ivx, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ivy, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ivz, -1.0); CHKERRQ(ierr);
	ierr = VecSet(md->ip,  -1.0); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(md->fs->DA_X,   md->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(md->fs->DA_Y,   md->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(md->fs->DA_Z,   md->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(md->fs->DA_CEN, md->ip,  &ip);   CHKERRQ(ierr);

	//=======================================================
	// compute interlaced global numbering of the local nodes
	//=======================================================
/*
	if     (md->type == _MONOLITHIC_)  { stv = md->fs->dof.st;  stp = dof->st + dof->lnv; }
	else if(md->type == _BLOCK_)       { stv = md->fs->dof.stv; stp = dof->stp;           }

	//---------
	// X-points
	//---------
	ierr = DMDAGetCorners(dof->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivx[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Y-points
	//---------
	ierr = DMDAGetCorners(dof->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivy[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// Z-points
	//---------
	ierr = DMDAGetCorners(dof->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ivz[k][j][i] = (PetscScalar)stv; stv++;
	}
	END_STD_LOOP

	//---------
	// P-points
	//---------
	ierr = DMDAGetCorners(dof->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		ip[k][j][i] = (PetscScalar)stp; stp++;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(dof->DA_X,   dof->ivx, &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_Y,   dof->ivy, &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_Z,   dof->ivz, &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dof->DA_CEN, dof->ip,  &ip);   CHKERRQ(ierr);

	// get ghost point indices
	LOCAL_TO_LOCAL(dof->DA_X,   dof->ivx)
	LOCAL_TO_LOCAL(dof->DA_Y,   dof->ivy)
	LOCAL_TO_LOCAL(dof->DA_Z,   dof->ivz)
	LOCAL_TO_LOCAL(dof->DA_CEN, dof->ip)

	// store index mode
	dof->idxmod = idxmod;
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode MatDataCoarsen(MatData *coarse, MatData *fine)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode MatDataRestricParam(MatData *coarse, MatData *fine)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

PetscErrorCode MatDataRestricBC(MatData *coarse, MatData *fine)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------



/*
PetscErrorCode MGLevelCreate(MGLevel *lvl, MGLevel *fine, FDSTAG *fs, BCCtx *bc)
{
	PetscInt ln=0, lnfine=0;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(!fine)
	{
		// setup top level
		lvl->fs     = fs;
		lvl->bcvx   = bc->bcvx;
		lvl->bcvy   = bc->bcvy;
		lvl->bcvz   = bc->bcvz;
		lvl->bcp    = bc->bcp;
		lvl->R      = NULL;
		lvl->P      = NULL;
	}
	else
	{
		// allocate staggered grid context
		ierr = PetscMalloc(sizeof(FDSTAG), &lvl->fs); CHKERRQ(ierr);

		// clear memory
		ierr = PetscMemzero(lvl->fs, sizeof(FDSTAG)); CHKERRQ(ierr);

		// coarsen staggered grid
		ierr = FDSTAGCoarsen(lvl->fs, fine->fs); CHKERRQ(ierr);

		// create restricted boundary condition vectors
		ierr = DMCreateLocalVector(lvl->fs->DA_X,   &lvl->bcvx);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->fs->DA_Y,   &lvl->bcvy);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->fs->DA_Z,   &lvl->bcvz);  CHKERRQ(ierr);
		ierr = DMCreateLocalVector(lvl->fs->DA_CEN, &lvl->bcp);   CHKERRQ(ierr);

		// compute index arrays
		ierr = DOFIndexCompute(&lvl->fs->dof, fine->fs->dof.idxmod); CHKERRQ(ierr);

		// get matrix sizes
		if     (lvl->fs->dof.idxmod == IDXCOUPLED)   { ln = lvl->fs->dof.ln;  lnfine = fine->fs->dof.ln;  }
		else if(lvl->fs->dof.idxmod == IDXUNCOUPLED) { ln = lvl->fs->dof.lnv; lnfine = fine->fs->dof.lnv; }

		// preallocate restriction & prolongation matrices
		// WARNING! CONSTANT SIZE PREALLOCATION
		// ADD VARIABLE PREALLOCATION FOR THESE MATRICES
		ierr = MatAIJCreate(ln,     lnfine, 12, NULL, 4, NULL, &lvl->R); CHKERRQ(ierr);
		ierr = MatAIJCreate(lnfine, ln,     8,  NULL, 7, NULL, &lvl->P); CHKERRQ(ierr);
	}

	// create viscosity vectors
	ierr = DMCreateLocalVector (lvl->fs->DA_CEN, &lvl->eta);   CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(lvl->fs->DA_XY,  &lvl->etaxy); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(lvl->fs->DA_XZ,  &lvl->etaxz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(lvl->fs->DA_YZ,  &lvl->etayz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}



//---------------------------------------------------------------------------
PetscErrorCode PMatCreate(PMat *p_pm, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	//========================================================================
	// create preconditioner matrix context
	//========================================================================

	PMat pm;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PMat), &pm); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pm, sizeof(p_PMat)); CHKERRQ(ierr);

	// set context
	pm->jr = jr;

	// read options
	ierr = PMatSetFromOptions(pm); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------

PetscErrorCode PMatDestroy(PMat pm)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = pm->Destroy(pm); CHKERRQ(ierr);
	ierr = PetscFree(pm);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

PetscErrorCode DOFIndexCompute(DOFIndex *dof, idxtype idxmod)
{

}

PetscErrorCode MGLevelRestrictParam(MGLevel *lvl, MGLevel *fine)
{
	// restrict parameters from fine to coarse grid

	PetscInt    mnx, mny, mnz;
	PetscInt    I, J, K, Im1, Jm1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar sum, ***eta, ***etaxy, ***etaxz, ***etayz, ***feta;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// initialize index bounds in fine grid
	mnx = fine->fs->dsx.tnods - 1;
	mny = fine->fs->dsy.tnods - 1;
	mnz = fine->fs->dsz.tnods - 1;

	// exchange ghost points in fine grid
	LOCAL_TO_LOCAL(fine->fs->DA_CEN, fine->eta)

	// access viscosity vector in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	// initialize ghost points in coarse grid
	ierr = VecZeroEntries(lvl->eta); CHKERRQ(ierr);

	// access viscosity vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->fs->DA_CEN, lvl->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_XY,  lvl->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_XZ,  lvl->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_YZ,  lvl->etayz, &etayz); CHKERRQ(ierr);

	//---------------------------
	// cell centers (coarse grid)
	//---------------------------
	ierr = DMDAGetCorners(lvl->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

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

	ierr = DMDAGetCorners (lvl->fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners (lvl->fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners (lvl->fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->eta, &feta); CHKERRQ(ierr);

	// restore access (coarse grid)
	ierr = DMDAVecRestoreArray(lvl->fs->DA_CEN, lvl->eta,   &eta);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_XY,  lvl->etaxy, &etaxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_XZ,  lvl->etaxz, &etaxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_YZ,  lvl->etayz, &etayz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------

PetscErrorCode MGLevelRestrictBC(MGLevel *lvl, MGLevel *fine, PetscBool no_restric_bc)
{
	// restrict boundary condition vectors from fine to coarse level

	PetscInt    I, J, K;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx,   ***ivy,   ***ivz,   ***ip;
	PetscScalar ***fbcvx, ***fbcvy, ***fbcvz, ***fbcp;
	PetscScalar ***cbcvx, ***cbcvy, ***cbcvz, ***cbcp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// mark all variables unconstrained
	ierr = VecSet(lvl->bcvx, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvy, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcvz, DBL_MAX); CHKERRQ(ierr);
	ierr = VecSet(lvl->bcp,  DBL_MAX); CHKERRQ(ierr);

	// check activation
	if(no_restric_bc == PETSC_TRUE) PetscFunctionReturn(0);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->fs->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->fs->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->fs->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->fs->dof.ip,  &ip);    CHKERRQ(ierr);

	// access boundary condition vectors in fine grid
	ierr = DMDAVecGetArray(fine->fs->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->fs->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	// access boundary condition vectors in coarse grid
	ierr = DMDAVecGetArray(lvl->fs->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(lvl->fs->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	ierr = DMDAGetCorners(lvl->fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(lvl->fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(lvl->fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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

	if(lvl->fs->dof.idxmod == IDXCOUPLED)
	{
		//-----------------------
		// P-points (coarse grid)
		//-----------------------
		ierr = DMDAGetCorners(lvl->fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

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
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->fs->dof.ivx, &ivx);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->fs->dof.ivy, &ivy);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->fs->dof.ivz, &ivz);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->fs->dof.ip,  &ip);    CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fine->fs->DA_X,   fine->bcvx,    &fbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Y,   fine->bcvy,    &fbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_Z,   fine->bcvz,    &fbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->fs->DA_CEN, fine->bcp,     &fbcp);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(lvl->fs->DA_X,    lvl->bcvx,     &cbcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_Y,    lvl->bcvy,     &cbcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_Z,    lvl->bcvz,     &cbcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(lvl->fs->DA_CEN,  lvl->bcp,      &cbcp);  CHKERRQ(ierr);

	// exchange ghost point constraints
	LOCAL_TO_LOCAL(lvl->fs->DA_X,   lvl->bcvx)
	LOCAL_TO_LOCAL(lvl->fs->DA_Y,   lvl->bcvy)
	LOCAL_TO_LOCAL(lvl->fs->DA_Z,   lvl->bcvz)
	LOCAL_TO_LOCAL(lvl->fs->DA_CEN, lvl->bcp)

	PetscFunctionReturn(0);
}
*/
