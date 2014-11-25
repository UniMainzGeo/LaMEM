//---------------------------------------------------------------------------
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#ifndef __multigrid_h__
#define __multigrid_h__
//---------------------------------------------------------------------------

// Galerkin multigrid level data structure

typedef struct
{
	DM        DA_CEN;            // central points array
	DM        DA_X, DA_Y, DA_Z;  // face points arrays
	DOFIndex  dof;               // indexing vectors
	Mat       R, P;              // restriction & prolongation operators (not set on finest grid)

	// ******** fine level ************
	//     |                   ^
	//     R-matrix            |
	//     |                   P-matrix
	//     v                   |
	// ******** this level ************

} MGLevel;

//---------------------------------------------------------------------------

PetscErrorCode MGLevelCreate(MGLevel *lvl, MGLevel *fine, FDSTAG *fs);

PetscErrorCode MGLevelDestroy(MGLevel *lvl);

PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelSetupProlong(MGLevel *lvl, MGLevel *fine);

// PetscErrorCode MGLevelAllocRestrict(MGLevel *lvl, MGLevel *fine);

// PetscErrorCode MGLevelAllocProlong(MGLevel *lvl, MGLevel *fine);

//---------------------------------------------------------------------------

typedef struct
{
	// PETSc level numbering (inverse w.r.t. coarsening sequence):
	// 0   - coarse grid
	// n-1 - fine grid
	// R & P matrices connect with coarse level (not set on coarsest level).
	// Coarsening step yields coarse grid operator. Own operator is prescribed.

	// LaMEM level numbering (natural w.r.t. coarsening sequence):
	// 0   - fine grid
	// n-1 - coarse grid
	// R & P matrices connect with fine level (not set on finest grid).
	// Coarsening step yields own operator. Fine level operator is prescribed.

	PetscInt  nlvl; // number of levels
	MGLevel  *lvls; // multigrid levles

	PC        pc;   // internal preconditioner context

	FDSTAG   *fs;   // finest level grid
	BCCtx    *bc;   // finest level boundary conditions

} MG;

//---------------------------------------------------------------------------

PetscErrorCode MGCreate(MG *mg, FDSTAG *fs, BCCtx *bc);

PetscErrorCode MGDestroy(MG *mg);

PetscErrorCode MGSetup(MG *mg, Mat A);

PetscErrorCode MGApply(PC pc, Vec x, Vec y);

//PetscErrorCode MGSetDiagOnLevels(MG *mg);

PetscErrorCode MGDumpMat(MG *mg);

PetscErrorCode MGGetNumLevels(MG *mg);

//---------------------------------------------------------------------------
#endif
