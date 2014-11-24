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
	Mat       R, P;              // restriction & prolongation matrices

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

//---------------------------------------------------------------------------

typedef struct
{
	// PETSc multigrid level numbering:
	// 0   - coarse grid
	// n-1 - fine grid
	// n   - number of levels

	// LaMEM multigrid level numbering:
	// 0   - fine grid
	// n-1 - coarse grid
	// n   - number of levels


	PetscInt  nlvl; // number of levels
	MGLevel  *lvls; // multigrid levles (except finest)

	PC        pc;   // internal preconditioner context

	FDSTAG   *fs;   // finest level grid
	BCCtx    *bc;   // finest level boundary conditions


} MG;

//---------------------------------------------------------------------------

PetscErrorCode MGCreate(MG *mg, FDSTAG *fs, BCCtx *bc);

PetscErrorCode MGDestroy(MG *mg);

PetscErrorCode MGSetup(MG *mg, Mat A);

PetscErrorCode MGApply(PC pc, Vec x, Vec y);

PetscErrorCode MGSetDiagOnLevels(MG *mg);

PetscErrorCode MGDumpMat(MG *mg);

//---------------------------------------------------------------------------

PetscErrorCode CheckMGRestrict(FDSTAG *fs, PetscInt *pncors);

// PetscErrorCode PreAllocRestrictStep(Mat R, FDSTAG *cors, FDSTAG *fine, BCCtx *bccors, idxtype idxmod);

// PetscErrorCode PreAllocProlongStep(Mat P, FDSTAG *fine, FDSTAG *cors, BCCtx *bcfine, idxtype idxmod);

PetscErrorCode SetupRestrictStep(Mat R, MGLevel *lvl, MGLevel *fine);

PetscErrorCode SetupProlongStep(Mat P, MGLevel *lvl, MGLevel *fine);


//---------------------------------------------------------------------------
#endif
