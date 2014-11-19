//---------------------------------------------------------------------------
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#ifndef __multigrid_h__
#define __multigrid_h__
//---------------------------------------------------------------------------

typedef struct
{
	DOFIndex *dof;
	DM       DA_CEN;
	DM       DA_X,  DA_Y,  DA_Z;
	Mat      R;      // restriction matrix for every level (except coarsest)
	Mat      P;      // prolongation matrix for every level (except finest)

	PetscBool fine;   // finest level flag
	PetscBool coarse; // coarsest level flag

} MGLevel;

//---------------------------------------------------------------------------

PetscErrorCode MGLevelCreateFine(MGLevel *lvl, FDSTAG *fs);

PetscErrorCode MGLevelCreateCoarse(MGLevel *fine, MGLevel *coarse);


//---------------------------------------------------------------------------
// Galerkin multigrid level data structure

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

	PetscInt  ncors;  // number of coarsening steps (grids)
//	Mat      *R;      // restriction matrices for every level (except coarsest)
//	Mat      *P;      // prolongation matrices for every level (except finest)
//	FDSTAG   *mgfs;   // staggered grid for every level (except finest)
//	BCCtx    *mgbc;   // boundary condition contexts for every level (except finest)
	PC        pc;     // internal preconditioner context
	FDSTAG   *fs;     // finest level grid
	BCCtx    *bc;     // finest level boundary conditions
	idxtype   idxmod; // indexing mode

	MGLevel  *lvl;



} MG;

//---------------------------------------------------------------------------

PetscErrorCode MGCreate(MG *mg, FDSTAG *fs, BCCtx *bc, idxtype idxmod);

PetscErrorCode MGDestroy(MG *mg);

PetscErrorCode MGSetup(MG *mg, Mat A);

PetscErrorCode MGApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

PetscErrorCode MGSetDiagOnLevels(MG *mg);

PetscErrorCode MGDumpMat(MG *mg);

PetscErrorCode SetupRestrictStep(Mat R, FDSTAG *cors, FDSTAG *fine, BCCtx *bccors, idxtype idxmod);

PetscErrorCode SetupProlongStep(Mat P, FDSTAG *fine, FDSTAG *cors, BCCtx *bcfine, idxtype idxmod);

PetscErrorCode CheckMGRestrict(FDSTAG *fs, PetscInt *_ncors);

//---------------------------------------------------------------------------
#endif
