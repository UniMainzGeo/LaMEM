//---------------------------------------------------------------------------
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#ifndef __multigrid_h__
#define __multigrid_h__
//---------------------------------------------------------------------------
// Galerkin multigrid level data structure
typedef struct
{
	// PETSc multigrid level numbering:
	// 0   - coarse grid
	// n-1 - fine grid
	// n   - number of levels

	PetscInt  ncors; // number of coarsening steps (grids)
	Mat      *R;     // restriction matrices for every level (except coarsest)
	Mat      *P;     // prolongation matrices for every level (except finest)
	FDSTAG   *mgfs;  // staggered grid for every level (except finest)
	BCCtx    *mgbc;  // boundary condition contexts for every level (except finest)
	FDSTAG   *fs;    // finest level grid
	BCCtx    *bc;    // finest level boundary conditions

} MGCtx;
//---------------------------------------------------------------------------

PetscErrorCode MGCheckGrid(FDSTAG *fs, PetscInt *_ncors);

PetscErrorCode MGCtxCreate(MGCtx *mg, FDSTAG *fs, BCCtx *bc, PC pc, idxtype idxmod);

PetscErrorCode MGCtxDestroy(MGCtx *mg);

PetscErrorCode MGCtxSetup(MGCtx *mg, idxtype idxmod);

PetscErrorCode MGCtxSetDiagOnLevels(MGCtx *mg, PC pcmg);

PetscErrorCode SetupRestrictStep(Mat R, FDSTAG *cors, FDSTAG *fine, BCCtx *bccors, idxtype idxmod);

PetscErrorCode SetupProlongStep(Mat P, FDSTAG *fine, FDSTAG *cors, BCCtx *bcfine, idxtype idxmod);

//---------------------------------------------------------------------------
#endif
