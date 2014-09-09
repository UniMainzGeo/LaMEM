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
	FDSTAG   *mgfs;  // staggered grid for every level (except finest)
	Mat      *R;     // restriction matrices for every level (except coarsest)
	Mat      *P;     // prolongation matrices for every level (except finest)
	BCCtx    *mgbc;  // boundary condition contexts for every level (except finest)
//	PC        pc;    // multigrid preconditioner


} MGCtx;
//---------------------------------------------------------------------------

PetscErrorCode MGCheckGrid(FDSTAG *fs, PetscInt *ncels);

PetscErrorCode MGCtxCreate(MGCtx *mg, FDSTAG *fs, PC pc, idxtype idxmod);

PetscErrorCode MGCtxDestroy(MGCtx *mg);

PetscErrorCode MGCtxSetup(MGCtx *mg, FDSTAG *fs, BCCtx *bc, idxtype idxmod);

PetscErrorCode MGCtxSetDiagOnLevels(MGCtx *mg, PC pc);

PetscErrorCode SetupRestrictStep(Mat R, FDSTAG *cors, FDSTAG *fine, BCCtx *bccors, idxtype idxmod);

PetscErrorCode SetupProlongStep(Mat P, FDSTAG *fine, FDSTAG *cors, BCCtx *bcfine, idxtype idxmod);

//---------------------------------------------------------------------------
#endif
