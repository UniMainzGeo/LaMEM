//---------------------------------------------------------------------------
//.......................   LINEAR SOLVER ROUTINES   ........................
//---------------------------------------------------------------------------
#ifndef __lsolve_h__
#define __lsolve_h__
//---------------------------------------------------------------------------
// Stokes block preconditioner type
typedef enum
{
	STOKES_DIRECT,          // direct solver for artificial compressibility
	STOKES_POWELL_HESTENES, // iterative solver for artificial compressibility
	STOKES_SCHUR,           // pressure Schur complement reduction
	STOKES_FIELDSPLIT,      // block factorization preconditioner
	STOKES_GALERKIN_MG,     // coupled multigrid with Galerkin coarsening
	// ... add more
} StokesPCType;
//---------------------------------------------------------------------------
// velocity block preconditioner type (only relevant for SCHUR & FIELDSPLIT)
typedef enum
{
	VEL_DIRECT,      // direct solver
	VEL_GALERKIN_MG, // coupled multigrid
	// ... add more
} VelPCType;

//---------------------------------------------------------------------------


// Block equation system matrix
typedef struct
{
	Mat           Avv, Avp;     // velocity matrices
	Mat           Apv;          // pressure matrices
	Vec           kIM;          // penalized inverse pressure mass matrix

	Vec           wv,  wp;      // velocity & pressure work vectors
	IS            isv, isp;     // index sets
	VecScatter    vsv, vsp;     // vector scatter contexts

	Mat           A;            // block matrix (MATNEST) or monolithic
	Mat           P;            // block matrix (MATNEST) or monolithic


//	PetscScalar   nrmVV, nrmVP; // norms of velocity matrices
//	PetscScalar   nrmPV, nrmPP; // norms of pressure matrices
//	PetscScalar   nrmRHSV;      // norm of velocity right-hand-side
//	PetscScalar   nrmRHSP;      // norm of pressure right-hand-side
//	PetscScalar   rtolV, atolV; // velocity tolerances (relative & absolute)
//	PetscScalar   rtolP, atolP; // pressure tolerances (relative & absolute)
	PC            pc;           // field-split preconditioner

} BlockMat;
//---------------------------------------------------------------------------

PetscErrorCode BlockMatCreate(BlockMat *bmat, FDSTAG *fs, Vec b);

PetscErrorCode BlockMatDestroy(BlockMat *bmat);

PetscErrorCode BlockMatCompute(BlockMat *bmat, FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx);

PetscErrorCode BlockMatClearSubMat(BlockMat *bmat);

PetscErrorCode BlockMatBlockToMonolithic(BlockMat *bmat, Vec b);

PetscErrorCode BlockMatMonolithicToBlock(BlockMat *bmat, Vec x);

//---------------------------------------------------------------------------

PetscErrorCode PowellHestenes(BlockMat *bmat, Vec r, Vec x);


// block stop test
//PetscErrorCode KSPBlockStopTest(KSP ksp, PetscInt n, PetscScalar rnorm, KSPConvergedReason *reason, void *mctx);

//---------------------------------------------------------------------------

// fieldsplit preconditioner
//PetscErrorCode ApplyFieldSplit(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------
#endif
