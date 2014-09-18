//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __lsove_h__
#define __lsove_h__
//---------------------------------------------------------------------------

// Stokes block preconditioner type
typedef enum
{
	STOKES_AL,  // augmented Lagrangian
	STOKES_BF,  // block factorization
	STOKES_MG,  // Galerkin multigrid
	STOKES_USER // user defined

} PCStokesType;

//---------------------------------------------------------------------------

// velocity preconditioner type
typedef enum
{
	VEL_MG,  // Galerkin multigrid
	VEL_USER // user-defined

} PCVelType;

//---------------------------------------------------------------------------

// Augmented Lagrangian preconditioner context
typedef struct
{
	Mat         Avv, Avp; // velocity sub-matrices
	Mat         Apv, App; // pressure sub-matrices
	Vec         M;        // diagonal penalty matrix
	PetscScalar pgamma;   // penalty parameter
	KSP         ksp;      // internal Krylov solver context

} PCStokesALCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesALCreate(PCStokesALCtx *alctx);

PetscErrorCode PCStokesALDestroy(PCStokesALCtx *alctx);

PetscErrorCode PCStokesALSetup(PCStokesALCtx *alctx);

PetscErrorCode PCStokesALApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// Block Factorization preconditioner context
typedef struct
{
	Mat   Avv, Avp; // velocity sub-matrices
	Mat   Apv, App; // pressure sub-matrices
	IS    isv, isp; // block index sets
	Mat   S;        // Schur complement preconditioner
	Mat   P;        // block matrix
	MGCtx mg;       // velocity multigrid context
	PC    pc;       // internal preconditioner context

} PCStokesBFCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesBFDestroy(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesBFSetup(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesBFApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
typedef struct
{
	Mat   P;      // monolithic matrix
	Mat   InvEta; // inverse viscosity matrix (Jacobian compensation)
	MGCtx mg;     // multigrid context
	PC    pc;     // internal preconditioner context

} PCStokesMGCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesMGDestroy(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesMGSetup(PCStokesALCtx *bfctx);

PetscErrorCode PCStokesMGApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// scatter block vectors to monolithic format & reverse
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode);

//---------------------------------------------------------------------------

/*
PetscErrorCode BlockMatCreate(BlockMat *bmat, FDSTAG *fs, Vec b);

PetscErrorCode BlockMatDestroy(BlockMat *bmat);

PetscErrorCode BlockMatCompute(BlockMat *bmat, FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx);

PetscErrorCode BlockMatClearSubMat(BlockMat *bmat);

PetscErrorCode BlockMatBlockToMonolithic(BlockMat *bmat, Vec b);

PetscErrorCode BlockMatMonolithicToBlock(BlockMat *bmat, Vec x);

//---------------------------------------------------------------------------

PetscErrorCode PowellHestenes(BlockMat *bmat, Vec r, Vec x);
*/

// block stop test
//PetscErrorCode KSPBlockStopTest(KSP ksp, PetscInt n, PetscScalar rnorm, KSPConvergedReason *reason, void *mctx);

//---------------------------------------------------------------------------

// fieldsplit preconditioner
//PetscErrorCode ApplyFieldSplit(PC pc, Vec x, Vec y);






//---------------------------------------------------------------------------
#endif
