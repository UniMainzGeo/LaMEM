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
	STOKES_AMG, // coupled algebraic multigrid
	STOKES_USER // user defined

} PCStokesType;

//---------------------------------------------------------------------------

// velocity preconditioner type
typedef enum
{
	VEL_MG,  // Galerkin multigrid
	VEL_AMG, // algebraic multigrid
	VEL_USER // user-defined

} PCVelType;

//---------------------------------------------------------------------------

// solution context
typedef struct
{
	FDSTAG    *fs;
	BCCtx     *bc;
	JacResCtx *jrctx;

} SolCtx;

//---------------------------------------------------------------------------

void SolCtxSet(SolCtx *sc, FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx);

//---------------------------------------------------------------------------

// Augmented Lagrangian preconditioner context
typedef struct
{
	SolCtx      *sc;     // solution context
	BMat         P;      // block preconditiong matrix
	Mat          M;      // diagonal penalty matrix
	PetscScalar  pgamma; // penalty parameter
	KSP          ksp;    // internal Krylov solver context
	Vec          f, h;   // residual blocks
	Vec          u, p;   // solution blocks

} PCStokesALCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesALCreate(PCStokesALCtx *al, SolCtx *sc);

PetscErrorCode PCStokesALDestroy(PCStokesALCtx *al);

PetscErrorCode PCStokesALSetup(PCStokesALCtx *al);

PetscErrorCode PCStokesALApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// Block Factorization preconditioner context
typedef struct
{
	SolCtx   *sc;    // solution context
	BMat      P;     // block preconditiong matrix
	Mat       M;     // pressure Schur complement preconditioner
	IS        isv;   // velocity index set
	IS        isp;   // pressure index set
	PC        pc;    // block preconditioner context
	PCVelType vtype; // velocity solver type
	MGCtx     mg;    // multigrid context
	PC        pcmg;  // velocity preconditioner context

} PCStokesBFCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokesBFCtx *bf, SolCtx *sc);

PetscErrorCode PCStokesBFDestroy(PCStokesBFCtx *bf);

PetscErrorCode PCStokesBFSetup(PCStokesBFCtx *bf);

PetscErrorCode PCStokesBFApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
typedef struct
{
	SolCtx *sc; // solution context
	Mat     P;    // monolithic matrix
	Mat     M;    // inverse viscosity matrix (Jacobian compensation)
	MGCtx   mg;   // multigrid context
	PC      pc;   // coupled preconditioner context

} PCStokesMGCtx;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokesMGCtx *mg, SolCtx *sc);

PetscErrorCode PCStokesMGDestroy(PCStokesMGCtx *mg);

PetscErrorCode PCStokesMGSetup(PCStokesMGCtx *mg);

PetscErrorCode PCStokesMGApply(PC pc, Vec x, Vec y);

//---------------------------------------------------------------------------

// scatter block vectors to monolithic format & reverse
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode);

//---------------------------------------------------------------------------
#endif
