//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __lsolve_h__
#define __lsolve_h__
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

typedef struct _p_PCStokes *PCStokes;

typedef struct _p_PCStokes
{
	JacRes *jr;
	void   *data; // type-specific context

	// internal operations
	PetscErrorCode (*Create) (PCStokes pc);
	PetscErrorCode (*Setup)  (PCStokes pc);
	PetscErrorCode (*Destroy)(PCStokes pc);
	PetscErrorCode (*Apply)  (Mat P, Vec x, Vec y);
	PetscErrorCode (*Picard) (Mat J, Vec x, Vec y);

}p_PCStokes;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesCreate(
	PCStokes      pc,
	JacRes       *jr,
	PCStokesType  ptype);

PetscErrorCode PCStokesSetup(PCStokes pc);

PetscErrorCode PCStokesDestroy(PCStokes pc);

//---------------------------------------------------------------------------

// Augmented Lagrangian preconditioner context
typedef struct
{
	BMat         P;      // block preconditiong matrix
	Mat          M;      // diagonal penalty matrix
	PetscScalar  pgamma; // penalty parameter
	KSP          ksp;    // internal Krylov solver context
	Vec          f, h;   // residual blocks
	Vec          u, p;   // solution blocks

} PCStokesAL;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesALCreate(PCStokes pc);

PetscErrorCode PCStokesALDestroy(PCStokes pc);

PetscErrorCode PCStokesALSetup(PCStokes pc);

PetscErrorCode PCStokesALApply(Mat P, Vec x, Vec y);

PetscErrorCode PCStokesALPicard(Mat J, Vec x, Vec y);

//---------------------------------------------------------------------------

// Block Factorization preconditioner context
typedef struct
{
	BMat       P;     // block preconditiong matrix
	Mat        M;     // pressure Schur complement preconditioner
	IS         isv;   // velocity index set
	IS         isp;   // pressure index set
	PC         pc;    // internal preconditioner context
	PCVelType  vtype; // multigrid velocity solver flag
	MGCtx      vctx;  // velocity multigrid context
	PC         vpc;   // velocity preconditioner context

} PCStokesBF;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFDestroy(PCStokes pc);

PetscErrorCode PCStokesBFSetup(PCStokes pc);

PetscErrorCode PCStokesBFApply(Mat P, Vec x, Vec y);

PetscErrorCode PCStokesBFPicard(Mat J, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
typedef struct
{
	Mat   P;   // monolithic matrix
	Mat   M;   // inverse viscosity matrix (Jacobian compensation)
	PC    pc;  // internal preconditioner context
	MGCtx ctx; // multigrid context
	Vec   w;   // work vector for computing Jacobian action

} PCStokesMG;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokes pc);

PetscErrorCode PCStokesMGDestroy(PCStokes pc);

PetscErrorCode PCStokesMGSetup(PCStokes pc);

PetscErrorCode PCStokesMGApply(Mat P, Vec x, Vec y);

PetscErrorCode PCStokesMGPicard(Mat J, Vec x, Vec y);

//---------------------------------------------------------------------------

// scatter block vectors to monolithic format & reverse
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode);

//---------------------------------------------------------------------------
#endif
