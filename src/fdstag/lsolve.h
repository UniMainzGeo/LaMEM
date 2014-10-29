//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __lsolve_h__
#define __lsolve_h__
//---------------------------------------------------------------------------
// Stokes preconditioner type
typedef enum
{
	_STOKES_BF_,  // block factorization
	_STOKES_MG_,  // Galerkin multigrid
	_STOKES_USER_ // user-defined

} PCStokesType;

//---------------------------------------------------------------------------

// velocity block preconditioner type (bf only)
typedef enum
{
	_VEL_MG_,  // Galerkin multigrid
	_VEL_USER_ // user-defined

} PCVelType;
//---------------------------------------------------------------------------

typedef struct _p_PCStokes *PCStokes;

typedef struct _p_PCStokes
{
	PCStokesType  type;
	PMat          pm;   // preconditioner matrix
	void         *data; // type-specific context

	// operations
	PetscErrorCode (*Create)  (PCStokes pc);
	PetscErrorCode (*Setup)   (PCStokes pc);
	PetscErrorCode (*Destroy) (PCStokes pc);
	PetscErrorCode (*Apply)   (Mat P, Vec x, Vec y);

} p_PCStokes;

// PCStokes - pointer to an opaque structure (to be used in declarations)
// sizeof(p_PCStokes) - size of the opaque structure

//---------------------------------------------------------------------------

PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm);

PetscErrorCode PCStokesSetFromOptions(PCStokes pc);

PetscErrorCode PCStokesSetup(PCStokes pc);

PetscErrorCode PCStokesDestroy(PCStokes pc);

//---------------------------------------------------------------------------

// Block Factorization preconditioner context
typedef struct
{
	PCVelType  vtype; // velocity solver type
	KSP        vksp;  // velocity solver
	MG         vmg;   // velocity multigrid context

} PCStokesBF;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc);

PetscErrorCode PCStokesBFDestroy(PCStokes pc);

PetscErrorCode PCStokesBFSetup(PCStokes pc);

PetscErrorCode PCStokesBFApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
typedef struct
{
	MG mg; // coupled multigrid context

} PCStokesMG;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokes pc);

PetscErrorCode PCStokesMGDestroy(PCStokes pc);

PetscErrorCode PCStokesMGSetup(PCStokes pc);

PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// User-defined
typedef struct
{
	PC pc; // general preconditioner object

} PCStokesUser;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesUserCreate(PCStokes pc);

PetscErrorCode PCStokesUserDestroy(PCStokes pc);

PetscErrorCode PCStokesUserSetup(PCStokes pc);

PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------
#endif
