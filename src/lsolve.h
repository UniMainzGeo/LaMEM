/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __lsolve_h__
#define __lsolve_h__
//---------------------------------------------------------------------------

// Stokes preconditioner type
enum PCStokesType
{
	_STOKES_BF_,  // block factorization
	_STOKES_MG_,  // Galerkin multigrid
	_STOKES_USER_ // user-defined

};

//---------------------------------------------------------------------------

// block factorization type
enum PCBFType
{
	_UPPER_,  // upper triangular factorization
	_LOWER_   // lower triangular factorization

};

//---------------------------------------------------------------------------

// velocity block preconditioner type (bf only)
enum PCVelType
{
	_VEL_MG_,  // Galerkin multigrid
	_VEL_USER_ // user-defined

};

//---------------------------------------------------------------------------

typedef struct _p_PCStokes *PCStokes;

typedef struct _p_PCStokes
{
	PCStokesType pctype;     // preconditioner type
	PCBFType     ftype;      // factorization type
	PCVelType    vtype;      // velocity solver type
	PetscScalar  pgamma;     // penalty parameter
	PetscInt     buildwBFBT; // flag to build wbfbt matrix
	PetscInt     buildCvv;   // flag to build clean velocity sub-matix
	PMat         pm;         // preconditioner matrix
	MatData      *md;        // assembly context
	void         *data;      // type-specific context

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

// Block factorization preconditioner context
struct PCStokesBF
{
	PCBFType    ftype; // factorization type
	PCVelType   vtype; // velocity solver type
	KSP         vksp;  // velocity solver
	MG          vmg;   // velocity multigrid context
	KSP 	    pksp;  // pressure solver
};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFDestroy(PCStokes pc);

PetscErrorCode PCStokesBFSetup(PCStokes pc);

PetscErrorCode PCStokesBFApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
struct PCStokesMG
{
	MG mg; // coupled multigrid context
};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokes pc);

PetscErrorCode PCStokesMGDestroy(PCStokes pc);

PetscErrorCode PCStokesMGSetup(PCStokes pc);

PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// User-defined
struct PCStokesUser
{
	PC pc; // general preconditioner object
};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesUserCreate(PCStokes pc);

PetscErrorCode PCStokesUserDestroy(PCStokes pc);

PetscErrorCode PCStokesUserSetup(PCStokes pc);

PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

#endif
