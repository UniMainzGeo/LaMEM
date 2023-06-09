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

// Stokes preconditioner type
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
struct PCStokesBF
{
	PCVelType vtype; // velocity solver type
	KSP       vksp;  // velocity solver
	MG        vmg;   // velocity multigrid context
	PCBFType  type;  // factorization type

};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc);

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
	PC pc;       // general preconditioner object
	IS isv, isp; // velocity and pressure index sets

} ;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesUserCreate(PCStokes pc);

PetscErrorCode PCStokesUserAttachIS(PCStokes pc);

PetscErrorCode PCStokesUserDestroy(PCStokes pc);

PetscErrorCode PCStokesUserSetup(PCStokes pc);

PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

#endif
