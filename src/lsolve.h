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
//......................   PRECONDITIONER ROUTINES   ........................
//---------------------------------------------------------------------------
#ifndef __lsolve_h__
#define __lsolve_h__
//---------------------------------------------------------------------------
// Picard operator storage type
enum PicardType
{
	_PICARD_ASSEMBLED_,
	_PICARD_MAT_FREE_
};
//---------------------------------------------------------------------------
// preconditioner type
enum PCStokesType
{
	_STOKES_BF_,  // block factorization
	_STOKES_MG_,  // coupled Galerkin multigrid
	_STOKES_USER_ // user-defined
};
//---------------------------------------------------------------------------
// block factorization type
enum PCBFType
{
	_BF_UPPER_,  // upper triangular factorization
	_BF_LOWER_   // lower triangular factorization
};
//---------------------------------------------------------------------------
// velocity block preconditioner type
enum PCVelType
{
	_VEL_MG_,   // Galerkin multigrid
	_VEL_USER_  // user-defined
};
//---------------------------------------------------------------------------
// Schur complement preconditioner type
enum PCSchurType
{
	_SCHUR_INV_ETA_, // inverse viscosity
	_SCHUR_WBFBT_    // wBFBT
};
//--------------------------------------------------------------------------

struct PCParam
{
	PicardType   ps_type; // Picard operator storage type
	PCStokesType pc_type; // preconditioner type
	PCBFType     bf_type; // block factorization type
	PCVelType    vs_type; // velocity solver type
	PCSchurType  sp_type; // Schur preconditioner type
	PetscScalar  pgamma;  // penalty parameter
};

//--------------------------------------------------------------------------

PetscErrorCode PCParamSetFromOptions(PCParam *p);

//--------------------------------------------------------------------------

struct PCData
{
	PCParam param; // parameters and options
	void    *data; // type-specific data
};

//---------------------------------------------------------------------------

PetscErrorCode PCDataCreate(PCData *pc, JacRes *jr, Mat J, Mat P);

PetscErrorCode PCDataDestroy(PCData *pc);

PetscErrorCode PCDataSetup(PCData *pc, JacRes *jr);

//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------

struct PCDataMG
{
	PCParam  *param; // parameters and options
	PMatMono  pm;    // monolithic matrix
	MatData   md;    // matrix assembly context
	MG        mg;    // coupled multigrid context
};

PetscErrorCode PCDataMGCreate(PCDataMG *pc, PCParam *param, JacRes *jr, Mat J, Mat P);

PetscErrorCode PCDataMGDestroy(PCDataMG *pc);

PetscErrorCode PCDataMGSetup(PCDataMG *pc, JacRes *jr);

PetscErrorCode PCDataMGApply(Mat P, Vec r, Vec x);

//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------

struct PCDataBF
{
	PCParam   *param; // parameters and options
	PMatBlock  pm;    // block matrix
	MatData    md;    // matrix assembly context
	MG         vmg;   // velocity multigrid context
	KSP        vksp;  // velocity solver
	KSP        pksp;  // pressure solver
};

PetscErrorCode PCDataBFCreate(PCDataBF *pc, PCParam *param, JacRes *jr, Mat J, Mat P);

PetscErrorCode PCDataBFDestroy(PCDataBF *pc);

PetscErrorCode PCDataBFSetup(PCDataBF *pc, JacRes *jr);

PetscErrorCode PCDataBFApply(Mat P, Vec r, Vec x);

PetscErrorCode PCDataBFBTApply(PCDataBF *pc, Vec x, Vec y);


//---------------------------------------------------------------------------
//............................. USER-DEFINED ................................
//---------------------------------------------------------------------------

struct PCDataUser
{
	PCParam  *param; // parameters and options
	PMatMono  pm;    // monolithic matrix
	MatData   md;    // matrix assembly context
	PC        pc;    // general preconditioner object
};

PetscErrorCode PCDataUserCreate(PCDataUser *pc, PCParam *param, JacRes *jr, Mat J, Mat P);

PetscErrorCode PCDataUserDestroy(PCDataUser *pc);

PetscErrorCode PCDataUserSetup(PCDataUser *pc, JacRes *jr);

PetscErrorCode PCDataUserApply(Mat P, Vec r, Vec x);

//--------------------------------------------------------------------------
#endif
