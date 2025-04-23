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

// Picard operator storage type
enum PicardType
{
	_PICARD_ASSEMBLED_,
	_PICARD_MAT_FREE_
};

//---------------------------------------------------------------------------

// preconditioner matrix type
enum PMatType
{
	_MONOLITHIC_,
	_BLOCK_
};

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

// user-defined preconditioner context
struct PCStokesUser
{
	PMatMono *pm; // monolithic matrix
	PC        pc; // general preconditioner object
};

//--------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
struct PCStokesMG
{
	PMatMono *pm; // monolithic matrix
	MG        mg; // coupled multigrid context
};

//--------------------------------------------------------------------------

// block factorization preconditioner context
struct PCStokesBF
{
	PMatBlock *pm; // block matrix
	MG         vmg;   // velocity multigrid context
	KSP        vksp;  // velocity solver
	KSP        pksp;  // pressure solver
};

//--------------------------------------------------------------------------

struct PCData
{
	PMatType     pm_type; // preconditioner matrix type
	PicardType   ps_type; // Picard operator storage type
	PCStokesType pc_type; // preconditioner type
	PCBFType     bf_type; // block factorization type
	PCVelType    vs_type; // velocity solver type
	PCSchurType  sp_type; // Schur preconditioner type
	PetscScalar  pgamma;  // penalty parameter

	MatData      *md;        // matrix assembly context
	void         *data;      // matrix

	PMatMono *Pm;
	PMatMono *Pb;


	/*

		PetscErrorCode (*Create)  (PMat pm);
		PetscErrorCode (*Assemble)(PMat pm);
		PetscErrorCode (*Destroy) (PMat pm);
		PetscErrorCode (*Picard)  (Mat J, Vec x, Vec y);
		PetscErrorCode (*Apply)   (Mat P, Vec x, Vec y);
		*/
};

//---------------------------------------------------------------------------

PetscErrorCode PCDataSetFromOptions(PCData *pc);

PetscErrorCode PCDataCreate(PCData *pc, JacRes *jr);

PetscErrorCode PCDataDestroy(PCData *pc);

PetscErrorCode PCDataSetup(PCData *pc);

PetscErrorCode PCDataSetApply (PCData *pc, Mat P);

PetscErrorCode PCDataSetPicard(PCData *pc, Mat J);


/*

PetscErrorCode PCStokesMGCreate(PCStokesMG *pc);

PetscErrorCode PCStokesMGDestroy(PCStokesMG *pc);

PetscErrorCode PCStokesMGSetup(PCStokesMG *pc);

PetscErrorCode PCStokesMGApply(Mat P, Vec x, Vec y);

//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFDestroy(PCStokes pc);

PetscErrorCode PCStokesBFSetup(PCStokes pc);

PetscErrorCode PCStokesBFApply(Mat JP, Vec x, Vec y);


//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

PetscErrorCode PCStokesUserCreate(PCStokes pc);

PetscErrorCode PCStokesUserDestroy(PCStokes pc);

PetscErrorCode PCStokesUserSetup(PCStokes pc);

PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

*/

#endif
