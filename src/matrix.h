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
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#ifndef __matrix_h__
#define __matrix_h__
//---------------------------------------------------------------------------

struct JacRes;
struct DOFIndex;

//---------------------------------------------------------------------------

PetscErrorCode MatAIJCreate(PetscInt m, PetscInt n, PetscInt d_nz,
	const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[], Mat *P);

PetscErrorCode MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P);

PetscErrorCode MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag);

PetscErrorCode MatAIJSetNullSpace(Mat P, DOFIndex *dof);

//---------------------------------------------------------------------------

// preconditioning matrix storage format
enum PMatType
{
	_MONOLITHIC_,
	_BLOCK_
};

//---------------------------------------------------------------------------

typedef struct p_PMat *PMat;

struct p_PMat
{
	JacRes     *jr;     // assembly context
	void       *data;   // type-specific context
	PMatType    type;   // matrix type
	PetscScalar pgamma; // penalty parameter

	// operations
	PetscErrorCode (*Create)  (PMat pm);
	PetscErrorCode (*Assemble)(PMat pm);
	PetscErrorCode (*Destroy) (PMat pm);
	PetscErrorCode (*Picard)  (Mat J, Vec x, Vec y);
};

// PMat - pointer to an opaque structure (to be used in declarations)
// sizeof(p_PMat) - size of the opaque structure

//---------------------------------------------------------------------------

PetscErrorCode PMatSetFromOptions(PMat pm);

PetscErrorCode PMatCreate(PMat *p_pm, JacRes *jr);

PetscErrorCode PMatAssemble(PMat pm);

PetscErrorCode PMatDestroy(PMat pm);

//---------------------------------------------------------------------------
//.........................   MONOLITHIC MATRIX   ...........................
//---------------------------------------------------------------------------

struct PMatMono
{
	Mat A; // monolithic matrix
	Mat M; // penalty terms compensation matrix
	Vec w; // work vector for computing Jacobian action
};

PetscErrorCode PMatMonoCreate(PMat pm);

PetscErrorCode PMatMonoAssemble(PMat pm);

PetscErrorCode PMatMonoDestroy(PMat pm);

PetscErrorCode PMatMonoPicard(Mat J, Vec x, Vec y);

//---------------------------------------------------------------------------
//...........................   BLOCK MATRIX   ..............................
//---------------------------------------------------------------------------

struct PMatBlock
{
	Mat Avv, Avp;  // velocity sub-matrices
	Mat Apv, App;  // pressure sub-matrices
	Mat iS;        // inverse of pressure Schur complement preconditioner
	Mat Cvv;       // clean velocity sub-matix

	Vec rv, rp;   // residual blocks
	Vec xv, xp;   // solution blocks
	Vec wv, wp;   // work vectors

	// wBFBT data
	PetscBool wbfbt; // wbfbt preconditioner flag
	DM        DA_P;  // cell-based grid
	Mat       K;     // Schur complement preconditioner matrix
	Mat       C;     // diagonal viscosity weighting matrix
	Vec       w;     // working vector in velocity space
};

//---------------------------------------------------------------------------

PetscErrorCode PMatBlockSetFromOptions(PMat pm);

PetscErrorCode PMatBlockCreate(PMat pm);

PetscErrorCode PMatBlockAssemble(PMat pm);

PetscErrorCode PMatBlockDestroy(PMat pm);

PetscErrorCode PMatBlockPicard(Mat J, Vec x, Vec y);

//---------------------------------------------------------------------------
// SERVICE FUNCTIONS
//---------------------------------------------------------------------------

// compute cell stiffness matrix with deviatoric projection
void getStiffMat(
	PetscScalar eta, PetscScalar diag,
	PetscScalar *v,  PetscScalar *cf,
	PetscScalar dx,  PetscScalar dy,   PetscScalar dz,
	PetscScalar fdx, PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx, PetscScalar bdy,  PetscScalar bdz);

void addDensGradStabil(
	PetscScalar fssa, PetscScalar *v,
	PetscScalar rho,  PetscScalar dt,   PetscScalar *grav,
	PetscScalar fdx,  PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx,  PetscScalar bdy,  PetscScalar bdz);

// compute velocity Schur complement
void getVelSchur(PetscScalar v[], PetscScalar d[], PetscScalar g[]);

// extract sub-matrices from stiffness matrix
void getSubMat(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[]);

// apply two-point constraints on the ghost nodes
void getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[]);

// constrain local matrix
void constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[]);

//---------------------------------------------------------------------------

// scatter block vectors to monolithic format & reverse
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// check locality / existence of the global DOF
#define CHECK_DOF(ind, start, num, nd, no) { if(ind != -1) { if(ind >= start && ind < start + num) nd++; else no++; } }

// set pressure two-point constraint
#define SET_PRES_TPC(bc, i, j, k, ind, lim, cf) { cf = 1.0; if(ind == lim && bc[k][j][i] != DBL_MAX) cf = 2.0; }

// compute scaled shear stress stencil for constrained internal velocity case
#define RESCALE_STENCIL(rescal, d, df, db, cf, cb, dr) { dr = 0.0; if(rescal) { if(cf == DBL_MAX) { dr += df; } if(cb == DBL_MAX) { dr += db; } } if(dr) { d = dr/2.0; } }

//---------------------------------------------------------------------------

#endif
