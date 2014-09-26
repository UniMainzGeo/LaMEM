//---------------------------------------------------------------------------
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#ifndef __matrix_h__
#define __matrix_h__
//---------------------------------------------------------------------------

PetscErrorCode PMatSetDiag(Mat P, PetscInt start, PetscInt ln, PetscScalar d, InsertMode mode);

PetscErrorCode PMatCreate(PetscInt m, PetscInt n, PetscInt d_nz,
	const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[], Mat *P);

PetscErrorCode PMatAssemble(Mat P, PetscInt numRows, const PetscInt rows[]);

//---------------------------------------------------------------------------

typedef struct
{
	Mat A;        // block matrix
	Mat Avv, Avp; // velocity sub-matrices
	Mat Apv, App; // pressure sub-matrices

} BMat;

//---------------------------------------------------------------------------

PetscErrorCode BMatCreate(BMat *bmat,
	PetscInt  lnv,       PetscInt  lnp,
	PetscInt *Avv_d_nnz, PetscInt *Avv_o_nnz,
	PetscInt *Avp_d_nnz, PetscInt *Avp_o_nnz,
	PetscInt *Apv_d_nnz, PetscInt *Apv_o_nnz);

PetscErrorCode BMatDestroy(BMat *bmat);

PetscErrorCode BMatClearSubMat(BMat *bmat);

PetscErrorCode BMatAssemble(BMat *bmat, BCCtx *bc);

//---------------------------------------------------------------------------

PetscErrorCode PMatCreateMonolithic(
	FDSTAG *fs,
	Mat    *P,
	Mat    *M);

PetscErrorCode PMatAssembleMonolithic(
	JacRes *jr,
	Mat     P,
	Mat     M);

//---------------------------------------------------------------------------

PetscErrorCode PMatCreateBlock(
	FDSTAG *fs,
	BMat   *P,
	Mat    *M);

PetscErrorCode PMatAssembleBlock(
	JacRes      *jr,
	BMat        *P,
	Mat          M,
	PetscScalar  pgamma);

//---------------------------------------------------------------------------

// apply two-point constraints on the ghost nodes
PetscErrorCode getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[]);

// constrain local matrix
PetscErrorCode constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[]);

// get velocity Schur complement from cell center local matrix, extract divergence & gradient sub-matrices
PetscErrorCode getVelSchurComp(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[], PetscScalar k);

// extract sub-matrices from stiffeness matrix
PetscErrorCode getSubMats(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[]);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// check locality / existence of the global DOF
#define CHECK_DOF(ind, start, num, nd, no) { if(ind != -1) { if(ind >= start && ind < start + num) nd++; else no++; } }

//---------------------------------------------------------------------------

#endif
