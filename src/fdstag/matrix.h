//---------------------------------------------------------------------------
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#ifndef __matrix_h__
#define __matrix_h__
//---------------------------------------------------------------------------

PetscErrorCode PMatCreate(PetscInt m, PetscInt n, PetscInt d_nz,
	const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[], Mat *P);

PetscErrorCode PMatAssemble(Mat P, PetscInt numRows, const PetscInt rows[]);

//---------------------------------------------------------------------------

PetscErrorCode PMatCreateMonolithic(FDSTAG  *fs, Mat *A, Mat *InvEta);

PetscErrorCode PMatAssembleMonolithic(FDSTAG  *fs, BCCtx *bc, JacResCtx *jrctx, Mat A, Mat InvEta, PetscBool precond);

//---------------------------------------------------------------------------

PetscErrorCode PMatCreateBlock(FDSTAG *fs, BlockMat *bmat);

PetscErrorCode PMatAssembleBlock(FDSTAG  *fs, BCCtx *bc, JacResCtx *jrctx, BlockMat *bmat, PetscScalar pgamma);

//---------------------------------------------------------------------------

// apply two-point constraints on the ghost nodes
PetscErrorCode getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[]);

// constrain local matrix
PetscErrorCode constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[]);

// get velocity Schur complement from cell center local matrix, extract divergence & gradient sub-matrices
PetscErrorCode getVelSchurComp(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[], PetscScalar k);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// check locality / existence of the global DOF
#define CHECK_DOF(ind, start, num, nd, no) { if(ind != -1) { if(ind >= start && ind < start + num) nd++; else no++; } }

//---------------------------------------------------------------------------

#endif
