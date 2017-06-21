/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   matrix.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//.....................   PRECONDITIONING MATRICES   ........................
//---------------------------------------------------------------------------
#ifndef __matrix_h__
#define __matrix_h__
//---------------------------------------------------------------------------

struct JacRes;

// WARNING! Add MatSetNearNullSpace for all matrix types

//---------------------------------------------------------------------------

PetscErrorCode MatAIJCreate(PetscInt m, PetscInt n, PetscInt d_nz,
	const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[], Mat *P);

PetscErrorCode MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P);

PetscErrorCode MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag);

//---------------------------------------------------------------------------
// preconditioning matrix storage format
enum PMatType
{
	_MONOLITHIC_,
	_BLOCK_

} ;
//---------------------------------------------------------------------------

typedef struct _p_PMat *PMat;

typedef struct _p_PMat
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

	// get cell stiffness matrix
	void (*getStiffMat)(
		PetscScalar,  PetscScalar,
		PetscScalar*, PetscScalar*,
		PetscScalar,  PetscScalar, PetscScalar,
		PetscScalar,  PetscScalar, PetscScalar,
		PetscScalar,  PetscScalar, PetscScalar);

} p_PMat;

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

PetscErrorCode PMatMonoPicard(Mat J, Vec x, Vec y);

PetscErrorCode PMatMonoDestroy(PMat pm);

//---------------------------------------------------------------------------
//...........................   BLOCK MATRIX   ..............................
//---------------------------------------------------------------------------

struct PMatBlock
{
	Mat Avv, Avp; // velocity sub-matrices
	Mat Apv, App; // pressure sub-matrices
	Mat iS;       // inverse of Schur complement preconditioner

	Vec rv, rp;   // residual blocks
	Vec xv, xp;   // solution blocks
	Vec wv, wp;   // work vectors

};

//---------------------------------------------------------------------------

PetscErrorCode PMatBlockCreate(PMat pm);

PetscErrorCode PMatBlockAssemble(PMat pm);

PetscErrorCode PMatBlockPicardClean(Mat J, Vec x, Vec y);

PetscErrorCode PMatBlockPicardSchur(Mat J, Vec x, Vec y);

PetscErrorCode PMatBlockDestroy(PMat pm);

//---------------------------------------------------------------------------
// SERVICE FUNCTIONS
//---------------------------------------------------------------------------

// compute cell stiffness matrix with deviatoric projection
void getStiffMatDevProj(
	PetscScalar eta, PetscScalar diag,
	PetscScalar *v,  PetscScalar *cf,
	PetscScalar dx,  PetscScalar dy,   PetscScalar dz,
	PetscScalar fdx, PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx, PetscScalar bdy,  PetscScalar bdz);

// compute cell stiffness matrix without deviatoric projection
void getStiffMatClean(
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

//---------------------------------------------------------------------------

#endif
