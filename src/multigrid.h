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
 **    filename:   multigrid.h
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
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#ifndef __multigrid_h__
#define __multigrid_h__

//---------------------------------------------------------------------------

struct BCCtx;
struct FDSTAG;
struct JacRes;

//---------------------------------------------------------------------------

// Galerkin multigrid level data structure

struct MGLevel
{
	// Constrained DOF stores parent DOF index in the boundary condition vector.
	// Parent DOF index is the only nonzero that is set in the row of R-matrix
	// and column of P-matrix to impose the constraints in a coarse grid operator
	// automatically. The finest grid uses standard boundary condition vectors.

	DM        DA_CEN;                // central points array
	DM        DA_X, DA_Y, DA_Z;      // face points arrays
	DOFIndex  dof;                   // indexing vectors
	Vec       bcvx, bcvy, bcvz, bcp; // restricted boundary condition vectors
	Vec       eta, etax, etay, etaz; // viscosity vectors
	Mat       R, P;                  // restriction & prolongation operators (not set on finest grid)


	// ******** fine level ************
	//     |                   ^
	//     R-matrix            |
	//     |                   P-matrix
	//     v                   |
	// ******** this level ************

} ;

//---------------------------------------------------------------------------

PetscErrorCode MGLevelCreate(MGLevel *lvl, MGLevel *fine, FDSTAG *fs, BCCtx *bc);

PetscErrorCode MGLevelDestroy(MGLevel *lvl);

PetscErrorCode MGLevelInitEta(MGLevel *lvl, JacRes *jr);

PetscErrorCode MGLevelRestrictEta(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelAverageEta(MGLevel *lvl);

PetscErrorCode MGLevelRestrictBC(MGLevel *lvl, MGLevel *fine, PetscBool restric_bc);

PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelSetupProlong(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelAllocRestrict(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelAllocProlong(MGLevel *lvl, MGLevel *fine);

//---------------------------------------------------------------------------

// setup row of restriction matrix
void getRowRestrict(PetscBool scale,
	PetscScalar parent, PetscInt n, PetscInt idx[], PetscScalar bc[],
	PetscScalar v[], PetscScalar vs[], PetscScalar eta_fine[], PetscScalar eta_crs);

// setup row of prolongation matrix
void getRowProlong(PetscBool scale,
	PetscInt parent, PetscScalar parent_bc, PetscInt n, PetscScalar bc[],
	PetscScalar v[], PetscScalar vs[], PetscScalar eta_crs[], PetscScalar eta_fine);

//---------------------------------------------------------------------------

struct MG
{
	// PETSc level numbering (inverse w.r.t. coarsening sequence):
	// 0   - coarse grid
	// n-1 - fine grid
	// R & P matrices connect with coarser level (i.e. not set on coarsest grid).
	// Coarsening step yields coarse grid operator. Own operator is prescribed.

	// LaMEM level numbering (natural w.r.t. coarsening sequence):
	// 0   - fine grid
	// n-1 - coarse grid
	// R & P matrices connect with finer level (i.e. not set on finest grid).
	// Coarsening step yields own operator. Fine level operator is prescribed.

	PetscInt  nlvl; // number of levels
	MGLevel  *lvls; // multigrid levles

	PC        pc;   // internal preconditioner context
	JacRes   *jr;   // finest level context

	PetscBool crs_setup;     // coarse solver setup flag
	PetscBool no_restric_bc; // boundary constraint restriction deactivation flag

};

//---------------------------------------------------------------------------

PetscErrorCode MGCreate(MG *mg, JacRes *jr);

PetscErrorCode MGDestroy(MG *mg);

PetscErrorCode MGSetupCoarse(MG *mg, Mat A);

PetscErrorCode MGSetup(MG *mg, Mat A);

PetscErrorCode MGApply(PC pc, Vec x, Vec y);

PetscErrorCode MGDumpMat(MG *mg);

PetscErrorCode MGGetNumLevels(MG *mg);

//---------------------------------------------------------------------------
#endif
