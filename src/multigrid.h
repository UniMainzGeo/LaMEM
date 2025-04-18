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

	FDSTAG    *fs;                       // staggered grid
	Vec        bcvx, bcvy, bcvz, bcp;    // restricted boundary condition vectors
	Vec        eta, etaxy, etaxz, etayz; // restricted viscosity vectors
	Mat        R, P;                     // restriction & prolongation operators (not set on finest grid)


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

PetscErrorCode MGLevelRestrictBC(MGLevel *lvl, MGLevel *fine, PetscBool no_restric_bc);

PetscErrorCode MGLevelSetupRestrict(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelSetupProlong(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelAllocRestrict(MGLevel *lvl, MGLevel *fine);

PetscErrorCode MGLevelAllocProlong(MGLevel *lvl, MGLevel *fine);

//---------------------------------------------------------------------------

// setup row of restriction matrix
void getRowRestrict(
		PetscScalar parent,
		PetscInt    n,
		PetscInt    idx[],
		PetscScalar bc[],
		PetscScalar v[],
		PetscScalar vs[]);

// setup row of prolongation matrix
void getRowProlong(
		PetscInt    parent,
		PetscScalar parent_bc,
		PetscInt    n,
		PetscScalar bc[],
		PetscScalar v[],
		PetscScalar vs[]);

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
