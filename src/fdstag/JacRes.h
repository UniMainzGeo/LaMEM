//---------------------------------------------------------------------------
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#ifndef __JacRes_h__
#define __JacRes_h__
//---------------------------------------------------------------------------
// FDSTAG Jacobian and residual evaluation context
typedef struct
{
	// coupled solution & residual vectors
	Vec gsol, gres; // global
//	Vec lsol, lres; // local (ghosted)

	// velocity	components
	Vec gvx,  gvy, gvz;  // global
	Vec lvx,  lvy, lvz;  // local (ghosted)

	// momentum residual components
	Vec gfx,  gfy, gfz;  // global
	Vec lfx,  lfy, lfz;  // local (ghosted)

	// strain-rate components
	Vec gdxx, gdyy, gdzz, gdxy, gdxz, gdyz; // global
	Vec ldxx, ldyy, ldzz, ldxy, ldxz, ldyz; // local (ghosted)

	// pressure & temperature
	Vec gp,  gT; // global
	Vec lp,  lT; // local (ghosted)

	// global continuity & energy residuals
	Vec gc, ge;

	// heat conductivity on cells, and shear heating terms on edges
//	Vec gk, ghxy, ghxz, ghyz; // global
//	Vec lk, lhxy, lhxz, lhyz; // local (ghosted)

	VecScatter g2lctx; // global to local scatter context

	// solution variables
	SolVarCell  *svCell;   // cell centers
	SolVarEdge  *svXYEdge; // XY edges
	SolVarEdge  *svXZEdge; // XZ edges
	SolVarEdge  *svYZEdge; // YZ edges
	PetscScalar *svBuff;   // storage for phRat

	// phase parameters
	PetscInt     numPhases; // number phases
	Material_t  *phases;    // phase parameters
	PetscInt     numSoft;   // number material softening laws
	Soft_t      *matSoft;   // material softening law parameters
	MatParLim    matLim;    // phase parameters limiters

	// global parameters
	PetscScalar dt;          // time step
	PetscScalar grav[SPDIM]; // global gravity components

	// scaling
	Scaling scal;

} JacResCtx;
//---------------------------------------------------------------------------

// create residual & Jacobian evaluation context
PetscErrorCode FDSTAGCreateJacResCtx(
	FDSTAG    *fs,
	JacResCtx *jrctx,
	PetscInt   numPhases,
	PetscInt   numSoft);

// destroy residual & Jacobian evaluation context
PetscErrorCode FDSTAGDestroyJacResCtx(JacResCtx *jrctx);

// compute effective inverse elastic viscosity
PetscErrorCode FDSTAGetI2Gdt(FDSTAG *fs, JacResCtx *jrctx);

// evaluate effective strain rate components in basic nodes
PetscErrorCode FDSTAGetEffStrainRate(FDSTAG *fs, JacResCtx *jrctx);

// compute nonlinear residual vectors
PetscErrorCode FDSTAGetResidual(FDSTAG *fs, JacResCtx *jrctx);

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode FDSTAGCopySol(FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx, Vec x);

// copy residuals to global vector
PetscErrorCode FDSTAGCopyRes(FDSTAG *fs, BCCtx *bc, JacResCtx *jrctx, Vec f);

//---------------------------------------------------------------------------

// initialize material parameter limits
PetscErrorCode SetMatParLim(MatParLim *matLim, UserContext *usr);

//---------------------------------------------------------------------------


/*
PetscErrorCode FDSTAGScatterSol(FDSTAG *fs, JacResCtx *jrctx);

PetscErrorCode FDSTAGScatterRes(FDSTAG *fs, JacResCtx *jrctx);

PetscErrorCode FDSTAGCreateScatter(FDSTAG *fs, JacResCtx *jrctx);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

// check existence of the global DOF
#define CHECK_DOF_INTERNAL(ind, start, num, gidx, lidx) { if(ind != -1) { gidx[num] = ind; lidx[num] = start; num++; } }
*/
//---------------------------------------------------------------------------
#endif


