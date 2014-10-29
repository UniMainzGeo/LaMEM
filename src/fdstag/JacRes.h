//---------------------------------------------------------------------------
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#ifndef __JacRes_h__
#define __JacRes_h__
//---------------------------------------------------------------------------
// FDSTAG Jacobian and residual evaluation context
typedef struct
{
	// external handles
	FDSTAG *fs;  // staggered-grid layout
	BCCtx  *cbc; // boundary condition context (coupled)
	BCCtx  *ubc; // boundary condition context (uncoupled)

	// coupled solution & residual vectors
	Vec gsol, gres; // global
//	Vec lsol, lres; // local (ghosted)

	// velocity	components
	Vec gvx,  gvy, gvz;  // global
	Vec lvx,  lvy, lvz;  // local (ghosted)

	// momentum residual components
	Vec gfx,  gfy, gfz;  // global
	Vec lfx,  lfy, lfz;  // local (ghosted)

	// strain-rate components (also used as buffer vectors)
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

//	VecScatter g2lctx; // global to local scatter context

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

} JacRes;
//---------------------------------------------------------------------------

// create residual & Jacobian evaluation context
PetscErrorCode JacResCreate(
	JacRes   *jr,
	FDSTAG   *fs,
	BCCtx    *cbc,
	BCCtx    *ubc,
	PetscInt  numPhases,
	PetscInt  numSoft);

// destroy residual & Jacobian evaluation context
PetscErrorCode JacResDestroy(JacRes *jr);

// compute effective inverse elastic viscosity
PetscErrorCode JacResGetI2Gdt(JacRes *jr);

// evaluate effective strain rate components in basic nodes
PetscErrorCode JacResGetEffStrainRate(JacRes *jr);

// compute components of vorticity vector
PetscErrorCode JacResGetVorticity(JacRes *jr);

// compute nonlinear residual vectors
PetscErrorCode JacResGetResidual(JacRes *jr);

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopySol(JacRes *jr, Vec x);

// copy residuals to global vector
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f);

PetscErrorCode JacResViewRes(JacRes *jr);

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


