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
static inline PetscInt constrEdgeNode(
	PetscInt    ix[],
	PetscInt    iy[],
	PetscScalar bx[],
	PetscScalar by[],
	PetscScalar vx[],
	PetscScalar vy[],
	PetscScalar dx,
	PetscScalar dy)
{
	PetscInt    cx = -1, cy = -1, px = -1, py = -1;
	PetscScalar epsb = 0.0, eps = 0.0, cfx = 0.0, cfy = 0.0;

	// epsb - is a boundary value for the strain rate component.
	// For the free surface it should be zero.
	// Alternatively, if a nonzero boundary stress is required,
	// this strain rate should be computed from the constitutive model
	// (probably nonlinear).
	// Currently only the free surface is implemented, so epsb is set to zero.

	// eps - is a cumulative strain rate from internal and Dirichlet ghost nodes.
	// It is used for enforcing the Neumann constraints.

	// cfx, cfy - are the strain-rate terms pre-multipliers used in the expansion
	// of Neumann constraints. They assume values of -1 or 1, depending on which
	// node along each direction is constrained, first or second (correspondingly).

	// determine primary internal nodes & constrained ghost nodes
	if(ix[0] == -1) { cx = 0; px = 1; cfx = -1.0; }
	if(ix[1] == -1) { cx = 1; px = 0; cfx =  1.0; }
	if(iy[0] == -1) { cy = 0; py = 1; cfy = -1.0; }
	if(iy[1] == -1) { cy = 1; py = 0; cfy =  1.0; }

	// sort out the internal nodes
	if(cx == -1 && cy == -1) return 0;

	// compute strain rates from internal nodes
	if(cx == -1) eps += (vx[1] - vx[0])/dy/2.0;
	if(cy == -1) eps += (vy[1] - vy[0])/dx/2.0;

	// expand Dirichlet constraints, update strain rates from Dirichlet nodes
	if(cx != -1 && bx[cx] != DBL_MAX) { vx[cx] = 2.0*bx[cx] - vx[px]; eps += (vx[1] - vx[0])/dy/2.0; }
	if(cy != -1 && by[cy] != DBL_MAX) { vy[cy] = 2.0*by[cy] - vy[py]; eps += (vy[1] - vy[0])/dx/2.0; }

	// expand Neumann constraints
	if(cx != -1 && bx[cx] == DBL_MAX) { vx[cx] = vx[px] + 2.0*dy*cfx*(epsb - eps); }
	if(cy != -1 && by[cy] == DBL_MAX) { vy[cy] = vy[py] + 2.0*dx*cfy*(epsb - eps); }

	return 1;
}

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


