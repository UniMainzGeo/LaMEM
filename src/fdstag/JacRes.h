//---------------------------------------------------------------------------
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#ifndef __JacRes_h__
#define __JacRes_h__

//---------------------------------------------------------------------------
// * replace setting time parameters consistently in the entire code

//---------------------------------------------------------------------------

// FDSTAG Jacobian and residual evaluation context
typedef struct
{
	// external handles
	FDSTAG  *fs;  // staggered-grid layout
	BCCtx   *bc;  // boundary condition context

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
	Vec ldxx, ldyy, ldzz, ldxy, ldxz, ldyz; // local (ghosted)
	Vec                   gdxy, gdxz, gdyz; // global
	// (ADVInterpMarkToEdge & ADVInterpFieldToMark is the only
	//  couple of functions where global vectors (gdxy, gdxz, gdyz) are used.
	//  Get a fuck rid of this ugly averaging between markers & edges!
	//  In ADVInterpFieldToMark it's easy.
	//  In ADVInterpMarkToEdge it's impossible because of assembly operation.
	//  Really really really need to switch to ghost marker approach!
	//  Also to get communication pattern independent of number of phases.

	// pressure & temperature
	Vec gp,  gT; // global
	Vec lp,  lT; // local (ghosted)

	// global continuity & energy residuals
	Vec gc, ge;

	// corner buffer
	Vec lbcor; // local (ghosted)

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
	PetscInt     numPhases;              // number phases
	Material_t   phases[max_num_phases]; // phase parameters
	PetscInt     numSoft;                // number material softening laws
	Soft_t       matSoft[max_num_soft];  // material softening law parameters
	MatParLim    matLim;                 // phase parameters limiters

	// parameters & controls
	Scaling     scal;        // scaling
	TSSol       ts;          // time-stepping parameters
	PetscScalar grav[SPDIM]; // global gravity components
	PetscScalar FSSA;        // density gradient penalty parameter
	//                          (a.k.a. free-surface-stabilization-algorithm)
	PetscScalar gtol;        // geometry tolerance

	PetscScalar pShift;      // pressure shift for plasticity model and output
	PetscBool   pShiftAct;   // pressure shift activation flag

} JacRes;
//---------------------------------------------------------------------------

PetscErrorCode JacResClear(JacRes *jr);

PetscErrorCode JacResSetFromOptions(JacRes *jr);

// create residual & Jacobian evaluation context
PetscErrorCode JacResCreate(
	JacRes   *jr,
	FDSTAG   *fs,
	BCCtx    *bc);

// destroy residual & Jacobian evaluation context
PetscErrorCode JacResDestroy(JacRes *jr);

// initialize and setup scaling object, perform scaling
PetscErrorCode JacResInitScale(JacRes *jr, UserCtx *usr);

// compute effective inverse elastic viscosity
PetscErrorCode JacResGetI2Gdt(JacRes *jr);

// get average pressure near the top surface
PetscErrorCode JacResGetPressShift(JacRes *jr);

// evaluate effective strain rate components in basic nodes
PetscErrorCode JacResGetEffStrainRate(JacRes *jr);

// compute components of vorticity vector
PetscErrorCode JacResGetVorticity(JacRes *jr);

// compute nonlinear residual vectors
PetscErrorCode JacResGetResidual(JacRes *jr);

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopySol(JacRes *jr, Vec x);

// copy residuals from local to global vectors, enforce boundary constraints
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f);

// copy momentum residuals from global to local vectors for output
PetscErrorCode JacResCopyMomentumRes(JacRes *jr, Vec f);

// copy continuity residuals from global to local vectors for output
PetscErrorCode JacResCopyContinuityRes(JacRes *jr, Vec f);

PetscErrorCode JacResViewRes(JacRes *jr);

PetscScalar JacResGetTime(JacRes *jr);

PetscInt JacResGetStep(JacRes *jr);

PetscErrorCode JacResGetCourantStep(JacRes *jr);

PetscErrorCode JacResInitTemp(JacRes *jr);

//---------------------------------------------------------------------------

// get maximum inverse time step on local domain
PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax);

//---------------------------------------------------------------------------

// initialize material parameter limits
PetscErrorCode SetMatParLim(MatParLim *matLim, UserCtx *usr);

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


