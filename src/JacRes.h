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
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#ifndef __JacRes_h__
#define __JacRes_h__
//---------------------------------------------------------------------------

struct FB;
struct Scaling;
struct TSSol;
struct FDSTAG;
struct FreeSurf;
struct BCCtx;
struct DBMat;
struct DBPropDike;
struct Dike;
struct Tensor2RN;
struct PData;
struct AdvCtx;

//---------------------------------------------------------------------------
//.....................   Deviatoric solution variables   ...................
//---------------------------------------------------------------------------

struct SolVarDev
{
	PetscScalar  eta;    // total effective viscosity
	PetscScalar  eta_st; // stabilization viscosity
	PetscScalar  I2Gdt;  // inverse elastic parameter (1/2G/dt)
	PetscScalar  Hr;     // shear heating term contribution
	PetscScalar  APS;    // accumulated plastic strain
	PetscScalar  PSR;    // plastic strain-rate contribution

};

//---------------------------------------------------------------------------
//.....................   Volumetric solution variables   ...................
//---------------------------------------------------------------------------

struct SolVarBulk
{
	PetscScalar  theta;  // volumetric strain rate
	PetscScalar  rho;    // strain- & temperature-dependent density
	PetscScalar  IKdt;   // inverse bulk elastic parameter (1/K/dt)
	PetscScalar  alpha;  // effective thermal expansion
	PetscScalar  Tn;     // history temperature
	PetscScalar  pn;     // history pressure
	PetscScalar  rho_pf; // fluid density from phase diagram
	PetscScalar  mf;     // melt fraction from phase diagram
	PetscScalar  phi;    // PSD angle
	PetscScalar  Ha ;    // Adiabatic heating
	PetscScalar  cond ;  // conductivity

};

//---------------------------------------------------------------------------
//........................   Cell solution variables   ......................
//---------------------------------------------------------------------------

struct SolVarCell
{
	SolVarDev    svDev;         // deviatoric variables
	SolVarBulk   svBulk;        // volumetric variables
	PetscScalar  sxx, syy, szz; // deviatoric stress
	PetscScalar  hxx, hyy, hzz; // history stress (elastic)
	PetscScalar  dxx, dyy, dzz; // total deviatoric strain rate
	PetscScalar *phRat;         // phase ratios in the control volume
	PetscInt     FreeSurf;      // indicates whether the control volume contains the internal free surface
	PetscScalar  U[3];          // total displacement
	PetscScalar  ATS;           // accumulated total strain
	PetscScalar  eta_cr;        // creep viscosity
	PetscScalar  DIIdif;        // relative diffusion creep strain rate
	PetscScalar  DIIdis;        // relative dislocation creep strain rate
	PetscScalar  DIIprl;        // relative Peierls creep strain rate
	PetscScalar  DIIfk;         // relative Frank-Kamenetzky creep strain rate
	PetscScalar  DIIpl;         // relative plastic strain rate
	PetscScalar  yield;         // average yield stress in control volume

};

//---------------------------------------------------------------------------
//........................   Edge solution variables   ......................
//---------------------------------------------------------------------------

struct SolVarEdge
{
	SolVarDev    svDev; // deviatoric variables
	PetscScalar  s;     // xy, xz, yz deviatoric stress components
	PetscScalar  h;     // xy, xz, yz history stress components (elastic)
	PetscScalar  d;     // xy, xz, yz total deviatoric strain rate components
	PetscScalar  ws;    // normalization for distance-dependent interpolation
	PetscScalar *phRat; // phase ratios in the control volume

};

//---------------------------------------------------------------------------
//...................   Runtime parameters and controls .....................
//---------------------------------------------------------------------------

// Ground water level type
enum GWLevelType
{
	_GW_NONE_,   // don't compute pore pressure
	_GW_TOP_,    // top of the domain
	_GW_SURF_,   // free surface
	_GW_LEVEL_   // fixed level

};

// solution access mode
enum AccessMode
{
	_interp_,    // assigned edges and corners for interpolation
	_no_interp_, // no interpolation is required
};

struct Controls
{
	PetscScalar grav[3];       // global gravity components
	PetscScalar FSSA;          // free surface stabilization parameter [0 - 1]
	PetscScalar shearHeatEff;  // shear heating efficiency parameter [0 - 1]
	PetscScalar biot;          // Biot pressure parameter [0 - 1]

	PetscScalar AdiabHeat;      // Adiabatic Heating efficiency
	PetscInt    actTemp;        // temperature diffusion activation flag
	PetscInt    actExp;         // thermal expansion activation flag
	PetscInt    actSteadyTemp;  // steady-state temperature initial guess flag
	PetscScalar steadyTempStep; // time for (quasi-)steady-state temperature initial guess
	PetscInt    steadyNumStep;  // number of steps for (quasi-)steady-state temperature initial guess
	PetscInt    actHeatRech;    // heat recharge setting
	PetscInt    initLithPres;   // set initial pressure to lithostatic pressure
	PetscInt    initGuess;      // initial guess activation flag
	PetscInt    pLithoVisc;     // use lithostatic pressure for creep laws
	PetscInt    pLithoPlast;    // use lithostatic pressure for plasticity
	PetscInt    pLimPlast;      // limit pressure at first iteration for plasticity
	PetscScalar pShift;         // shift the pressure by a constant value while evaluating plasticity & for output
	PetscInt    pShiftAct;      // pressure shift activation flag (zero pressure in the top cell layer)
	PetscInt    printNorms;     // priny norms of velocity/pressure/temperature?

	PetscScalar eta_min;        // minimum viscosity
	PetscScalar eta_max;        // maximum viscosity
	PetscScalar eta_ref;        // reference viscosity (initial guess)
	PetscScalar TRef;           // reference temperature
	PetscScalar Rugc;           // universal gas constant
	PetscScalar minCh;          // minimum cohesion
	PetscScalar minFr;          // minimum friction
	PetscScalar tauUlt;         // ultimate yield stress

	PetscScalar rho_fluid;      // fluid density
	GWLevelType gwType;         // type of ground water level (none, top, surf, level)
	PetscScalar gwLevel;        // fixed ground water level

	PetscInt    getPermea;      // effective permeability computation activation flag
	PetscInt    rescal;         // stencil rescaling flag (for interval constraints)

	PetscScalar mfmax;          // maximum melt fraction affecting viscosity reduction

	PetscInt    lmaxit;         // maximum number of local rheology iterations
	PetscScalar lrtol;          // local rheology iterations relative tolerance
	PetscInt    Phasetrans;     // Flag to activate phase transition routines
	PetscInt    Passive_Tracer; // Flag to activate passive tracer routine
	PetscScalar Adiabatic_gr;   // Adiabatic gradient

	PetscInt    actDike;        // Flag to activate dike, additional term on RHS of divergence
	PetscInt    useTk;          // activation flag for using temperature-dependent conductivity
	PetscInt    dikeHeat;       // activation flag for using Behn & Ito heat source in dike
};

//---------------------------------------------------------------------------
//.............. FDSTAG Jacobian and residual evaluation context ............
//---------------------------------------------------------------------------

struct JacRes
{
	// external handles
	Scaling    *scal;   // scaling
	TSSol      *ts;     // time-stepping parameters
	FDSTAG     *fs;     // staggered-grid layout
	FreeSurf   *surf;   // free surface
	BCCtx      *bc;     // boundary condition context
	DBPropDike *dbdike; // dike database
	DBMat      *dbm;    // material database
	PData      *Pd;      // Phase diagram

	// parameters and controls
	Controls ctrl;

	// coupled solution & residual vectors
	Vec gsol, gres; // global

	// solution variables
	SolVarCell  *svCell;   // cell centers
	SolVarEdge  *svXYEdge; // XY edges
	SolVarEdge  *svXZEdge; // XZ edges
	SolVarEdge  *svYZEdge; // YZ edges
	PetscScalar *svBuff;   // storage for phRat

	//=========
	// pressure
	//=========

	Vec lp_lith; // lithostatic pressure
	Vec lp_pore; // pore pressure

	//=======================
	// temperature parameters
	//=======================

	Vec gT;   // temperature solution vector (global)
	DM  DA_T; // temperature cell-centered grid with star stencil
	Mat Att;  // temperature preconditioner matrix
	Vec dT;   // temperature increment (global)
	Vec ge;   // energy residual (global)
	KSP tksp; // temperature diffusion solver

	// reference energy residual norm for automatic tolerance setting
	PetscScalar ts_ksp_ref_norm;

	//==========================
	// 2D integration primitives
	//==========================
	DM DA_CELL_2D; // 2D cell center grid

};
//---------------------------------------------------------------------------

// create residual & Jacobian evaluation context
PetscErrorCode JacResCreate(JacRes *jr, FB *fb);

PetscErrorCode JacResCreateData(JacRes *jr);

PetscErrorCode JacResReadRestart(JacRes *jr, FILE *fp);

PetscErrorCode JacResWriteRestart(JacRes *jr, FILE *fp);

// destroy residual & Jacobian evaluation context
PetscErrorCode JacResDestroy(JacRes *jr);

// compute effective inverse elastic parameter
PetscErrorCode JacResGetI2Gdt(JacRes *jr);

// form residual vector
PetscErrorCode JacResFormResidual(JacRes *jr, Vec x, Vec f);

// access solution vectors, set two-point constraints, prepare for interpolation (optionally)
PetscErrorCode JacResGetSolution(JacRes *jr, Vec x, Vec *lvx, Vec *lvy, Vec *lvz, Vec *lp, Vec *lT, AccessMode mode);

// restore solution vectors
PetscErrorCode JacResRestoreSolution(JacRes *jr, Vec *lvx, Vec *lvy, Vec *lvz, Vec *lp, Vec *lT);

// access current velocity
PetscErrorCode JacResGetVel(JacRes *jr, Vec x, Vec lvx, Vec lvy, Vec lvz);

// access current pressure
PetscErrorCode JacResGetPres(JacRes *jr, Vec x, Vec lp);

// get average pressure near the top surface
PetscErrorCode JacResGetPressShift(JacRes *jr, Vec lp);

// evaluate effective strain rate components in basic nodes
PetscErrorCode JacResGetEffStrainRate(JacRes *jr,
                                      Vec lvx,  Vec lvy,  Vec lvz,
                                      Vec ldxx, Vec ldyy, Vec ldzz,
                                      Vec ldxy, Vec ldxz, Vec ldyz);

// compute velocity gradients for output
PetscErrorCode JacResGetVelGrad(JacRes *jr,
                                Vec lvx,   Vec lvy,   Vec lvz,
                                Vec dvxdx, Vec dvxdy, Vec dvxdz,
                                Vec dvydx, Vec dvydy, Vec dvydz,
                                Vec dvzdx, Vec dvzdy, Vec dvzdz);

// compute components of vorticity vector
PetscErrorCode JacResGetVorticity(JacRes *jr,
                                  Vec lvx,  Vec lvy,  Vec lvz,
                                  Vec ldxy, Vec ldxz, Vec ldyz);

// initialize pressure
PetscErrorCode JacResInitPres(JacRes *jr);

// initialize pressure to lithostatic pressure
PetscErrorCode JacResInitLithPres(JacRes *jr, AdvCtx *actx, TSSol *ts);

// assemble residual
PetscErrorCode JacResAssembleRes(JacRes *jr, Vec f, Vec lfx, Vec lfy, Vec lfz, Vec gc);

// print residual
PetscErrorCode JacResViewRes(JacRes *jr);

//---------------------------------------------------------------------------

// compute maximum horizontal compressive stress (SHmax) orientation
PetscErrorCode JacResGetSHmax(JacRes *jr, Vec cx, Vec cy);

// compute maximum horizontal extension rate (EHmax) orientation
PetscErrorCode JacResGetEHmax(JacRes *jr, Vec cx, Vec cy);

//---------------------------------------------------------------------------
// Effective permeability functions
//---------------------------------------------------------------------------

PetscErrorCode JacResGetPermea(JacRes *jr, PetscInt bgPhase, PetscInt step, char *outfile);

//---------------------------------------------------------------------------
//......................   TEMPERATURE FUNCTIONS   ..........................
//---------------------------------------------------------------------------

PetscErrorCode JacResGetTempParam(
    JacRes      *jr,
    PetscScalar *phRat,
    PetscScalar *k_,      // conductivity
    PetscScalar *rho_Cp_, // volumetric heat capacity
    PetscScalar *rho_A_,  // volumetric radiogenic heat
    PetscScalar Tc,       // temperature of cell
    PetscScalar y_c,
    PetscInt J);          // coordinate of cell

// check whether thermal material parameters are properly defined
PetscErrorCode JacResCheckTempParam(JacRes *jr);

// setup temperature parameters
PetscErrorCode JacResCreateTempParam(JacRes *jr);

// destroy temperature parameters
PetscErrorCode JacResDestroyTempParam(JacRes *jr);

// initialize temperature from markers
PetscErrorCode JacResInitTemp(JacRes *jr);

// access current temperature
PetscErrorCode JacResGetTemp(JacRes *jr, Vec lT);

// compute temperature residual vector
PetscErrorCode JacResGetTempRes(JacRes *jr, PetscScalar dt);

// assemble temperature preconditioner matrix
PetscErrorCode JacResGetTempMat(JacRes *jr, PetscScalar dt);

//---------------------------------------------------------------------------
//......................   INTEGRATION FUNCTIONS   ..........................
//---------------------------------------------------------------------------

// compute overpressure field in the cell centers
PetscErrorCode JacResGetOverPressure(JacRes *jr, Vec p);

// compute lithostatic pressure in the cell centers
PetscErrorCode JacResGetLithoStaticPressure(JacRes *jr);

// compute pore pressure from phase properties and lithostatic stress
PetscErrorCode JacResGetPorePressure(JacRes *jr);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define SET_TPC(bc, a, k, j, i, pmdof) { \
    if(bc[k][j][i] == DBL_MAX) a[k][j][i] = pmdof; \
    else                       a[k][j][i] = 2.0*bc[k][j][i] - pmdof; }

//---------------------------------------------------------------------------
#endif


