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
 **    filename:   JacRes.h
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
//...................   FDSTAG JACOBIAN AND RESIDUAL  .......................
//---------------------------------------------------------------------------
#ifndef __JacRes_h__
#define __JacRes_h__

struct FB;
struct Scaling;
struct TSSol;
struct FDSTAG;
struct FreeSurf;
struct BCCtx;
struct DBMat;
struct DBPropDike;
struct Tensor2RN;
struct PData;
struct AdvCtx;
//struct ConstEqCtx;

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

struct Controls
{
	PetscScalar grav[3];       // global gravity components
	PetscScalar FSSA;          // free surface stabilization parameter [0 - 1]
	PetscInt    FSSA_allVel;   // Use all velocity components for FSSA?
	PetscScalar shearHeatEff;  // shear heating efficiency parameter [0 - 1]
	PetscScalar biot;          // Biot pressure parameter [0 - 1]

	PetscInt    AdiabHeat;		// Adiabatic Heating flag
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
	PetscInt    printNorms;		// priny norms of velocity/pressure/temperature?

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

  PetscInt    useTk;     // activation flag for using temperature-dependent conductivity

  PetscInt  dikeHeat;   // activation flag for using Behn & Ito heat source in dike
};

//---------------------------------------------------------------------------
//.............. FDSTAG Jacobian and residual evaluation context ............
//---------------------------------------------------------------------------

struct JacRes
{
	// external handles
	Scaling  *scal;  // scaling
	TSSol    *ts;    // time-stepping parameters
	FDSTAG   *fs;    // staggered-grid layout
	FreeSurf *surf;  // free surface
	BCCtx    *bc;    // boundary condition context
  DBPropDike *dbdike; // dike database
	DBMat    *dbm;   // material database
  //  ConstEqCtx *ctx;
  
	// parameters and controls
	Controls ctrl;

	// coupled solution & residual vectors
	Vec gsol, gres; // global

	// velocity	components
	Vec gvx,  gvy, gvz;  // global
	Vec lvx,  lvy, lvz;  // local (ghosted)
	Vec dvxdx,dvxdy, dvxdz,dvydx,dvydy,dvydz,dvzdx,dvzdy,dvzdz;  // velocity tensor components

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

	// For almost all the purposes only one center-based array is necessary instead of three
	// for example - strain rate contributions from centers can be stored in one array

	// IN GENERAL GET RID OF BUFFER VECTORS, USE LOCAL DMGetLocalVector (GET_INIT_LOCAL_VECTOR)

	// pressure
	Vec gp;      // global
	Vec lp;      // local (ghosted)
	Vec lp_lith; // lithostatic pressure
	Vec lp_pore; // pore pressure

	// continuity residual
	Vec gc; // global

	// corner buffer
	Vec lbcor; // local (ghosted)

	// solution variables
	SolVarCell  *svCell;   // cell centers
	SolVarEdge  *svXYEdge; // XY edges
	SolVarEdge  *svXZEdge; // XZ edges
	SolVarEdge  *svYZEdge; // YZ edges
	PetscScalar *svBuff;   // storage for phRat
	PetscScalar  mean_p;  // average lithostatic pressure

	// Phase diagram
	PData       *Pd;

	// Adjoint field based gradients
	Vec          lgradfield;
	Vec          phi; // PSD context

	//=======================
	// temperature parameters
	//=======================

	Vec lT;   // temperature (box stencil, active even without diffusion)
	DM  DA_T; // temperature cell-centered grid with star stencil
	Mat Att;  // temperature preconditioner matrix
	Vec dT;   // temperature increment (global)
	Vec ge;   // energy residual (global)
	KSP tksp; // temperature diffusion solver

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

// form residual vector
PetscErrorCode JacResFormResidual(JacRes *jr, Vec x, Vec f);

// compute effective inverse elastic parameter
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

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopyVel(JacRes *jr, Vec x);

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopyPres(JacRes *jr, Vec x);

// initialize pressure
PetscErrorCode JacResInitPres(JacRes *jr);

// initialize pressure to lithostatic pressure
PetscErrorCode JacResInitLithPres(JacRes *jr, AdvCtx *actx);

// copy residuals from local to global vectors, enforce boundary constraints
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f);

// copy momentum residuals from global to local vectors for output
PetscErrorCode JacResCopyMomentumRes(JacRes *jr, Vec f);

// copy continuity residuals from global to local vectors for output
PetscErrorCode JacResCopyContinuityRes(JacRes *jr, Vec f);

PetscErrorCode JacResViewRes(JacRes *jr);

//---------------------------------------------------------------------------

// compute velocity gradient and normalized velocities at cell center
PetscErrorCode getGradientVel(
	FDSTAG *fs, PetscScalar ***lvx, PetscScalar ***lvy, PetscScalar ***lvz,
	PetscInt i, PetscInt j, PetscInt k, PetscInt sx, PetscInt sy, PetscInt sz,
	Tensor2RN *L, PetscScalar *vel, PetscScalar *pvnrm);

//---------------------------------------------------------------------------

// compute maximum horizontal compressive stress (SHmax) orientation
PetscErrorCode JacResGetSHmax(JacRes *jr);

// compute maximum horizontal extension rate (EHmax) orientation
PetscErrorCode JacResGetEHmax(JacRes *jr);

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
	PetscScalar Tc);      // temperature of cell 

// check whether thermal material parameters are properly defined
PetscErrorCode JacResCheckTempParam(JacRes *jr);

// setup temperature parameters
PetscErrorCode JacResCreateTempParam(JacRes *jr);

// destroy temperature parameters
PetscErrorCode JacResDestroyTempParam(JacRes *jr);

// initialize temperature from markers
PetscErrorCode JacResInitTemp(JacRes *jr);

// correct temperature for diffusion (Newton update)
PetscErrorCode JacResUpdateTemp(JacRes *jr);

// apply temperature two-point constraints
PetscErrorCode JacResApplyTempBC(JacRes *jr);

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

#define SET_EDGE_CORNER(n, a, K, J, I, k, j, i, pmdof) \
	a[K][J][I] = a[k][j][I] + a[k][J][i] + a[K][j][i] - 2.0*pmdof;

//---------------------------------------------------------------------------
#endif


