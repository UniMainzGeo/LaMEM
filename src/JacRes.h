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
struct Tensor2RN;
struct PData;

//---------------------------------------------------------------------------
//.....................   Deviatoric solution variables   ...................
//---------------------------------------------------------------------------

struct SolVarDev
{
	PetscScalar  DII;   // effective strain rate
	PetscScalar  eta;   // effective tangent viscosity
	PetscScalar  I2Gdt; // inverse elastic parameter (1/2G/dt)
	PetscScalar  Hr;    // shear heating term (partial)
	PetscScalar  DIIpl; // plastic strain rate
	PetscScalar  APS;   // accumulated plastic strain
	PetscScalar  PSR;   // plastic strain-rate contribution
	PetscScalar  dEta;  // dEta/dDII derivative (Jacobian)
	PetscScalar  fr;    // effective friction coefficient (Jacobian)
	PetscScalar  yield; // average yield stress in control volume
	PetscScalar  mf;    // melt fraction
};

//---------------------------------------------------------------------------
//.....................   Volumetric solution variables   ...................
//---------------------------------------------------------------------------

struct SolVarBulk
{
	PetscScalar  theta; // volumetric strain rate
	PetscScalar  rho;   // strain- & temperature-dependent density
	PetscScalar  IKdt;  // inverse bulk elastic parameter (1/K/dt)
	PetscScalar  alpha; // effective thermal expansion
	PetscScalar  Tn;    // history temperature
	PetscScalar  pn;    // history pressure
	PetscScalar  rho_pd;// Density from phase diagram
	PetscScalar  rho_pf;// Fluid Density from phase diagram
	PetscScalar  mf;    // Melt fraction from phase diagram

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
	PetscScalar  eta_creep;     // effective creep viscosity (output)
	PetscScalar  eta_vp;        // viscoplastic viscosity (output)
	PetscScalar  U[3];          // displacement
	PetscScalar  es;            // nadai strain (octahedral shear strain)
	PetscScalar  nu;            // lodes ratio
	PetscScalar  uxx, uyy, uzz; // stretch tensor
    PetscScalar  FSA[3];        // major strain axis
    PetscScalar  tr;			// FSA trend
    PetscScalar  dp;			// FSA trend

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
    PetscScalar  u;     // xy, xz, yz stretch tensor components

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
	PetscScalar shearHeatEff;  // shear heating efficiency parameter [0 - 1]
	PetscScalar biot;          // Biot pressure parameter [0 - 1]

	PetscInt    actTemp;	   // temperature diffusion activation flag
	PetscInt    actExp;	       // thermal expansion activation flag
	PetscInt    actSteadyTemp; // steady-state temperature initial guess flag
	PetscInt    pShiftAct;     // pressure shift activation flag (zero pressure in the top cell layer)
	PetscScalar pShift;        // pressure shift for plasticity model and output
	PetscInt    initGuess;     // initial guess activation flag
	PetscInt    pLithoVisc;    // use lithostatic pressure for creep laws
	PetscInt    pLithoPlast;   // use lithostatic pressure for plasticity
	PetscInt    pLimPlast;     // limit pressure at first iteration for plasticity
	PetscInt    jac_mat_free;  // matrix-free analytical Jacobian activation flag
	PetscInt    quasiHarmAvg;  // use quasi-harmonic averaging regularization for plasticity

	PetscScalar eta_min;       // minimum viscosity
	PetscScalar inv_eta_max;   // inverse of maximum viscosity
	PetscScalar eta_ref;       // reference viscosity (initial guess)
	PetscScalar TRef;          // reference temperature
	PetscScalar Rugc;          // universal gas constant
	PetscScalar DII_ref;       // background (reference) strain-rate
	PetscScalar minCh;         // minimum cohesion
	PetscScalar minFr;         // minimum friction
	PetscScalar tauUlt;        // ultimate yield stress

	PetscScalar cf_eta_min;    // visco-plastic regularization parameter (plasticity)
	PetscScalar n_pw;          // power-law regularization parameter (plasticity)

	PetscScalar rho_fluid;     // fluid density
	GWLevelType gwType;        // type of ground water level (none, top, surf, level)
	PetscScalar gwLevel;       // fixed ground water level

	PetscInt    getPermea;     // effective permeability computation activation flag
	PetscInt    rescal;        // stensil rescaling flag (for interval constraints)

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
	DBMat    *dbm;   // material database

	// parameters and controls
	Controls ctrl;

	// coupled solution & residual vectors
	Vec gsol, gres; // global

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

	// Phase diagram
	PData       *Pd;

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

// copy residuals from local to global vectors, enforce boundary constraints
PetscErrorCode JacResCopyRes(JacRes *jr, Vec f);

// copy momentum residuals from global to local vectors for output
PetscErrorCode JacResCopyMomentumRes(JacRes *jr, Vec f);

// copy continuity residuals from global to local vectors for output
PetscErrorCode JacResCopyContinuityRes(JacRes *jr, Vec f);

PetscErrorCode JacResViewRes(JacRes *jr);

//---------------------------------------------------------------------------
// Infinite Strain Axis (ISA) computation functions
//---------------------------------------------------------------------------

// compute velocity gradient and normalized velocities at cell center
PetscErrorCode getGradientVel(
	FDSTAG *fs, PetscScalar ***lvx, PetscScalar ***lvy, PetscScalar ***lvz,
	PetscInt i, PetscInt j, PetscInt k, PetscInt sx, PetscInt sy, PetscInt sz,
	Tensor2RN *L, PetscScalar *vel, PetscScalar *pvnrm);

// compute Infinite Strain Axis (ISA)
PetscErrorCode JacResGetISA(JacRes *jr);

// compute Grain Orientation Lag (GOL) parameter
PetscErrorCode JacResGetGOL(JacRes *jr);

//---------------------------------------------------------------------------

// compute maximum horizontal compressive stress (SHmax) orientation
PetscErrorCode JacResGetSHmax(JacRes *jr);

// compute maximum horizontal extension rate (EHmax) orientation
PetscErrorCode JacResGetEHmax(JacRes *jr);

//---------------------------------------------------------------------------
// Effective permeability functions
//---------------------------------------------------------------------------

PetscErrorCode JacResGetPermea(JacRes *jr, PetscInt bgPhase, PetscInt step);

//---------------------------------------------------------------------------
//......................   TEMPERATURE FUNCTIONS   ..........................
//---------------------------------------------------------------------------

PetscErrorCode JacResGetTempParam(
	JacRes      *jr,
	PetscScalar *phRat,
	PetscScalar *k_,      // conductivity
	PetscScalar *rho_Cp_, // volumetric heat capacity
	PetscScalar *rho_A_); // volumetric radiogenic heat

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


