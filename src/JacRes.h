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
//---------------------------------------------------------------------------
// max number of phases
#define max_num_phases 32

// max number of soft laws
#define max_num_soft   10

//---------------------------------------------------------------------------

// FDSTAG Jacobian and residual evaluation context
typedef struct
{
	// external handles
	Scaling  *scal; // scaling
	TSSol    *ts;   // time-stepping parameters
	FDSTAG   *fs;   // staggered-grid layout
	BCCtx    *bc;   // boundary condition context

	// parameters & controls
	PetscScalar grav[SPDIM]; // global gravity components
	PetscScalar FSSA;        // density gradient penalty parameter
	PetscScalar gtol;	     // geometry tolerance
	PetscInt    actTemp;	 // temperature diffusion activation flag
	PetscInt    pShiftAct;   // pressure shift activation flag (zero pressure in the top cell layer)

	// phase parameters
	PetscInt     numPhases;              // number phases
	Material_t   phases[max_num_phases]; // phase parameters
	PetscInt     numSoft;                // number material softening laws
	Soft_t       matSoft[max_num_soft];  // material softening law parameters
	MatParLim    matLim;                 // phase parameters limiters

	// external and runtime parameters
	PetscInt    AirPhase;    // air phase number
	PetscScalar avg_topo;    // average topography (a copy from free surface)
	PetscScalar pShift;      // pressure shift for plasticity model and output

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
	Vec gp;        // global
	Vec lp;        // local (ghosted)
	Vec lp_lithos; // lithostatic pressure

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

} JacRes;
//---------------------------------------------------------------------------

// create residual & Jacobian evaluation context
PetscErrorCode JacResCreate(JacRes *jr, FB *fb);

PetscErrorCode JacResCreateData(JacRes *jr);

PetscErrorCode JacResReadRestart(JacRes *jr, FILE *fp);

PetscErrorCode JacResWriteRestart(JacRes *jr, FILE *fp);

// destroy residual & Jacobian evaluation context
PetscErrorCode JacResDestroy(JacRes *jr);

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

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopyVel(JacRes *jr, Vec x);

// copy solution from global to local vectors, enforce boundary constraints
PetscErrorCode JacResCopyPres(JacRes *jr, Vec x);

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

PetscErrorCode JacResSetVelRotation(JacRes *jr);

//---------------------------------------------------------------------------

// get maximum inverse time step on local domain
PetscErrorCode getMaxInvStep1DLocal(Discret1D *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax);

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
PetscErrorCode JacResGetTempRes(JacRes *jr);

// assemble temperature preconditioner matrix
PetscErrorCode JacResGetTempMat(JacRes *jr);

//---------------------------------------------------------------------------
//......................   INTEGRATION FUNCTIONS   ..........................
//---------------------------------------------------------------------------

// compute overpressure field in the cell centers
PetscErrorCode JacResGetOverPressure(JacRes *jr, Vec p);

// compute lithostatic pressure in the cell centers
PetscErrorCode JacResGetLithoStaticPressure(JacRes *jr);

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


