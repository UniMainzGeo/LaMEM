//---------------------------------------------------------------------------
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#ifndef __constEq_h__
#define __constEq_h__
//---------------------------------------------------------------------------

// nonlinear constitutive equations evaluation context
typedef struct
{
	PetscScalar  DII;   // effective strain-rate
	PetscScalar  A_els; // elasticity constant
	PetscScalar  A_dif; // diffusion constant
	PetscScalar  A_dis; // dislocation constant
	PetscScalar  N_dis; // dislocation exponent
	PetscScalar  A_prl; // Peierls constant
	PetscScalar  N_prl; // Peierls exponent
	PetscScalar  taupl; // plastic yield stress
	PetscBool    cfsol; // closed-form solution flag

} ConstEqCtx;

//---------------------------------------------------------------------------

// setup nonlinear constitutive equation evaluation context
// evaluate dependence on constant parameters (pressure, temperature)
PetscErrorCode ConstEqCtxSetup(
	ConstEqCtx  *ctx,  // evaluation context
	Material_t  *mat,  // phase parameters
	MatParLim   *lim,  // phase parameters limits
	PetscScalar  DII,  // effective strain-rate
	PetscScalar  APS,  // accumulated plastic strain
	PetscScalar  dt,   // time step
	PetscScalar  p,    // pressure
	PetscScalar  T);    // temperature

// compute residual of the visco-elastic constitutive equation
PetscScalar GetConsEqRes(PetscScalar eta, void *pctx);

// solve effective viscosity from nonlinear visco-elastic constitutive equations
PetscErrorCode GetEffVisc(
	ConstEqCtx  *ctx,
	MatParLim   *lim,
	PetscScalar *eta,
	PetscScalar *DIIpl);

// apply strain softening to a parameter (friction, cohesion)
PetscScalar ApplyStrainSoft(Soft_t *sl, PetscScalar APS, PetscScalar par);

// compute inverse elastic viscosity in control volume
PetscScalar GetI2Gdt(
	PetscInt     numPhases,
	Material_t  *phases,
	PetscScalar *phRat,
	PetscScalar  dt);

// Evaluate deviatoric constitutive equations in control volume
PetscErrorCode DevConstEq(
	SolVarDev   *svDev,     // solution variables
	PetscInt     numPhases, // number phases
	Material_t  *phases,    // phase parameters
	PetscScalar *phRat,     // phase ratios
	MatParLim   *lim,       // phase parameters limits
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T);        // temperature

// Evaluate volumetric constitutive equations in control volume
PetscErrorCode VolConstEq(
	SolVarBulk  *svBulk,    // solution variables
	PetscInt     numPhases, // number phases
	Material_t  *phases,    // phase parameters
	PetscScalar *phRat,     // phase ratios
	MatParLim   *lim,       // phase parameters limits
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T);        // temperature

// compute stress, plastic strain-rate and shear heating term on cell
PetscErrorCode GetStressCell(
		SolVarCell  *svCell, // solution variables
		MatParLim   *lim,    // phase parameters limits
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz);   // ...

// compute stress, plastic strain-rate and shear heating term on edge
PetscErrorCode GetStressEdge(
	SolVarEdge  *svEdge, // solution variables
	MatParLim   *lim,    // phase parameters limits
	PetscScalar  d);     // effective shear strain rate component

//---------------------------------------------------------------------------
#endif
