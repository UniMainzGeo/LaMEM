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
 **    filename:   constEq.h
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
	PetscScalar  fr;    // effective friction coefficient

} ConstEqCtx;

//---------------------------------------------------------------------------

// setup nonlinear constitutive equation evaluation context
// evaluate dependence on constant parameters (pressure, temperature)
PetscErrorCode ConstEqCtxSetup(
	ConstEqCtx  *ctx,      // evaluation context
	Material_t  *mat,      // phase parameters
	Controls    *ctrl,     // parameters and controls
	PetscScalar  DII,      // effective strain-rate
	PetscScalar  APS,      // accumulated plastic strain
	PetscScalar  dt,       // time step
	PetscScalar  p,        // pressure
	PetscScalar  p_lithos, // lithostatic pressure
	PetscScalar  p_pore,   // pore pressure
	PetscScalar  T);       // temperature

// compute residual of the visco-elastic constitutive equation
PetscScalar GetConsEqRes(PetscScalar eta, void *pctx);

PetscErrorCode GetEffVisc(
	ConstEqCtx  *ctx,
	Controls    *ctrl,        // parameters and controls
	PetscScalar *eta_total,
	PetscScalar *eta_creep,
	PetscScalar *eta_vp,
	PetscScalar *DIIpl,
	PetscScalar *dEta,
	PetscScalar *fr);

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
	PetscScalar *eta_creep, // creep viscosity (for output)
	PetscScalar *eta_vp,    // viscoplastic viscosity (for output)
	PetscInt     numPhases, // number phases
	Material_t  *phases,    // phase parameters
	PetscScalar *phRat,     // phase ratios
	Controls    *ctrl,       // parameters and controls
	PetscScalar  p_lithos,  // lithostatic pressure
	PetscScalar  p_pore,    // pore pressure
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T);        // temperature

// Evaluate volumetric constitutive equations in control volume
PetscErrorCode VolConstEq(
	SolVarBulk  *svBulk,    // solution variables
	PetscInt     numPhases, // number phases
	Material_t  *phases,    // phase parameters
	PetscScalar *phRat,     // phase ratios
	Controls    *ctrl,      // parameters and controls
	PetscScalar  depth,     // depth for depth-dependent density model
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T);        // temperature

// compute stress, plastic strain-rate and shear heating term on cell
PetscErrorCode GetStressCell(
		SolVarCell  *svCell, // solution variables
		Controls    *ctrl,   // parameters and controls
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz);   // ...

// compute stress, plastic strain-rate and shear heating term on edge
PetscErrorCode GetStressEdge(
	SolVarEdge  *svEdge, // solution variables
	Controls    *ctrl,    // parameters and controls
	PetscScalar  d);     // effective shear strain rate component

//---------------------------------------------------------------------------
#endif
