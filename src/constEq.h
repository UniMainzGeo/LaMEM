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

struct Material_t;
struct Soft_t;
struct Controls;
struct SolVarDev;
struct SolVarBulk;
struct SolVarCell;
struct SolVarEdge;
struct PData;

//---------------------------------------------------------------------------

// constitutive equations evaluation context
struct ConstEqCtx
{
	// database parameters
	PetscInt     numPhases; // number phases
	Material_t  *phases;    // phase parameters
	Soft_t      *soft;      // material softening laws
	Controls    *ctrl;      // parameters and controls
	PData       *pd;        // phase diagram data
	PetscScalar  dt;        // time step
	PetscScalar  stats[3];  // total number of [starts, fails, iterations]

	// control volume parameters
	PetscScalar *phRat;  // phase ratios in the control volume
	PetscScalar  p;      // pressure
	PetscScalar  p_lith; // lithostatic pressure
	PetscScalar  p_pore; // pore pressure
	PetscScalar  depth;  // depth for depth-dependent density model
	PetscScalar  T;      // temperature
	PetscScalar  DII;    // effective strain rate
	PetscScalar  APS;    // accumulated plastic strain

	// phase parameters
	PetscScalar  A_els;  // elasticity constant
	PetscScalar  A_dif;  // diffusion constant
	PetscScalar  A_max;  // upper bound constant
	PetscScalar  A_dis;  // dislocation constant
	PetscScalar  N_dis;  // dislocation exponent
	PetscScalar  A_prl;  // Peierls constant
	PetscScalar  N_prl;  // Peierls exponent
	PetscScalar  taupl;  // plastic yield stress

	// control volume results
	PetscScalar  eta;    // effective viscosity
	PetscScalar  eta_cr; // creep viscosity
	PetscScalar  eta_vp; // visco-plastic viscosity
	PetscScalar  DIIdif; // diffusion creep strain rate
	PetscScalar  DIIdis; // dislocation creep strain rate
	PetscScalar  DIIprl; // Peierls creep strain rate
	PetscScalar  DIIpl;  // plastic strain rate
	PetscScalar  yield;  // yield stress

};

//---------------------------------------------------------------------------
// evaluate volumetric constitutive equations in control volume
PetscErrorCode volConstEq(SolVarBulk *svBulk, ConstEqCtx *ctx);

// evaluate deviatoric constitutive equations in control volume
PetscErrorCode devConstEq(SolVarDev *svDev, ConstEqCtx *ctx);

// compute phase viscosities and strain rate partitioning
PetscErrorCode getPhaseVisc(ConstEqCtx *ctx, PetscInt ID);

// setup nonlinear constitutive equation evaluation context
PetscErrorCode setUpPhase(ConstEqCtx *ctx, PetscInt ID);

// compute stress, plastic strain-rate and shear heating term on cell
PetscErrorCode getStressCell(
		SolVarCell  *svCell, // solution variables
		ConstEqCtx  *ctx,    // evaluation context
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz);   // ...

// compute stress, plastic strain-rate and shear heating term on edge
PetscErrorCode getStressEdge(
		SolVarEdge  *svEdge, // solution variables
		ConstEqCtx  *ctx,    // evaluation context
		PetscScalar  d);     // effective shear strain rate component

// compute residual of the visco-elastic constitutive equation
PetscScalar getConsEqRes(PetscScalar eta, void *pctx);

// apply strain softening to a parameter (friction, cohesion)
PetscScalar applyStrainSoft(
		Soft_t      *soft, // material softening laws
		PetscInt     ID,   // softening law ID
		PetscScalar  APS,  // accumulated plastic strain
		PetscScalar  par); // softening parameter

// compute inverse elastic parameter in control volume
PetscScalar getI2Gdt(
		PetscInt     numPhases, // number phases
		Material_t  *phases,    // phase parameters
		PetscScalar *phRat,     // phase ratios in the control volume
		PetscScalar  dt);       // time step

//---------------------------------------------------------------------------
//.............................. PHASE DIAGRAM  .............................
//---------------------------------------------------------------------------

PetscErrorCode setDataPhaseDiagram(
		PData       *pd,
		PetscScalar  p,
		PetscScalar  T,
		char         pdn[]);

//---------------------------------------------------------------------------
#endif
