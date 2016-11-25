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
 **    filename:   solVar.h
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
//.................   SOLUTION AND HISTORY VARIABLES   ......................
//---------------------------------------------------------------------------
#ifndef __solVar_h__
#define __solVar_h__

//---------------------------------------------------------------------------
//....    Non-Symmetric second rank tensor (gradient & rotation tensors) ....
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar xx, xy, xz;
	PetscScalar yx, yy, yz;
	PetscScalar zx, zy, zz;

} Tensor2RN;

//---------------------------------------------------------------------------
//.......   Symmetric second rank tensor (stress & strain tensors)   ........
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar xx, xy, xz;
	PetscScalar     yy, yz;
	PetscScalar         zz;

} Tensor2RS;

//---------------------------------------------------------------------------
//............   Material marker (history variables advection)   ............
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt    phase; // phase identifier
	PetscScalar X[3];  // global coordinates
	PetscScalar p;     // pressure
	PetscScalar T;     // temperature
	PetscScalar APS;   // accumulated plastic strain
	Tensor2RS   S;     // deviatoric stress
	PetscScalar U[3];  // displacement

} Marker;

//---------------------------------------------------------------------------
//.....................   Deviatoric solution variables   ...................
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar  DII;   // effective strain rate
	PetscScalar  eta;   // effective tangent viscosity
	PetscScalar  I2Gdt; // inverse elastic viscosity (1/2G/dt)
	PetscScalar  Hr;    // shear heating term (partial)
	PetscScalar  DIIpl; // plastic strain rate
	PetscScalar  APS;   // accumulated plastic strain
	PetscScalar  PSR;   // plastic strain-rate contribution
	PetscScalar  dEta;  // dEta/dDII derivative (Jacobian)
	PetscScalar  fr;    // effective friction coefficient (Jacobian)

} SolVarDev;

//---------------------------------------------------------------------------
//.....................   Volumetric solution variables   ...................
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar  theta; // volumetric strain rate
	PetscScalar  rho;   // strain- & temperature-dependent density
	PetscScalar  IKdt;  // inverse bulk elastic viscosity (1/K/dt)
	PetscScalar  alpha; // effective thermal expansion
	PetscScalar  Tn;    // history temperature
	PetscScalar  pn;    // history pressure

} SolVarBulk;

//---------------------------------------------------------------------------
//........................   Cell solution variables   ......................
//---------------------------------------------------------------------------

typedef struct
{
	SolVarDev    svDev;         		// deviatoric variables
	SolVarBulk   svBulk;        		// volumetric variables
	PetscScalar  sxx, syy, szz; 		// deviatoric stress
	PetscScalar  hxx, hyy, hzz; 		// history stress (elastic)
	PetscScalar  dxx, dyy, dzz; 		// total deviatoric strain rate
	PetscScalar *phRat;         		// phase ratios in the control volume
	PetscScalar  eta_creep;     		// effective creep viscosity (output)
	PetscScalar  eta_viscoplastic;     	// viscoplastic viscosity (output)
	PetscScalar  U[3];          		// displacement

} SolVarCell;

//---------------------------------------------------------------------------
//........................   Edge solution variables   ......................
//---------------------------------------------------------------------------

typedef struct
{
	SolVarDev    svDev; // deviatoric variables
	PetscScalar  s;     // xy, xz, yz deviatoric stress components
	PetscScalar  h;     // xy, xz, yz history stress components (elastic)
	PetscScalar  d;     // xy, xz, yz total deviatoric strain rate components
	PetscScalar  ws;    // normalization for distance-dependent interpolation
	PetscScalar *phRat; // phase ratios in the control volume

} SolVarEdge;

//---------------------------------------------------------------------------
//.......................   Softening Law Parameters  .......................
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt    ID;   // softening law ID
	PetscScalar APS1; // begin of softening APS
	PetscScalar APS2; // end of softening APS
	PetscScalar A;    // reduction ratio

} Soft_t;

//---------------------------------------------------------------------------
//......................   Material parameter table   .......................
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt     ID;      // material ID
	// density parameters
	PetscScalar  rho;     // reference density                          [kg/m^3]
	PetscScalar  rho_n;   // depth-dependent density model parameter    [ ]
	PetscScalar  rho_c;   // depth-dependent density model parameter    [1/m]
	PetscScalar  beta;    // pressure-dependent density model parameter [1/Pa]
	// elasticity parameters
	PetscScalar  K;       // bulk modulus                               [Pa]
	PetscScalar  Kp;      // pressure dependence parameter              [ ]
	PetscScalar  G;       // shear modulus                              [Pa]
	// diffusion creep parameters
	PetscScalar  Bd;      // pre-exponential constant                   [1/Pa/s]
	PetscScalar  Ed;      // activation energy                          [J/mol]
	PetscScalar  Vd;      // activation volume                          [m^3/mol]
	// dislocation creep parameters
	PetscScalar  Bn;      // pre-exponential constant                   [1/Pa^n/s]
	PetscScalar  n;       // power law exponent                         [ ]
	PetscScalar  En;      // activation energy                          [J/mol]
	PetscScalar  Vn;      // activation volume                          [m^3/mol]
	// Peierls creep parameters
	PetscScalar  Bp;      // pre-exponential constant                   [1/s]
	PetscScalar  Ep;      // activation energy                          [J/mol]
	PetscScalar  Vp;      // activation volume                          [m^3/mol]
	PetscScalar  taup;    // scaling stress                             [Pa]
	PetscScalar  gamma;   // approximation parameter                    [ ]
	PetscScalar  q;       // stress-dependence parameter                [ ]
	// plasticity parameters
	PetscScalar  fr;      // friction angle                             [deg]
	PetscScalar  ch;      // cohesion                                   [Pa]
	Soft_t      *frSoft;  // friction softening law parameters
	Soft_t      *chSoft;  // cohesion softening law parameters
	// thermal parameters
	PetscScalar  alpha;   // thermal expansivity                        [1/K]
	PetscScalar  Cp;      // cpecific heat (capacity)                   [J/kg/K]
	PetscScalar  k;       // thermal conductivity                       [W/m/k]
	PetscScalar  A;       // radiogenic heat production                 [W/kg]

} Material_t;

//---------------------------------------------------------------------------
//.....................   Material parameter limiters .......................
//---------------------------------------------------------------------------

typedef struct
{
	// viscosity limits
	PetscScalar eta_min;
	PetscScalar inv_eta_max;
	PetscScalar eta_ref;     // reference viscosity (initial guess)
	// reference temperature
	PetscScalar TRef;
	// universal gas constant
	PetscScalar Rugc;
	// background (reference) strain-rate
	PetscScalar DII_ref;
	// plasticity parameters limits
	PetscScalar minCh;  // minimum cohesion
	PetscScalar minFr;  // maximum friction
	PetscScalar tauUlt; // ultimate yield stress
	// fluid density for depth-dependent density model
	PetscScalar rho_fluid;
	// thermo-mechanical coupling controls
	PetscScalar shearHeatEff; // shear heating efficiency parameter [0 - 1]
	// rheology controls
	PetscInt    quasiHarmAvg; // quasi-harmonic averaging regularization flag (plasticity)
	PetscScalar cf_eta_min;   // visco-plastic regularization parameter (plasticity)
	PetscScalar n_pw;         // power-law regularization parameter (plasticity)
	// pressure parameters
	PetscInt    pLithoVisc;   // use lithostatic pressure for creep laws
	PetscInt    pLithoPlast;  // use lithostatic pressure for plasticity
	PetscInt    pLimPlast;    // limit pressure at first iteration for plasticity
	// direction to the North for stress orientation
	// counter-clockwise positive measured from x-axis
	PetscScalar theta_north;
	// print warning messages
	PetscInt    warn;
	// matrix-free closed-form jacobian
	PetscInt    jac_mat_free;
	// initial guess computation flag
	PetscInt    initGuess;

} MatParLim;
//---------------------------------------------------------------------------
#endif
