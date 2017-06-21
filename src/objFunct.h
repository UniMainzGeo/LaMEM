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
 **    filename:   objFunct.h
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
//.....................   OBJECTIVE FUNCTION ROUTINES   .....................
//---------------------------------------------------------------------------
#ifndef __objFunct_h__
#define __objFunct_h__
//---------------------------------------------------------------------------

// maximum number of obervational types
#define _max_num_obs_ 7
#define _max_len_name_ 8

//-----------------------------------------------------------------------------

struct FB;
struct FreeSurf;

//-----------------------------------------------------------------------------
// model parameter type enumeration

#define _max_num_ModParam_type_ 30
#define _max_num_MatParam_type_ 30

enum PTypes// List of model parameter types (30)
{
	// -- material model parameter types --
	_RHO0_, _RHON_, _RHOC_,                             // density
	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
	_ETA0_,	_E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
	_BP_, _TAUP_, _GAMMA_, _Q_, _EP_, _VP_,             // Peierls creep
	_SHEAR_, _BULK_, _KP_,                              // elasticity
	_COHESION_, _FRICTION_, _CHSOFTID_, _FRSOFTID_,     // plasticity (Drucker-Prager)
	_ALPHA_, _CP_, _K_, _A_,                             // energy
	_MFR_                                               // melt fraction
	// -- others --
	// ... (geometry, pushing box , etc. ...)

};

/*
const char *PTypesName[] ={
		// -- material model parameter types --
		"rho0","rho_n","rho_c",                        // density
		"eta","Bd","Ed","Vd",                          // newtonian linear diffiusion
		"eta0","e0","Bn","n","En","Vn",                // power-law (dislocation) creep
		"Bp","taup","gamma","q","Ep","Vp",             // Peierls creep
		"shear","bulk","Kp",                           // elasticity
		"cohesion","friction","chSoftID","frSoftID",   // plasticity (Drucker Prager)
		"alpha","cp","k","A"                           // energy

		// -- others --
		// ... (geometry, pushing box , etc. ...)
};
*/

//-----------------------------------------------------------------------------
// Structure that holds inversion parameters
struct ModParam
{
	PetscInt         use;  // 0 = NO 1 = Tobi's inversion 2 = Compute gradients 3 = full inversion 4 = save this forward simulation as comparison simulation
	PetscInt         mdN;  // number of model parameters
	PetscInt         mID;  // current model number
	PetscInt        *phs;  // model phase number
	PetscInt        *typ;  // model parameter type
	PetscScalar     *val;  // model value
	PetscScalar     *grd;  // gradient value
	PetscScalar      mfit; // misfit value for current model parameters

	// Variables additionally needed for the adjoint TAO solver
	Vec              xini;      	// Comparison velocity field for adjoint inversion
	Vec              P;				// vector containing parameters
	Vec              fcconv;        // Vector containing all f/fini values to track convergence
	PetscInt         Ab;    		// Use adjoint bounds (only works with Tao)?
	PetscInt         Tao;    		// Use Tao?
	PetscInt         Adv;      		// Advect the point?
	PetscInt         count;			// iteration counter
	PetscInt         mdI;    		// number of indices
	PetscInt         Ap;        	// 1 = several indices ; 2 = whole domain ; 3 = surface
	PetscInt         reg;       	// 1 = Tikhonov regularization of the adjoint cost function ; 2 = total variation regularization (TV)
	PetscScalar      mfitini; 		// initial misfit value for current model parameters
	PetscScalar      tol; 		    // tolerance for F/Fini after which code has converged
	PetscScalar      factor1;   	// factor to multiply the gradients (should be set such that the highest gradient scales around 1/100 of its parameter ; only used without tao)
	PetscScalar      factor2;   	// factor that increases the convergence velocity (this value is added to itself after every succesful gradient descent ; only used without tao)
	PetscScalar      maxfactor2;	// limit on the factor (only used without tao)
	PetscScalar     *Ax;			// X-coordinates of comparison points
	PetscScalar     *Ay;			// Y-coordinates of comparison points
	PetscScalar     *Az;  			// Z-coordinates of comparison points
	PetscScalar     *Av;			// Velocity components of comparison points
	PetscScalar     *W;        		// Array of weights for the regularization
};

// observation type
enum ObsType
{
	_VELX_,            // 1: horizontal velocity (x) at the surface
	_VELY_,            // 2: horizontal velocity (y) at the surface
	_VELZ_,            // 3: vertical velocity (z) at the surface
	_TOPO_,            // 4: surface topography
	_BOUG_,            // 5: bouguer anomaly
	_ISA_,             // 6: orientation of isa (<-> sks-seismic anisotropy)
	_SHMAX_            // 7: orientation of SHmax

};

//---------------------------------------------------------------------------
//........................ Objective function object ........................
//---------------------------------------------------------------------------
struct ObjFunct
{
	FreeSurf     *surf;                 // free surface object
	char         *infile;               // input file name
	PetscBool    CompMfit;              // Compute misfit?
	PetscInt     otUse[_max_num_obs_+1];// array of boolean USED flags
	PetscInt     otN;                   // number of USED observation types
	PetscInt     ocN;                   // total number of observational constraints
	PetscScalar  err[_max_num_obs_];    // array containing individual sums of errors
	PetscScalar  errtot;                // total error
	Vec          obs[_max_num_obs_];    // vectors containing the observations
	Vec          qul[_max_num_obs_];    // vectors containing quality info (quality/sigma)^2, where quality (0..1)

	// missing ...
	// (data) covariance matrix

};
//---------------------------------------------------------------------------

// create objective function object
PetscErrorCode ObjFunctCreate(ObjFunct *objf, ModParam *IOparam, FreeSurf *surf, FB *fb);

// destroy object
PetscErrorCode ObjFunctDestroy(ObjFunct *objf);

// read command line options
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[], FB *fb);

// compute error
PetscErrorCode ObjFunctCompErr(ObjFunct *objf);

// compute weighted Least square error for surface vectors
PetscErrorCode VecErrSurf(Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal);

//---------------------------------------------------------------------------

#endif
