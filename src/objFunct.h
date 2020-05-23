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
#include "parsing.h"	// filebuffer
#include "adjoint.h"    // defines the global variables _MAX_PAR_ and _MAX_OBS_, which we need here

struct FB;
struct FreeSurf;

//-----------------------------------------------------------------------------
enum PTypes// List of model parameter types (30)
{
	// -- material model parameter types --
	_RHO0_, _RHON_, _RHOC_,                             // density
	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
	_ETA0_,	_E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
	_BP_, _TAUP_, _GAMMA_, _Q_, _EP_, _VP_,             // Peierls creep
	_SHEAR_, _BULK_, _KP_,                              // elasticity
	_COHESION_, _FRICTION_, _CHSOFTID_, _FRSOFTID_,     // plasticity (Drucker-Prager)
	_ALPHA_, _CP_, _K_, _A_                             // energy
	// -- others --
	// ... (geometry, pushing box , etc. ...)

};

enum InvTypes// List of inversion types
{
	_none_, _inversion_, _adjointgradients_, _gradientdescent_, _syntheticforwardrun_, 
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
	PetscInt         use;                               // Choose one of InvTypes
	PetscInt         mdN;                               // number of model parameters
	PetscInt         mID;                               // current model number
	char 			 type_name[_MAX_PAR_][_str_len_];   // stores the name of the adjoint parameters
	PetscBool 		 FD_gradient[_MAX_PAR_];			// Compute gradient via (brute force) finite differences, or with
    PetscInt      	 phs[_MAX_PAR_];                    // phase of the parameter
	PetscScalar      grd[_MAX_PAR_];                    // gradient value
	PetscScalar     *val;                               // model value
	PetscScalar      mfit;                              // misfit value for current model parameters
    DBMat            dbm_modified;                      // holds the (modified) LaMEM material database
	FB 				*fb;								// holds a copy of the filebuffer	

	// Variables additionally needed for the adjoint TAO solver
	Vec              xini;      	                    // Comparison velocity field for adjoint inversion
	Vec              P;				                    // vector containing parameters
	Vec              fcconv;                            // Vector containing all f/fini values to track convergence
	PetscInt         Ab;    		                    // Use adjoint bounds (only works with Tao)?
	PetscInt         Tao;    		                    // Use Tao?
	PetscInt         Adv;      		                    // Advect the point?
	PetscInt         count;			                    // iteration counter
	PetscInt         SCF;                               // Scale cost function?
	PetscInt         mdI;    		                    // number of indices
	PetscInt         Ap;        	                    // 1 = several indices ; 2 = whole domain ; 3 = surface
	PetscInt         FS;                                // 1 = pointwise gradient
	PetscInt         Gr;                                // 1 = Grad w.r.t solution; 0 = Grad w.r.t to cost function
	PetscInt         OFdef;                             // Objective function defined by hand?
	PetscInt         maxit;                             // maximum number of inverse iteration
	PetscInt         maxitLS;                           // maximum number of backtracking
	PetscScalar      Scale_Grad;                        // scale parameter update with initial gradient?
	PetscScalar      mfitini;   	                    // initial misfit value for current model parameters
	PetscScalar      tol; 	   	                        // tolerance for F/Fini after which code has converged
	PetscScalar      facLS;      	                    // factor in the line search that multiplies current line search parameter if GD update was successful (increases convergence speed)
	PetscScalar      facB;      	                    // backtrack factor that multiplies current line search parameter if GD update was not successful
	PetscScalar      factor2array[51];                  // factor that increases the convergence velocity (this value is added to itself after every successful gradient descent ; only used without tao)
	PetscScalar      maxfac;	                        // limit on the factor (only used without tao)
	PetscScalar      vel_scale;                         // normalization of the observation (currently mean; classically the variance)
	PetscScalar      DII_ref;                           // SUPER UNNECESSARY but DII is otherwise not accesible
	PetscScalar      Coord[3];		                    // Temp Coordinates of comparison points
	PetscScalar      Ax[_MAX_OBS_];                     // X-coordinates of comparison points
	PetscScalar      Ay[_MAX_OBS_];	                    // Y-coordinates of comparison points
	PetscScalar      Az[_MAX_OBS_];                     // Z-coordinates of comparison points
	PetscScalar      Ae[_MAX_OBS_];                     // Velocity target value of comparison points
	PetscInt         Av[_MAX_OBS_];	                    // Velocity components [x/y/z] of comparison points
	PetscScalar      Avel_num[_MAX_OBS_];             	// Numerically computed velocity at the comparison points
	PetscBool        Apoint_on_proc[_MAX_OBS_];         // Is the observation point on the current processor or not (simplified printing)?

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
