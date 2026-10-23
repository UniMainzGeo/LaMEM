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
//.....................   OBJECTIVE FUNCTION ROUTINES   .....................
//---------------------------------------------------------------------------
#ifndef __objFunct_h__
#define __objFunct_h__
//---------------------------------------------------------------------------

struct FB;

//-----------------------------------------------------------------------------
enum PTypes// List of model parameter types (30)
{
	// -- material model parameter types --
	_RHO0_, _RHON_, _RHOC_,                             // density
	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
	_ETA0_, _E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
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

//-----------------------------------------------------------------------------
// Structure that holds inversion parameters
struct ModParam
{
	InvTypes         use;                               // Choose one of InvTypes
	PetscInt         mdN;                               // number of model parameters
	PetscInt         mID;                               // current model number
	char             type_name[_MAX_PAR_][_str_len_];   // stores the name of the adjoint parameters
	PetscInt         FD_gradient[_MAX_PAR_];            // Compute gradient via (brute force) finite differences
	PetscScalar      FD_eps[_MAX_PAR_];                 // Perturbation for this parameter if using (brute force) FD; can override default value
	PetscInt         phs[_MAX_PAR_];                    // phase of the parameter
	PetscInt         par_log10[_MAX_PAR_];              // Is the value indicated log10(par)
	PetscScalar      grd[_MAX_PAR_];                    // gradient value
	PetscBool        BruteForce_FD;                     // indicate whether we compute Brute force FD or not

	PetscScalar     *val;                               // model value
	PetscScalar      mfit, mfitCenter;                  // misfit value for current model parameters
	DBMat            dbm_modified;                      // holds the (modified) LaMEM material database
	FB              *fb;                                // holds a copy of the filebuffer
	PetscScalar      ReferenceDensity;                  // Reference density (perturbations are computed w.r.t. this value)

	// Variables additionally needed for the adjoint TAO solver
	Vec              xini;                              // Comparison velocity field for adjoint inversion
	Vec              P;                                 // vector containing parameters
	Vec              fcconv;                            // Vector containing all f/fini values to track convergence
	PetscInt         Ab;                                // Use adjoint bounds (only works with Tao)?
	PetscInt         Tao;                               // Use TAO?
	PetscInt         ScalLaws;                          // Print scaling laws?
	PetscInt         SetInitAdjParam;                   // Use AdjointParameterStart/End InitialGuess as initial parameter or not?
	PetscInt         Adv;                               // Advect the point?
	PetscInt         count;                             // iteration counter
	PetscInt         SCF;                               // Scale cost function?
	PetscInt         mdI;                               // number of indices
	PetscInt         Ap;                                // 1 = several indices ; 2 = whole domain ; 3 = surface
	PetscInt         Gr;                                // 1 = Grad w.r.t solution; 0 = Grad w.r.t to cost function
	PetscInt         OFdef;                             // Objective function defined by hand?
	PetscInt         maxit;                             // maximum number of inverse iteration
	PetscInt         maxitLS;                           // maximum number of backtracking
	PetscInt         MfitType;                          // What observable? 0 = Vel; 1 = PSD
	PetscScalar      Scale_Grad;                        // scale parameter update with initial gradient?
	PetscScalar      mfitini;                           // initial misfit value for current model parameters
	PetscScalar      tol;                               // tolerance for F/Fini after which code has converged
	PetscScalar      facLS;                             // factor in the line search that multiplies current line search parameter if GD update was successful (increases convergence speed)
	PetscScalar      facB;                              // backtrack factor that multiplies current line search parameter if GD update was not successful
	PetscScalar      factor2array[51];                  // factor that increases the convergence velocity (this value is added to itself after every successful gradient descent ; only used without tao)
	PetscScalar      maxfac;                            // limit on the factor (only used without tao)
	PetscScalar      vel_scale;                         // normalization of the observation (currently mean; classically the variance)
	PetscScalar      Coord[3];                          // Temp Coordinates of comparison points
	PetscScalar      Ax[_MAX_OBS_];                     // X-coordinates of comparison points
	PetscScalar      Ay[_MAX_OBS_];                     // Y-coordinates of comparison points
	PetscScalar      Az[_MAX_OBS_];                     // Z-coordinates of comparison points
	PetscScalar      Ae[_MAX_OBS_];                     // Velocity target value of comparison points
	PetscInt         Av[_MAX_OBS_];                     // Velocity components [x/y/z] of comparison points
	char             ObsName[_MAX_OBS_][5];             // Type of the current observation we compare with
	PetscScalar      Avel_num[_MAX_OBS_];               // Numerically computed velocity at the comparison points
	PetscBool        Apoint_on_proc[_MAX_OBS_];         // Is the observation point on the current processor or not (simplified printing)?
	char             ScalLawFilename[_str_len_];        // Name of scaling law file
};

//---------------------------------------------------------------------------

#endif
