// This is a framework code for the adjoint inversion optimization
//---------------------------------------------------------------------------
// COMPUTATION OF ADJOINT INVERSION
//---------------------------------------------------------------------------
// RECIPE:
// Objective function    F(x,x(p)) = (1/2)*[P*(x-x_ini)' * P*(x-x_ini)]      // p = parameter ; x = converged solution ; xini = comparison solution (same size as jr->gsol) ; P = Projection vector containing the proportions of solution influence
// Derivative I          dF/dx     = p*x-P*x_ini
// Adjoint operation     psi       = J^-T * dF/dx           // J = converged Jacobain matrix
// Derivative II         dr/dp     = [r(p+h) - r(p)]/h      // finite difference approximation of derivative of residual r vs parameter
// Gradients             dF/dp     = -psi^T * dr/dp
//
// ------------------------------------------------------------
// USAGE:
// In this file you need to define:
// IOparam.use	                   		= 1	                                // 0 = NO 1 = Free for other inversion methods 2 = Compute gradients 3 = full inversion 4 = save this forward simulation as comparison simulation
// IOparam.Tao                     		= 0                                 // 0 = Self written gradient descent using (BFGS) Hessian 1 = Use PETSCs TAO (LMVM/BLMVM)
// IOparam.mdN                     		= 4                                 // Number of parameters (same as lengths of phs & P )
// IOparam.mdI                     		= 1                                 // Number of indices (same as lengths of Ax & Ay & Az )
// IOparam.Ab                      		= 0                                 // 0 = No usage of bounds 1 = use bounds
// IOparam.Ap                      		= 1                                 // 1 = several indices ; 2 = the whole domain (will exclude mdI & Ax & Ay & Az) ; 3 = the point with maximum velocity
// IOparam.Ax					   		= {1.2}                             // Array (same length as Ay, Az) containing the x coordinates of the point where you want to compute the gradients
// IOparam.Ay					   		= {0.6}                             // Array (same length as Ax, Az) containing the y coordinates of the point where you want to compute the gradients
// IOparam.Az					   		= {0.4}                             // Array (same length as Ax, Ay) containing the z coordinates of the point where you want to compute the gradients
// IOparam.phs                     		= {1 2 1 2}                         // Array (same length as mdN) containing the phase of the parameter
// IOparam.P                       		= {_RHO0_, _RHO0_, _ETA_,_ETA_}     // Array (same length as mdN) containing the parameter corresponding to the phase
// IOparam.Av                      		= {3}                               // Array (same length as Ax, Ay, Az, mdI) containing the related velocity direction in which to compute the gradient
// IOparam.Adv                     		= 1                                 // Should the point be advected?
// IOparam.tol 							= 1e-3;    							// tolerance for F/Fini after which code has converged
// (optional - GD)  Lb                  = {1 0 0 0}                         // Array (same length as mdN) containing the lower bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
// (optional - GD)  Ub                  = {5 8 8 8}                         // Array (same length as mdN) containing the upper bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
// (optional - tao) IOparam.factor1     = 1e1 								// factor to multiply the gradients (should be set such that the highest gradient scales around 1/100 of its parameter)
// (optional - tao) IOparam.factor2 	= 1.5  								// factor that increases the convergence velocity (this value is added to itself after every succesful gradient descent)
// (optional - tao) IOparam.maxfactor2 	= 100 								// limit on the factor2
//
// + provide the characteristic (nondimensional scaling) values that you defined in your input file (to guarantee consistency)
//
// --> COMPILE laMEM with this file as LaMEM.c
// --> EXECUTE your input file (*.dat) (f.e. mpiexec -n 2 /home/user/lamem/bin/opt/LaMEM -ParamFile Input.dat)
//
// ------------------------------------------------------------
// Possible parameters:
//	_RHO0_, _RHON_, _RHOC_,                             // density
//	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
//	_ETA0_,	_E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
//	_BP_, _TAUP_, _GAMMA_, _Q_, _EP_, _VP_,             // Peierls creep
//	_SHEAR_, _BULK_, _KP_,                              // elasticity
//	_COHESION_, _FRICTION_, _CHSOFTID_, _FRSOFTID_,     // plasticity (Drucker-Prager)
//	_ALPHA_, _CP_, _K_, _A_                             // energy
//
// ------------------------------------------------------------
// Possible Velocities:  (Vx   Vy   Vz
//                       (1    2    3)
//
// ------------------------------------------------------------
// LINEAR SOLVER:
// You can control the behaviour of the KSP object for the adjoint with the prefix "as_" (probably the same as "js_")
//
// ------------------------------------------------------------
// FULL INVERSION REMARKS:
// 1) In case you want to perform the full adjoint inversion (ComputeAdjointGradients = 3) make sure that you have a comparison file with a petsc vector the same size as jr->gsol
//    and called Forward_Solution.bin. Most likely you want to run a forward simulation and then change the parameters to do so just run your forward model with ComputeAdjointGradients = 4
//    which automatically saves this file. Perturb the values within the input script and solve again with ComputeAdjointGradients = 3.
//
// 2) You can control the behaviour of the TAO object for the adjoint with the prefix "tao_" (example: '-tao_type lmvm' ; '-tao_fatol 1e-15' ; '-tao_converged_reason')
//
// 3) If not using tao, play around with the factors: factor 1 should be set such that the smallest gradient multiplied by this factor is around 1/100 of its paramter.
//    Size of factor 2 controls how fast the step size will increase (gives faster convergence but you might overstep couple of minimas).
//
// ------------------------------------------------------------
// IMPORTANT REMARKS:
// 1) Since the Adjoint needs the Jacobian matrix for computing the gradients it's crucial to make sure that
//    you compute the Jacobian matrix in the timesteps where you want to compute the gradients
//    (f.e. a linear problem would need a low value for -snes_atol [1e-20] and a low max iteration count -snes_max_it [2] to
//    guarantee the computation of the Jacobian + the option '-snes_type ksponly')
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "paraViewOutBin.h"
#include "paraViewOutSurf.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "multigrid.h"
#include "advect.h"
#include "marker.h"
#include "paraViewOutMark.h"
#include "input.h"
#include "matProps.h"
#include "objFunct.h"
#include "AVDView.h"
#include "break.h"
#include "parsing.h"
#include "adjoint.h"
//---------------------------------------------------------------------------
static char help[] = "Adjoint inversion computation .\n\n";
//--------------------------------------------------------------------------
extern PetscErrorCode PCCreate_SemiRedundant(PC);
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{

	// TODO: 1) Currently nothing

	PetscErrorCode 	ierr;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0, help); CHKERRQ(ierr);

	ModParam        IOparam;
	Scaling         scal;
	PetscScalar     *gradar,*Ubar,*Lbar, *Par, *fcconvar, F;
	Vec             val, Ub, Lb, grad, P;
	PetscInt        i;

	// INITIALIZE
	IOparam.use = 3;   		// 0 = no inversion ; 1 = Tobis inversion ; 2 = only compute adjoint gradients ; 3 = 'full' adjoint inversion with TAO ; 4 = assume this as a forward simulation and save the solution
	IOparam.Tao = 0;   		// Use TAO?
	IOparam.mdN = 3;   		// Number of parameters
	IOparam.mdI = 1;   		// Number of indices
	IOparam.Ab  = 0;   		// Apply bounds?
	IOparam.Ap  = 1;   		// 1 = several indices ; 2 = the whole domain ; 3 = surface
	IOparam.reg = 0;   		// 1 = tikhonov regularization of the cost function (TN) 2 = total variation regularization (TV)
	IOparam.Adv = 0;   		// 1 = advect the point
	IOparam.tol = 1e-15;    // tolerance for F/Fini after which code has converged
	IOparam.factor1 = 1e9; 	// factor to multiply the gradients (should be set such that the highest gradient scales around 1/100 of its parameter ; only used without tao)
	IOparam.factor2 = 1.5;  // factor that increases the convergence velocity (this value is added to itself after every succesful gradient descent ; only used without tao)
	IOparam.maxfactor2 = 100; // limit on the factor (only used without tao)

	IOparam.count = 1;  // iteration counter for the initial cost function

	// VECTORS
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &Lb);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &Ub);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &val);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &P);    CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &grad); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, 1501       , PETSC_DETERMINE, &IOparam.fcconv);   CHKERRQ(ierr);   // 1500 is the maximum inversion iterations that are accepted

	// TEMPORARY VARIABLES
	PetscInt		phsar[IOparam.mdN];
	PetscInt      	typar[IOparam.mdN];
	PetscScalar     W[IOparam.mdN];
	PetscScalar     Ax[IOparam.mdI];
	PetscScalar     Ay[IOparam.mdI];
	PetscScalar     Az[IOparam.mdI];
	PetscScalar     Av[IOparam.mdI];

	// CHARACTERISTIC VALUES (Bounds need to be normalized)
	PetscScalar 	length,viscosity,temperature,stress,time,force,acceleration,mass;


	/*// RAYLEIGH  Boris'Repro
	scal.utype  = _NONE_;//_GEO_;   // _NONE_ or _SI_ or _GEO_

	length      = 1;
	viscosity   = 1;
	temperature = 1;
	stress      = 1;

	time         = viscosity/stress;
	force        = stress*length*length;
	acceleration = length/time/time;
	mass         = force/acceleration;

	scal.inp_mass        = mass;
	scal.inp_time        = time;
	scal.inp_length      = length;
	scal.inp_temperature = temperature;
	scal.inp_force       = force;

	ierr = ScalingCreate(&scal);
*/



	// DETACHEMTN FOLDING NAIARA
	scal.utype  = _GEO_;   // _NONE_ or _SI_ or _GEO_

	length      = 1e3;
	viscosity   = 1e20;
	temperature = 1;
	stress      = 1e6;

	time         = viscosity/stress;
	force        = stress*length*length;
	acceleration = length/time/time;
	mass         = force/acceleration;

	scal.inp_mass        = mass;
	scal.inp_time        = time;
	scal.inp_length      = length;
	scal.inp_temperature = temperature;
	scal.inp_force       = force;

	ierr = ScalingCreate(&scal);


	/*
	// SUBDUCTION
	scal.utype  = _GEO_;   // _NONE_ or _SI_ or _GEO_

	length      = 1e5;
	viscosity   = 1e20;
	temperature = 1;
	stress      = 1e6;

	time         = viscosity/stress;
	force        = stress*length*length;
	acceleration = length/time/time;
	mass         = force/acceleration;

	scal.inp_mass        = mass;
	scal.inp_time        = time;
	scal.inp_length      = length;
	scal.inp_temperature = temperature;
	scal.inp_force       = force;

	ierr = ScalingCreate(&scal);
	*/


	/*////////////////////////////////////////
	//             RAYLEIGH  BORIS REPRO Hihalf             //
	////////////////////////////////////////
	// PHASES
	phsar[0] = 1;
	phsar[1] = 1;
	phsar[2] = 2;
	phsar[3] = 2;
	phsar[4] = 1;
	phsar[5] = 2;
	IOparam.phs = phsar;

	// X-COORDINATE
	Ax[0] = 4.45	/scal.length;    // 4.45
	IOparam.Ax = Ax;

	// Y-COORDINATE
	Ay[0] = 0		/scal.length;
	IOparam.Ay = Ay;

	// Z-COORDINATE
	Az[0] = 0.5	/scal.length;   // 0.5
	IOparam.Az = Az;

	// VELOCITY COMPONENT
	Av[0] = 3;
	// Av[3] = 3;
	IOparam.Av = Av;

	// PARAMETER TYPES
	typar[0] = _RHO0_;
	typar[1] = _ETA0_;
	typar[2] = _RHO0_;
	typar[3] = _ETA0_;
	typar[4] = _N_;
	typar[5] = _N_;
	IOparam.typ = typar;

	// PARAMETER VALUES (only taken into account if IOparam.use = 3)
	VecGetArray(P,&Par);
	Par[0] = 0;    // 1 (initial solution)
	Par[1] = 1;        // 1 (initial solution)
	Par[2] = 1;     // 0 (initial solution)
	Par[3] = 0.1;    // 100 (initial solution)
	Par[4] = 1.3;    // 1.1 (initial solution)
	Par[4] = 1.3;    // 1 (initial solution)
	VecRestoreArray(P,&Par);

	// UPPER BOUND
	VecGetArray(Ub,&Ubar);
	Ubar[0] = 5;
	Ubar[1] = 5;
	Ubar[2] = 1000;
	Ubar[3] = 2;
	Ubar[4] = 2;
	VecRestoreArray(Ub,&Ubar);

	// LOWER BOUND
	VecGetArray(Lb,&Lbar);
	Lbar[0] = 0;
	Lbar[1] = 0;
	Lbar[2] = 0;
	Lbar[3] = 0;
	Lbar[4] = 0;
	VecRestoreArray(Lb,&Lbar);

	// INITIALIZE GRADIENTS
	VecGetArray(grad,&gradar);
	gradar[0] = 0;
	gradar[1] = 0;
	gradar[2] = 0;
	gradar[3] = 0;
	gradar[4] = 0;
	gradar[5] = 0;
	IOparam.grd = gradar;
	VecRestoreArray(grad,&gradar);

	// TIKHONOV WEIGHTS
	W[0] = 1;
	W[1] = 1;
	// gradar[2] = 0;
	// gradar[3] = 0;
	IOparam.W = W;
*/

	/*////////////////////////////////////////
		//             RAYLEIGH  BORIS REPRO Hismall  (linear)           //
		////////////////////////////////////////
		// PHASES
		phsar[0] = 1;
		phsar[1] = 1;
		phsar[2] = 2;
		phsar[3] = 2;
		IOparam.phs = phsar;

		// X-COORDINATE
		Ax[0] = 2.15	/scal.length;    // 2.15
		//Ax[1] = 1.6	/scal.length;    // 4.45
		//Ax[2] = 2.16	/scal.length;    // 4.45
		//Ax[3] = 2.3	/scal.length;    // 4.45
		//Ax[4] = 4.44	/scal.length;    // 4.45
		//Ax[5] = 5.31	/scal.length;    // 4.45
		//Ax[6] = 6.41	/scal.length;    // 4.45
		//Ax[7] = 8.34	/scal.length;    // 4.45
		IOparam.Ax = Ax;

		// Y-COORDINATE
		Ay[0] = 0		/scal.length;
		//Ay[1] = 0		/scal.length;
		//Ay[2] = 0		/scal.length;
		//Ay[3] = 0		/scal.length;
		//Ay[4] = 0		/scal.length;
		//Ay[5] = 0		/scal.length;
		//Ay[6] = 0		/scal.length;
		//Ay[7] = 0		/scal.length;
		IOparam.Ay = Ay;

		// Z-COORDINATE
		Az[0] = 0.0897293	/scal.length;   // 0.5     // 0.105569
		//Az[1] = 0.1	/scal.length;   // 0.5
		//Az[2] = 0.1	/scal.length;   // 0.5
		//Az[3] = 0.1	/scal.length;   // 0.5
		//Az[4] = 0.1	/scal.length;   // 0.5
		//Az[5] = 0.1	/scal.length;   // 0.5
		//Az[6] = 0.1	/scal.length;   // 0.5
		//Az[7] = 0.1	/scal.length;   // 0.5
		IOparam.Az = Az;


		// VELOCITY COMPONENT
		Av[0] = 3;
		//Av[1] = 3;
		//Av[2] = 3;
		//Av[3] = 3;
		//Av[4] = 3;
		//Av[5] = 3;
		//Av[6] = 3;
		//Av[7] = 3;
		IOparam.Av = Av;

		// PARAMETER TYPES
		typar[0] = _RHO0_;
		typar[1] = _ETA_;
		typar[2] = _RHO0_;
		typar[3] = _ETA_;
		IOparam.typ = typar;

		// PARAMETER VALUES (only taken into account if IOparam.use = 3)
		VecGetArray(P,&Par);
		Par[0] = 0;    // 1 (initial solution)
		Par[1] = 1;        // 1 (initial solution)
		Par[2] = 1;     // 0 (initial solution)
		Par[3] = 0.1;    // 100 (initial solution)
		VecRestoreArray(P,&Par);

		// UPPER BOUND
		VecGetArray(Ub,&Ubar);
		Ubar[0] = 5;
		Ubar[1] = 5;
		Ubar[2] = 1000;
		Ubar[3] = 2;
		VecRestoreArray(Ub,&Ubar);

		// LOWER BOUND
		VecGetArray(Lb,&Lbar);
		Lbar[0] = 0;
		Lbar[1] = 0;
		Lbar[2] = 0;
		Lbar[3] = 0;
		VecRestoreArray(Lb,&Lbar);

		// INITIALIZE GRADIENTS
		VecGetArray(grad,&gradar);
		gradar[0] = 0;
		gradar[1] = 0;
		gradar[2] = 0;
		gradar[3] = 0;
		IOparam.grd = gradar;
		VecRestoreArray(grad,&gradar);

		// TIKHONOV WEIGHTS
		W[0] = 1;
		W[1] = 1;
		// gradar[2] = 0;
		// gradar[3] = 0;
		IOparam.W = W;
*/

	////////////////////////////////////////
	//              NAIARA                  //
	////////////////////////////////////////
	// PHASES
	phsar[0] = 1;
	phsar[1] = 2;
	phsar[2] = 3;
	//phsar[3] = 1;
	//phsar[4] = 2;
	//phsar[5] = 3;
	//phsar[6] = 1;
	//phsar[7] = 2;
	//phsar[8] = 3;
	IOparam.phs = phsar;

	// PARAMETER TYPES
	typar[0] = _RHO0_;
	typar[1] = _RHO0_;
	typar[2] = _RHO0_;
	//typar[3] = _ETA0_;
	//typar[4] = _ETA0_;
	//typar[5] = _ETA0_;
	//typar[6] = _N_;
	//typar[7] = _N_;
	//typar[8] = _N_;
	IOparam.typ = typar;

	// PARAMETER VALUES (only taken into account if IOparam.use = 3)
	VecGetArray(P,&Par);
	Par[0] = 1000;    // /scal.density;    //  2900(initial solution)
	Par[1] = 2200;    // /scal.viscosity;;        //  2200(initial solution)
	Par[2] = 2600;   // /scal.viscosity;;        //  1e23(initial solution)      1e19
	//Par[3] = 1e23;   // /scal.density;        //  1e18(initial solution)         1e24
	//Par[4] = 1e21;   // /scal.density;        //  (initial solution)
	//Par[5] = 1e20;   // /scal.viscosity;;        //  (initial solution)
	//Par[6] = 1.2;   // /scal.viscosity;;        //  (initial solution)
	//Par[7] = 1.5;   // /scal.viscosity;;        //  (initial solution)
	//Par[8] = 1.2;   // /scal.viscosity;;        //  (initial solution)
	VecRestoreArray(P,&Par);

	// INITIALIZE GRADIENTS
	VecGetArray(grad,&gradar);
	gradar[0] = 0;
	gradar[1] = 0;
	gradar[2] = 0;
	//gradar[3] = 0;
	//gradar[4] = 0;
	//gradar[5] = 0;
	//gradar[6] = 0;
	//gradar[7] = 0;
	//gradar[8] = 0;
	IOparam.grd = gradar;
	VecRestoreArray(grad,&gradar);

	// UPPER BOUND
	VecGetArray(Ub,&Ubar);
	Ubar[0] = 5000;
	Ubar[1] = 5000;
	Ubar[2] = 5000;
	VecRestoreArray(Ub,&Ubar);

	// LOWER BOUND
	VecGetArray(Lb,&Lbar);
	Lbar[0] = 1000;
	Lbar[1] = 1000;
	Lbar[2] = 1000;
	VecRestoreArray(Lb,&Lbar);

	// X-COORDINATE
	Ax[0] = 20	/scal.length;
	IOparam.Ax = Ax;

	// Y-COORDINATE
	Ay[0] = 0	/scal.length;
	IOparam.Ay = Ay;

	// Z-COORDINATE
	Az[0] = 3	    /scal.length;
	IOparam.Az = Az;

	// VELOCITY COMPONENT
	Av[0] = 3;
	IOparam.Av = Av;

	// TIKHONOV WEIGHTS
	W[0] = 0.000001;
	W[1] = 0.000001;
	W[2] = 0.000001;
	// gradar[2] = 0;
	// gradar[3] = 0;
	IOparam.W = W;

	/*
	////////////////////////////////////////
	//              SUBDUCTION            //
	////////////////////////////////////////
	// PHASES
	phsar[0] = 1;
	phsar[1] = 2;
	phsar[2] = 3;
	phsar[3] = 4;
	phsar[4] = 5;
	phsar[5] = 6;
	phsar[6] = 2;
	phsar[7] = 3;
	phsar[8] = 4;
	phsar[9] = 5;
	phsar[10] = 6;
	phsar[11] = 1;
	phsar[12] = 2;
	phsar[13] = 3;
	phsar[14] = 4;
	phsar[15] = 5;
	phsar[16] = 6;
	IOparam.phs = phsar;

	// PARAMETER TYPES
	typar[0] = _EN_;
	typar[1] = _EN_;
	typar[2] = _EN_;
	typar[3] = _EN_;
	typar[4] = _EN_;
	typar[5] = _EN_;
	typar[6] = _N_;
	typar[7] = _N_;
	typar[8] = _N_;
	typar[9] = _N_;
	typar[10] = _N_;
	typar[11] = _RHO0_;
	typar[12] = _RHO0_;
	typar[13] = _RHO0_;
	typar[14] = _RHO0_;
	typar[15] = _RHO0_;
	typar[16] = _RHO0_;
	IOparam.typ = typar;

	// PARAMETER VALUES (only taken into account if IOparam.use = 3)
	VecGetArray(P,&Par);
	Par[0] = 350e3;    // /scal.density;    //  2900(initial solution)
	Par[1] = 350e3;    // /scal.viscosity;;        //  2200(initial solution)
	Par[2] = 350e3;    // /scal.density;        //  2500(initial solution)
	Par[3] = 350e3;   // /scal.viscosity;;        //  1e23(initial solution)      1e19
	Par[4] = 350e3;   // /scal.viscosity;;        //  1e23(initial solution)      1e19
	Par[5] = 350e3;   // /scal.viscosity;;        //  1e23(initial solution)      1e19
	VecRestoreArray(P,&Par);

	// INITIALIZE GRADIENTS
	VecGetArray(grad,&gradar);
	gradar[0] = 0;
	gradar[1] = 0;
	gradar[2] = 0;
	gradar[3] = 0;
	gradar[4] = 0;
	gradar[5] = 0;
	gradar[6] = 0;
	gradar[7] = 0;
	gradar[8] = 0;
	gradar[9] = 0;
	gradar[10] = 0;
	gradar[11] = 0;
	gradar[12] = 0;
	gradar[13] = 0;
	gradar[14] = 0;
	gradar[15] = 0;
	gradar[16] = 0;
	IOparam.grd = gradar;
	VecRestoreArray(grad,&gradar);

	// UPPER BOUND
	VecGetArray(Ub,&Ubar);
	Ubar[0] = 800e3;
	Ubar[1] = 800e3;
	Ubar[2] = 800e3;
	Ubar[3] = 800e3;
	Ubar[4] = 800e3;
	Ubar[5] = 800e3;
	VecRestoreArray(Ub,&Ubar);

	// LOWER BOUND
	VecGetArray(Lb,&Lbar);
	Lbar[0] = 100e3;
	Lbar[1] = 100e3;
	Lbar[2] = 100e3;
	Lbar[3] = 100e3;
	Lbar[4] = 100e3;
	Lbar[5] = 100e3;
	VecRestoreArray(Lb,&Lbar);

	// X-COORDINATE
	Ax[0] = 400	/scal.length;
	IOparam.Ax = Ax;

	// Y-COORDINATE
	Ay[0] = 99	/scal.length;
	IOparam.Ay = Ay;

	// Z-COORDINATE
	Az[0] = 300	    /scal.length;
	IOparam.Az = Az;

	// VELOCITY COMPONENT
	Av[0] = 3;
	IOparam.Av = Av;
	*/



	/*// NORMALIZE
	VecGetArray(P,&Par);
	Par[0] = (Par[0]-2000)/(3200-2000);
	Par[1] = (Par[1]-2000)/(3200-2000);
	Par[2] = (Par[2]-2000)/(3200-2000);
	Par[3] = (Par[3]-2000)/(3200-2000);
	Par[4] = (Par[4]-2000)/(3200-2000);
	Par[5] = (Par[5]-2000)/(3200-2000);
	Par[6] = (Par[6]-2000)/(3200-2000);
	Par[7] = (Par[7]-2000)/(3200-2000);
	Par[8] = (Par[8]-2000)/(3200-2000);
	Par[9] = (Par[9]-1e10)/(1e20-1e10);
	Par[10] = (Par[10]-1e10)/(1e20-1e10);
	VecRestoreArray(P,&Par);*/
	/*// DENORMALIZE
	VecGetArray(P,&Par);
	Par[0] = (Par[0]*(3200-2000))+2000;
	Par[1] = (Par[1]*(3200-2000))+2000;
	Par[2] = (Par[2]*(3200-2000))+2000;
	Par[3] = (Par[3]*(3200-2000))+2000;
	Par[4] = (Par[4]*(3200-2000))+2000;
	Par[5] = (Par[5]*(3200-2000))+2000;
	Par[6] = (Par[6]*(3200-2000))+2000;
	Par[7] = (Par[7]*(3200-2000))+2000;
	Par[8] = (Par[8]*(3200-2000))+2000;
	Par[9] = (Par[9]*(1e20-1e10))+1e10;
	Par[10] = (Par[10]*(1e20-1e10))+1e10;
	VecRestoreArray(P,&Par);*/

	////
	ierr = PCRegister("pc_semiredundant",PCCreate_SemiRedundant);

	//===============
	// SOLVE ADJOINT
	//===============
	// only compute the adjoint gradients or simply forward code
 	if(IOparam.use == 2 || IOparam.use == 0)
 	{
 		VecDuplicate(P,&IOparam.P);
 		VecCopy(P,IOparam.P);

 		// call LaMEM main library function
 		ierr = LaMEMLib(&IOparam); CHKERRQ(ierr);
 	}
 	// compute 'full' adjoint inversion
 	else if(IOparam.use == 3)
 	{

 		VecDuplicate(P,&IOparam.P);
 		VecCopy(P,IOparam.P);

 		// if tao is used try the LMVM/BLMVM algorithms
 		if(IOparam.Tao == 1)
 		{
 	 		Tao tao;

 	 		ierr = TaoCreate(PETSC_COMM_WORLD,&tao);

 	 	 	// 1.Check if bounds are available and sets them
 	 		if (IOparam.Ab == 1)
 	 		{
 	 	 	 	ierr = TaoSetVariableBounds(tao,Lb,Ub);	 								CHKERRQ(ierr);
 	 	 	 	ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr);
 	 		}
 	 		else
 	 		{
 	 			ierr = TaoSetType(tao,TAOLMVM);CHKERRQ(ierr);                    // TAOLMVM, TAOBLMVM, TAOBMRM (bad), TAOCG (all 4), TAOTRON (might crash but is fast)
 	 		}

 	 		// 2. Set up Tao
 	 	 	ierr = TaoSetObjectiveAndGradientRoutine(tao, AdjointOptimisationTAO, &IOparam);	 	CHKERRQ(ierr);
 	 	 	ierr = TaoSetInitialVector(tao,P);	 									CHKERRQ(ierr);
 	 	 	ierr = TaoSetTolerances(tao,1e-30,1e-30,1e-30,1e-30,1e-30);	CHKERRQ(ierr);
 	 	 	ierr = TaoSetFunctionLowerBound(tao,0.0001);CHKERRQ(ierr);
 	 	 	ierr = TaoSetFromOptions(tao);	 										CHKERRQ(ierr);

 	 	 	// 3. Solve Tao & view result
 	 	 	ierr = TaoSolve(tao);	 												CHKERRQ(ierr);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 	 	 	TaoView(tao,PETSC_VIEWER_STDOUT_WORLD);
 	 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");

 	 	 	// 4. Clean
 	 	 	ierr = TaoDestroy(&tao);
 		}
 		else
 		// Without TAO try line search tuned gradient descent
 		{
 			ierr = AdjointOptimisation(P, F, grad, &IOparam);
 		}

 		PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------\n");
 		PetscPrintf(PETSC_COMM_WORLD,"*         INVERSION RESULT SUMMARY       *\n");
 		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n");
 		PetscPrintf(PETSC_COMM_WORLD,"Number of inversion iterations: %d\n",IOparam.count);
 		PetscPrintf(PETSC_COMM_WORLD,"F/Fini:\n");
 		VecGetArray(IOparam.fcconv,&fcconvar);
 		for(i=0;i<IOparam.count;i++)
 		{
 			PetscPrintf(PETSC_COMM_WORLD,"%8.15f\n",fcconvar[i]);
 		}
 		VecRestoreArray(IOparam.fcconv,&fcconvar);
 		PetscPrintf(PETSC_COMM_WORLD,"\nFinal cost function:\n");
 		PetscPrintf(PETSC_COMM_WORLD,"%8.15f\n",IOparam.mfit);
 		PetscPrintf(PETSC_COMM_WORLD,"\nFinal Parameters:\n");
		VecGetArray(IOparam.P,&Par);
		for(i=0;i<IOparam.mdN;i++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%.12f\n",Par[i]);
		}
		VecRestoreArray(IOparam.P,&Par);
 		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n\n");

 	}
 	// this is a forward simulation that we want to save as comparison solution
 	else if(IOparam.use == 4)
 	{
 		// call LaMEM main library function
 		ierr = LaMEMLib(&IOparam); CHKERRQ(ierr);

 		// Save output
 		PetscViewer     viewer;
 		PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
 		PetscViewerSetType(viewer,PETSCVIEWERBINARY);
 		PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
 		PetscViewerFileSetName(viewer,"Forward_Solution.bin");
 	 	VecView(IOparam.xini,viewer);
 	 	PetscViewerDestroy(&viewer);
 	 	PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n        Forward Solution succesfully saved\n------------------------------------------\n");
 	}

	ierr = PetscMemzero(&IOparam, sizeof(ModParam)); CHKERRQ(ierr);

	// cleanup PETSC
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointOptimisation"
PetscErrorCode AdjointOptimisation(Vec P, PetscScalar F, Vec grad, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	PetscInt 		i, j, a, LScount;
	PetscScalar 	*Par, *Paroldar, *gradoldar, *Scalar, *gradar, tempdot1, *rtemp, *dgradtemp, *dPtemp, *temptemp1, *fcconvar, b;
	PetscScalar   	val, Fold;
	ModParam    	*IOparam;
	IOparam     	= (ModParam*)ctx;
	Vec         	dP,dgrad,Pold,gradold,r,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10;
	double 	    	H[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp1[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp2[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp3[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp4[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp5[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp6[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp7[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp8[IOparam->mdN-1][IOparam->mdN-1];
	double			Htemp9[IOparam->mdN-1][IOparam->mdN-1];

	// get parameter values
	VecDuplicate(IOparam->P,&dP);
	VecDuplicate(IOparam->P,&Pold);
	VecDuplicate(grad,&gradold);
	VecDuplicate(grad,&dgrad);
	VecDuplicate(grad,&r);
	VecDuplicate(grad,&temp1);
	VecDuplicate(grad,&temp2);
	VecDuplicate(grad,&temp3);
	VecDuplicate(grad,&temp4);
	VecDuplicate(grad,&temp5);
	VecDuplicate(grad,&temp6);
	VecDuplicate(grad,&temp7);
	VecDuplicate(grad,&temp8);
	VecDuplicate(grad,&temp9);
	VecDuplicate(grad,&temp10);
	VecCopy(P,IOparam->P);

	// Initialize parameters
	a = -1;
	b = 1;   // initial value for factor2

	// Initialize cost functions
	F = 1e100;
	Fold = 1e100;

	while(F>IOparam->tol)
	{
		// Give the updated values to the code
		VecCopy(P,IOparam->P);

		// Reset line search counter
		LScount = 1;

		// call LaMEM main library function
		ierr = LaMEMLib(IOparam); CHKERRQ(ierr);

		// Save intial cost function & create initial Hessian
		if(IOparam->count==1)
		{
			IOparam->mfitini = IOparam->mfit;
			for(i=0;i<IOparam->mdN;i++)
			{
				H[i][i] = 1;
			}
		}

		// Save cost function
		F = IOparam->mfit;

		// If cost function in this timestep is larger then before perform bisection line search
		while(F>Fold)
		{
			b = 1;
			PetscPrintf(PETSC_COMM_WORLD,"\n- - - - - - - - - - - - - - - - - - - - - - - - - - - \n");
			PetscPrintf(PETSC_COMM_WORLD,"              LINE SEARCH IT %d                       \n",LScount);

			VecGetArray(P,&Par);
			VecGetArray(Pold,&Paroldar);
			VecGetArray(dP,&dPtemp);
			for(i=0;i<IOparam->mdN;i++)
			{
				for(j=0;j<IOparam->mdN;j++)
				{
					dPtemp[i] = dPtemp[i]/2;
				}
			}

			// Update parameter
			for(i=0;i<IOparam->mdN;i++)
			{
				Par[i] = Paroldar[i] + dPtemp[i];
			}
			VecRestoreArray(P,&Par);
			VecRestoreArray(Pold,&Paroldar);
			VecRestoreArray(dP,&dPtemp);

			// Give the updated values to the code
			VecCopy(P,IOparam->P);

			// call LaMEM main library function
			ierr = LaMEMLib(IOparam); CHKERRQ(ierr);

			F = IOparam->mfit;

			LScount+=1;
			if(LScount>10)
			{
				PetscPrintf(PETSC_COMM_WORLD,"******************************************************\n");
				PetscPrintf(PETSC_COMM_WORLD,"*             YOUR SOLUTION DIVERGED                 *\n");
				PetscPrintf(PETSC_COMM_WORLD,"******************************************************\n\n");

				// Return parameters for final output
				VecCopy(P,IOparam->P);

				PetscFunctionReturn(0);
			}
		}

		// Zero out perturbation
		VecDuplicate(IOparam->P,&dP);

		// restore parameter values
		VecDuplicate(IOparam->P,&P);
		VecCopy(IOparam->P,P);

		VecGetArray(grad,&gradar);
		for(j = 0; j < IOparam->mdN; j++)
		{
			gradar[j] = IOparam->grd[j];
		}
		VecRestoreArray(grad,&gradar);

/*
// -----------------------------------------------------------------//
//              BFGS HESSIAN                                        //
// -----------------------------------------------------------------//

		VecCopy(Pold,dP);
		VecGetArray(Pold,&Paroldar);
		VecGetArray(P,&Par);
		VecGetArray(dP,&dPtemp);
		for(i=0;i<IOparam->mdN;i++)
		{
			dPtemp[i] = Par[i]-Paroldar[i];
			if(dPtemp[i]<1e-16)
			{
				dPtemp[i] = 1;
			}
		}
		VecRestoreArray(dP,&dPtemp);
		VecRestoreArray(Pold,&Paroldar);
		VecRestoreArray(P,&Par);

		VecCopy(gradold,dgrad);
		VecGetArray(gradold,&gradoldar);
		VecGetArray(grad,&gradar);
		VecGetArray(dgrad,&dgradtemp);
		for(i=0;i<IOparam->mdN;i++)
		{
			dgradtemp[i] = gradar[i]-gradoldar[i];
			if(dgradtemp[i]<1e-16)
			{
				dgradtemp[i] = 1;
			}
		}
		VecRestoreArray(dgrad,&dgradtemp);
		VecRestoreArray(gradold,&gradoldar);
		VecRestoreArray(grad,&gradar);

		// Right part of minus sign
		VecGetArray(dgrad,&dgradtemp);
		VecGetArray(temp1,&temptemp1);
		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				temptemp1[i] = temptemp1[i] + H[i][j]*dgradtemp[i];
			}
		}
		VecRestoreArray(temp1,&temptemp1);
		VecRestoreArray(dgrad,&dgradtemp);

		ierr =  VecDot(temp1,dgrad,&tempdot1);
		ierr =  VecSet(temp3,tempdot1);
		ierr =  VecPointwiseDivide(temp4,temp1,temp3);

		// Left part of minus sign
		ierr = VecPointwiseMult(temp5,dP,dgrad);
		ierr = VecPointwiseDivide(temp6,dP,temp5);

		//r =
		ierr =  VecCopy(temp4,r);
		ierr =  VecAYPX(r,-1,temp6);

		VecView(dP,PETSC_VIEWER_STDOUT_WORLD);


			///////////////////////////////////////////////////////////////////////////////////
			//                        r
			//         temp6                                 temp4(r)
			//                                                           tempdot1(temp3)
			//            temp5                       temp1         temp1
			// r = (dp./(dp.*dgrad')) - ((NUM.Adjoint.H*dgrad')/(dgrad*NUM.Adjoint.H*dgrad'));
			///////////////////////////////////////////////////////////////////////////////////


		VecGetArray(r,&rtemp);
		VecGetArray(dP,&dPtemp);
		VecGetArray(dgrad,&dgradtemp);
		VecGetArray(temp1,&temptemp1);
		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				Htemp1[i][j] = dPtemp[i]*dPtemp[j];
				Htemp2[i][j] = dPtemp[i]*dgradtemp[j];
				Htemp3[i][j] = rtemp[i]*rtemp[j];
				Htemp4[i][j] = temptemp1[i]*temptemp1[j];
			}
		}
		VecRestoreArray(r,&rtemp);
		VecRestoreArray(dP,&dPtemp);
		VecRestoreArray(dgrad,&dgradtemp);
		VecRestoreArray(temp1,&temptemp1);

		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				Htemp5[i][j] = Htemp4[i][j]/tempdot1;
				Htemp6[i][j] = Htemp1[i][j]/Htemp2[i][j];
				Htemp7[i][j] = Htemp3[i][j]*tempdot1;
			}
		}

		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				Htemp8[i][j] = Htemp5[i][j]+Htemp6[i][j];
			}
		}

		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				Htemp9[i][j] = Htemp8[i][j]+Htemp7[i][j];
			}
		}

		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				H[i][j] = H[i][j]-Htemp9[i][j];
			}
		}

			///////////////////////////////////////////////////////////////////////////////////
			// Hessian inverse
			//                            H
			//                                                                                                                                          Htemp9
			//                                                                                                            Htemp8
			//                                                                            Htemp5                                     Htemp6                                Htemp7
			//                                                   Htemp4
			//										temp1						temp1					tempdot1				Htemp1		Htemp2				tempdot1			      Htemp3
			//NUM.Adjoint.H = NUM.Adjoint.H - ((NUM.Adjoint.H*dgrad'*(NUM.Adjoint.H*dgrad')')/(dgrad*NUM.Adjoint.H*dgrad')) + ((dp*dp')  /  (dp'*dgrad')) + (dgrad*NUM.Adjoint.H*dgrad'  *(r*r'));
			///////////////////////////////////////////////////////////////////////////////////

		PetscPrintf(PETSC_COMM_WORLD,"\nHESSIAN MATRIX:\n");
		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				PetscPrintf(PETSC_COMM_WORLD,"%3.6f  ",H[i][j]);
			}
			PetscPrintf(PETSC_COMM_WORLD,"\n");
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n");
*/


		PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"%d. IT INVERSION RESULT: line search its = %d ; F / FINI = %.12f\n\n",IOparam->count,LScount-1,IOparam->mfit/IOparam->mfitini);
		PetscPrintf(PETSC_COMM_WORLD,"FOLD = %.20f \n   F = %.20f\n\n",Fold,F);

		// BEFORE UPDATING the par vector store the old gradient & Parameter vector
		VecCopy(P,Pold);
		VecCopy(grad,gradold);
		Fold = F;

		VecGetArray(grad,&gradar);
		VecGetArray(P,&Par);
		VecGetArray(dP,&dPtemp);
		for(i=0;i<IOparam->mdN;i++)
		{
			for(j=0;j<IOparam->mdN;j++)
			{
				dPtemp[i] = dPtemp[i] + (H[i][j]*gradar[i]);
			}
		}

		VecRestoreArray(dP,&dPtemp);
		PetscReal min;
		ierr =  VecMin(dP,NULL,&min);
		if(min<0)
		{
			min = -min;
		}
		VecGetArray(dP,&dPtemp);

		PetscPrintf(PETSC_COMM_WORLD,"Factor2 = %2.1f\n\n",b);
		for(i=0;i<IOparam->mdN;i++)
		{
			dPtemp[i] = b*dPtemp[i]*a;
		}

		b = b+IOparam->factor2;
		if(b>IOparam->maxfactor2)
		{
			b = IOparam->maxfactor2;
		}

		// Display the current state of the parameters
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%D. Diff parameter value = %.12f\n",j+1,dPtemp[j]);
		}

		PetscPrintf(PETSC_COMM_WORLD,"\n");

		// Update parameter
		for(i=0;i<IOparam->mdN;i++)
		{
			Par[i] = Par[i] + dPtemp[i];
		}
		VecRestoreArray(grad,&gradar);
		VecRestoreArray(P,&Par);
		VecRestoreArray(dP,&dPtemp);

		// Display the current state of the parameters
		VecGetArray(P,&Par);
		for(j = 0; j < IOparam->mdN; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD,"%D. Parameter value = %.12f\n",j+1,Par[j]);
		}
		VecRestoreArray(P,&Par);

		PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------\n\n");

		// Give the updated values to the code  (actually unfortunately necessary here and at the top of this function - need to rearrange that)
		VecCopy(P,IOparam->P);

		VecGetArray(IOparam->fcconv,&fcconvar);
		fcconvar[IOparam->count] = IOparam->mfit/IOparam->mfitini;
		VecRestoreArray(IOparam->fcconv,&fcconvar);

		// count
		IOparam->count += 1;
		if(IOparam->count>1500)
		{
			PetscPrintf(PETSC_COMM_WORLD,"\n\n\nEXCEEDED 1500 FUNCTION EVALUATIONS (consider changing inversion options)\n\n\n");
			PetscFunctionReturn(0);
		}

	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "AdjointOptimisationTAO"
PetscErrorCode AdjointOptimisationTAO(Tao tao, Vec P, PetscReal *F, Vec grad, void *ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize
	PetscInt j;
	PetscScalar *Par, *gradar;
	PetscReal   norm;
	ModParam    *IOparam;
	IOparam     = (ModParam*)ctx;

	// get parameter values
	VecDuplicate(P,&IOparam->P);
	VecCopy(P,IOparam->P);

	// call LaMEM main library function
	ierr = LaMEMLib(IOparam); CHKERRQ(ierr);

	// restore parameter values
	VecDuplicate(IOparam->P,&P);
	VecCopy(IOparam->P,P);

	// Store the gradient & misfit
	VecGetArray(grad,&gradar);
	for(j = 0; j < IOparam->mdN; j++)
	{
		gradar[j] = IOparam->grd[j];
	}
	VecRestoreArray(grad,&gradar);

	*F = IOparam->mfit;

	// Save intial cost function
	if(IOparam->count==1){IOparam->mfitini = IOparam->mfit;}

	// Display the current state of the parameters
	VecGetArray(P,&Par);
	for(j = 0; j < IOparam->mdN; j++)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%D. Parameter value = %.12f\n",j+1,Par[j]);
	}
	VecRestoreArray(P,&Par);

	// Relative cost function
	PetscPrintf(PETSC_COMM_WORLD,"mfit / mfit0 = %.12f\n------------------------------------------\n\n",IOparam->mfit/IOparam->mfitini);

	// count
	IOparam->count += 1;
	if(IOparam->count>1500)
	{
		PetscPrintf(PETSC_COMM_WORLD,"\n\n\nEXCEEDED 1500 FUNCTION EVALUATIONS (consider changing inversion options)\n\n\n");
		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}


