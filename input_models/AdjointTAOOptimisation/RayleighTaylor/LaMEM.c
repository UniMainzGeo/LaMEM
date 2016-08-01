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
// IOparam.use	                   = 1	                               // 0 = NO 1 = Tobi's inversion 2 = Compute gradients 3 = full inversion 4 = save this forward simulation as comparison simulation
// IOparam.mdN                     = 4                                 // Number of parameters (same as lengths of phs & P )
// IOparam.mdI                     = 1                                 // Number of indices (same as lengths of Ax & Ay & Az )
// IOparam.Ab                      = 0                                 // 0 = No usage of bounds 1 = use bounds
// IOparam.Ap                      = 1                                 // 1 = several indices ; 2 = the whole domain (will exclude mdI & Ax & Ay & Az)
// IOparam.Ax					   = {1.2}                             // Array (same length as Ay, Az) containing the x coordinates of the point where you want to compute the gradients
// IOparam.Ay					   = {0.6}                             // Array (same length as Ax, Az) containing the y coordinates of the point where you want to compute the gradients
// IOparam.Az					   = {0.4}                             // Array (same length as Ax, Ay) containing the z coordinates of the point where you want to compute the gradients
// IOparam.phs                     = {1 2 1 2}                         // Array (same length as mdN) containing the phase of the parameter
// IOparam.P                       = {_RHO0_, _RHO0_, _ETA_,_ETA_}     // Array (same length as mdN) containing the parameter corresponding to the phase
// IOparam.Av                      = {3}                               // Array (same length as Ax, Ay, Az, mdI) containing the related velocity direction in which to compute the gradient
// (optional) Lb                   = {1 0 0 0}                         // Array (same length as mdN) containing the lower bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
// (optional) Ub                   = {5 8 8 8}                         // Array (same length as mdN) containing the upper bounds for the parameters (only taken into account if 'IOparam.Ab  = 1;')
//
// + provide the characteristic values that you defined in your input file (to guarantee consistency)
//
// --> COMPILE laMEM with this file as LaMEM.c
// --> EXECUTE your input file (*.dat) (f.e. mpiexec -n 2 /home/greuber/Codes/lamem/bin/opt/LaMEM -tao_fatol 1e-15 -ParamFile Input.dat)
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

	// TODO: 1) Include whole phases as indices
	//       2) Eventually put the perturbation as par=par+(perturb*par)

	PetscErrorCode 	ierr;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0, help); CHKERRQ(ierr);

	ModParam        IOparam;
	Scaling         scal;
	PetscScalar     *gradar,*Ubar,*Lbar, *Par;
	Vec             val, Ub, Lb, grad, P;

	// INITIALIZE
	IOparam.use = 3;   // 0 = no inversion ; 1 = Tobis inversion ; 2 = only compute adjoint gradients ; 3 = 'full' adjoint inversion with TAO ; 4 = assume this as a forward simulation and save the solution
	IOparam.mdN = 2;   // Number of parameters
	IOparam.mdI = 1;   // Number of indices
	IOparam.Ab  = 1;   // Apply bounds?
	IOparam.Ap  = 1;   // 1 = several indices ; 2 = the whole domain
	IOparam.reg = 0;   // 1 = tikhonov regularization of the cost function (TN) 2 = total variation regularization (TV)

	IOparam.count = 1;  // iteration counter for the initial cost function

	// VECTORS
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &Lb);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &Ub);   CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &val);  CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &P);    CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, IOparam.mdN, PETSC_DETERMINE, &grad); CHKERRQ(ierr);

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

	// RAYLEIGH
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

	////////////////////////////////////////
	//             RAYLEIGH               //
	////////////////////////////////////////
	// PHASES
	phsar[0] = 1;
	phsar[1] = 2;
	IOparam.phs = phsar;

	// X-COORDINATE
	Ax[0] = 1.49	/scal.length;
	IOparam.Ax = Ax;

	// Y-COORDINATE
	Ay[0] = 0		/scal.length;
	IOparam.Ay = Ay;

	// Z-COORDINATE
	Az[0] = 0.5099  /scal.length;
	IOparam.Az = Az;

	// VELOCITY COMPONENT
	Av[0] = 3;
	IOparam.Av = Av;

	// PARAMETER TYPES
	typar[0] = _RHO0_;
	typar[1] = _ETA_;
	IOparam.typ = typar;

	// PARAMETER VALUES (only taken into account if IOparam.use = 3)
	VecGetArray(P,&Par);
	Par[0] = 0.00000001;    // 3 (initial solution)
	Par[1] = 0.001;        // 5 (initial solution)
	VecRestoreArray(P,&Par);

	// UPPER BOUND
	VecGetArray(Ub,&Ubar);
	Ubar[0] = 10;
	Ubar[1] = 10;
	VecRestoreArray(Ub,&Ubar);

	// LOWER BOUND
	VecGetArray(Lb,&Lbar);
	Lbar[0] = 0;
	Lbar[1] = 0;
	VecRestoreArray(Lb,&Lbar);

	// INITIALIZE GRADIENTS
	VecGetArray(grad,&gradar);
	gradar[0] = 0;
	gradar[1] = 0;
	IOparam.grd = gradar;
	VecRestoreArray(grad,&gradar);

	// TIKHONOV WEIGHTS
	W[0] = 1;
	W[1] = 1;
	IOparam.W = W;


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
 	 	ierr = TaoSetObjectiveAndGradientRoutine(tao, AdjointOptimisation, &IOparam);	 	CHKERRQ(ierr);
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
PetscErrorCode AdjointOptimisation(Tao tao, Vec P, PetscReal *F, Vec grad, void *ctx)
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

	PetscFunctionReturn(0);
}
