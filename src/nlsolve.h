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
//.....................   NONLINEAR SOLVER ROUTINES   .......................
//---------------------------------------------------------------------------
#ifndef __nlsolve_h__
#define __nlsolve_h__
//---------------------------------------------------------------------------
/*

http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESLineSearchSetFromOptions.html
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESLineSearchPreCheckPicard.html#SNESLineSearchPreCheckPicard

	-snes_linesearch_type <type> 	- basic, bt, l2, cp, shell
	-snes_linesearch_order <order> 	- 1, 2, 3. Most types only support certain orders (bt supports 2 or 3)
	-snes_linesearch_norms 	- Turn on/off the linesearch norms for the basic linesearch type
	-snes_linesearch_minlambda 	- The minimum step length
	-snes_linesearch_maxstep 	- The maximum step size
	-snes_linesearch_rtol 	- Relative tolerance for iterative line searches
	-snes_linesearch_atol 	- Absolute tolerance for iterative line searches
	-snes_linesearch_ltol 	- Change in lambda tolerance for iterative line searches
	-snes_linesearch_max_it 	- The number of iterations for iterative line searches
	-snes_linesearch_monitor 	- Print progress of line searches
	-snes_linesearch_damping 	- The linesearch damping parameter
	-snes_linesearch_keeplambda 	- Keep the previous search length as the initial guess.
	-snes_linesearch_precheck_picard 	- Use precheck that speeds up convergence of picard method
	-snes_linesearch_precheck_picard_angle 	- Angle used in picard precheck method

http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateMFFD.html#MatCreateMFFD

-mat_mffd_type 	- wp or ds (see MATMFFD_WP or MATMFFD_DS)
-mat_mffd_err 	- square root of estimated relative error in function evaluation
-mat_mffd_period -how often h is recomputed, defaults to 1, everytime

-mat_mffd_err <error_rel> 	- Sets error_rel
-mat_mffd_unim <umin> 	- Sets umin (for default PETSc routine that computes h only)
-mat_mffd_check_positivity	-

 */

//---------------------------------------------------------------------------

// Jacobian type
enum JacType
{
	//===================
	// assembled matrices
	//===================
	_PICARD_,   // constant effective coefficients approximation (viscosity, conductivity, stress)
//	_FDCOLOR_,  // finite difference coloring approximation with full sparsity pattern
//	_ANALYTIC_, // analytic Jacobian with full sparsity pattern
//	_APPROX_,   // analytic Jacobian truncated to Picard sparsity pattern (possibly with diagonal compensation)
//	_FDAPPROX_, // finite difference coloring approximation truncated to Picard sparsity pattern
	//============
	// matrix-free
	//============
	_MFFD_ // built-in finite difference approximation

};

//---------------------------------------------------------------------------
struct NLSol
{
	Mat       J;      // Jacobian matrix
	Mat       P;      // preconditioner
	Mat       MFFD;   // matrix-free finite difference Jacobian
	PCStokes  pc;     // Stokes preconditioner

	JacType     jtype;    // actual type of Jacobian operator
	PetscInt    it;       // iteration counter
	PetscInt    it_Nwt;   // newton iteration counter
	PetscScalar refRes;   // reference residual norm
	PetscInt    nPicIt;   // number of Picard iteraions before switch to Newton
	PetscScalar rtolPic;  // relative Picard residual reduction tolerance
	PetscInt    nNwtIt;   // number of Newton iterations before switch to Picard
	PetscScalar rtolNwt;  // Newton divergence tolerance
} ;

//---------------------------------------------------------------------------

PetscErrorCode NLSolClear(NLSol *nl);

PetscErrorCode NLSolCreate(NLSol *nl, PCStokes pc, SNES *p_snes);

PetscErrorCode NLSolDestroy(NLSol *nl);

//---------------------------------------------------------------------------

// compute residual vector
PetscErrorCode FormResidual(SNES snes, Vec x, Vec f, void *ctx);

// compute Jacobian matrix and preconditioner
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx);

//---------------------------------------------------------------------------

PetscErrorCode JacApplyMFFD(Mat A, Vec x, Vec y);

//---------------------------------------------------------------------------

PetscErrorCode SNESPrintConvergedReason(SNES snes, 	PetscLogDouble t_beg);

PetscErrorCode SNESCoupledTest(
	SNES                snes,
	PetscInt            it,
	PetscReal           xnorm,
	PetscReal           gnorm,
	PetscReal           f,
	SNESConvergedReason *reason,
	void                *cctx);

//---------------------------------------------------------------------------

PetscErrorCode DisplaySpecifiedSolverOptions(PCStokes pc, SNES snes);

//---------------------------------------------------------------------------

#endif
