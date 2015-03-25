/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
 */

#include "LaMEM.h"
#include "Solvers.h"


//==========================================================================================================
// Perform Powell-Hestenes iterations, with a diagonal inverse PP matrix as preconditioner.
#undef __FUNCT__
#define __FUNCT__ "StokesSolve_PowellHestenes"
PetscErrorCode StokesSolve_PowellHestenes(Mat VV_MAT, Mat VP_MAT, Mat PV_MAT, Mat PP_MAT,
		Vec U, Vec P, Vec f, Vec h)
{
	PetscErrorCode ierr;
	Mat 			invPP_times_PV, VP_times_invPP_times_PV, cp_VV_MAT;
	Vec 			diagPP, f_new, dP, Div, invPP_h, PV_U;
	PetscScalar 	kappa, MaximumDivergence, div_max;
	KSP 			ksp_v;
	PetscInt 		PH_MaxInnerIterations, iter;
	PetscLogDouble	t0,t1;
	PetscBool 		assemble_inverse_mass_matrix, flg;

	PetscFunctionBegin;

	// initialize default parameters
	kappa                        = 1e4;
	PH_MaxInnerIterations        = 50;
	MaximumDivergence            = 1e-8;
	assemble_inverse_mass_matrix = PETSC_FALSE;

	// read options from command line
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-ph_kappa",                     &kappa,                        &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt   (PETSC_NULL, "-ph_max_it",                    &PH_MaxInnerIterations,        &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(PETSC_NULL, "-ph_tol",                       &MaximumDivergence,            &flg); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool  (PETSC_NULL, "-assemble_inverse_mass_matrix", &assemble_inverse_mass_matrix, &flg); CHKERRQ(ierr);

	// print actual options used
	PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes parameter: kappa  = %1.4e \n", kappa );
	PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes parameter: max_it = %lld \n", (LLD)PH_MaxInnerIterations );
	PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes parameter: tol    = %1.4e \n", MaximumDivergence );
	if (assemble_inverse_mass_matrix) {
		PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes parameter: using assembled inverse mass matrix \n");
	} else {
		PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes parameter: using Bt.inv[diag(M)].B \n");
	}

	PetscTime(&t0);

	if (!assemble_inverse_mass_matrix)
	{	ierr = MatGetVecs( PP_MAT, &diagPP, PETSC_NULL ); CHKERRQ(ierr);
		ierr = MatGetDiagonal(PP_MAT, diagPP); CHKERRQ(ierr);
		ierr = VecReciprocal(diagPP); CHKERRQ(ierr);      //  1/diagPP  as an approximation of inv(PP)
		ierr = VecScale( diagPP, kappa ); CHKERRQ(ierr);

		// A better option would be:
		// VecPointwiseMult( diagPP, diagPP, kappa );
		// where kappa is vector with the correct viscosity scaling

		ierr = MatDuplicate( PV_MAT, MAT_COPY_VALUES, &invPP_times_PV ); CHKERRQ(ierr);
		ierr = MatDiagonalScale( invPP_times_PV, diagPP, PETSC_NULL ); CHKERRQ(ierr);
		ierr = MatMatMult( VP_MAT, invPP_times_PV, MAT_INITIAL_MATRIX, 1.2, &VP_times_invPP_times_PV ); CHKERRQ(ierr);
	}
	else
	{
		ierr = MatScale(PP_MAT,kappa); CHKERRQ(ierr);
		ierr = MatMatMult( PP_MAT, PV_MAT, MAT_INITIAL_MATRIX, 1.2, &invPP_times_PV ); CHKERRQ(ierr);
		ierr = MatMatMult( VP_MAT, invPP_times_PV, MAT_INITIAL_MATRIX, 1.2, &VP_times_invPP_times_PV ); CHKERRQ(ierr);
	}

	// Add-hoc solution of the nonzero pattern inconsistency problem.
	// We just duplicate the VV_MAT input matrix. Doubles the memory, but not an issue
	// since this solver is not intended for large problems.

	MatDuplicate(VV_MAT, MAT_COPY_VALUES, &cp_VV_MAT);

	// This slows down stuff incredibly (think about constructing this matrix outside!)
	ierr = MatAXPY(cp_VV_MAT,1.0,VP_times_invPP_times_PV,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr); // VV_new = VV + kappa*VP*invPP*PV

	PetscTime(&t1);
	PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes setup: time = %e \n", t1-t0 );

	// (3) Create vectors
	ierr = VecDuplicate( f,  &f_new ); CHKERRQ(ierr);
	ierr = VecDuplicate( P,  &dP    ); CHKERRQ(ierr);
	ierr = VecDuplicate( P,  &Div   ); CHKERRQ(ierr);
	ierr = VecDuplicate( P,  &PV_U  ); CHKERRQ(ierr);
	ierr = VecDuplicate( P,  &invPP_h ); CHKERRQ(ierr);

	// (4)	Initialize solution and pressure vectors
	ierr = VecSet( U       ,0.0); CHKERRQ(ierr);
	ierr = VecSet( P       ,0.0); CHKERRQ(ierr);
	ierr = VecSet( f_new   ,0.0); CHKERRQ(ierr);
	ierr = VecSet( Div     ,0.0); CHKERRQ(ierr);

	// (5) Create KSP solver environment
	PetscTime(&t0);

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_v); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp_v, cp_VV_MAT, cp_VV_MAT); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp_v,"vs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_v); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp_v);CHKERRQ(ierr);

	PetscTime(&t1);
	PetscPrintf( PETSC_COMM_WORLD, "# Powell-Hestenes solver setup: time = %e \n", t1-t0 );

	// (6) Perform Powell-Hesteness iterations
	iter 	= 1;
	div_max = 1e10;
	while (div_max > MaximumDivergence){
		// f_new = f + VP*(diagPP*h - P)
		if (!assemble_inverse_mass_matrix) {
			ierr = VecPointwiseMult( invPP_h, diagPP, h); CHKERRQ(ierr);  // invPP_h = diagPP*h
		} else {
			ierr = MatMult(PP_MAT,h, invPP_h);CHKERRQ(ierr);
		}

		ierr = VecAXPY( invPP_h, -1.0, P); CHKERRQ(ierr);					// invPP_h = invPP_h - P
		ierr = MatMult( VP_MAT, invPP_h, f_new ); CHKERRQ(ierr);			// f_new = VP*invPP_h
		ierr = VecAXPY( f_new, 1.0, f); CHKERRQ(ierr);						// f_new = f_new + f

		// Solve Vel=VV_new\f_new
		ierr = KSPSolve(ksp_v,f_new,U); CHKERRQ(ierr);

		// dP = kappa*invPP*PV*Vel
		ierr = MatMult( PV_MAT, U, PV_U); CHKERRQ(ierr); 		// PV_U = PV*U
		ierr = VecAXPY( PV_U, -1.0, h); CHKERRQ(ierr); 			// PV_U = PV_U -h

		if (!assemble_inverse_mass_matrix) {
			ierr = VecPointwiseMult( dP, diagPP, PV_U); CHKERRQ(ierr); 	// dP = invPP*PV_U
		} else {
			ierr = MatMult(PP_MAT,PV_U, dP);CHKERRQ(ierr);
		}

		// P = P + dP
		ierr = VecAXPY(P, 1.0, dP); CHKERRQ(ierr);

		// Compute divergence
		ierr = MatMult(PV_MAT,U,Div); CHKERRQ(ierr);
		ierr = VecAXPY( Div, -1.0, h ); CHKERRQ(ierr); // Div <- Div - h

		ierr = VecNorm(Div,NORM_INFINITY,&div_max); CHKERRQ(ierr);

		{
			PetscScalar ru2,rp2;

			ierr = VecCopy(f,f_new);CHKERRQ(ierr);
			ierr = VecScale(f_new,-1.0);CHKERRQ(ierr);

			ierr = MatMultAdd(cp_VV_MAT,U,f_new,f_new);CHKERRQ(ierr); // [A + kappa.Bt.inv(M).B] U
			// remove contribution from kappa.Bt.inv(M).B
			ierr = VecScale(U,-1.0);CHKERRQ(ierr);
			ierr = MatMultAdd(VP_times_invPP_times_PV,U,f_new,f_new);CHKERRQ(ierr); // [- kappa.Bt.inv(M).B ]U
			ierr = VecScale(U,-1.0);CHKERRQ(ierr);

			ierr = MatMultAdd(VP_MAT,P,f_new,f_new);CHKERRQ(ierr);
			ierr = VecNorm(f_new,NORM_2,&ru2);CHKERRQ(ierr);
			ierr = VecNorm(Div,NORM_2,&rp2);CHKERRQ(ierr);

			// print progress
			PetscPrintf(PETSC_COMM_WORLD, "        Powell-Hesteness iteration %4D, |(ru,rp)| = %1.10e %1.10e : max(Div) = %1.10e \n",iter, ru2,rp2,div_max);

		}

		iter = iter + 1;

		// Emergency exit
		if (iter > PH_MaxInnerIterations){
			div_max = 0;
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Powell-Hesteness iterations did not converge; Check your setup and kappa.");
			break;
		}
	}

	// restore PP
	if (assemble_inverse_mass_matrix) {
		ierr = MatScale(PP_MAT,1.0/kappa); CHKERRQ(ierr);
	}

	// Cleaning up
	if (!assemble_inverse_mass_matrix) {
		ierr = VecDestroy(&diagPP); CHKERRQ(ierr);
	}
	ierr = VecDestroy(&PV_U); CHKERRQ(ierr);
	ierr = VecDestroy(&invPP_h); CHKERRQ(ierr);
	ierr = MatDestroy(&invPP_times_PV); CHKERRQ(ierr);
	ierr = MatDestroy(&VP_times_invPP_times_PV); CHKERRQ(ierr);
	ierr = VecDestroy(&f_new); CHKERRQ(ierr);
	ierr = VecDestroy(&dP); CHKERRQ(ierr);
	ierr = VecDestroy(&Div); CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp_v); CHKERRQ(ierr);

	ierr = MatDestroy(&cp_VV_MAT); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//==========================================================================================================
/*
 * Schur complement reduction
 *
 *  / A00  A01\   / x1 \     / b1 \
 * |           | |      | = |      | <=> S x2 = f_hat
 *  \ A10  A11/   \ x2 /     \ b2 /
 *
 * where:
 *
 * S = A11 - A10 inv(A00) A01 - Schur complement matrix
 *
 * f_hat = b2 - A10 inv(A00) b1 - reduced rhs
 *
 */
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "StokesSolve_SCR"
PetscErrorCode StokesSolve_SCR(Mat A11, Mat A12, Mat A21, Mat A22, Mat pc_mat_schur,
		Vec x1, Vec x2, Vec b1, Vec b2, UserContext *user)
{
	PetscErrorCode  ierr;
	Mat             S;
	KSP             ksp_scr, ksp_v;
	Vec             f_hat, tmp, h;
	PetscLogDouble  t0, t1;
	PetscBool       flg;

	// ksp_scr - Schur complement system solver (database prefix is "scr_")
	// ksp_v   - velocity block solver          (database prefix is "vs_")

	PetscFunctionBegin;

	if (b2==PETSC_NULL)
	{	ierr = VecDuplicate(x2, &h);CHKERRQ(ierr);
	ierr = VecZeroEntries(h);CHKERRQ(ierr);
	}
	else h = b2;

	// scale Schur preconditioning matrix  (PETSc FIELDSPLIT uses S = -A21.inv(A11).A12 + A22)
	ierr = MatScale(pc_mat_schur, -1.0); CHKERRQ(ierr);

	// hack needed in 3.4.2
	ierr = KSPMatRegisterAll(); CHKERRQ(ierr);

	// create Schur complement matrix
	ierr = MatCreateSchurComplement(A11, A11, A12, A21, A22, &S); CHKERRQ(ierr);

	// setup velocity solver
	ierr = MatSchurComplementGetKSP(S, &ksp_v); CHKERRQ(ierr);
	// Attach DM to the velocity solver to enable geometric multigrid
	if(user->VelocitySolver == 2)
	{	ierr = KSPSetDM(ksp_v, user->DA_Vel); CHKERRQ(ierr);
	ierr = KSPSetDMActive(ksp_v,PETSC_FALSE); CHKERRQ(ierr);
	}
	ierr = KSPSetOptionsPrefix(ksp_v,"vs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_v); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp_v);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp_v, PETSC_TRUE); CHKERRQ(ierr);


	// Tell us if there is a zero initial guess or not
	ierr = KSPGetInitialGuessNonzero(ksp_v, &flg);
	if (flg){
		//    PetscPrintf( PETSC_COMM_WORLD, "# SCR: using nonzero initial guess for ksp_v\n");
	}
	else{
		//   PetscPrintf( PETSC_COMM_WORLD, "# SCR: using zero initial guess for ksp_v\n");
	}

	// compute f_hat
	ierr = VecDuplicate(h, &f_hat);CHKERRQ(ierr);
	ierr = VecDuplicate(x1, &tmp);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_v, b1,tmp);CHKERRQ(ierr); // tmp = inv(A11).b1
	ierr = VecScale(tmp, -1.0);CHKERRQ(ierr);
	ierr = VecCopy(h, f_hat);CHKERRQ(ierr); // f_hat = h
	ierr = MatMultAdd(A21, tmp, f_hat, f_hat);CHKERRQ(ierr); // f_hat = f_hat + tmp = b2 - A21.inv(A11).b1

	// create and setup Schur complement solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_scr);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp_scr,"scr_");CHKERRQ(ierr);

	if (pc_mat_schur == PETSC_NULL) {
		ierr = KSPSetOperators(ksp_scr, S, S); CHKERRQ(ierr);
	} else {
		ierr = KSPSetOperators(ksp_scr, S, pc_mat_schur);CHKERRQ(ierr);
	}
	ierr = KSPSetFromOptions(ksp_scr); CHKERRQ(ierr);

	// solve for p
	ierr = KSPSetInitialGuessNonzero(ksp_scr, PETSC_TRUE); CHKERRQ(ierr);
	PetscTime(&t0);



	ierr = KSPSolve(ksp_scr, f_hat, x2);CHKERRQ(ierr);
	PetscTime(&t1);
	PetscPrintf( PETSC_COMM_WORLD, "# Time: p solve = %2.2e \n", t1-t0 );

	// obtain solution for u
	ierr = MatMult(A12, x2, tmp);CHKERRQ(ierr);
	ierr = VecAYPX(tmp, -1.0, b1);CHKERRQ(ierr); // t <- -t + b1
	PetscTime(&t0);
	ierr = KSPSolve(ksp_v, tmp, x1);CHKERRQ(ierr);
	PetscTime(&t1);
	PetscPrintf( PETSC_COMM_WORLD, "# Time: u solve = %2.2e \n", t1-t0 );

	// clean up
	ierr = KSPDestroy(&ksp_scr);CHKERRQ(ierr);
	ierr = VecDestroy(&f_hat);CHKERRQ(ierr);
	ierr = VecDestroy(&tmp);CHKERRQ(ierr);
	ierr = MatDestroy(&S);CHKERRQ(ierr);
	if (b2==PETSC_NULL) {
		ierr = VecDestroy(&h);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "CopyNestToBlock"
PetscErrorCode CopyNestToBlock(DM da_nest, Vec vec_nest, Vec v_0, Vec v_1)
{
	// function copies contents of nested vector into separate blocks
	Vec v_tmp[2];
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = DMCompositeGetAccess(da_nest, vec_nest, &v_tmp[0], &v_tmp[1]);     CHKERRQ(ierr);
	ierr = VecCopy(v_tmp[0], v_0);                                            CHKERRQ(ierr);
	ierr = VecCopy(v_tmp[1], v_1);                                            CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(da_nest, vec_nest, &v_tmp[0], &v_tmp[1]); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "CopyBlockToNest"
PetscErrorCode CopyBlockToNest(DM da_nest, Vec vec_nest, Vec v_0, Vec v_1)
{
	// function assembles nested vector from separate blocks
	Vec v_tmp[2];
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = DMCompositeGetAccess(da_nest, vec_nest, &v_tmp[0], &v_tmp[1]);     CHKERRQ(ierr);
	ierr = VecCopy(v_0, v_tmp[0]);                                            CHKERRQ(ierr);
	ierr = VecCopy(v_1, v_tmp[1]);                                            CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(da_nest, vec_nest, &v_tmp[0], &v_tmp[1]); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//==========================================================================================================
// MatNest + Fieldsplit
#undef __FUNCT__
#define __FUNCT__ "StokesSolve_FC"
PetscErrorCode StokesSolve_FC(Mat A11, Mat A12, Mat A21, Mat A22, Mat pc_mat_schur,
		Vec x1, Vec x2, Vec b1, Vec b2, UserContext *user)
{
	PetscErrorCode ierr;
	Mat stokes[4], A, PP_MAT, PP_MAT_ZERO = PETSC_NULL;
	PetscInt nv_s, nv_e, np_s, np_e, nv, np, nsplits, Verbose=0;
	IS  isu, isp;
	Vec x, b, P_lithostatic, rhs_correction;
	KSP ksp, *sub_ksp;
	PC  pc;
	PetscBool use_stokes_res, use_stokes_monitor, use_def_schur_prec, flg, EliminateLithostaticPressure;
	PCSide side;
	DM	DM_Stokes;
	PetscLogDouble     cputime_start, cputime_end;

	PetscFunctionBegin;

	PetscOptionsGetInt(PETSC_NULL ,"-Verbose",	&Verbose	, PETSC_NULL); //default 0; (1: show timing info)

	//================
	// PRESSURE MATRIX
	//================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating PP matrix ... ");
	}
	if(A22 == PETSC_NULL)
	{	ierr = MatDuplicate(user->PP_MAT, MAT_DO_NOT_COPY_VALUES, &PP_MAT_ZERO); CHKERRQ(ierr);
	ierr = MatZeroEntries(PP_MAT_ZERO);                                     CHKERRQ(ierr);
	PP_MAT = PP_MAT_ZERO;
	}
	else PP_MAT = A22;
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}


	//=============
	// COMPOSITE DM
	//=============
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating Composite DM ... ");
	}

	// create DMComposite that contains both the Velocity and Pressure DA's
	ierr =	DMCompositeCreate(PETSC_COMM_WORLD, &DM_Stokes); CHKERRQ(ierr);
	ierr =	DMCompositeAddDM(DM_Stokes, user->DA_Vel);		 CHKERRQ(ierr);	// Velocity
	ierr =	DMCompositeAddDM(DM_Stokes, user->DA_Pres);		 CHKERRQ(ierr);	// Pressure

	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}
	//=======
	// MATRIX
	//=======




	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating Nested A matrix ... ");
	}
	// create nested matrix
	stokes[0] = A11;
	stokes[1] = A12;
	stokes[2] = A21;
	stokes[3] = PP_MAT;

	ierr = MatCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, 2, PETSC_NULL, (const Mat*)stokes, &A); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}

	//====================
	// ELIMINATE (LITHOSTATIC) PRESSURE FROM INITIAL SYSTEM OF EQUATIONS
	//====================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Eliminating lithostatic pressure if required ... ");
	}
	EliminateLithostaticPressure = PETSC_FALSE;
	PetscOptionsGetBool( PETSC_NULL,"-EliminateLithostaticPressure",&EliminateLithostaticPressure,PETSC_FALSE );
	if (EliminateLithostaticPressure){

		// initialization stuff
		ierr = VecDuplicate(x2,&P_lithostatic);         CHKERRQ(ierr);  // create vector
		ierr = VecCopy(x2,P_lithostatic);               CHKERRQ(ierr);  // copy values
		ierr = VecDuplicate(b1,&rhs_correction);        CHKERRQ(ierr);
		ierr = VecZeroEntries(rhs_correction);          CHKERRQ(ierr);  // ensure it's initialized to zero

		// RHS correction = A11*P_lithostatic
		ierr = MatMult(A12, P_lithostatic, rhs_correction);

		// Add to rhs of velocity   f = f-rhs_correction
		ierr = VecAXPY(b1, -1.0 ,rhs_correction);       CHKERRQ(ierr);

		// The initial pressure vector must be zeroed out
		ierr = VecZeroEntries(x2);                      CHKERRQ(ierr);

	}
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}





	//====================
	// RHS & INITIAL GUESS
	//====================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating RHS and initial guess vectors ... ");
	}
	// generate temporary nested solution vectors from DMComposite object
	ierr = DMCreateGlobalVector(DM_Stokes, &b); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(DM_Stokes, &x); CHKERRQ(ierr);

	// copy rhs and initial guess to nested vectors
	ierr = CopyBlockToNest(DM_Stokes, b, b1, b2);
	ierr = CopyBlockToNest(DM_Stokes, x, x1, x2);
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}
	//=========================
	// KSP FULLY COUPLED SOLVER
	//=========================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating KSP Solver ... ");
	}
	// create krylov solver for fully coupled system, attach DMComposite & prefix
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);                CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, A, A);                       CHKERRQ(ierr);
	ierr = KSPSetDM(ksp, DM_Stokes);                         CHKERRQ(ierr);
	ierr = KSPSetDMActive(ksp, PETSC_FALSE);                 CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp,"fc_");                   CHKERRQ(ierr);
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}

	//==========================
	// FIELDSPLIT PRECONDITIONER
	//==========================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Creating Fieldsplit and SCHUR Preconditioner ... ");
	}

	// access preconditioner object
	ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

	// define block upper triangular preconditioner
	ierr = PCSetType(pc, PCFIELDSPLIT);                                      CHKERRQ(ierr);
	ierr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);                      CHKERRQ(ierr);
	ierr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER); CHKERRQ(ierr);

	// create velocity and pressure index sets
	ierr = VecGetOwnershipRange(x1, &nv_s, &nv_e); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(x2, &np_s, &np_e); CHKERRQ(ierr);
	nv = nv_e - nv_s;
	np = np_e - np_s;
	ierr = ISCreateStride(PETSC_COMM_WORLD, nv, nv_s+np_s,   1, &isu); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, np, nv_s+np_s+nv,1, &isp); CHKERRQ(ierr);

	// pass index sets to define fieldsplit
	ierr = PCFieldSplitSetIS(pc, PETSC_NULL, isu);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc, PETSC_NULL, isp);CHKERRQ(ierr);

	// set custom schur complement preconditioner (if requested)
	use_def_schur_prec = PETSC_FALSE;
	PetscOptionsGetBool(PETSC_NULL, "-use_default_schur_preconditioner", &use_def_schur_prec, &flg);
	if(use_def_schur_prec == PETSC_FALSE)
	{	// scale schur preconditioning matrix  (PETSc FIELDSPLIT uses S = -A21.inv(A11).A12 + A22)
		ierr = MatScale(pc_mat_schur, -1.0); CHKERRQ(ierr);
		ierr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, pc_mat_schur);
	}

	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}

	//==================
	// FINALIZE KSP & PC
	//==================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Finalize KSP and PC setting up ... ");
	}
	// read additional command line options
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

	// force setup to enable fetching sub-ksp
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	// fetch sub-ksp objects for velocity & pressure splits
	ierr = PCFieldSplitGetSubKSP(pc, &nsplits, &sub_ksp); CHKERRQ(ierr);
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}

	//================
	// VELOCITY SOLVER
	//================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Setting up velocity solver ... ");
	}
	// Attach DM to the velocity solver to enable geometric multigrid
	if(user->VelocitySolver == 2)
	{	ierr = KSPSetDM(sub_ksp[0], user->DA_Vel);      CHKERRQ(ierr);
	ierr = KSPSetDMActive(sub_ksp[0], PETSC_FALSE); CHKERRQ(ierr);
	}

	ierr = KSPSetOptionsPrefix(sub_ksp[0],"vs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(sub_ksp[0]);         CHKERRQ(ierr);
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, " almost ... ");
	}
	ierr = KSPSetUp(sub_ksp[0]);                  CHKERRQ(ierr);
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}



	/* Set Nullspace for VV matrix

	{


			MatNullSpace           matnull;
			Vec             		vel_coords;


			ierr = DMDAGetCoordinates(user->DA_Vel,&vel_coords);   		CHKERRQ(ierr);
			ierr = MatNullSpaceCreateRigidBody(vel_coords,&matnull);   	CHKERRQ(ierr);
//			ierr = MatSetNearNullSpace(A11,matnull);  		 	CHKERRQ(ierr);

			ierr = KSPSetNullSpace(sub_ksp[0], matnull); CHKERRQ(ierr);

			//MatNullSpaceView(matnull,PetscViewer viewer)

			ierr = MatNullSpaceDestroy(&matnull);  						CHKERRQ(ierr);

		}
	 */


	//===================================
	// SCHUR COMPLEMENT (PRESSURE) SOLVER
	//===================================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Setting up pressure solver ... ");
	}
	ierr = KSPSetOptionsPrefix(sub_ksp[1],"ps_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(sub_ksp[1]);         CHKERRQ(ierr);
	ierr = KSPSetUp(sub_ksp[1]);                  CHKERRQ(ierr);
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done [%g] s \n ", cputime_end-cputime_start);
	}
	//============================
	// STOPPING CONDITON & MONITOR
	//============================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Defining stopping condition and monitor ... \n");
	}
	// define a block based stopping condition
	use_stokes_res = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL, "-use_stokes_residual", &use_stokes_res, &flg); CHKERRQ(ierr);
	if(use_stokes_res == PETSC_TRUE)
	{
		ierr = KSPGetPCSide(ksp, &side); CHKERRQ(ierr);
		if(side == PC_LEFT)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Using preconditioned Stokes residual to define the FC stopping condition \n");
		}
		if(side == PC_RIGHT)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Using unpreconditioned (true) Stokes residual to define the FC stopping condition \n");
		}
		ierr = Petsc_DMComposite_BlockStoppingConditionConfiguration(ksp); CHKERRQ(ierr);
	}

	// set residual monitors
	use_stokes_monitor = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL, "-use_stokes_monitor", &use_stokes_monitor, &flg ); CHKERRQ(ierr);
	if(use_stokes_monitor == PETSC_TRUE )
	{
		PetscPrintf( PETSC_COMM_WORLD, "StokesResidual: Monitor activated \n");
		ierr = KSPMonitorSet(ksp, Petsc_DMComposite_BlockMonitor, 0, 0); CHKERRQ(ierr);
	}

	// Tell us if there is a zero initial guess or not
	ierr = KSPGetInitialGuessNonzero(ksp, &flg);
	if (flg){
		PetscPrintf( PETSC_COMM_WORLD, "# FC: using nonzero initial guess \n");
	}
	else{
		PetscPrintf( PETSC_COMM_WORLD, "# FC: using zero initial guess \n");
	}
	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done with stopping condition and monitor [%g] s \n ", cputime_end-cputime_start);
	}

	//=====================
	// SOLVE COUPLED SYSTEM
	//=====================
	if (Verbose>0){
		PetscTime(&cputime_start);
		PetscPrintf(PETSC_COMM_WORLD, "Solving system of equations ... \n");
	}

	ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

	// copy solution blocks from nested vector
	ierr = CopyNestToBlock(DM_Stokes, x, x1, x2); CHKERRQ(ierr);

	if (Verbose>0){
		PetscTime(&cputime_end);
		PetscPrintf(PETSC_COMM_WORLD, " done solving system of equations [%g] s \n ", cputime_end-cputime_start);
	}


	//===================
	// ADD LITHOSTATIC COMPONENT BACK TO PRESSURE
	//===================
	if (EliminateLithostaticPressure ){
		PetscScalar norm_p;

		VecNorm(x2, NORM_2 ,&norm_p);
		PetscPrintf(PETSC_COMM_WORLD,"Test pressure norm before adding lithostatic component: norm2 = %f \n", norm_p);

		ierr = VecAXPY(x2, 1.0 ,P_lithostatic);     CHKERRQ(ierr);


		VecNorm(x2, NORM_2 ,&norm_p);
		PetscPrintf(PETSC_COMM_WORLD,"Test pressure norm after adding lithostatic component: norm2 = %f \n", norm_p);


		VecNorm(P_lithostatic, NORM_2 ,&norm_p);
		PetscPrintf(PETSC_COMM_WORLD," lithostatic component: norm2 = %f \n", norm_p);


	}




	//========
	// CLEANUP
	//========

	// delete temporary objects
	ierr = KSPDestroy(&ksp);      CHKERRQ(ierr);
	ierr = ISDestroy(&isu);       CHKERRQ(ierr);
	ierr = ISDestroy(&isp);       CHKERRQ(ierr);
	ierr = VecDestroy(&x);        CHKERRQ(ierr);
	ierr = VecDestroy(&b);        CHKERRQ(ierr);
	ierr = MatDestroy(&A);        CHKERRQ(ierr);
	ierr = DMDestroy(&DM_Stokes); CHKERRQ(ierr);
	ierr = PetscFree(sub_ksp);    CHKERRQ(ierr);

	if (EliminateLithostaticPressure){
		ierr = VecDestroy(&P_lithostatic);   CHKERRQ(ierr);
		ierr = VecDestroy(&rhs_correction);  CHKERRQ(ierr);
	}

	if(PP_MAT_ZERO != PETSC_NULL) { ierr = MatDestroy(&PP_MAT_ZERO); CHKERRQ(ierr); }

	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "Petsc_DMComposite_BlockMonitor"
PetscErrorCode Petsc_DMComposite_BlockMonitor(KSP ksp,PetscInt n,PetscScalar rnorm, void *mctx)
{
	PetscErrorCode		ierr;
	PetscInt			b, nDM;
	Vec					Br,r[5],v,w;
	DM					DM_Stokes;
	DMType              DM_Stokes_type;
	Mat					A;
	PetscBool			is_dmcomposite=PETSC_FALSE;
	PetscFunctionBegin;

	if(rnorm) rnorm = 0.0;
	if(mctx)  mctx  = NULL;

	// Compute residual vector for ALL blocks [velocity & pressure in case of Stokes ]
	ierr =	KSPGetOperators( ksp, &A, NULL);		CHKERRQ(ierr);
	ierr =	MatGetVecs( A, &w, &v );				CHKERRQ(ierr);
	ierr =	KSPBuildResidual( ksp, v, w, &Br );		CHKERRQ(ierr);

	// get DM attached to ksp
	ierr =  KSPGetDM(ksp, &DM_Stokes); CHKERRQ(ierr);
	// check that DMType is a DMCOMPOSITE
	ierr = DMGetType(DM_Stokes, &DM_Stokes_type); CHKERRQ(ierr);
	is_dmcomposite = (PetscBool)(!strcmp(DM_Stokes_type, DMCOMPOSITE));


	if (is_dmcomposite == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"The routine PetscMatNestFS_BlockMonitor only works for KSP that have a DMCOMPOSITE attached to it.");
	}

	/* Fetch vectors for every sub-block in DMCOMPOSITE */
	ierr =	DMCompositeGetNumberDM(DM_Stokes,&nDM);						CHKERRQ(ierr);
	if		(nDM==1){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0]);								CHKERRQ(ierr);
	}
	else if (nDM==2){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0], &r[1]);						CHKERRQ(ierr);
	}
	else if (nDM==3){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0], &r[1], &r[2]);					CHKERRQ(ierr);
	}
	else if (nDM==4){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0], &r[1], &r[2], &r[3]);			CHKERRQ(ierr);
	}
	else if (nDM==5){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0], &r[1], &r[2], &r[3], &r[4]);	CHKERRQ(ierr);
	}
	else {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"The routine PetscMatNestFS_BlockMonitor can currently only handle %i blocks.",nDM);
	}

	/* Output to screen */
	if (n == 0 && ((PetscObject)ksp)->prefix) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norms for %s solve.\n",((PetscObject)ksp)->prefix);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  # KSP Iterations      /  L2_norm             scaled_L2_norm      Linf_norm           /  ... for each %3D blocks  [%s] \n", nDM, ((PetscObject)ksp)->prefix ); CHKERRQ(ierr);
	}

	/* In case we are solving a Stokes system, we can dump the solution vectors to disk at every intermediate step [for debugging only] */



	ierr = PetscPrintf(PETSC_COMM_WORLD,"%3D KSP Residual norms  ",n);CHKERRQ(ierr);
	for( b=0; b<nDM; b++ ) {
		PetscScalar L2_norm, scaled_L2_norm, Linf_norm;
		PetscInt length;

		ierr = VecGetSize( r[b], &length ); CHKERRQ(ierr);
		ierr = VecNorm( r[b], NORM_2, &L2_norm ); CHKERRQ(ierr);
		ierr = VecNorm( r[b], NORM_INFINITY, &Linf_norm ); CHKERRQ(ierr);
		scaled_L2_norm = ( L2_norm ) / sqrt( (PetscScalar)length );


		ierr = PetscPrintf(PETSC_COMM_WORLD, "/  %14.12e  %14.12e  %14.12e  ",L2_norm, scaled_L2_norm, Linf_norm );CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"/  [%s] \n", ((PetscObject)ksp)->prefix );CHKERRQ(ierr);

	/* Cleaning up */
	if		(nDM==1){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0]);								CHKERRQ(ierr);
	}
	else if (nDM==2){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0], &r[1]);						CHKERRQ(ierr);
	}
	else if (nDM==3){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0], &r[1], &r[2]);					CHKERRQ(ierr);
	}
	else if (nDM==4){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0], &r[1], &r[2], &r[3]);			CHKERRQ(ierr);
	}
	else if (nDM==5){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0], &r[1], &r[2], &r[3], &r[4]);	CHKERRQ(ierr);
	}

	ierr = VecDestroy(&v);CHKERRQ(ierr);
	ierr = VecDestroy(&w);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//==========================================================================================================
// Stops when |r_i|_2 / length(r_i) < abs_tolerance
#undef __FUNCT__
#define __FUNCT__ "KSPStokes_DMComposite_BlockConvergenceTest"
PetscErrorCode KSPStokes_DMComposite_BlockConvergenceTest(KSP ksp,PetscInt n,PetscScalar rnorm,KSPConvergedReason *reason,void *mctx)
{
	PetscErrorCode		ierr;
	PetscScalar			norms[10];
	PetscInt			sizes[10];
	PetscInt			b, nDM;
	Vec					r[2], Br;
	Vec					w, v;
	PetscBool			rel_crit, abs_crit,rel_add_crit,abs_add_crit,stp_crit[2];
	PetscInt			stp_reason[2][2];
	Mat					A;
	DM					DM_Stokes;
	DMType              DM_Stokes_type;
	PetscBool			is_dmcomposite=PETSC_FALSE;

	if(rnorm) rnorm = 0.0;
	KSPBlockStoppingConditionCtx *ctx = (KSPBlockStoppingConditionCtx*)mctx;

	// Compute residual vector for ALL blocks [velocity & pressure in case of Stokes ]
	ierr = KSPGetOperators(ksp, &A, NULL);    CHKERRQ(ierr);
	ierr = MatGetVecs(A, &w, &v);              CHKERRQ(ierr);
	ierr = KSPBuildResidual(ksp, v, w, &Br);  CHKERRQ(ierr);

	// get DM attached to ksp
	ierr =  KSPGetDM(ksp, &DM_Stokes); CHKERRQ(ierr);

	// check that DMType is a DMCOMPOSITE
	ierr = DMGetType(DM_Stokes, &DM_Stokes_type); CHKERRQ(ierr);
	is_dmcomposite = (PetscBool)(!strcmp(DM_Stokes_type, DMCOMPOSITE));

	if (is_dmcomposite == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"The routine KSPStokes_DMComposite_BlockConvergenceTest only works for KSP that have a DMCOMPOSITE attached to it.");
	}

	/* Fetch vectors for every sub-block in DMCOMPOSITE */
	ierr =	DMCompositeGetNumberDM(DM_Stokes,&nDM);										CHKERRQ(ierr);
	if		(nDM==1){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0]);								CHKERRQ(ierr);
	}
	else if (nDM==2){
		ierr =	DMCompositeGetAccess(DM_Stokes,Br,&r[0], &r[1]);						CHKERRQ(ierr);
	}
	else {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"The routine KSPStokes_DMComposite_BlockConvergenceTest can currently only handle %i blocks.",nDM);
	}

	if(ctx->initial_nonzero_norm==PETSC_NULL) {
		ierr = PetscMalloc( sizeof(PetscScalar)*(size_t)nDM, &ctx->initial_nonzero_norm );CHKERRQ(ierr);
	}

	/* Set the residual */
	for( b=0; b<nDM; b++ ) {

		ierr = VecGetSize( r[b], &sizes[b] ); CHKERRQ(ierr);

		/* Linf */
		if( ctx->norm_type == 10 ) {
			ierr		=	VecNorm( r[b], NORM_INFINITY, &norms[b] ); CHKERRQ(ierr);
		}
		/* L2 */
		else if( ctx->norm_type == 20 ) {
			ierr		=	VecNorm( r[b], NORM_2, &norms[b] ); CHKERRQ(ierr);
		}
		/* scaled L2 */
		else if( ctx->norm_type == 30 ) {
			ierr		=	VecNorm( r[b], NORM_2, &norms[b] ); CHKERRQ(ierr);
			norms[b]	=	norms[b]/sqrt(sizes[b]);
		}
		else {
			SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown norm option provided");
		}


		if (n==1) {
			ctx->initial_nonzero_norm[b] = norms[b];
		}
	}

	/* Cleanup */
	if		(nDM==1){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0]);								CHKERRQ(ierr);
	}
	else if (nDM==2){
		ierr =	DMCompositeRestoreAccess(DM_Stokes,Br,&r[0], &r[1]);						CHKERRQ(ierr);
	}

	ierr = VecDestroy( &v );		CHKERRQ(ierr);
	ierr = VecDestroy( &w );		CHKERRQ(ierr);

	// --- Check convergence of KSP ---
	stp_crit[0]		=	PETSC_FALSE;
	stp_crit[1]		=	PETSC_FALSE;
	rel_crit 		= 	PETSC_FALSE;
	abs_crit 		= 	PETSC_FALSE;
	rel_add_crit	= 	PETSC_FALSE;
	abs_add_crit	= 	PETSC_FALSE;

	for( b=0; b<nDM; b++ ) {

		// check relative tolerance
		if(n>=1){
			if( norms[b]/ctx->initial_nonzero_norm[b] < ctx->fc_ksp_rtol[b] ) {
				if((norms[b] < ctx->fc_ksp_min_atol[b]) && (n >= ctx->fc_ksp_min_it)){
					rel_crit 		 = PETSC_TRUE;
					rel_add_crit  	 = PETSC_TRUE;
					stp_reason[b][0] = 1;
				}
				else{
					rel_crit 		 = PETSC_TRUE;
					rel_add_crit	 = PETSC_FALSE;
					stp_reason[b][0] = 0;
				}
			}else{
				rel_crit 			= PETSC_FALSE;
				rel_add_crit  		= PETSC_FALSE;
				stp_reason[b][0] 	= 0;
			}
		}

		// check absolute tolerance
		if((norms[b] < ctx->fc_ksp_atol[b]) && (n >= ctx->fc_ksp_min_it)) {
			abs_crit 		 = PETSC_TRUE;
			abs_add_crit  	 = PETSC_TRUE;
			stp_reason[b][1] = 1;
		}
		else{
			abs_crit 		 = PETSC_FALSE;
			abs_add_crit  	 = PETSC_FALSE;
			stp_reason[b][1] = 0;
		}

		// combined decision
		if((abs_crit && abs_add_crit) || (rel_crit && rel_add_crit)){
			stp_crit[b]		 =	PETSC_TRUE;
		}
		else{
			stp_crit[b]		 =	PETSC_FALSE;
		}
	}

	if(stp_crit[0] && stp_crit[1] ){
		// stop iterations
		*reason = KSP_CONVERGED_HAPPY_BREAKDOWN;

		// print reason
		if(stp_reason[0][0]==1 && stp_reason[0][1]==0) PetscPrintf( PETSC_COMM_WORLD, "Velocity convergence reason: RTOL \n");
		if(stp_reason[0][0]==0 && stp_reason[0][1]==1) PetscPrintf( PETSC_COMM_WORLD, "Velocity convergence reason: ATOL \n");
		if(stp_reason[0][0]==1 && stp_reason[0][1]==1) PetscPrintf( PETSC_COMM_WORLD, "Velocity convergence reason: ATOL + RTOL\n");

		if(stp_reason[1][0]==1 && stp_reason[1][1]==0) PetscPrintf( PETSC_COMM_WORLD, "Pressure convergence reason: RTOL \n");
		if(stp_reason[1][0]==0 && stp_reason[1][1]==1) PetscPrintf( PETSC_COMM_WORLD, "Pressure convergence reason: ATOL \n");
		if(stp_reason[1][0]==1 && stp_reason[1][1]==1) PetscPrintf( PETSC_COMM_WORLD, "Pressure convergence reason: ATOL + RTOL\n");

	}
	else{
		// continue iterations
		*reason = KSP_CONVERGED_ITERATING;
		PetscFunctionReturn(0);
	}


	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "Petsc_DMComposite_BlockStoppingConditionDestroy"
PetscErrorCode Petsc_DMComposite_BlockStoppingConditionDestroy(void *data)
{
	PetscErrorCode ierr;
	KSPBlockStoppingConditionCtx *ctx = (KSPBlockStoppingConditionCtx*)data;
	if(ctx->initial_nonzero_norm) { ierr = PetscFree(ctx->initial_nonzero_norm);CHKERRQ(ierr); }
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "Petsc_DMComposite_BlockStoppingConditionConfiguration"
PetscErrorCode Petsc_DMComposite_BlockStoppingConditionConfiguration( KSP ksp )
{
	PetscErrorCode ierr;
	PetscBool      flg, is_block = PETSC_FALSE;
	PetscInt       norm_int;
	Mat            A;
	MatType        A_type;
	KSPBlockStoppingConditionCtx *ctx;
	char opt_string[PETSC_MAX_PATH_LEN];

	// check operators are blocks
	ierr     = KSPGetOperators( ksp, &A, NULL); CHKERRQ(ierr);
	ierr     = MatGetType(A, &A_type); CHKERRQ(ierr);
	is_block = (PetscBool)(!strcmp(A_type, MATNEST));

	if( is_block == PETSC_FALSE ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_SUP, "Only valid for Mat type \"nest\"" );
	}

	/* look for norm option choice */
	/* default */
	norm_int = 20;
	ierr = PetscOptionsGetString(PETSC_NULL, "-use_stokes_norm", opt_string, PETSC_MAX_PATH_LEN, &flg); CHKERRQ(ierr);
	if      (!strcmp(opt_string, "Linf"))     norm_int = 10;
	else if (!strcmp(opt_string, "L2"))       norm_int = 20;
	else if (!strcmp(opt_string, "scaledL2")) norm_int = 30;
	else    PetscPrintf( PETSC_COMM_WORLD, "Incorrect value specified for -use_stokes_norm option (%s). Using L2 instead\n", opt_string);

	/* print norm type */
	if     ( norm_int==10 ) { PetscPrintf( PETSC_COMM_WORLD, "StokesResidual: Using L_inf \n"); }
	else if( norm_int==20)  { PetscPrintf( PETSC_COMM_WORLD, "StokesResidual: Using L_2 \n"); }
	else if( norm_int==30)  { PetscPrintf( PETSC_COMM_WORLD, "StokesResidual: Using scaled L_2 \n"); }

	ierr = PetscMalloc( sizeof(KSPBlockStoppingConditionCtx), &ctx );CHKERRQ(ierr);
	ctx->norm_type = norm_int;
	ctx->initial_nonzero_norm = PETSC_NULL;

	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_v_atol",&ctx->fc_ksp_atol[0],PETSC_NULL );			CHKERRQ(ierr);
	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_v_rtol",&ctx->fc_ksp_rtol[0],PETSC_NULL );			CHKERRQ(ierr);
	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_v_min_atol",&ctx->fc_ksp_min_atol[0],PETSC_NULL );	CHKERRQ(ierr);
	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_p_atol",&ctx->fc_ksp_atol[1],PETSC_NULL );			CHKERRQ(ierr);
	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_p_rtol",&ctx->fc_ksp_rtol[1],PETSC_NULL );			CHKERRQ(ierr);
	ierr = PetscOptionsGetReal( PETSC_NULL,"-fc_ksp_p_min_atol",&ctx->fc_ksp_min_atol[1],PETSC_NULL );	CHKERRQ(ierr);
	ierr = PetscOptionsGetInt( PETSC_NULL,"-fc_ksp_min_it",&ctx->fc_ksp_min_it,PETSC_NULL );			CHKERRQ(ierr);
	ierr = PetscPrintf( PETSC_COMM_WORLD, "Stopping conditions: \n");
	//ierr = PetscPrintf( PETSC_COMM_WORLD, "Velocity: RTOL: %e	ATOL: %e\n",ctx->fc_ksp_rtol[0],ctx->fc_ksp_atol[0]);	CHKERRQ(ierr);
	//ierr = PetscPrintf( PETSC_COMM_WORLD, "Pressure: RTOL: %e	ATOL: %e\n",ctx->fc_ksp_rtol[1],ctx->fc_ksp_atol[1]);	CHKERRQ(ierr);
	//ierr = PetscPrintf( PETSC_COMM_WORLD, "# Additional criteria: minATOL(v): %e	minATOL(p): %efc_ksp_min_it %lld\n",ctx->fc_ksp_minatol[0],ctx->fc_ksp_minatol[1],(LLD)ctx->fc_ksp_min_it);	CHKERRQ(ierr);

	ierr = PetscPrintf( PETSC_COMM_WORLD, "Velocity: RTOL: %e	ATOL: %e	minATOL: %e		fc_ksp_min_it %lld\n",ctx->fc_ksp_rtol[0],ctx->fc_ksp_atol[0],ctx->fc_ksp_min_atol[0],(LLD)ctx->fc_ksp_min_it);	CHKERRQ(ierr);
	ierr = PetscPrintf( PETSC_COMM_WORLD, "Pressure: RTOL: %e	ATOL: %e	minATOL: %e		fc_ksp_min_it %lld\n",ctx->fc_ksp_rtol[1],ctx->fc_ksp_atol[1],ctx->fc_ksp_min_atol[1],(LLD)ctx->fc_ksp_min_it);	CHKERRQ(ierr);

	ierr = KSPSetConvergenceTest( ksp, KSPStokes_DMComposite_BlockConvergenceTest, (void*)ctx, Petsc_DMComposite_BlockStoppingConditionDestroy ); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "VelSolverTest"
PetscErrorCode VelSolverTest(Mat A, Vec x, Vec b, UserContext *user)
{
	// test performance of single velocity solve with zero initial guess
	KSP            ksp;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	// create Krylov solver for the velocity system, attach DM & prefix
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);  CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, A, A);        CHKERRQ(ierr);
	if(user->VelocitySolver == 2)
	{	ierr = KSPSetDM(ksp, user->DA_Vel);       CHKERRQ(ierr);
		ierr = KSPSetDMActive(ksp, PETSC_FALSE); CHKERRQ(ierr);
	}
	ierr = KSPSetOptionsPrefix(ksp,"vs_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);         CHKERRQ(ierr);
	ierr = KSPSetUp(ksp);                   CHKERRQ(ierr);

	// solve velocity subsystem
	ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================


