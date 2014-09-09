
#ifndef __Solvers_h__
#define __Solvers_h__

typedef struct {
	PetscInt norm_type; /* 10, 20, 30, L_inf, L2, L2/sqrt(L) */
	PetscBool use_relative_norm;
	PetscScalar *initial_nonzero_norm;
	PetscScalar fc_ksp_rtol[2],fc_ksp_atol[2],fc_ksp_min_atol[2];
	PetscInt	fc_ksp_min_it;
} KSPBlockStoppingConditionCtx;

//==========================================================================================================

// Solver 1
PetscErrorCode StokesSolve_PowellHestenes( Mat VV_MAT, Mat VP_MAT, Mat PV_MAT, Mat PP_MAT,
		Vec Velocity, Vec Pressure, Vec f, Vec h);

// Solver 2
PetscErrorCode StokesSolve_SCR(Mat A11, Mat A12, Mat A21, Mat A22, Mat pc_mat_schur,
		Vec x1, Vec x2, Vec b1, Vec b2, UserContext *user);

// Solver 3
PetscErrorCode StokesSolve_FC(Mat A11, Mat A12, Mat A21, Mat A22, Mat pc_mat_schur,
		Vec x1, Vec x2, Vec b1, Vec b2, UserContext *user);

// Velocity solver test
PetscErrorCode VelSolverTest(Mat A, Vec x, Vec b, UserContext *user);

//==========================================================================================================

PetscErrorCode Petsc_DMComposite_BlockMonitor(KSP ksp,PetscInt n,PetscScalar rnorm,void *mctx);

PetscErrorCode KSPStokes_DMComposite_BlockConvergenceTest(KSP ksp,PetscInt n,PetscScalar rnorm,KSPConvergedReason *reason,void *mctx);

PetscErrorCode Petsc_DMComposite_BlockStoppingConditionDestroy(void *data);

PetscErrorCode Petsc_DMComposite_BlockStoppingConditionConfiguration( KSP ksp );

PetscErrorCode CopyNestToBlock(DM da_nest, Vec vec_nest, Vec v_0, Vec v_1);

PetscErrorCode CopyBlockToNest(DM da_nest, Vec vec_nest, Vec v_0, Vec v_1);

//==========================================================================================================

#endif
