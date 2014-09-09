
#ifndef __Assembly_FDSTAG_h__
#define __Assembly_FDSTAG_h__

PetscErrorCode DAGetGlobalIndex(DM da,PetscInt i,PetscInt j,PetscInt k,PetscInt d,PetscInt *gidx);
PetscErrorCode FDSTAG_GetStencilValues(DM dav,DM dap,PetscInt vsL,MatStencil vs[],PetscInt psL,MatStencil ps[],PetscInt vgidx[],PetscInt pgidx[]);

PetscErrorCode ComputeStiffnessMatrixes_FDSTAG( DM da, DM da_pres,
		Mat VV_MAT, 		Mat PP_MAT, 	Mat PV_MAT, 	Mat VP_MAT,
		UserContext *user, 	PetscScalar dt, Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec pv_rhs_push, Vec ViscosityScaling, Mat approx_S );

PetscErrorCode FDSTAG_FDStencil_Incompressibility(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[], PetscScalar dy_vec[], PetscScalar dz_vec[],
		MatStencil  row_pp[], PetscScalar v_pp[],
		MatStencil  col_pv[], PetscScalar v_pv[]);

PetscErrorCode FDSTAG_FDStencil_ForceBalance1(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XY[4],		PetscScalar Eta_XZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar		v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[]);

PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance1(PetscInt j, PetscInt k, PetscInt nel_y, PetscInt nel_z,
		PetscScalar v_vv[], UserContext *user);

PetscErrorCode FDSTAG_FDStencil_ForceBalance2(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XY[4],		PetscScalar Eta_YZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar	v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[]);

PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance2(PetscInt i, PetscInt k, PetscInt nel_x, PetscInt nel_z,
		PetscScalar v_vv[], UserContext *user);

PetscErrorCode FDSTAG_FDStencil_ForceBalance3(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XZ[4],		PetscScalar Eta_YZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar	v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[], 		PetscScalar	FreeSurface_Fraction);

PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance3(PetscInt i, PetscInt j, PetscInt nel_x, PetscInt nel_y,
		PetscScalar v_vv[], PetscScalar v_vp[], PetscInt FreeSurfaceCells[], UserContext *user, PetscBool eliminate_stickyair_from_system);

PetscErrorCode FDSTAG_InterpolateVelocityPressureToCornerpoints(UserContext *user, Vec sol_vel);

PetscErrorCode ComputeRHS_FDSTAG(DM da, Vec b, UserContext *user);

PetscErrorCode SetBoundaryConditionsRHS_FDSTAG( DM da, UserContext *user, Vec b, PetscInt SetNonZeroValuesFlag);

PetscErrorCode FDSTAG_SetInitialViscosityFields(UserContext *user);

PetscErrorCode FDSTAG_SetInitialSolutionVector(UserContext *user, Vec sol_vel, Vec Pressure);

PetscErrorCode ComputeStiffnessMatrixRHSTemperature_FDSTAG(DM da_temp, DM da, Mat T_MAT, UserContext *user,
		Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt );

#endif
