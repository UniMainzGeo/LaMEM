
#ifndef __Assembly_h__
#define __Assembly_h__

PetscErrorCode ComputeVelocityLocal2Global( const PetscInt _edof, MatStencil *row,MatStencil *col, PetscInt local2global[], const PetscInt i ,const PetscInt j , const PetscInt k, Mat MATRIX );

PetscErrorCode ComputePressureLocal2Global( const PetscInt, PetscInt, PetscInt, PetscInt, PetscInt PRES_local2global[], Mat);

PetscErrorCode SetValuesRHS( Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k );

PetscErrorCode ComputeStiffnessMatrixes( LaMEMVelPressureDA C, DM da, DM da_pres,
		Mat VV_MAT, 		Mat PP_MAT, 	Mat PV_MAT, 	Mat VP_MAT,
		UserContext user, 	PetscScalar dt, Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec rhs_add_pBC, Vec pv_rhs_push, Vec ViscosityScaling, Mat approx_S, PetscInt *remesh );


PetscErrorCode ComputeStiffnessMatrixes_FEM( LaMEMVelPressureDA C, DM da, DM da_pres,
		Mat VV_MAT, 		Mat PP_MAT, 	Mat PV_MAT, 	Mat VP_MAT,
		UserContext *user, 	Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec rhs_add_pBC, Vec ViscosityScaling, Mat approx_S, PetscInt *remesh );

PetscErrorCode ComputeRHS( LaMEMVelPressureDA C, DM da, DM da_p, Vec b, Vec Pressure, UserContext user);

PetscErrorCode ComputeRHS_FEM( LaMEMVelPressureDA C, DM da, DM da_p, Vec b, Vec Pressure, UserContext *user);

PetscErrorCode SetBC_ElementStiffnessMatrix_Symmetric( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], double F_BC[],double P_BC[],
		DMDACoor3d coord_elem[],
		PetscInt iel_x, PetscInt iel_y, PetscInt iel_z,
		PetscInt nel_x, PetscInt nel_y, PetscInt nel_z, sBC BoundaryConditions, UserContext *user);

PetscErrorCode Add_BC_ToLocalStiffness(LaMEMVelPressureDA C,double VV[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], double F_BC[], double P_BC[],
		PetscInt nodes[], PetscInt ndof, PetscInt SetDOF,
		double BoundaryVelocity[MAX_nnel_1D*MAX_nnel_1D][3] );

PetscErrorCode SetBoundaryConditionsRHS( DM da, UserContext *user, Vec b, PetscInt SetNonZeroValuesFlag );

PetscErrorCode ComputeElementStiffnessMatrixes( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], 		double PP[MAX_npres][MAX_npres], double PP_prec[MAX_npres][MAX_npres],
		double IntPoint[MAX_ndim][MAX_ngp_vel],	double IntWeight[],	double *mu,
		DMDACoor3d coord_elem[], PetscScalar *DiagonalRatio );

PetscErrorCode WriteElementStiffnessMatrixToDisk( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], 			double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], 					double V_RHS[], 			double PP[MAX_npres][MAX_npres],
		double mu_average,
		double rho_g, 								DMDACoor3d coord_elem[], 		PetscInt number,
		PetscInt ielx, 								PetscInt iely, 					PetscInt ielz);


#endif
