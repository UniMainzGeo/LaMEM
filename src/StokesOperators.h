
#ifndef __StokesOperators_h__
#define __StokesOperators_h__


void ComputeKinmtx( const PetscInt _nnel, const PetscInt _edof, double KinMtx[6][MAX_edof], double **dhdPhys, double **dhdsVel, double **InvJacob );

void ComputeMaterialMatrix( const double mu, double MatMtx[6][6] );

void Multiply_Kin_Times_Mat_Times_Kin( const PetscInt _edof, double VV_add[MAX_edof][MAX_edof], double Kinmtx[6][MAX_edof], double Matmtx[6][6] );

PetscErrorCode ComputeInversePP( double PP[MAX_npres][MAX_npres], double InversePP[MAX_npres][MAX_npres] );

void VP_invPP_PV( const PetscInt _edof, const PetscInt _npres, double VV_correction[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres], double PP[MAX_npres][MAX_npres], double PV[MAX_npres][MAX_edof] );

PetscErrorCode ComputeDivergence( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Divergence, PetscScalar *MaximumDivergence);

PetscErrorCode ComputeGradient( LaMEMVelPressureDA C, DM da, DM da_pres, Vec Divergence, Vec b, Vec b_local, UserContext *user);

/* matrix free constructors */
PetscErrorCode MatCreate_MFOperatorDivergence( LaMEMVelPressureDA C, DM da_vel,Vec vel, DM da_pres,Vec pres, Mat *A );
PetscErrorCode MatCreate_MFOperatorGradient( LaMEMVelPressureDA C, DM da_vel,Vec vel, DM da_pres,Vec pres, UserContext *data, Mat *A );

void test_MFOperatorDivergence( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Divergence, PetscScalar *MaximumDivergence );
void test_MFOperatorGradient( LaMEMVelPressureDA C, DM da_vel, DM da_pres, Vec pressure, Vec vel, Vec vel_local, UserContext *user );

void LaMEMSetPressueElemDataMemoryFromArray(
		PressureElemDynamic *data,
		const PetscInt i, const PetscInt j, const PetscInt k,
		const PetscInt block_size,
		PetscScalar ***list );

PetscErrorCode MatDestroy_MFOperatorDivergence( Mat A );
PetscErrorCode MatMult_MFOperatorDivergence( Mat A, Vec x, Vec y );
PetscErrorCode MatDestroy_MFOperatorGradient( Mat A );
PetscErrorCode MatMult_MFOperatorGradient( Mat A, Vec x, Vec y );


#endif
