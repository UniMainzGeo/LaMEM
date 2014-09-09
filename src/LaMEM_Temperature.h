
#ifndef __LaMEM_Temperature_h__
#define __LaMEM_Temperature_h__


void ComputeKinmtxTemp( const PetscInt _nnel, const PetscInt _edof_temp, double **KinMtx, double **dhdPhys, double **dhdsVel,  double **InvJacob );

void ComputeMaterialMatrixTemp( const double k, double MatMtx[3][3] );

void Multiply_Kin_Times_Mat_Times_Kin_Temp( const PetscInt _edof_temp, double TT_add[MAX_edof_temp][MAX_edof_temp], double **Kinmtx, double Matmtx[3][3] );

void Multiply_ShapeT_Times_ShapeT( const PetscInt _edof_temp, double ShapeTemp[], double TT_mass[MAX_edof_temp][MAX_edof_temp] );

PetscErrorCode CreateStencilInGlobalStiffnessTemp( const PetscInt _edof_temp, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );

PetscErrorCode GetTemperatureElement( PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );

PetscErrorCode SetValuesRHS_Temp( PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );

PetscErrorCode ComputeStiffnessMatrixRHSTemperature( LaMEMVelPressureDA C, DM da_temp, DM da, Mat J, UserContext user,
		Vec Temp_local, Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt );

PetscErrorCode ComputeStiffnessMatrixRHSTemperature_FEM( LaMEMVelPressureDA C, DM da_temp, DM da, Mat J, UserContext *user,
		Vec Temp_local, Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt );

PetscErrorCode FindNearestNode( PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );

PetscErrorCode ParticleTemperatureToNodes( LaMEMVelPressureDA C, UserContext *user, DM da_temp, Vec Temp_local, Vec Temp );

PetscErrorCode TemperatureDiffNodesToParticles( LaMEMVelPressureDA C, UserContext *user, DM da_temp, Vec Temp, Vec Temp_old );

/*==========================================================================================================*/
/* Get global row number, given i,j and k */
static inline void GetGlobalIndexTemp( PetscInt nx, PetscInt ny, PetscInt i, PetscInt j, PetscInt k, PetscInt dof, PetscInt *ind )
{
	PetscInt totdof = 1;

	*ind = totdof*(k*(nx*ny) + j*nx + i) + dof;
}
/*==========================================================================================================*/


#endif
