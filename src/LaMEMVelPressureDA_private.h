

#ifndef __LaMEMVelPressureDA_private_h__
#define __LaMEMVelPressureDA_private_h__

// below is the good example of how man can overcomplicate things
// use C++ if you want proper classes, or just use functions & classical C

struct _p_LaMEMVelPressureDA {
	DAVPElementType type;
	char *type_name;
	char *type_description;

	// where are comments, by the way?
	
	PetscInt nnel;
	PetscInt ngp_vel; // what does this variable mean, in particular? asks naively-minded new user of LaMEM  
	PetscInt nintp_1D;
	PetscInt ElementType;
	PetscInt npres;
	PetscInt edof;
	PetscInt edof_temp;
	PetscInt nnode_el_1D;
	PetscInt nnel_1D;
	PetscInt ndim;
	
	PetscInt boundary_nodes_left[9];
	PetscInt boundary_nodes_right[9];
	PetscInt boundary_nodes_front[9];
	PetscInt boundary_nodes_back[9];
	PetscInt boundary_nodes_lower[9];
	PetscInt boundary_nodes_upper[9];
	
	/* get element index */
	void (*fp_ComputeElementIndex)( LaMEMVelPressureDA, Particles, PetscInt*, PetscInt*, PetscInt* );
	
	/* Material.c */ // can be replaced with a interpolation type function
	
	/* Temperature.c */
	void (*fp_CreateStencilInGlobalStiffnessTemp)( LaMEMVelPressureDA, MatStencil*,MatStencil*,const PetscInt,const PetscInt, const PetscInt );
	void (*fp_GetTemperatureElement)(LaMEMVelPressureDA, PetscScalar***, PetscScalar*, PetscInt, PetscInt, PetscInt);
	void (*fp_SetValuesRHS_Temp)( LaMEMVelPressureDA, PetscScalar***, PetscScalar*, PetscInt, PetscInt, PetscInt);
	void (*fp_FindNearestNode)( LaMEMVelPressureDA, PetscScalar, PetscScalar, PetscScalar, PetscInt*, PetscInt*, PetscInt*);
	
	/* Assembly.c */
	void (*fp_ComputeVelocityLocal2Global)( LaMEMVelPressureDA, MatStencil*,MatStencil*,PetscInt*, const PetscInt,const PetscInt, const PetscInt, Mat);
	void (*fp_CreateStencilInGlobalRhs)( LaMEMVelPressureDA, PetscInt*,const PetscInt,const PetscInt, const PetscInt, const PetscInt, const PetscInt, const PetscInt );
	void (*fp_SetValuesRHS)( LaMEMVelPressureDA, Field ***, PetscScalar*, PetscInt, PetscInt, PetscInt );
	
	/* Elements.c */
	PetscErrorCode (*fp_ComputeShapeFunctionVelocity)( LaMEMVelPressureDA, double*, double**, const double* );
	void (*fp_ComputeShapeFunctionPressure)( LaMEMVelPressureDA, double*, const DMDACoor3d, double* );
	PetscErrorCode (*fp_GetVelocityElement)( LaMEMVelPressureDA, Field***, PetscScalar*, PetscInt, PetscInt, PetscInt );
	PetscErrorCode (*fp_GetElementCoords)( LaMEMVelPressureDA, DMDACoor3d*, DMDACoor3d***, PetscInt,PetscInt, PetscInt, PetscInt );
};
//--------------------------------------------------------------------------------------------------------------------------------

void ComputeElementIndex_FDSTAG( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz );
void CreateStencilInGlobalStiffnessTemp_FDSTAG( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );
void GetTemperatureElement_FDSTAG( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );
void SetValuesRHS_Temp_FDSTAG( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );
void FindNearestNode_FDSTAG( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );
void ComputeVelocityLocal2Global_FDSTAG(LaMEMVelPressureDA C, MatStencil *row, MatStencil *col, PetscInt local2global[], const PetscInt i ,const PetscInt j , const PetscInt k, Mat MATRIX );
void CreateStencilInGlobalRhs_FDSTAG( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz );
void SetValuesRHS_FDSTAG( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode ComputeShapeFunctionVelocity_FDSTAG( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] );
void ComputeShapeFunctionPressure_FDSTAG( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );
PetscErrorCode GetVelocityElement_FDSTAG( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode GetElementCoords_FDSTAG( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag );
PetscErrorCode LaMEMVelPressureDACreate_FDSTAG( LaMEMVelPressureDA C );

void ComputeElementIndex_Q1P0( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz );
void CreateStencilInGlobalStiffnessTemp_Q1P0( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );
void GetTemperatureElement_Q1P0( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );
void SetValuesRHS_Temp_Q1P0( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );
void FindNearestNode_Q1P0( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );
void ComputeVelocityLocal2Global_Q1P0(LaMEMVelPressureDA C,MatStencil *row, MatStencil *col,PetscInt local2global[],const PetscInt i ,const PetscInt j , const PetscInt k,Mat MATRIX );
void CreateStencilInGlobalRhs_Q1P0( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz );
void SetValuesRHS_Q1P0( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode ComputeShapeFunctionVelocity_Q1P0( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] );
void ComputeShapeFunctionPressure_Q1P0( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );
PetscErrorCode GetVelocityElement_Q1P0( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode GetElementCoords_Q1P0( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag );
PetscErrorCode LaMEMVelPressureDACreate_Q1P0( LaMEMVelPressureDA C );

void ComputeElementIndex_Q1Q1( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz );
void CreateStencilInGlobalStiffnessTemp_Q1Q1( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );
void GetTemperatureElement_Q1Q1( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );
void SetValuesRHS_Temp_Q1Q1( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );
void FindNearestNode_Q1Q1( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );
void ComputeVelocityLocal2Global_Q1Q1(LaMEMVelPressureDA C,MatStencil *row, MatStencil *col,PetscInt local2global[],const PetscInt i ,const PetscInt j , const PetscInt k,Mat MATRIX );
void CreateStencilInGlobalRhs_Q1Q1( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz );
void SetValuesRHS_Q1Q1( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode ComputeShapeFunctionVelocity_Q1Q1( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] );
void ComputeShapeFunctionPressure_Q1Q1( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );
PetscErrorCode GetVelocityElement_Q1Q1( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode GetElementCoords_Q1Q1( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag );
PetscErrorCode LaMEMVelPressureDACreate_Q1Q1( LaMEMVelPressureDA C );

void ComputeElementIndex_Q2PM1( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz );
void CreateStencilInGlobalStiffnessTemp_Q2PM1( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );
void GetTemperatureElement_Q2PM1( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );
void SetValuesRHS_Temp_Q2PM1( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );
void FindNearestNode_Q2PM1( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );
void ComputeVelocityLocal2Global_Q2PM1( LaMEMVelPressureDA C,MatStencil *row, MatStencil *col, PetscInt local2global[],const PetscInt i ,const PetscInt j , const PetscInt k, Mat MATRIX );
void CreateStencilInGlobalRhs_Q2PM1( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz );
void SetValuesRHS_Q2PM1( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode ComputeShapeFunctionVelocity_Q2PM1( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] );
void ComputeShapeFunctionPressure_Q2PM1_GLOBAL( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );
void ComputeShapeFunctionPressure_Q2PM1_LOCAL( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );
PetscErrorCode GetVelocityElement_Q2PM1( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k );
PetscErrorCode GetElementCoords_Q2PM1( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag );
PetscErrorCode LaMEMVelPressureDACreate_Q2PM1( LaMEMVelPressureDA C );

//--------------------------------------------------------------------------------------------------------------------------------

#endif

