
#ifndef __Utils_h__
#define __Utils_h__

PetscErrorCode DACompareStructures( DM da1, DM da2, PetscBool *flg );

PetscErrorCode LaMEMMod( PetscInt x, PetscInt y, PetscInt *result );

PetscInt Mod( PetscInt x,PetscInt y );

PetscBool pnpoly(PetscInt nvert, PetscScalar *vertx, PetscScalar *verty, PetscScalar testx, PetscScalar testy);

PetscScalar InterpolateWithin2DLinearElement(PetscScalar *p_x,PetscScalar *p_y,PetscScalar p_z[4], PetscScalar x, PetscScalar y);

PetscErrorCode NaturalCoords_Linear2D_Element( PetscScalar *eta, PetscScalar *xsi,PetscScalar x_real, PetscScalar y_real, PetscScalar ECOORD_x[4], PetscScalar ECOORD_y[4]);

PetscErrorCode DAGetProcessorSubset_VerticalDirection(DM da, PetscInt *NumProcs_Z, MPI_Comm *comm);

PetscErrorCode LaMEMReadInputFile( UserContext *user );

PetscErrorCode InitializeCode( UserContext *user );

PetscErrorCode  StencilToLocalNumbering(Mat,PetscInt,const MatStencil[],PetscInt[]);

PetscErrorCode ComputeGlobalProperties( DM da, UserContext *user, PetscInt itime, Vec Velocity, LaMEMVelPressureDA C );

PetscErrorCode DASetGhostedCoordinates(DM da,Vec c);

PetscErrorCode WriteStiffnessMatrixToDisk( DM da, Mat A11, Mat A12, Mat A21, Mat A22, Mat approx_S, Vec f, Vec h, Vec Sol_vel, Vec Pressure, UserContext *user);

PetscErrorCode WriteTemperatureStiffnessMatrixToDisk( DM da_temp, Mat T_MAT, Vec rhs_temp, Vec Temp, UserContext *user);

PetscErrorCode LaMEMCreate2dArray( const PetscInt M, const PetscInt N, PetscScalar ***_A2, PetscScalar **_A );

PetscErrorCode LaMEMDestroy2dArray(PetscScalar ***_A2, PetscScalar **_A );

PetscErrorCode GetMaterialPropertiesFromCommandLine(UserContext *user);

PetscErrorCode GetLithColumnFromCommandLine(UserContext *user);
//---------------------------------------------------------------------------
PetscErrorCode DMDAView(const char * name, DM da);

PetscErrorCode DMDAViewVTK(const char * filename, DM da);

PetscErrorCode DMDAGetProcessorRank(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_col);

PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n);

PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n);

//---------------------------------------------------------------------------

PetscErrorCode DebugSave2Bin_GlobalVec(Vec vec_data,const char *filename,PetscInt itime, const char *folder);

PetscErrorCode DebugSave2Bin_GlobalMat(Mat mat_data,const char *filename,PetscInt itime, const char *folder);

//---------------------------------------------------------------------------
PetscErrorCode DMDACreate2dFrom3d(DM da,PetscInt gp,DMDADirection dir,DM *da_plane_out,MPI_Comm *PLANE_COMM);
PetscErrorCode DMDADestroy2dFrom3d(DM da_plane, MPI_Comm PLANE_COMM);

PetscErrorCode DMExtractGlobalVec2dFromGlobalVec3d(DM da,Vec gvec_da,PetscInt gp,DM da_plane,MPI_Comm PLANE_COMM,Vec *gvec_plane);
PetscErrorCode DMDestroyGlobalVec2dFromGlobalVec3d(Vec gvec_plane, MPI_Comm PLANE_COMM);


PetscErrorCode SaveProcessorPartitioning(UserContext *user);

//---------------------------------------------------------------------------

PetscErrorCode  LaMEMInitializeMaterialProperties( UserContext *user );

//---------------------------------------------------------------------------
// Get global row number, given i,j and k
static inline void GetGlobalIndex( PetscInt nx, PetscInt ny, PetscInt i, PetscInt j, PetscInt k, PetscInt dof, PetscInt *ind)
{
	PetscInt totdof = 3;
	*ind = totdof*(k*(nx*ny) + j*nx + i) + dof;
}
//---------------------------------------------------------------------------
PetscErrorCode ReadMeshSegDir(
	FILE        *fp,
	const char  *name,
	PetscScalar  beg,
	PetscScalar  end,
	PetscInt    *tncels,
	MeshSegInp  *msi);
//---------------------------------------------------------------------------

#endif

