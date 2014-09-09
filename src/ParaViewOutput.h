
#ifndef __LaMEM_ParaViewOutput_h__
#define __LaMEM_ParaViewOutput_h__


typedef struct {
	PetscBool vtk_ascii;
	PetscBool mu,rho,G; /* scaled by non-dimensional paramaters */
	PetscBool numParticles;
	PetscBool pressure, temperature;
	PetscBool devStrainrateInv;
	PetscBool devStressInv;
	PetscBool stress,strainrate;
	PetscBool strain,plasticStrain;
	PetscBool phase, n,C,phi,k,Cp,Q,alpha,FK;
} LaMEMView_QuadratureFields;




PetscErrorCode DAView3DPVTS(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length);
PetscErrorCode DAView2DPVTS(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length);

PetscErrorCode DMView_3DVTK_StructuredGrid( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[],PetscScalar scaling_length );
PetscErrorCode DMView_3DVTK_StructuredGrid_3ComponentVector( DM da, Vec FIELD, const char file_prefix[], const char DirectoryName[], const char vectorfield_name[],PetscScalar scaling_length);
PetscErrorCode DMView_2DVTK_StructuredGrid_Topo( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[], Vec SurfaceTopography_Vx, Vec SurfaceTopography_Vy, Vec SurfaceTopography_Vz, PetscScalar scaling_length );
PetscErrorCode DMView_2DVTK_StructuredGrid_Topo_Erosion( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[], PetscScalar scaling_length );


PetscErrorCode ParaviewPVDOpen(const char pvdfilename[]);
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[]);

PetscErrorCode DMViewVTK_write_PieceExtend( FILE *vtk_fp, PetscInt indent_level, DM da, const char local_file_prefix[] );
PetscErrorCode LaMEM_DMView_3DVTK_StructuredGrid_QuadPoints(LaMEMView_QuadratureFields *view, DM DA_Materials_fine, DM DA_Processors,DM DA_Quadrature, Vec Materials_fine, const char file_prefix[], PetscInt ngp_vel_1D, UserContext *user, const char DirectoryName[] );

PetscErrorCode  LaMEM_CreateOutputDirectory(const char *DirectoryName);
PetscErrorCode LaMEM_SetQuadraturePointCoords_to_QuadratureDA(DM DA_Quadrature, PetscInt ngp_vel_1D, DM DA_Materials_fine, Vec Materials_fine );
PetscErrorCode LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DM DA_Quadrature, PetscInt ngp_vel_1D, DM DA_Materials_fine, Vec Materials_fine, Vec DataArray, UserContext *user, PetscInt DataType );
PetscErrorCode LaMEM_DAView_3DVTK_StructuredGrid_QuadPoints(LaMEMView_QuadratureFields *view, DM DA_Materials_fine, DM DA_Processors, DM DA_Quadrature, Vec Materials_fine, const char file_prefix[], PetscInt ngp_vel_1D, UserContext *user, const char DirectoryName[]);

PetscErrorCode LaMEMView_QuadratureFieldsInit( LaMEMView_QuadratureFields *view );
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTU( LaMEMView_QuadratureFields *view, UserContext *user, DM DA_Processors_fine, DM DA_Materials_fine, Vec Materials_fine, PetscInt ngp_vel, const char NAME[] );
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTS( LaMEMView_QuadratureFields *view,  UserContext *user, DM DA_Processors_fine, DM DA_Materials_fine, Vec Material_fine, PetscInt ngp_vel_1D, PetscInt itime, const char DirectoryName[] );
PetscErrorCode LaMEMViewQuadraturePoints_Meta3DPVTS( LaMEMView_QuadratureFields *view, UserContext *user, const char file_prefix[], const char local_file_prefix[], DM da, const char DirectoryName[]  );


PetscErrorCode WriteTemperatureOutputFile_VTS(UserContext *user,PetscInt itime, Vec Temp, const char DirectoryName[]);
PetscErrorCode WriteVelocityOutputFile_VTS(UserContext *user,PetscInt itime, Vec sol, const char DirectoryName[]);
PetscErrorCode WriteTopographyOutputFile_VTS(UserContext *user,PetscInt itime, const char DirectoryName[]);
PetscErrorCode WritePhasesOutputFile_VTS(LaMEMVelPressureDA C, UserContext *user,PetscInt itime, const char DirectoryName[]);

PetscErrorCode DMView_2DVTK_PStructuredGrid( DM da, const char file_prefix[], const char local_file_prefix[], UserContext *user );
PetscErrorCode DMView_3DVTK_PStructuredGrid( DM da, const char file_prefix[], const char local_file_prefix[], UserContext *user );
PetscErrorCode DMView2DPVTS_Topo(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length);
PetscErrorCode DMView3DPVTS(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length);
PetscErrorCode DMViewVTK_write_PieceExtend_Topo( FILE *vtk_fp, PetscInt indent_level, DM da, const char local_file_prefix[] );
PetscErrorCode LaMEMViewQuadraturePoints_Meta3DPVTU( LaMEMView_QuadratureFields *view, const char file_prefix[], const char local_file_prefix[]);
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTU_appended( LaMEMView_QuadratureFields *view, UserContext *user, PetscInt xm, PetscInt ym, PetscInt zm, DM DA_Materials_fine, Vec Materials_fine, PetscInt ngp_vel, const char file_prefix[] );



#endif
