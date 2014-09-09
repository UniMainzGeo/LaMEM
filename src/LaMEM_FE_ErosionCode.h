
#ifndef __LaMEM_FD_ErosionCode_h__
#define __LaMEM_FD_ErosionCode_h__

/*------------------------STRUCTURES-----------------------------*/

// Structure that contains all information we need in the code
typedef struct {
	PetscInt 		M,N,nt,output_freq;
    PetscInt        nelx, nely, n_int,fill_lake;
	PetscInt		BC, nbre_river;
    PetscScalar     location_river[100];
    PetscScalar     rain_river; 
    PetscInt        mode_river;
	PetscScalar	 	lx,ly, Hmax, noise, dt, rain, k0, c, hack;
	DM				da_elevation, da_element, da_elevation_output;
	Vec				Elevation, Elevation_init, Elevation_el, Streamflow_el, properties;

} FE_Erosion_UserData;


// structure containing data for the elements
typedef struct {
    PetscScalar H_el, k_el, Q_el, s_el;
} FE_ElementInfo;


typedef struct  {
    PetscScalar H;
    PetscInt num;
    PetscInt low;
} NodeData;

/*-----------------------------------------------------------------*/

// Functions for the erosion code
PetscErrorCode FE_ErosionCode_TectonicTimestep( UserContext *user);
PetscErrorCode TransferLaMEMDataToAndInitialize_FD_ErosionCode( UserContext *user, FE_Erosion_UserData *Erosion_UserData, PetscScalar dt);
PetscErrorCode FE_Erosion_UpdateTopography_DA( FE_Erosion_UserData *Erosion_UserData);

PetscErrorCode Solve_NonlinearDiffusion( FE_Erosion_UserData *Erosion_UserData, PetscInt bcdof[], PetscScalar bcval[]);
PetscErrorCode GetShapeFunction( PetscScalar coord_int0 ,PetscScalar coord_int1 , PetscScalar Ni[], PetscScalar dNi[2][4]);
PetscErrorCode ConstructDetJacobian( PetscScalar dNi[2][4], PetscScalar coords[], PetscScalar *det_J);
PetscErrorCode ConstructdNdx( PetscScalar dNi[2][4], PetscScalar coords[], PetscScalar dNdx[2][4]);
PetscErrorCode Evaluate_Integrale( PetscInt n_int, const PetscScalar k_el, const PetscScalar s_el, PetscScalar coord_el[], PetscScalar MM[], PetscScalar KM[], PetscScalar F[]);
PetscErrorCode DAGetElementCorners( DM da, PetscInt *sx, PetscInt *sy,PetscInt *mx, PetscInt *my,PetscInt *mz);
PetscErrorCode DAGetLocalElementSize( DM da, PetscInt *mxl, PetscInt *myl,PetscInt *mzl);
PetscErrorCode ConstructRow(const PetscInt m, const PetscInt istart, const PetscInt jstart, const PetscInt ei, const PetscInt ej, const PetscInt iglobal[], PetscInt row[]);
PetscErrorCode Assemble_matrixes_vector(FE_Erosion_UserData *Erosion_UserData, Mat KLG, Mat KRG, Vec FG);
PetscErrorCode GetElementCoord(PetscScalar coord_el[], PetscInt i, PetscInt j, DMDACoor2d **coors);
PetscErrorCode DAGetElementOwnershipRanges( DM da, PetscInt **_lx, PetscInt **_ly );
PetscErrorCode  Set_Diffusion_Coefficient( FE_Erosion_UserData *Erosion_UserData);
PetscErrorCode  NodeDataComparison(const void *a, const void *b);
PetscErrorCode  NodeDataSort(PetscInt L,PetscScalar Hvec[],PetscInt low[],PetscInt indices[]);
PetscErrorCode  lamem_min(PetscScalar*,PetscInt,PetscInt*,PetscScalar*);
PetscErrorCode  low_neigh(PetscScalar Hvec[], PetscInt low [], PetscInt M, PetscInt N);
PetscErrorCode  area_node_reggrid(PetscScalar area[], PetscScalar lx, PetscScalar ly, PetscInt M, PetscInt N);
PetscErrorCode stream_flow_D8(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar Qvec[], PetscScalar ly, PetscScalar Hvec[], PetscInt fill_lake);
PetscErrorCode  GetAverageNodes(PetscScalar table[], PetscInt i, PetscInt j, PetscScalar **matrix, PetscScalar *mean);
PetscErrorCode low_neigh_grad(PetscScalar H[], PetscInt low [], PetscInt M, PetscInt N, PetscScalar lx, PetscScalar ly);
PetscErrorCode GetInputParameters( FE_Erosion_UserData *Erosion_UserData);
PetscErrorCode InitializeMesh( FE_Erosion_UserData *Erosion_UserData);
PetscErrorCode FE_Erosion_UpdateTopography_DA( FE_Erosion_UserData *Erosion_UserData);
PetscErrorCode DefineBC(FE_Erosion_UserData *Erosion_UserData, PetscInt *bcdof, PetscScalar *bcval);
PetscErrorCode computeBC(PetscInt N, PetscInt M, PetscScalar *H, PetscInt *bcdof, PetscScalar *bcval, PetscInt BC);
PetscErrorCode TransferLaMEMDataToAndInitialize_FE_ErosionCode( UserContext *user, FE_Erosion_UserData *Erosion_UserData, PetscScalar dt_erosion);
PetscErrorCode SaveInitialErosionSurface( UserContext *user, const char *FileName);
PetscErrorCode LoadInitialErosionSurface( UserContext *user);
PetscErrorCode NodeNeighbour_LowestHeight(PetscInt Mx,PetscInt My,PetscScalar *H,PetscInt low[]);
PetscErrorCode Solve_linearDiffusion( FE_Erosion_UserData *Erosion_UserData, PetscInt bcdof[], PetscScalar bcval[]);
PetscErrorCode Assemble_matrixes_vector_smoother(FE_Erosion_UserData *Erosion_UserData, Mat KLG, Mat KRG, Vec FG);
PetscErrorCode stream_flow_D8_location_river(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar Qvec[], PetscScalar lx,PetscScalar ly, PetscScalar Hvec[], PetscScalar location_river[], PetscInt nbre_river);
PetscErrorCode choise_int(PetscScalar arrondi_sup, PetscScalar arrondi_inf, PetscScalar result, PetscInt *result_int);
PetscErrorCode stream_flow(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar rain_river, PetscScalar Qvec[], PetscScalar lx,PetscScalar ly, PetscScalar Hvec[], PetscScalar location_river[], PetscInt nbre_river, PetscInt fill_lake);

#endif
