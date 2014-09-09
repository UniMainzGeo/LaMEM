
#ifndef __Mesh_h__
#define __Mesh_h__

PetscErrorCode AdvectGrid( DM da, Vec Velocity, DM da_coord, UserContext *user);

PetscErrorCode ComputeNodesPerProcessor( DM da, PetscInt *lx, PetscInt *ly,PetscInt *lz, const DAVPElementType element_type );

PetscErrorCode SaveInitialMesh( UserContext *user, DM da_nodes,const char *FileName );

PetscErrorCode AdvectAndUpdate_InternalFreeSurface( UserContext *user, Vec Velocity );

PetscErrorCode ReadMeshFromFile( DM da, UserContext *user );

PetscErrorCode GenerateMeshFine( DM da, UserContext *user );

PetscErrorCode ComputeNeighbors( DM da, PetscInt NeighborCPU[3][3][3] );

PetscErrorCode ApplySedimentationToInternalFreeSurface( UserContext *user);

PetscErrorCode ApplyInfinitelyFastErosionToInternalFreeSurface( UserContext *user);

PetscErrorCode ApplySerialErosionCodeToInternalFreeSurface(UserContext *user );

PetscErrorCode InterpolateErosionSurfaceToInternalFreeSurface(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography);

PetscErrorCode AdvectInternalErosionSurfaceOnRankZeroWithTectonicVelocitiesAndReinterpolateToRegularMesh(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec FreeSurface_Vx, Vec FreeSurface_Vy, Vec FreeSurface_Vz );

PetscErrorCode ApplyErosion( DM da_coord, UserContext *user );

PetscErrorCode RemeshGrid( DM da, UserContext *user );

PetscErrorCode UpdateSurfaceAndBottomTopography_FEM( UserContext *user );

PetscErrorCode ScatterTopographyData_to_DA_Surface( DM da, Vec GlobalTopography, PetscInt SurfaceTopography);

PetscErrorCode DAUpdatedGhostedCoordinates( DM da );

PetscErrorCode DASetCoordinatesFromLocalVector( DM da, Vec local_coords );

PetscBool PetscCompareScalar(PetscScalar f1, PetscScalar f2);

PetscErrorCode DeformFDSTAGMeshWithBackgroundStrainrate( UserContext *user );

PetscErrorCode SetUniformCoordinates_FDSTAG( UserContext *user );

PetscErrorCode CopyInternalFreeSurfaceToRankZero(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography, Vec FreeSurface_Vx, Vec FreeSurface_Vy, Vec FreeSurface_Vz );

PetscErrorCode CopyErodedSurfaceFromRankZeroToAllOtherRanks(UserContext *user, DM DMDA_InternalFreeSurfaceOnRankZero, Vec InternalFreeSurfaceTopography );

double round_double( double d, PetscInt places );

PetscErrorCode DARoundCoordinates( DM da, const PetscInt places );

PetscErrorCode RemeshGrid_old( DM da, UserContext *user );

#endif

