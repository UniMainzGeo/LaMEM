/*
 * GetSurfaceVelocity.h
 *
 *  Created on: 13.02.2012
 *      Author: tobibaumann
 */

#ifndef GETSURFACEVELOCITY_H_
#define GETSURFACEVELOCITY_H_

// free slip box surface
PetscErrorCode GetSurfaceVelocity(UserContext *user,DM da_nodes,Vec gvec_Vel);
PetscErrorCode GetMisfit_SurfaceVelocity(UserContext *user,Vec lvec_Vx, Vec lvec_Vy,MPI_Comm SURF_COMM);
PetscErrorCode Sum_Vectors(MPI_Comm G_COMM ,Vec *lvec_Vall,PetscInt vec_length);
PetscErrorCode VecSeq2VecMpi(PetscMPIInt rank,Vec lvec,Vec *gvec);
PetscErrorCode VecAbsMax(Vec gvec,PetscScalar *max);

// internal free surface
PetscErrorCode GetSurfaceVelocity_InternalFreeSurface(UserContext *user,PetscInt itime);
PetscErrorCode GetMisfit_SurfaceVelocity_InternalFreeSurface(UserContext *user,Vec gvec_Vx, Vec gvec_Vy,Vec gvec_Vz,MPI_Comm PLANE_COMM);

#endif /* GETSURFACEVELOCITY_H_ */
