/*
 * GetGravityField.h
 *
 *  Created on: 12.02.2012
 *      Author: tobibaumann
 */
#ifndef GETGRAVITYFIELD_H_
#define GETGRAVITYFIELD_H_


typedef struct {
	PetscScalar 	xs,xm,ys,ym,zs,zm;
	PetscScalar		dx,dy,dz,dxh,dyh,dzh,cellvol;
} grid3D_struct;

typedef struct {
	PetscInt		nx,ny;
	PetscScalar 	xs,xm,ys,ym,zs,zm;
	PetscScalar		dx,dy;
	PetscScalar		x,y,z;
} survey_struct;



PetscErrorCode GetAiryIsostasy(UserContext *user, const PetscInt ngp_vel,grid3D_struct	grid3D, PetscInt itime);
PetscErrorCode GetMisfit_IsostaticTopography(UserContext *user,Vec gvec_Tiso, MPI_Comm PLANE_COMM);
PetscErrorCode GetGravityField(UserContext *user, const PetscInt _ngp_vel,PetscInt itime);
PetscErrorCode GetGravityEffectNumerical(PetscScalar dV,PetscInt num_int,PetscScalar x_survey,PetscScalar y_survey ,PetscScalar z_survey,PetscScalar **C ,PetscScalar *result);
PetscErrorCode GetGravityEffectAnalytical(PetscScalar x_survey,PetscScalar y_survey ,PetscScalar z_survey,PetscScalar **cornervec,PetscScalar *result);
PetscErrorCode GetGravityEffectAnalytical2(PetscScalar *x,PetscScalar *y,PetscScalar *z,PetscScalar *gsum);
PetscErrorCode GetMisfit_GravityField(UserContext *user, Vec gvec_survey_dg);
PetscErrorCode SaveGravityField2VTK(UserContext *user, Vec lvec_dg, Vec lvec_coords,PetscInt itime);


#endif /* GETGRAVITYFIELD_H_ */
