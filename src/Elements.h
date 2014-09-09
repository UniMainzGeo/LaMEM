
#ifndef __Elements_h__
#define __Elements_h__

PetscErrorCode ComputeShapeFunctionVelocity( double ShapeVel[], double **dhdsVel, const double Point[] );

PetscErrorCode ComputeShapeFunctionPressure( double ShapeP[], const DMDACoor3d CoordIntp, double Point[] );

PetscErrorCode ComputeJacobianElement( const PetscInt _nnel, double **dhdsVel, const DMDACoor3d coord_elem[], double **Jacob,double **InvJacob, double *DetJacob );

PetscErrorCode ComputeCoordIntp( const PetscInt _nnel, const double ShapeVel[], const  DMDACoor3d coord_elem[], DMDACoor3d *CoordIntp_out );

PetscErrorCode GetVelocityElement( Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k );

PetscErrorCode GetPressureElement( P_array ***pressure, PetscScalar P_element[], PetscInt i, PetscInt j, PetscInt k );

PetscErrorCode GetElementCoords( DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag);

PetscErrorCode ElementSpacing( DMDACoor3d *coord_elem, PetscScalar *dx, PetscScalar *dy, PetscScalar *dz);

PetscErrorCode FDSTAG_ComputeSpacing(DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, UserContext *user,
		PetscScalar dx_NormalVelocity[7], PetscScalar dy_NormalVelocity[7], PetscScalar dz_NormalVelocity[7],
		PetscScalar dx_P[2], PetscScalar dy_P[2], PetscScalar dz_P[2],  PetscScalar *z_center);

PetscErrorCode CorrectElementCoordsForPeriodicity( DMDACoor3d *coord_elem, PetscInt ix, PetscInt iy, UserContext *user, PetscInt Flag);


PetscErrorCode FDSTAG_ExtractMaterialParameters(PetscScalar ***viscosity_center, PetscScalar ***viscosity_XY, PetscScalar ***viscosity_YZ, PetscScalar ***viscosity_XZ,
		PetscScalar ***density_center, PetscScalar ***PhaseProportionsAir_Center,  PetscScalar ***LocalSurfaceTopography, PetscInt ielx,PetscInt iely, PetscInt ielz, PetscInt zs_FreeSurface,
		PetscScalar dx_P[2], PetscScalar dy_P[2], PetscScalar dz_P[2], UserContext *user,
		PetscScalar ViscosityCenter[7], PetscScalar Viscosity_XY[4], PetscScalar Viscosity_YZ[4],PetscScalar Viscosity_XZ[4], PetscScalar dRho_dxdydz[3], PetscInt FreeSurfaceCells[7],
		PetscScalar	*z_FreeSurface );



#endif
