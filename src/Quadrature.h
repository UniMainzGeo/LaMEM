
#ifndef __Quadrature_h__
#define __Quadrature_h__

PetscErrorCode IntegrationPoints( const PetscInt _nnel, const PetscInt _nintp_1D, double IntPoint[3][MAX_ngp_vel], double IntWeight[] );

PetscErrorCode IntPointProperties( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Pressure,
		DM da_temp, Vec Temp, UserContext *user, PetscScalar dt, PetscInt ComputeFull );

#endif

