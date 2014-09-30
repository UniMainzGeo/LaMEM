//-----------------------------------------------------------------------------

#ifndef __LaMEMLib_FDSTAG_private_h__
#define __LaMEMLib_FDSTAG_private_h__

//-----------------------------------------------------------------------------

// PetscErrorCode ViewLinearSolveResidual(UserContext *user);

PetscErrorCode CalculateMisfitValues(
	UserContext        *user,
	LaMEMVelPressureDA  C,
	PetscInt            itime,
	PetscScalar        *LaMEM_OutputParameters);

PetscErrorCode CalculateTimeStep(UserContext *user, PetscInt itime);

PetscErrorCode CheckVelocityError(UserContext *user);

//-----------------------------------------------------------------------------
#endif
