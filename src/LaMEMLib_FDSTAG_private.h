//-----------------------------------------------------------------------------

#ifndef __LaMEMLib_FDSTAG_private_h__
#define __LaMEMLib_FDSTAG_private_h__

//-----------------------------------------------------------------------------

PetscErrorCode CreateSolutionVectors(UserContext *user);

PetscErrorCode DestroySolutionObjects(UserContext *user, LaMEMVelPressureDA *C);

// PetscErrorCode ViewLinearSolveResidual(UserContext *user);

PetscErrorCode CalculateMisfitValues(
	UserContext        *user,
	LaMEMVelPressureDA  C,
	PetscInt            itime,
	PetscScalar        *LaMEM_OutputParameters);

PetscErrorCode CalculateTimeStep(UserContext *user, PetscInt itime);

PetscErrorCode CheckVelocityError(UserContext *user);

PetscErrorCode FDSTAGCompPrecond(
		Mat VV_MAT, Mat VP_MAT, Mat PV_MAT, Mat PP_MAT, Mat approx_S,
		Vec ViscosityScaling, PetscScalar *nrmVV, UserContext *user);

//-----------------------------------------------------------------------------
#endif
