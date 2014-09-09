
#ifndef __LaMEM_Initialize_h__
#define __LaMEM_Initialize_h__


PetscErrorCode  LaMEM_Initialize_StokesSolver				( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C);

PetscErrorCode  LaMEM_Initialize_StokesSolver_FEM			( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C);

PetscErrorCode  LaMEM_Initialize_StokesSolver_FDSTAG		( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C);

PetscErrorCode  LaMEM_Initialize_TemperatureSolver			( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C);

PetscErrorCode  LaMEM_Initialize_TemperatureSolver_FEM		( UserContext *user, DAVPElementType vpt_element_type, LaMEMVelPressureDA C);

PetscErrorCode  LaMEM_Initialize_TemperatureSolver_FDSTAG	( UserContext *user, DAVPElementType vpt_element_type);

PetscErrorCode  LaMEM_Initialize_CheckInput					( void );

PetscErrorCode InitializeInternalErosionSurfaceOnRankZero	( UserContext *user);

PetscErrorCode SetSinusoidalPerturbation(PetscScalar SinusoidalFreeSurfaceAmplitude, UserContext *user);

PetscErrorCode FDSTAG_Preallocation_A11FDSTAG(DM daproc,DM da,DM da_pres,Mat A, PetscInt nvel_rowStart, PetscInt nvel_rowEnd);

PetscErrorCode FDSTAG_Preallocation_A12FDSTAG(DM daproc,DM da,DM da_pres,Mat A,
		PetscInt nvel_rowStart,  PetscInt nvel_rowEnd, PetscInt npres_rowStart, PetscInt npres_rowEnd);

PetscErrorCode FDSTAG_Preallocation_A21FDSTAG(DM daproc,DM da,DM da_pres,Mat A,
	PetscInt nvel_rowStart,  PetscInt nvel_rowEnd, PetscInt npres_rowStart, PetscInt npres_rowEnd);

#endif
