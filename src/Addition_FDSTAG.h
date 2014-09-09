
#ifndef __Addition_FDSTAG_h__
#define __Addition_FDSTAG_h__

PetscErrorCode AddPushingToModelPV_MAT(UserContext *user,Mat PV_MAT,Mat PV_MAT_push,Vec pv_rhs_push);
PetscErrorCode AddPushingToModel(UserContext *user, Mat VV_MAT,Mat VP_MAT,Vec rhs,Vec rhs_p, Vec pv_rhs_push);
PetscErrorCode GetGlobalIndicesForPushing(UserContext *user, PetscInt *rowidx_array, PetscInt *numRowsx, PetscInt *rowidy_array, PetscInt *numRowsy, PetscBool flg);
PetscErrorCode ModifyStiffnessMatrixForPushing(UserContext *user,Mat VV_MAT,Mat VP_MAT,Vec rhs,PetscInt *rowidx_array, PetscInt numRowsx,PetscInt *rowidy_array, PetscInt numRowsy,PetscInt ind_change);
PetscErrorCode ModifySolutionsVectorForPushing2(UserContext *user,Vec x_push,PetscInt *rowidx_array,PetscInt numRowsx,PetscInt *rowidy_array,PetscInt numRowsy,PetscInt ind_change);
PetscErrorCode ZeroColumnsPV_MAT(UserContext *user,Vec pv_rhs_push,Mat PV_MAT,Mat PV_MAT_push,PetscInt *rowidx_array,PetscInt numRowsx,PetscInt *rowidy_array,PetscInt numRowsy,PetscInt ind_change);
PetscErrorCode AdvectThePushingBlockCoordinates(UserContext *user, PetscInt ind_change);

#endif
