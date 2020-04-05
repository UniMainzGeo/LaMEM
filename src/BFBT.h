#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_


PetscErrorCode JacResGetViscRes(JacRes *jr, PetscScalar dt);

PetscErrorCode JacResGetViscMat(JacRes *jr, PetscScalar dt);

PetscErrorCode LumpMatrixToVector(PMat pm);

#endif /* SRC_BFBT_H_ */
