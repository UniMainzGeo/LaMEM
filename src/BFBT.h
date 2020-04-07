#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_


PetscErrorCode JacResGetViscRes(JacRes *jr, PetscScalar dt);

PetscErrorCode JacResGetViscMat(PMat pm);

PetscErrorCode LumpMatrixToVector(PMat pm);

#endif /* SRC_BFBT_H_ */
