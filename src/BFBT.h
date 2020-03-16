#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_


PetscErrorCode JacResGetViscRes(JacRes *jr, PetscScalar dt);

PetscErrorCode JacResGetViscMat(JacRes *jr, PetscScalar dt);



#endif /* SRC_BFBT_H_ */
