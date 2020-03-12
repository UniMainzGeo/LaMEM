#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_









PetscErrorCode JacResGetViscParam(JacRes *jr, PetscScalar *phRat, PetscScalar *visc);

PetscErrorCode JacResGetViscMat(JacRes *jr, PetscScalar dt);



#endif /* SRC_BFBT_H_ */
