#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_


PetscErrorCode JacResGetViscMat(PMat pm);

PetscErrorCode CopyViscosityToScalingVector(Vec a, Vec b, Vec c, Vec ScalingVec);

#endif /* SRC_BFBT_H_ */
