/*
 * BFBT.h
 *
 *  Created on: 04.03.2020
 *      Author: daniel
 */

#ifndef SRC_BFBT_H_
#define SRC_BFBT_H_


PetscErrorCode GetViscMat(PMat pm);

PetscErrorCode CreateViscMat(PMat pm);

PetscErrorCode CopyViscosityToScalingVector(Vec a, Vec b, Vec c, Vec ScalingVec);


#endif /* SRC_BFBT_H_ */
