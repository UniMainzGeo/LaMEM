/*
 * phase_transition.h
 *
 *  Created on: Apr 20, 2020
 *      Author: piccolo
 */

#ifndef melt_parametrization_h_
#define melt_parametrization_h_
struct ConstEqCtx;
struct Scaling;
struct DBMat;
struct Material_t;

// read phase transition law
PetscScalar  Compute_Melt_Fraction(PetscScalar   P, PetscScalar T ,Material_t *phase, ConstEqCtx *ctx);



#endif /* MELT_PARAMETRIZATION_H_ */
