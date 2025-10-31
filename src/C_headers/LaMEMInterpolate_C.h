#ifndef LAMEMINTERPOLATE_C_H
#define LAMEMINTERPOLATE_C_H

#include <petsc.h>


#ifdef __cplusplus
extern "C" {
#endif

// InterpFlags structure wrapper
typedef struct {
    PetscBool use_bound;
    PetscBool update;
} LaMEMInterpFlags;

// Thin wrappers around existing interpolation functions
PetscErrorCode LaMEMInterpolate_XFaceCorner(void *fs, Vec XFace, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_YFaceCorner(void *fs, Vec YFace, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_ZFaceCorner(void *fs, Vec ZFace, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_CenterCorner(void *fs, Vec Center, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_XYEdgeCorner(void *fs, Vec XYEdge, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_XZEdgeCorner(void *fs, Vec XZEdge, Vec Corner, LaMEMInterpFlags iflag);
PetscErrorCode LaMEMInterpolate_YZEdgeCorner(void *fs, Vec YZEdge, Vec Corner, LaMEMInterpFlags iflag);

#ifdef __cplusplus
}
#endif

#endif