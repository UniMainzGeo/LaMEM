#include "LaMEMInterpolate_C.h"
#include "../LaMEM.h"
#include "../interpolate.h"
#include "../fdstag.h"

extern "C" {

// Convert C struct to C++ struct
static InterpFlags convertInterpFlags(LaMEMInterpFlags c_flags) {
    InterpFlags cpp_flags;
    cpp_flags.use_bound = c_flags.use_bound;
    cpp_flags.update = c_flags.update;
    return cpp_flags;
}

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMInterpolate_XFaceCorner(void *fs, Vec XFace, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpXFaceCorner((FDSTAG*)fs, XFace, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_YFaceCorner(void *fs, Vec YFace, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpYFaceCorner((FDSTAG*)fs, YFace, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_ZFaceCorner(void *fs, Vec ZFace, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpZFaceCorner((FDSTAG*)fs, ZFace, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_CenterCorner(void *fs, Vec Center, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpCenterCorner((FDSTAG*)fs, Center, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_XYEdgeCorner(void *fs, Vec XYEdge, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpXYEdgeCorner((FDSTAG*)fs, XYEdge, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_XZEdgeCorner(void *fs, Vec XZEdge, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpXZEdgeCorner((FDSTAG*)fs, XZEdge, Corner, cpp_flags);
}

PetscErrorCode LaMEMInterpolate_YZEdgeCorner(void *fs, Vec YZEdge, Vec Corner, LaMEMInterpFlags iflag) {
    InterpFlags cpp_flags = convertInterpFlags(iflag);
    return InterpYZEdgeCorner((FDSTAG*)fs, YZEdge, Corner, cpp_flags);
}

}