//---------------------------------------------------------------------------
//.............   FDSTAG OUTPUT VECTOR INTERPOLATION ROUTINES   .............
//---------------------------------------------------------------------------
#ifndef __interpolate_h__
#define __interpolate_h__

//---------------------------------------------------------------------------

// Interpolation flags

typedef struct
{
	PetscBool update;    // update vs. overwrite target vector
	PetscBool use_bound; // use boundary ghost points for interpolation

} InterpFlags;

//---------------------------------------------------------------------------
// Interpolation functions:
//
// x-face  -> corner   FDSTAGInterpXFaceCorner
// y-face  -> corner   FDSTAGInterpYFaceCorner
// z-face  -> corner   FDSTAGInterpZFaceCorner
// center  -> corner   FDSTAGInterpCenterCorner
// xy-edge -> corner   FDSTAGInterpXYEdgeCorner
// xz-edge -> corner   FDSTAGInterpXZEdgeCorner
// yz-edge -> corner   FDSTAGInterpYZEdgeCorner
//
// All functions perform distance-based interpolation.
// All functions assume input vectors in local format.
// Entire boundary & ghost points are supposed to be initialized beforehand.
// Boundary ghost points are assumed to exist for:
//    - face arrays: in tangential directions (i.e. y & z for x)
//    - central arrays: in all directions
// No boundary ghost points are assumed to exist for edge arrays.

//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInterpXFaceCorner (FDSTAG *fs, Vec XFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpYFaceCorner (FDSTAG *fs, Vec YFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpZFaceCorner (FDSTAG *fs, Vec ZFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpCenterCorner(FDSTAG *fs, Vec Center, Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpXYEdgeCorner(FDSTAG *fs, Vec XYEdge, Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpXZEdgeCorner(FDSTAG *fs, Vec XZEdge, Vec Corner, InterpFlags iflag);

PetscErrorCode FDSTAGInterpYZEdgeCorner(FDSTAG *fs, Vec YZEdge, Vec Corner, InterpFlags iflag);

//---------------------------------------------------------------------------
#endif
