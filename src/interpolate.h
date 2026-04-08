/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//.............   FDSTAG OUTPUT VECTOR INTERPOLATION ROUTINES   .............
//---------------------------------------------------------------------------
#ifndef __interpolate_h__
#define __interpolate_h__

//---------------------------------------------------------------------------

struct FDSTAG;

//---------------------------------------------------------------------------
// Edge and corner initialization functions:

PetscErrorCode SetEdgeCornerCenter(FDSTAG *fs, Vec Center);

PetscErrorCode SetEdgeCornerXFace (FDSTAG *fs, Vec XFace);

PetscErrorCode SetEdgeCornerYFace (FDSTAG *fs, Vec YFace);

PetscErrorCode SetEdgeCornerZFace (FDSTAG *fs, Vec ZFace);

//---------------------------------------------------------------------------

// Interpolation flags

struct InterpFlags
{
 public:

	PetscInt update;    // update vs. overwrite target vector
	PetscInt use_bound; // use boundary ghost points for interpolation
};

//---------------------------------------------------------------------------
// Interpolation functions:
//
// x-face  -> corner   InterpXFaceCorner
// y-face  -> corner   InterpYFaceCorner
// z-face  -> corner   InterpZFaceCorner
// center  -> corner   InterpCenterCorner
// xy-edge -> corner   InterpXYEdgeCorner
// xz-edge -> corner   InterpXZEdgeCorner
// yz-edge -> corner   InterpYZEdgeCorner
//
// All functions perform distance-based interpolation.
// All functions assume input vectors in local format.
// Entire boundary & ghost points are supposed to be initialized beforehand.
// Boundary ghost points are assumed to be initialized for:
//    - face arrays: in tangential directions (i.e. y & z for x)
//    - central arrays: in all directions
// No boundary ghost points are assumed to exist for edge arrays.

//---------------------------------------------------------------------------

PetscErrorCode InterpXFaceCorner (FDSTAG *fs, Vec XFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode InterpYFaceCorner (FDSTAG *fs, Vec YFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode InterpZFaceCorner (FDSTAG *fs, Vec ZFace,  Vec Corner, InterpFlags iflag);

PetscErrorCode InterpCenterCorner(FDSTAG *fs, Vec Center, Vec Corner, InterpFlags iflag);

PetscErrorCode InterpXYEdgeCorner(FDSTAG *fs, Vec XYEdge, Vec Corner, InterpFlags iflag);

PetscErrorCode InterpXZEdgeCorner(FDSTAG *fs, Vec XZEdge, Vec Corner, InterpFlags iflag);

PetscErrorCode InterpYZEdgeCorner(FDSTAG *fs, Vec YZEdge, Vec Corner, InterpFlags iflag);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define SET_EDGE_CORNER(a, K, J, I, k, j, i, pmdof) \
	a[K][J][I] = a[k][j][I] + a[k][J][i] + a[K][j][i] - 2.0*pmdof;

//---------------------------------------------------------------------------
#endif
