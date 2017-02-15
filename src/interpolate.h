/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   interpolate.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//.............   FDSTAG OUTPUT VECTOR INTERPOLATION ROUTINES   .............
//---------------------------------------------------------------------------
#ifndef __interpolate_h__
#define __interpolate_h__

//---------------------------------------------------------------------------

// Interpolation flags

struct InterpFlags
{
	PetscBool update;    // update vs. overwrite target vector
	PetscBool use_bound; // use boundary ghost points for interpolation

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
// Boundary ghost points are assumed to exist for:
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
#endif
