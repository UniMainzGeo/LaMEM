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
 **    filename:   tools.h
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
// ........................... UTILITY FUNCTIONS ............................
//---------------------------------------------------------------------------
#ifndef __tools_h__
#define __tools_h__
//---------------------------------------------------------------------------

/* $Id: tools.h 5681 2015-02-20 20:57:42Z ltbaumann $ */

typedef enum
{
	_NOT_FOUND_ERROR_,
	_NOT_FOUND_EXIT_

} exitType;

//---------------------------------------------------------------------------
//  basic statistic functions
//---------------------------------------------------------------------------

PetscScalar getArthMean(PetscScalar *data, PetscInt n);

PetscScalar getVar(PetscScalar *data, PetscInt n);

PetscScalar getStdv(PetscScalar *data, PetscInt n);

//---------------------------------------------------------------------------
// read arrays from PETSC options database with error checking
//---------------------------------------------------------------------------

PetscErrorCode GetScalDataItemCheckScale(
	const char  ident[],
	const char  name[],
	exitType    extp,
	PetscInt    n,
	PetscScalar *a,
	PetscScalar amin,
	PetscScalar amax,
	PetscScalar scal);

PetscErrorCode GetIntDataItemCheck(
	const char  ident[],
	const char  name[],
	exitType    extp,
	PetscInt    n,
	PetscInt    *a,
	PetscInt    amin,
	PetscInt    amax);

PetscErrorCode DMDAGetProcessorRank(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_col);

PetscErrorCode makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n);

PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n);

PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n);

//---------------------------------------------------------------------------
// checks whether processor has a zero rank in the communicator
PetscInt ISRankZero(MPI_Comm comm);

// check whether communicator is parallel (has more than one rank)
PetscInt ISParallel(MPI_Comm comm);

PetscErrorCode LaMEMCreateOutputDirectory(const char *DirectoryName);

//---------------------------------------------------------------------------

void polygon_box(
	PetscInt    *pnv,    // number of polygon vertices (can be modified)
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar  rtol,   // relative tolerance
	PetscScalar *atol,   // absolute tolerance
	PetscScalar *box);   // bounding box of a polygon

void in_polygon(
	PetscInt     np,     // number of test points
	PetscScalar *pcoord, // coordinates of test points
	PetscInt     nv,     // number of polygon vertices
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar *box,    // bounding box of a polygon (optimization)
	PetscScalar  atol,   // absolute tolerance
	PetscInt    *in);    // point location flags (1-inside, 0-outside)

//---------------------------------------------------------------------------

static inline void RotDispPoint2D(PetscScalar Xa[], PetscScalar Xb[], PetscScalar costh, PetscScalar sinth, PetscScalar xa[], PetscScalar xb[])
{
	PetscScalar r[2];

	// get radius vector
	r[0] = xa[0] - Xa[0];
	r[1] = xa[1] - Xa[1];

	// rotate & translate
	xb[0] = costh*r[0] - sinth*r[1] + Xb[0];
	xb[1] = sinth*r[0] + costh*r[1] + Xb[1];
}

//---------------------------------------------------------------------------
#endif
