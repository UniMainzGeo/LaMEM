/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This sofware was developed at:
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

//---------------------------------------------------------------------------
#endif
