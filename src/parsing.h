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
 **    filename:   parsing.h
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

/*
 *  Originally developed by Dave A. May
 *  Copyright 2011 Geophysical Fluid Dynamics. All rights reserved.
 *
 */

//---------------------------------------------------------------------------
//...................   Input file parsing routines   .......................
//---------------------------------------------------------------------------
#ifndef __parsing_h__
#define __parsing_h__
//-----------------------------------------------------------------------------

typedef enum
{
	_REQUIRED_,
	_OPTIONAL_

} ParamType;

//-----------------------------------------------------------------------------
// Input file buffer
//-----------------------------------------------------------------------------

typedef struct
{
	//=====================================================================
	//
	// input file buffer
	//
	// supports reading information from data blocks of similar type
	// e.g. material properties or softening law parameters
	// reading PETSc options is disabled in block mode
	//
	// corresponding access modes are:
	//    * flat  (parse entire file)
	//    * block (parse current data block, defined by line ranges)
	//
	// All strings reserve two null characters in the end to detect overrun
	//
	//=====================================================================

	PetscInt   nchar;   // number of characters
	PetscInt   nlines;  // number of lines
	char      *buff;    // file buffer
	char     **line;    // pointers to lines
	PetscInt   nblocks; // number of data blocks
	PetscInt   blockID; // ID of block to be retrieved
	PetscInt  *blBeg;   // starting lines of blocks
	PetscInt  *blEnd;   // ending lines of blocks

} FB;

//-----------------------------------------------------------------------------

PetscErrorCode FBLoad(FB **pfb);

PetscErrorCode FBDestroy(FB **pfb);

PetscErrorCode FBFindBlocks(FB *fb, ParamType ptype, const char *keybeg, const char *keyend);

PetscErrorCode FBFreeBlocks(FB *fb);

PetscErrorCode FBGetRanges(FB *fb, PetscInt *lnbeg, PetscInt *lnend);

PetscErrorCode FBGetIntArray(
		FB         *fb,
		const char *key,
		PetscInt   *nvalues,
		PetscInt   *values,
		PetscInt    num,
		PetscBool  *found);

PetscErrorCode FBGetScalarArray(
		FB          *fb,
		const char  *key,
		PetscInt    *nvalues,
		PetscScalar *values,
		PetscInt     num,
		PetscBool   *found);

PetscErrorCode FBGetString(
		FB         *fb,
		const char *key,
		char       *str,    // output string
		size_t      len,    // full size of output string
		PetscBool  *found);

//-----------------------------------------------------------------------------
// Wrappers
//-----------------------------------------------------------------------------

PetscErrorCode getIntParam(
		FB         *fb,
		ParamType   ptype,
		const char *key,
		PetscInt   *val,
		PetscInt    num,
		PetscInt    maxval);

PetscErrorCode getScalarParam(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		PetscScalar *val,
		PetscInt     num,
		PetscScalar  scal);

// string is initialized with default value, if available, otherwise set to zero
PetscErrorCode getStringParam(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		char        *str,         // output string
		size_t       len,         // full size of output string
		const char  *_default_);  // default value (optional)

//-----------------------------------------------------------------------------
// PETSc options parsing functions
//-----------------------------------------------------------------------------

PetscErrorCode PetscOptionsReadFromFile(FB *fb);

PetscErrorCode PetscOptionsReadRestart(FILE *fp);

PetscErrorCode PetscOptionsWriteRestart(FILE *fp);

PetscErrorCode  PetscOptionsGetCheckString(
	const char   name[],
	char         string[],
	size_t       len,
	PetscBool   *set);

//-----------------------------------------------------------------------------
#endif
