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
 **    filename:   LaMEM.h
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
//.........................   Main include file   ...........................
//---------------------------------------------------------------------------

#ifndef __LaMEM_h__
#define __LaMEM_h__

//-----------------------------------------------------------------------------
//   PREFERABLE VARIABLES
//
//   PetscInt    - for all indices                    (can be int or long long int)
//   PetscScalar - for all floating point variables   (can be float or double)
//   float       - for reduced size output
//   size_t      - for all sizes offsets & counters   (unsigned long long int)
//   PetscMPIInt - for passing service integer parameters to MPI functions (int)
//   MPIU_SCALAR - appropriate MPI Data Type for sending/receiving PetsScalar
//   MPIU_INT    - appropriate MPI Data Type for sending/receiving PetscInt
//
//-----------------------------------------------------------------------------
// space dimension
#define SPDIM 3

// cast macros
#define LLD long long int

// use this to enable asprintf
#ifndef _GNU_SOURCE
	#define _GNU_SOURCE
#endif

// Identify gcc compiler.
// Take care that other compilers may also define __GNUC__ macro.

#undef GCC_COMPILER

#if defined (__GNUC__) && !defined (__INTEL_COMPILER)

	#define GCC_COMPILER

#endif

//-----------------------------------------------------------------------------
// EXTERNAL INCLUDES
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <sys/stat.h>
#include <errno.h>

// Get rid of unnecessary PETSc-induced warnings when using gcc compiler
// This is supposed to work with -Wall -Wconversion -Wcast-qual -Wshadow options

#ifdef GCC_COMPILER

	#pragma GCC diagnostic ignored "-Wconversion"
	#pragma GCC diagnostic ignored "-Wsign-conversion"
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	#pragma GCC diagnostic ignored "-Wcast-qual"

	#include <petsc.h>

	#pragma GCC diagnostic warning "-Wconversion"
	#pragma GCC diagnostic warning "-Wsign-conversion"
	#pragma GCC diagnostic warning "-Wunused-parameter"
	#pragma GCC diagnostic warning "-Wcast-qual"

#else

	#include <petsc.h>

#endif

//-----------------------------------------------------------------------------
// UDEFINE COMPLEX UNIT MACRO
//-----------------------------------------------------------------------------
#ifdef I
	#undef I
#endif

//-----------------------------------------------------------------------------
// PROTOTYPES
//-----------------------------------------------------------------------------

// LaMEM library main function

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode LaMEMLibMain(void *param);

#ifdef __cplusplus
}
#endif

#include "types.h"


//-----------------------------------------------------------------------------
#endif
