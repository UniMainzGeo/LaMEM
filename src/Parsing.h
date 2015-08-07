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
 **    filename:   nlsolve.c
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
//...................   Input file parsing routines   .......................
//---------------------------------------------------------------------------
#ifndef __Parsing_h__
#define __Parsing_h__
//-----------------------------------------------------------------------------

#define MAX_LINE_LEN 1024
#define _TRUE 1
#define _FALSE 0

//-----------------------------------------------------------------------------

void parse_GetInt( FILE *fp, const char key[], PetscInt *value, PetscInt *found );

void parse_GetDouble( FILE *fp, const char key[], double *value, PetscInt *found );

void parse_GetString( FILE *fp, const char key[], char value[], PetscInt max_L, PetscInt *found );

void parse_GetDoubleArray( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt *found );

void parse_GetIntArray( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt *found );

void parse_GetDoubleAllInstances( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt max_L, PetscInt *found );

void parse_GetIntAllInstances( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt max_L, PetscInt *found );

PetscInt is_comment_line( const char line[] );

void strip( char line[] );

char *trim(char *str);

PetscInt key_matches( const char key[], char line[] );

void strip_all_whitespace( char str[], char str2[] );

void trim_past_comment( char line[] );

void strip_L( char str[] );

PetscInt material_key_matches( const char key[], char line[] );

PetscErrorCode PetscOptionsReadFromFile(const char filename[],PetscInt *found, PetscInt PrintOut );

//-----------------------------------------------------------------------------
#endif
