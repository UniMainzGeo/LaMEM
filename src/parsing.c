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
 **    filename:   parsing.c
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

#include "LaMEM.h"
#include "parsing.h"

// define separators
const char comments[] = "#!|";
#define NUL '\0'

//-----------------------------------------------------------------------------
char *trim(char *str)
{
	// Courtesy of http://gd.tuwien.ac.at/languages/c/cref-mleslie/CONTRIB/SNIP/trim.c

	char *ibuf = str, *obuf = str;
	PetscInt i = 0, cnt = 0;


	//  Trap NULL

	if (str)
	{
		// Remove leading spaces (from RMLEAD.C)

		for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf)
				;
		if (str != ibuf)
				memmove(str, ibuf, (size_t)(ibuf - str));

		//  Collapse embedded spaces (from LV1WS.C)

		while (*ibuf)
		{
			if (isspace(*ibuf) && cnt)
					ibuf++;
			else
			{
				if (!isspace(*ibuf))
						cnt = 0;
				else
				{
					*ibuf = ' ';
					cnt = 1;
				}
				obuf[i++] = *ibuf++;
			}
		}
		obuf[i] = NUL;

		//  Remove trailing spaces (from RMTRAIL.C)

		while (--i >= 0)
		{
			if (!isspace(obuf[i]))
					break;
		}
		obuf[++i] = NUL;
	}
	return str;
}
//-----------------------------------------------------------------------------
PetscInt is_comment_line( const char line[] )
{
	PetscInt i,L;

	L = (PetscInt)strlen(comments);
	for( i=0; i<L; i++ ) {
		if( line[0] == comments[i] ) {
			return _TRUE;
		}
	}

	return _FALSE;
}
//-----------------------------------------------------------------------------
void strip( char line[] )
{
	// removes all characters before and including the = sign

	PetscInt i,line_L;

	line_L = (PetscInt)strlen(line);

	for( i=0; i<line_L; i++ ) {
		if( line[i] != '=' ) {
			line[i] = ' ';
		}
		else {
			line[i] = ' ';
			return;
		}
	}
}
//-----------------------------------------------------------------------------
void strip_L( char str[] )
{
	PetscInt count, i = 0;
	char *p = str;
	PetscInt L;

	while(*p) {
		if( (*p == ' ') || (*p == '\t') ) {
			i++;
			p++;
		}
		else { break; }
	}
	count = i;

	L = (PetscInt)strlen( str );
	memmove( &str[0], &str[count], sizeof(char)*(size_t)(L-count) );
	str[L-count] = 0;
}
//-----------------------------------------------------------------------------
void trim_past_comment( char line[] )
{
	// removes all characters after and including the first comment character encountered

	PetscInt i,j,k,L,line_L;

	line_L = (PetscInt)strlen(line);
	L = (PetscInt)strlen(comments);

	for( i=0; i<line_L; i++ ) {
		for( j=0; j<L; j++ ) {
			if( line[i] == comments[j] ) {
				// prune everything past
				for( k=i; k<line_L; k++ ) {
					line[k] = ' ';
				}
				// remove white space
				trim(line);
				return;
			}
		}
	}
}
//-----------------------------------------------------------------------------
void strip_all_whitespace( char str[], char str2[] )
{
	PetscInt i = 0;
	char *p = str;
	while(*p)
	{
		if( (*p != ' ') && (*p != '\t') )
				str2[i++] = *p;
		p++;
	}
	str2[i] = 0;
}
//-----------------------------------------------------------------------------
PetscInt key_matches( const char key[], char line[] )
{
	PetscInt key_L, line_L, i;
	PetscInt c, cmp;
	char LINE[10000];


	strip_all_whitespace(line, LINE);

	c = 0;
	line_L = (PetscInt)strlen( LINE );
	key_L = (PetscInt)strlen( key );

	for( i=0; i<line_L; i++ ) {
		if( LINE[i] == '=' ) {
			c = i;
			break;
		}
	}
	// blank line or a comment line
	if( c== 0 ) {
		return _FALSE;
	}

	// The line length should not be smaller than the key
	if (c<key_L){
		return _FALSE;
	}

	cmp = strncmp(key,LINE, (size_t)(c) );

	if( cmp == 0 ) {
		return _TRUE;
	}
	else {
		return _FALSE;
	}
}
//-----------------------------------------------------------------------------
PetscInt material_key_matches( const char key[], char line[] )
{
	PetscInt line_L, i;
	PetscInt c, cmp;
	char LINE[10000];

	strip_all_whitespace(line, LINE);

	c = 0;
	line_L = (PetscInt)strlen( LINE );

	for( i=0; i<line_L; i++ ) {
		if( LINE[i] == '>' ) {
			c = i;
			break;
		}
	}
	// blank line or a comment line
	if( c== 0 ) {
		return _FALSE;
	}

	cmp = strncmp(key,LINE, (size_t)(c) );

	if( cmp == 0 ) {
		return _TRUE;
	}
	else {
		return _FALSE;
	}
}
//-----------------------------------------------------------------------------
// search file look for key in the file
void parse_GetInt( FILE *fp, const char key[], PetscInt *value, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match, int_val;

	// init
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		int_val = (PetscInt)strtol( line, NULL, 0 );

		*value = int_val;
		*found = _TRUE;
		return;
	}
}
//-----------------------------------------------------------------------------
void parse_GetIntAllInstances( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	PetscInt int_val;
	PetscInt count;

	// init
	*found = _FALSE;
	count = 0;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		int_val = (PetscInt)strtol( line, NULL, 0 );

		if( max_L == count ) {
			printf("parse_GetIntAllInstances: Error, input array is not large enough to hold result \n");
			return;
		}
		values[count] = int_val;


		count++;

		*found = _TRUE;
	}

	*nvalues = count;
}
//-----------------------------------------------------------------------------
void parse_GetDouble( FILE *fp, const char key[], double *value, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;

	// init
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		double_val = strtod( line, NULL );

		*value = double_val;
		*found = _TRUE;
		return;
	}
}
//-----------------------------------------------------------------------------
void parse_GetDoubleAllInstances( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;
	PetscInt count;

	// init
	*found = _FALSE;
	count = 0;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		double_val = strtod( line, NULL );

		if( max_L == count ) {
			printf("parse_GetDoubleAllInstances: Error, input array is not large enough to hold result \n");
			return;
		}
		values[count] = double_val;

		count++;

		*found = _TRUE;
	}

	if( *found == _TRUE ) {
		PetscInt i;

		printf("Key(%s) -- Values( %e", key, values[0] );
		for( i=1; i<count; i++ ) {
			printf(", %e", values[i] );
		}printf(" )\n");
	}

	*nvalues = count;
}
//-----------------------------------------------------------------------------
void parse_GetDoubleArray( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;
	char *_line;
	PetscInt count;

	// init
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		count = 0;
		_line = line;
		for(;;) {
			char * endp;

			double_val = strtod(_line, &endp);
			values[count] = double_val;

			if( endp == _line ) {
				break;
			}

			_line = endp;

			count++;
		}

		*nvalues = count;
		*found = _TRUE;
		return;
	}
}
//-----------------------------------------------------------------------------
void parse_GetIntArray( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	PetscInt int_val;
	char *_line;
	PetscInt count;

	// init
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);

		count = 0;
		_line = line;
		for(;;) {
			char * endp;

			int_val = (PetscInt)strtol(_line, &endp,0);
			values[count] = int_val;

			if( endp == _line ) {
				break;
			}

			_line = endp;

			count++;
		}

		*nvalues = count;
		*found = _TRUE;
		return;
	}
}
//-----------------------------------------------------------------------------
void parse_GetString( FILE *fp, const char key[], char value[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match, line_L;
	char LINE[MAX_LINE_LEN];

	// init
	*found = _FALSE;

	memset( value, 0, sizeof(char)*(size_t)max_L );

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// strip word and equal sign
		strip(line);
		strip_all_whitespace(line, LINE);

		trim_past_comment(LINE);
		line_L = (PetscInt)strlen( LINE );

		// one byte for terminating null character should be reserved
		if( line_L > max_L-1 ) {
			printf("parse_GetString: Error, input string is not large enough to hold result \n");
			return;
		}

		strncpy( value, LINE, (size_t)line_L );

		*found = _TRUE;
		return;
	}
}
//-----------------------------------------------------------------------------
PetscErrorCode PetscOptionsReadFromFile(
		const char filename[],
		PetscInt *found, PetscInt PrintOut )
{
	PetscErrorCode ierr;
	char key[] = "<PetscOptionsStart>";
	char key_end[] = "<PetscOptionsEnd>";
	char line[MAX_LINE_LEN], next[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match, end_match;
	FILE *fp;
	char *all_options;

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED  
	ierr = PetscOptionsGetAll(NULL,  &all_options );  CHKERRQ(ierr); /* copy all command line args */

	PetscPrintf(PETSC_COMM_WORLD," Parsing input file : %s \n", filename);

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		printf("LaMEM_parse_input: Cannot open input file %s \n", filename );
		exit(1);
	}

	// init
	*found = _FALSE;

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// if line is blank
		if( strlen(line) == 0 ) {	continue;	}

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = material_key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		// scan until the end-match symbol is found
		end_match = material_key_matches( key_end, line );

		while( end_match == _FALSE ) {
			fgets( next, MAX_LINE_LEN-1, fp );

			// get rid of white space
			trim(next);

			// if line is blank
			if( strlen(next) == 0 ) {	continue;	}

			// is first character a comment ?
			comment = is_comment_line( next );
			if( comment == _TRUE ) {   continue;  }

			end_match = material_key_matches( key_end, next );
			if( end_match == _FALSE ) {
				trim_past_comment(next);

				if (PrintOut==1){
					PetscPrintf(PETSC_COMM_WORLD,"# Adding PetscOption: %s \n", next );
				}

				// Add to petsc options
				ierr = PetscOptionsInsertString(NULL, next); CHKERRQ(ierr);

				continue;
			}

			*found = _TRUE;
		}
	}

	fclose(fp);

	// force command line args in to override defaults
	ierr = PetscOptionsInsertString(NULL, all_options); CHKERRQ(ierr);

	ierr = PetscFree(all_options);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
