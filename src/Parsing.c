/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Parsing.h"

// forward slash was removed to enable reading paths from string parameters
// const char comments[] = "#!/|";
   const char comments[] = "#!|";


/*
Courtesy of
http://gd.tuwien.ac.at/languages/c/cref-mleslie/CONTRIB/SNIP/trim.c
*/
#define NUL '\0'
char *trim(char *str)
{
	char *ibuf = str, *obuf = str;
	PetscInt i = 0, cnt = 0;

	/*
	**  Trap NULL
	*/

	if (str)
	{
		/*
		**  Remove leading spaces (from RMLEAD.C)
		*/

		for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf)
				;
		if (str != ibuf)
				memmove(str, ibuf, (size_t)(ibuf - str));

		/*
		**  Collapse embedded spaces (from LV1WS.C)
		*/

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

		/*
		**  Remove trailing spaces (from RMTRAIL.C)
		*/

		while (--i >= 0)
		{
			if (!isspace(obuf[i]))
					break;
		}
		obuf[++i] = NUL;
	}
	return str;
}

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

/* removes all characters before and including the = sign */
void strip( char line[] )
{
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


/* removes all characters after and including the first comment character encountered */
void trim_past_comment( char line[] )
{
	PetscInt i,j,k,L,line_L;

	line_L = (PetscInt)strlen(line);
	L = (PetscInt)strlen(comments);

	for( i=0; i<line_L; i++ ) {
		for( j=0; j<L; j++ ) {
			if( line[i] == comments[j] ) {
				/* prune everything past */
				for( k=i; k<line_L; k++ ) {
					line[k] = ' ';
				}
				/* remove white space */
				trim(line);
				return;
			}
		}
	}

}

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
	/* blank line or a comment line */
	if( c== 0 ) {
		return _FALSE;
	}

	/* The linelength should not be smaller than the key */
	if (c<key_L){
		return _FALSE;
	}

	cmp = strncmp(key,LINE, (size_t)(c) );

	//cmp = strncmp(key,LINE,key_L );

	//printf("%s : %s --  %lld \n", key, LINE, cmp );

	if( cmp == 0 ) {
	//	printf("NAME = %s   KEY=%s, keylength=%lld, c=%lld \n", LINE, key, key_L,c );
		return _TRUE;
	}
	else {
		return _FALSE;
	}
}


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
	/* blank line or a comment line */
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





/* search file look for key in the file */
void parse_GetInt( FILE *fp, const char key[], PetscInt *value, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
//	PetscBool flg;
	PetscInt match, int_val;

	/* init */
	*found = _FALSE;

	/* reset to start of file */
	//fseek( fp, 0, SEEK_SET);
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }


//		printf("Searching for %s \n", key );
//		printf("  Found %s \n", line );

		/* strip word and equal sign */
		strip(line);
//		printf("  stripped %s \n", line );


		int_val = (PetscInt)strtol( line, NULL, 0 );
//		printf("    %lld \n", int_val );
#ifdef _LOG_PARSER
		printf("Key(%s) -- Value(%lld) \n", key, (LLD)int_val );
#endif
		*value = int_val;
		*found = _TRUE;
		return;
	}

//	Check for command line override
//	asprintf( &option, "-%s", key );
//	PetscOptionsGetInt( PETSC_NULL, option, &value, &flg );
//	free( key );
}

void parse_GetIntAllInstances( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	PetscInt int_val;
	PetscInt count;

	/* init */
	*found = _FALSE;
	count = 0;

	/* reset to start of file */
	rewind( fp );


	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }


		/* strip word and equal sign */
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

#ifdef _LOG_PARSER
	if( *found == _TRUE ) {
		PetscInt i;

		printf("Key(%s) -- Values( %lld", key, (LLD)values[0] );
		for( i=1; i<count; i++ ) {
			printf(", %lld", (LLD)values[i] );
		}printf(" )\n");
	}
#endif

	*nvalues = count;

}



void parse_GetDouble( FILE *fp, const char key[], double *value, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;

	/* init */
	*found = _FALSE;

	/* reset to start of file */
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }
#ifdef _LOG_PARSER
			printf("Line(%s) key=%s \n", line,key );
#endif

		/* strip word and equal sign */
		strip(line);

		double_val = strtod( line, NULL );
#ifdef _LOG_PARSER

		printf("Key(%s) -- Value(%e) \n", key, double_val );
#endif
		*value = double_val;
		*found = _TRUE;
		return;
	}

}

void parse_GetDoubleAllInstances( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;
	PetscInt count;

	/* init */
	*found = _FALSE;
	count = 0;

	/* reset to start of file */
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }


		/* strip word and equal sign */
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



void parse_GetDoubleArray( FILE *fp, const char key[], PetscInt *nvalues, double values[], PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	double double_val;
	char *_line;
	PetscInt count;


	/* init */
	*found = _FALSE;

	/* reset to start of file */
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		/* strip word and equal sign */
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
#ifdef _LOG_PARSER
		printf("Key(%s) -- Values( %e", key, values[0] );
		PetscInt i;
		for( i=1; i<count; i++ ) {
			printf(", %e", values[i] );
		}printf(" )\n");
#endif

		*nvalues = count;
		*found = _TRUE;
		return;
	}

}

void parse_GetIntArray( FILE *fp, const char key[], PetscInt *nvalues, PetscInt values[], PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match;
	PetscInt int_val;
	char *_line;
	PetscInt count;


	/* init */
	*found = _FALSE;

	/* reset to start of file */
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		/* strip word and equal sign */
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
#ifdef _LOG_PARSER
		printf("Key(%s) -- Values( %e", key, values[0] );
		PetscInt i;
		for( i=1; i<count; i++ ) {
			printf(", %e", values[i] );
		}printf(" )\n");
#endif

		*nvalues = count;
		*found = _TRUE;
		return;
	}

}

void parse_GetString( FILE *fp, const char key[], char value[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match, line_L;
	char LINE[MAX_LINE_LEN];

	/* init */
	*found = _FALSE;

	memset( value, 0, sizeof(char)*(size_t)max_L );

	/* reset to start of file */
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		/* strip word and equal sign */
		strip(line);
		strip_all_whitespace(line, LINE);

		trim_past_comment(LINE);
		line_L = (PetscInt)strlen( LINE );
		if( line_L > max_L ) {
			printf("parse_GetString: Error, input string is not large enough to hold result \n");
			return;
		}

		strncpy( value, LINE, (size_t)line_L );
#ifdef _LOG_PARSER
		printf("Key(%s) -- Value(%s) \n", key, value );
#endif
		*found = _TRUE;
		return;
	}
}


void MaterialTypeParser_FindDouble( const char key[],
		FILE* fp, const long start, const long end, double default_val, double *val, PetscInt *found )
{

	*val = default_val;

	Material_helper_FindDouble( key, fp, start, end, val, found );

#ifdef _LOG_PARSER
	if( essential == _TRUE ) {
		if(!*found) {
			printf("\t  -- Cound not locate essential parameter '%s' \n",key);
		}
		else {
			printf("\t  ++ Found essential parameter '%s'=%e \n", key, *val );
		}
	}
	else {
		if(!*found) {
			printf("\t  -- Cound not locate non-essential parameter. Using default '%s'=%e \n", key, *val );
		}
		else {
			printf("\t  ++ Found non-essential parameter '%s'=%e \n", key, *val );
		}
	}
#endif

}
void MaterialTypeParser_FindInt( const char key[],
		FILE* fp, const long start, const long end, PetscInt default_val, PetscInt *val, PetscInt *found )
{

	*val = default_val;

	Material_helper_FindInt( key, fp, start, end, val, found );
#ifdef _LOG_PARSER
	if( essential == _TRUE ) {
		if(!*found) {
			printf("\t  -- Cound not locate essential parameter '%s' \n",key);
		}
		else {
			printf("\t  ++ Found essential parameter '%s'=%lld \n", key, (LLD)(*val) );
		}
	}
	else {
		if(!*found) {
			printf("\t  -- Cound not locate non-essential parameter. Using default '%s'=%lld \n", key, (LLD)(*val) );
		}
		else {
			printf("\t  ++ Found non-essential parameter '%s'=%lld \n", key, (LLD)(*val) );
		}
	}
#endif
}









void MaterialReadFromFile(
		Material M, const char filename[],
		void (*fp_parse_material_type)(Material M, const char*,const char*,FILE*,const long,const long,void**,PetscInt*),
		PetscInt *found )
{
	char key[] = "<MaterialStart>";
	char key_end[] = "<MaterialEnd>";
	char line[MAX_LINE_LEN], next[MAX_LINE_LEN], NAME[MAX_LINE_LEN];
	PetscInt comment;
	PetscInt match, end_match;
	PetscInt matfield_match;
	long start_mat, end_mat;
	PetscInt type_found;
	FILE *fp;
	char *phase_name, *attr_name, *type_name;
	void *data;
	PetscInt phase_id_val;

	end_mat = 0;

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		printf("LaMEM_parse_input: Cannot open input file %s \n", filename );
		exit(1);
	}



	/* init */
	*found = _FALSE;

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = material_key_matches( key, line );
		if( match == _FALSE ) {   continue;   }


		/* mark file position - beginning of this material */
		start_mat = ftell( fp );

		phase_name   = NULL;
		phase_id_val = -1;
		attr_name    = NULL;
		type_name    = NULL;

		/* find end of this material in file and find viscosity law */
		end_match = material_key_matches( key_end, line );
		while( end_match == _FALSE ) {
			fgets( next, MAX_LINE_LEN-1, fp );

			/* get rid of white space */
			trim(next);

			/* is first character a comment ? */
			comment = is_comment_line( next );
			if( comment == _TRUE ) {   continue;  }


		//	printf("CONTENTS LINE: %s \n", next );

			/* find phase, attribute, type */
			matfield_match = key_matches( "phase", next );
			if( matfield_match == _TRUE ) {
				strip(next);
				trim_past_comment(next);
				strip_all_whitespace(next,NAME);
				asprintf( &phase_name, "%s", NAME );
				continue;
			}

			matfield_match = key_matches( "phase_id", next );
			if( matfield_match == _TRUE ) {
				strip(next);
				phase_id_val = (PetscInt)strtol( next, NULL, 0 );
				continue;
			}

			matfield_match = key_matches( "attribute", next );
			if( matfield_match == _TRUE ) {
				strip(next);
				trim_past_comment(next);
				strip_all_whitespace(next,NAME);
				asprintf( &attr_name, "%s", NAME );
				continue;
			}

			matfield_match = key_matches( "type", next );
			if( matfield_match == _TRUE ) {
				strip(next);
				trim_past_comment(next);
				strip_all_whitespace(next,NAME);
				asprintf( &type_name, "%s", NAME );
				continue;
			}


			end_match = material_key_matches( key_end, next );
			if( end_match == _FALSE ) continue;
//			else { end_mat = ftell( fp ); }
			end_mat = ftell( fp );

		//	printf("LINE: %s \n", next );
			*found = _TRUE;
		}


		/*
		If the phase name is omitted but the phase_id is provided, create a phase name.
		*/
		if( (phase_name==NULL) && (phase_id_val!=-1) ) {
			asprintf( &phase_name, "phase_%lld", (LLD)phase_id_val );
		}


		/* jump to start of this material */
		//printf("*** Found { %s, %s, %s } ***\n", phase_name, attr_name, type_name );

		if(phase_name==NULL) {
			printf("Error(MaterialReadFromFile): Could not determine value of 'phase' \n");
			continue;
		}
		if( phase_id_val == -1 ) {
			printf("Error(MaterialReadFromFile): Could not determine value of 'phase_id' \n");
			continue;
		}
		if(attr_name==NULL)  {
			printf("Error(MaterialReadFromFile): Could not determine value of 'attribute' \n");
			continue;
		}
		if( type_name==NULL) {
			printf("Error(MaterialReadFromFile): Could not determine value of 'type' \n");
			continue;
		}

		/* look for data in file */
		fp_parse_material_type(M, attr_name, type_name, fp, start_mat, end_mat, &data, &type_found );
		if(type_found==_TRUE)  {
			MaterialAddPhase( M, phase_name, phase_id_val);
			MaterialAddAttribute( M, phase_name, attr_name);

			MaterialAddAttributeType( M, phase_name, attr_name, type_name, data);
		}


		if(phase_name!=NULL){ free(phase_name); }
		if(attr_name!=NULL){ free(attr_name); }
		if(type_name!=NULL){ free(type_name); }

	}


	fclose(fp);
}

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
	ierr = PetscOptionsGetAll( &all_options );  CHKERRQ(ierr); /* copy all command line args */

	PetscPrintf(PETSC_COMM_WORLD," Parsing input file : %s \n", filename);

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		printf("LaMEM_parse_input: Cannot open input file %s \n", filename );
		exit(1);
	}

	/* init */
	*found = _FALSE;

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* get rid of white space */
		trim(line);

		/* if line is blank */
		if( strlen(line) == 0 ) {	continue;	}

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = material_key_matches( key, line );
		if( match == _FALSE ) {   continue;   }

		/* scan until the end-match symbol is found */
		end_match = material_key_matches( key_end, line );

		while( end_match == _FALSE ) {
			fgets( next, MAX_LINE_LEN-1, fp );

			/* get rid of white space */
			trim(next);

			/* if line is blank */
			if( strlen(next) == 0 ) {	continue;	}

			/* is first character a comment ? */
			comment = is_comment_line( next );
			if( comment == _TRUE ) {   continue;  }


			end_match = material_key_matches( key_end, next );
			if( end_match == _FALSE ) {
				trim_past_comment(next);

				if (PrintOut==1){
					PetscPrintf(PETSC_COMM_WORLD,"# Adding PetscOption: %s \n", next );
				}
				// Add to petsc options
				ierr = PetscOptionsInsertString(next); CHKERRQ(ierr);

				continue;
			}

			*found = _TRUE;
		}



	}


	fclose(fp);

	ierr = PetscOptionsInsertString( all_options ); CHKERRQ(ierr); /* force command line args in to override defaults */

	ierr = PetscFree(all_options);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

