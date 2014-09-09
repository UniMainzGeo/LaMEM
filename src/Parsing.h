
#ifndef __Parsing_h__
#define __Parsing_h__


#define MAX_LINE_LEN 5000
#define _TRUE 1
#define _FALSE 0
//#define _LOG_PARSER




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
void MaterialTypeParser_FindDouble( const char key[], FILE* fp, const long start, const long end, double default_val, double *val, PetscInt *found );
void MaterialTypeParser_FindInt( const char key[],FILE* fp, const long start, const long end, PetscInt default_val, PetscInt *val, PetscInt *found );

PetscErrorCode PetscOptionsReadFromFile(const char filename[],PetscInt *found, PetscInt PrintOut );

void strip_L( char str[] );
PetscInt material_key_matches( const char key[], char line[] );



#endif
