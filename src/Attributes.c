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

void MaterialCreate( Material *m )
{
	Material material;

	material = (Material)malloc( sizeof(struct _p_Material) );

	material->n_phases = 0;
	material->phase = NULL;

	material->n_attributes_MAX = 0;
	material->n_types_MAX = 0;
	material->size_ptr = sizeof(PetscInt*);
	material->material_phase_table = NULL;

	char_ptr_stack_create(&material->a);

	*m = material;
}

void MaterialDestroy( Material M )
{
	PetscInt p, a, t;


	for( p=0; p<M->n_phases; p++ ) {
		Phase phase = M->phase[p];

		free(phase->name); // get rid of memory leak

		for( a=0; a<phase->n_attributes; a++ ) {
			Attribute attr = phase->attribute[a];

			free(attr->name); // get rid of memory leak

			for( t=0; t<attr->n_types; t++ ) {
				AttributeType type = attr->type[t];

				free(type->name); // get rid of memory leak

				if( type->data != NULL ) {
					free( type->data );
				}
				free( type );
			}
			free( attr->type);

			free( attr );
		}
		free( phase->attribute );
		free( phase );
	}
	free( M->phase );


	if( M->material_phase_table != NULL ) {
		void ****map = M->material_phase_table;

		for( p=0; p<M->n_phases; p++ ) {
			for( a=0; a<M->n_attributes_MAX; a++ ) {
				free( map[p][a] );
			}
			free( map[p] );
		}
		free( map );
	}


	// squeeze my shit out
	char_ptr_stack_delete(&M->a);

	free( M );

}




void _PhaseInit( Phase p )
{
	p->name         = NULL;
	p->phase_number = 0;
	p->P_id         = 0;
	p->n_attributes = 0;
	p->attribute    = NULL;
}
void _AttributeInit( Attribute a )
{
	a->name         = NULL;
	a->A_id           = 0;
	a->n_types = 0;
	a->type    = NULL;
}
void _AttributeTypeInit( AttributeType t )
{
	t->name         = NULL;
	t->T_id           = 0;
	t->data         = NULL;
}


void MaterialAddPhase( Material m, const char phase_name[], const PetscInt phase_number)
{
	PetscInt p;
	Phase this_phase;


	/* check phase name doesn't exist */
	for( p=0; p<m->n_phases; p++ ) {
		if( strcmp(m->phase[p]->name,phase_name) == 0 ) {
#ifdef _LOG_PARSER
			printf("MaterialAddPhase: Phase '%s' already exists. Ignoring... \n", phase_name );
#endif
			return;
		}
	}

	/* check phase number doesn't exist */
	for( p=0; p<m->n_phases; p++ ) {
		if( m->phase[p]->phase_number == phase_number ) {
#ifdef _LOG_PARSER
			printf("MaterialAddPhase: Phase '%s' already exists with the phase_number '%lld'. Ignoring... \n", m->phase[p]->name, (LLD)phase_number );
#endif
			return;
		}
	}


	/* check if first phase */
	if( m->phase == NULL ) {
		m->phase = (Phase*)malloc( sizeof(Phase) );
	}
	else {
		Phase *tmp_p;
		tmp_p = (Phase*)realloc( m->phase, sizeof(Phase)*(size_t)(m->n_phases+1) );
		m->phase = tmp_p;
	}

	this_phase = (Phase)malloc( sizeof(struct _p_Phase) );
	_PhaseInit(this_phase);

	asprintf( &this_phase->name, "%s", phase_name );
	this_phase->phase_number = phase_number;
	this_phase->P_id = m->n_phases;

	m->phase[ m->n_phases ] = this_phase;

	m->n_phases++;
}

void _MaterialGetPhase( Material m, const char phase_name[], Phase *phase )
{
	PetscInt p;

	*phase = NULL;
	for( p=0; p<m->n_phases; p++ ) {
		if( strcmp(m->phase[p]->name,phase_name) == 0 ) {
			*phase = m->phase[p];
			return;
		}
	}
}

void _PhaseGetAttribute( Phase phase, const char attr_name[], Attribute *attr )
{
	PetscInt a;

	*attr = NULL;
	for( a=0; a<phase->n_attributes; a++ ) {
		if( strcmp(phase->attribute[a]->name,attr_name) == 0 ) {
			*attr = phase->attribute[a];
			return;
		}
	}
}

void _AttributeGetType( Attribute attr, const char type_name[], AttributeType *type )
{
	PetscInt t;

	*type = NULL;
	for( t=0; t<attr->n_types; t++ ) {
		if( strcmp(attr->type[t]->name,type_name) == 0 ) {
			*type = attr->type[t];
			return;
		}
	}
}



void MaterialAddAttribute( Material m, const char phase_name[], const char attr_name[])
{
	Phase phase;
	Attribute attr, this_attr;

	/* check phase exists */
	_MaterialGetPhase( m, phase_name, &phase );
	if( phase == NULL ) {
#ifdef _LOG_PARSER
		printf("MaterialAddAttribute: Phase '%s' does not exist... \n", phase_name );
#endif
		return;
	}


	/* check attribute doesn't exist */
	_PhaseGetAttribute( phase, attr_name, &attr );
	if( attr != NULL ) {
#ifdef _LOG_PARSER
		printf("MaterialAddPhase: Phase '%s' already exists. Ignoring... \n", phase_name );
#endif
		return;
	}

	/* check if first attribute */
	if( phase->attribute == NULL ) {
		phase->attribute = (Attribute*)malloc( sizeof(Attribute) );
	}
	else {
		Attribute *tmp_a;
		tmp_a = (Attribute*)realloc( phase->attribute, sizeof(Attribute)*(size_t)(phase->n_attributes+1) );
		phase->attribute = tmp_a;
	}

	this_attr = (Attribute)malloc( sizeof(struct _p_Attribute) );
	_AttributeInit(this_attr);

	asprintf( &this_attr->name, "%s", attr_name );
	this_attr->A_id = phase->n_attributes;

	phase->attribute[ phase->n_attributes ] = this_attr;

	phase->n_attributes++;
}

void MaterialAddAttributeType( Material m, const char phase_name[], const char attr_name[], const char type_name[], void *data)
{
	Phase phase;
	Attribute attr;
	AttributeType type, this_type;


	/* check phase exists */
	_MaterialGetPhase( m, phase_name, &phase );
	if( phase == NULL ) {
#ifdef _LOG_PARSER
		printf("MaterialAddAttributeType: Phase '%s' does not exist... \n", phase_name );
#endif
		return;
	}
	/* check attribute exists */
	_PhaseGetAttribute( phase, attr_name, &attr );
	if( attr == NULL ) {
#ifdef _LOG_PARSER
		printf("MaterialAddAttributeType: Attribute '%s' does not exist... \n", attr_name );
#endif
		return;
	}
	/* check type doesn't exist */
	_AttributeGetType( attr, type_name, &type );
	if( type != NULL ) {
#ifdef _LOG_PARSER
		printf("MaterialAddAttributeType: Type '%s' already exists. Ignoring... \n", type_name );
#endif
		return;
	}


	/* check if first type */
	if( attr->type == NULL ) {
		attr->type = (AttributeType*)malloc( sizeof(AttributeType) );
	}
	else {
		AttributeType *tmp_t;
		tmp_t = (AttributeType*)realloc( attr->type, sizeof(AttributeType)*(size_t)(attr->n_types+1) );
		attr->type = tmp_t;
	}

	this_type = (AttributeType)malloc( sizeof(struct _p_AttributeType) );
	_AttributeTypeInit(this_type);

	asprintf( &this_type->name, "%s", type_name );
	this_type->T_id = attr->n_types;
	this_type->data = data;

	attr->type[ attr->n_types ] = this_type;


	attr->n_types++;

}

void MaterialQueryGetIndices( Material m, const char phase_name[], const char attr_name[], const char type_name[], PetscInt *p_idx, PetscInt *a_idx, PetscInt *t_idx )
{
	PetscInt p,a,t;


	if( p_idx != NULL ) { *p_idx = -1; }
	if( a_idx != NULL ) { *a_idx = -1; }
	if( t_idx != NULL ) { *t_idx = -1; }

	for( p=0; p<m->n_phases; p++ ) {
		Phase phase = m->phase[p];

		if( strcmp(phase->name,phase_name)==0 ) {
			if( p_idx != NULL ) { *p_idx = p; }

			for( a=0; a<phase->n_attributes; a++ ) {
				Attribute attr = phase->attribute[a];

				if( strcmp(attr->name,attr_name)==0 ) {
					if( a_idx != NULL ) { *a_idx = a; }

					for( t=0; t<attr->n_types; t++ ) {
						AttributeType type = attr->type[t];

						if( strcmp(type->name,type_name)==0 ) {
							if( t_idx != NULL ) { *t_idx = t; }
						}
					}
				}
			}
		}
	}
}

void _MaterialCreateTable( Material m )
{
	PetscInt n_phases_MAX;
	PetscInt p,a,t;
	void ****map;


	if( m->material_phase_table != NULL ) { return; }

	/* NO ERROR CHECKING */
	/* check table exists, create if necessary */
	if( m->material_phase_table == NULL ) {
		/* get sizes */
		n_phases_MAX = m->n_phases;
		m->n_attributes_MAX = 0;
		m->n_types_MAX = 0;
		for( p=0; p<m->n_phases; p++ ) {
			Phase phase = m->phase[p];

			if( phase->n_attributes > m->n_attributes_MAX ) {
				m->n_attributes_MAX = phase->n_attributes;
			}

			for( a=0; a<phase->n_attributes; a++ ) {
				Attribute attr = phase->attribute[a];

				if( attr->n_types > m->n_types_MAX ) {
					m->n_types_MAX = attr->n_types;
				}
			}

		}


		/* allocate space */
		map = (void****)malloc( m->size_ptr * (size_t)n_phases_MAX );
		for( p=0; p<n_phases_MAX; p++ ) {
			map[p] = (void***)malloc( m->size_ptr * (size_t)m->n_attributes_MAX );
			for( a=0; a<m->n_attributes_MAX; a++ ) {
				map[p][a] = (void**)malloc( m->size_ptr * (size_t)m->n_types_MAX );
				for( t=0; t<m->n_types_MAX; t++ ) {
					/* init to NULL */
					map[p][a][t] = NULL;
				}
			}
		}
	}


	/* now fill it out */
	for( p=0; p<m->n_phases; p++ ) {
		Phase phase = m->phase[p];

		for( a=0; a<phase->n_attributes; a++ ) {
			Attribute attr = phase->attribute[a];

			for( t=0; t<attr->n_types; t++ ) {
				AttributeType type = attr->type[t];

				/* init to NULL */
				map[p][a][t] = (void*)type->data;
			}
		}
	}


	m->material_phase_table = (void****)map;

}

void MaterialGetPhaseAttributeType( Material m, const PetscInt phase_idx, const PetscInt attr_idx, const PetscInt type_idx, void **data )
{
	if( m->material_phase_table == NULL ) {
		_MaterialCreateTable(m);
	}
	*data = NULL;

	if( phase_idx >= m->n_phases ) { return; }
	if( attr_idx >= m->n_attributes_MAX ) { return; }
	if( type_idx >= m->n_types_MAX ) { return; }

	if( phase_idx < 0 ) { return; }
	if( attr_idx < 0 ) { return; }
	if( type_idx < 0 ) { return; }

	(*data) = m->material_phase_table[phase_idx][attr_idx][type_idx];
}

void MaterialGetAllPhases( Material M, PetscInt *n, Phase **phases )
{
	if(n!=NULL){ *n = M->n_phases; }
	if(phases!=NULL) { *phases = M->phase; }
}

void MaterialGetPhaseByName( Material m, const char name[], PetscInt *P_id, Phase *phase )
{
	PetscInt p;

	if(phase!=NULL) { *phase = NULL; }
	if(P_id!=NULL)  { *P_id = -1; }

	for( p=0; p<m->n_phases; p++ ) {
		if( strcmp(m->phase[p]->name,name) == 0 ) {

			if(phase!=NULL) { *phase= m->phase[p]; }
			if(P_id!=NULL)  { *P_id = m->phase[p]->P_id; }

			return;
		}
	}
}
void MaterialGetPhaseByID( Material m, PetscInt phase_id, PetscInt *P_id, Phase *phase )
{
	PetscInt p;

	if(phase!=NULL) { *phase = NULL; }
	if(P_id!=NULL)  { *P_id = -1; }

	for( p=0; p<m->n_phases; p++ ) {
		if( m->phase[p]->phase_number == phase_id ) {

			if(phase!=NULL) { *phase = m->phase[p]; }
			if(P_id!=NULL) { *P_id = m->phase[p]->P_id; }

			return;
		}
	}
}
void PhaseGetAllAttributes( Phase p, PetscInt *n, Attribute **attrs )
{
	if(n!=NULL){ *n = p->n_attributes; }
	if(attrs!=NULL) { *attrs = p->attribute; }
}
void PhaseGetAttributeByName( Phase p, const char name[], PetscInt *A_id, Attribute *attr )
{
	PetscInt a;

	if(attr!=NULL) { *attr = NULL; }
	if(A_id!=NULL)  { *A_id = -1; }

	for( a=0; a<p->n_attributes; a++ ) {
		if( strcmp(p->attribute[a]->name,name) == 0 ) {

			if(attr!=NULL) { *attr = p->attribute[a]; }
			if(A_id!=NULL)  { *A_id = p->attribute[a]->A_id; }

			return;
		}
	}
}
void AttributeGetAllTypes( Attribute a, PetscInt *n, AttributeType **types )
{
	if(n!=NULL){ *n = a->n_types; }
	if(types!=NULL) { *types = a->type; }
}
void AttributeGetTypeByName( Attribute a, const char name[], PetscInt *T_id, AttributeType *type )
{
	PetscInt t;

	if(type!=NULL) { *type = NULL; }
	if(T_id!=NULL)  { *T_id = -1; }

	for( t=0; t<a->n_types; t++ ) {
		if( strcmp(a->type[t]->name,name) == 0 ) {

			if(type!=NULL) { *type = a->type[t]; }
			if(T_id!=NULL)  { *T_id = a->type[t]->T_id; }

			return;
		}
	}
}
void AttributeTypeGetData( AttributeType t, void **data )
{
	*data = t->data;
}

void PhaseCompare( Phase p, const char name[], PetscInt *truth )
{
	*truth = _FALSE;

	if( strcmp(p->name,name) == 0 ) {
		*truth = _TRUE;
	}
}
void AttributeCompare( Attribute a, const char name[], PetscInt *truth )
{
	*truth = _FALSE;

	if( strcmp(a->name,name) == 0 ) {
		*truth = _TRUE;
	}
}
void AttributeTypeCompare( AttributeType t, const char name[], PetscInt *truth )
{
	*truth = _FALSE;

	if( strcmp(t->name,name) == 0 ) {
		*truth = _TRUE;
	}
}

void AttributeExist( Phase p, const char name[], PetscInt *truth )
{
	*truth = _FALSE;
	Attribute a,*attrs;
	PetscInt 	AA, n_attrs;

	PhaseGetAllAttributes( p, &n_attrs, &attrs );
	*truth = _FALSE;
	for( AA=0; AA<n_attrs; AA++ ) {
		a 			= 	attrs[AA];
		if( strcmp(a->name,name) == 0 ) {
			*truth = _TRUE;
		}
	}
}



void MaterialUpdateData( Material m,  PetscInt p, PetscInt a, PetscInt t, void *data)
{
	Phase 			phase = m->phase[p];
	Attribute 		attribute = phase->attribute[a];
	AttributeType 	type = attribute->type[t];

	type->data = data;

}




void MaterialView( Material m )
{
	PetscInt p,a,t;

	printf("Material:\n");
	for( p=0; p<m->n_phases; p++ ) {
		Phase phase = m->phase[p];

		printf("\tPhase[%lld]: %s, phase_id = %lld\n",(LLD)(phase->P_id), phase->name, (LLD)(phase->phase_number) );
		for( a=0; a<phase->n_attributes; a++ ) {
			Attribute attribute = phase->attribute[a];

			printf("\t\tAttribute[%lld]: %s\n",(LLD)(attribute->A_id), attribute->name );
			for( t=0; t<attribute->n_types; t++ ) {
				AttributeType type = attribute->type[t];

				printf("\t\t\tAttributeType[%lld]: %s\n",(LLD)(type->T_id), type->name );
				printf("\t\t\t  data: %p\n", type->data );
			}
		}

	}
}


void Material_helper_FindInt( const char key[], FILE* fp, const long start, const long end, PetscInt *val, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt match;
	long curr_pos;
	PetscInt comment;


	*found = _FALSE;

	fseek( fp, start, SEEK_SET );
	curr_pos = start;
	while( curr_pos != end ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _TRUE ) {
			strip(line);
			*val = (PetscInt)strtol( line, NULL, 0 );
			*found  = _TRUE;
			break;
		}

		curr_pos = ftell( fp );
	}
}
void Material_helper_FindDouble( const char key[], FILE* fp, const long start, const long end, double *val, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt match;
	long curr_pos;
	PetscInt comment;

	*found = _FALSE;

	fseek( fp, start, SEEK_SET );
	curr_pos = start;
	while( curr_pos != end ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( key, line );
		if( match == _TRUE ) {
			strip(line);
			*val = strtod( line, NULL );
			*found  = _TRUE;
			break;
		}

		curr_pos = ftell( fp );
	}
}

void Material_helper_CheckType( const char sought_type[], FILE* fp, const long start, const long end, PetscInt *found )
{
	char line[MAX_LINE_LEN], LINE[MAX_LINE_LEN];
	PetscInt match;
	long curr_pos;
	PetscInt comment;
	PetscInt sought_type_L, LINE_L;

	*found = _FALSE;
	sought_type_L = (PetscInt)strlen(sought_type);

	fseek( fp, start, SEEK_SET );
	curr_pos = start;
	while( curr_pos != end ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		/* is first character a comment ? */
		comment = is_comment_line( line );
		if( comment == _TRUE ) {   continue;  }

		match = key_matches( "type", line );
		if( match == _TRUE ) {
			/* strip word and equal sign */
			strip(line);
			strip_all_whitespace(line, LINE);
			trim_past_comment(LINE);
			LINE_L = (PetscInt)strlen(LINE);

			if( LINE_L != sought_type_L ) {
				continue;
			}

			if( strncmp( sought_type, LINE, (size_t)sought_type_L ) == 0 ) {
				*found  = _TRUE;
				break;
			}
		}

		curr_pos = ftell( fp );
	}
}

void test_map( void )
{
	Material M;
	double a, b, c;
	void *ptr;

	a = 1.1;
	b = 2.2;
	c = 3.3;

	MaterialCreate( &M );

	MaterialAddPhase( M, "mantle", 0);
	MaterialAddAttribute( M, "mantle", "VISCOUS");
	MaterialAddAttributeType( M, "mantle", "VISCOUS", "constant_a", (void*)&a);
	MaterialAddAttributeType( M, "mantle", "VISCOUS", "constant_b", (void*)&b);
	MaterialAddAttributeType( M, "mantle", "VISCOUS", "constant_c", (void*)&c);

	MaterialAddAttribute( M, "mantle", "ELASTIC");
	MaterialAddAttributeType( M, "mantle", "ELASTIC", "constant_b", (void*)&b);

	MaterialAddAttribute( M, "mantle", "PLASTIC");

	MaterialAddPhase( M, "crust", 2);
	MaterialAddAttribute( M, "crust", "PLASTIC");
	MaterialAddAttributeType( M, "crust", "PLASTIC", "constant_b", (void*)&b);
	MaterialAddAttributeType( M, "crust", "PLASTIC", "constant_c", (void*)&c);
	MaterialAddAttribute( M, "crust", "ELASTIC");

	MaterialAddPhase( M, "water", 1);
	MaterialAddAttribute( M, "water", "VISCOUS");

	MaterialAddAttribute( M, "crap", "VISCOUS");

	MaterialView(M);


	MaterialGetPhaseAttributeType( M, 0,0,0, &ptr );
	printf("ptr = %p \n", ptr );
	MaterialGetPhaseAttributeType( M, 0,0,1, &ptr );
	printf("ptr = %p \n", ptr );
	MaterialGetPhaseAttributeType( M, 0,0,2, &ptr );
	printf("ptr = %p \n", ptr );

	MaterialGetPhaseAttributeType( M, 0,1,0, &ptr );
	printf("ptr = %p \n", ptr );


	MaterialGetPhaseAttributeType( M, 1,0,0, &ptr );
	printf("ptr = %p \n", ptr );
	MaterialGetPhaseAttributeType( M, 1,0,1, &ptr );
	printf("ptr = %p \n", ptr );
	MaterialGetPhaseAttributeType( M, 1,1,0, &ptr );
	printf("ptr = %p \n", ptr );

	MaterialGetPhaseAttributeType( M, 3,0,0, &ptr );
	printf("ptr = %p \n", ptr );


}

/* material type helper functions */
void MaterialTypeCreate(
		Material M,
		const char attr_name[], const char type_name[],
		size_t size,
		void **type_data )
{
	void *data;
	struct _p_MaterialTypeBase *base;

	data = (void*)malloc( size );
	memset( data, 0, size );

	base = (struct _p_MaterialTypeBase*)data;

	asprintf( &base->attribute_name, "%s", attr_name );
	asprintf( &base->type_name, "%s", type_name );

	*type_data = data;

	// fill my shit in
	char_ptr_stack_insert( &M->a, base->attribute_name);
	char_ptr_stack_insert( &M->a, base->type_name);

}

void MaterialCompareAttributeTypeNames(
		const char attribute_a[], const char attribute_b[],
		const char type_a[], const char type_b[], PetscInt *match )
{
	*match = _FALSE;

	if( strcmp(attribute_a,attribute_b) == 0 ) {
		if( strcmp(type_a,type_b) == 0 ) {
			*match = _TRUE;
	}}
}


/*=================================================================================================================================*/
void ConstantViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "constant";
	const size_t  type_size             = sizeof(struct Material_ConstantElasticity);
	struct Material_ConstantViscosity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
//	ConstantViscosityCreate( &type_data );
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "eta0", fp, start_mat, end_mat, 1.0, &type_data->eta0, &value_found );

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetConstantViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *eta0 )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "constant";
	struct Material_ConstantViscosity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( eta0 != NULL ) {   *eta0 = data->eta0;   }

}
/*=================================================================================================================================*/



/*=================================================================================================================================*/
void PowerLawViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "power_law";
	const size_t  type_size             = sizeof(struct Material_PowerLawViscosity);
	struct Material_PowerLawViscosity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif


	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "eta0", fp, start_mat,  end_mat,  1.0, &type_data->eta0, &value_found );
	MaterialTypeParser_FindDouble( "N_exp", fp, start_mat, end_mat, 1.0, &type_data->N_exp, &value_found );
	MaterialTypeParser_FindDouble( "e0", fp, start_mat,    end_mat, 1.0, &type_data->e0,   &value_found );

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetPowerLawViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *eta0, double *N_exp, double *e0 )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "power_law";
	struct Material_PowerLawViscosity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( eta0  != NULL  ){   *eta0 = data->eta0;   }
	if( N_exp != NULL ){   *N_exp = data->N_exp; }
	if( e0    != NULL ){   *e0    = data->e0; }


}
/*=================================================================================================================================*/



/*=================================================================================================================================*/
void TempDepViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "temp_dep";
	const size_t  type_size          = sizeof(struct Material_TempDepViscosity);
	struct Material_TempDepViscosity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif


	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "PreExpFactor", fp, start_mat,  end_mat,  1.0, &type_data->PreExpFactor, &value_found );
	MaterialTypeParser_FindDouble( "N_exp", fp, start_mat, end_mat, 1.0, &type_data->N_exp, &value_found );
	MaterialTypeParser_FindDouble( "e0", fp, start_mat,    end_mat, 1.0, &type_data->e0,   &value_found );
	MaterialTypeParser_FindDouble( "ActivationEnergy", fp, start_mat,    end_mat, 1.0, &type_data->ActivationEnergy,   &value_found );

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetTempDepViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *PreExpFactor, double *N_exp, double *e0, double *ActivationEnergy )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "temp_dep";
	struct Material_TempDepViscosity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( PreExpFactor  != NULL  ){   *PreExpFactor = data->PreExpFactor;   }
	if( N_exp != NULL ){   *N_exp = data->N_exp; }
	if( e0    != NULL ){   *e0    = data->e0; }
	if( ActivationEnergy    != NULL ){   *ActivationEnergy    = data->ActivationEnergy; }


}
/*=================================================================================================================================*/



/*=================================================================================================================================*/
void TempDepNoPowerLawViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "tempdep_nopowerlaw";
	const size_t  type_size          = sizeof(struct Material_TempDepNoPowerLawViscosity);
	struct Material_TempDepNoPowerLawViscosity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif


	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "PreExpFactor", fp, start_mat,  end_mat,  1.0, &type_data->PreExpFactor, &value_found );
	MaterialTypeParser_FindDouble( "N_exp", fp, start_mat, end_mat, 1.0, &type_data->N_exp, &value_found );
	MaterialTypeParser_FindDouble( "e0", fp, start_mat,    end_mat, 1.0, &type_data->e0,   &value_found );
	MaterialTypeParser_FindDouble( "ActivationEnergy", fp, start_mat,    end_mat, 1.0, &type_data->ActivationEnergy,   &value_found );

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetTempDepNoPowerLawViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *PreExpFactor, double *N_exp, double *e0, double *ActivationEnergy )
{
	const char THIS_ATTRIBUTE_NAME[] = "VISCOSITY";
	const char THIS_TYPE_NAME[]      = "tempdep_nopowerlaw";
	struct Material_TempDepNoPowerLawViscosity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( PreExpFactor  != NULL  ){   *PreExpFactor = data->PreExpFactor;   }
	if( N_exp != NULL ){   *N_exp = data->N_exp; }
	if( e0    != NULL ){   *e0    = data->e0; }
	if( ActivationEnergy    != NULL ){   *ActivationEnergy    = data->ActivationEnergy; }


}
/*=================================================================================================================================*/



/*=================================================================================================================================*/
void ConstantElasticityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "ELASTICITY";
	const char    THIS_TYPE_NAME[]      = "constant";
	const size_t  type_size             = sizeof(struct Material_ConstantElasticity);
	/* auxillary params */
	struct Material_ConstantElasticity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "shear", fp, start_mat, end_mat, 1.0, &type_data->shear, &value_found );
	MaterialTypeParser_FindDouble( "bulk",  fp, start_mat, end_mat, 1.0, &type_data->bulk,  &value_found );

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetConstantElasticityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *shear, double *bulk )
{
	const char THIS_ATTRIBUTE_NAME[] = "ELASTICITY";
	const char THIS_TYPE_NAME[]      = "constant";
	struct Material_ConstantElasticity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( shear != NULL ) {   *shear = data->shear;   }
	if( bulk != NULL ) {    *bulk = data->bulk;     }

}
/*=================================================================================================================================*/

/*=================================================================================================================================*/
/* temperature-dependent density */
void TemperatureDependentDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char    THIS_TYPE_NAME[]      = "temperature_dependent";
	const size_t  type_size             = sizeof(struct Material_TemperatureDependentDensity);
	/* auxillary params */
	struct Material_TemperatureDependentDensity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }


#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "rho0", fp, start_mat, end_mat, 1.0, &type_data->rho0, &value_found );	 	// essential param
	MaterialTypeParser_FindDouble( "alpha",  fp, start_mat, end_mat, 0.0, &type_data->alpha,  &value_found );  	// non-essential
	MaterialTypeParser_FindDouble( "T0",  fp, start_mat, end_mat, 0.0, &type_data->T0,  &value_found );  		// non-essential


	/* set output */
	*data = (void*)type_data;
}

void MaterialGetTemperatureDependentDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0, double *alpha, double *T0 )
{
	const char THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char THIS_TYPE_NAME[]      = "temperature_dependent";
	struct Material_TemperatureDependentDensity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( rho0 != NULL ) {   	*rho0 	= data->rho0;   }
	if( alpha != NULL ) {  	*alpha 	= data->alpha;  }
	if( T0 != NULL ) {    	*T0 	= data->T0;     }

}
/*=================================================================================================================================*/

/*=================================================================================================================================*/
/* convection density (for use in convection simulations) */
void ConvectionDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char    THIS_TYPE_NAME[]      = "convection";
	const size_t  type_size             = sizeof(struct Material_TemperatureDependentDensity);
	/* auxillary params */
	struct Material_ConvectionDensity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }


#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "Ra",  fp, start_mat, end_mat, 0.0, &type_data->Ra,  &value_found );  		// essential
	MaterialTypeParser_FindDouble( "rho0", fp, start_mat, end_mat, 1.0, &type_data->rho0, &value_found );	 	// non-essential param
	MaterialTypeParser_FindDouble( "T0",  fp, start_mat, end_mat, 0.0, &type_data->T0,  &value_found );  		// non-essential


	/* set output */
	*data = (void*)type_data;
}

void MaterialGetConvectionDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0, double *Ra, double *T0 )
{
	const char THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char THIS_TYPE_NAME[]      = "convection";
	struct Material_ConvectionDensity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( rho0 != NULL ) {   	*rho0 	= data->rho0;   }
	if( Ra != NULL ) {  	*Ra 	= data->Ra;  }
	if( T0 != NULL ) {    	*T0 	= data->T0;     }

}
/*=================================================================================================================================*/
/*=================================================================================================================================*/
/* Artificial temperature dependent density */
void ArtificialTemperatureDependentDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char    THIS_TYPE_NAME[]      = "artificial_temperature_dependent";
	const size_t  type_size             = sizeof(struct Material_ArtificialTemperatureDependentDensity);
	/* auxillary params */
	struct Material_ArtificialTemperatureDependentDensity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }


#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "rho0", fp, start_mat, end_mat, 1.0, &type_data->rho0, &value_found );	 	// essential param

	/* set output */
	*data = (void*)type_data;
}

void MaterialGetArtificialTemperatureDependentDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0)
{
	const char THIS_ATTRIBUTE_NAME[] = "DENSITY";
	const char THIS_TYPE_NAME[]      = "artificial_temperature_dependent";
	struct Material_ArtificialTemperatureDependentDensity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( rho0 != NULL ) {   	*rho0 	= data->rho0;   }


}
/*=================================================================================================================================*/
/*=================================================================================================================================*/
/* Get constant energy parameter */
void ConstantEnergyReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "ENERGY";
	const char    THIS_TYPE_NAME[]      = "constant";
	const size_t  type_size             = sizeof(struct Material_ConstantEnergy);
	/* auxillary params */
	struct Material_ConstantEnergy *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "ThermalConductivity", 	fp, start_mat, end_mat, 1.0, &type_data->ThermalConductivity, 	&value_found );	 	// essential param
	MaterialTypeParser_FindDouble( "HeatCapacity",  		fp, start_mat, end_mat, 0.0, &type_data->HeatCapacity,  				&value_found );  	// essential
	MaterialTypeParser_FindDouble( "RadioactiveHeat",  		fp, start_mat, end_mat, 0.0, &type_data->RadioactiveHeat,  				&value_found );  	// non-essential
	MaterialTypeParser_FindDouble( "ShearHeating",  		fp, start_mat, end_mat, 0.0, &type_data->ShearHeating, 				 &value_found );  	// non-essential


	/* set output */
	*data = (void*)type_data;
}

void MaterialGetConstantEnergyParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *ThermalConductivity,
		double *HeatCapacity, double *RadioactiveHeat, double *ShearHeating )
{
	const char THIS_ATTRIBUTE_NAME[] = "ENERGY";
	const char THIS_TYPE_NAME[]      = "constant";
	struct Material_ConstantEnergy *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( ThermalConductivity != NULL ) 	{ 	*ThermalConductivity 	= data->ThermalConductivity;   	}
	if( HeatCapacity != NULL ) 			{  	*HeatCapacity 			= data->HeatCapacity;  			}
	if( RadioactiveHeat != NULL ) 		{    	*RadioactiveHeat 	= data->RadioactiveHeat;     	}
	if( ShearHeating != NULL ) 			{    	*ShearHeating 		= data->ShearHeating;     		}


}
/*=================================================================================================================================*/

/*=================================================================================================================================*/
void ConstantPlasticityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found )
{
	/* parameters to define type */
	const char    THIS_ATTRIBUTE_NAME[] = "PLASTICITY";
	const char    THIS_TYPE_NAME[]      = "DruckerPrager";
	const size_t  type_size             = sizeof(struct Material_ConstantPlasticity);
	/* auxillary params */
	struct Material_ConstantPlasticity *type_data;
	PetscInt value_found;

	/* initialise */
	*found = _FALSE;
	if( data != NULL ) { *data = NULL; }

	/* check names of attibutes and types match */
	MaterialCompareAttributeTypeNames( attribute,THIS_ATTRIBUTE_NAME, type,THIS_TYPE_NAME, found );
	if( *found == _FALSE ) { return; }

#ifdef _LOG_PARSER
	printf("\tParsing Attribute(%s) => Type(%s)\n", attribute, type );
#endif

	/* allocate space for structure */
	MaterialTypeCreate(M, THIS_ATTRIBUTE_NAME, THIS_TYPE_NAME, type_size, (void**)&type_data );

	/* parse contents */
	MaterialTypeParser_FindDouble( "Cohesion", fp, start_mat, end_mat, 1.0, &type_data->Cohesion, &value_found );
	MaterialTypeParser_FindDouble( "FrictionAngle",  fp, start_mat, end_mat, 0.0, &type_data->FrictionAngle,  &value_found );

	MaterialTypeParser_FindDouble( "Weakening_PlasticStrain_Begin", fp, start_mat, end_mat, 0.0, &type_data->Weakening_PlasticStrain_Begin,&value_found );
	MaterialTypeParser_FindDouble( "Weakening_PlasticStrain_End",   fp, start_mat, end_mat, 0.0, &type_data->Weakening_PlasticStrain_End,  &value_found );
	MaterialTypeParser_FindDouble( "CohesionAfterWeakening",  		fp, start_mat, end_mat, type_data->Cohesion, 		&type_data->CohesionAfterWeakening,  					&value_found );
	MaterialTypeParser_FindDouble( "FrictionAngleAfterWeakening",   fp, start_mat, end_mat, type_data->FrictionAngle, 	&type_data->FrictionAngleAfterWeakening,  &value_found );


	/* set output */
	*data = (void*)type_data;
}

void MaterialGetConstantPlasticityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *Cohesion, double *FrictionAngle,
		double *Weakening_PlasticStrain_Begin, double *Weakening_PlasticStrain_End, double *CohesionAfterWeakening, double *FrictionAngleAfterWeakening)
{
	const char THIS_ATTRIBUTE_NAME[] = "PLASTICITY";
	const char THIS_TYPE_NAME[]      = "DruckerPrager";
	struct Material_ConstantPlasticity *data;
	PetscInt found;


	MaterialGetPhaseAttributeType( M, p_idx,a_idx,t_idx, (void**)&data );
	if( data == NULL ) { return; }

	/* compare attribute and type names */
	MaterialCompareAttributeTypeNames( data->attribute_name,THIS_ATTRIBUTE_NAME, data->type_name,THIS_TYPE_NAME, &found );
	if( found == _FALSE ) { return; }

	if( Cohesion != NULL ) 						{    *Cohesion = data->Cohesion;   }
	if( FrictionAngle != NULL ) 				{    *FrictionAngle = data->FrictionAngle;     }
	if( Weakening_PlasticStrain_Begin != NULL ) {    *Weakening_PlasticStrain_Begin = data->Weakening_PlasticStrain_Begin;     }
	if( Weakening_PlasticStrain_End != NULL ) 	{    *Weakening_PlasticStrain_End = data->Weakening_PlasticStrain_End;     }
	if( CohesionAfterWeakening != NULL ) 		{    *CohesionAfterWeakening = data->CohesionAfterWeakening;     }
	if( FrictionAngleAfterWeakening != NULL ) 	{    *FrictionAngleAfterWeakening = data->FrictionAngleAfterWeakening;     }



}
/*=================================================================================================================================*/

//-----------------------------------------------------------------------------
void char_ptr_stack_create(char_ptr_stack *a)
{
	// initialize stack
	a->cap  = 128;
	a->sz   = 0;
	a->data = (char**) malloc((size_t)a->cap * sizeof(char*));
}
//-----------------------------------------------------------------------------
void char_ptr_stack_delete(char_ptr_stack *a)
{
	// delete stack and its contents
	PetscInt i;
	for(i = 0; i < a->sz; i++)
		free(a->data[i]);
	free(a->data);
}
//-----------------------------------------------------------------------------
void char_ptr_stack_resize(char_ptr_stack *a)
{
	PetscInt i;
	// double capacity of the stack
	a->cap *= 2;
	char **tmp = (char**) malloc((size_t)a->cap * sizeof(char*));
	// copy current contents
	for(i = 0; i < a->sz; i++)
		tmp[i] = a->data[i];
	// clear memory & store new pointer
	free(a->data);
	a->data = tmp;
}
//-----------------------------------------------------------------------------
void char_ptr_stack_insert(char_ptr_stack *a, char *c)
{
	// increase capacity if necessary
	if(a->sz == a->cap)
		char_ptr_stack_resize(a);
	// add new pointer
	a->data[a->sz++] = c;
}
//-----------------------------------------------------------------------------

