#ifndef __Attributes_h__
#define __Attributes_h__

//-----------------------------------------------------------------------------
typedef struct
{
	PetscInt cap;  // current capacity
	PetscInt sz;   // current size
	char ** data; // data storage
} char_ptr_stack;
//-----------------------------------------------------------------------------
void char_ptr_stack_create(char_ptr_stack *a);
void char_ptr_stack_delete(char_ptr_stack *a);
void char_ptr_stack_resize(char_ptr_stack *a);
void char_ptr_stack_insert(char_ptr_stack *a, char *c);
//-----------------------------------------------------------------------------

/* material Attributes type helper functions */
typedef struct _p_AttributeType* AttributeType;
struct _p_AttributeType {
	char *name;
	PetscInt T_id;
	void *data;
};


typedef struct _p_Attribute* Attribute;
struct _p_Attribute {
	char *name;
	PetscInt A_id;
	PetscInt n_types;
	AttributeType *type;
};

typedef struct _p_Phase* Phase;
struct _p_Phase {
	char *name;
	PetscInt phase_number;
	PetscInt P_id;
	PetscInt n_attributes;
	Attribute *attribute;
};

typedef struct _p_Material* Material;
struct _p_Material {
	PetscInt n_phases;
	Phase *phase;

	PetscInt n_attributes_MAX;
	PetscInt n_types_MAX;
	size_t size_ptr;
	void ****material_phase_table;  /* PHASE, ATTRIBUTE, TYPE*/

	char_ptr_stack a;
};

#define _MaterialTypeBase \
	char *attribute_name; \
	char *type_name;
struct _p_MaterialTypeBase { _MaterialTypeBase };

void MaterialReadFromFile(
		Material M, const char filename[],
		void (*fp_parse_material_type)(Material M, const char*,const char*,FILE*,const long,const long,void**,PetscInt*),
		PetscInt *found );

void MaterialCreate( Material *m );
void MaterialDestroy( Material M );
void MaterialView( Material m );
void MaterialGetAllPhases( Material M, PetscInt *n, Phase **phases );
void PhaseGetAllAttributes( Phase p, PetscInt *n, Attribute **attrs );
void AttributeGetAllTypes( Attribute a, PetscInt *n, AttributeType **types );
void AttributeCompare( Attribute a, const char name[], PetscInt *truth );
void AttributeTypeCompare( AttributeType t, const char name[], PetscInt *truth );
void Material_helper_FindInt( const char key[], FILE* fp, const long start, const long end, PetscInt *val, PetscInt *found );
void Material_helper_FindDouble( const char key[], FILE* fp, const long start, const long end, double *val, PetscInt *found );
void MaterialAddPhase( Material m, const char phase_name[], const PetscInt phase_number);
void MaterialAddAttribute( Material m, const char phase_name[], const char attr_name[]);
void MaterialAddAttributeType( Material m, const char phase_name[], const char attr_name[], const char type_name[], void *data);
void AttributeExist( Phase p, const char name[], PetscInt *truth );
void MaterialGetPhaseAttributeType( Material m, const PetscInt phase_idx, const PetscInt attr_idx, const PetscInt type_idx, void **data );

/* constant viscosity definition -----------------------------------*/
struct Material_ConstantViscosity {
	_MaterialTypeBase
	double eta0;
};

void ConstantViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetConstantViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *eta0 );
/*-------------------------------------------------------------------*/


/* power-law viscosity definition -----------------------------------*/
struct Material_PowerLawViscosity {
	_MaterialTypeBase
	double eta0, N_exp, e0;
};

void PowerLawViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetPowerLawViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *eta0, double *N_exp, double *e0 );
/*-------------------------------------------------------------------*/


/* temperature-dependent viscosity definition -----------------------------------*/
struct Material_TempDepViscosity {
	_MaterialTypeBase
	double PreExpFactor, N_exp, e0, ActivationEnergy;
};

void TempDepViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetTempDepViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *PreExpFactor, double *N_exp, double *e0, double *ActivationEnergy );
/*-------------------------------------------------------------------*/


/* temperature-dependent without power-law viscosity definition -----------------------------------*/
struct Material_TempDepNoPowerLawViscosity {
	_MaterialTypeBase
	double PreExpFactor, N_exp, e0, ActivationEnergy;
};

void TempDepNoPowerLawViscosityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetTempDepNoPowerLawViscosityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *PreExpFactor, double *N_exp, double *e0, double *ActivationEnergy );
/*-------------------------------------------------------------------*/


/* constant elasticity definition -----------------------------------*/
struct Material_ConstantElasticity {
	_MaterialTypeBase
	double shear, bulk;
};

void ConstantElasticityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetConstantElasticityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *shear, double *bulk );
/*-------------------------------------------------------------------*/


/* temperature dependent density definition -------------------------*/
struct Material_TemperatureDependentDensity {
	_MaterialTypeBase
	double rho0, alpha, T0;
};

void TemperatureDependentDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetTemperatureDependentDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0, double *alpha, double *T0 );
/*-------------------------------------------------------------------*/

/* convection dependent density definition -------------------------*/
struct Material_ConvectionDensity {
	_MaterialTypeBase
	double rho0, Ra, T0;
};

void ConvectionDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetConvectionDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0, double *Ra, double *T0 );
/*-------------------------------------------------------------------*/

/* artificial temperature dependent density definition --------------*/
struct Material_ArtificialTemperatureDependentDensity {
	_MaterialTypeBase
	double rho0;
};

void ArtificialTemperatureDependentDensityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetArtificialTemperatureDependentDensityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *rho0 );
/*-------------------------------------------------------------------*/

/* constant energy parameters -------------------------*/
struct Material_ConstantEnergy {
	_MaterialTypeBase
	double ThermalConductivity, HeatCapacity, RadioactiveHeat, ShearHeating;
};

void ConstantEnergyReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetConstantEnergyParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *ThermalConductivity,
		double *HeatCapacity, double *RadioactiveHeat, double *ShearHeating );
/*-------------------------------------------------------------------*/

/* constant plasticity definition -----------------------------------*/
struct Material_ConstantPlasticity {
	_MaterialTypeBase
	double Cohesion, FrictionAngle, Weakening_PlasticStrain_Begin, Weakening_PlasticStrain_End, CohesionAfterWeakening, FrictionAngleAfterWeakening;
};

void ConstantPlasticityReadFromFile(Material M,
		const char attribute[], const char type[],
		FILE *fp, const long start_mat, const long end_mat,
		void **data, PetscInt *found );

void MaterialGetConstantPlasticityParams( Material M, const PetscInt p_idx, const PetscInt a_idx, const PetscInt t_idx, double *Cohesion, double *FrictionAngle,
			double *Weakening_PlasticStrain_Begin, double *Weakening_PlasticStrain_End, double *CohesionAfterWeakening, double *FrictionAngleAfterWeakening );
/*-------------------------------------------------------------------*/

void _PhaseInit( Phase p );
void _AttributeInit( Attribute a );
void _AttributeTypeInit( AttributeType t );
void _MaterialGetPhase( Material m, const char phase_name[], Phase *phase );
void _PhaseGetAttribute( Phase phase, const char attr_name[], Attribute *attr );
void _AttributeGetType( Attribute attr, const char type_name[], AttributeType *type );
void MaterialQueryGetIndices( Material m, const char phase_name[], const char attr_name[], const char type_name[], PetscInt *p_idx, PetscInt *a_idx, PetscInt *t_idx );
void MaterialGetPhaseByName( Material m, const char name[], PetscInt *P_id, Phase *phase );
void MaterialGetPhaseByID( Material m, PetscInt phase_id, PetscInt *P_id, Phase *phase );
void PhaseGetAttributeByName( Phase p, const char name[], PetscInt *A_id, Attribute *attr );
void AttributeGetTypeByName( Attribute a, const char name[], PetscInt *T_id, AttributeType *type );
void AttributeTypeGetData( AttributeType t, void **data );
void PhaseCompare( Phase p, const char name[], PetscInt *truth );
void MaterialUpdateData( Material m,  PetscInt p, PetscInt a, PetscInt t, void *data);
void MaterialTypeCreate(Material M, const char attr_name[], const char type_name[], size_t size, void **type_data );
void MaterialCompareAttributeTypeNames(const char attribute_a[], const char attribute_b[], const char type_a[], const char type_b[], PetscInt *match );
void Material_helper_CheckType( const char sought_type[], FILE* fp, const long start, const long end, PetscInt *found );
void _MaterialCreateTable( Material m );
void test_map( void );

/*-------------------------------------------------------------------*/
#endif
