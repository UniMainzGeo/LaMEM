//---------------------------------------------------------------------------
//.....................   OBJECTIVE FUNCTION ROUTINES   .....................
//---------------------------------------------------------------------------
#ifndef __objFunct_h__
#define __objFunct_h__
//---------------------------------------------------------------------------
// maximum number of obervational types
#define _max_num_obs_ 7
#define _max_len_name_ 8


typedef enum 	// observation type
{
	_VELX_,            // 1: horizontal velocity (x) at the surface
	_VELY_,            // 2: horizontal velocity (y) at the surface
	_VELZ_,            // 3: vertical velocity (z) at the surface
	_TOPO_,            // 4: surface topography
	_BOUG_,            // 5: bouguer anomaly
	_ISA_,             // 6: orientation of isa (<-> sks-seismic anisotropy)
	_SHMAX_            // 7: orientation of SHmax
} ObsType;

//---------------------------------------------------------------------------
//........................ Objective function object ........................
//---------------------------------------------------------------------------
typedef struct
{
	FreeSurf     *surf;                 // free surface object
	char         *infile;               // input file name
	PetscBool    CompMfit;              // Compute misfit?
	PetscInt     otUse[_max_num_obs_+1];// array of boolean USED flags
	PetscInt     otN;                   // number of USED observation types
	PetscInt     ocN;                   // total number of observational constraints
	PetscScalar  err[_max_num_obs_];    // array containing individual sums of errors
	PetscScalar  errtot;                // total error
	Vec          obs[_max_num_obs_];    // vectors containing the observations
	Vec          qul[_max_num_obs_];    // vectors containing quality info (quality/sigma)^2, where quality (0..1)

	// missing ...
	// (data) covariance matrix

} ObjFunct;
//---------------------------------------------------------------------------

// destroy object
PetscErrorCode ObjFunctDestroy(ObjFunct *objf);

// create objective function object
PetscErrorCode ObjFunctCreate(ObjFunct *objf, FreeSurf *surf);

// read command line options
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[]);

// compute error
PetscErrorCode ObjFunctCompErr(ObjFunct *objf);

// compute weighted Least square error for surface vectors
PetscErrorCode VecErrSurf(Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal);

//---------------------------------------------------------------------------

#endif
