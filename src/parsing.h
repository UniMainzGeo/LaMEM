/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...................   Input file parsing routines   .......................
//---------------------------------------------------------------------------
#ifndef __parsing_h__
#define __parsing_h__

//-----------------------------------------------------------------------------

struct FDSTAG;

//-----------------------------------------------------------------------------

enum ParamType
{
	_REQUIRED_,
	_OPTIONAL_
};

//-----------------------------------------------------------------------------
// Input file buffer
//-----------------------------------------------------------------------------

struct FB
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
	char      *fbuf;    // file buffer

	char      *lbuf;    // line buffer

	PetscInt   nfLines; // number of flat lines
	char     **pfLines; // pointers to flat lines

	PetscInt   nbLines; // number of block lines
	char     **pbLines; // pointers to block lines

	PetscInt   nblocks; // number of data blocks
	PetscInt   blockID; // ID of block to be retrieved
	PetscInt  *blBeg;   // starting lines of blocks
	PetscInt  *blEnd;   // ending lines of blocks

    PetscInt   ID;      // ID of the current phase or softening law 
};

//-----------------------------------------------------------------------------

PetscErrorCode FBLoad(FB **pfb, PetscBool DisplayOutput, char *restartFileName = NULL);

PetscErrorCode FBDestroy(FB **pfb);

PetscErrorCode FBParseBuffer(FB *fb);

PetscErrorCode FBFindBlocks(FB *fb, ParamType ptype, const char *keybeg, const char *keyend);

PetscErrorCode FBFreeBlocks(FB *fb);

char ** FBGetLineRanges(FB *fb, PetscInt *lnbeg, PetscInt *lnend);

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
		const char  *_default_);  // default value (optional)

//-----------------------------------------------------------------------------
// PETSc options parsing functions
//-----------------------------------------------------------------------------

PetscErrorCode PetscOptionsReadFromFile(FB *fb, PetscBool DisplayOutput);

PetscErrorCode PetscOptionsReadRestart(FILE *fp);

PetscErrorCode PetscOptionsWriteRestart(FILE *fp);

PetscErrorCode  PetscOptionsGetCheckString(
	const char   key[],
	char         str[],
	PetscBool   *set);

//-----------------------------------------------------------------------------
// Default solver options
//-----------------------------------------------------------------------------

PetscErrorCode solverOptionsReadFromFile(FB *fb);

PetscErrorCode solverOptionsSetRequired();

PetscErrorCode set_tolerances(const char *prefix, PetscScalar tolerances[3]);

PetscErrorCode set_default_smoother(
		const char *smoother_type,
		char       *smoother_ksp,
		char       *smoother_pc);

PetscErrorCode set_smoother_options(
		const char *prefix,
		const char *smoother_ksp,
		const char *smoother_pc,
		PetscScalar smoother_damping,
		PetscScalar smoother_omega,
		PetscInt    smoother_num_sweeps,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells);

PetscErrorCode set_subdomain_options(
		const char *prefix,
		const char *smoother_pc,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells);

PetscErrorCode get_num_mg_levles(
		FDSTAG  *fs,
		PetscInt &num_mg_levels);

PetscErrorCode get_coarse_reduction_factor(
		FDSTAG  *fs,
		PetscInt num_mg_levels,
		PetscInt coarse_cells_per_cpu,
		PetscInt &reduction_factor);



PetscErrorCode set_integer_option(const char *key, const PetscInt val, const char *prefix = NULL);

PetscErrorCode set_scalar_option(const char *key, const PetscScalar val, const char *prefix = NULL);

PetscErrorCode set_string_option(const char *key, const char *val, const char *prefix = NULL);

PetscErrorCode set_empty_option(const char *key, const char *prefix = NULL);

//-----------------------------------------------------------------------------
#endif
