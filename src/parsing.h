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

enum StokesSolverType
{
	_DIRECT_STOKES_,
	_MULTIGRID_STOKES_,
	_BLOCK_STOKES_,
	_wBFBT_STOKES_
};

enum CoarseSolverType
{
	_DIRECT_COARSE_,
	_HYPRE_COARSE_,
	_ASM_COARSE_,
	_BJACOBI_COARSE_
};

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
// Set default solver options
//-----------------------------------------------------------------------------

PetscErrorCode solverOptionsReadFromFile(FB *fb);

PetscErrorCode solverOptionsSetRequired(FB *fb);

//-----------------------------------------------------------------------------
#endif
