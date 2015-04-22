/* Definitions for use with LaMEM */
/* $Id$ */

#ifndef __LaMEM_h__
#define __LaMEM_h__

//-----------------------------------------------------------------------------
// DEFINITIONS
//-----------------------------------------------------------------------------

// temporary global variable
#define ELEMENT_Q1P0		2
#define ELEMENT_Q1Q1		3
#define ELEMENT_Q2P1		1
#define ELEMENT_FDSTAG		4
#define ELEMENT_Q2P1_LOCAL	22
#define ELEMENT_Q2P1_GLOBAL	33

// set maximums based on Q2 storage requirements
#define MAX_nnel 27
#define MAX_ngp_vel 27
#define MAX_nintp_1D 3
#define MAX_npres 8
#define MAX_edof 81
#define MAX_edof_temp 27
#define MAX_nnode_el_1D 2
#define MAX_nnel_1D 3
#define MAX_ndim 3

#ifdef PARTICLES
#define particle_buffer 1.25
#define particle_props 28
#endif

#define max_num_phases 32 // max no of phases
#define max_num_soft   10 // max no of soft laws
#define MaxNumCPU      524288
#define MaxNumSteps    1000000

// maximum number of mesh segments in every direction
#define MaxNumMeshSegs 10

#define NumMaterialPropsElem 33
#define NumTimeDepData 	 (8+6*max_num_phases + 2*6)

// cast macros
#define LLD long long int

// space dimension
#define SPDIM 3

// use this to enable asprintf
#define _GNU_SOURCE

// Identify gcc compiler.
// Take care that other compilers may also define __GNUC__ macro.

#undef GCC_COMPILER

#if defined (__GNUC__) && !defined (__INTEL_COMPILER)

	#define GCC_COMPILER

#endif

#define MAX_NAME_LEN 64
#define MAX_PATH_LEN 256

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
	#pragma GCC diagnostic ignored "-Wshadow"

	#include <petsc.h>
	#include <petsc-private/matimpl.h>

	#pragma GCC diagnostic warning "-Wconversion"
	#pragma GCC diagnostic warning "-Wsign-conversion"
	#pragma GCC diagnostic warning "-Wunused-parameter"
	#pragma GCC diagnostic warning "-Wcast-qual"
	#pragma GCC diagnostic warning "-Wshadow"

#else

	#include <petsc.h>
	#include <petsc-private/matimpl.h>

#endif

//-----------------------------------------------------------------------------
// INTERNAL TYPE DEFINITIONS
//-----------------------------------------------------------------------------

#include "Attributes.h"

// used only in FDSTAG Canonical
#include "fdstagTypes.h"

// LaMEM Legacy
#include "LaMEM_Types.h"
#include "LaMEMVelPressureDA.h"
#include "LaMEMVelPressureDA_private.h"

//-----------------------------------------------------------------------------
// PROTOTYPES
//-----------------------------------------------------------------------------

// LaMEM library main function
PetscErrorCode LaMEMLib(PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id);
PetscErrorCode LaMEMLib_FDSTAG(void *echange_ctx);
PetscErrorCode LaMEMLib_Legacy(PetscBool InputParamFile, const char *ParamFile, PetscScalar *LaMEM_OutputParameters, PetscInt *mpi_group_id);

//-----------------------------------------------------------------------------
// GLOBAL VARIABLES
//-----------------------------------------------------------------------------

extern PetscInt __ELEMENT_TYPE__;
extern PetscInt __Q2_TYPE__;

//-----------------------------------------------------------------------------
// MACROS
//-----------------------------------------------------------------------------

#define EMERGENCY_EXIT(message){SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-> LaMEM error encountered - %s", message);}

// invert the value of PETSC_BOOL variable
#define PETSC_NEGATE(a) ((a == PETSC_TRUE ) ? PETSC_FALSE : PETSC_TRUE)

#define LAMEM_FREE(a) if(a) { free(a); a = NULL; }

#define LAMEM_CHECKEQ(a, b, rtol, atol) (PetscAbsScalar((a)-(b)) <= rtol*(PetscAbsScalar(a) + PetscAbsScalar(b)) + atol)

#define IS_POWER_OF_TWO(x) ((x) && !((x) & ((x) - 1)))

//-----------------------------------------------------------------------------
// INLINE FUNCTIONS
//-----------------------------------------------------------------------------

// this function returns global rank of processor in DMDA
static inline PetscMPIInt getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p)
{
	if (i < 0 || i >= m || j < 0 || j >= n || k < 0 || k >= p) return -1;
	return (PetscMPIInt)(i + j*m + k*m*n);
}

// this function computes local ranks of processor in DMDA
static inline void getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n)
{
	(*k) =  rank/(m*n);
	(*j) = (rank - (*k)*m*n)/m;
	(*i) =  rank - (*k)*m*n - (*j)*m;
}

// Modified bisection algorithm (ltbaumann 210113)
// Returns index i of the closet gridpoint x_i of an arbitrary value x
// px - 1D-grid coordinates
// L  - first index
// R  - last index
static inline PetscInt Bisection(PetscScalar *px, PetscInt L, PetscInt R, PetscScalar x)
{
	PetscInt M;
	while((R-L) > 1)
	{	M = (L+R)/2;
		if(px[M] <= x) L=M;
		if(px[M] >= x) R=M;
	}
	if(PetscAbsScalar(px[L]-x) <= PetscAbsScalar(px[R]-x)) return(L);
	else                                                   return(R);
}
//-----------------------------------------------------------------------------
static inline PetscErrorCode sfexp(PetscScalar x, PetscScalar *y)
{
	// y = e^x with checking range errors
	errno = 0;
	(*y) = exp(x);
	if(errno == EDOM)   { PetscPrintf(PETSC_COMM_WORLD,"Domain Error!\n"); return(-1); }
	if(errno == ERANGE) { PetscPrintf(PETSC_COMM_WORLD,"Range Error!\n");  return(-1); }
	return(0);
/*
	return (-1);
	feclearexcept(FE_ALL_EXCEPT);
	fe = fetestexcept (FE_ALL_EXCEPT);
	if (fe & FE_DIVBYZERO) puts ("FE_DIVBYZERO");
	if (fe & FE_INEXACT)   puts ("FE_INEXACT");
	if (fe & FE_INVALID)   puts ("FE_INVALID");
	if (fe & FE_OVERFLOW)  puts ("FE_OVERFLOW");
	if (fe & FE_UNDERFLOW) puts ("FE_UNDERFLOW");
*/
}
//---------------------------------------------------------------------------
static inline PetscErrorCode sfpow(PetscScalar a, PetscScalar x, PetscScalar *y)
{
	// y = a^x with checking range errors
	errno = 0;
	(*y) = pow(a, x);
	if(errno == EDOM)   { PetscPrintf(PETSC_COMM_WORLD,"Domain Error!\n"); return(-1); }
	if(errno == ERANGE) { PetscPrintf(PETSC_COMM_WORLD,"Range Error!\n");  return(-1); }
	return (0);
}
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
#endif
