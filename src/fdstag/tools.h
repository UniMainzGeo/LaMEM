//---------------------------------------------------------------------------
// ........................... UTILITY FUNCTIONS ............................
//---------------------------------------------------------------------------
#ifndef __tools_h__
#define __tools_h__
//---------------------------------------------------------------------------

/* $Id: tools.h 5681 2015-02-20 20:57:42Z ltbaumann $ */

//---------------------------------------------------------------------------
//  basic statistic functions
//---------------------------------------------------------------------------

PetscScalar getArthMean(PetscScalar *data, PetscInt n);

PetscScalar getVar(PetscScalar *data, PetscInt n);

PetscScalar getStdv(PetscScalar *data, PetscInt n);

//---------------------------------------------------------------------------
// read arrays from PETSC options database with error checking
//---------------------------------------------------------------------------

PetscErrorCode GetScalArrayCheckScale(
	const char  ident[],
	const char  name[],
	PetscInt    n,
	PetscScalar a[],
	PetscScalar amin,
	PetscScalar amax,
	PetscScalar scal);

PetscErrorCode GetIntArrayCheck(
	const char  ident[],
	const char  name[],
	PetscInt    n,
	PetscInt    a[],
	PetscInt    amin,
	PetscInt    amax);

//---------------------------------------------------------------------------
#endif
