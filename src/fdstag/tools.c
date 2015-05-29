//---------------------------------------------------------------------------
// ........................... UTILITY FUNCTIONS ............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tools.h"
//---------------------------------------------------------------------------
//  basic statistic functions
//---------------------------------------------------------------------------
PetscScalar getArthMean(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar sum = 0.0;

    for (k=0; k<n; k++)
        sum += data[k];
    return (sum/(PetscScalar)n);
}
//---------------------------------------------------------------------------
PetscScalar getVar(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar mean = getArthMean(data,n);
    PetscScalar temp = 0.0;
//	PetscPrintf(PETSC_COMM_WORLD,"mean=%g \n",mean);
    for (k=0; k<n; k++)
        temp += (mean-data[k])*(mean-data[k]);
    return (temp/(PetscScalar)n);
}
//---------------------------------------------------------------------------
PetscScalar getStdv(PetscScalar *data, PetscInt n)
{
    return sqrt(getVar(data,n));
}
//---------------------------------------------------------------------------
// read arrays from PETSC options database with error checking
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetScalDataItemCheckScale"
PetscErrorCode GetScalDataItemCheckScale(
	const char  ident[],
	const char  name[],
	exitType    extp,
	PetscInt    n,
	PetscScalar *a,
	PetscScalar amin,
	PetscScalar amax,
	PetscScalar scal)
{
	PetscInt  i;
	PetscBool found;
	PetscInt  nmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// skip trivial case
	if(!n) PetscFunctionReturn(0);

	// set expected number of elements
	nmax = n;

	// read array
	ierr = PetscOptionsGetRealArray(NULL, ident,  a, &nmax, &found); CHKERRQ(ierr);

	// check data item exists
	if(found != PETSC_TRUE)
	{
		if(extp == _NOT_FOUND_ERROR_)
		{
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Data item \"%s\" is not found\n", name);
		}
		else if(extp == _NOT_FOUND_EXIT_)
		{
			PetscFunctionReturn(0);
		}
	}

	// check correct number of elements is provided
	if(nmax != n)
	{
		SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "Wrong number of elements in data item \"%s\" , actual: %lld, expected: %lld\n", name, (LLD)nmax, (LLD)n);
	}

	// check ranges
	if(amin && amax)
	{
		for(i = 0; i < n; i++)
		{
			if(a[i] < amin || a[i] > amax)
			{
				SETERRQ5(PETSC_COMM_SELF, PETSC_ERR_USER, "Data item \"%s\" is out of bound, actual: %e, range: [%e - %e], element %lld\n",
					name, a[i], amin, amax, (LLD)i);
			}
		}
	}

	// scale
	if(scal)
	{
		for(i = 0; i < n; i++) a[i] /= scal;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetIntDataItemCheck"
PetscErrorCode GetIntDataItemCheck(
	const char  ident[],
	const char  name[],
	exitType    extp,
	PetscInt    n,
	PetscInt    *a,
	PetscInt    amin,
	PetscInt    amax)
{
	PetscInt  i;
	PetscBool found;
	PetscInt  nmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// skip trivial case
	if(!n) PetscFunctionReturn(0);

	// set expected number of elements
	nmax = n;

	// read array
	ierr = PetscOptionsGetIntArray(NULL, ident,  a, &nmax, &found); CHKERRQ(ierr);

	// check data item exists
	if(found != PETSC_TRUE)
	{
		if(extp == _NOT_FOUND_ERROR_)
		{
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Data item \"%s\" is not found\n", name);
		}
		else if(extp == _NOT_FOUND_EXIT_)
		{
			PetscFunctionReturn(0);
		}
	}

	// check correct number of elements is provided
	if(nmax != n)
	{
		SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_USER, "Wrong number of elements in data item \"%s\", actual: %lld, expected: %lld\n", name, (LLD)nmax, (LLD)n);
	}

	// check ranges
	if(amin && amax)
	{
		for(i = 0; i < n; i++)
		{
			if(a[i] < amin || a[i] > amax)
			{
				SETERRQ5(PETSC_COMM_SELF, PETSC_ERR_USER, "Data item \"%s\" is out of bound, actual: %lld, range: [%lld - %lld], element %lld\n",
					name, (LLD)a[i], (LLD)amin, (LLD)amax, (LLD)i);
			}
		}
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
