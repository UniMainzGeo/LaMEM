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
#define __FUNCT__ "GetScalArrayCheckScale"
PetscErrorCode GetScalArrayCheckScale(
	const char  ident[],
	const char  name[],
	PetscInt    n,
	PetscScalar a[],
	PetscScalar amin,
	PetscScalar amax,
	PetscScalar scal)
{
	PetscInt  i;
	PetscBool set;
	PetscInt  nmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// skip trivial case
	if(!n) PetscFunctionReturn(0);

	// set expected number of elements
	nmax = n;

	// read array
	ierr = PetscOptionsGetRealArray(NULL, ident,  a, &nmax, &set); CHKERRQ(ierr);

	// check array exists
	if(set != PETSC_TRUE)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Array \"%s\" is not found\n", name);
	}

	// check correct number of elements is provided
	if(nmax != n)
	{
		SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_USER, "Wrong number of elements in array \"%s\" , actual: %lld, expected: %lld\n", (LLD)nmax, (LLD)n);
	}

	// check ranges
	if(amin && amax)
	{
		for(i = 0; i < n; i++)
		{
			if(a[i] < amin || a[i] > amax)
			{
				SETERRQ5(PETSC_COMM_SELF, PETSC_ERR_USER, "Element %lld of array \"%s\" is out of bound, actual: %e, range: [%e - %e]\n",
					(LLD)i, name, a[i], amin, amax);
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
#define __FUNCT__ "GetIntArrayCheck"
PetscErrorCode GetIntArrayCheck(
	const char  ident[],
	const char  name[],
	PetscInt    n,
	PetscInt    a[],
	PetscInt    amin,
	PetscInt    amax)
{
	PetscInt  i;
	PetscBool set;
	PetscInt  nmax;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// skip trivial case
	if(!n) PetscFunctionReturn(0);

	// set expected number of elements
	nmax = n;

	// read array
	ierr = PetscOptionsGetIntArray(NULL, ident,  a, &nmax, &set); CHKERRQ(ierr);

	// check array exists
	if(set != PETSC_TRUE)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Array \"%s\" is undefined\n", name);
	}

	// check correct number of elements is provided
	if(nmax != n)
	{
		SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_USER, "Wrong number of elements in array \"%s\" , actual: %lld, expected: %lld\n", (LLD)nmax, (LLD)n);
	}

	// check ranges
	if(amin && amax)
	{
		for(i = 0; i < n; i++)
		{
			if(a[i] < amin || a[i] > amax)
			{
				SETERRQ5(PETSC_COMM_SELF, PETSC_ERR_USER, "Element %lld of array \"%s\" is out of bound, actual: %lld, range: [%lld - %lld]\n",
					(LLD)i, name, (LLD)a[i], (LLD)amin, (LLD)amax);
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
