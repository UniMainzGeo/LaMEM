/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   tools.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

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


//---------------------------------------------------------------------------


#undef __FUNCT__
#define __FUNCT__ "DMDAGetProcessorRank"
PetscErrorCode DMDAGetProcessorRank(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_col)
{
	PetscMPIInt	rank;
	PetscInt	m, n, p, i, j, k, colind;
	PetscFunctionBegin;
	// get MPI processor rank
	MPI_Comm_rank(((PetscObject)da)->comm, &rank);
	// get number processors in each coordinate direction
	DMDAGetInfo(da, 0, 0, 0, 0, &m, &n, &p, 0, 0, 0, 0, 0, 0);
	// determine i-j-k coordinates of processor
	// x-index runs first (i), then y (j), followed by z (k)
	getLocalRank(&i, &j, &k, rank, m, n);
	// compute index of x-y column of processors
	// (same rule as above for x and y coordinates)
	colind = i + j*m;
	// assign output
	if(rank_x)   (*rank_x)   = i;
	if(rank_y)   (*rank_y)   = j;
	if(rank_z)   (*rank_z)   = k;
	if(rank_col) (*rank_col) = colind;
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "makeMPIIntArray"
PetscErrorCode makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n)
{
	PetscMPIInt    *tmp;
	size_t          sz;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;
	// compute size in bytes
	sz = (size_t)n*sizeof(PetscMPIInt);
	// allocate space
	ierr = PetscMalloc(sz, &tmp); CHKERRQ(ierr);
	// initialize memory from vector (if required)
	if(init) { ierr = PetscMemcpy(tmp, init, sz); CHKERRQ(ierr); }
	// or just clear memory
	else { ierr = PetscMemzero(tmp, sz); CHKERRQ(ierr); }
	// return pointer
	*arr = tmp;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "makeIntArray"
PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n)
{
	PetscInt       *tmp;
	size_t          sz;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;
	// compute size in bytes
	sz = (size_t)n*sizeof(PetscInt);
	// allocate space
	ierr = PetscMalloc(sz, &tmp); CHKERRQ(ierr);
	// initialize memory from vector (if required)
	if(init) { ierr = PetscMemcpy(tmp, init, sz); CHKERRQ(ierr); }
	// or just clear memory
	else { ierr = PetscMemzero(tmp, sz); CHKERRQ(ierr); }
	// return pointer
	*arr = tmp;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "makeScalArray"
PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n)
{
	PetscScalar    *tmp;
	size_t          sz;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;
	// compute size in bytes
	sz = (size_t)n*sizeof(PetscScalar);
	// allocate space
	ierr = PetscMalloc(sz, &tmp); CHKERRQ(ierr);
	// initialize memory from vector (if required)
	if(init) { ierr = PetscMemcpy(tmp, init, sz); CHKERRQ(ierr); }
	// or just clear memory
	else { ierr = PetscMemzero(tmp, sz); CHKERRQ(ierr); }
	// return pointer
	*arr = tmp;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//==========================================================================================================
PetscInt ISRankZero(MPI_Comm comm)
{
	PetscMPIInt rank;

	MPI_Comm_rank(comm, &rank);

	return (rank == 0);
}
//==========================================================================================================
PetscInt ISParallel(MPI_Comm comm)
{
	PetscMPIInt size;

	MPI_Comm_size(comm, &size);

	return (size > 1);
}
//==========================================================================================================
// Creates an output directory
#undef __FUNCT__
#define __FUNCT__ "LaMEMCreateOutputDirectory"
PetscErrorCode LaMEMCreateOutputDirectory(const char *DirectoryName)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// generate a new directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		if(mkdir(DirectoryName, S_IRWXU))
		{
			PetscPrintf(PETSC_COMM_WORLD," Writing to existing directory %s \n", DirectoryName);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD," Created new directory %s \n", DirectoryName);
		}
	}

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================

