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
 **    filename:   tools.h
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
#ifndef __tools_h__
#define __tools_h__

//---------------------------------------------------------------------------
//  basic statistic functions
//---------------------------------------------------------------------------

PetscScalar getArthMean(PetscScalar *data, PetscInt n);

PetscScalar getVar(PetscScalar *data, PetscInt n);

PetscScalar getStdv(PetscScalar *data, PetscInt n);

PetscErrorCode DMDAGetProcessorRank(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_col);

PetscErrorCode makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n);

PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n);

PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n);

//---------------------------------------------------------------------------
// checks whether processor has a zero rank in the communicator
PetscInt ISRankZero(MPI_Comm comm);

// check whether communicator is parallel (has more than one rank)
PetscInt ISParallel(MPI_Comm comm);

PetscErrorCode LaMEMCreateOutputDirectory(const char *DirectoryName);

//---------------------------------------------------------------------------

#define LAMEM_CHECKEQ(a, b, rtol, atol) (PetscAbsScalar((a)-(b)) <= rtol*(PetscAbsScalar(a) + PetscAbsScalar(b)) + atol)

#define IS_POWER_OF_TWO(x) ((x) && !((x) & ((x) - 1)))

//---------------------------------------------------------------------------

void polygon_box(
	PetscInt    *pnv,    // number of polygon vertices (can be modified)
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar  rtol,   // relative tolerance
	PetscScalar *atol,   // absolute tolerance
	PetscScalar *box);   // bounding box of a polygon

void in_polygon(
	PetscInt     np,     // number of test points
	PetscScalar *pcoord, // coordinates of test points
	PetscInt     nv,     // number of polygon vertices
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar *box,    // bounding box of a polygon (optimization)
	PetscScalar  atol,   // absolute tolerance
	PetscInt    *in);    // point location flags (1-inside, 0-outside)

//---------------------------------------------------------------------------

static inline void RotDispPoint2D(PetscScalar Xa[], PetscScalar Xb[], PetscScalar costh, PetscScalar sinth, PetscScalar xa[], PetscScalar xb[])
{
	PetscScalar r[2];

	// get radius vector
	r[0] = xa[0] - Xa[0];
	r[1] = xa[1] - Xa[1];

	// rotate & translate
	xb[0] = costh*r[0] - sinth*r[1] + Xb[0];
	xb[1] = sinth*r[0] + costh*r[1] + Xb[1];
}
//---------------------------------------------------------------------------

static inline PetscScalar ARCCOS(PetscScalar x)
{
	if(x >  1.0 - DBL_EPSILON) x =  1.0 - DBL_EPSILON;
	if(x < -1.0 + DBL_EPSILON) x = -1.0 + DBL_EPSILON;

	return acos(x);
}

//---------------------------------------------------------------------------

static inline PetscScalar ODDROOT(PetscScalar x, PetscScalar a)
{

	if(x < 0.0) return -pow(-x, a);
	else        return  pow( x, a);

}

//---------------------------------------------------------------------------

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

//---------------------------------------------------------------------------
#endif
