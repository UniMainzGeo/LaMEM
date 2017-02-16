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
// Read and write global vectors
//---------------------------------------------------------------------------

PetscErrorCode VecReadRestart (Vec x, FILE *fp);

PetscErrorCode VecWriteRestart(Vec x, FILE *fp);

//---------------------------------------------------------------------------
//  basic statistic functions
//---------------------------------------------------------------------------

PetscScalar getArthMean(PetscScalar *data, PetscInt n);

PetscScalar getVar(PetscScalar *data, PetscInt n);

PetscScalar getStdv(PetscScalar *data, PetscInt n);

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

#define CHECKEQ(a, b, rtol, atol) (PetscAbsScalar((a)-(b)) <= rtol*(PetscAbsScalar(a) + PetscAbsScalar(b)) + atol)

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
//---------------------------------------------------------------------------
// this function computes local ranks of processor in DMDA
static inline void getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n)
{
	(*k) =  rank/(m*n);
	(*j) = (rank - (*k)*m*n)/m;
	(*i) =  rank - (*k)*m*n - (*j)*m;
}
//-----------------------------------------------------------------------------
// calculate displacements from counts, return number elements
PetscInt calcDisp(PetscInt n, PetscInt *counts, PetscInt *displ);

// rewind displacements after using them as access iterators
void rewinDisp(PetscInt n, PetscInt *displ);

//---------------------------------------------------------------------------
// key-value sort using standard functions (scalar-index)
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar key;
	PetscInt    val;

} Pair;

// comparison function for sorting key-value pairs
int comp_key_val(const void * a, const void * b);

PetscErrorCode sort_key_val(PetscScalar *a, PetscInt *idx, PetscInt n);

//---------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------

// compute pointers from counts, return total count
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]);

// rewind pointers after using them as access iterators
void rewindPtr(PetscInt n, PetscInt ptr[]);

// compute phase ratio array
PetscErrorCode getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum);

// find ID of the cell containing point (call this function for local point only!)
static inline PetscInt FindPointInCell(
	PetscScalar *px, // node coordinates
	PetscInt     L,  // index of the leftmost node
	PetscInt     R,  // index of the rightmost node
	PetscScalar  x)  // point coordinate
{
	// get initial guess assuming uniform grid
	PetscInt M = L + (PetscInt)((x-px[L])/((px[R]-px[L])/(PetscScalar)(R-L)));

	if(M == R) return R-1;

	if(px[M]   <= x) L=M;
	if(px[M+1] >= x) R=M+1;

	while((R-L) > 1)
	{
		M = (L+R)/2;
		if(px[M] <= x) L=M;
		if(px[M] >= x) R=M;

	}
	return(L);
}
//-----------------------------------------------------------------------------
static inline PetscScalar InterpLin3D(
	PetscScalar ***lv,
	PetscInt    i,
	PetscInt    j,
	PetscInt    k,
	PetscInt    sx,
	PetscInt    sy,
	PetscInt    sz,
	PetscScalar xp,
	PetscScalar yp,
	PetscScalar zp,
	PetscScalar *cx,
	PetscScalar *cy,
	PetscScalar *cz)
{
	PetscScalar xb, yb, zb, xe, ye, ze, v;

	// get relative coordinates
	xe = (xp - cx[i])/(cx[i+1] - cx[i]); xb = 1.0 - xe;
	ye = (yp - cy[j])/(cy[j+1] - cy[j]); yb = 1.0 - ye;
	ze = (zp - cz[k])/(cz[k+1] - cz[k]); zb = 1.0 - ze;

	// interpolate & return result
	v =
	lv[sz+k  ][sy+j  ][sx+i  ]*xb*yb*zb +
	lv[sz+k  ][sy+j  ][sx+i+1]*xe*yb*zb +
	lv[sz+k  ][sy+j+1][sx+i  ]*xb*ye*zb +
	lv[sz+k  ][sy+j+1][sx+i+1]*xe*ye*zb +
	lv[sz+k+1][sy+j  ][sx+i  ]*xb*yb*ze +
	lv[sz+k+1][sy+j  ][sx+i+1]*xe*yb*ze +
	lv[sz+k+1][sy+j+1][sx+i  ]*xb*ye*ze +
	lv[sz+k+1][sy+j+1][sx+i+1]*xe*ye*ze;

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar InterpLin2D(
	PetscScalar ***lv,
	PetscInt    i,
	PetscInt    j,
	PetscInt    L,
	PetscInt    sx,
	PetscInt    sy,
	PetscScalar xp,
	PetscScalar yp,
	PetscScalar *cx,
	PetscScalar *cy)
{
	PetscScalar xb, yb, xe, ye, v;

	// get relative coordinates
	xe = (xp - cx[i])/(cx[i+1] - cx[i]); xb = 1.0 - xe;
	ye = (yp - cy[j])/(cy[j+1] - cy[j]); yb = 1.0 - ye;

	// interpolate & return result
	v =
	lv[L][sy+j  ][sx+i  ]*xb*yb +
	lv[L][sy+j  ][sx+i+1]*xe*yb +
	lv[L][sy+j+1][sx+i  ]*xb*ye +
	lv[L][sy+j+1][sx+i+1]*xe*ye;

	return v;
}
//-----------------------------------------------------------------------------


#endif
