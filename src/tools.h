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
// Printing functions
//---------------------------------------------------------------------------

void PrintStart(PetscLogDouble *t_beg, const char *msg, const char *filename);

void PrintDone(PetscLogDouble t_beg);

void PrintStep(PetscInt step);

//---------------------------------------------------------------------------
// Read and write global vectors
//---------------------------------------------------------------------------

PetscErrorCode VecReadRestart (Vec x, FILE *fp);

PetscErrorCode VecWriteRestart(Vec x, FILE *fp);

//---------------------------------------------------------------------------
// Basic statistic functions
//---------------------------------------------------------------------------

PetscScalar getArthMean(PetscScalar *data, PetscInt n);

PetscScalar getVar(PetscScalar *data, PetscInt n);

PetscScalar getStdv(PetscScalar *data, PetscInt n);

PetscErrorCode makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n);

PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n);

PetscErrorCode clearIntArray(PetscInt *arr, const PetscInt n);

PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n);

//---------------------------------------------------------------------------
// Rank checking functions
//---------------------------------------------------------------------------

// checks whether processor has a zero rank in the communicator
PetscInt ISRankZero(MPI_Comm comm);

// check whether communicator is parallel (has more than one rank)
PetscInt ISParallel(MPI_Comm comm);

// get global rank of processor in DMDA
static inline PetscMPIInt getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p)
{
	if (i < 0 || i >= m || j < 0 || j >= n || k < 0 || k >= p) return -1;
	return (PetscMPIInt)(i + j*m + k*m*n);
}

// get local ranks of processor in DMDA
static inline void getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n)
{
	(*k) =  rank/(m*n);
	(*j) = (rank - (*k)*m*n)/m;
	(*i) =  rank - (*k)*m*n - (*j)*m;
}

//---------------------------------------------------------------------------
// Directory management functions
//---------------------------------------------------------------------------

PetscErrorCode DirMake(const char *name);

PetscErrorCode DirRemove(const char *name);

PetscErrorCode DirRename(const char *old_name, const char *new_name);

PetscErrorCode DirCheck(const char *name, PetscInt *exists);

//---------------------------------------------------------------------------
// Numerical functions
//---------------------------------------------------------------------------

#define CHECKEQ(a, b, rtol, atol) (PetscAbsScalar((a)-(b)) <= (rtol)*(PetscAbsScalar(a) + PetscAbsScalar(b)) + (atol))

#define IS_POWER_OF_TWO(x) ((x) && !((x) & ((x) - 1)))

static inline PetscScalar ARCCOS(PetscScalar x)
{
	if(x >  1.0 - DBL_EPSILON) x =  1.0 - DBL_EPSILON;
	if(x < -1.0 + DBL_EPSILON) x = -1.0 + DBL_EPSILON;

	return acos(x);
}

static inline PetscScalar ODDROOT(PetscScalar x, PetscScalar a)
{
	if(x < 0.0) return -pow(-x, a);
	else        return  pow( x, a);
}

//---------------------------------------------------------------------------
// Polygon location functions
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
// Polygon stretching functions
//---------------------------------------------------------------------------

// generate linearly interpolated values
void linSpace(
	PetscScalar  min,
	PetscScalar  max,
	PetscInt     N,
	PetscScalar *outVec);

// interpolate stretch parameters for all polygons
void interpStretch(
	PetscScalar *Sx,
    PetscScalar *Sy,
    PetscInt     numCtrlPoly,
    PetscInt    *CtrlPoly,
    PetscInt     numPoly,
    PetscScalar *SxAll,
    PetscScalar *SyAll);

// find center of mass of polygon
void findCenterMass(
	PetscScalar *coords,
	PetscInt     nN,
	PetscScalar &x_cen,
	PetscScalar &y_cen);

// stretch Polygon
void stretchPolygon(
	PetscScalar *coords,
	PetscInt nN,
	PetscScalar Sx,
	PetscScalar Sy);

//---------------------------------------------------------------------------
// indexing functions
//---------------------------------------------------------------------------

// compute pointers from counts, return total count
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]);

// rewind pointers after using them as access iterators
void rewindPtr(PetscInt n, PetscInt ptr[]);

//-----------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------

// compute phase ratio array
PetscErrorCode getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum);

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

// bisection algorithm for scalar nonlinear equation
PetscInt solveBisect(
		PetscScalar a,
		PetscScalar b,
		PetscScalar tol,
		PetscScalar maxit,
		PetscScalar &x,
		PetscInt    &it,
		PetscScalar (*f) (PetscScalar x, void *pctx),
		void *pctx);

//---------------------------------------------------------------------------
// Interpolation functions
//---------------------------------------------------------------------------

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
