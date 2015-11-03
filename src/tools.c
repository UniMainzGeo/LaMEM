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
//---------------------------------------------------------------------------
//
// Fast detection points inside a polygonal region
//
// Originally written as a MATLAB mexFunction by:
//
// A. David Redish      email: adr@nsma.arizona.edu
// Guillaume Jacquenot  email: guillaume.jacquenot@gmail.com
//
// Modified to be callable directly from C by:
//
// Anton A. Popov
//
// Output is simplified to skip distinguishing points IN and ON the polygon
//
//---------------------------------------------------------------------------
void polygon_box(
	PetscInt    *pnC,
	PetscScalar *cx,
	PetscScalar *cy,
	PetscScalar *box)
{
	PetscInt    iC, nC;
	PetscScalar ax, bx, ay, by;
	PetscScalar xmin, xmax, ymin, ymax;

	nC = (*pnC);
	ax = cx[0];
	ay = cy[0];
	bx = cx[nC - 1];
	by = cy[nC - 1];

	// decrease number of points if the polygon is closed
	if(ax == bx && ay == by) nC--;

	// calculate bounding box of a polygon
	xmin = xmax = cx[0];
	ymin = ymax = cy[0];

	for(iC = 0; iC < nC; iC++)
	{
		if(cx[iC] < xmin) xmin = cx[iC];
		if(cx[iC] > xmax) xmax = cx[iC];
		if(cy[iC] < ymin) ymin = cy[iC];
		if(cy[iC] > ymax) ymax = cy[iC];
	}

	box[0] = xmin;
	box[1] = xmax;
	box[2] = ymin;
	box[3] = ymax;

	(*pnC) = nC;
}
//---------------------------------------------------------------------------
void in_polygon(
	PetscInt     nP,   // number of points
	PetscScalar *px,   // x-coordinates of points
	PetscScalar *py,   // y-coordinates of points
	PetscInt     nC,   // number of polygon vertices
	PetscScalar *cx,   // x-coordinates of polygon vertices
	PetscScalar *cy,   // y-coordinates of polygon vertices
	PetscScalar *box,  // bounding box of a polygon (optimization)
	PetscScalar  gtol, // geometry tolerance
	PetscInt    *in)   // point location flags (1-inside, 0-outside)
{
	PetscInt    iP, iC, ind;
	PetscInt    points_on, points_in;
	PetscScalar xmin, xmax, ymin, ymax;
	PetscScalar ax, bx, ay, by;
	PetscScalar nIntersect, intersecty, tmp;

	// get bounding box
	xmin = box[0];
	xmax = box[1];
	ymin = box[2];
	ymax = box[3];

	// test whether each point is in polygon
	for(iP = 0; iP < nP; iP++)
	{
		// assume point is outside
		in[iP] = 0;

		// check bounding box
		if(px[iP] < xmin) continue;
		if(px[iP] > xmax) continue;
		if(py[iP] < ymin) continue;
		if(py[iP] > ymax) continue;

		// count the number of intersections
		nIntersect = 0.0;
		points_on  = 0;

		for(iC = 0; iC < nC; iC++)
		{
			// does the line PQ intersect the line AB?
			if(iC == nC-1)
			{
				ax = cx[nC-1];
				ay = cy[nC-1];
				bx = cx[0];
				by = cy[0];
			}
			else
			{
				ax = cx[iC];
				ay = cy[iC];
				bx = cx[iC+1];
				by = cy[iC+1];
			}

			if(ax == bx)
			{
				// Vertical points
				if(px[iP] == ax)
				{
					// ensure order correct
					if(ay > by)
					{
						tmp = ay; ay = by; by = tmp;
					}
					if(py[iP] >= ay && py[iP] <= by)
					{
						points_on  = 1;
						nIntersect = 0.0;
						break;
					}
				}
			}
			else
			{
				// Non Vertical points
				if( px[iP] < MIN(ax, bx) || MAX(ax, bx) < px[iP]) continue;

				intersecty = ay + (px[iP] - ax)/(bx - ax)*(by - ay);

				if(fabs(intersecty - py[iP]) < gtol)
				{
					points_on  = 1;
					nIntersect = 0.0;
					break;
				}
				else if(intersecty < py[iP] && (ax == px[iP] || bx == px[iP]))
				{
					if(ax == px[iP])
					{
						if(iC == 0)
						{
							ind = nC-1;
						}
						else
						{
							ind = iC-1;
						}
						if(MIN(bx, cx[ind]) < px[iP] && px[iP] < MAX(bx, cx[ind]))
						{
							nIntersect += 1.0;
						}
					}
				}
				else if (intersecty < py[iP])
				{
					nIntersect += 1.0;
				}
			}
		}

		// check if the contour polygon is closed
		points_in = (PetscInt)(nIntersect - 2.0*floor(nIntersect/2.0));
		in[iP]    = MAX(points_on, points_in);
	}
}
//---------------------------------------------------------------------------
/*

	// polygon
	PetscInt    i;
	PetscInt    nC   = 5;
	PetscScalar cx[] = { 0.0, 2.0, 1.0, 0.0, 0.0 };
	PetscScalar cy[] = { 0.0, 0.0, 1.0, 1.0, 0.0 };

	PetscScalar box[4];

	polygon_box(&nC, cx, cy, box);

	printf("Bounding box: xmin=%f, xmax=%f, ymin=%f, ymax=%f\n", box[0], box[1], box[2], box[3]);

	printf("Number of vertices: %d\n", nC);

	// points

	PetscInt nP = 4;

	PetscScalar px[] = { 0.5, 0.0, 2.0, 1.6 };
	PetscScalar py[] = { 0.5, 0.5, 0.0, 0.5 };

	// test
	PetscInt in[3];

	in_polygon(nP, px, py, nC, cx, cy, box, 1e-10, in);

	for(i = 0; i < nP; i++)
	{
		if(in[i]) printf("Point %d [x=%f, y=%f] is inside polygon \n", i, px[i], py[i]);
		else      printf("Point %d [x=%f, y=%f] is outside polygon \n", i, px[i], py[i]);
	}

*/
