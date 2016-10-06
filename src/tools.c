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
	ierr = PetscOptionsGetRealArray(NULL, NULL, ident,  a, &nmax, &found); CHKERRQ(ierr);

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
	ierr = PetscOptionsGetIntArray(NULL, NULL, ident,  a, &nmax, &found); CHKERRQ(ierr);

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
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
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
// Fast detection points inside a polygonal region.
//
// Originally written as a MATLAB mexFunction by:
// A. David Redish      email: adr@nsma.arizona.edu
// Guillaume Jacquenot  email: guillaume.jacquenot@gmail.com
//
// Modified to be callable directly from C by:
// Anton A. Popov       email: popov@uni-mainz.de
//
// Output is simplified to skip distinguishing points IN and ON the polygon.
// Coordinates storage format is changed from strided to interleaved.
// Bounding box computation is separated from the main test function.
//
//---------------------------------------------------------------------------
void polygon_box(
	PetscInt    *pnv,    // number of polygon vertices (can be modified)
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar  rtol,   // relative tolerance
	PetscScalar *atol,   // absolute tolerance
	PetscScalar *box)    // bounding box of a polygon
{
	PetscInt    iv, nv;
	PetscScalar ax, bx, ay, by, xv, yv;
	PetscScalar xmin, xmax, ymin, ymax, dmin;

	nv = (*pnv);
	ax = vcoord[0];
	ay = vcoord[1];
	bx = vcoord[2*(nv-1)  ];
	by = vcoord[2*(nv-1)+1];

	// decrease number of points if the polygon is closed
	if(ax == bx && ay == by) nv--;

	// calculate bounding box of a polygon
	xmin = xmax = vcoord[0];
	ymin = ymax = vcoord[1];

	for(iv = 0; iv < nv; iv++)
	{
		// get vertex coordinates
		xv = vcoord[2*iv  ];
		yv = vcoord[2*iv+1];

		if(xv < xmin) xmin = xv;
		if(xv > xmax) xmax = xv;
		if(yv < ymin) ymin = yv;
		if(yv > ymax) ymax = yv;
	}

	box[0] = xmin;
	box[1] = xmax;
	box[2] = ymin;
	box[3] = ymax;

	// get smallest extent of the polygon
	dmin = xmax-xmin; if(ymax-ymin < dmin) dmin = ymax-ymin;

	// store parameters
	(*atol) = rtol*dmin; // absolute tolerance
	(*pnv)  = nv;        // number of vertices
}
//---------------------------------------------------------------------------
void in_polygon(
	PetscInt     np,     // number of test points
	PetscScalar *pcoord, // coordinates of test points
	PetscInt     nv,     // number of polygon vertices
	PetscScalar *vcoord, // coordinates of polygon vertices
	PetscScalar *box,    // bounding box of a polygon (optimization)
	PetscScalar  atol,   // absolute tolerance
	PetscInt    *in)     // point location flags (1-inside, 0-outside)
{
	PetscInt    ip, iv, ind;
	PetscInt    point_on, point_in;
	PetscScalar ax, bx, ay, by;
	PetscScalar nIntersect, intersecty, tmp;
	PetscScalar xmin, xmax, ymin, ymax, xp, yp, xvind;

	// get bounding box
	xmin = box[0];
	xmax = box[1];
	ymin = box[2];
	ymax = box[3];

	// test whether each point is in polygon
	for(ip = 0; ip < np; ip++)
	{
		// assume point is outside
		in[ip] = 0;

		// get point coordinates
		xp = pcoord[2*ip    ];
		yp = pcoord[2*ip + 1];

		// check bounding box
		if(xp < xmin) continue;
		if(xp > xmax) continue;
		if(yp < ymin) continue;
		if(yp > ymax) continue;

		// count the number of intersections
		nIntersect = 0.0;
		point_on   = 0;

		for(iv = 0; iv < nv; iv++)
		{
			// does the line PQ intersect the line AB?
			if(iv == nv-1)
			{
				ax = vcoord[2*(nv-1)  ];
				ay = vcoord[2*(nv-1)+1];
				bx = vcoord[0         ];
				by = vcoord[1         ];
			}
			else
			{
				ax = vcoord[2*iv      ];
				ay = vcoord[2*iv+1    ];
				bx = vcoord[2*(iv+1)  ];
				by = vcoord[2*(iv+1)+1];
			}

			if(ax == bx)
			{
				// vertical points
				if(xp == ax)
				{
					// ensure order correct
					if(ay > by)
					{
						tmp = ay; ay = by; by = tmp;
					}
					if(yp >= ay && yp <= by)
					{
						point_on   = 1;
						nIntersect = 0.0;
						break;
					}
				}
			}
			else
			{
				// non-vertical points
				if(xp < MIN(ax, bx) || MAX(ax, bx) < xp) continue;

				intersecty = ay + (xp - ax)/(bx - ax)*(by - ay);

				if(fabs(intersecty - yp) < atol)
				{
					point_on   = 1;
					nIntersect = 0.0;
					break;
				}
				else if(intersecty < yp && (ax == xp || bx == xp))
				{
					if(ax == xp)
					{
						if(iv == 0)
						{
							ind = nv-1;
						}
						else
						{
							ind = iv-1;
						}

						xvind = vcoord[2*ind];

						if(MIN(bx, xvind) < xp && xp < MAX(bx, xvind))
						{
							nIntersect += 1.0;
						}
					}
				}
				else if (intersecty < yp)
				{
					nIntersect += 1.0;
				}
			}
		}

		// check if the contour polygon is closed
		point_in = (PetscInt)(nIntersect - 2.0*floor(nIntersect/2.0));
		in[ip]   = MAX(point_on, point_in);
	}
}
//---------------------------------------------------------------------------
