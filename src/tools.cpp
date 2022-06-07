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
#include <unistd.h>

//---------------------------------------------------------------------------
// Printing functions
//---------------------------------------------------------------------------
void PrintStart(PetscLogDouble *t_beg, const char *msg, const char *filename)
{
	PetscTime(t_beg);

	if(filename)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s file(s) <%s> ... ", msg, filename);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s ... ", msg);
	}
}
//---------------------------------------------------------------------------
void PrintDone(PetscLogDouble t_beg)
{
	PetscLogDouble t_end;

	MPI_Barrier(PETSC_COMM_WORLD);

	PetscTime(&t_end);

	PetscPrintf(PETSC_COMM_WORLD,"done (%g sec)\n", t_end - t_beg);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
}
//---------------------------------------------------------------------------
void PrintStep(PetscInt step)
{
	char *number, *p;

	char line[] = "==========================================================================";
	char left[]  = " STEP ";
	char right[] = " ";
	asprintf(&number, "%d", step);

	p = line + (strlen(line) - strlen(left) - strlen(number) - strlen(right))/2;

	memcpy(p, left,   strlen(left));   p += strlen(left);
	memcpy(p, number, strlen(number)); p += strlen(number);
	memcpy(p, right,  strlen(right));

	free(number);

	PetscPrintf(PETSC_COMM_WORLD,"%s\n", line);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecReadRestart"
PetscErrorCode VecReadRestart(Vec x, FILE *fp)
{
	PetscInt     size;
	PetscScalar *xarr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecGetLocalSize(x, &size); CHKERRQ(ierr);

	// get vector array
	ierr = VecGetArray(x, &xarr); CHKERRQ(ierr);

	// write to file
	fread(xarr, sizeof(PetscScalar), (size_t)size, fp);

	// restore vector array
	ierr = VecRestoreArray(x, &xarr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecWriteRestart"
PetscErrorCode VecWriteRestart(Vec x, FILE *fp)
{
	PetscInt     size;
	PetscScalar *xarr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecGetLocalSize(x, &size); CHKERRQ(ierr);

	// get vector array
	ierr = VecGetArray(x, &xarr); CHKERRQ(ierr);

	// write to file
	fwrite(xarr, sizeof(PetscScalar), (size_t)size, fp);

	// restore vector array
	ierr = VecRestoreArray(x, &xarr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
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
#define __FUNCT__ "clearIntArray"
PetscErrorCode clearIntArray(PetscInt *arr, const PetscInt n)
{
	size_t          sz;
	PetscErrorCode 	ierr;

	PetscFunctionBegin;

	// compute size in bytes
	sz = (size_t)n*sizeof(PetscInt);

	// clear memory
	ierr = PetscMemzero(arr, sz); CHKERRQ(ierr);

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
PetscInt ISRankZero(MPI_Comm comm)
{
	PetscMPIInt rank;

	MPI_Comm_rank(comm, &rank);

	return (rank == 0);
}
//---------------------------------------------------------------------------
PetscInt ISParallel(MPI_Comm comm)
{
	PetscMPIInt size;

	MPI_Comm_size(comm, &size);

	return (size > 1);
}
//---------------------------------------------------------------------------
PetscMPIInt getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p)
{
	// get global rank of processor in DMDA

	if(i < 0 || i >= m || j < 0 || j >= n || k < 0 || k >= p) return -1;

	return (PetscMPIInt)(i + j*m + k*m*n);
}
//---------------------------------------------------------------------------
PetscMPIInt getGlobalRankPeriodic(
		PetscInt i,  PetscInt j,  PetscInt k,
		PetscInt m,  PetscInt n,  PetscInt p,
		PetscInt pi, PetscInt pj, PetscInt pk)
{
	// get global rank of processor in DMDA
	if(pi) { if(i < 0) { i = m-1; } if(i >= m) {i = 0; } }
	if(pj) { if(j < 0) { j = n-1; } if(j >= n) {j = 0; } }
	if(pk) { if(k < 0) { k = p-1; } if(k >= p) {k = 0; } }

	if(i < 0 || i >= m || j < 0 || j >= n || k < 0 || k >= p) return -1;

	return (PetscMPIInt)(i + j*m + k*m*n);
}
//---------------------------------------------------------------------------
void getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n)
{
	// get local ranks of processor in DMDA

	(*k) =  rank/(m*n);
	(*j) = (rank - (*k)*m*n)/m;
	(*i) =  rank - (*k)*m*n - (*j)*m;

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DirMake"
PetscErrorCode DirMake(const char *name)
{
	int status;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create a new directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{

#ifdef _WIN32
		// call this on windows machines
		status = mkdir(name);
#else
		// standard access pattern drwxr-xr-x
		status = mkdir(name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif 

		if(status && errno != EEXIST)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to create directory %s", name);
		}
	}

	// synchronize
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DirRemove"
PetscErrorCode DirRemove(const char *name)
{
	int status;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// synchronize
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	// remove directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		status = rmdir(name);

		if(status)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to remove directory %s", name);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DirRename"
PetscErrorCode DirRename(const char *old_name, const char *new_name)
{
	int status;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// synchronize
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	// rename directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		status = rename(old_name, new_name);

		if(status)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to rename directory %s into %s", old_name, new_name);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DirCheck"
PetscErrorCode DirCheck(const char *name, PetscInt *exists)
{
	struct stat s;
	int         status;
	PetscInt    check;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		status = stat(name, &s);

		// check whether file exists and is a directory
		check = (!status && S_ISDIR(s.st_mode));
	}

	// synchronize
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Bcast(&check, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}

	(*exists) = check;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
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
				if(xp < PetscMin(ax, bx) || PetscMax(ax, bx) < xp) continue;

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

						if(PetscMin(bx, xvind) < xp && xp < PetscMax(bx, xvind))
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
		in[ip]   = PetscMax(point_on, point_in);
	}
}
//---------------------------------------------------------------------------
void linSpace(
	PetscScalar  min,
	PetscScalar  max,
	PetscInt     N,
	PetscScalar *outVec)
{
    PetscScalar delta = (max-min)/((PetscScalar)(N-1));

    for(PetscInt i = 0; i < N; i++)
    {
        outVec[i] = min + ((PetscScalar)i)*delta;
    }
}
//---------------------------------------------------------------------------
void interpStretch(
	PetscScalar *Sx,
    PetscScalar *Sy,
    PetscInt     numCtrlPoly,
    PetscInt    *CtrlPoly,
    PetscInt     numPoly,
    PetscScalar *SxAll,
    PetscScalar *SyAll)
{
    PetscInt     i,j,k;
	PetscInt     maxSize = 0;
    
    // enter first value
    SxAll[CtrlPoly[0]] = Sx[0];
    SyAll[CtrlPoly[0]] = Sy[0];

    // get maximum size for SxNew and SyNew
    for (i=1; i < numCtrlPoly; i++)
    {  
        if (CtrlPoly[i]-CtrlPoly[i-1]+1 > maxSize)
        {
            maxSize = CtrlPoly[i]-CtrlPoly[i-1]+1;
        }
    }

    // interpolate the rest
    for (i=1; i < numCtrlPoly; i++)
    {
        PetscScalar SxNew[maxSize];
        PetscScalar SyNew[maxSize];
        linSpace(Sx[i-1],Sx[i],CtrlPoly[i] - CtrlPoly[i-1] + 1, SxNew);
        linSpace(Sy[i-1],Sy[i],CtrlPoly[i] - CtrlPoly[i-1] + 1, SyNew);
        for (j = CtrlPoly[i-1], k = 0; j < CtrlPoly[i]+1; j++, k++)
        {
            SxAll[j] = SxNew[k];
            SyAll[j] = SyNew[k];
        }
    }

    // extend stretch to first and last polygon
    if (CtrlPoly[0] != 0)
    {
        for (i=0; i < CtrlPoly[0]; i++)
        {
            SxAll[i] = Sx[0];
            SyAll[i] = Sy[0];
        }
    }
    if (CtrlPoly[numCtrlPoly-1] != numCtrlPoly-1)
    {
        for (i=CtrlPoly[numCtrlPoly-1]; i < numPoly; i++)
        {
            SxAll[i] = Sx[numCtrlPoly-1];
            SyAll[i] = Sy[numCtrlPoly-1];
        }
    }

    // no stretchfactor should be smaller than 0.1
    for (i=0; i < numPoly; i++)
    {
    	if (SxAll[i] < 0.1) {SxAll[i] = 0.1;}
    	if (SyAll[i] < 0.1) {SyAll[i] = 0.1;}
    }
}
//---------------------------------------------------------------------------
void findCenterMass(
	PetscScalar *coords,
	PetscInt     nN,
	PetscScalar &x_cen,
	PetscScalar &y_cen)
{
    PetscInt    i;
    PetscScalar sumX = 0;
    PetscScalar sumY = 0;
    PetscScalar A = 0;
    PetscScalar xc = 0; 
    PetscScalar yc = 0;
    PetscScalar meanX, meanY;
    PetscScalar xp[nN], yp[nN], a[nN];

	// split coords into x and y
	PetscScalar x[nN];
	PetscScalar y[nN];
	for (i=0; i < nN; i++)
	{
		x[i] = coords[2*i];
		y[i] = coords[2*i+1];
	}

    // get mean of all nodes
    for (i=0; i < nN; i++)
    {
        sumX += x[i];
        sumY += y[i];
    }
    meanX = sumX/nN;
    meanY = sumY/nN;

    // shift nodes by mean for better precision
    for (i=0; i < nN; i++)
    {
        x[i] -= meanX;
        y[i] -= meanY;
    }

    // summations for CCW boundary
    for (i=0; i < nN; i++)
    {
        if (i == nN-1) 
        {
            xp[i] = x[0];
            yp[i] = y[0];
        }
        else
        {
            xp[i] = x[i+1];
            yp[i] = y[i+1];
        }
        a[i]   =  x[i]*yp[i] - xp[i]*y[i];
        xc    += (x[i]+xp[i])*a[i];
        yc    += (y[i]+yp[i])*a[i];
        A     +=  a[i];
    }
    A = A/2;
    xc = xc/(6*A);
    yc = yc/(6*A);

    // shift back to original coords
    x_cen = xc + meanX;
    y_cen = yc + meanY;
}
//---------------------------------------------------------------------------
// stretch Polygon
void stretchPolygon(
	PetscScalar *coords,
	PetscInt nN,
	PetscScalar Sx,
	PetscScalar Sy)
{
    PetscScalar x_cen, y_cen;    
    PetscInt    iNode;
    PetscScalar x,y;

	// find center of mass
    findCenterMass(coords,nN,x_cen,y_cen);

    // loop over nodes
    for (iNode = 0; iNode < nN; iNode++)
    {
        // transform to relative coordinates
        x        = coords[2*iNode]   - x_cen;
        y        = coords[2*iNode+1] - y_cen;

        // stretch
        x        = x * Sx;
        y        = y * Sy;

        // transform back to absolute corrdinates
        coords[2*iNode]   = x + x_cen;
        coords[2*iNode+1] = y + y_cen;
    }
}
//---------------------------------------------------------------------------
// indexing functions
//---------------------------------------------------------------------------
PetscInt getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[])
{
	// compute pointers from counts, return total count
	// if counts and pointers are stored in the same array, counts are overwritten

	PetscInt i, t, tcnt = 0;

	for(i = 0; i < n; i++)
	{
		t      = counts[i];
		ptr[i] = tcnt;
		tcnt  += t;
	}
	return tcnt;
}
//---------------------------------------------------------------------------
void rewindPtr(PetscInt n, PetscInt ptr[])
{
	// rewind pointers after using them as access iterators

	PetscInt i, prev = 0, next;

	for(i = 0; i < n; i++)
	{
		next   = ptr[i];
		ptr[i] = prev;
		prev   = next;
	}
}
//---------------------------------------------------------------------------
/*
	PetscInt arr[18];
	PetscInt counts[] = {3, 5, 7, 2, 1, 0};
	PetscInt tnum = getPtrCnt(5, counts, counts);
	counts[5] = tnum;
	for(int i = 0; i < 5; i++) { for(int j = counts[i]; j < counts[i+1]; j++) arr[j] = i; }
	for(int i = 0; i < 18; i++) counts[arr[i]]++;
	rewindPtr(5, counts);
*/
//---------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getPhaseRatio"
PetscErrorCode getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum)
{
	// compute phase ratio array

	PetscInt    i;
	PetscScalar sum = 0.0;

	PetscFunctionBegin;

	for(i = 0; i < n; i++) sum  += v[i];

	if(sum == 0.0)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, " Empty control volume");
	}

	for(i = 0; i < n; i++) v[i] /= sum;

	(*rsum) = sum;

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
// bisection algorithm for scalar nonlinear equation
PetscInt solveBisect(
		PetscScalar a,
		PetscScalar b,
		PetscScalar tol,
		PetscScalar maxit,
		PetscScalar &x,
		PetscInt    &it,
		PetscScalar (*f) (PetscScalar x, void *pctx),
		void *pctx)
{
	PetscScalar fa, fx;

	// initialize
	x  = a;
	it = 1;

	// get residual of left bound (initial guess)
	fa = f(a, pctx);

	// check whether closed-form solution exists
    if(PetscAbsScalar(fa) <= tol)
	{
    	// return convergence flag
    	return 1;
    }

	do
	{	// get new iterate
	    x = (a + b)/2.0;

	    // get new residual
	    fx = f(x, pctx);

	    // update interval
	    if(fa*fx < 0.0)
	    {
	    	b = x;
	    }
	    else
	    {
	    	a = x; fa = fx;
	    }

	    // update iteration count
	    it++;

	} while(PetscAbsScalar(fx) > tol && it < maxit);

	// return convergence flag
	return PetscAbsScalar(fx) <= tol;
}
//-----------------------------------------------------------------------------
