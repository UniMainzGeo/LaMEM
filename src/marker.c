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
 **    filename:   marker.c
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
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "marker.h"
#include "parsing.h"
#include "advect.h"
#include "fdstag.h"
#include "scaling.h"
#include "JacRes.h"
#include "phase.h"
#include "tools.h"

/*
#START_DOC#
#END_DOC#
*/
//---------------------------------------------------------------------------
// initialization function header
#define INIT_FUNCTION_HEADER \
	FDSTAG      *fs; \
	Scaling     *scal; \
	PetscInt     imark; \
	PetscScalar  lx, ly, lz, bx, by, bz, ex, ey, ez, cf, *X; \
	PetscErrorCode ierr; \
	PetscFunctionBegin; \
	fs   = advect->fs; \
	scal = advect->scal; \
	cf   = scal->length; \
	ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr); \
	lx = ex - bx; \
	ly = ey - by; \
	lz = ez - bz;
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInit"
PetscErrorCode ADVMarkInit(AdvCtx *actx, FB *fb)
{
	FDSTAG   *fs;
	PetscInt  nmarkx, nmarky, nmarkz, nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	PetscPrintf(PETSC_COMM_WORLD," Starting marker initialization routine\n");

	// allocate storage for uniform distribution
	if(actx->msetup != _FILES_)
	{
		// get local number of markers
		nmarkx  = fs->dsx.ncels * actx->NumPartX;
		nmarky  = fs->dsy.ncels * actx->NumPartY;
		nmarkz  = fs->dsz.ncels * actx->NumPartZ;
		nummark = nmarkx*nmarky*nmarkz;

		// allocate storage
		ierr = ADVReAllocStorage(actx, nummark); CHKERRQ(ierr);

		// set number of markers
		actx->nummark = nummark;
	}

	// initialize coordinates, add random noise
	if(actx->msetup != _FILES_
	&& actx->msetup != _POLYGONS_)
	{
		ierr = ADVMarkInitCoord(actx); CHKERRQ(ierr);
	}

	// display info
	PetscPrintf(PETSC_COMM_WORLD," Marker setup employed [msetup] : ");

	// initialize marker phase, temperature, etc.
	if     (actx->msetup == _GEOM_)       { PetscPrintf(PETSC_COMM_WORLD,"%s\n","geom");     ierr = ADVMarkInitGeom    (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _FILES_)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","files");    ierr = ADVMarkInitFiles   (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _POLYGONS_)   { PetscPrintf(PETSC_COMM_WORLD,"%s\n","polygons"); ierr = ADVMarkInitPolygons(actx, fb); CHKERRQ(ierr); }

	// optionally read temperature from file
	if(actx->msetup != _FILES_)
	{
		ierr = ADVMarkSetTempFromFile(actx, fb); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD," Finished marker initialization routine\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitCoord"
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx)
{
	// initializes coordinates and adds random noise if required for hard-coded setups

	FDSTAG      *fs;
	PetscScalar  x, y, z, dx, dy, dz;
	PetscInt     i, j, k, ii, jj, kk;
	PetscInt     imark;
	PetscRandom  rctx;
	PetscScalar  cf_rand;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	if(actx->randNoise)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Adding random noise to marker distribution \n");

		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);
	}

	// marker counter
	imark = 0;

	// create uniform distribution of markers/cell for variable grid
	for(k = 0; k < fs->dsz.ncels; k++)
	{
		// spacing of particles
		dz = (fs->dsz.ncoor[k+1]-fs->dsz.ncoor[k])/(PetscScalar)actx->NumPartZ;
		for(j = 0; j < fs->dsy.ncels; j++)
		{
			// spacing of particles
			dy = (fs->dsy.ncoor[j+1]-fs->dsy.ncoor[j])/(PetscScalar)actx->NumPartY;
			for(i = 0; i < fs->dsx.ncels; i++)
			{
				// spacing of particles
				dx = (fs->dsx.ncoor[i+1]-fs->dsx.ncoor[i])/(PetscScalar)actx->NumPartX;

				// loop over markers in cells
				for (kk = 0; kk < actx->NumPartZ; kk++)
				{
					if(kk == 0) z = fs->dsz.ncoor[k] + dz*0.5;
					else        z = fs->dsz.ncoor[k] + dz*0.5 + (PetscScalar)kk*dz;

					for (jj = 0; jj < actx->NumPartY; jj++)
					{
						if(jj == 0) y = fs->dsy.ncoor[j] + dy*0.5;
						else        y = fs->dsy.ncoor[j] + dy*0.5 + (PetscScalar)jj*dy;

						for(ii = 0; ii < actx->NumPartX; ii++)
						{
							if(ii == 0) x = fs->dsx.ncoor[i] + dx*0.5;
							else        x = fs->dsx.ncoor[i] + dx*0.5 + (PetscScalar)ii*dx;

							// set marker coordinates
							actx->markers[imark].X[0] = x;
							actx->markers[imark].X[1] = y;
							actx->markers[imark].X[2] = z;

							if(actx->randNoise)
							{
								// add random noise
								// decrease/increase amount of noise by changing A in: (cf_rand-0.5)*dx/A
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[0] += (cf_rand - 0.5)*dx/1;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[1] += (cf_rand - 0.5)*dy/1;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[2] += (cf_rand - 0.5)*dz/1;
							}

							// increment local counter
							imark++;
						}
					}
				}
			}
		}
	}

	// destroy random context
	if(actx->randNoise)
	{
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSave"
PetscErrorCode ADVMarkSave(AdvCtx *actx)
{
	int          fd;
	PetscInt     imark;
	Marker      *P;
	PetscViewer  view_out;
	char        *filename, path[_STR_LEN_];
	PetscScalar *markbuf, *markptr, header, chLen, chTemp, Tshift, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!actx->saveMark) PetscFunctionReturn(0);

	// access context
	chLen  = actx->jr->scal->length;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;

	PetscPrintf(PETSC_COMM_WORLD,"Saving markers in parallel to files: %s.xxxxxxxx.dat \n", actx->saveFile);

	// extract directory path
	strcpy(path, actx->saveFile); (*strrchr(path, '/')) = '\0';

	// create directory
	ierr = DirMake(path); CHKERRQ(ierr);

	// compile file name
	asprintf(&filename, "%s.%1.8lld.dat", actx->saveFile, (LLD)actx->iproc);

	// open file for binary output
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &view_out); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_out, &fd);                                CHKERRQ(ierr);

	// initialize file header for MATLAB compatibility
	header = -1;

	// create write buffer
	ierr = PetscMalloc((size_t)(5*actx->nummark)*sizeof(PetscScalar), &markbuf); CHKERRQ(ierr);

	// copy data from storage into buffer
	for(imark = 0, markptr = markbuf; imark < actx->nummark; imark++, markptr += 5)
	{
		P          =              &actx->markers[imark];
		markptr[0] =              P->X[0]*chLen;
		markptr[1] =              P->X[1]*chLen;
		markptr[2] =              P->X[2]*chLen;
		markptr[3] = (PetscScalar)P->phase;
		markptr[4] =              P->T*chTemp - Tshift;
	}

	// write binary output
	s_nummark = (PetscScalar)actx->nummark;
	ierr = PetscBinaryWrite(fd, &header,    1,               PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fd, &s_nummark, 1,               PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fd, markbuf,    5*actx->nummark, PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);

	// destroy file handle & file name
	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
	free(filename);

	// destroy buffer
	ierr = PetscFree(markbuf); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished saving markers in parallel \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkCheckMarkers"
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx)
{
	// check initial marker distribution
	FDSTAG      *fs;
	PetscScalar *X;
	PetscInt     error;
	PetscScalar  bx, by, bz;
	PetscScalar  ex, ey, ez;
	PetscInt     *numMarkCell, rbuf[4], sbuf[4];
	PetscInt     i, maxid, NumInvalidPhase, numNonLocal, numEmpty, numWrong, refMarkCell;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// get maximum phase ID
	maxid = actx->dbm->numPhases - 1;

	// get reference number of markers per cell
	refMarkCell = actx->NumPartX*actx->NumPartY*actx->NumPartZ;

	// get local coordinate bounds
	ierr = FDSTAGGetLocalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	// allocate marker counter array
	ierr = makeIntArray(&numMarkCell, NULL, fs->nCells); CHKERRQ(ierr);

	// clear error flag
	error = 0;

	// count markers with invalid phase ID & non-local markers
	NumInvalidPhase = 0;
	numNonLocal     = 0;

	for(i = 0; i < actx->nummark; i++)
	{
		// marker should have a valid phase ID
		if(actx->markers[i].phase > maxid) NumInvalidPhase++;

		// get marker coordinates
		X = actx->markers[i].X;

		// marker must be local (check bounding box)
		if(X[0] < bx || X[0] > ex
		|| X[1] < by || X[1] > ey
		|| X[2] < bz || X[2] > ez) numNonLocal++;

		// count number of markers in the cells
		numMarkCell[actx->cellnum[i]]++;
	}

	// count empty & sparse cells
	numEmpty = 0;
	numWrong = 0;

	for(i = 0; i < fs->nCells; i++)
	{
		if(numMarkCell[i] == 0)           numEmpty++;
		if(numMarkCell[i] != refMarkCell) numWrong++;
	}

	// clear
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	// get global figures
	if(actx->nproc != 1)
	{
		sbuf[0] = NumInvalidPhase;
		sbuf[1] = numNonLocal;
		sbuf[2] = numEmpty;
		sbuf[3] = numWrong;

		ierr = MPI_Allreduce(sbuf, rbuf, 4, MPIU_INT, MPI_SUM, actx->icomm); CHKERRQ(ierr);

		NumInvalidPhase = rbuf[0];
		numNonLocal     = rbuf[1];
		numEmpty        = rbuf[2];
		numWrong        = rbuf[3];
	}

	// print diagnostics
	if(NumInvalidPhase)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of markers with invalid phase ID: %lld\n", (LLD)NumInvalidPhase); CHKERRQ(ierr);
		error = 1;
	}

	if(numNonLocal)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of non-local markers: %lld\n", (LLD)numNonLocal); CHKERRQ(ierr);
		error = 1;
	}

	if(numEmpty)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of exactly empty cells: %lld\n", (LLD)numEmpty); CHKERRQ(ierr);
		error = 1;
	}

	if(numWrong)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of cells with incorrect number of markers (nmark_x*nmark_y*nmark_z): %lld\n", (LLD)numWrong); CHKERRQ(ierr);
		error = 1;
	}

	if(error)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Problems with initial marker distribution (see the above message)");
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSetTempFromFile"
PetscErrorCode ADVMarkSetTempFromFile(AdvCtx *actx, FB *fb)
{
	FDSTAG       *fs;
	int           fd;
	Marker       *P;
	PetscViewer   view_in;
	char          filename[_STR_LEN_];
	PetscScalar   header[2], dim[3];
	PetscInt      Fsize, imark, nummark, nmarkx, nmarky, nmarkz;
	PetscScalar   DX, DY, DZ, bx, by, bz, ex, ey, ez;
	PetscScalar   xp, yp, zp, Xc, Yc, Zc, xpL, ypL, zpL;
	PetscScalar  *Temp;
	PetscInt      Ix, Iy, Iz;
	PetscScalar   chTemp, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "temp_file", filename, NULL); CHKERRQ(ierr);

	// check whether file is provided
	if(!strlen(filename)) PetscFunctionReturn(0);

	// access context
	fs     = actx->fs;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;

	PetscPrintf(PETSC_COMM_WORLD," Loading temperature redundantly from file: %s \n", filename);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header & size of file
	ierr = PetscBinaryRead(fd, &header, 2, PETSC_SCALAR); CHKERRQ(ierr);
	Fsize = (PetscInt)(header[1])-3;

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)Fsize*sizeof(PetscScalar), &Temp); CHKERRQ(ierr);

	// read entire file
	ierr = PetscBinaryRead(fd, &dim, 3,     PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, Temp, Fsize, PETSC_SCALAR); CHKERRQ(ierr);

	// get mesh extents
	ierr = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

	// get grid spacing
	DX = (ex - bx)/(dim[0] - 1.0);
	DY = (ey - by)/(dim[1] - 1.0);
	DZ = (ez - bz)/(dim[2] - 1.0);

	// get local number of markers
	nmarkx  = fs->dsx.ncels * actx->NumPartX;
	nmarky  = fs->dsy.ncels * actx->NumPartY;
	nmarkz  = fs->dsz.ncels * actx->NumPartZ;
	nummark = nmarkx*nmarky*nmarkz;

	PetscInt nx, ny;
	nx = (PetscInt)dim[0];
	ny = (PetscInt)dim[1];

	for(imark = 0; imark < nummark; imark++)
	{
		// get current marker
		P = &actx->markers[imark];

		// get global marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// index of the lower left corner of the element (of the temperature grid) in which the particle is
		Ix = (PetscInt)floor((xp - bx)/DX);
		Iy = (PetscInt)floor((yp - by)/DY);
		Iz = (PetscInt)floor((zp - bz)/DZ);

		// coordinate of the first corner (lower left deepest)
		Xc = bx + (PetscScalar)Ix*DX;
		Yc = by + (PetscScalar)Iy*DY;
		Zc = bz + (PetscScalar)Iz*DZ;

		// Local coordinate of the particle inside a temperature element
		xpL = (xp - Xc)/DX;
		ypL = (yp - Yc)/DY;
		zpL = (zp - Zc)/DZ;

		// Interpolate value on the particle using trilinear shape functions
		P->T = ((
		(1.0-xpL) * (1.0-ypL) * (1.0-zpL) * Temp[Iz    *nx*ny + Iy     * nx + Ix   ] +
		 xpL      * (1.0-ypL) * (1.0-zpL) * Temp[Iz    *nx*ny + Iy     * nx + Ix+1 ] +
		 xpL      *  ypL      * (1.0-zpL) * Temp[Iz    *nx*ny + (Iy+1) * nx + Ix+1 ] +
		(1.0-xpL) *  ypL      * (1.0-zpL) * Temp[Iz    *nx*ny + (Iy+1) * nx + Ix   ] +
		(1.0-xpL) * (1.0-ypL) *  zpL      * Temp[(Iz+1)*nx*ny + Iy     * nx + Ix   ] +
		 xpL      * (1.0-ypL) *  zpL      * Temp[(Iz+1)*nx*ny + Iy     * nx + Ix+1 ] +
		 xpL      *  ypL      *  zpL      * Temp[(Iz+1)*nx*ny + (Iy+1) * nx + Ix+1 ] +
		(1.0-xpL) *  ypL      *  zpL      * Temp[(Iz+1)*nx*ny + (Iy+1) * nx + Ix   ] ) + Tshift)/chTemp;
	}

	// clear memory
	PetscFree(Temp);
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
// Specific initialization routines
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFiles"
PetscErrorCode ADVMarkInitFiles(AdvCtx *actx, FB *fb)
{
	int          fd;
	Marker      *P;
	PetscViewer  view_in;
	char        *filename, file[_STR_LEN_];
	PetscScalar *markbuf, *markptr, header, chTemp, chLen, Tshift, s_nummark;
	PetscInt     imark, nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "mark_load_file", file, "./markers/mdb"); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"Loading markers in parallel from files: %s.xxxxxxxx.dat \n", file);

	// compile input file name with extension
	asprintf(&filename, "%s.%1.8lld.dat", file, (LLD)actx->iproc);

	// open file
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd);                               CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header
	ierr = PetscBinaryRead(fd, &header, 1, PETSC_SCALAR); CHKERRQ(ierr);

	// read number of local of markers
	ierr = PetscBinaryRead(fd, &s_nummark, 1, PETSC_SCALAR); CHKERRQ(ierr);
	nummark = (PetscInt)s_nummark;

	// allocate marker storage
	ierr = ADVReAllocStorage(actx, nummark); CHKERRQ(ierr);

	// set number of markers
	actx->nummark = nummark;

	// allocate marker buffer
	ierr = PetscMalloc((size_t)(5*actx->nummark)*sizeof(PetscScalar), &markbuf); CHKERRQ(ierr);

	// read markers into buffer
	ierr = PetscBinaryRead(fd, markbuf, 5*actx->nummark, PETSC_SCALAR); CHKERRQ(ierr);

	// destroy file handle & file name
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
	free(filename);

	// get characteristic length & temperature
	chLen  = actx->jr->scal->length;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;

	// copy buffer to marker storage
	for(imark = 0, markptr = markbuf; imark < actx->nummark; imark++, markptr += 5)
	{
		P        =           &actx->markers[imark];
		P->X[0]  =           markptr[0]/chLen;
		P->X[1]  =           markptr[1]/chLen;
		P->X[2]  =           markptr[2]/chLen;
		P->phase = (PetscInt)markptr[3];
		P->T     =          (markptr[4] + Tshift)/chTemp;
	}

	// free marker buffer
	ierr = PetscFree(markbuf); CHKERRQ(ierr);

	// wait until all processors finished reading markers
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished Loading markers in parallel \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitGeom"
PetscErrorCode ADVMarkInitGeom(AdvCtx *actx, FB *fb)
{
	Marker      *P;
	PetscScalar chLen;
	PetscInt    jj, iter, ngeom, imark, maxPhaseID;
	GeomPrim    geom[_max_geom_], *sphere, *box, *hex, *layer;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	iter       = 0;
	ngeom      = 0;
	maxPhaseID = actx->dbm->numPhases - 1;
	chLen      = actx->jr->scal->length;

	// clear storage
	ierr = PetscMemzero(geom, sizeof(GeomPrim)*(size_t)_max_geom_); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"Reading geometric primitives: \n\n");

	//========
	// SPHERES
	//========

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<SphereStart>", "<SphereEnd>"); CHKERRQ(ierr);

	ngeom += fb->nblocks;

	if(ngeom > _max_geom_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many geometric primitives! Max allowed: %lld", (LLD)_max_geom_);
	}

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		sphere = &geom[iter++];

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &sphere->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "center",  sphere->center, 3, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "radius", &sphere->radius, 1, chLen);      CHKERRQ(ierr);

		sphere->setPhase = setPhaseSphere;

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//======
	// BOXES
	//======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<BoxStart>", "<BoxEnd>"); CHKERRQ(ierr);

	ngeom += fb->nblocks;

	if(ngeom > _max_geom_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many geometric primitives! Max allowed: %lld", (LLD)_max_geom_);
	}

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		box = &geom[iter++];

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &box->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bounds",  box->bounds, 6, chLen);      CHKERRQ(ierr);

		box->setPhase = setPhaseBox;

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//======
	// HEXES
	//======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<HexStart>", "<HexEnd>"); CHKERRQ(ierr);

	ngeom += fb->nblocks;

	if(ngeom > _max_geom_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many geometric primitives! Max allowed: %lld", (LLD)_max_geom_);
	}

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		hex = &geom[iter++];

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &hex->phase, 1,  maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "coord",   hex->coord, 24, chLen);      CHKERRQ(ierr);

		hex->setPhase = setPhaseHex;

		// compute bounding box
		HexGetBoundingBox(hex->coord, hex->bounds);

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);


	//=======
	// LAYERS
	//=======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<LayerStart>", "<LayerEnd>"); CHKERRQ(ierr);

	ngeom += fb->nblocks;

	if(ngeom > _max_geom_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many geometric primitives! Max allowed: %lld", (LLD)_max_geom_);
	}

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		layer = &geom[iter++];

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &layer->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "top",    &layer->top,    1, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bottom", &layer->bot,    1, chLen);      CHKERRQ(ierr);

		layer->setPhase = setPhaseLayer;

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//==============
	// ASSIGN PHASES
	//==============

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		//set default
		P->phase = actx->bgPhase;

		// override from geometric primitives
		for(jj = 0; jj < ngeom; jj++)
		{
			if(geom[jj].setPhase(&geom[jj], P)) break;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitPolygons"
PetscErrorCode ADVMarkInitPolygons(AdvCtx *actx, FB *fb)
{
	// REDUNDANTLY loads a file with 2D-polygons that coincide with the marker planes
	// each processor uses the full polygonal shapes to find assign phase ids to local markers

	FDSTAG        *fs;
	int            fd;
	PetscViewer    view_in;
	char           filename[_STR_LEN_];
	PetscScalar    header[2];
	PetscInt       tstart[3], tend[3], nmark[3], nidx[3], nidxmax;
	PetscInt       k, n, kvol, Fcount, Fsize, VolN, Nmax, Lmax, kpoly;
	Polygon2D      Poly;
	PetscInt      *polyin;
	PetscInt      *idx;
	PetscScalar   *X,*PolyLen,*PolyIdx,*PolyFile;
	PetscInt       imark, imarkx, imarky, imarkz, icellx, icelly, icellz;
	PetscScalar    dx, dy, dz, x, y, z;
	PetscScalar    chLen;
	PetscLogDouble t0, t1;
	char           normalDir[4] = {"xyz"};
	PetscRandom    rctx;
	PetscScalar    cf_rand;
	PetscInt       nPoly;
	PetscScalar    atol;
	PetscScalar    box[4];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "poly_file", filename, "./input/poly.dat"); CHKERRQ(ierr);

	// initialize
	fs = actx->fs;
	x = y = z = dx = dy = dz = 0.0;
	chLen = actx->jr->scal->length;

	// initialize the random number generator
	if(actx->randNoise)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Adding random noise to marker distribution \n");

		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);
	}

	// --- initialize markers ---
	
	// marker counter
	imark  = 0;
	icellx = 0;
	icelly = 0;
	icellz = 0;

	// initialize makers in a processor wise manner
	for(imarkz = 0; imarkz < fs->dsz.ncels*actx->NumPartZ; imarkz++)
	{
		if(!(imarkz%actx->NumPartZ))
		{
			dz = (fs->dsz.ncoor[icellz+1] - fs->dsz.ncoor[icellz]) / (PetscScalar) (actx->NumPartZ);
			z  = fs->dsz.ncoor[icellz] + 0.5*dz;
			icellz++;
		}
		else
		{
			z += dz;
		}
		icelly = 0;

		for(imarky = 0; imarky < fs->dsy.ncels*actx->NumPartY; imarky++)
		{
			if(!(imarky%actx->NumPartY))
			{
				dy = (fs->dsy.ncoor[icelly+1] - fs->dsy.ncoor[icelly]) / (PetscScalar) (actx->NumPartY);
				y  = fs->dsy.ncoor[icelly] + 0.5*dy;
				icelly++;
			}
			else
			{
				y += dy;
			}
			icellx = 0;

			for(imarkx = 0; imarkx < fs->dsx.ncels*actx->NumPartX; imarkx++)
			{
				if(!(imarkx%actx->NumPartX))
				{
					dx = (fs->dsx.ncoor[icellx+1] - fs->dsx.ncoor[icellx]) / (PetscScalar) (actx->NumPartX);
					x  = fs->dsx.ncoor[icellx] + 0.5*dx;
					icellx++;
				}
				else
				{
					x += dx;
				}

				// set marker coordinates
				actx->markers[imark].X[0] = x;
				actx->markers[imark].X[1] = y;
				actx->markers[imark].X[2] = z;
				
				
				if(actx->randNoise)
				{
					// add random noise
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[0] += (cf_rand-0.5)*dx/( (PetscScalar) actx->NumPartX);
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[1] += (cf_rand-0.5)*dy/( (PetscScalar) actx->NumPartY);
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[2] += (cf_rand-0.5)*dz/( (PetscScalar) actx->NumPartZ);
				}

				// increment local counter
				imark++;
			}
		}
	}

	// --- local grid/marker info ---

	// get first global index of marker plane
	tstart[0] = fs->dsx.pstart * actx->NumPartX;
	tstart[1] = fs->dsy.pstart * actx->NumPartY;
	tstart[2] = fs->dsz.pstart * actx->NumPartZ;

	// get local number of markers per direction
	nmark[0]  = fs->dsx.ncels * actx->NumPartX;
	nmark[1]  = fs->dsy.ncels * actx->NumPartY;
	nmark[2]  = fs->dsz.ncels * actx->NumPartZ;

	// get last global index of marker plane
	for(k = 0; k < 3; k++)
	{
		tend[k] = tstart[k] + nmark[k] - 1;
	}

	// how many markers on the marker plane ?
	nidx[0] = nmark[1] * nmark[2]; nidxmax = nidx[0];
	nidx[1] = nmark[0] * nmark[2]; if (nidx[1] > nidxmax) nidxmax = nidx[1];
	nidx[2] = nmark[0] * nmark[1]; if (nidx[2] > nidxmax) nidxmax = nidx[2];

	// read file
	PetscPrintf(PETSC_COMM_WORLD," Loading polygons redundantly from file: %s \n", filename);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd);                               CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header & size of file
	ierr = PetscBinaryRead(fd, &header, 2, PETSC_SCALAR); CHKERRQ(ierr);
	Fsize = (PetscInt)(header[1]);

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)Fsize  *sizeof(PetscScalar),&PolyFile); CHKERRQ(ierr);
	Fcount = 0;

	// read entire file 
	ierr = PetscBinaryRead(fd, PolyFile, Fsize, PETSC_SCALAR); CHKERRQ(ierr);

	// read number of volumes
	VolN = (PetscInt)(PolyFile[Fcount]); Fcount++;
	Nmax = (PetscInt)(PolyFile[Fcount]); Fcount++;
	Lmax = (PetscInt)(PolyFile[Fcount]); Fcount++;

    // allocate space for index array & the coordinates of the largest polygon
	ierr = PetscMalloc((size_t)Nmax  *sizeof(PetscScalar),&PolyLen); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)Nmax  *sizeof(PetscScalar),&PolyIdx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)Lmax*2*sizeof(PetscScalar),&Poly.X);  CHKERRQ(ierr);

	// allocate temporary arrays
	ierr = PetscMalloc((size_t)nidxmax*sizeof(PetscInt),&idx);     CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nidxmax*sizeof(PetscBool),&polyin); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nidxmax*2*sizeof(PetscScalar),&X);  CHKERRQ(ierr);

	// --- loop over all volumes ---
	for(kvol = 0; kvol < VolN; kvol++)
	{
		PetscTime(&t0);		

		// read volume header
		Poly.dir   = (PetscInt)(PolyFile[Fcount]); Fcount++; // normal vector of polygon plane
		Poly.phase = (PetscInt)(PolyFile[Fcount]); Fcount++; // phase that polygon defines
		Poly.type  = (PetscInt)(PolyFile[Fcount]); Fcount++; // type of assigning the phases
		Poly.num   = (PetscInt)(PolyFile[Fcount]); Fcount++; // number of polygon slices defining the volume
		Poly.nmark = 0;

		// define axes the span the polygon plane
		if (Poly.dir==0)
		{
			Poly.ax[0] = 1; Poly.ax[1] = 2;
		}
		else if (Poly.dir==1)
		{
			Poly.ax[0] = 0; Poly.ax[1] = 2;
		}
		else if (Poly.dir==2)
		{
			Poly.ax[0] = 0; Poly.ax[1] = 1;
		}
		else
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "The 'Dir' argument is wrong; should be 0, 1 or 2.");
		}

		// get lengths of polygons (PetscScalar !)
		for(kpoly = 0; kpoly < Poly.num; kpoly++)
		{
			PolyIdx[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// get lengths of polygons (PetscScalar !)
		for (kpoly=0; kpoly<Poly.num;kpoly++)
		{
			PolyLen[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// --- loop through all slices ---
		for(kpoly = 0; kpoly < Poly.num; kpoly++)
		{
			// read polygon
			Poly.len  = (PetscInt)(PolyLen[kpoly]);
			Poly.gidx = (PetscInt)(PolyIdx[kpoly]);
			Poly.lidx = (PetscInt)(PolyIdx[kpoly])-tstart[Poly.dir];

			// check if slice is part of local proc
			if(Poly.gidx >= tstart[Poly.dir] && Poly.gidx <= tend[Poly.dir])
			{
				// read polygon
				for (n=0; n<Poly.len*2;n++)
				{
					Poly.X[n] = PolyFile[Fcount]; Fcount++;
				}

				// get local markers that locate on polygon plane
				ADVMarkSecIdx(actx, Poly.dir, Poly.lidx, idx);

				for(k = 0; k < nidx[Poly.dir]; k++)
				{
					X[k*2]   = actx->markers[idx[k]].X[Poly.ax[0]] * chLen;
					X[k*2+1] = actx->markers[idx[k]].X[Poly.ax[1]] * chLen;
				}

				// get bounding box of a polygon
				nPoly = Poly.len;

				polygon_box(&nPoly, Poly.X, 1e-12, &atol, box);

				in_polygon(nidx[Poly.dir], X, nPoly, Poly.X, box, atol, polyin);

				// set marker phase
				for(k = 0; k < nidx[Poly.dir]; k++)
				{
					if(polyin[k])
					{
						if(Poly.type == 1) // additive
						{
							actx->markers[idx[k]].phase += Poly.phase;
						}
						else if(Poly.type == 2) // grid additive
						{
							if(actx->markers[idx[k]].phase % 2 == 1) // avoid adding twice when contours are over imposed (e.g. at grid intersection)
							{
								actx->markers[idx[k]].phase += Poly.phase;
							}
						}
						else // overwriting
						{
							actx->markers[idx[k]].phase = Poly.phase;
						}
						Poly.nmark++;
					}
				}
			}
			else
			{
				// increase counter of the buffer
				Fcount += Poly.len*2;
			}
		}

		PetscTime(&t1);

		PetscPrintf(PETSC_COMM_WORLD,"[Rank 0] Created vol %lld/%lld [%g sec]: phase %lld, type %lld, %lld slices, %c-normal-dir; found %lld markers \n",(LLD)kvol+1,(LLD)VolN, t1-t0, (LLD)Poly.phase, (LLD)Poly.type,(LLD)Poly.num, normalDir[Poly.dir], (LLD)Poly.nmark);
	}

	// free
	PetscFree(idx);
	PetscFree(polyin);
	PetscFree(X);
	PetscFree(PolyIdx);
	PetscFree(PolyLen);
	PetscFree(Poly.X);
	PetscFree(PolyFile);
	
	if(actx->randNoise)
	{
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);

	// wait until all processors finished reading markers
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished setting markers with polygons\n");

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSecIdx"
void ADVMarkSecIdx(AdvCtx *actx, PetscInt dir, PetscInt Islice, PetscInt *idx)
{
	FDSTAG   *fs;
	PetscInt i,ix,iy,iz,nmarkx,nmarky,nmarkz;
	PetscInt d,c;
	
	// get fdstag info
	fs = actx->fs;

	// get local number of markers
	nmarkx  = fs->dsx.ncels * actx->NumPartX;
	nmarky  = fs->dsy.ncels * actx->NumPartY;
	nmarkz  = fs->dsz.ncels * actx->NumPartZ;

	if(dir == 0) // yz plane
	{
		d = 0;
		c = Islice;
		for(iz=0; iz<nmarkz; iz++)
		{
			for (iy=0; iy<nmarky; iy++)
			{
				idx[d] =c;
				c += nmarkx;
				d++;
			}
		}
	}
	else if(dir == 1) // xz plane
	{
		d = 0;
		c = Islice *nmarkx;
		for(iz=0; iz<nmarkz;iz++ )
		{
			for(ix=0; ix<nmarkx;ix++)
			{
				idx[d] = c;
				c++;
				d++;
			}
			c += nmarkx*nmarky-nmarkx;
		}
	}
	else if(dir == 2) // xy plane
	{
		d = 0;
		for(i=0; i<(nmarkx*nmarky);i++)
		{
			idx[d] = i + (Islice*nmarkx*nmarky);
			d++;
		}
	}

	return;
}
//---------------------------------------------------------------------------
// geometric primitives functions
//---------------------------------------------------------------------------
PetscInt setPhaseSphere(GeomPrim *sphere, Marker *P)
{
	PetscScalar dx, dy, dz;

	dx = P->X[0] - sphere->center[0];
	dy = P->X[1] - sphere->center[1];
	dz = P->X[2] - sphere->center[2];

	if(sqrt(dx*dx + dy*dy + dz*dz) <= sphere->radius)
	{
		P->phase = sphere->phase;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
PetscInt setPhaseBox(GeomPrim *box, Marker *P)
{
	if(P->X[0] >= box->bounds[0] && P->X[0] <= box->bounds[1]
	&& P->X[1] >= box->bounds[2] && P->X[1] <= box->bounds[3]
	&& P->X[2] >= box->bounds[4] && P->X[2] <= box->bounds[5])
	{
		P->phase = box->phase;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
PetscInt setPhaseLayer(GeomPrim *layer, Marker *P)
{
	if(P->X[2] >= layer->bot && P->X[2] <= layer->top)
	{
		P->phase = layer->phase;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
PetscInt setPhaseHex(GeomPrim *hex, Marker *P)
{
	PetscInt    i;
	PetscScalar tol = 1e-6;

	// cell tetrahedrization
	PetscInt tet [] =
	{
		0, 1, 2, 5, // 0
		2, 6, 5, 7, // 1
		0, 3, 2, 7, // 2
		0, 5, 4, 7, // 3
		0, 2, 5, 7, // 4
	};

	// check bounding box
	if(P->X[0] >= hex->bounds[0] && P->X[0] <= hex->bounds[1]
	&& P->X[1] >= hex->bounds[2] && P->X[1] <= hex->bounds[3]
	&& P->X[2] >= hex->bounds[4] && P->X[2] <= hex->bounds[5])
	{
		// check tetrahedrons
		for(i = 0; i < 5; i++)
		{
			if(TetPointTest(hex->coord, tet + 4*i, P->X, tol))
			{
				P->phase = hex->phase;
				return 1;
			}
		}
	}

	return 0;
}
//---------------------------------------------------------------------------
void HexGetBoundingBox(
		PetscScalar *coord,  // hex coordinates
		PetscScalar *bounds) // bounding box
{
	PetscInt     i;
	PetscScalar *x;

	bounds[0] = bounds[1] = coord[0];
	bounds[2] = bounds[3] = coord[1];
	bounds[4] = bounds[5] = coord[2];

	// compute bounding box
	for(i = 1; i < 8; i++)
	{
		x = coord + 3*i;

		if(bounds[0] > x[0]) bounds[0] = x[0];
		if(bounds[1] < x[0]) bounds[1] = x[0];
		if(bounds[2] > x[1]) bounds[2] = x[1];
		if(bounds[3] < x[1]) bounds[3] = x[1];
		if(bounds[4] > x[2]) bounds[4] = x[2];
		if(bounds[5] < x[2]) bounds[5] = x[2];
	}
}
//---------------------------------------------------------------------------
PetscInt TetPointTest(
		PetscScalar *coord, // tetrahedron coordinates
		PetscInt    *ii,    // corner indices
		PetscScalar *xp,    // point coordinate
		PetscScalar  tol)   // relative tolerance
{
	// macro for computing 3x3 matrix determinant
	#define DET PetscAbsScalar(a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31))

	PetscInt     j1, j2, j3, j4;
	PetscScalar  x, y, z, r, s, t, q, d;
	PetscScalar  a11, a12, a13, a21, a22, a23, a31, a32, a33;
	PetscScalar  x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

	// assign point coordinates
	x = xp[0];
	y = xp[1];
	z = xp[2];

	// assign nodal coordinates
	j1 = 3*ii[0]; x1 = coord[j1]; y1 = coord[j1+1]; z1 = coord[j1+2]; // node 1
	j2 = 3*ii[1]; x2 = coord[j2]; y2 = coord[j2+1]; z2 = coord[j2+2]; // node 2
	j3 = 3*ii[2]; x3 = coord[j3]; y3 = coord[j3+1]; z3 = coord[j3+2]; // node 3
	j4 = 3*ii[3]; x4 = coord[j4]; y4 = coord[j4+1]; z4 = coord[j4+2]; // node 4

	// compute total volume
	a11 = x2 - x1;   a12 = x3 - x1;   a13 = x4 - x1;
	a21 = y2 - y1;   a22 = y3 - y1;   a23 = y4 - y1;
	a31 = z2 - z1;   a32 = z3 - z1;   a33 = z4 - z1;
	d   = DET;

	// compute sub-volume (p2-p3-p4-p) (N1)
	a11 = x2 - x;    a12 = x3 - x;    a13 = x4 - x;
	a21 = y2 - y;    a22 = y3 - y;    a23 = y4 - y;
	a31 = z2 - z;    a32 = z3 - z;    a33 = z4 - z;
	r   = DET;

	// compute sub-volume (p1-p3-p4-p) (N2)
	a11 = x1 - x;
	a21 = y1 - y;
	a31 = z1 - z;
	s   = DET;

	// compute sub-volume (p1-p2-p4-p) (N3)
	a12 = x2 - x;
	a22 = y2 - y;
	a32 = z2 - z;
	t   = DET;

	// compute sub-volume (p1-p2-p3-p) (N4)
	a13 = x3 - x;
	a23 = y3 - y;
	a33 = z3 - z;
	q   = DET;

	// point test
	if(r + s + t + q > d*(1.0 + tol)) return 0;

	return 1;
}
//---------------------------------------------------------------------------
