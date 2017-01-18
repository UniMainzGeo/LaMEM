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
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "advect.h"
#include "marker.h"
#include "paraViewOutBin.h"
#include "paraViewOutSurf.h"
#include "paraViewOutMark.h"
#include "paraViewOutAVD.h"
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

	// initialize variables for marker control
	actx->nmin = (PetscInt) (actx->NumPartX*actx->NumPartY*actx->NumPartZ/2); // min # of markers/cell -50%
	actx->nmax = (PetscInt) (actx->NumPartX*actx->NumPartY*actx->NumPartZ*3); // max # of markers/cell 300%
	actx->avdx = actx->NumPartX * 3;
	actx->avdy = actx->NumPartY * 3;
	actx->avdz = actx->NumPartZ * 3;

	// initialize coordinates and add random noise if required for hard-coded setups
	if(actx->msetup != _FILES_
	&& actx->msetup != _POLYGONS_)
	{
		ierr = ADVMarkInitCoord(actx); CHKERRQ(ierr);
	}

	// display info
	PetscPrintf(PETSC_COMM_WORLD," Marker setup employed [msetup] : ");

	// initialize marker phase, temperature, etc.
	if     (actx->msetup == _FILES_)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","files");      ierr = ADVMarkInitFiles     (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _POLYGONS_)   { PetscPrintf(PETSC_COMM_WORLD,"%s\n","polygons");   ierr = ADVMarkInitPolygons  (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _DIAPIR_)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","diapir");     ierr = ADVMarkInitDiapir    (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _FBLOCK_)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","block");      ierr = ADVMarkInitBlock     (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _SUBDUCTION_) { PetscPrintf(PETSC_COMM_WORLD,"%s\n","subduction"); ierr = ADVMarkInitSubduction(actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _FOLDING_)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","folding");    ierr = ADVMarkInitFolding   (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _DETACHMENT_) { PetscPrintf(PETSC_COMM_WORLD,"%s\n","detachment"); ierr = ADVMarkInitDetachment(actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _SLAB_)       { PetscPrintf(PETSC_COMM_WORLD,"%s\n","slab");       ierr = ADVMarkInitSlab      (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _SPHERES_)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","spheres");    ierr = ADVMarkInitSpheres   (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _BANDS_)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","bands");      ierr = ADVMarkInitBands     (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _DOMES_)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","domes");      ierr = ADVMarkInitDomes     (actx, fb); CHKERRQ(ierr); }

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
	char        *filename;
	PetscViewer  view_out;
	PetscScalar *markbuf, *markptr, header, chLen, chTemp, Tshift, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!actx->saveMark) PetscFunctionReturn(0);

	// access context
	chLen  = actx->jr->scal->length;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;

	PetscPrintf(PETSC_COMM_WORLD," Saving markers in parallel to files: ./%s/%s.xxx.dat \n", actx->savePath, actx->saveName);

	// create directory
	ierr = LaMEMCreateOutputDirectory(actx->savePath); CHKERRQ(ierr);

	// compile file name
	asprintf(&filename, "./%s/%s.%1.8lld.dat", actx->savePath, actx->saveName, (LLD)actx->iproc);

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
	PetscScalar  xs, ys, zs;
	PetscScalar  xe, ye, ze;
	PetscInt     *numMarkCell, rbuf[4], sbuf[4];
	PetscInt     i, maxid, NumInvalidPhase, numNonLocal, numEmpty, numSparse;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// get maximum Phase
	maxid = actx->jr->numPhases - 1;

	// get local mesh sizes
	GET_DOMAIN_BOUNDS(xs, xe, fs->dsx)
	GET_DOMAIN_BOUNDS(ys, ye, fs->dsy)
	GET_DOMAIN_BOUNDS(zs, ze, fs->dsz)

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
		if(X[0] < xs || X[0] > xe
		|| X[1] < ys || X[1] > ye
		|| X[2] < zs || X[2] > ze) numNonLocal++;

		// count number of markers in the cells
		numMarkCell[actx->cellnum[i]]++;
	}

	// count empty & sparse cells
	numEmpty  = 0;
	numSparse = 0;

	for(i = 0; i < fs->nCells; i++)
	{
		if(numMarkCell[i] == 0) numEmpty++;
		if(numMarkCell[i] <  8) numSparse++;
	}

	// get global figures
	if(actx->nproc != 1)
	{
		sbuf[0] = NumInvalidPhase;
		sbuf[1] = numNonLocal;
		sbuf[2] = numEmpty;
		sbuf[3] = numSparse;

		ierr = MPI_Allreduce(sbuf, rbuf, 4, MPIU_INT, MPI_SUM, actx->icomm); CHKERRQ(ierr);

		NumInvalidPhase = rbuf[0];
		numNonLocal     = rbuf[1];
		numEmpty        = rbuf[2];
		numSparse       = rbuf[3];
	}

	// print diagnostics
	if(NumInvalidPhase)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of markers that have invalid phase ID: %lld\n", (LLD)NumInvalidPhase); CHKERRQ(ierr);
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

	if(numSparse)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "WARNING! Number of cells with less than 8 markers: %lld.\n", (LLD)numSparse); CHKERRQ(ierr);
	}

	if(error)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Problems with initial marker distribution (see the above message)");
	}

	// clear
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

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
	char          filename[MAX_PATH_LEN];
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
	ierr = PetscMemzero(filename, sizeof(char)*MAX_PATH_LEN); CHKERRQ(ierr);

	ierr = getStringParam(fb, _OPTIONAL_, "temp_file", filename, MAX_PATH_LEN); CHKERRQ(ierr);

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
	char        *filename, name[MAX_NAME_LEN], path[MAX_PATH_LEN];
	PetscScalar *markbuf, *markptr, header, chTemp, chLen, Tshift, s_nummark;
	PetscInt     imark, nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name & path
	ierr = PetscMemzero(name, sizeof(char)*MAX_NAME_LEN); CHKERRQ(ierr);
	ierr = PetscMemzero(path, sizeof(char)*MAX_PATH_LEN); CHKERRQ(ierr);

	ierr = getStringParam(fb, _REQUIRED_, "mark_load_name", name, MAX_NAME_LEN); CHKERRQ(ierr);
	ierr = getStringParam(fb, _REQUIRED_, "mark_load_path", path, MAX_PATH_LEN); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Loading markers in parallel from files: ./%s/%s.xxx.dat \n", path, name);

	// compile input file name
	asprintf(&filename, "./%s/%s.%1.8lld.dat", path, name, (LLD)actx->iproc);

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
#define __FUNCT__ "ADVMarkInitDiapir"
PetscErrorCode ADVMarkInitDiapir(AdvCtx *actx, FB *fb)
{
	PetscInt imark;

	PetscFunctionBegin;

	// Diapir setup - not done - should have some perturbation
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD," DIAPIR SETUP\n");
	PetscPrintf(PETSC_COMM_WORLD," Setup Parameters: Setup.Diapir_Hi = [%g]\n",user->Setup_Diapir_Hi);

	for(imark = 0; imark < actx->nummark; imark++)
	{
		// phase
		if(actx->markers[imark].X[2] >= user->Setup_Diapir_Hi) actx->markers[imark].phase = 1;
		else                                                   actx->markers[imark].phase = 0;

		// temperature
		actx->markers[imark].T = 0.0;
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitBlock"
PetscErrorCode ADVMarkInitBlock(AdvCtx *actx, FB *fb)
{
	// falling block setup

	PetscInt    imark, nel_x, nel_y, nel_z;
	PetscScalar dx,dy,dz;
	PetscScalar bleft, bright, bfront, bback, bbottom, btop;
	PetscScalar blx, bly, blz;
	PetscBool   b2D, b2Dy;

	PetscErrorCode ierr;
	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  FALLING BLOCK SETUP \n");

	// number of elements on finest resolution
	nel_x = user->nel_x;
	nel_y = user->nel_y;
	nel_z = user->nel_z;

	// spacing
	dx = user->W/((PetscScalar)nel_x);
	dy = user->L/((PetscScalar)nel_y);
	dz = user->H/((PetscScalar)nel_z);

	// block dimensions
	blx = 0.5*(PetscScalar)nel_x*dx;
	bly = 0.5*(PetscScalar)nel_y*dy;
	blz = 0.5*(PetscScalar)nel_z*dz;

	bleft   = 0.25*(PetscScalar)nel_x*dx + user->x_left ; bright = bleft   + blx; // left and right side of block
	bfront  = 0.25*(PetscScalar)nel_y*dy + user->y_front; bback  = bfront  + bly; // front and back side of block
	bbottom = 0.25*(PetscScalar)nel_z*dz + user->z_bot  ; btop   = bbottom + blz; // bottom and top side of block

	PetscPrintf(PETSC_COMM_WORLD,"      Block coordinates: [left,right]=[%g,%g]; \n",bleft,bright);
    PetscPrintf(PETSC_COMM_WORLD,"                         [front,back]=[%g,%g]; \n",bfront,bback);
    PetscPrintf(PETSC_COMM_WORLD,"                         [bottom,top]=[%g,%g]; \n",bbottom,btop);
    

	b2D  = PETSC_FALSE;
	b2Dy = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL, NULL, "-Block_2D" , &b2D , NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL, NULL, "-Block_2Dy", &b2Dy, NULL); CHKERRQ(ierr);

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 1.0 - actx->markers[imark].X[2];

		if((actx->markers[imark].X[2] > bbottom)
		&& (actx->markers[imark].X[2] < btop)
		&& (actx->markers[imark].X[0] > bleft)
		&& (actx->markers[imark].X[0] < bright)
		&& (b2Dy == PETSC_FALSE))
		{
			if(b2D == PETSC_TRUE)
			{
				// 2D block
				actx->markers[imark].phase = 1;
				actx->markers[imark].T     = 1.0;
			}
			else
			{
				if((actx->markers[imark].X[1] > bfront)
				&& (actx->markers[imark].X[1] < bback))
				{
					// 3D block
					actx->markers[imark].phase = 1;
					actx->markers[imark].T     = 1.0;
				}
			}
		}
		if(b2Dy == PETSC_TRUE)
		{
			// create 2D block in y-direction
			if((actx->markers[imark].X[2] > bbottom)
			&& (actx->markers[imark].X[2] < btop)
			&& (actx->markers[imark].X[1] > bfront)
			&& (actx->markers[imark].X[1] < bback))
			{
				// 2D block
				actx->markers[imark].phase = 1;
				actx->markers[imark].T     = 1.0;
			}
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSubduction"
PetscErrorCode ADVMarkInitSubduction(AdvCtx *actx, FB *fb)
{
	// subduction setup with air

	PetscScalar H_air, SlabThickness, SlabWidth, SlabLength;
	PetscScalar DistanceFromLeft, SlabMaxSubdDepth,Slab_Xright, dz;
	PetscScalar Air_thickness, Slab_ThicknessFactor, Air_ThicknessFactor, Slab_WidthFactor;
	PetscScalar chLen;
	PetscInt    imark, nz;

	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  SUBDUCTION WITH STICKY AIR SETUP \n");

	// TOTAL number of CELLS and spacing in z-direction
	nz = user->nel_z;
	dz = user->H/((PetscScalar)nz);

	// read some input parameters
	Slab_ThicknessFactor = 0.1; // in % of the total H
	PetscOptionsGetReal(NULL, NULL ,"-Slab_ThicknessFactor", &Slab_ThicknessFactor, NULL);

	Air_ThicknessFactor  = 0.1; // in % of the total H
	PetscOptionsGetReal(NULL, NULL ,"-Air_ThicknessFactor" , &Air_ThicknessFactor , NULL);

	Slab_WidthFactor     = 1.0;
	PetscOptionsGetReal(NULL, NULL ,"-Slab_WidthFactor"    , &Slab_WidthFactor    , NULL);

	// setup parameters
	SlabThickness = (PetscScalar) (round(user->H/dz)*dz*Slab_ThicknessFactor);
	Air_thickness = (PetscScalar) (round(user->H/dz)*dz*Air_ThicknessFactor );
	H_air         = (PetscScalar) (user->z_bot + user->H - Air_thickness    );

	SlabWidth        = Slab_WidthFactor*user->L;
	SlabLength       = 0.4 * (user->W);
	DistanceFromLeft = (user->W) * 0.1;
	SlabMaxSubdDepth = 2 * SlabThickness;
	Slab_Xright      = DistanceFromLeft + SlabLength-(user->H-H_air);

	// get characteristic length
	chLen  = actx->jr->scal.length;

	PetscPrintf(PETSC_COMM_WORLD," Setup Parameters: SlabThickness = [%g km] AirThickness = [%g km] H = [%g km] SlabWidthFactor = [%g] \n",SlabThickness*chLen*0.001, Air_thickness*chLen*0.001, user->H*chLen*0.001, Slab_WidthFactor);

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		if((actx->markers[imark].X[2] > (H_air - SlabThickness))
		&& (actx->markers[imark].X[2] <  H_air)
		&& (actx->markers[imark].X[0] > (user->x_left  + DistanceFromLeft))
		&& (actx->markers[imark].X[0] < (user->x_left  + DistanceFromLeft + SlabLength))
		&& (actx->markers[imark].X[1] < (user->y_front + SlabWidth)))
		{
			actx->markers[imark].phase = 1; // slab flat part
		}

		if((actx->markers[imark].X[2] < (user->H - (actx->markers[imark].X[0]-Slab_Xright)))
		&& (actx->markers[imark].X[2] > (user->H-2*SlabThickness - (actx->markers[imark].X[0]-Slab_Xright)))
		&& (actx->markers[imark].X[1] < (user->y_front+SlabWidth))
		&& (actx->markers[imark].X[2] > (H_air-SlabMaxSubdDepth)))
		{
			actx->markers[imark].phase = 1; // slab inclined part
		}

		if (actx->markers[imark].X[2] > H_air)
		{
			actx->markers[imark].phase = 2; // air
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFolding"
PetscErrorCode ADVMarkInitFolding(AdvCtx *actx, FB *fb)
{
	// multilayer folding setup (Zagros)

	PetscInt    i, imark;
	PetscBool   flg, DisplayInfo, FoldingSetupDimensionalUnits;
	PetscScalar zbot[10], ztop[10], BottomStepPerturbation_Amplitude;
	PetscScalar chLen;

	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  MULTILAYER FOLDING SETUP \n");

	// initialize arrays for fractions
	for(i = 0;i < 10; i++) { zbot[i] = 0.0; ztop[i] = 0.0; }

	// get input data
	DisplayInfo = PETSC_TRUE;

    PetscOptionsHasName(NULL, NULL,"-FoldingSetupDimensionalUnits"      ,&FoldingSetupDimensionalUnits);
    if (DisplayInfo && FoldingSetupDimensionalUnits) PetscPrintf(PETSC_COMM_WORLD,"   Assuming that the folding setup is given in dimensional units \n");
    
	PetscOptionsGetReal(NULL, NULL,"-Layer1_bottom"      ,&zbot[0], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer1_top"         ,&ztop[0], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=1 from z=[%g,%g]\n",zbot[0],ztop[0]);

	PetscOptionsGetReal(NULL, NULL,"-Layer2_bottom"      ,&zbot[1], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer2_top"         ,&ztop[1], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=2 from z=[%g,%g]\n",zbot[1],ztop[1]);

	PetscOptionsGetReal(NULL, NULL,"-Layer3_bottom"      ,&zbot[2], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer3_top"         ,&ztop[2], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=3 from z=[%g,%g]\n",zbot[2],ztop[2]);

	PetscOptionsGetReal(NULL, NULL,"-Layer4_bottom"      ,&zbot[3], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer4_top"         ,&ztop[3], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=4 from z=[%g,%g]\n",zbot[3],ztop[3]);

	PetscOptionsGetReal(NULL, NULL,"-Layer5_bottom"      ,&zbot[4], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer5_top"         ,&ztop[4], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=5 from z=[%g,%g]\n",zbot[4],ztop[4]);

	PetscOptionsGetReal(NULL, NULL,"-Layer6_bottom"      ,&zbot[5], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer6_top"         ,&ztop[5], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=6 from z=[%g,%g]\n",zbot[5],ztop[5]);

	PetscOptionsGetReal(NULL, NULL,"-Layer7_bottom"      ,&zbot[6], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer7_top"         ,&ztop[6], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=7 from z=[%g,%g]\n",zbot[6],ztop[6]);

	PetscOptionsGetReal(NULL, NULL,"-Layer8_bottom"      ,&zbot[7], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer8_top"         ,&ztop[7], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=8 from z=[%g,%g]\n",zbot[7],ztop[7]);

	PetscOptionsGetReal(NULL, NULL,"-Layer9_bottom"      ,&zbot[8], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer9_top"         ,&ztop[8], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=9 from z=[%g,%g]\n",zbot[8],ztop[8]);

	PetscOptionsGetReal(NULL, NULL,"-Layer10_bottom"      ,&zbot[9], NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer10_top"         ,&ztop[9], &flg);
	if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase=10 from z=[%g,%g]\n",zbot[9],ztop[9]);

	// get characteristic length
	chLen  = actx->jr->scal.length;

	if (FoldingSetupDimensionalUnits)
	{
		for(i = 0; i < 10; i++) { zbot[i] = zbot[i]/chLen; ztop[i] = ztop[i]/chLen;}
	}
	else
	{
		// transform fractions into domain coordinates
		for(i = 0; i < 10; i++) { zbot[i] = zbot[i]*user->H; ztop[i] = ztop[i]*user->H;}
	}

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		for(i = 0; i < 10; i++)
		{
			if((actx->markers[imark].X[2] >= zbot[i])
			&& (actx->markers[imark].X[2] < ztop[i]))
			{
				actx->markers[imark].phase = i+1;
			}
		}
	}

	BottomStepPerturbation_Amplitude = 0.1;
	PetscOptionsGetReal(NULL, NULL,"-BottomStepPerturbation_Amplitude"         ,&BottomStepPerturbation_Amplitude, &flg);
	if (flg) PetscPrintf(PETSC_COMM_WORLD,"   Using bottom step like perturbation with amplitude %f  \n", BottomStepPerturbation_Amplitude);

	if (FoldingSetupDimensionalUnits) BottomStepPerturbation_Amplitude = BottomStepPerturbation_Amplitude/chLen;

	if (DisplayInfo && flg){
		PetscScalar BottomStepPerturbation_X, BottomStepPerturbation_Y;

		// add a step-like pertubation at the bottom layer
		BottomStepPerturbation_X = 0;
		PetscOptionsGetReal(NULL, NULL,"-BottomStepPerturbation_X"         ,&BottomStepPerturbation_X, &flg);
		if (FoldingSetupDimensionalUnits){
			BottomStepPerturbation_X = BottomStepPerturbation_X/chLen;
		}

		BottomStepPerturbation_Y = user->y_front + user->L;
		PetscOptionsGetReal(NULL, NULL,"-BottomStepPerturbation_Y"         ,&BottomStepPerturbation_Y, &flg);
		if (FoldingSetupDimensionalUnits & flg){
			BottomStepPerturbation_Y = BottomStepPerturbation_Y/chLen;
		}


		// loop over local markers
		for(imark = 0; imark < actx->nummark; imark++)
		{
			if (actx->markers[imark].phase > 0){
				if ((actx->markers[imark].X[0]<BottomStepPerturbation_X) & (actx->markers[imark].X[1]<BottomStepPerturbation_Y)){
					if ( actx->markers[imark].X[2]< (BottomStepPerturbation_Amplitude+zbot[0]) ){
						actx->markers[imark].phase = 0;
					}
				}
			}
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitDetachment"
PetscErrorCode ADVMarkInitDetachment(AdvCtx *actx, FB *fb)
{
	// 1-layer over detachment with two linear shaped perturbation (Grasemann & Schmalholz 2012)
	// one perturbation is fixed in position, the other one's position is varied with -Heterogeneity_Offset

	PetscInt    imark;
	PetscScalar zbot = 0.0, ztop = 0.0;
	PetscScalar Het_L, Het_W, Het_A, Het_Off;
	PetscScalar chLen;
	PetscBool   flg, DisplayInfo;

	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  ONE-LAYER OVER DETACHMENT WITH 2 LINEAR PERTURBATIONS SETUP \n");

	DisplayInfo = PETSC_TRUE;
	PetscOptionsGetReal(NULL, NULL,"-Heterogeneity_L"      , &Het_L  , NULL);
	PetscOptionsGetReal(NULL, NULL,"-Heterogeneity_W"      , &Het_W  , NULL);
	PetscOptionsGetReal(NULL, NULL,"-Heterogeneity_A"      , &Het_A  , NULL);
	PetscOptionsGetReal(NULL, NULL,"-Heterogeneity_Offset" , &Het_Off, NULL);

	PetscOptionsGetReal(NULL, NULL,"-Layer1_bottom", &zbot, NULL);
	PetscOptionsGetReal(NULL, NULL,"-Layer1_top"   , &ztop, &flg);
	if(DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD," - Adding Layer with Phase = 1 from z=[%g,%g]\n", zbot, ztop);

	// transform fractions into domain coordinates
	zbot = zbot*user->H; ztop = ztop*user->H;

	// get characteristic length
	chLen  = actx->jr->scal.length;

	// non-dimensionalization
	Het_L   = Het_L  /chLen;
	Het_W   = Het_W  /chLen;
	Het_A   = Het_A  /chLen;
	Het_Off = Het_Off/chLen;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		// one layer properties
		if((actx->markers[imark].X[2] >= zbot)
		&& (actx->markers[imark].X[2] <  ztop))
		{
			actx->markers[imark].phase = 1;
		}

		// seed fixed
		if((actx->markers[imark].X[1] >= (user->L*0.5 - Het_L*0.5))
		&& (actx->markers[imark].X[1] <= (user->L*0.5 + Het_L*0.5))
		&& (actx->markers[imark].X[0] <= (Het_W))
		&& (actx->markers[imark].X[2] <= (zbot + Het_A)))
		{
			actx->markers[imark].phase = 0;
		}

		// seed with offset
		if((actx->markers[imark].X[1] >= (Het_Off + user->L*0.5 - Het_L*0.5))
		&& (actx->markers[imark].X[1] <= (Het_Off + user->L*0.5 + Het_L*0.5))
		&& (actx->markers[imark].X[0] >= (user->W - Het_W))
		&& (actx->markers[imark].X[2] <= (zbot    + Het_A)))
		{
			actx->markers[imark].phase = 0;
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSlab"
PetscErrorCode ADVMarkInitSlab(AdvCtx *actx, FB *fb)
{
	// slab detachment (Thieulot et al. 2014)
	// domain parameters should be given in params file because they are needed to create the grid
	// W = 1000.0e3; L = 500.0e3; H = 700.0e3; ThicknessAir = 40.0e3;

	PetscInt imark;

	// input parameters - chosen as function of L (y-dir)
	PetscScalar hslab = 0.16; // thickness of slab L = 80 km
	PetscScalar lslab = 0.5;  // length of slab    L = 250 km
	PetscScalar dslab = 0.5;  // depth of slab     L = 250 km
	PetscScalar dair  = 0.08; // thickness air     L = 40 km
	PetscScalar chLen_km;

	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  SLAB DETACHMENT SETUP \n");

	// ACHTUNG!
	// get characteristic length in km
	if (actx->jr->scal.utype == _SI_) chLen_km  = actx->jr->scal.length/1000.0;
	else                              chLen_km  = actx->jr->scal.length;

	// non-dimensionalization
	hslab = hslab*user->L;
	lslab = lslab*user->L;
	dslab = dslab*user->L;
	dair  = dair *user->L;

	PetscPrintf(PETSC_COMM_WORLD," Setup Parameters: W = [%g km] L = [%g km] H = [%g km] SlabDimensions = [%g km,%g km,%g km] DepthSlab = [%g km] ThicknessAir = [%g km] \n",user->W*chLen_km,user->L*chLen_km,user->H*chLen_km, hslab*chLen_km,lslab*chLen_km,dslab*chLen_km,dslab*chLen_km,dair*chLen_km);

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		// initialize
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		// flat slab
		if(actx->markers[imark].X[2] >= (user->H-hslab - dair)
		&& actx->markers[imark].X[2] <= (user->H-dair))
		{
			actx->markers[imark].phase = 1;
		}
		else
		{	// vertical slab
			if(actx->markers[imark].X[0] >= (user->W/2.0 - hslab/2.0)
			&& actx->markers[imark].X[0] <= (user->W/2.0 + hslab/2.0)
			&& actx->markers[imark].X[1] <= lslab
			&& actx->markers[imark].X[2] >= (user->H - hslab - dslab - dair)
			&& actx->markers[imark].X[2] <= (user->H - hslab - dair))
			{
				actx->markers[imark].phase = 1;
			}
		}

		// sticky air layer
		if(actx->markers[imark].X[2] >= (user->H - dair))
		{
			actx->markers[imark].phase = 2;
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSpheres"
PetscErrorCode ADVMarkInitSpheres(AdvCtx *actx, FB *fb)
{
	// multiple falling spheres

	PetscInt    nsphere = 10, ind, imark;
	PetscScalar rsphere = 0.1;
	PetscScalar x = 0.0, y = 0.0, z = 0.0;
	PetscScalar xc[20], yc[20], zc[20];

	PetscFunctionBegin;
/*
	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  MULTIPLE FALLING SPHERES SETUP \n");

	// get data
	PetscOptionsGetInt (NULL, NULL, "-NumSpheres",   &nsphere, NULL);
	PetscOptionsGetReal(NULL, NULL, "-SphereRadius", &rsphere, NULL);

	// distribution - set for now such that the spheres do not touch
	xc[0] = 0.25; yc[0] = 0.15; zc[0] = 0.30;
	xc[1] = 0.15; yc[1] = 0.65; zc[1] = 0.30;
	xc[2] = 0.25; yc[2] = 0.55; zc[2] = 0.70;
	xc[3] = 0.35; yc[3] = 0.25; zc[3] = 0.70;
	xc[4] = 0.45; yc[4] = 0.35; zc[4] = 0.30;
	xc[5] = 0.35; yc[5] = 0.65; zc[5] = 0.50;
	xc[6] = 0.65; yc[6] = 0.15; zc[6] = 0.50;
	xc[7] = 0.55; yc[7] = 0.45; zc[7] = 0.50;
	xc[8] = 0.65; yc[8] = 0.35; zc[8] = 0.70;
	xc[9] = 0.55; yc[9] = 0.65; zc[9] = 0.30;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		// initialize
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		// spheres
		for(ind = 0; ind < nsphere; ind++)
		{
			x = actx->markers[imark].X[0] - xc[ind];
			y = actx->markers[imark].X[1] - yc[ind];
			z = actx->markers[imark].X[2] - zc[ind];

			if(x*x + y*y + z*z <= rsphere*rsphere)
			{
				actx->markers[imark].phase = 1;
			}
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitBands"
PetscErrorCode ADVMarkInitBands(AdvCtx *actx, FB *fb)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;
/*
	scal   = actx->jr->scal.length;
	Tshift = actx->jr->scal.Tshift;

	// set default values
	H      = 10.0;               // layer thickness               [km]
	Hb     = 2.0;                // basal layer thickness         [km]
	size   = 1.0;                // inclusion size in XZ-plane    [km]
	offset = 0.0;                // inclusion offset along X axis [km]
	length = (user->L*scal)/2.0; // inclusion length along Y axis [km]

	// get layer thickness
	ierr = PetscOptionsGetReal(NULL, NULL, "-H_layer",  &H,  NULL);  CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-H_bottom", &Hb, NULL);  CHKERRQ(ierr);

	// get inclusion parameters
	ierr = PetscOptionsHasName(NULL, NULL, "-Inclusion_use", &use_inc); CHKERRQ(ierr);
	if(use_inc == PETSC_TRUE)
	{
		ierr = PetscOptionsGetReal(NULL, NULL, "-Inclusion_size",   &size,   NULL); CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL, NULL, "-Inclusion_offset", &offset, NULL); CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL, NULL, "-Inclusion_length", &length, NULL); CHKERRQ(ierr);
	}

	// scale
	H      /= scal;
	Hb     /= scal;
	size   /= scal;
	offset /= scal;
	length /= scal;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		// assign phase in layers
		if     (z >  Hb + H)            P->phase = 2; // air
		else if(z >= Hb && z <= Hb + H) P->phase = 1; // matrix
		else                            P->phase = 0; // basal layer

		// check inclusion
		if(use_inc == PETSC_TRUE && z >= Hb && z <= Hb + size
		&& ((y <= length         && x >= -(offset + size)/2.0 && x <= -(offset - size)/2.0)
		||  (y >= user->L-length && x >=  (offset - size)/2.0 && x <=  (offset + size)/2.0))) P->phase = 0;

		// assign temperature
		P->T = Tshift;
	}
*/
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitDomes"
PetscErrorCode ADVMarkInitDomes(AdvCtx *actx, FB *fb)
{
	// water phase    -> 0
	// basement phase -> 1
	// salt phase     -> 2

	Marker     *P;
	PetscInt    imark, anhydrite_phase;
	PetscScalar surf_level, anhydrite_bot, anhydrite_top, x, z, scal, base_top;

	PetscErrorCode ierr;
	PetscFunctionBegin;
/*
	scal = actx->jr->scal.length; OutBufPut3DVecComp

	// get anhydrite layer parameters (by default salt)
	anhydrite_bot   = user->z_bot;
	anhydrite_top   = user->z_bot;
	anhydrite_phase = 2;

	ierr = PetscOptionsGetReal(NULL, NULL, "-domes_anhydrite_bot",   &anhydrite_bot,   NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-domes_anhydrite_top",   &anhydrite_top,   NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt (NULL, NULL, "-domes_anhydrite_phase", &anhydrite_phase, NULL); CHKERRQ(ierr);

	// scale
	anhydrite_bot /= scal;
	anhydrite_top /= scal;

	surf_level = actx->jr->avg_topo;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		z = P->X[2];

		// water
		P->phase = 0;

		// salt
		if(z <= surf_level) P->phase = 2;

		// anhydrite
		if(z < anhydrite_top && z > anhydrite_bot) P->phase = anhydrite_phase;

		// basement
		if(x < 0.0) base_top = user->z_bot - (surf_level-user->z_bot)*(2.0*x/user->W);
		else        base_top = user->z_bot + (surf_level-user->z_bot)*(2.0*x/user->W);

		if(z <= base_top) P->phase = 1;

		// assign temperature
		P->T = 0.0;
	}
*/
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
	char           filename[MAX_PATH_LEN];
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
	ierr = PetscMemzero(filename, sizeof(char)*MAX_NAME_LEN); CHKERRQ(ierr);

	ierr = getStringParam(fb, _REQUIRED_, "poly_file", filename, MAX_PATH_LEN); CHKERRQ(ierr);

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

	// optionally read temperature from file
	ierr = ADVMarkSetTempFromFile(actx, fb); CHKERRQ(ierr);

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
