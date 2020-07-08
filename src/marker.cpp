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
#include "bc.h"
#include "surf.h"

/*
#START_DOC#
#END_DOC#
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInit"
PetscErrorCode ADVMarkInit(AdvCtx *actx, FB *fb)
{
	FDSTAG    *fs;
	PetscInt  nmarkx, nmarky, nmarkz, nummark;
	PetscBool LoadPhaseDiagrams;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

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

	// initialize markers
	if     (actx->msetup == _GEOM_)       { ierr = ADVMarkInitGeom    (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _FILES_)      { ierr = ADVMarkInitFiles   (actx, fb); CHKERRQ(ierr); }
	else if(actx->msetup == _POLYGONS_)   { ierr = ADVMarkInitPolygons(actx, fb); CHKERRQ(ierr); }

	// set temperature (optional methods)

	// linear gradient
	ierr = ADVMarkSetTempGrad(actx); CHKERRQ(ierr);

	// from file
	ierr = ADVMarkSetTempFile(actx, fb); CHKERRQ(ierr);

	// phase-based
	ierr = ADVMarkSetTempPhase(actx); CHKERRQ(ierr);

	// Load phase diagrams for the phases where it is required + interpolate the reference density for the first timestep
	LoadPhaseDiagrams = PETSC_FALSE;
	
	for(PetscInt i = 0; i < actx->jr->dbm->numPhases; i++)
	{
		if(actx->jr->dbm->phases[i].pdAct)
		{
			LoadPhaseDiagrams = PETSC_TRUE;
		}
	}

	if(LoadPhaseDiagrams)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Phase Diagrams: \n");
	}	

	for(PetscInt i=0; i<actx->jr->dbm->numPhases; i++)
	{
		if(actx->jr->dbm->phases[i].pdAct)
		{
			PetscPrintf(PETSC_COMM_WORLD,"   Phase %i,  ", i);

			ierr = LoadPhaseDiagram(actx, actx->jr->dbm->phases, i); CHKERRQ(ierr);
		}
	}

	if(LoadPhaseDiagrams)
	{
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}


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
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkPerturb"
PetscErrorCode ADVMarkPerturb(AdvCtx *actx)
{
	FDSTAG      *fs;
	PetscScalar *X;
	PetscInt     i, ID, I, J, K, nx,ny;
	PetscScalar  dx,dy,dz;
	PetscRandom  rctx;
	PetscScalar  cf_rand;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// return if not set
	if(!actx->randNoiseGP) PetscFunctionReturn(0);

	PetscPrintf(PETSC_COMM_WORLD,"Apply random noise to markers after initialization\n");

	fs = actx->fs;

	// get random number context
	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);

	// get number of cells
	nx = fs->dsx.ncels;
	ny = fs->dsy.ncels;

	// loop over all local particles
	for(i = 0; i < actx->nummark; i++)
	{
		// get marker coordinates
		X = actx->markers[i].X;

		// get consecutive index of the host cell
		ID = actx->cellnum[i];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get subgrid cell widths
		dx = SIZE_CELL(I, 0, fs->dsx)/(PetscScalar)actx->NumPartX;
		dy = SIZE_CELL(J, 0, fs->dsy)/(PetscScalar)actx->NumPartY;
		dz = SIZE_CELL(K, 0, fs->dsz)/(PetscScalar)actx->NumPartZ;

		// Perturb marker location
		ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
		X[0] += (cf_rand - 0.5)*dx;
		ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
		X[1] += (cf_rand - 0.5)*dy;
		ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
		X[2] += (cf_rand - 0.5)*dz;
	}

	// destroy random context
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSave"
PetscErrorCode ADVMarkSave(AdvCtx *actx)
{
	int            fd;
	PetscInt       imark;
	Marker         *P;
	PetscViewer    view_out;
	PetscLogDouble t;
	char           *filename, path[_str_len_];
	PetscScalar    *markbuf, *markptr, header, chLen, chTemp, Tshift, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(actx->advect == ADV_NONE) PetscFunctionReturn(0);

	if(!actx->saveMark) PetscFunctionReturn(0);

	PrintStart(&t, "Saving markers in parallel to", actx->saveFile);

	// access context
	chLen  = actx->jr->scal->length;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;


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

	PrintDone(t);

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
#define __FUNCT__ "ADVMarkSetTempGrad"
PetscErrorCode ADVMarkSetTempGrad(AdvCtx *actx)
{
	// initialize temperature on markers based on linear gradient
	FDSTAG      *fs;
	BCCtx       *bc;
	Marker      *P;
	PetscInt     imark, nummark;
	PetscScalar  dTdz, zbot, ztop, zp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	bc      = actx->jr->bc;
	fs      = actx->fs;
	nummark = actx->nummark;

	// return if not set
	if(!bc->initTemp) PetscFunctionReturn(0);

	// get grid coordinate bounds in z-direction
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, &ztop); CHKERRQ(ierr);

	// override top boundary with free surface level
	if(actx->surf->UseFreeSurf)
	{
		ztop = actx->surf->InitLevel;
	}

	// get temperature gradient in z-direction
	dTdz = (bc->Ttop - bc->Tbot)/(ztop - zbot);

	// set temperature based on temperature gradient
	for(imark = 0; imark < nummark; imark++)
	{
		// get current marker
		P = &actx->markers[imark];

		// get global marker coordinates
		zp = P->X[2];

		// set temperature
		if(zp > ztop) P->T = bc->Ttop;
		else          P->T = bc->Tbot + dTdz*(zp - zbot);
	}

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSetTempPhase"
PetscErrorCode ADVMarkSetTempPhase(AdvCtx *actx)
{
	// initialize temperature on markers based on phase temperature
	Material_t  *phases;
	Marker      *P;
	PetscInt     i, n, phase_set, imark, nummark;
	PetscScalar  phase_temp[_max_num_phases_];

	PetscFunctionBegin;

	n         = actx->dbm->numPhases;
	phases    = actx->dbm->phases;
	nummark   = actx->nummark;
	phase_set = 0;

	// set temperature based on phase
	for(i = 0; i < n; i++)
	{
		if(phases[i].T) { phase_temp[i] = phases[i].T; phase_set = 1; }
		else              phase_temp[i] = 0.0;
	}

	// check activation
	if(!phase_set) PetscFunctionReturn(0);

	for(imark = 0; imark < nummark; imark++)
	{
		// get current marker
		P = &actx->markers[imark];

		// assign phase temperature to markers, if initial phase temperature is set
		if(phase_temp[P->phase]) P->T = phase_temp[P->phase];
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSetTempFile"
PetscErrorCode ADVMarkSetTempFile(AdvCtx *actx, FB *fb)
{
	FDSTAG         *fs;
	int            fd;
	Marker         *P;
	PetscViewer    view_in;
	PetscLogDouble t;
	char           filename[_str_len_];
	PetscScalar    header, dim[3];
	PetscInt       GridSize, nx, ny, nz, imark, nummark, nmarkx, nmarky, nmarkz;
	PetscScalar    DX, DY, DZ, bx, by, bz, ex, ey, ez;
	PetscScalar    xp, yp, zp, Xc, Yc, Zc, xpL, ypL, zpL;
	PetscScalar    *Temp;
	PetscInt       Ix, Iy, Iz;
	PetscScalar    chTemp, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "temp_file", filename, NULL); CHKERRQ(ierr);

	// check whether file is provided
	if(!strlen(filename)) PetscFunctionReturn(0);

	PrintStart(&t, "Loading temperature redundantly from", filename);

	// access context
	fs     = actx->fs;
	chTemp = actx->jr->scal->temperature;
	Tshift = actx->jr->scal->Tshift;

	// open and read the file
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header
	ierr = PetscBinaryRead(fd, &header, 1, PETSC_SCALAR); CHKERRQ(ierr);

	// read grid dimensions
	ierr = PetscBinaryRead(fd, &dim, 3,  PETSC_SCALAR); CHKERRQ(ierr);

	// compute grid size
	nx = (PetscInt)dim[0];
	ny = (PetscInt)dim[1];
	nz = (PetscInt)dim[2];
	GridSize = nx * ny * nz;

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)GridSize*sizeof(PetscScalar), &Temp); CHKERRQ(ierr);

	// read entire file
	ierr = PetscBinaryRead(fd, Temp, GridSize, PETSC_SCALAR); CHKERRQ(ierr);

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

	PrintDone(t);

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSetTempVector"
PetscErrorCode ADVMarkSetTempVector(AdvCtx *actx)
{
	FDSTAG         *fs;
	JacRes         *jr;
	Marker         *P;
	PetscInt       sx, sy, sz, nx, ny, jj, ID, I, J, K, II, JJ, KK, AirPhase;
	PetscScalar    *ccx, *ccy, *ccz, ***lT;
	PetscScalar    xc, yc, zc, xp, yp, zp, Ttop;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs   =  actx->fs;
	jr   =  actx->jr;

	// get air phase
	AirPhase = -1;
	Ttop     =  0.0;

	if(actx->surf->UseFreeSurf)
	{
		AirPhase = actx->surf->AirPhase;
		Ttop     = actx->jr->bc->Ttop;
	}

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// cell coordinates
	ccx = fs->dsx.ccoor;
	ccy = fs->dsy.ccoor;
	ccz = fs->dsz.ccoor;

	// access temperature vector
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	// scan all markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get coordinates of cell center
		xc = ccx[I];
		yc = ccy[J];
		zc = ccz[K];

		// map marker on the cells of center grids
		if(xp > xc) { II = I; } else { II = I-1; }
		if(yp > yc) { JJ = J; } else { JJ = J-1; }
		if(zp > zc) { KK = K; } else { KK = K-1; }

		// interpolate temperature on the marker
		P->T = InterpLin3D(lT, II, JJ, KK,  sx, sy, sz, xp, yp, zp, ccx, ccy, ccz);

		// override temperature of air phase
		if(AirPhase != -1 && P->phase == AirPhase) P->T = Ttop;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Specific initialization routines
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFiles"
PetscErrorCode ADVMarkInitFiles(AdvCtx *actx, FB *fb)
{
	int            fd;
	Marker         *P;
	PetscViewer    view_in;
	PetscLogDouble t;
	char           *filename, file[_str_len_];
	PetscScalar    *markbuf, *markptr, header, chTemp, chLen, Tshift, s_nummark;
	PetscInt       imark, nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "mark_load_file", file, "./markers/mdb"); CHKERRQ(ierr);

	PrintStart(&t, "Loading markers in parallel from", file);

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

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitGeom"
PetscErrorCode ADVMarkInitGeom(AdvCtx *actx, FB *fb)
{
	Marker         *P;
	PetscLogDouble  t;
	PetscScalar     chLen, chTime, chTemp;
	char            TemperatureStructure[_str_len_];
	PetscInt        jj, ngeom, imark, maxPhaseID;
	GeomPrim        geom[_max_geom_], *pgeom[_max_geom_], *sphere, *box, *hex, *layer, *cylinder;

	// map container to sort primitives in the order of appearance
	map<PetscInt, GeomPrim*> cgeom;
	map<PetscInt, GeomPrim*>::iterator it, ie;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ngeom      = 0;
	maxPhaseID = actx->dbm->numPhases - 1;
	chLen      = actx->jr->scal->length;
	chTime     = actx->jr->scal->time;
	chTemp     = actx->jr->scal->temperature;

	// clear storage
	ierr = PetscMemzero(geom,  sizeof(GeomPrim) *(size_t)_max_geom_); CHKERRQ(ierr);
	ierr = PetscMemzero(pgeom, sizeof(GeomPrim*)*(size_t)_max_geom_); CHKERRQ(ierr);

	PrintStart(&t, "Reading geometric primitives", NULL);

	//=======
	// LAYERS
	//=======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<LayerStart>", "<LayerEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		fb->ID  = jj;								// allows command-line parsing
		GET_GEOM(layer, geom, ngeom, _max_geom_);

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &layer->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "top",    &layer->top,    1, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bottom", &layer->bot,    1, chLen);      CHKERRQ(ierr);

		// optional sinusoidal perturbation of layer interface:
		//  (adds amplitude*sin(2*pi/wavelength*x) to the interface)
		layer->cosine = 0;
		ierr = getIntParam   (fb, _OPTIONAL_, "cosine",  &layer->cosine,  1, maxPhaseID); CHKERRQ(ierr);
		if (layer->cosine==1){
			ierr = getScalarParam   (fb, _REQUIRED_, "wavelength",  &layer->wavelength,  1, maxPhaseID); CHKERRQ(ierr);
			ierr = getScalarParam   (fb, _REQUIRED_, "amplitude",   &layer->amplitude,   1, maxPhaseID); CHKERRQ(ierr);
		}

		// random noise
		layer->rand_amplitude = 0.0;
		ierr = getScalarParam   (fb, _OPTIONAL_, "rand_ampl",  &layer->rand_amplitude,  1, maxPhaseID); CHKERRQ(ierr);

		// Optional temperature options:
		layer->setTemp = 0;
		ierr = getStringParam(fb, _OPTIONAL_, "Temperature",     TemperatureStructure,       NULL ); CHKERRQ(ierr);
		if 		(!strcmp(TemperatureStructure, "constant"))	    {layer->setTemp=1;}
		else if (!strcmp(TemperatureStructure, "linear"))	    {layer->setTemp=2;}
		else if (!strcmp(TemperatureStructure, "halfspace"))    {layer->setTemp=3;}
		
		// Depending on temperature options, get required input parameters
		if (layer->setTemp==1){
			ierr = getScalarParam(fb, _REQUIRED_, "cstTemp", 	&layer->cstTemp, 1, 1);     CHKERRQ(ierr); 
		
			// take potential shift C->K into account	
			layer->cstTemp = (layer->cstTemp +  actx->jr->scal->Tshift)/actx->jr->scal->temperature; 		
		}
		if (layer->setTemp>1){
			ierr = getScalarParam(fb, _REQUIRED_, "topTemp", 	&layer->topTemp, 1, 1);     CHKERRQ(ierr); 
			ierr = getScalarParam(fb, _REQUIRED_, "botTemp", 	&layer->botTemp, 1, 1);     CHKERRQ(ierr); 

			// take potential shift C->K into account	
			layer->topTemp = (layer->topTemp +  actx->jr->scal->Tshift)/actx->jr->scal->temperature; 		
			layer->botTemp = (layer->botTemp +  actx->jr->scal->Tshift)/actx->jr->scal->temperature;
		}
		if (layer->setTemp==3){
			ierr 			= getScalarParam(fb, _REQUIRED_, "thermalAge", &layer->thermalAge, 1, chTime); CHKERRQ(ierr); 
			layer->kappa    = 1e-6/( (actx->jr->scal->length_si)*(actx->jr->scal->length_si)/(actx->jr->scal->time_si)); // thermal diffusivity in m2/s	
		}

		layer->setPhase = setPhaseLayer;

		cgeom.insert(make_pair(fb->blBeg[fb->blockID++], layer));
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//========
	// SPHERES
	//========

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<SphereStart>", "<SphereEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		fb->ID  = jj;								// allows command-line parsing
		
		GET_GEOM(sphere, geom, ngeom, _max_geom_);
		
		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &sphere->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "radius", &sphere->radius, 1, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "center",  sphere->center, 3, chLen);      CHKERRQ(ierr);

		// Optional temperature options:
		sphere->setTemp = 0;
		ierr = getStringParam(fb, _OPTIONAL_, "Temperature",     TemperatureStructure,       NULL ); CHKERRQ(ierr);
		if 		(!strcmp(TemperatureStructure, "constant"))	    {sphere->setTemp=1;}
		
		// Depending on temperature options, get required input parameters
		if (sphere->setTemp==1){
			ierr = getScalarParam(fb, _REQUIRED_, "cstTemp", 	&sphere->cstTemp, 1, 1);     CHKERRQ(ierr); 
		
			// take potential shift C->K into account	
			sphere->cstTemp = (sphere->cstTemp +  actx->jr->scal->Tshift)/actx->jr->scal->temperature; 		
		}
		
		sphere->setPhase = setPhaseSphere;

		cgeom.insert(make_pair(fb->blBeg[fb->blockID++], sphere));
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//======
	// BOXES
	//======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<BoxStart>", "<BoxEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		fb->ID  = jj;								// allows command-line parsing
		GET_GEOM(box, geom, ngeom, _max_geom_);

		box->setTemp = 0;	//default is no	
		ierr = getIntParam   (fb, _REQUIRED_, "phase",  	&box->phase,   1, maxPhaseID); 	CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "bounds",  	 box->bounds,  6, chLen);      	CHKERRQ(ierr);
		box->bot = box->bounds[4]; box->top = box->bounds[5];

		// Optional temperature options:
		box->setTemp = 0;
		ierr = getStringParam(fb, _OPTIONAL_, "Temperature",        TemperatureStructure,       NULL );          CHKERRQ(ierr);
		if 		(!strcmp(TemperatureStructure, "constant"))	    {box->setTemp=1;}
		else if (!strcmp(TemperatureStructure, "linear"))	    {box->setTemp=2;}
		else if (!strcmp(TemperatureStructure, "halfspace"))    {box->setTemp=3;}
		
		// Depending on temperature options, get required input parameters
		if (box->setTemp==1){
			ierr = getScalarParam(fb, _REQUIRED_, "cstTemp", 	&box->cstTemp, 1, chTemp);     CHKERRQ(ierr); 
		}
		if (box->setTemp>1){
			ierr = getScalarParam(fb, _REQUIRED_, "topTemp", 	&box->topTemp, 1, chTemp);     CHKERRQ(ierr); 
			ierr = getScalarParam(fb, _REQUIRED_, "botTemp", 	&box->botTemp, 1, chTemp);     CHKERRQ(ierr); 
		}
		if (box->setTemp==3){
			
			ierr = getScalarParam(fb, _REQUIRED_, "thermalAge", &box->thermalAge, 1, chTime); CHKERRQ(ierr); 

			box->kappa      = 1e-6/( (actx->jr->scal->length_si)*(actx->jr->scal->length_si)/(actx->jr->scal->time_si)); // thermal diffusivity in m2/s	
		}
		
		box->setPhase = setPhaseBox;

		cgeom.insert(make_pair(fb->blBeg[fb->blockID++], box));
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//======
	// HEXES
	//======

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<HexStart>", "<HexEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		fb->ID  = jj;								// allows command-line parsing
		GET_GEOM(hex, geom, ngeom, _max_geom_);

		ierr = getIntParam   (fb, _REQUIRED_, "phase",  &hex->phase, 1,  maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "coord",   hex->coord, 24, chLen);      CHKERRQ(ierr);

		// compute bounding box
		HexGetBoundingBox(hex->coord, hex->bounds);

		hex->setPhase = setPhaseHex;

		cgeom.insert(make_pair(fb->blBeg[fb->blockID++], hex));
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//==========
	// CYLINDERS
	//==========

	ierr = FBFindBlocks(fb, _OPTIONAL_, "<CylinderStart>", "<CylinderEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		fb->ID  = jj;								// allows command-line parsing
		GET_GEOM(cylinder, geom, ngeom, _max_geom_);

		ierr = getIntParam   (fb, _REQUIRED_, "phase",   &cylinder->phase,  1, maxPhaseID); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "radius",  &cylinder->radius, 1, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "base",     cylinder->base,   3, chLen);      CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "cap",      cylinder->cap,    3, chLen);      CHKERRQ(ierr);

		// Optional temperature options:
		cylinder->setTemp = 0;
		ierr = getStringParam(fb, _OPTIONAL_, "Temperature",     TemperatureStructure,       NULL ); CHKERRQ(ierr);
		if 		(!strcmp(TemperatureStructure, "constant"))	    {cylinder->setTemp=1;}
		
		// Depending on temperature options, get required input parameters
		if (cylinder->setTemp==1){
			ierr = getScalarParam(fb, _REQUIRED_, "cstTemp", 	&cylinder->cstTemp, 1, 1);     CHKERRQ(ierr); 
		
			// take potential shift C->K into account	
			cylinder->cstTemp = (cylinder->cstTemp +  actx->jr->scal->Tshift)/actx->jr->scal->temperature; 		
		}

		cylinder->setPhase = setPhaseCylinder;

		cgeom.insert(make_pair(fb->blBeg[fb->blockID++], cylinder));
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// store pointers to primitives in the order of appearance in the file
	for(it = cgeom.begin(), ie = cgeom.end(), ngeom = 0; it != ie; it++)
	{
		pgeom[ngeom++] = it->second;
	}

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
			pgeom[jj]->setPhase(pgeom[jj], P);
		}
	}

	PrintDone(t);

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
	char           filename[_str_len_];
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
	PetscLogDouble t;
	PetscRandom    rctx;
	PetscScalar    cf_rand;
	PetscInt       nPoly;
	PetscScalar    atol;
	PetscScalar    box[4];
	CtrlP          CtrlPoly;
	PetscInt       VolID, nCP;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get file name
	ierr = getStringParam(fb, _OPTIONAL_, "poly_file", filename, "./input/poly.dat"); CHKERRQ(ierr);

	PrintStart(&t, "Loading polygons redundantly from", filename);

	// initialize
	fs = actx->fs;
	x = y = z = dx = dy = dz = 0.0;
	chLen = actx->jr->scal->length;

	// initialize the random number generator
	if(actx->randNoise)
	{
		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);
	}

	//===========================
	// --- initialize markers ---
	//===========================
	
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

	//===============================
	// --- local grid/marker info ---
	//===============================

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

	// read geometry variations
	ierr = ADVMarkReadCtrlPoly(fb, &CtrlPoly, VolID, nCP); CHKERRQ(ierr);

	// --- loop over all volumes ---
	for(kvol = 0; kvol < VolN; kvol++)
	{
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

		// get position of polygons (PetscScalar !)
		for(kpoly = 0; kpoly < Poly.num; kpoly++)
		{
			PolyIdx[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// get lengths of polygons (PetscScalar !)
		for (kpoly=0; kpoly<Poly.num;kpoly++)
		{
			PolyLen[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// interpolate stretch parameters
		PetscScalar SxAll[Poly.num];
		PetscScalar SyAll[Poly.num];
		if (kvol == VolID)
		{
			PetscPrintf(PETSC_COMM_WORLD,"\nVarying volume %d (phase: %d, type: %d) \n", VolID, Poly.phase, Poly.type);
			
			// shift index of control polys by 1 to be in line with c indexing
			PetscInt    i;
    		for (i=0; i < nCP; ++i)
    		{
				// also check if control polygon is out of bounds
				if (CtrlPoly.Pos[i] > Poly.num)
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Control Polygon out of bounds. Volume only has %d polygons", Poly.num);
				}
				PetscPrintf(PETSC_COMM_WORLD,"CtrlPoly %d: Pos: %d, Sx: %.6f, Sy: %.6f \n",i+1,CtrlPoly.Pos[i],CtrlPoly.Sx[i],CtrlPoly.Sy[i]);
				CtrlPoly.Pos[i] = CtrlPoly.Pos[i] - 1;
    		}

			// interpolate stretch parameters
			interpStretch(CtrlPoly.Sx,CtrlPoly.Sy,nCP,CtrlPoly.Pos,Poly.num,SxAll,SyAll);
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

				// vary Polygon geometry
				if (kvol == VolID)
				{
					stretchPolygon(Poly.X,Poly.len,SxAll[kpoly],SyAll[kpoly]);
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

	PrintDone(t);

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkReadCtrlPoly"
PetscErrorCode ADVMarkReadCtrlPoly(FB *fb, CtrlP *CtrlPoly, PetscInt &VolID, PetscInt &nCP)
{
	PetscInt       jj;
	
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// find blocks
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<vG_ControlPolyStart>", "<vG_ControlPolyEnd>"); CHKERRQ(ierr);
	nCP  = fb->nblocks;

	// check number of control polygons
	if (nCP > _max_ctrl_poly_)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "%d exceeds maximum number of control polygons (%d) \n",nCP,_max_ctrl_poly_);
	}

	// loop over blocks
	for(jj = 0; jj < nCP; jj++)
	{
		ierr = getIntParam   (fb, _REQUIRED_, "PolyID",  &CtrlPoly->ID[jj],    1, 0);   CHKERRQ(ierr);
		ierr = getIntParam   (fb, _REQUIRED_, "VolID",   &CtrlPoly->VolID[jj], 1, 0);   CHKERRQ(ierr);
		ierr = getIntParam   (fb, _REQUIRED_, "PolyPos", &CtrlPoly->Pos[jj],   1, 0);   CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "Sx",      &CtrlPoly->Sx[jj],    1, 1);   CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "Sy",      &CtrlPoly->Sy[jj],    1, 1);   CHKERRQ(ierr);

		if (CtrlPoly->VolID[jj] != CtrlPoly->VolID[0])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "All control polygons should have the same volume ID \n");
		}

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	if (nCP > 0)
	{
		VolID = CtrlPoly->VolID[0];
	}

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
// get the density from a phase diagram
#undef __FUNCT__
#define __FUNCT__ "LoadPhaseDiagram"
PetscErrorCode LoadPhaseDiagram(AdvCtx *actx, Material_t  *phases, PetscInt i)
{
	FILE          *fp;
    PetscInt       i_pd,j,ij,lineStart,n,found, NumberOfPhaseDiagramProperties;
    PetscScalar    fl[2];
    char           buf[1000],name[_str_len_];
    PData         *pd;
    Scaling       *scal;
   
	PetscFunctionBegin;

	scal = actx->jr->scal;
	pd   = actx->jr->Pd;

	found = 0;
	// Get the next empty row in the buffer
	for(j=0; j<_max_num_pd_; j++)
	{
		if(!pd->rho_pdns[0][j])
		{
			found 	= 1;
			i_pd 	= j;
			break;
		}
		else
		{
			found 	= 1;
			// Check if we have this diagram already in the buffer
			for(ij=0; ij<_pd_name_sz_; ij++)
			{
				if((pd->rho_pdns[ij][j] != phases[i].pdn[ij]))
				{
					found = 0;
					break;
				}
			}
			if(found == 1)
			{
				// We already loaded that diagram so no need to do anything here except setting the flags for the melt
				sprintf(name,"%s.in",phases[i].pdn);  // is this ever used?
				fp=fopen(phases[i].pdf,"r");
				for(j=0;j<1;j++)
				{
					if(j==0)
					{
						fscanf(fp, "%lf,",&fl[0]);
					}
				}
				if(fl[0] == 4 || fl[0] == 5)
				{
					phases[i].pdAct = 1;
				}
				fclose(fp);
				PetscFunctionReturn(0);
			}
		}
	}

	if(found == 0)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Phase diagram buffer too small!\n\n");
		PetscFunctionReturn(0);
	}

	// Create the name
	sprintf(name,"%s.in",phases[i].pdn);

	lineStart = 50;    // 50 lines are reserved for the header in the phase diagram

	fp=fopen(phases[i].pdf,"r");
	if (fp==NULL)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "No such phase diagram: %s\n",name);
	}

	// Read header
	for(j=0;j<lineStart;j++)
	{
		if(j==0)
		{
			fscanf(fp, "%i,",&pd->numProps[i_pd]);
		}
		else
		{
			fgets(buf, 1000, fp);
		}
	}

	// Read important phase diagram info about the pressure & temperature range of the diagram
	fscanf(fp, "%lf,",&pd->minT[i_pd]);														// minimum T of diagram [in Kelvin]
	pd->minT[i_pd] 			=	pd->minT[i_pd]/scal->temperature;							// non-dimensionalize
	fscanf(fp, "%lf,",&pd->dT[i_pd]);														// Temperature increment
	pd->dT[i_pd] 			=	pd->dT[i_pd]/scal->temperature;								// non-dimensionalize
	fscanf(fp, "%i,",&pd->nT[i_pd]);														// # of temperature points in diagram 
	pd->maxT[i_pd] 	 		=	pd->minT[i_pd] + (PetscScalar)(pd->nT[i_pd])*pd->dT[i_pd];	// maximum T of diagram
	fscanf(fp, "%lf,",&pd->minP[i_pd]);														// minimum P of diagram [in bar]
	pd->minP[i_pd] 			=	(pd->minP[i_pd]*1e5)/scal->stress_si;						// non-dimensionalize
	fscanf(fp, "%lf,",&pd->dP[i_pd]);														// Pressure increment
	pd->dP[i_pd] 			=	(pd->dP[i_pd]*1e5)/scal->stress_si;							// non-dimensionalize
	fscanf(fp, "%i,",&pd->nP[i_pd]);														// # of pressure points in diagram 
	pd->maxP[i_pd] 	 		=	pd->minP[i_pd] + (PetscScalar)(pd->nP[i_pd])*pd->dP[i_pd];	// maximum P of diagram
	
	n = pd->nT[i_pd]*pd->nP[i_pd]; // number of points

	/*
	Check what data is available:
	1 column = rho fluid [kg/m3]
	2 column = melt fraction []
	3 column = density [kg/m3]
	4 column = T [K]
	5 column = P [b]
	*/

	NumberOfPhaseDiagramProperties = pd->numProps[i_pd];
	if (NumberOfPhaseDiagramProperties == 3)  // density
	{
		fscanf(fp,"%lf %lf %lf,",&pd->rho_v[0][i_pd],&fl[0],&fl[1]);
		pd->rho_v[0][i_pd] /= scal->density;
		for (j=1; j<n; j++)
		{
			fscanf(fp, "%lf %lf %lf,",&pd->rho_v[j][i_pd],&fl[0],&fl[1]);
			pd->rho_v[j][i_pd] /= scal->density;
		}
	}
	else if(NumberOfPhaseDiagramProperties == 4)   // density + mf
	{
		fscanf(fp, "%lf %lf %lf %lf,",&pd->Me_v[0][i_pd],&pd->rho_v[0][i_pd],&fl[0],&fl[1]);
		pd->rho_v[0][i_pd] /= scal->density;

		for (j=1; j<n; j++)
		{
			fscanf(fp, "%lf %lf %lf %lf,",&pd->Me_v[j][i_pd],&pd->rho_v[j][i_pd],&fl[0],&fl[1]);
			pd->rho_v[j][i_pd] /= scal->density;
		}
	}
	else if(NumberOfPhaseDiagramProperties == 5)   // density + mf + density_fluid
	{
		fscanf(fp, "%lf %lf %lf %lf %lf,",&pd->rho_f_v[0][i_pd],&pd->Me_v[0][i_pd],&pd->rho_v[0][i_pd],&fl[0],&fl[1]);
		pd->rho_v[0][i_pd] /= scal->density;
		pd->rho_f_v[0][i_pd] /= scal->density;

		for (j=1; j<n; j++)
		{
			fscanf(fp, "%lf %lf %lf %lf %lf,",&pd->rho_f_v[j][i_pd],&pd->Me_v[j][i_pd],&pd->rho_v[j][i_pd],&fl[0],&fl[1]);
			pd->rho_v[j][i_pd] /= scal->density;
			pd->rho_f_v[j][i_pd] /= scal->density;
		}
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"Unknown phase diagram data!\n");
		PetscFunctionReturn(0);
	}

	// Interpolate the name
	for(j=0; j<_pd_name_sz_; j++)
	{
		pd->rho_pdns[j][i_pd] = phases[i].pdn[j];
	}
	fclose(fp);

	// Uncomment to debug values
	// PetscPrintf(PETSC_COMM_WORLD,"RHO = %.20f ; scal = %lf\n 2 = %lf\n  3 = %lf\n 3m = %lf\n  4 = %.20f ; scal = %lf\n 5 = %lf\n 6 = %lf\n 6m = %lf\n n = %i ; scal = %lf\n",pd->rho_v[20000][0], scal.temperature,pd->rho_pdval[1][i_pd],pd->rho_pdval[2][i_pd],pd->rho_pdval[3][i_pd],pd->rho_pdval[4][i_pd], scal.stress_si,pd->rho_pdval[5][i_pd],pd->rho_pdval[6][i_pd],pd->rho_pdval[7][i_pd],n, scal.density);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// geometric primitives functions
//---------------------------------------------------------------------------
void setPhaseSphere(GeomPrim *sphere, Marker *P)
{
	PetscScalar dx, dy, dz;

	dx = P->X[0] - sphere->center[0];
	dy = P->X[1] - sphere->center[1];
	dz = P->X[2] - sphere->center[2];

	if(sqrt(dx*dx + dy*dy + dz*dz) <= sphere->radius)
	{
		P->phase = sphere->phase;
		if (sphere->setTemp>0)
		{	
			PetscScalar T=0;
			computeTemperature(sphere, P, &T);

			P->T = T; 			// set Temperature
		}
	}
}
//---------------------------------------------------------------------------
void setPhaseBox(GeomPrim *box, Marker *P)
{
	if(P->X[0] >= box->bounds[0] && P->X[0] <= box->bounds[1]
	&& P->X[1] >= box->bounds[2] && P->X[1] <= box->bounds[3]
	&& P->X[2] >= box->bounds[4] && P->X[2] <= box->bounds[5])
	{
		P->phase = box->phase;
		if (box->setTemp>0)
		{	
			PetscScalar T=0;
			computeTemperature(box, P, &T);

			P->T = T; 			// set Temperature
		}

	}
}
//---------------------------------------------------------------------------
void setPhaseLayer(GeomPrim *layer, Marker *P)
{
	PetscScalar bot, top,pert,pert_random;


	bot = layer->bot; 
	top = layer->top;
	if (layer->cosine==1){
		// Add sinusoidal perturbation
		pert 	= 	-layer->amplitude*PetscCosScalar(2*PETSC_PI/layer->wavelength*P->X[0]);	
		bot 	= 	bot + pert;
		top 	= 	top + pert;
	}

	// add random noise
	pert_random 	= (rand()/PetscScalar(RAND_MAX)-0.5)*layer->rand_amplitude;
	bot 			= 	bot + pert_random;
	top 			= 	top + pert_random;

	if(P->X[2] >= bot && P->X[2] <= top)
	{
		P->phase = layer->phase;
		if (layer->setTemp>0)
		{	
			PetscScalar T=0;
			computeTemperature(layer, P, &T);

			P->T = T; 			// set Temperature
		}
	}
}
//---------------------------------------------------------------------------
void setPhaseHex(GeomPrim *hex, Marker *P)
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
				return;
			}
		}
	}
}
//---------------------------------------------------------------------------
void setPhaseCylinder(GeomPrim *cylinder, Marker *P)
{
	PetscScalar px, py, pz, ax, ay, az, dx, dy, dz, t;

	// get vector between a test point and cylinder base
	px = P->X[0] - cylinder->base[0];
	py = P->X[1] - cylinder->base[1];
	pz = P->X[2] - cylinder->base[2];

	// get cylinder axis vector
	ax = cylinder->cap[0] - cylinder->base[0];
	ay = cylinder->cap[1] - cylinder->base[1];
	az = cylinder->cap[2] - cylinder->base[2];

	// find normalized parametric coordinate of a point-axis projection
	t = (ax*px + ay*py + az*pz)/(ax*ax + ay*ay + az*az);

	// find distance vector between point and axis
	dx = px - t*ax;
	dy = py - t*ay;
	dz = pz - t*az;

	// check cylinder
	if(t >= 0.0 && t <= 1.0 && sqrt(dx*dx + dy*dy + dz*dz) <= cylinder->radius)
	{
		P->phase = cylinder->phase;

		if (cylinder->setTemp>0)
		{	
			PetscScalar T=0;
			computeTemperature(cylinder, P, &T);

			P->T = T; 			// set Temperature
		}
	}
}
//---------------------------------------------------------------------------
// geometric primitives temperature functions
//---------------------------------------------------------------------------
void computeTemperature(GeomPrim *geom, Marker *P, PetscScalar *T)
{
	// computes the temperature at the point based on the top of the geometric object
	
	if(geom->setTemp == 1)
	{
		// constant temperature
		(*T) = geom->cstTemp;
	}
	else if (geom->setTemp==2)
	{
		// linear temperature between top & bottom
		PetscScalar z_top, z_bot, z, T_top, T_bot;
		
		z_top = geom->top;
		z_bot = geom->bot;
		T_top = geom->topTemp;
		T_bot = geom->botTemp;
		z     = P->X[2];
		(*T)  = (z-z_top)*(T_top - T_bot)/(z_top-z_bot) + T_top; // linear gradient between top & bottom
		

	}
	else if (geom->setTemp==3)
	{
		// Half space cooling profile
		PetscScalar z_top, z, T_top, T_bot, kappa, thermalAge;

		z_top      = geom->top;
		T_top      = geom->topTemp;
		T_bot      = geom->botTemp;
		thermalAge = geom->thermalAge;
		z          = PetscAbs(P->X[2]-z_top);
		kappa      = geom->kappa;
		(*T)       = (T_bot-T_top)*erf(z/2.0/sqrt(kappa*thermalAge)) + T_top;
	}
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
