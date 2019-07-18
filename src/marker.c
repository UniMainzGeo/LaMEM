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
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
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
#include "AVDView.h"
#include "break.h"

/*
#START_DOC#
#END_DOC#
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInit"
PetscErrorCode ADVMarkInit(AdvCtx *actx, UserCtx *user)
{
	FDSTAG   *fs;
	PetscInt  nmarkx, nmarky, nmarkz, nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	PetscPrintf(PETSC_COMM_WORLD," Starting marker initialization routine\n");

	// restart of simulation - markers read from file
	if (user->restart==1) user->msetup = RESTART;

	// allocate storage for uniform distribution
	if(user->msetup != PARALLEL
	&& user->msetup != RESTART)
	{
		// get local number of markers
		nmarkx  = fs->dsx.ncels * user->NumPartX;
		nmarky  = fs->dsy.ncels * user->NumPartY;
		nmarkz  = fs->dsz.ncels * user->NumPartZ;
		nummark = nmarkx*nmarky*nmarkz;

		// allocate storage
		ierr = ADVReAllocStorage(actx, nummark); CHKERRQ(ierr);

		// set number of markers
		actx->nummark = nummark;
	}

	// initialize variables for marker control
	actx->nmin = (PetscInt) (user->NumPartX*user->NumPartY*user->NumPartZ/2); // min no. of markers/cell -50%
	actx->nmax = (PetscInt) (user->NumPartX*user->NumPartY*user->NumPartZ*3);   // max no. of markers/cell 300%
	actx->avdx = user->NumPartX * 3;
	actx->avdy = user->NumPartY * 3;
	actx->avdz = user->NumPartZ * 3;

	// read options from command line
	PetscOptionsGetInt(NULL, NULL ,"-markers_min",   &actx->nmin,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-markers_max",   &actx->nmax,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-markers_avdx",  &actx->avdx,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-markers_avdy",  &actx->avdy,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-markers_avdz",  &actx->avdz,   NULL);

	// initialize coordinates and add random noise if required for hard-coded setups
	if(user->msetup != PARALLEL
	&& user->msetup != REDUNDANT
	&& user->msetup != RESTART
	&& user->msetup != POLYGONS)
	{
		ierr = ADVMarkInitCoord(actx, user); CHKERRQ(ierr);
	}

	// display info on-screen
	PetscPrintf(PETSC_COMM_WORLD," Marker setup employed [msetup] : ");


	// initialize marker phase, temperature, etc.
	if     (user->msetup == PARALLEL)   { PetscPrintf(PETSC_COMM_WORLD,"%s\n","parallel");        ierr = ADVMarkInitFileParallel (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == REDUNDANT)  { PetscPrintf(PETSC_COMM_WORLD,"%s\n","redundant");       ierr = ADVMarkInitFileRedundant(actx, user); CHKERRQ(ierr); }
	else if(user->msetup == POLYGONS)   { PetscPrintf(PETSC_COMM_WORLD,"%s\n","polygons");        ierr = ADVMarkInitFilePolygons (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == DIAPIR)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","diapir");          ierr = ADVMarkInitDiapir       (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == BLOCK)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","block");           ierr = ADVMarkInitBlock        (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == SUBDUCTION) { PetscPrintf(PETSC_COMM_WORLD,"%s\n","subduction");      ierr = ADVMarkInitSubduction   (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == FOLDING)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","folding");         ierr = ADVMarkInitFolding      (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == DETACHMENT) { PetscPrintf(PETSC_COMM_WORLD,"%s\n","detachment");      ierr = ADVMarkInitDetachment   (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == SLAB)       { PetscPrintf(PETSC_COMM_WORLD,"%s\n","slab");            ierr = ADVMarkInitSlab         (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == SPHERES)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","spheres");         ierr = ADVMarkInitSpheres      (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == BANDS)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","bands");           ierr = ADVMarkInitBands        (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == PIPES)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","pipes");           ierr = ADVMarkInitPipes        (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == GEOTH)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","geoth");           ierr = ADVMarkInitGeoth        (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == FAULT)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","fault");           ierr = ADVMarkInitFault        (actx, user); CHKERRQ(ierr); }
	//else if(user->msetup == ROZHKO)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","rozhko");          ierr = ADVMarkInitRozhkoComplex       (actx, user); CHKERRQ(ierr); }
	//else if(user->msetup == ROZHKO)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","rozhko");          ierr = ADVMarkInitRozhko       (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == KM8)     { PetscPrintf(PETSC_COMM_WORLD,"%s\n","km8");          ierr = ADVMarkInitKM8       (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == CON)        { PetscPrintf(PETSC_COMM_WORLD,"%s\n","con");             ierr = ADVMarkInitCon       (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == PREFRAC)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","prefrac");         ierr = ADVMarkInitPrefrac       (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == DOMES)      { PetscPrintf(PETSC_COMM_WORLD,"%s\n","domes");           ierr = ADVMarkInitDomes        (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == ROTATION)   { PetscPrintf(PETSC_COMM_WORLD,"%s\n","rotation");        ierr = ADVMarkInitRotation     (actx, user); CHKERRQ(ierr); }
	else if(user->msetup == RESTART)    { PetscPrintf(PETSC_COMM_WORLD,"%s\n","restart");         ierr = BreakReadMark           (actx      ); CHKERRQ(ierr); }
	else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER," *** Incorrect option for initialization of markers");


	// compute host cells for all the markers
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	// check marker distribution
	ierr = ADVMarkCheckMarkers(actx); CHKERRQ(ierr);

	// project initial history from markers to grid
	ierr = ADVProjHistMarkToGrid(actx); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished marker initialization routine\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitCoord"
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, UserCtx *user)
{
	//=============================================================
	// Initializes coordinates and adds random noise if required for hard-coded setups
	// HANDLES VARIABLE MESH SPACING!
	// WARNING! Random noise only for internal setups
	//==============================================================

	// generate coordinates of uniformly distributed markers within a cell
	FDSTAG      *fs;
	PetscScalar  x, y, z, dx, dy, dz;
	PetscInt     i, j, k, ii, jj, kk;
	PetscInt     imark;
	PetscRandom  rctx;
	PetscScalar  cf_rand;
	PetscBool    AddRandomNoiseParticles;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// random noise
	AddRandomNoiseParticles = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL, NULL,"-AddRandomNoiseParticles", &AddRandomNoiseParticles, NULL); CHKERRQ(ierr);

	if(AddRandomNoiseParticles==PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD, " Adding random noise to marker distribution \n");

	// initialize the random number generator
	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

	// marker counter
	imark = 0;

	// create uniform distribution of markers/cell for variable grid
	for(k = 0; k < fs->dsz.ncels; k++)
	{
		// spacing of particles
		dz = (fs->dsz.ncoor[k+1]-fs->dsz.ncoor[k])/(PetscScalar)user->NumPartZ;
		for(j = 0; j < fs->dsy.ncels; j++)
		{
			// spacing of particles
			dy = (fs->dsy.ncoor[j+1]-fs->dsy.ncoor[j])/(PetscScalar)user->NumPartY;
			for(i = 0; i < fs->dsx.ncels; i++)
			{
				// spacing of particles
				dx = (fs->dsx.ncoor[i+1]-fs->dsx.ncoor[i])/(PetscScalar)user->NumPartX;

				// loop over markers in cells
				for (kk = 0; kk < user->NumPartZ; kk++)
				{
					if (kk == 0) z = fs->dsz.ncoor[k] + dz*0.5;
					else         z = fs->dsz.ncoor[k] + dz*0.5 + (PetscScalar)kk*dz;

					for (jj = 0; jj < user->NumPartY; jj++)
					{
						if (jj == 0) y = fs->dsy.ncoor[j] + dy*0.5;
						else         y = fs->dsy.ncoor[j] + dy*0.5 + (PetscScalar)jj*dy;

						for (ii = 0; ii < user->NumPartX; ii++)
						{
							if (ii == 0) x = fs->dsx.ncoor[i] + dx*0.5;
							else         x = fs->dsx.ncoor[i] + dx*0.5 + (PetscScalar)ii*dx;

							// set marker coordinates
							actx->markers[imark].X[0] = x;
							actx->markers[imark].X[1] = y;
							actx->markers[imark].X[2] = z;

							if(AddRandomNoiseParticles)
							{
								// add random noise
								// decrease/increase amount of noise by changing A in: (cf_rand-0.5)*dx/A
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[0] += (cf_rand-0.5)*dx/1;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[1] += (cf_rand-0.5)*dy/1;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								actx->markers[imark].X[2] += (cf_rand-0.5)*dz/1;
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
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSave"
PetscErrorCode ADVMarkSave(AdvCtx *actx, UserCtx *user)
{
	int          fd;
	PetscInt     imark;
	Marker      *P;
	char        *SaveFileName;
	PetscViewer  view_out;
	PetscScalar *markbuf, *markptr, header, chLen, chTemp, Tshift, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!user->SaveParticles) PetscFunctionReturn(0);

	PetscPrintf(PETSC_COMM_WORLD," Saving markers in parallel to files: ./%s/%s.xxx.out \n",
		user->SaveInitialParticlesDirectory, user->ParticleFilename);

	// initialize file header for MATLAB compatibility
	header = -1;

	// create write buffer
	ierr = PetscMalloc((size_t)(5*actx->nummark)*sizeof(PetscScalar), &markbuf); CHKERRQ(ierr);

	chLen  = actx->jr->scal.length;
	chTemp = actx->jr->scal.temperature;

	// temperature shift
	Tshift = actx->jr->scal.Tshift;

	// copy data from storage into buffer
	for(imark = 0, markptr = markbuf; imark < actx->nummark; imark++, markptr += 5)
	{
		P          =              &actx->markers[imark];
		markptr[0] =              P->X[0]*chLen;
		markptr[1] =              P->X[1]*chLen;
		markptr[2] =              P->X[2]*chLen;
		markptr[3] = (PetscScalar)P->phase;
		markptr[4] =             (P->T - Tshift)*chTemp;
	}

	// create directory
	ierr = LaMEMCreateOutputDirectory(user->SaveInitialParticlesDirectory); CHKERRQ(ierr);

	// compile file name
	asprintf(&SaveFileName, "./%s/%s.%lld.out",
		user->SaveInitialParticlesDirectory,
		user->ParticleFilename, (LLD)actx->iproc);

	// open file for binary output
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, SaveFileName, FILE_MODE_WRITE, &view_out); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_out, &fd);                                    CHKERRQ(ierr);

	// write binary output
	s_nummark = (PetscScalar)actx->nummark;
	ierr = PetscBinaryWrite(fd, &header,    1,               PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fd, &s_nummark, 1,               PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fd, markbuf,    5*actx->nummark, PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);

	// destroy file handle & file name
	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
	free(SaveFileName);

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
	PetscBool    error;
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
	error = PETSC_FALSE;

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
		error = PETSC_TRUE;
	}

	if(numNonLocal)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of non-local markers: %lld\n", (LLD)numNonLocal); CHKERRQ(ierr);
		error = PETSC_TRUE;
	}

	if(numEmpty)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of exactly empty cells: %lld\n", (LLD)numEmpty); CHKERRQ(ierr);
		error = PETSC_TRUE;
	}

	if(numSparse)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "WARNING! Number of cells with less than 8 markers: %lld.\n", (LLD)numSparse); CHKERRQ(ierr);
	}

	if(error == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Problems with initial marker distribution (see the above message)");
	}

	// clear
	ierr = PetscFree(numMarkCell); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Specific initialization routines
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFileParallel"
PetscErrorCode ADVMarkInitFileParallel(AdvCtx *actx, UserCtx *user)
{
	//=============================================================
	// WARNING! Random noise only for internal setups!!
	// Markers from file are assumed to have the coordinates set
	//==============================================================

	int          fd;
	Marker      *P;
	PetscViewer  view_in;
	char        *LoadFileName;
	PetscScalar *markbuf, *markptr, header, chTemp, chLen, Tshift, s_nummark;
	PetscInt     imark, nummark, NumPropsPerMarker;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD," Loading markers in parallel from files: ./%s/%s.xxx.out \n",
		user->LoadInitialParticlesDirectory, user->ParticleFilename);

	// compile input file name
	asprintf(&LoadFileName, "./%s/%s.%lld.out",
		user->LoadInitialParticlesDirectory,
		user->ParticleFilename, (LLD)actx->iproc);

	// open file
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, LoadFileName, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd);                                   CHKERRQ(ierr);

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
	free(LoadFileName);

	// get characteristic length & temperature
	chLen  = actx->jr->scal.length;
	chTemp = actx->jr->scal.temperature;

	// temperature shift
	Tshift = actx->jr->scal.Tshift;

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
	// WARNING! NOT SURE WHETHER THIS IS REALLY NECESSARY
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished Loading markers in parallel \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFileRedundant"
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, UserCtx *user)
{
	//=============================================================
	// HANDLES VARIABLE MESH SPACING!
	// WARNING! Random noise only for internal setups!!
	// Markers from file are assumed to have the coordinates set
	//==============================================================

	// REDUNDANTLY loads the X,Y,Z-coord, phase and temperature data from a file
	// each processor uses its own part after loading, and ignores the rest
	// THEREFORE! not to be used for large markers files - use parallel read for that!
	// coord and temperature are assumed to be dimensional

	FDSTAG       *fs;
	int           fd;
	Marker       *P;
	PetscViewer   view_in;
	char         *LoadFileName;
	PetscScalar   maxtemp, header, chLen, chTemp, Tshift;
	PetscScalar   x, y, z, xs[3], xe[3], info[3], *phase, *temp, *xcoor, *ycoor, *zcoor;
	PetscInt      nmarkx, nmarky, nmarkz, nmark;
	PetscInt      tmarkx, tmarky, tmarkz, tmark;
	PetscInt      i, imark, maxphase;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;

	// compile input file name
	asprintf(&LoadFileName, "./%s/%s.dat",
		user->LoadInitialParticlesDirectory,
		user->ParticleFilename);

	PetscPrintf(PETSC_COMM_WORLD," Loading markers redundantly from file: %s \n", LoadFileName);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, LoadFileName, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header
	ierr = PetscBinaryRead(fd, &header, 1, PETSC_SCALAR); CHKERRQ(ierr);

	// read number of markers in each direction
	ierr = PetscBinaryRead(fd, info, 3, PETSC_SCALAR); CHKERRQ(ierr);

	nmarkz = (PetscInt)(info[0]);
	nmarky = (PetscInt)(info[1]);
	nmarkx = (PetscInt)(info[2]);
	nmark  = nmarkx*nmarky*nmarkz;

	// compute total number of markers
	tmarkx = fs->dsx.tcels*user->NumPartX;
	tmarky = fs->dsy.tcels*user->NumPartY;
	tmarkz = fs->dsz.tcels*user->NumPartZ;
	tmark  = tmarkx*tmarky*tmarkz;

	// check whether file contains the expected number of markers
	if(tmark != nmark)
	{
		SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER," The file does not contain the expected number of markers [expected = %lld vs. present = %lld]", (LLD)tmark, (LLD)nmark);
	}

	// create buffer for x,y,z, phase & temperature
	ierr = PetscMalloc((size_t)nmark*sizeof(PetscScalar), &xcoor); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nmark*sizeof(PetscScalar), &ycoor); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nmark*sizeof(PetscScalar), &zcoor); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nmark*sizeof(PetscScalar), &phase); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nmark*sizeof(PetscScalar), &temp);  CHKERRQ(ierr);

	// read phase and temperature into buffer
	ierr = PetscBinaryRead(fd, xcoor, nmark, PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, ycoor, nmark, PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, zcoor, nmark, PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, phase, nmark, PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, temp , nmark, PETSC_SCALAR); CHKERRQ(ierr);

	// destroy file-handle & file name
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
	free(LoadFileName);

	// fill markers array
	imark    = 0;
	maxphase = 0;
	maxtemp  = 0.0;

	// coordinates of the local domain
	xs[0] = fs->dsx.ncoor[0];
	xs[1] = fs->dsy.ncoor[0];
	xs[2] = fs->dsz.ncoor[0];

	xe[0] = fs->dsx.ncoor[fs->dsx.ncels];
	xe[1] = fs->dsy.ncoor[fs->dsy.ncels];
	xe[2] = fs->dsz.ncoor[fs->dsz.ncels];

	// get characteristic length & temperature
	chLen  = actx->jr->scal.length;
	chTemp = actx->jr->scal.temperature;

	// temperature shift
	Tshift = actx->jr->scal.Tshift;

	// loop over all markers and put them in the local domain
	for(i = 0; i < nmark; i++)
	{
		x  = xcoor[i]/chLen;
		y  = ycoor[i]/chLen;
		z  = zcoor[i]/chLen;

		// truncate non-local markers
		if(x >= xs[0] && x <= xe[0]
		&& y >= xs[1] && y <= xe[1]
		&& z >= xs[2] && z <= xe[2])
		{
			P        = &actx->markers[imark];
			P->X[0]  = x;
			P->X[1]  = y;
			P->X[2]  = z;
			P->phase = (PetscInt)phase[i];
			P->T     = (temp[i] + Tshift)/chTemp;

			// increment local counter
			imark++;
		}

		// calculate max no of phases
		if((PetscInt)phase[i] > maxphase)
		{
			maxphase = (PetscInt) phase[i];
		}

		// calculate max temperature
		if(temp[i] > maxtemp)
		{
			maxtemp = temp[i];
		}
	}

	// error checking
	if((maxphase + 1) != actx->jr->numPhases)
	{
		SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER," No. of detected phases %lld does not correspond to the no. of phases given in parameters file %lld", (LLD)maxphase, (LLD)actx->jr->numPhases);
	}

	PetscPrintf(PETSC_COMM_WORLD," Statistics markers: MaxPhase = %lld, MaxTemp = %g, NumberPhases = %lld \n", (LLD)maxphase, maxtemp, (LLD)actx->jr->numPhases);

	// free specific allocated memory
	ierr = PetscFree(xcoor); CHKERRQ(ierr);
	ierr = PetscFree(ycoor); CHKERRQ(ierr);
	ierr = PetscFree(zcoor); CHKERRQ(ierr);
	ierr = PetscFree(phase); CHKERRQ(ierr);
	ierr = PetscFree(temp);  CHKERRQ(ierr);

	// wait until all processors finished reading markers
	// WARNING! NOT SURE WHETHER THIS IS REALLY NECESSARY
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished Loading markers redundantly \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitDiapir"
PetscErrorCode ADVMarkInitDiapir(AdvCtx *actx, UserCtx *user)
{
	PetscInt imark;

	PetscFunctionBegin;

	// Diapir setup - not done - should have some perturbation

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitBlock"
PetscErrorCode ADVMarkInitBlock(AdvCtx *actx, UserCtx *user)
{
	// falling block setup

	PetscInt    imark, nel_x, nel_y, nel_z;
	PetscScalar dx,dy,dz;
	PetscScalar bleft, bright, bfront, bback, bbottom, btop;
	PetscScalar blx, bly, blz;
	PetscBool   b2D, b2Dy;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		actx->markers[imark].p     = 0.0;

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSubduction"
PetscErrorCode ADVMarkInitSubduction(AdvCtx *actx, UserCtx *user)
{
	// subduction setup with air

	PetscScalar H_air, SlabThickness, SlabWidth, SlabLength;
	PetscScalar DistanceFromLeft, SlabMaxSubdDepth,Slab_Xright, dz;
	PetscScalar Air_thickness, Slab_ThicknessFactor, Air_ThicknessFactor, Slab_WidthFactor;
	PetscScalar chLen;
	PetscInt    imark, nz;

	PetscFunctionBegin;

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFolding"
PetscErrorCode ADVMarkInitFolding(AdvCtx *actx, UserCtx *user)
{
	// multilayer folding setup (Zagros)

	PetscInt    i, imark;
	PetscBool   flg, DisplayInfo, FoldingSetupDimensionalUnits;
	PetscScalar zbot[10], ztop[10], BottomStepPerturbation_Amplitude;
	PetscScalar chLen;

	PetscFunctionBegin;

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
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitDetachment"
PetscErrorCode ADVMarkInitDetachment(AdvCtx *actx, UserCtx *user)
{
	// 1-layer over detachment with two linear shaped perturbation (Grasemann & Schmalholz 2012)
	// one perturbation is fixed in position, the other one's position is varied with -Heterogeneity_Offset

	PetscInt    imark;
	PetscScalar zbot = 0.0, ztop = 0.0;
	PetscScalar Het_L, Het_W, Het_A, Het_Off;
	PetscScalar chLen;
	PetscBool   flg, DisplayInfo;

	PetscFunctionBegin;

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSlab"
PetscErrorCode ADVMarkInitSlab(AdvCtx *actx, UserCtx *user)
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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitSpheres"
PetscErrorCode ADVMarkInitSpheres(AdvCtx *actx, UserCtx *user)
{
	// multiple falling spheres

	PetscInt    nsphere = 10, ind, imark;
	PetscScalar rsphere = 0.1;
	PetscScalar x = 0.0, y = 0.0, z = 0.0;
	PetscScalar xc[20], yc[20], zc[20];

	if(user) user = NULL;

 	PetscFunctionBegin;

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitBands"
PetscErrorCode ADVMarkInitBands(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		if     (z >  Hb + H)            {
			P->phase = 2; // air
		}
		else if(z >= Hb && z <= Hb + H) {
			P->phase = 1; // matrix
		}
		else                            {
			P->phase = 0; // basal layer
		}

		// check inclusion
		if(use_inc == PETSC_TRUE && z >= Hb && z <= Hb + size
		&& ((y <= length         && x >= -(offset + size)/1.0 && x <= -(offset - size)/1.0)
		||  (y >= user->L-length && x >=  (offset - size)/1.0 && x <=  (offset + size)/1.0))){
			P->phase = 0;
		}

		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitRozhko"
PetscErrorCode ADVMarkInitRozhko(AdvCtx *actx, UserCtx *user)
{
	/*Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		if     (z >  Hb + H)            {
			P->phase = 3; // air
		}
		else if(z >= Hb && z <= Hb + H) {
			P->phase = 2; // matrix
		}
		else                            {
			P->phase = 1; // basal layer
		}

		// check inclusion
		//if(use_inc == PETSC_TRUE && z >= Hb && z <= Hb + length
		//&& (x < -2.0*size/3.0 || x >= size/3.0)){

		if(use_inc == PETSC_TRUE && z >= Hb && z <= Hb + length
					&& (x < -4000.0/151.0/2.0 || x >= 4000.0/151.0/2.0)){
			P->phase = 1;
		}
		//if(use_inc == PETSC_TRUE && z <= Hb + length/2){
		//			P->phase = 1;
		//		}

		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);*/



	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		if     (z >  Hb + H)            {
			P->phase = 2; // air
		}
		else if(z >= Hb && z <= Hb + H) {
			P->phase = 1; // matrix
		}
		else                            {
			P->phase = 2; // basal layer
		}

		// check inclusion
		if(use_inc == PETSC_TRUE && z >= Hb && z <= Hb + size
		&& ((y <= length         && x >= -(offset + 2*size)/2.0 && x <= -(offset - 2*size)/2.0)
		||  (y >= user->L-length && x >=  (offset - 2*size)/2.0 && x <=  (offset + 2*size)/2.0))){
			P->phase = 1;
		}


		if (x >= ((user->W+user->x_left)) - Hb || x <= (-(user->W+user->x_left)) + Hb) {
			P->phase = 2; // low permeability
		}

		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);



}
//---------------------------------------------------------------------------

/*//---------------------------------------------------------------------------
// Para ComplexSetup5
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitRozhkoComplex"
PetscErrorCode ADVMarkInitRozhkoComplex(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift, P1x, P1z, P2x, P2z, faultWidth, Fshift;

	// initialize the random number generator
	PetscRandom   rctx;
	PetscScalar   cf_rand;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);

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

		P->phase = 1;

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		// Fault
		faultWidth = 50.0;
		P1x = 0.0; //5.0*user->W/8.0 - user->x_left;
		P1z = 0.0;
		P2x = 2100.0; //user->W-user->x_left;
		P2z = 3000.0 ; //user->H/2.0;
		if  (  z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
		{
			P->phase = 6; // core
		}
		else if  ( (  z >= faultWidth +  P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < 3.0*faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
				|| (  z >= -2.0*faultWidth +  P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) ))
		{
			P->phase = 7; // 7 damaged area of the fault
		}
		else {


			if  (  z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) ) Fshift = 200.0;
			else Fshift = 0.0;
			// assign phase in layers
			if     (z <  user->H/3.0 - Fshift)
			{
				// generate random phase
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 1;// granite
				else                 P->phase = 8;
			}
			else if(z < user->H/2.0 - Fshift) {
				// generate random phase
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 2;// impermeable
				else                 P->phase = 8;
			}
			else if(z < 4.0*user->H/6.0  - Fshift) {
				// generate random phase
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 3; // aquifer
				else                 P->phase = 8;
			}
			else if(z < 5.0*user->H/6.0  - Fshift) {
				// generate random phase
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 4; // sandstone
				else                  P->phase = 8;
			}
			else {
				// generate random phase
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 5; // sand
				else                  P->phase = 8;
			}
		}

		//pipes
		if ((x<-1050 && x>= -1100 && z>=1500) || (x>=-950 && x<= -900 && z>=1500) ){
			//P->phase = 9;
		}
		//pipes
		if (x>=1050 && x<= 950 && z>=550)  {
			//P->phase = 3;
		}

		// check inclusion
		if(use_inc == PETSC_TRUE && z <= size/2
		&& ((y <= length         && x >= -(offset + size)/1.0 && x <= -(offset - size)/1.0)
		||  (y >= user->L-length && x >=  (offset - size)/1.0 && x <=  (offset + size)/1.0))){
			//P->phase = 0;
		}


		// assign temperature
		P->T = Tshift;
	}

	// destroy random context
	ierr = PetscRandomDestroy(&rctx);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitRozhkoComplex"
PetscErrorCode ADVMarkInitRozhkoComplex(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift, P1x, P1z, P2x, P2z, faultWidth, Fshift;

	// initialize the random number generator
	PetscRandom   rctx;
	PetscScalar   cf_rand;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);

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

		P->phase = 1;

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		// Fault
		faultWidth = 50.0;
		P1x = 0.0; //5.0*user->W/8.0 - user->x_left;
		P1z = 0.0;
		P2x = 2100.0; //user->W-user->x_left;
		P2z = 3000.0 ; //user->H/2.0;
		if  (  z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
		{
			P->phase = 6; // core
		}
		else if  ( (  z >= faultWidth +  P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < 3.0*faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
				|| (  z >= -2.0*faultWidth +  P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) ))
		{
			P->phase = 7; // 7 damaged area of the fault
		}
		else {


			if  (  z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) ) Fshift = 200.0;
			else Fshift = 0.0;
			// assign phase in layers
			if     (z <  user->H/3.0 - Fshift)
			{
				/* generate random phase*/
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 1;// granite
				else                 P->phase = 8;
			}
			else if(z < user->H/2.0 - Fshift) {
				/* generate random phase*/
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 2;// impermeable
				else                 P->phase = 8;
			}
			else if(z < 4.0*user->H/6.0  - Fshift) {
				/* generate random phase*/
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 3; // aquifer
				else                 P->phase = 8;
			}
			else if(z < 5.0*user->H/6.0  - Fshift) {
				/* generate random phase*/
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 4; // sandstone
				else                  P->phase = 8;
			}
			else {
				/* generate random phase*/
				ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				if (cf_rand <= 0.75 ) P->phase = 5; // sand
				else                  P->phase = 8;
			}
		}

		//pipes
		if ((x<-1050 && x>= -1100 && z>=1500) || (x>=-950 && x<= -900 && z>=1500) ){
			//P->phase = 9;
		}
		//pipes
		if (x>=1050 && x<= 950 && z>=550)  {
			//P->phase = 3;
		}

		// check inclusion
		if(use_inc == PETSC_TRUE && z <= size/2
		&& ((y <= length         && x >= -(offset + size)/1.0 && x <= -(offset - size)/1.0)
		||  (y >= user->L-length && x >=  (offset - size)/1.0 && x <=  (offset + size)/1.0))){
			//P->phase = 0;
		}


		// assign temperature
		P->T = Tshift;
	}

	// destroy random context
	ierr = PetscRandomDestroy(&rctx);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitKM8"
PetscErrorCode ADVMarkInitKM8(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, x, y, z, scal, Tshift, P1x, P1z, P2x, P2z, faultWidth, Fshift;
	PetscScalar x1,y1,z1,x2,y2,z2,x3,y3,z3, aux, aux2, aux3, aux4;

	// initialize the random number generator
	PetscRandom   rctx;
	PetscScalar   cf_rand;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);           CHKERRQ(ierr);

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

		P->phase = 1;

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		if (x<2900.0 && x>-2900.0 && y<2900.0 && y>-2900.0 && z>100.0 && z< 3500.0)
		{

			//pre existing fractures
			aux2=10e6;
			x1=-2000.0; y1=0.0; z1=960.0;
			x2=-1900.0; y2=0.0; z2=860.0;
			x3=-1900.0; y3=10.0; z3=860.0;
			aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
			if ( aux  < aux2 && aux  > -aux2)
			{
				P->phase = 2;
			}
		}


        /*// faults
		aux2=10e6;

		if (x<2900.0 && x>-2900.0 && y<2900.0 && y>-2900.0 && z>100.0 && z< 3500.0)
		{

			//P2, P7, P8 (KM-ML)
			x1=2400.0; y1=0.0; z1=960.0;
			x2=500.0; y2=-2500.0; z2=960.0;
			x3=500.0; y3=0.0; z3=2400.0;
			aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
			if ( aux  < aux2 && aux  > -aux2)
			{
				P->phase = 3;
			}


			//P4, P7, P9 (KMWB)
			if (z<1900.0)
			{
			x1=-2550.0; y1=0.0; z1=960.0;
			x2=500.0; y2=-2600.0; z2=960.0;
			x3=-2600.0; y3=0.0; z3=1345.0;
			//P13, P7, P9
			x1=-2150.0; y1=0.0; z1=0.0;
			x2=500.0; y2=-2600.0; z2=960.0;
			x3=-2700.0; y3=0.0; z3=1345.0;
			aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
			if ( aux  < aux2/1.0 && aux  > -aux2/1.0)
			{
				P->phase = 4;
			}
			}

			//P1, P5, P10 (kme)
			if (z<1237.0)
			{
			x1=1200.0; y1=0.0; z1=960.0;
			x2=0.0; y2=1500.0; z2=960.0;
			x3=700.0; y3=0.0; z3=340.0;
			aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
			if ( aux  < aux2/1.0 && aux  > -aux2/1.0)
			{
				P->phase = 5;
			}
			}

			//P3, P11, P12 (kmw)
			if (z<1900.0)
			{
			x1=-1000.0; y1=0.0; z1=960.0;
			x2=-200.0; y2=-2200.0; z2=960.0;
			x3=-1700.0; y3=0.0; z3=1874.0;
			aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
			if ( aux  < aux2/1.0 && aux  > -aux2/1.0)
			{
				P->phase = 6;
			}
			}



			// assign phase in layers
			if  (   (z <  1050.0) && (z>1020.0))
			{
				//P->phase = 2;
			}
		}
		//Aquifers
		if (   (    (x-400.0)*(x-400.0)  +  (z-(-4000.0))*(z-(-4000.0))  < 6050.0*6050.0)
				&& ((x-400.0)*(x-400.0)  +  (z-(-4000.0))*(z-(-4000.0))  > 6000.0*6000.0)
                && (x>-1700.0) && (x<600.0)
		)
		{

			//P3, P11, P12 (kmw)

						x1=-1000.0; y1=0.0; z1=960.0;
						x2=-200.0; y2=-2200.0; z2=960.0;
						x3=-1700.0; y3=0.0; z3=1874.0;
						aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
						if ( aux  < aux2/1.0 )
						{


							//P2, P7, P8 (KM-ML)
							x1=2400.0; y1=0.0; z1=960.0;
							x2=500.0; y2=-2500.0; z2=960.0;
							x3=500.0; y3=0.0; z3=2400.0;
							aux=(x-x1)*(y2-y1)*(z3-z1)+(y-y1)*(z2-z1)*(x3-x1)+(x2-x1)*(y3-y1)*(z-z1)-(x3-x1)*(y2-y1)*(z-z1)-(y-y1)*(x2-x1)*(z3-z1)-(z2-z1)*(y3-y1)*(x-x1);
							if ( aux  > -aux2/5.0)
							{
								P->phase = 2;
							}

						}


			//P->phase = 2;
		}
		*/

	}

	// destroy random context
	ierr = PetscRandomDestroy(&rctx);    CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitCon"
PetscErrorCode ADVMarkInitCon(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscBool   use_inc;
	PetscScalar H, Hb, size, offset, length, height, x, y, z, scal, Tshift, radius, center, aux;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal   = actx->jr->scal.length;
	Tshift = actx->jr->scal.Tshift;

	// set default values
	H      = 10.0;               // layer thickness               [km]
	Hb     = 2.0;                // basal layer thickness         [km]
	size   = 1.0;                // inclusion size in XZ-plane    [km]
	offset = 0.0;                // inclusion offset along X axis [km]
	length = (user->L*scal)/2.0; // inclusion length along Y axis [km]
	height = (user->H*scal);

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
	height /= scal;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		// assign phase in layers
		P->phase = 2; // air
		//if     ((x*x + (y-length)*(y-length) + z*z) > (height *height))
		aux = 1.1;
		radius =  (z-height/aux)*(z-height/aux)*1.5;
		if     ( ( (z< height/aux)   &&   (  x*x + (y-length)*(y-length)  < radius ) ) ) { //|| (z < height/8.0 )) {
			P->phase = 1; // basal layer
		}

		//Cilinder
		center = user->W/2.0;
		radius = user->W/4.0;

		if     (  (x-center)*(x-center) + (y-length)*(y-length)  < radius*radius   ) { //|| (z < height/8.0 )) {
					P->phase = 2; // basal layer
		}

		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitGeoth"
PetscErrorCode ADVMarkInitGeoth(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscScalar x, y, z;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		y = P->X[1];
		z = P->X[2];

		//if     (    (x >= user->x_left+(user->W/40.0) && x<= user->x_left+2.0*(user->W/40.0) && z>user->H/4.0) || (x <= user->W+user->x_left-(user->W/40.0) && x>= user->W+user->x_left-2.0*(user->W/40.0) && z>3.0*user->H/4.0)   )
		if     (     (x <= user->W+user->x_left-2.0*(user->W/40.0) && x>= user->W+user->x_left-4.0*(user->W/40.0) && z>3.0*user->H/4.0)   )
		{
			P->phase = 2;
		}
		else if ( (z>user->H/5.0) && (z<4.0*user->H/5.0))
		{
			P->phase = 1;
		}
		else
		{
			P->phase = 0;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFault"
PetscErrorCode ADVMarkInitFault(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscScalar P1x, P1z, P2x, P2z, faultWidth, x, z, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	Tshift = actx->jr->scal.Tshift;

	// set default values

	// fault from point P1 to point P2

	/*P1x = (user->W+2.0*user->x_left)/2.0  - (user->W / 5.0) ;
	P1z = user->H/4.0;

	P2x = P1x + 2.0 * (user->W / 5.0);
	P2z = P1z + user->H/2.0;*/

	P1x = -20000.0;
	P1z = 5000.0;

	P2x = -10000;
	P2z = 15000.0;

	faultWidth = user->H / 100.0;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		z = P->X[2];

		// assign phase in layers
		if     (   x >  -10000.0 && x < 10000.0 &&  z >= 5000.0 + 1.0/2.0*(x+10000.0) && z <= 5200.0 + 1.0/2.0*(x+10000.0)     )
		//if     (  x >  P1x && x <  P2x  && z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
		{
			P->phase = 1; //fault
		}
		else
		{
			P->phase = 0;
		}
		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitPipes"
PetscErrorCode ADVMarkInitPipes(AdvCtx *actx, UserCtx *user)
{

	Marker     *P;
	PetscInt    imark;
	PetscScalar P1x, P1z, P2x, P2z, P3x, P3z, P4x, P4z, P5x, P5z, pipeWidth, x, z, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	Tshift = actx->jr->scal.Tshift;

	// set default values

	// pipe by points P1->P2->P3->P4

	P1x = user->x_left+(user->W/5.0);
	P1z = 5.0*user->H/8.0;

	P2x = P1x;
	P2z = user->H/3.0;

	P3x = P1x + (user->W/5.0);
	P3z = P2z;

	P4x =  user->x_left + 4.0*user->W/5.0;
	P4z = 3.0*user->H/4.0;

	P5x = P4x;
	P5z = 5.0*user->H/6.0;

	pipeWidth = user->H / 80.0;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		z = P->X[2];

		if (z<=2.0*user->H/4.0 || z>=3.0*user->H/4.0)
		{
			P->phase = 0;
		}
		else
		{
			P->phase = 1;
		}
		if     (  (x >=  P1x && x <=  P1x+pipeWidth  && z <= P1z && z >= P2z )    || ( x>= P2x && x<= P3x && z >= P2z && z <= P2z +pipeWidth)    ||
			(  x >  P3x && x <  P4x  && z >= P3z + (x-P3x)/(P4x-P3x)*(P4z-P3z) && z < pipeWidth + P3z + (x-P3x)/(P4x-P3x)*(P4z-P3z ) ) ||
			(x >=  P4x && x <=  P4x+pipeWidth  && z >= P4z && z <= P5z ) )
		{
			P->phase = 2;
		}

		// assign temperature
		P->T = Tshift;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitPrefrac"
PetscErrorCode ADVMarkInitPrefrac(AdvCtx *actx, UserCtx *user)
{
	Marker     *P;
	PetscInt    imark;
	PetscScalar P1x, P1z, P2x, P2z, P3x, P3z, P4x, P4z, P5x, P5z, P6x, P6z, faultWidth, x, z, Tshift;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	Tshift = actx->jr->scal.Tshift;

	// set default values

	// fault from point P1 to point P2
	P1x = 0.0; //(user->W+2.0*user->x_left)/3.0  - (user->W / 5.0) ;
	P1z = user->H/4.0;

	P2x = user->W / 2.0 - user->W / 5.0; //P1x + 2.0 * (user->W / 5.0);
	P2z = P1z + user->H/2.0;

	// fault from point P3 to point P4
	P3x = 0.0; //P1x-user->W / 5.0;
	P3z = 3000.0; //user->H/4.0;

	P4x = P2x-user->W / 5.0;
	P4z = P1z + user->H/2.0;

	// fault from point P3 to point P4
	P3x = -300.0;	//0.0; //P1x-user->W / 5.0;
	P3z = 2300.0; //user->H/4.0;

	P4x = 1300.0; 	//P2x-user->W / 5.0;
	P4z = 3600.0; 	//P1z + user->H/2.0;

	P5x = -700.0; //P1x-user->W / 5.0;
	P5z = 2800.0; //user->H/4.0;

	P6x = 500.0;
	P6z = 3800.0;

	faultWidth = user->H / 500.0;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		P = &actx->markers[imark];

		// get coordinates
		x = P->X[0];
		z = P->X[2];

		// assign phase in layers
		//if     (  x >  P1x && x <  P2x  && z >= P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z) && z < faultWidth + P1z + (x-P1x)/(P2x-P1x)*(P2z-P1z ) )
		//{
		//	P->phase = 1; //fault
		//}
		//else
			if (  x >=  P3x && x <  P4x  && z >= P3z + (x-P3x)/(P4x-P3x)*(P4z-P3z) && z < faultWidth + P3z + (x-P3x)/(P4x-P3x)*(P4z-P3z ) )
		{
			P->phase = 2; //fault
		}
		else if (  x >=  P5x && x <  P6x  && z >= P5z + (x-P5x)/(P6x-P5x)*(P6z-P5z) && z < faultWidth + P5z + (x-P5x)/(P6x-P5x)*(P6z-P5z ) )
		{
			P->phase = 2; //fault
		}
		else
		{
			P->phase = 1;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitDomes"
PetscErrorCode ADVMarkInitDomes(AdvCtx *actx, UserCtx *user)
{
	// water phase    -> 0
	// basement phase -> 1
	// salt phase     -> 2

	Marker     *P;
	PetscInt    imark, anhydrite_phase;
	PetscScalar surf_level, anhydrite_bot, anhydrite_top, x, z, scal, base_top;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = actx->jr->scal.length;

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

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitRotation"
PetscErrorCode ADVMarkInitRotation(AdvCtx *actx, UserCtx *user)
{
	// falling block setup

	PetscInt    imark;
	PetscScalar x, z, rad, dist;

	PetscFunctionBegin;

	// print info
	PetscPrintf(PETSC_COMM_WORLD,"  ROTATION SETUP \n");

	// init skip solver
	user->SkipStokesSolver = PETSC_TRUE;

	// circle dimensions
	rad  =  0.15;
	dist = -0.25;

	// loop over local markers
	for(imark = 0; imark < actx->nummark; imark++)
	{
		actx->markers[imark].phase = 0;
		actx->markers[imark].T     = 0.0;

		x = actx->markers[imark].X[0];
		z = actx->markers[imark].X[2];

		if(x*x + (z-dist)*(z-dist)<=rad*rad)
		{
			// 2D circle
			actx->markers[imark].phase = 1;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkInitFilePolygons"
PetscErrorCode ADVMarkInitFilePolygons(AdvCtx *actx, UserCtx *user)
{
	// REDUNDANTLY loads a file with 2D-polygons that coincide with the marker planes
	// each processor uses the full polygonal shapes to find assign phase ids to local markers

	FDSTAG       *fs;
	int           fd;	
	PetscViewer   view_in;
	char         *LoadFileName;
	PetscScalar   header[2];
	PetscInt      tstart[3],tend[3], nmark[3], nidx[3], nidxmax;
	PetscInt      k,n,kvol,Fcount,Fsize,VolN,Nmax,Lmax,kpoly;
	Polygon2D     Poly;
	PetscBool     AddRandomNoise, ReducedOutput_Polygons;
	PetscInt     *polyin;
	PetscInt     *idx;
	PetscScalar  *X,*PolyLen,*PolyIdx,*PolyFile;
	PetscInt      imark,imarkx,imarky,imarkz,icellx,icelly,icellz;
	PetscScalar   dx=0.0,dy=0.0,dz=0.0,x=0.0,y=0.0,z=0.0;
	PetscScalar   chLen;
	PetscLogDouble t0,t1;
	char          normalDir[4] = {"xyz"};
	PetscRandom   rctx;
	PetscScalar   cf_rand;
	PetscInt      nPoly;
	PetscScalar   atol;
	PetscScalar   box[4];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get marker context
	fs = actx->fs;

	// random noise
	AddRandomNoise = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL, NULL,"-AddRandomNoiseParticles" , &AddRandomNoise , NULL); CHKERRQ(ierr);
	if(AddRandomNoise) PetscPrintf(PETSC_COMM_WORLD, " Adding random noise to marker distribution \n");


	ReducedOutput_Polygons = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL, NULL,"-ReducedOutput_Polygons" , &ReducedOutput_Polygons , NULL); CHKERRQ(ierr);
	if(ReducedOutput_Polygons) PetscPrintf(PETSC_COMM_WORLD, " Reduced output for Polygons activated \n");


	// initialize the random number generator
	ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);
	
	// set characteristic length
	chLen = actx->jr->scal.length;

	// --- initialize markers ---
	
	// marker counter
	imark  = 0;
	icellx = 0;
	icelly = 0;
	icellz = 0;

	// initialize makers in a processor wise manner
	for(imarkz = 0; imarkz < fs->dsz.ncels*user->NumPartZ; imarkz++)
	{
		if (!(imarkz%user->NumPartZ))
		{
			dz = (fs->dsz.ncoor[icellz+1] - fs->dsz.ncoor[icellz]) / (PetscScalar) (user->NumPartZ);
			z  = fs->dsz.ncoor[icellz] + 0.5*dz;
			icellz++;
		}
		else
		{
			z += dz;
		}
		icelly = 0;

		for(imarky = 0; imarky < fs->dsy.ncels*user->NumPartY; imarky++)
		{
			if (!(imarky%user->NumPartY))
			{
				dy = (fs->dsy.ncoor[icelly+1] - fs->dsy.ncoor[icelly]) / (PetscScalar) (user->NumPartY);
				y  = fs->dsy.ncoor[icelly] + 0.5*dy;
				icelly++;
			}
			else
			{
				y += dy;
			}
			icellx = 0;

			for(imarkx = 0; imarkx < fs->dsx.ncels*user->NumPartX; imarkx++)
			{
				if (!(imarkx%user->NumPartX))
				{
					dx = (fs->dsx.ncoor[icellx+1] - fs->dsx.ncoor[icellx]) / (PetscScalar) (user->NumPartX);
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
				
				
				if(AddRandomNoise)
				{
					// add random noise
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[0] += (cf_rand-0.5)*dx/( (PetscScalar) user->NumPartX);
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[1] += (cf_rand-0.5)*dy/( (PetscScalar) user->NumPartY);
					ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
					actx->markers[imark].X[2] += (cf_rand-0.5)*dz/( (PetscScalar) user->NumPartZ);
				}

				// increment local counter
				imark++;
			}
		}
	}

	// --- local grid/marker info ---

	// get first global index of marker plane
	tstart[0] = fs->dsx.pstart * user->NumPartX;
	tstart[1] = fs->dsy.pstart * user->NumPartY;
	tstart[2] = fs->dsz.pstart * user->NumPartZ;

	// get local number of markers per direction
	nmark[0]  = fs->dsx.ncels * user->NumPartX;
	nmark[1]  = fs->dsy.ncels * user->NumPartY;
	nmark[2]  = fs->dsz.ncels * user->NumPartZ;

	// get last global index of marker plane
	for (k=0;k<3;k++)
	{
		tend[k]   = tstart[k] + nmark[k] - 1;
	}

	// how many markers on the marker plane ?
	nidx[0] = nmark[1] * nmark[2]; nidxmax = nidx[0];
	nidx[1] = nmark[0] * nmark[2]; if (nidx[1] > nidxmax) nidxmax = nidx[1];
	nidx[2] = nmark[0] * nmark[1]; if (nidx[2] > nidxmax) nidxmax = nidx[2];

	// compile input file name
	asprintf(&LoadFileName, "./%s/%s",
		user->LoadInitialParticlesDirectory,
		user->ParticleFilename);

	PetscPrintf(PETSC_COMM_WORLD," Loading polygons redundantly from file: %s \n", LoadFileName);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, LoadFileName, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header & size of file
	ierr = PetscBinaryRead(fd, &header, 2, PETSC_SCALAR); CHKERRQ(ierr);
	Fsize = (PetscInt)(header[1]);

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)Fsize  *sizeof(PetscScalar),&PolyFile); CHKERRQ(ierr);
	Fcount=0;

	// read entire file 
	ierr = PetscBinaryRead(fd, PolyFile, Fsize, PETSC_SCALAR); CHKERRQ(ierr);

	// read number of volumes
	VolN = (PetscInt)(PolyFile[Fcount]); Fcount++;
	Nmax = (PetscInt)(PolyFile[Fcount]); Fcount++;
	Lmax = (PetscInt)(PolyFile[Fcount]); Fcount++;

    // allocate space for index array & the coordinates of the largest polygon
	ierr = PetscMalloc((size_t)Nmax  *sizeof(PetscScalar),&PolyLen); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)Nmax  *sizeof(PetscScalar),&PolyIdx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)Lmax*2*sizeof(PetscScalar),&Poly.X); CHKERRQ(ierr);

	// allocate temporary arrays
	ierr = PetscMalloc((size_t)nidxmax*sizeof(PetscInt),&idx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nidxmax*sizeof(PetscBool),&polyin); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)nidxmax*2*sizeof(PetscScalar),&X); CHKERRQ(ierr);


	// --- loop over all volumes ---
	for (kvol=0; kvol<VolN; kvol++)
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
		for (kpoly=0; kpoly<Poly.num;kpoly++)
		{
			PolyIdx[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// get lengths of polygons (PetscScalar !)
		for (kpoly=0; kpoly<Poly.num;kpoly++)
		{
			PolyLen[kpoly] = PolyFile[Fcount]; Fcount++;
		}

		// --- loop through all slices ---
		for (kpoly=0; kpoly<Poly.num;kpoly++)
		{
			// read polygon
			Poly.len  = (PetscInt)(PolyLen[kpoly]);
			Poly.gidx = (PetscInt)(PolyIdx[kpoly]);
			Poly.lidx = (PetscInt)(PolyIdx[kpoly])-tstart[Poly.dir];

			// check if slice is part of local proc
			if (Poly.gidx  >= tstart[Poly.dir] && Poly.gidx <= tend[Poly.dir])
			{
				// read polygon
				for (n=0; n<Poly.len*2;n++)
				{
					Poly.X[n] = PolyFile[Fcount]; Fcount++;
				}

				// get local markers that locate on polygon plane
				ADVMarkSecIdx(actx,user,Poly.dir,Poly.lidx, idx);
				for (k=0;k<nidx[Poly.dir];k++)
				{
					X[k*2]   = actx->markers[idx[k]].X[Poly.ax[0]] * chLen;
					X[k*2+1] = actx->markers[idx[k]].X[Poly.ax[1]] * chLen;
				}

				// get bounding box of a polygon
				nPoly = Poly.len;

				polygon_box(&nPoly, Poly.X, 1e-12, &atol, box);

				in_polygon(nidx[Poly.dir], X, nPoly, Poly.X, box, atol, polyin);

				// set marker phase
				for (k=0;k<nidx[Poly.dir];k++)
				{
					if (polyin[k])
					{
						if (Poly.type == 1) // additive
						{
							actx->markers[idx[k]].phase += Poly.phase;
						}
						else if (Poly.type == 2) // grid additive
						{
							if ( actx->markers[idx[k]].phase % 2 == 1 ) // avoid adding twice when contours are over imposed (e.g. at grid intersection)
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

		if (!ReducedOutput_Polygons){
			PetscPrintf(PETSC_COMM_WORLD,"[Rank 0] Created vol %lld/%lld [%g sec]: phase %lld, type %lld, %lld slices, %c-normal-dir; found %lld markers \n",(LLD)kvol+1,(LLD)VolN, t1-t0, (LLD)Poly.phase, (LLD)Poly.type,(LLD)Poly.num, normalDir[Poly.dir], (LLD)Poly.nmark);
		}
	}

	// Set temperature from file if a Temperature file is specified in the input
	if(strcmp(user->TemperatureFilename,"noTemperatureFileName")!=0)
	{
		ierr = ADVMarkSetTempFromFile(actx,user);
		CHKERRQ(ierr);
	}
	
	// free
	PetscFree(idx);
	PetscFree(polyin);
	PetscFree(X);

	PetscFree(PolyIdx);
	PetscFree(PolyLen);
	PetscFree(Poly.X);

	PetscFree(PolyFile);
	
	// destroy random context
	ierr = PetscRandomDestroy(&rctx);    CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
	free(LoadFileName);

	// wait until all processors finished reading markers
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished setting markers with polygons\n");
	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSetTempFromFile"
PetscErrorCode ADVMarkSetTempFromFile(AdvCtx *actx, UserCtx *user)
{
	FDSTAG       *fs;
	int          fd;
	Marker       *P;
	PetscViewer  view_in;
	char         *LoadFileName;
	PetscScalar  header[2],dim[3];
	PetscInt     Fsize,  imark,nummark, nmarkx, nmarky, nmarkz;
	PetscScalar  DX,DY,DZ;
	PetscScalar  xp,yp,zp, Xc, Yc, Zc, xpL, ypL, zpL;
	PetscScalar  *Temp;
	//PetscScalar  *dim;
	PetscInt Ix,Iy,Iz;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	PetscScalar chTemp;
	
	chTemp = actx->jr->scal.temperature;
	fs = actx->fs;

	// create filename
	asprintf(&LoadFileName, "./%s/%s",
	user->LoadInitialParticlesDirectory,
	user->TemperatureFilename);

	PetscPrintf(PETSC_COMM_WORLD," Loading temperature redundantly from file: %s \n", LoadFileName);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, LoadFileName, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

	// read (and ignore) the silent undocumented file header & size of file
	ierr = PetscBinaryRead(fd, &header, 2, PETSC_SCALAR); CHKERRQ(ierr);
	Fsize = (PetscInt)(header[1])-3;

	// allocate space for entire file & initialize counter
	ierr = PetscMalloc((size_t)Fsize*sizeof(PetscScalar), &Temp); CHKERRQ(ierr);

	// read entire file
	ierr = PetscBinaryRead(fd, &dim, 3,     PETSC_SCALAR); CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd, Temp, Fsize, PETSC_SCALAR); CHKERRQ(ierr);

	// grid spacing
	DX = user->W/(dim[0] - 1.0);
	DY = user->L/(dim[1] - 1.0);
	DZ = user->H/(dim[2] - 1.0);

	// get local number of markers
	nmarkx  = fs->dsx.ncels * user->NumPartX;
	nmarky  = fs->dsy.ncels * user->NumPartY;
	nmarkz  = fs->dsz.ncels * user->NumPartZ;
	nummark = nmarkx*nmarky*nmarkz;
	
	PetscInt nx, ny;
	nx = (PetscInt)dim[0];
	ny = (PetscInt)dim[1];
	//nummark = 2;

	for(imark = 0; imark < nummark; imark++)
	{
		// get curent marker
		P = &actx->markers[imark];

		// get global marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// index of the lower left corner of the element (of the temperature grid) in which the particle is
		Ix = (PetscInt)floor((xp - user->x_left) /DX);
		Iy = (PetscInt)floor((yp - user->y_front)/DY);
		Iz = (PetscInt)floor((zp - user->z_bot)  /DZ);

		//T3D_IxIyIz = Temp[Iz*dim[0]*dim[1] + Iy*dim[0] + Ix ];

		// Coordinate of the first corner (lower left deepest)
		Xc = user->x_left + (PetscScalar)Ix*DX;
		Yc = user->y_front+ (PetscScalar)Iy*DY;
		Zc = user->z_bot  + (PetscScalar)Iz*DZ;

		// Local coordinate of the particle inside a temperature element
		xpL = (xp - Xc)/DX;
		ypL = (yp - Yc)/DY;
		zpL = (zp - Zc)/DZ;

		// Interpolate value on the particle using trilinear shape functions
		P->T = (
		(1.0-xpL) * (1.0-ypL) * (1.0-zpL) * Temp[Iz    *nx*ny + Iy     * nx + Ix   ] +
		  xpL   * (1.0-ypL) * (1.0-zpL) * Temp[Iz    *nx*ny + Iy     * nx + Ix+1 ] +
		  xpL   *   ypL   * (1.0-zpL) * Temp[Iz    *nx*ny + (Iy+1) * nx + Ix+1 ] +
		(1.0-xpL) *   ypL   * (1.0-zpL) * Temp[Iz    *nx*ny + (Iy+1) * nx + Ix   ] +
		(1.0-xpL) * (1.0-ypL) *   zpL   * Temp[(Iz+1)*nx*ny + Iy     * nx + Ix   ] +
		  xpL   * (1.0-ypL) *   zpL   * Temp[(Iz+1)*nx*ny + Iy     * nx + Ix+1 ] +
		  xpL   *   ypL   *   zpL   * Temp[(Iz+1)*nx*ny + (Iy+1) * nx + Ix+1 ] +
		(1.0-xpL) *   ypL   *   zpL   * Temp[(Iz+1)*nx*ny + (Iy+1) * nx + Ix   ] )/chTemp;
	}

	// Clear memory
	PetscFree(Temp);

	ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
	free(LoadFileName);

	PetscFunctionReturn(ierr);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ADVMarkSecIdx"
void ADVMarkSecIdx(AdvCtx *actx, UserCtx *user, PetscInt dir, PetscInt Islice, PetscInt *idx)
{
	FDSTAG   *fs;
	PetscInt i,ix,iy,iz,nmarkx,nmarky,nmarkz;
	PetscInt d,c;
	
	// get fdstag info
	fs = actx->fs;

	// get local number of markers
	nmarkx  = fs->dsx.ncels * user->NumPartX;
	nmarky  = fs->dsy.ncels * user->NumPartY;
	nmarkz  = fs->dsz.ncels * user->NumPartZ;

	if (dir == 0) // yz plane
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
	else if (dir == 1) // xz plane
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
	else if (dir == 2) // xy plane
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
/*
#undef __FUNCT__
#define __FUNCT__ "inpoly"
void inpoly(PolyCtx *polydat, PetscInt N, PetscScalar *X, PetscScalar *node, PetscInt Nnode, PetscBool *in, PetscBool *bnd)
{
    PetscInt *cnect;
    PetscScalar *Xtemp , *nodetemp;
    PetscBool *cn , *on;
    PetscScalar *x, *y;
    PetscInt *idx;
    PetscScalar tol, tol1 , eps = 2.2204e-016, temp , minix=10e20, maxix=-10e20, miniy=10e20, maxiy=-10e20, norminf=-10e20;
    PetscScalar X1, X2 , Y1, Y2, xmin, xmax , ymin, ymax, XX, YY , x2x1 , y2y1;
    PetscScalar lim , ub;
    PetscInt i , j, l , i2 , ind , k ;
    PetscInt n1 , n2 , lower , upper, start;
    PetscInt N1, Ncnect;

    // WARNING! UNUSED PARAMETER!
    if(polydat) polydat = NULL;

    // Retrieve allocated arrays
	Xtemp = polydat->Xtemp;
	nodetemp = polydat->nodetemp;
	cnect = polydat->cnect;
	cn = polydat->cn;
	on = polydat->on;
	x = polydat->x;
	y = polydat->y;
	idx = polydat->idx;

	// constants
	N1 = N-1;
	Ncnect = Nnode;

	// allocate and initialize temporary arrays

	PetscMalloc((size_t)N*2*sizeof(PetscScalar),&Xtemp);
	for(i = 0; i < 2*N; i++)
	{
		Xtemp[i] = X[i];
	}

	PetscMalloc((size_t)Nnode*2*sizeof(PetscScalar),&nodetemp);
	for(i = 0; i < 2*Nnode; i++)
	{
		nodetemp[i] = node[i];
	}

	PetscMalloc((size_t)Ncnect*2*sizeof(PetscInt),&cnect);
	for (i = 0; i < Ncnect - 1; i++)
	{
		i2            = 2*i;
		cnect[i2]     = i + 1;
		cnect[i2 + 1] = i + 2;
	}
	cnect[Nnode*2-2] = Nnode;
	cnect[Nnode*2-1] = 1.0;

    // Temporary vectors
	PetscMalloc((size_t)N*sizeof(PetscBool),&cn);
	PetscMalloc((size_t)N*sizeof(PetscBool),&on);
	PetscMalloc((size_t)N*sizeof(PetscScalar),&x);
	PetscMalloc((size_t)N*sizeof(PetscScalar),&y);
	PetscMalloc((size_t)N*sizeof(PetscInt),&idx);

	// unmodified code of Darren Engwirda & Sebastien Paris ---
    for (i = 0; i < 2*N ; i=i+2)
    {
        temp  = Xtemp[i];
        if(temp < minix)
        {
            minix = temp;
        }
        if(temp > maxix)
        {
            maxix = temp;
        }
        temp  = X[i + 1];
        if(temp < miniy)
        {
            miniy = temp;
        }
        if(temp > maxiy)
        {
            maxiy = temp;
        }
    }
    for (i = 0 ; i < 2*Nnode ; i = i + 2)
    {
        temp = fabs(node[i]);
        if(temp > norminf)
        {
            norminf = temp;
        }
        temp = fabs(node[i + 1]);
        if(temp > norminf)
        {
            norminf = temp;
        }
    }
    tol  = norminf*eps;
    lim  = norminf + tol;
	tol1 = tol + 1.0;

    if ((maxix - minix) > (maxiy - miniy))
    {
        for (i = 0; i < 2*N ; i=i+2)
        {
            temp             = Xtemp[i];
            Xtemp[i]         = Xtemp[i + 1];
            Xtemp[i + 1]     = temp;
        }
        for (i = 0 ; i < 2*Nnode ; i=i+2)
        {
            temp             = nodetemp[i];
            nodetemp[i]      = nodetemp[i + 1];
            nodetemp[i + 1]  = temp;
        }
    }
    for (i = 0 ; i< N ; i++)
    {
        cn[i]    = 0;
        on[i]    = 0;
        y[i]     = Xtemp[2*i + 1];
        idx[i]   = i;

        in[i]    = 0; // added this to initialize *in and *bnd
        bnd[i]   = 0;
    }

    qsindex( y , idx , 0 , N1 );
    for (i = 0 ; i < N ; i++)
    {
        x[i]     = Xtemp[2*idx[i]];
    }

    for (k = 0 ; k < Ncnect ; k++)
    {
        n1   = 2*(cnect[2*k]     - 1);
        n2   = 2*(cnect[2*k + 1] - 1);
        X1   = nodetemp[n1];
        Y1   = nodetemp[n1 + 1];
        X2   = nodetemp[n2];
        Y2   = nodetemp[n2 + 1];
		x2x1 = X2 - X1;
		y2y1 = Y2 - Y1;

        if (X1 > X2)
        {
            xmin = X2;
            xmax = X1;
        }
        else
        {
            xmin = X1;
            xmax = X2;
        }
        if (Y1 > Y2)
        {
            ymin = Y2;
            ymax = Y1;
        }
        else
        {
            ymin = Y1;
            ymax = Y2;
        }
        if (y[0] == ymin)
        {
            start = 0;
        }
        else if (y[N1] <= ymin)
        {
            start = N1;
        }
        else
        {
            lower = 0;
            upper = N1;
            start = ((lower+upper)/2);
            for (l = 0 ; l < N ; l++)
            {
                if (y[start] < ymin)
                {
                    lower = start;
                    start = ((lower+upper)/2);
                }
                else if (y[start-1]<ymin)
                {
                    break;
                }
                else
                {
                    upper = start;
                    start = ((lower+upper)/2);
                }
            }
            start--;
        }

		start = max(0 , start);
        for (j = start ; j < N ; j++)
        {
            YY = y[j];
            if (YY <= ymax)
            {
                if (YY >= ymin)
                {
                    XX = x[j];

                    if (XX >= xmin)
                    {
                        if (XX <= xmax)
                        {
                            on[j] = on[j] || ( fabs( (Y2 - YY)*(X1 - XX) - (Y1 - YY)*(X2 - XX) ) < tol );
                            if (YY < ymax)
                            {
                                ub = ( x2x1*(Y1 - YY) - y2y1*(X1 - XX) )/( (XX - lim)*y2y1 );
                                if ( (ub > -tol) && (ub < tol1 ) )
                                {
                                    cn[j] = !cn[j];
                                }
                            }
                        }
                    }
                    else if (YY < ymax)
                    {
                        cn[j] = !cn[j];
                    }
                    
                    // Code added by Arthur
                    // solves a bug, but may trigger another in the future (?)
                    else if (YY==ymax)
                    {
                        cn[j] = 1;
                    }
                    
                    // ---
                    
                }
            }
            else
            {
                break;
            }
        }
    }
    for(i = 0; i < N ; i++)
    {
        ind      = idx[i];
        in[ind]  = (cn[i] || on[i]);
        bnd[ind] = on[i];
    }

     PetscFree(cn);
     PetscFree(on);
     PetscFree(x);
     PetscFree(y);
     PetscFree(idx);
     PetscFree(Xtemp);
     PetscFree(nodetemp);
     PetscFree(cnect);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CreatePolyCtx"
PetscErrorCode CreatePolyCtx(PolyCtx *polydat, PetscInt N, PetscInt Nnode)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate and initialize temporary arrays
	ierr = PetscMalloc((size_t)N*2*sizeof(PetscScalar),&polydat->Xtemp);
	ierr = PetscMalloc((size_t)Nnode*2*sizeof(PetscScalar),&polydat->nodetemp);
	ierr = PetscMalloc((size_t)Nnode*2*sizeof(PetscScalar),&polydat->cnect);
	// Temporary vectors
	ierr = PetscMalloc((size_t)N*sizeof(PetscBool),&polydat->cn); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)N*sizeof(PetscBool),&polydat->on); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)N*sizeof(PetscScalar),&polydat->x); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)N*sizeof(PetscScalar),&polydat->y); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)N*sizeof(PetscInt),&polydat->idx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DestroyPolyCtx"
PetscErrorCode DestroyPolyCtx(PolyCtx polydat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscFree(polydat.cn); CHKERRQ(ierr);
	ierr = PetscFree(polydat.on); CHKERRQ(ierr);
	ierr = PetscFree(polydat.x); CHKERRQ(ierr);
	ierr = PetscFree(polydat.y); CHKERRQ(ierr);
	ierr = PetscFree(polydat.idx); CHKERRQ(ierr);
	ierr = PetscFree(polydat.Xtemp); CHKERRQ(ierr);
	ierr = PetscFree(polydat.nodetemp); CHKERRQ(ierr);
	ierr = PetscFree(polydat.cnect); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void qsindex (PetscScalar  *a, PetscInt *idx , PetscInt lo, PetscInt hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted

    PetscInt i=lo, j=hi , ind;
    PetscScalar x=a[(lo+hi)/2] , h;

    // partition
    do
    {
        while (a[i]<x) i++;
        while (a[j]>x) j--;
        if (i<=j)
        {
            h        = a[i];
			a[i]     = a[j];
			a[j]     = h;
			ind      = idx[i];
			idx[i] = idx[j];
			idx[j] = ind;
            i++;
			j--;
        }
    }
	while (i<=j);

    // recursion
    if (lo<j) qsindex(a , idx , lo , j);
    if (i<hi) qsindex(a , idx , i , hi);
}
//---------------------------------------------------------------------------
*/
