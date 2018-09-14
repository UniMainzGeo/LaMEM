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
 **    filename:   paraViewOutSurf.h
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
//..............   FREE SRUFACE PARAVIEW XML OUTPUT ROUTINES   ..............
//---------------------------------------------------------------------------
#ifndef __paraViewOutSurf_h__
#define __paraViewOutSurf_h__
//---------------------------------------------------------------------------
// maximum number of components in the output vector
#define _max_num_comp_surf_ 3

//---------------------------------------------------------------------------

struct FB;
struct FreeSurf;

//---------------------------------------------------------------------------
//................ ParaView free surface output driver object ...............
//---------------------------------------------------------------------------

struct PVSurf
{
	FreeSurf  *surf;               // free surface object
	char       outfile[_STR_LEN_]; // output file name
	float     *buff;               // direct output buffer
	long int   offset;             // pvd file offset
	PetscInt   outsurf;            // free surface output flag
	PetscInt   outpvd;             // pvd file output flag
	PetscInt   velocity;           // velocity output flag
	PetscInt   topography;         // surface topography output flag
	PetscInt   amplitude;          // topography amplitude output flag
	PetscInt   newcontinental;    // new continental crust produced
	PetscInt   newmafic;          // new mafic crust produced


};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVSurfCreate(PVSurf *pvsurf, FB *fb);

// create buffer array
PetscErrorCode PVSurfCreateData(PVSurf *pvsurf);

// destroy ParaView output driver
PetscErrorCode PVSurfDestroy(PVSurf *pvsurf);

// write all time-step output files to disk (PVD, PVTS, VTS)
PetscErrorCode PVSurfWriteTimeStep(PVSurf *pvsurf, const char *dirName, PetscScalar ttime);

// parallel output file .pvts
PetscErrorCode PVSurfWritePVTS(PVSurf *pvsurf, const char *dirName);

// sequential output file .vts
PetscErrorCode PVSurfWriteVTS(PVSurf *pvsurf, const char *dirName);




//---------------------------------------------------------------------------

void OutputBufferWrite(
	FILE     *fp,
	float    *buff,
	PetscInt  cn);

PetscErrorCode PVSurfWriteCoord(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteVel(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteTopo(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteAmplitude(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteNewContinental(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteNewMafic(PVSurf *pvsurf, FILE *fp);



//---------------------------------------------------------------------------
#endif
