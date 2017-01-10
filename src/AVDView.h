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
 **    filename:   AVDView.h
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

#ifndef __AVDView_h__
#define __AVDView_h__

//---------------------------------------------------------------------------

typedef struct
{
	AdvCtx    *actx;    // advection context
	char      *outfile; // output file name
	long int  offset;   // pvd file offset
	PetscInt  outavd;   // AVD output flag
	PetscInt  refine;   // Voronoi Diagram refinement factor
	PetscInt  outpvd;   // pvd file output flag

} PVAVD;

//---------------------------------------------------------------------------

PetscErrorCode PVAVDCreate(PVAVD *pvavd, AdvCtx *actx, const char *filename);

PetscErrorCode PVAVDDestroy(PVAVD *pvavd);

PetscErrorCode PVAVDReadFromOptions(PVAVD *pvavd);

PetscErrorCode PVAVDWriteTimeStep(PVAVD *pvavd, const char *dirName, PetscScalar ttime, PetscInt tindx);

PetscErrorCode PVAVDWritePVTR(PVAVD *pvavd, const char *dirName);

PetscErrorCode PVAVDWriteVTR(PVAVD *pvavd, const char *dirName);

//---------------------------------------------------------------------------
#endif
