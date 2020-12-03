/***    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
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
 **    filename:   paraViewOutMark.h
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
//..............   MARKER PARAVIEW XML OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __paraViewOutPassiveTracers_h__
#define __paraViewOutPassiveTracers_h__
//---------------------------------------------------------------------------
//................ ParaView marker output driver object .....................
//---------------------------------------------------------------------------

struct FB;
struct AdvCtx;

//---------------------------------------------------------------------------

struct PVPtr
{
	AdvCtx    *actx;              // advection context
	char      outfile[_str_len_+20]; // output file name
	long int  offset;             // pvd file offset
	PetscInt  outptr;             // marker output flag
	PetscInt  outpvd;             // pvd file output flag
	PetscInt  Temperature;
	PetscInt  Pressure;
	PetscInt  Phase;
	PetscInt  MeltFraction;
	PetscInt  ID;
	PetscInt  Active;
	PetscInt  Grid_mf;

};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVPtrCreate(PVPtr *pvptr, FB *fb);

// write all time-step output files to disk (PVD, PVTU, VTU)
PetscErrorCode PVPtrWriteTimeStep(PVPtr *pvptr, const char *dirName, PetscScalar ttime);

// .vtu marker output
PetscErrorCode PVPtrWriteVTU(PVPtr *pvptr, const char *dirName);

// .pvtu marker output
PetscErrorCode PVPtrWritePVTU(PVPtr *pvptr, const char *dirName);


#endif
