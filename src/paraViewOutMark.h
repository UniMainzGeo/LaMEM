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
#ifndef __paraViewOutMark_h__
#define __paraViewOutMark_h__
//---------------------------------------------------------------------------
//................ ParaView marker output driver object .....................
//---------------------------------------------------------------------------

struct FB;
struct AdvCtx;

//---------------------------------------------------------------------------

struct PVMark
{
	AdvCtx    *actx;              // advection context
	char      outfile[_STR_LEN_]; // output file name
	long int  offset;             // pvd file offset
	PetscInt  outmark;            // marker output flag
	PetscInt  outpvd;             // pvd file output flag

};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVMarkCreate(PVMark *pvmark, FB *fb);

// write all time-step output files to disk (PVD, PVTU, VTU)
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime, PetscInt tindx);

// .vtu marker output
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName);

// .pvtu marker output
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName);

//---------------------------------------------------------------------------

#endif
