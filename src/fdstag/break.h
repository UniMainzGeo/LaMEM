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
 **    filename:   break.h
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
//........................  BREAKPOINT ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __break_h__
#define __break_h__

// check if breakpoints exist and restart new simulation if not
PetscErrorCode BreakCheck    (UserCtx *user);

//---------------------------------------------------------------------------
// Write breakpoint files
//---------------------------------------------------------------------------
// grid, mark, gsol, gtopo, info are written together
PetscErrorCode BreakWrite(UserCtx *user,
						AdvCtx    *actx,
						FreeSurf  *surf,
						PVOut     *pvout,
						PVSurf    *pvsurf,
						PVMark    *pvmark,
						PVAVD     *pvavd,
						JacType    jtype);
//---------------------------------------------------------------------------
// Read breakpoint files
//---------------------------------------------------------------------------
// gsol, info are read together; grid, topo, mark are read separately
PetscErrorCode BreakRead(UserCtx *user,
						AdvCtx   *actx,
						PVOut    *pvout,
						PVSurf   *pvsurf,
						PVMark   *pvmark,
						PVAVD    *pvavd,
						JacType  *jtype);

PetscErrorCode BreakReadGrid (UserCtx *user, FDSTAG *fs);
PetscErrorCode BreakReadSurf (FDSTAG *fs, FreeSurf *surf);
PetscErrorCode BreakReadMark (AdvCtx *actx);

//---------------------------------------------------------------------------
// Read and write global vectors
//---------------------------------------------------------------------------
PetscErrorCode BreakWriteVec (FILE *fp, Vec x, PetscInt n);
PetscErrorCode BreakReadVec  (FILE *fp, Vec x, PetscInt n);

//---------------------------------------------------------------------------
// Read and write FDSTAG 1D structures
//---------------------------------------------------------------------------
void BreakWriteDiscret1D(FILE *fp, Discret1D *ds, MeshSeg1D *ms);
void BreakReadDiscret1D (FILE *fp, Discret1D *ds, MeshSeg1D *ms);
//---------------------------------------------------------------------------

#endif
