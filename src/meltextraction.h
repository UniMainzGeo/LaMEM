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
 **    filename:   adjoint.h
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
//Melt Extraction Routine
//---------------------------------------------------------------------------
#ifndef __meltextraction_h__
#define __meltextraction__
//---------------------------------------------------------------------------

PetscErrorCode MeltExtractionCreate(JacRes *jr);
PetscErrorCode MeltExtractionDestroy(JacRes *jr);
PetscErrorCode MeltExtractionSave(JacRes *jr, AdvCtx *actx);
PetscErrorCode MeltExtractionUpdate(JacRes *jr, AdvCtx *actx);
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx, PetscInt iphase);
PetscErrorCode MeltExtractionInterpMarkerBackToGrid(AdvCtx *actx);
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr,PetscInt iphase,PetscInt update, AdvCtx *actx);
PetscErrorCode MeltExtractionInject(JacRes *jr,AdvCtx *actx, PetscInt ID, PetscInt I, PetscInt J, PetscInt K, PetscScalar UP,PetscInt iphase, PetscInt newM, PetscInt n_in);
PetscErrorCode Moho_Tracking(FreeSurf *surf);
PetscErrorCode Extrusion_melt(FreeSurf *surf,PetscInt iphase,AdvCtx *actx);
PetscErrorCode ExchangeMassME(SolVarBulk *svBulk,PetscScalar dx,PetscScalar dy,PetscScalar dz, PetscScalar dt);
PetscErrorCode Compute_Thickness(JacRes *jr);
PetscErrorCode OutPutVolume(FreeSurf *surf,PetscInt iphase);
PetscErrorCode AverageP_T(JacRes *jr, PetscInt iphase);
#endif
