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



struct Melt_Extraction_t
{
	FreeSurf  *surf;
	Vec   gdMV, gdMVmerge,gdMoho1, gdMoho, Thickness;
	Vec   ldMV,ldMoho,Miphase,Miphase2;
	Vec   gdc, ldc;
	Vec   SurfMoho,LSurfMoho;
	Vec   ldvecmerge, dgmvvecmerge;
	Vec   DMin, DMax,LDMin,LDMax;
	Vec   MeltID0,MeltID1,MeltID2,MeltID3;
	Vec   MeltID0S,MeltID1S,MeltID2S,MeltID3S;
};





PetscErrorCode DBMatReadMeltExtraction_Par(DBMat *dbm, FB *fb);
PetscErrorCode MeltExtractionCreate(JacRes *jr);
PetscErrorCode MeltExtractionDestroy(JacRes *jr);
PetscErrorCode MeltExtractionSave(JacRes *jr, AdvCtx *actx);
PetscErrorCode MeltExtractionUpdate(JacRes *jr, AdvCtx *actx);
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx, PetscInt ID_ME);
PetscErrorCode MeltExtractionInterpMarkerBackToGrid(AdvCtx *actx);
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr,PetscInt ID_ME,PetscInt update, AdvCtx *actx);
PetscErrorCode Moho_Tracking(JacRes *jr);
PetscErrorCode Extrusion_melt(FreeSurf *surf,PetscInt ID_ME,AdvCtx *actx);
PetscErrorCode Compute_Thickness(JacRes *jr);
PetscErrorCode MeltExtractionPhaseRatio(AdvCtx *actx);
PetscErrorCode ReadMelt_Extraction(JacRes *jr, FILE *fp);
PetscErrorCode Melt_Extraction__WriteRestart(JacRes *jr, FILE *fp);
PetscScalar Compute_dM(PetscScalar mf, Melt_Ex_t *M_Ex_t, PetscScalar dt);
PetscScalar Compute_dMex_Marker(AdvCtx *actx,PetscInt ID,PetscInt iphase );
PetscErrorCode Compute_Comulative_Melt_Extracted(JacRes *jr, AdvCtx *actx,PetscInt ID_ME,  Melt_Ex_t *M_Ex_t, PetscInt update);
PetscErrorCode Set_to_zero_Volumetric_source(JacRes *jr);
PetscErrorCode Update_Volumetric_source(JacRes *jr);
PetscErrorCode Save_Pressure_Temperature_Melt_Extracted(JacRes *jr,PetscInt ID_ME);
PetscErrorCode Set_to_zero_Vector(JacRes *jr);








#endif
