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
 **    filename:   adjoint.c
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
 **    This routine:
 **         Andrea Piccolo
 **         Georg  Reuber
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/


#include "LaMEM.h"
#include "cvi.h"
#include "phase.h"
#include "fdstag.h"
#include "advect.h"
#include "tools.h"
#include "JacRes.h"
#include "scaling.h"
#include "parsing.h"
#include "surf.h"
#include "tssolve.h"
#include "bc.h"
#include "subgrid.h"
#include "AVD.h"
#include "parsing.h"
#include "meltextraction.h"
#include "constEq.h"



//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBMatReadMeltExtraction_Par"
PetscErrorCode DBMatReadMeltExtraction_Par(DBMat *dbm, FB *fb)
{
	// read softening law from file

	Melt_Ex_t   *melt_par;
	Scaling      *scal;
	char         Type_[_str_len_];
	PetscInt  ID;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm->scal;

	// softening law ID
	ierr = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbm->numMEPar-1); CHKERRQ(ierr);

	// get pointer to specified softening law
	melt_par = dbm->matMexT + ID;
	// check ID
	if(melt_par->ID != -1)
	{
		 SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate melt extraction transition law!");
	}

	// set ID
	melt_par->ID = ID;

	ierr = getStringParam(fb, _REQUIRED_, "Name", melt_par->Name,0);  CHKERRQ(ierr);
	if (!melt_par->Name)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the name of melt extraction law ", (LLD)ID);
	}

	ierr    =   getStringParam(fb, _REQUIRED_, "Type",Type_,NULL);  CHKERRQ(ierr);

	if(!strcmp(Type_,"_Constant_"))
	{
		melt_par->Type = _Constant_mf_;
	}
	else if(!strcmp(Type_,"_Constant_flux_"))
	{
		melt_par->Type = _Constant_flux_;

	}
	else
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct type of melt extraction ", (LLD)ID);

	}

	ierr = getScalarParam(fb, _REQUIRED_, "M_left",&melt_par->Mleft, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "M_trs", &melt_par->Mtrs, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "M_Max", &melt_par->Mmax, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "IR",    &melt_par->IR, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "T_Int", &melt_par->TInt, 1, 1.0); CHKERRQ(ierr);

	melt_par->TInt=(melt_par->TInt+scal->Tshift)/scal->temperature;


	ierr = getIntParam(fb, _REQUIRED_, "Phase_effusive", &melt_par->PhExt, 1, _max_num_phases_); CHKERRQ(ierr);
	ierr = getIntParam(fb, _REQUIRED_, "Phase_intrusive", &melt_par->PhInt, 1, _max_num_phases_); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "VolCor",&melt_par->VolCor, 1, 1.0); CHKERRQ(ierr);
	if(!melt_par->VolCor)
	{
		melt_par->VolCor = 0.85; // Alternatively one might a more sophisticated approach
	}


	ierr = getScalarParam(fb, _REQUIRED_, "DInt",&melt_par->DInt, 1, 1.0); CHKERRQ(ierr);
	if(melt_par->Type == _Constant_flux_)
	{
		melt_par->mf_cap_flux = -1; // default value
		ierr = getScalarParam(fb, _REQUIRED_, "timescale",&melt_par->timescale, 1, 1.0); CHKERRQ(ierr);
		ierr = getIntParam(fb, _OPTIONAL_, "mf_max_flux",&melt_par->mf_cap_flux, 1, 1); CHKERRQ(ierr);


		melt_par->timescale = melt_par->timescale/ scal->time;

	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD,"   Melt Extraction Law [%lld]                                      \n", (LLD)(melt_par->ID))                                             ;
	PetscPrintf(PETSC_COMM_WORLD,"   Name                                            :   %s          \n", melt_par->Name)                                                  ;
	PetscPrintf(PETSC_COMM_WORLD,"   Type                                            :   %s          \n", Type_)                                                           ;
	PetscPrintf(PETSC_COMM_WORLD,"   Mleft                                           :   %1.3f [n.d.]\n", melt_par->Mleft)                                                 ;
	PetscPrintf(PETSC_COMM_WORLD,"   Mtrs  (retained melt fraction)                  :   %1.3f [n.d.]\n", melt_par->Mtrs)                                                  ;
	PetscPrintf(PETSC_COMM_WORLD,"   MMax  (change phase melt extraction threshold)  :   %1.3f [n.d.]\n", melt_par->Mmax)                                                  ;
	PetscPrintf(PETSC_COMM_WORLD,"   DInt (normalized crustal depth intrusion        :   %1.3f [n.d.]\n", melt_par->DInt)                                                  ;
	if(melt_par->Type == _Constant_flux_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"   timescale                                       :   %1.3f [%s]  \n", melt_par->timescale*scal->time,scal->lbl_time)    ;
		PetscPrintf(PETSC_COMM_WORLD,"   volumetric_flux (Mleft/timescale)               :   %1.3f [1/%s]\n", melt_par->Mleft/(melt_par->timescale)*scal->time,scal->lbl_time) ;
	}
	PetscPrintf(PETSC_COMM_WORLD,"   T Intrusion                                     :   %1.0f [%s]\n",   melt_par->TInt*scal->temperature - scal->Tshift, scal->lbl_temperature)                           ;
	PetscPrintf(PETSC_COMM_WORLD,"   IR (proportion of intrusion)                    :   %1.3f [n.d.]\n", melt_par->IR)                           ;
	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");




	PetscFunctionReturn(0);
}




//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionCreate"
PetscErrorCode MeltExtractionCreate(JacRes *jr)
{ // First functions that is called
	Melt_Extraction_t *MEPar;

	PetscErrorCode ierr;
	PetscFunctionBegin;
	if(jr->ctrl.MeltExt==0) PetscFunctionReturn(0);

	MEPar = jr->MEPar;


	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &MEPar->gdc)                       ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &MEPar->gdMoho)                    ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &MEPar->Miphase)                   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &MEPar->Miphase2)                   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &MEPar->gdMoho1)                   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &MEPar->dgmvvecmerge)              ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &MEPar->Thickness)                 ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->ldc)                        ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &MEPar->ldvecmerge)                 ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &MEPar->ldMoho)                     ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &MEPar->DMin)                    ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &MEPar->DMax)                    ;     CHKERRQ(ierr);

	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &MEPar->LDMin)                    ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &MEPar->LDMax)                    ;     CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->SurfMoho)                   ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF, &MEPar->LSurfMoho)                   ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF, &jr->surf->lmagmathick)                     ;     CHKERRQ(ierr);




	if(jr->dbm->numMEPar==1)
	{
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID0S)                 ;     CHKERRQ(ierr);

		ierr = VecZeroEntries(MEPar->MeltID0S); CHKERRQ(ierr);


	}
	else if(jr->dbm->numMEPar==2)
	{
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID1)                 ;     CHKERRQ(ierr);

		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID1S)                 ;     CHKERRQ(ierr);

		ierr = VecZeroEntries(MEPar->MeltID0S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID1S); CHKERRQ(ierr);

	}
	else if(jr->dbm->numMEPar==3)
	{
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID1)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID2)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID1S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID2S)                 ;     CHKERRQ(ierr);

		ierr = VecZeroEntries(MEPar->MeltID0S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID1S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID2S); CHKERRQ(ierr);


	}
	else
	{
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID1)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID2)                 ;     CHKERRQ(ierr);
		ierr = DMCreateLocalVector(jr->fs->DA_CEN, &MEPar->MeltID3)                 ;     CHKERRQ(ierr);

		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID1S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID2S)                 ;     CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &MEPar->MeltID3S)                 ;     CHKERRQ(ierr);


		ierr = VecZeroEntries(MEPar->MeltID0S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID1S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID2S); CHKERRQ(ierr);
		ierr = VecZeroEntries(MEPar->MeltID3S); CHKERRQ(ierr);

	}






	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionDestroy"
PetscErrorCode MeltExtractionDestroy(JacRes *jr)
{// Last
	Melt_Extraction_t *MEPar;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if(jr->ctrl.MeltExt==0) PetscFunctionReturn(0);
	MEPar = jr->MEPar;

	ierr = VecDestroy(&MEPar->ldMoho)                           ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->gdMoho)                           ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->gdMoho1)                          ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->Miphase)                          ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->Miphase2)                         ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->dgmvvecmerge)                     ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->ldvecmerge)                       ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->Thickness)                        ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->ldc)                              ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->gdc)                              ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->DMin)                             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->DMax)                             ;     CHKERRQ(ierr);


	ierr = VecDestroy(&MEPar->LDMin)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->LDMax)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->SurfMoho)                   ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->LSurfMoho)                   ;     CHKERRQ(ierr);
	ierr = VecDestroy(&MEPar->lmagmathick)                   ;     CHKERRQ(ierr);






	if(jr->dbm->numMEPar==1)
		{
			ierr = VecDestroy(&MEPar->MeltID0)                 ;     CHKERRQ(ierr);


			ierr = VecDestroy(&MEPar->MeltID0S)                 ;     CHKERRQ(ierr);



		}
		else if(jr->dbm->numMEPar==2)
		{
			ierr = VecDestroy( &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID1)                 ;     CHKERRQ(ierr);

			ierr = VecDestroy(&MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID1S)                 ;     CHKERRQ(ierr);



		}
		else if(jr->dbm->numMEPar==3)
		{
			ierr = VecDestroy( &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID1)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID2)                 ;     CHKERRQ(ierr);

			ierr = VecDestroy(&MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID1S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID2S)                 ;     CHKERRQ(ierr);



		}
		else
		{
			ierr = VecDestroy( &MEPar->MeltID0)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID1)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID2)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy( &MEPar->MeltID3)                 ;     CHKERRQ(ierr);

			ierr = VecDestroy(&MEPar->MeltID0S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID1S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID2S)                 ;     CHKERRQ(ierr);
			ierr = VecDestroy(&MEPar->MeltID3S)                 ;     CHKERRQ(ierr);


		}

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionSave"
PetscErrorCode MeltExtractionSave(JacRes *jr,AdvCtx *actx)
{    /* Short explanation: Melt extraction works using different classes. Each classes work as if it were a softening law.
 	 	 LaMEM loop over the melt extraction parametrization and each phase that belongs to this classes co-partecipate to the total volume.
 	 	 After this stage the information is send to melt extraction exchange volume. In this function all the contribute are sum up, generating
 	 	 a redundant matrix containing the total volume extracted from each column. Then the volumetric deformation associated with the melt extraction
 	 	 is distribuited within a Dz which is equal (or proportional) to the effective thickness. Since the melt extraction is extremely complicated
 	 	 function and require in average 0.5 seconds, it is not possible to couple with non linear iteration. As a consequence this function is useful
 	 	 only for generating the volumetric source term, which is enforced in the residuum of the continuity equation.
	*/
	Melt_Ex_t      *M_Ex_t	;
	PetscInt       update;
	PetscInt       ID_ME;


	PetscErrorCode ierr;
	PetscFunctionBegin;


	if((jr->ctrl.initGuess == 1) | (jr->ctrl.MeltExt == 0)) PetscFunctionReturn(0);

	update = 0;

	ierr=Moho_Tracking(jr); CHKERRQ(ierr);

	ierr=Compute_Thickness(jr); CHKERRQ(ierr);

	ierr= Set_to_zero_Volumetric_source(jr); CHKERRQ(ierr);

	for(ID_ME=0;ID_ME<jr->dbm->numMEPar;ID_ME++)
	{

		// Access to the melt extraction context
		M_Ex_t = jr->dbm->matMexT+ID_ME;

		ierr = Set_to_zero_Vector(jr); CHKERRQ(ierr);

		ierr = Compute_Comulative_Melt_Extracted(jr,actx, ID_ME,M_Ex_t); CHKERRQ(ierr);

		ierr = MeltExtractionExchangeVolume(jr,ID_ME,update,actx);            CHKERRQ(ierr);

		ierr = Update_Volumetric_source(jr);                                  CHKERRQ(ierr);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "MeltExtractionUpdate"
PetscErrorCode MeltExtractionUpdate(JacRes *jr, AdvCtx *actx)
{
	/*This function must be called inside the free surface routine. Such that it could handle even the "melt" sedimentation by phase.
	 */
	PetscErrorCode ierr;
	PetscFunctionBegin;
	// Variable declaration
	Melt_Ex_t    		*M_Ex_t	;



	PetscInt       update;
	PetscInt       ID_ME;

	if((jr->ctrl.initGuess == 1) | (jr->ctrl.MeltExt == 0)) PetscFunctionReturn(0);


	update = 1;

	// Initialize & get array

	ierr = ADVInterpMarkToCell(actx);CHKERRQ(ierr);

	for(ID_ME=0;ID_ME<jr->dbm->numMEPar;ID_ME++)
		{

			// Access to the melt extraction context
			M_Ex_t = jr->dbm->matMexT+ID_ME;

			ierr = Set_to_zero_Vector(jr); CHKERRQ(ierr);

			ierr = Compute_Comulative_Melt_Extracted(jr,actx, ID_ME,M_Ex_t); CHKERRQ(ierr);

			ierr = Save_Pressure_Temperature_Melt_Extracted(jr,ID_ME);            CHKERRQ(ierr);

			ierr = MeltExtractionExchangeVolume(jr,ID_ME,update,actx);            CHKERRQ(ierr);

			ierr = MeltExtractionInterpMarker(actx, ID_ME) ;                      CHKERRQ(ierr);

			ierr = MeltExtractionPhaseRatio(actx);                                CHKERRQ(ierr);

			ierr = MeltExtractionInterpMarkerBackToGrid(actx);                   CHKERRQ(ierr);

		}

	PetscFunctionReturn(0);

}
//-------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionExchangeVolume"
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr, PetscInt ID_ME,PetscInt update,AdvCtx *actx)
{

	Melt_Ex_t    		*M_Ex_t	;
	Melt_Extraction_t 	*Mext	;
	FDSTAG         		*fs   	;
    Discret1D 			*dsz	;
    FreeSurf            *surf   ;


	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, L,K1,K2,in;
	Vec          global_volume ;
	PetscScalar  bz, ez;
	PetscScalar  IR, dx, dy, dz,Z,DZ;
	PetscScalar  ***Thickness, D, D1, ***MohoG,***DMin,***DMax,vol,Vol_Cor;
	PetscScalar  *vdgmvvec, *vdgmvvecmerge, ***vdgmvvecmerge2, ***vdgmvvec2, ***Mipbuff;

	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	M_Ex_t = jr->dbm->matMexT;
	Mext   = jr->MEPar;
	fs = jr->fs;
	dsz = &fs->dsz;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	IR = M_Ex_t[ID_ME].IR; // Amount of intrusion that has to be injected within the crust
	surf = jr->surf;
	Vol_Cor = M_Ex_t[ID_ME].VolCor;


	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
	// Initialize the necessary vector
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &global_volume); CHKERRQ(ierr);
	ierr = VecZeroEntries(global_volume) ; CHKERRQ(ierr);


	// Retrieve the vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, Mext->Miphase, &Mipbuff); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, global_volume, &vdgmvvec2); CHKERRQ(ierr);

	// scan all local cells
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		vdgmvvec2[L][j][i] += Vol_Cor*Mipbuff[k][j][i];
	}
	END_STD_LOOP

	// Restore the vector
	ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->Miphase , &Mipbuff)     ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, global_volume, &vdgmvvec2)     ; CHKERRQ(ierr);
	// Perform the integral
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(global_volume, &vdgmvvec); CHKERRQ(ierr);
		ierr = VecGetArray(Mext->dgmvvecmerge, &vdgmvvecmerge); CHKERRQ(ierr);
		ierr = MPI_Allreduce(vdgmvvec, vdgmvvecmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(global_volume, &vdgmvvec); CHKERRQ(ierr);
		ierr = VecRestoreArray(Mext->dgmvvecmerge, &vdgmvvecmerge); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(global_volume,Mext->dgmvvecmerge);  CHKERRQ(ierr);
	}



	// Access to the vectors


	ierr = DMDAVecGetArray(fs->DA_CEN, Mext->Miphase, &Mipbuff)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->dgmvvecmerge, &vdgmvvecmerge2)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,Mext->Thickness,&Thickness); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,Mext->ldMoho,&MohoG); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->LDMax, &DMax)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->LDMin, &DMin)   ; CHKERRQ(ierr);

	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);
	// Scan all the cells
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	START_PLANE_LOOP
	{
		dx = SIZE_CELL(i,sx,fs->dsx);
		dy = SIZE_CELL(j,sy,fs->dsy);
		if(IR>0.0)
		{

			DZ=0.0;
			D1 =Thickness[L][j][i];
			D = MohoG[L][j][i] + D1*M_Ex_t[ID_ME].DInt;
			DMin[L][j][i]=D+0.5*IR*vdgmvvecmerge2[L][j][i]/(dx*dy);
			DMax[L][j][i]=D-0.5*IR*vdgmvvecmerge2[L][j][i]/(dx*dy);
			vol = 0.0;

			if(vdgmvvecmerge2[L][j][i] < 0.0)
			{
				// The grid is assumed to BE uniform.

					vol=-IR*(vdgmvvecmerge2[L][j][i])/(DMax[L][j][i]-DMin[L][j][i]);

					// Are the extreme of the interval in the same processor?
						if(DMin[L][j][i] >= bz && DMin[L][j][i] < ez && DMax[L][j][i] >= bz && DMax[L][j][i] < ez)
						{

							ierr = Discret1DFindPoint(&fs->dsz, DMin[L][j][i], K1); CHKERRQ(ierr);
							ierr = Discret1DFindPoint(&fs->dsz, DMax[L][j][i], K2); CHKERRQ(ierr);

							// Are the two point the same point?
							if(K1==K2)
							{

								Mipbuff[sz+K1][j][i] += -IR*(vdgmvvecmerge2[L][j][i]);
							}
							else
							{
								for(in=sz+K1;in<=sz+K2;in++)
									{
									dz = SIZE_CELL(in,sz,fs->dsz);
									if(in==sz+K1)
									{
										Z=COORD_CELL(in, sz, fs->dsz);
										dz=(dz/2)+(Z-DMin[L][j][i]);
									}
									if(in==sz+K2)
									{
										Z=COORD_CELL(in, sz, fs->dsz);
										dz=(dz/2)-(Z-DMax[L][j][i]);
									}

									Mipbuff[in][j][i]+=vol*dz;
									DZ+=dz;
									}
							}

						}
						else if(DMin[L][j][i] >= bz && DMin[L][j][i] < ez )
						{

							ierr = Discret1DFindPoint(&fs->dsz, DMin[L][j][i], K1); CHKERRQ(ierr);
							for(in=sz+K1;in<sz+nz;in++)
							{
								dz = SIZE_CELL(in,sz,fs->dsz);
								if(in==sz+K1)
								{
									Z=COORD_CELL(in, sz, fs->dsz);
									dz=(dz/2)+(Z-DMin[L][j][i]);
								}
								Mipbuff[in][j][i]+=vol*dz;
							}
						}
						else if(DMax[L][j][i] >= bz && DMax[L][j][i] < ez)
						{

							ierr = Discret1DFindPoint(&fs->dsz, DMax[L][j][i], K2); CHKERRQ(ierr);
							for(in=sz;in<=sz+K2;in++)
							{
								dz = SIZE_CELL(in,sz,fs->dsz);
								if(in==sz+K2)
								{
									Z=COORD_CELL(in, sz, fs->dsz);
									dz=(dz/2)-(Z-DMax[L][j][i]);
								}
								Mipbuff[in][j][i]+=vol*dz;

							}
						}
						else if (bz >=DMin[L][j][i] && bz < DMax[L][j][i] &&ez >=DMin[L][j][i] && ez < DMax[L][j][i] )
						{

							for(in=sz;in<sz+nz;in++)
							{
							dz = SIZE_CELL(in,sz,fs->dsz);
							Mipbuff[in][j][i]+=vol*dz;
							}
				}
			}
			if(update==1)
				{
					// The first time that Exchange volume is called, has as input the mass. The second time
					// it takes into account the volume. In order to retrive the thickness of the "melt extracted"
					// you need to divide for the area.
					dx = SIZE_CELL(i,sx,fs->dsx);
					dy = SIZE_CELL(j,sy,fs->dsy);

					vdgmvvecmerge2[L][j][i]= -(vdgmvvecmerge2[L][j][i])/(dx*dy);

				}
		}
		else
		{
			if(update==1)
			{
				vdgmvvecmerge2[L][j][i]= -(vdgmvvecmerge2[L][j][i]/(dx*dy));

			}
		}


	}
	END_PLANE_LOOP
	// Restore Vectors
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->dgmvvecmerge, &vdgmvvecmerge2); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->Miphase , &Mipbuff) ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,Mext->ldMoho,&MohoG) ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,Mext->Thickness,&Thickness) ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->LDMax, &DMax)   ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->LDMin, &DMin)   ; CHKERRQ(ierr);

	// Restore global & not useful vector
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&global_volume); CHKERRQ(ierr);

	// Create the local vector


	if(update==1)
	{
		GLOBAL_TO_LOCAL(jr->DA_CELL_2D,Mext->dgmvvecmerge,Mext->ldvecmerge)

		ierr=Extrusion_melt(surf,ID_ME,actx); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------------//

#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarker"
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx, PetscInt ID_ME)
{
	//=======================================================================
	// interpolate increments of the history field of melt extraction to markers
	//=======================================================================

	Melt_Ex_t    		*M_Ex_t	;
	Melt_Extraction_t 	*Mext	;
	PData               *pd;
	FDSTAG *fs;
	JacRes *jr;
	Marker *P;
	Material_t *phases;
	PetscScalar  ***Dm_save,***DMin,***DMax;
	PetscInt nx, ny, sx, sy, sz;
	PetscInt jj, ID, I, J, K, newphase,Ph_int;
	PetscInt L;
	PetscScalar  DM,newME, T_Int;
	PetscScalar  mfeff,dM, ***p,***T;


	PetscErrorCode ierr;
	PetscFunctionBegin;
	fs = actx->fs;
	jr = actx->jr;
	pd = jr->Pd;
	Mext = jr->MEPar;
	M_Ex_t = jr->dbm->matMexT+ID_ME;
	phases = jr->dbm->phases;
	L = (PetscInt)fs->dsz.rank; // rank of the processor

	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	Ph_int  = M_Ex_t->PhInt;
	T_Int   = M_Ex_t->TInt;


	// access 3D layouts of local vectors


	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->LDMax, &DMax)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->LDMin, &DMin)   ; CHKERRQ(ierr);


	for(jj = 0; jj < actx->nummark; jj++)
		{
			// access next marker
			P = &actx->markers[jj];
			// get consecutive index of the host cell
			ID = actx->cellnum[jj];
			// expand I, J, K cell indices
			GET_CELL_IJK(ID, I, J, K, nx, ny)


			if(P->X[2]>=DMin[L][sy+J][sx+I] && P->X[2]<=DMax[L][sy+J][sx+I])
			{

				P->phase = Ph_int;
				P->T     = T_Int;
				P->APS   = 0.0;
				P->U[0]  = 0.0;
				P->U[1]  = 0.0;
				P->U[2]  = 0.0;
			}


		}


	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->LDMax, &DMax)   ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->LDMin, &DMin)   ; CHKERRQ(ierr);

	ierr =  MeltExtractionPhaseRatio(actx); CHKERRQ(ierr);



	if(ID_ME==0)
	{
		ierr = DMDAVecGetArray(fs->DA_CEN, Mext->MeltID0, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME == 1)
	{
		ierr = DMDAVecGetArray(fs->DA_CEN, Mext->MeltID1, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==2)
	{
		ierr = DMDAVecGetArray(fs->DA_CEN, Mext->MeltID2, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==3)
	{
		ierr = DMDAVecGetArray(fs->DA_CEN, Mext->MeltID3, &Dm_save); CHKERRQ(ierr);
	}
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp, &p); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &T); CHKERRQ(ierr);



	for(jj = 0; jj < actx->nummark; jj++)
	{
		dM = 0.0;

		DM = 0.0;
		// access next marker
		P = &actx->markers[jj];
		// get consecutive index of the host cell
		ID = actx->cellnum[jj];
		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		DM = Dm_save[sz+K][sy+J][sx+I];

		if(DM!=0.0)
		{

			if(phases[P->phase].ID_MELTEXT==ID_ME)
			{
				ierr =  setDataPhaseDiagram(pd, p[sz+K][sy+J][sx+I], T[sz+K][sy+J][sx+I], jr->dbm->phases[P->phase].pdn); CHKERRQ(ierr);

				mfeff = pd->mf - P->MExt;
				dM = Compute_dM(mfeff, M_Ex_t, jr->ts->dt);

				P->MExt +=dM;
				P->MTot +=dM;

				if(P->MExt >= M_Ex_t->Mmax)
				{
					newphase = jr->dbm->phases[P->phase].PhNext;
					newME   = P->MExt -M_Ex_t->Mmax;
					P->phase=newphase;
					P->MExt=newME;

				}
			}

			}
	}



	if(ID_ME==0)
	{
		ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->MeltID0, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME == 1)
	{
		ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->MeltID1, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==2)
	{
		ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->MeltID2, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==3)
	{
		ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->MeltID3, &Dm_save); CHKERRQ(ierr);
	}

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp, &p); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT, &T); CHKERRQ(ierr);


	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);



	PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Moho_Tracking"
PetscErrorCode Moho_Tracking(JacRes *jr)
{
	FDSTAG *fs;
	Discret1D *dsz;
	Melt_Extraction_t 	*Mext	;


	PetscInt i, j, k,numPhases,ii,cnt,gcnt,I1,I2,J1,J2,mx,my;
	PetscInt sx, sy, sz, nx, ny, nz, iter,L;
	PetscScalar bz, ez, bottom;
	PetscScalar ***Mohovec2,*Mohovec22,*MMerge ,MantP;
	Material_t *phases; // Phases
	PetscScalar *phRat;
	PetscScalar ***Moho_s,Moho_b[4],***Moho,buf;

	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	fs = jr->fs;
	Mext = jr->MEPar;
	dsz = &fs->dsz;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	phases = jr->dbm->phases;
	numPhases = jr->dbm->numPhases; // take the number of phases from dbm structures

	// Find the bottom of the domain
	ierr = FDSTAGGetGlobalBox(fs, 0, 0, &bottom, 0, 0, 0); CHKERRQ(ierr);
	// get local coordinate bounds
	ierr = VecSet   (Mext->gdMoho1,bottom);	 CHKERRQ(ierr);
	ierr = VecZeroEntries   (Mext->gdMoho);	 CHKERRQ(ierr);
	ierr = DMDAVecGetArray  (jr->DA_CELL_2D,Mext->gdMoho1, &Mohovec2);	CHKERRQ(ierr);
	// Std Loop
	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	iter = 0;
	START_STD_LOOP{
		phRat = jr->svCell[iter++].phRat;
		MantP=0.0;
		for(ii=0;ii<numPhases;ii++)
		{
			if(phRat[ii] && phases[ii].pMant==1) MantP+=phRat[ii];
		}
		if(MantP>0.0)
		{
			Mohovec2[L][j][i] = COORD_NODE(k, sz, fs->dsz);
		}
	}END_STD_LOOP

	ierr =DMDAVecRestoreArray(jr->DA_CELL_2D,Mext->gdMoho1,&Mohovec2); CHKERRQ(ierr);

	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(Mext->gdMoho1, &Mohovec22)           ; CHKERRQ(ierr);
		ierr = VecGetArray(Mext->gdMoho, &MMerge) ; CHKERRQ(ierr);
		ierr = MPI_Allreduce(Mohovec22, MMerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_MAX, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(Mext->gdMoho1, &Mohovec22); CHKERRQ(ierr);
		ierr = VecRestoreArray(Mext->gdMoho, &MMerge); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(Mext->gdMoho1,Mext->gdMoho);  CHKERRQ(ierr);
	}
	// Find the values at the top
	GLOBAL_TO_LOCAL(jr->DA_CELL_2D, Mext->gdMoho, Mext->ldMoho);


	ierr = VecZeroEntries(Mext->SurfMoho); CHKERRQ(ierr);
	ierr = VecZeroEntries(Mext->LSurfMoho); CHKERRQ(ierr);

	// retrieve average topography
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;


	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->ldMoho,  &Moho);  CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->surf->DA_SURF, Mext->SurfMoho,  &Moho_s);  CHKERRQ(ierr);
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	cnt=0;
		// scan all free surface local points
		START_PLANE_LOOP
		{
			// This function is highly based on the sedimentation routine in surf.cpp. There is no need to
			// use dt, since the thickness is already in its integral form by default.

			I1 = i;
			I2 = i-1;
			J1 = j;
			J2 = j-1;
			// check index bounds if ghost points are undefined
			if(I1 == mx) I1--;
			if(I2 == -1) I2++;
			if(J1 == my) J1--;
			if(J2 == -1) J2++;

			Moho_b[0] = Moho[L][J1][I1];
			Moho_b[1] = Moho[L][J1][I2];
			Moho_b[2] = Moho[L][J2][I1];
			Moho_b[3] = Moho[L][J2][I2];

			buf = (Moho_b[0] + Moho_b[1] + Moho_b[2] + Moho_b[3])/4.0;

			// store advected topography

			Moho_s[L][j][i]=buf;

		}
		END_PLANE_LOOP

		// restore access
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->ldMoho,  &Moho);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(jr->surf->DA_SURF, Mext->SurfMoho,  &Moho_s);  CHKERRQ(ierr);


		if(ISParallel(PETSC_COMM_WORLD))
		{
			ierr = MPI_Allreduce(&cnt, &gcnt, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
		}
		else
		{
			gcnt = cnt;
		}

		GLOBAL_TO_LOCAL(jr->surf->DA_SURF, Mext->SurfMoho, Mext->LSurfMoho);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Extrusion_melt"
PetscErrorCode Extrusion_melt(FreeSurf *surf,PetscInt ID_ME, AdvCtx *actx)
{
	JacRes *jr;
	FDSTAG *fs;
	Melt_Ex_t    		*M_Ex_t	;
	Melt_Extraction_t 	*Mext	;


	PetscScalar ***topo,***lmelt,***Dm_save;
	PetscScalar zbot, ztop, z, Melt[4],Vol[4],dx,dy ,Layer,Ext_fraction,Vol_t;
	PetscInt L, cnt, gcnt;
	PetscInt i, j, nx, ny, sx, sy, I1, I2, J1, J2, mx, my;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	jr = surf->jr;
	Mext = jr->MEPar;
	M_Ex_t= jr->dbm->matMexT;

	fs = jr->fs;
	L  = (PetscInt)fs->dsz.rank;
	surf->MeltExtraction = 1;

	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;

	// get z-coordinates of the top and bottom boundaries
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, &ztop); CHKERRQ(ierr);
	Ext_fraction = 0.0;
	// store the phase that is being sedimented
	surf->PhExt = M_Ex_t[ID_ME].PhExt;
	Ext_fraction = 1-M_Ex_t[ID_ME].IR;


	if(ID_ME==0)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, Mext->MeltID0S, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME == 1)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, Mext->MeltID1S, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==2)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, Mext->MeltID2S, &Dm_save); CHKERRQ(ierr);
	}
	else if(ID_ME==3)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, Mext->MeltID3S, &Dm_save); CHKERRQ(ierr);
	}




	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Mext->ldvecmerge,  &lmelt);  CHKERRQ(ierr);
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	cnt=0;


	// scan all free surface local points
	START_PLANE_LOOP
	{
		// This function is highly based on the sedimentation routine in surf.cpp. There is no need to
		// use dt, since the thickness is already in its integral form by default.
		dx = SIZE_CELL(i,sx,fs->dsx);
		dy = SIZE_CELL(j,sy,fs->dsy);

		I1 = i;
		I2 = i-1;
		J1 = j;
		J2 = j-1;
		// check index bounds if ghost points are undefined
		if(I1 == mx) I1--;
		if(I2 == -1) I2++;
		if(J1 == my) J1--;
		if(J2 == -1) J2++;

		Melt[0] = lmelt[L][J1][I1];
		Melt[1] = lmelt[L][J1][I2];
		Melt[2] = lmelt[L][J2][I1];
		Melt[3] = lmelt[L][J2][I2];

		Vol[0] = lmelt[L][J1][I1]*dx*dy;
		Vol[1] = lmelt[L][J1][I2]*dx*dy;
		Vol[2] = lmelt[L][J2][I1]*dx*dy;
		Vol[3] = lmelt[L][J2][I2]*dx*dy;


		Layer = (Melt[0] + Melt[1] + Melt[2] + Melt[3])/4.0;
		Vol_t = (Vol[0] + Vol[1] + Vol[2] + Vol[3])/4.0;
		// get topography
		z = topo[L][j][i];
		//PetscPrintf(PETSC_COMM_WORLD," Topo Z = %6f\n",z*jr->scal->length);

		// uniformly advect
		z += Ext_fraction*Layer;

		// check if internal free surface goes outside the model domain
		if(z > ztop) z = ztop;
		if(z < zbot) z = zbot;


		// store advected topography
		if(Ext_fraction >0.0 && M_Ex_t[ID_ME].IR<1.0)
		{
			topo[L][j][i] = z;
		}

		Dm_save[L][j][i] += Vol_t/(dx*dy);

	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Mext->ldvecmerge,  &lmelt);  CHKERRQ(ierr);

	if(ID_ME==0)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, Mext->MeltID0S, &Dm_save); CHKERRQ(ierr);


	}
	else if(ID_ME == 1)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, Mext->MeltID1S, &Dm_save); CHKERRQ(ierr);


	}
	else if(ID_ME==2)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, Mext->MeltID2S, &Dm_save); CHKERRQ(ierr);


	}
	else if(ID_ME==3)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, Mext->MeltID3S, &Dm_save); CHKERRQ(ierr);

	}




	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&cnt, &gcnt, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gcnt = cnt;
	}


	// compute ghosted version of the advected surface topography
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);



	// compute & store average topography
	ierr = FreeSurfGetAvgTopo(surf); CHKERRQ(ierr);

	// compute host cells for all the markers
	ierr = ADVMarkCrossFreeSurf(actx); CHKERRQ(ierr);

	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

	ierr =  MeltExtractionPhaseRatio(actx);

	// update arrays for marker-cell interaction
	ierr = FreeSurfGetAirPhaseRatio(surf); CHKERRQ(ierr);



	// print info

	surf->MeltExtraction = 0;

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Compute_Thickness"
PetscErrorCode Compute_Thickness(JacRes *jr)
{
	FreeSurf *surf;
	Melt_Extraction_t 	*Mext	;


	FDSTAG *fs;
	PetscScalar cz[4],***h,***ntopo,z,***mg;
	PetscInt i, j, sx, sy, nx, ny, L;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	surf=jr->surf;
	Mext = jr->MEPar;
	fs=jr->fs;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	// Initialize the global vector
	ierr = VecZeroEntries(Mext->Thickness); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &ntopo); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D,Mext->Thickness,&h); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D,Mext->ldMoho,&mg); CHKERRQ(ierr);
	// scan all local cells
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		// get topography at cell corners
		cz[0] = ntopo[L][j  ][i  ];
		cz[1] = ntopo[L][j  ][i+1];
		cz[2] = ntopo[L][j+1][i  ];
		cz[3] = ntopo[L][j+1][i+1];
			// get average cell height
		z=(cz[0] + cz[1] + cz[2] + cz[3])/4.0;
		h[L][j][i]=z-mg[L][j][i];
	}
	END_PLANE_LOOP

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &ntopo); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,Mext->Thickness,&h); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,Mext->ldMoho,&mg); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//-----------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarkerBackToGrid"
PetscErrorCode MeltExtractionInterpMarkerBackToGrid(AdvCtx *actx)
{
	FDSTAG *fs;
	JacRes *jr;
	Marker *P;
	SolVarCell *svCell;
	Controls *ctrl;
	PetscInt ID, I, J, K, II, JJ, KK;
	PetscInt ii, jj, numPhases;
	PetscInt nx, ny, sx, sy, sz, nCells;
	PetscScalar xp, yp, zp, xc, yc, zc, wxn, wyn, wzn, wxc, wyc, wzc, w = 0.0;
	PetscScalar UPXY, UPXZ, UPYZ;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;

	PetscErrorCode ierr;

	PetscFunctionBegin;
	fs = actx->fs;
	jr = actx->jr;
	numPhases = actx->dbm->numPhases;
	ctrl= &jr->ctrl;
	if(ctrl->initGuess | !ctrl->MeltExt) PetscFunctionReturn(0);
	// check marker phases
	ierr = ADVCheckMarkPhases(actx); CHKERRQ(ierr);

	//======
	// CELLS
	//======

	// marker-to-grid projection (cell nodes)

	// number of cells
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nCells = fs->nCells;

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
		svCell->svBulk.mfext_cur = 0.0;
		svCell->svBulk.Tot_MExt  = 0.0;
	}

	// scan ALL markers
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
		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables
		svCell->svBulk.mfext_cur += w*P->MExt;
		svCell->svBulk.Tot_MExt  += w*P->MTot;

	}
	ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);


	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		// normalize history variables
		svCell->svBulk.mfext_cur /= w;
		svCell->svBulk.Tot_MExt   /= w;
	}

	//======
	// EDGES
	//======

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// clear local vectors
	ierr = VecZeroEntries(jr->ldxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldyz); CHKERRQ(ierr);

	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// set interpolated fields to defaults
	UPXY = 1.0; UPXZ = 1.0; UPYZ = 1.0;

	// scan ALL markers
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
		xc = fs->dsx.ccoor[I];
		yc = fs->dsy.ccoor[J];
		zc = fs->dsz.ccoor[K];

		// map marker on the control volumes of edge nodes
		if(xp > xc) { II = I+1; } else { II = I; }
		if(yp > yc) { JJ = J+1; } else { JJ = J; }
		if(zp > zc) { KK = K+1; } else { KK = K; }

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get interpolation weights in node control volumes
		wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);
		wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);
		wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);

		UPXY = P->MTot;  UPXZ = P->MTot;  UPYZ = P->MTot;

		// update required fields from marker to edge nodes
		lxy[sz+K ][sy+JJ][sx+II] += wxn*wyn*wzc*UPXY;
		lxz[sz+KK][sy+J ][sx+II] += wxn*wyc*wzn*UPXZ;
		lyz[sz+KK][sy+JJ][sx+I ] += wxc*wyn*wzn*UPYZ;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// assemble global vectors
	LOCAL_TO_GLOBAL(fs->DA_XY, jr->ldxy, jr->gdxy)
	LOCAL_TO_GLOBAL(fs->DA_XZ, jr->ldxz, jr->gdxz)
	LOCAL_TO_GLOBAL(fs->DA_YZ, jr->ldyz, jr->gdyz)

	// access 1D layouts of global vectors
	ierr = VecGetArray(jr->gdxy, &gxy);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdxz, &gxz);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdyz, &gyz);  CHKERRQ(ierr);

	// copy (normalized) data to the residual context
	for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].svDev.mfext_cur = gxy[jj]/jr->svXYEdge[jj].ws;
	for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].svDev.mfext_cur = gxz[jj]/jr->svXZEdge[jj].ws;
	for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].svDev.mfext_cur = gyz[jj]/jr->svYZEdge[jj].ws;

	// restore access
	ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//===========================================================================================
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionPhaseRatio"
PetscErrorCode MeltExtractionPhaseRatio(AdvCtx *actx)
{
	FDSTAG *fs;
	JacRes *jr;
	Marker *P;
	SolVarCell *svCell;
	Controls *ctrl;
	PetscInt ID, I, J, K;
	PetscInt ii, jj, numPhases;
	PetscInt nx, ny, nCells;
	PetscScalar xp, yp, zp, wxc, wyc, wzc, w = 0.0;

	PetscErrorCode ierr;

	PetscFunctionBegin;
	fs = actx->fs;
	jr = actx->jr;
	numPhases = actx->dbm->numPhases;
	ctrl= &jr->ctrl;
	if(ctrl->initGuess | !ctrl->MeltExt) PetscFunctionReturn(0);
	// check marker phases
	ierr = ADVCheckMarkPhases(actx); CHKERRQ(ierr);

	//======
	// CELLS
	//======

	// marker-to-grid projection (cell nodes)

	// number of cells
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nCells = fs->nCells;

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
	}

	// scan ALL markers
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
		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables

	}

	// normalize interpolated values
		for(jj = 0; jj < nCells; jj++)
		{
			// access solution variable
			svCell = &jr->svCell[jj];

			// normalize phase ratios
			ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		}



	PetscFunctionReturn(0);
}

//---------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Melt_Extraction_WriteRestart"
PetscErrorCode Melt_Extraction__WriteRestart(JacRes *jr, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// write solution vectors

	if(jr->dbm->numMEPar==1)
	{
		ierr = VecWriteRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);

	}
	else if(jr->dbm->numMEPar==2)
	{

		ierr = VecWriteRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);


	}
	else if(jr->dbm->numMEPar==3)
	{
		ierr = VecWriteRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID2S, fp); CHKERRQ(ierr);

	}
	else
	{
		ierr = VecWriteRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID2S, fp); CHKERRQ(ierr);
		ierr = VecWriteRestart(jr->MEPar->MeltID3S, fp); CHKERRQ(ierr);

	}
	PetscFunctionReturn(0);
}

//=========================================================

#undef __FUNCT__
#define __FUNCT__ "Read_Melt_Extraction"
PetscErrorCode ReadMelt_Extraction(JacRes *jr, FILE *fp)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read solution vectors

	ierr = MeltExtractionCreate(jr);CHKERRQ(ierr);


	if(jr->dbm->numMEPar==1)
		{
			ierr = VecReadRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);

		}
		else if(jr->dbm->numMEPar==2)
		{

			ierr = VecReadRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);


		}
		else if(jr->dbm->numMEPar==3)
		{
			ierr = VecReadRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID2S, fp); CHKERRQ(ierr);

		}
		else
		{
			ierr = VecReadRestart(jr->MEPar->MeltID0S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID1S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID2S, fp); CHKERRQ(ierr);
			ierr = VecReadRestart(jr->MEPar->MeltID3S, fp); CHKERRQ(ierr);

		}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Compute_dM"
PetscScalar Compute_dM(PetscScalar mfeff, Melt_Ex_t *M_Ex_t, PetscScalar dt)
{
	PetscScalar dM,mf_temp;
	dM = 0.0;



	if(M_Ex_t->Type == _Constant_mf_)
	{
		if(mfeff>M_Ex_t->Mtrs)
		{
			dM = M_Ex_t->Mleft;
			mf_temp=mfeff-dM;
			if(mf_temp<M_Ex_t->Mtrs)
			{
				dM=mfeff-M_Ex_t->Mtrs;
			}
		}
		else
		{
			dM = 0.0;
		}
	}
	if(M_Ex_t->Type == _Constant_flux_)
	{
		PetscScalar dM_flux;
		dM_flux = M_Ex_t->Mleft/M_Ex_t->timescale;
		dM =dM_flux * dt;
		if(M_Ex_t->mf_cap_flux == 1)
		{
			if(dM>M_Ex_t->Mleft)
			{
				dM = M_Ex_t->Mleft;
			}
		}


		if(mfeff>M_Ex_t->Mtrs)
		{
			mf_temp = mfeff-dM;

			if(mf_temp<M_Ex_t->Mtrs)
			{
			dM=mfeff-M_Ex_t->Mtrs;
			}
		}
		else
		{
			dM=0.0;
		}
	}


	return dM;
}
//----------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Compute_dMex_Marker"
PetscScalar Compute_mfeff_Marker(AdvCtx *actx,PetscInt ID,PetscInt iphase, PetscScalar Pr, PetscScalar Tc )
{
	PetscInt         n,ipn,c=0,phase,*mark_id,id_m;
	PetscScalar      mfeff_b,mfeff;
	PData            *pd;

	PetscErrorCode ierr;


	pd  = actx->jr->Pd;
	n = actx->markstart[ID+1] - actx->markstart[ID];
	mark_id = actx->markind + actx->markstart[ID];


	mfeff = 0.0;
	for(ipn=0;ipn<n;ipn++)
	{
		id_m = mark_id[ipn];
		phase = actx->markers[id_m].phase;
		if(phase == iphase)
		{
			ierr = setDataPhaseDiagram(pd, Pr, Tc, actx->dbm->phases[iphase].pdn); CHKERRQ(ierr);

			mfeff += pd->mf-actx->markers[id_m].MExt;

			c ++;
		}
	}

	if(c>0)
	{
		mfeff_b = mfeff/c;
	}
	else
	{
		mfeff_b = 0.0;
	}
	return mfeff_b;
}
//--------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Compute_Comulative_Melt_Extracted"
PetscErrorCode Compute_Comulative_Melt_Extracted(JacRes *jr, AdvCtx *actx,PetscInt ID_ME,  Melt_Ex_t *M_Ex_t)
{
	FDSTAG         		*fs   	;
	PData          		*pd    	;
	DBMat        		*dbm	;
	Melt_Extraction_t 	*Mext	;
	Material_t     *phases;

	PetscScalar    ***Mipbuff;
	PetscInt       i, j, k, nx, ny, nz, sx, sy, sz, iter, iphase;
	PetscInt       numPhases;
	PetscScalar    *phRat;
	PetscScalar    ***p,pc;
	PetscScalar    ***T,Tc;
	PetscScalar    mfeff,dx,dy,dz,dM;
	PetscInt       ID;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	dbm    = jr->dbm;
	fs = jr->fs;
	numPhases = jr->dbm->numPhases;
	pd  = jr->Pd;
	phases =dbm->phases;
	Mext   =  jr->MEPar;

	iter = 0;

	GET_CELL_RANGE(nx, sx, fs->dsx);
	GET_CELL_RANGE(ny, sy, fs->dsy);
	GET_CELL_RANGE(nz, sz, fs->dsz);
	//


	ierr = DMDAVecGetArray(fs->DA_CEN, Mext->Miphase, &Mipbuff); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp, &p); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &T); CHKERRQ(ierr);


	START_STD_LOOP
	{
		dx = SIZE_CELL(i,sx,fs->dsx);
		dy = SIZE_CELL(j,sy,fs->dsy);
		dz = SIZE_CELL(k,sz,fs->dsz);



		phRat = actx->jr->svCell[iter++].phRat; // take phase ratio on the central node
		GET_CELL_ID(ID, i-sx, j-sy, k-sz, fs->dsx.ncels, fs->dsy.ncels)



		for(iphase=0;iphase<numPhases;iphase++)
		{

			if(phases[iphase].ID_MELTEXT == ID_ME)
			{

				if(phRat[iphase]>0.0)
				{

					if(phRat[iphase]>1.0) PetscPrintf(PETSC_COMM_SELF," phRat = %6f \n",phRat[iphase]);

					pc = p[k][j][i];
											// current temperature
					Tc = T[k][j][i];


					dM  	= 0.0;
					mfeff   = 0.0;
					// current pressure
					if(phases[iphase].pdAct==1)
					{
						// check if exist a bit of melt within the cell
						ierr = setDataPhaseDiagram(pd, pc, Tc, phases[iphase].pdn); CHKERRQ(ierr);

						if(pd->mf>0.0)
						{
							// compute the effective melt fraction within the cell, by computing the average mfeff for the all the particles whose phase belongs to the melt extraction law
							mfeff=Compute_mfeff_Marker(actx, ID,iphase,pc,Tc);

						}

						if(mfeff<0.0) mfeff=0.0;

					}
					else
					{
						// place holder for the
					}

					if(mfeff>M_Ex_t->Mtrs)
					{
						// compute the dM
						dM = Compute_dM(mfeff, M_Ex_t, jr->ts->dt);

						Mipbuff[k][j][i] += -phRat[iphase] * dM*dx*dy*dz;

					}
					else
					{
						Mipbuff[k][j][i] +=0.0;
					}

				}
			}
		}
	}END_STD_LOOP


	ierr = DMDAVecRestoreArray(fs->DA_CEN, Mext->Miphase, &Mipbuff); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,&p); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,&T); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(actx->jr->fs->DA_CEN,Mext->Miphase,Mext->ldc)

	PetscFunctionReturn(0);

}
//--------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_to_zero_Volumetric_source"
PetscErrorCode Set_to_zero_Volumetric_source(JacRes *jr)
{
	SolVarBulk     		*svBulk ;
	SolVarCell     		*svCell ;
	PetscInt i,j,k,iter,nx,ny,nz,sx,sy,sz;

	iter = 0;

	GET_CELL_RANGE(nx, sx, jr->fs->dsx);
	GET_CELL_RANGE(ny, sy, jr->fs->dsy);
	GET_CELL_RANGE(nz, sz, jr->fs->dsz);
	START_STD_LOOP{
		svCell = &jr->svCell[iter]; // take the central node based properties
		svBulk = &svCell->svBulk;
		svBulk->Vol_S=0.0;          // set to zero volumetric deformation associate with ME
		iter++;
	}END_STD_LOOP


	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Update_Volumetric_source"
PetscErrorCode Update_Volumetric_source(JacRes *jr)
{
	SolVarBulk            *svBulk ;
	SolVarCell            *svCell ;
	Melt_Extraction_t     *Mext;

	PetscScalar         ***Mipbuff;
	PetscScalar         dx,dy,dz;
	PetscInt i,j,k,iter,nx,ny,nz,sx,sy,sz;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	Mext = jr->MEPar;
	iter=0;
	// Here we compute the right-side term to plug in the residuum (JacRes). Mass_in is the initial mass, while dMass is representing the term that has to be
	// add to svBulk->Mass. This allow us to avoid to dirty JacRes
	ierr = DMDAVecGetArray(jr->fs->DA_CEN, Mext->Miphase, &Mipbuff); CHKERRQ(ierr);

	GET_CELL_RANGE(nx, sx, jr->fs->dsx);
	GET_CELL_RANGE(ny, sy, jr->fs->dsy);
	GET_CELL_RANGE(nz, sz, jr->fs->dsz);

	START_STD_LOOP{
			svCell = &jr->svCell[iter] ;// take the central node based properties
			svBulk = &svCell->svBulk ;
			dx = SIZE_CELL(i,sx,jr->fs->dsx);
			dy = SIZE_CELL(j,sy,jr->fs->dsy);
			dz = SIZE_CELL(k,sz,jr->fs->dsz);

			svBulk->Vol_S += (1-(dx*dy*dz)/(dx*dy*dz+Mipbuff[k][j][i]));


			iter++;
	}END_STD_LOOP

	ierr = DMDAVecRestoreArray(jr->fs->DA_CEN, Mext->Miphase, &Mipbuff) ; CHKERRQ(ierr);


	PetscFunctionReturn(0);

}
//-------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Save_Pressure_Temperature_Melt_Extracted"
PetscErrorCode Save_Pressure_Temperature_Melt_Extracted(JacRes *jr,PetscInt ID_ME)
{
	Melt_Extraction_t     *Mext;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	Mext = jr->MEPar;

	if(ID_ME==0)
		{
			ierr = VecZeroEntries(Mext->MeltID0);     CHKERRQ(ierr);
			ierr = VecCopy(Mext->ldc,Mext->MeltID0);  CHKERRQ(ierr);

		}
		else if(ID_ME == 1)
		{
			ierr = VecZeroEntries(Mext->MeltID1);     CHKERRQ(ierr);
			ierr = VecCopy(Mext->ldc,Mext->MeltID1);  CHKERRQ(ierr);

		}
		else if(ID_ME==2)
		{
			ierr = VecZeroEntries(Mext->MeltID2);     CHKERRQ(ierr);
			ierr = VecCopy(Mext->ldc,Mext->MeltID2);  CHKERRQ(ierr);

		}
		else if(ID_ME==3)
		{
			ierr = VecZeroEntries(Mext->MeltID3);     CHKERRQ(ierr);
			ierr = VecCopy(Mext->ldc,Mext->MeltID3);  CHKERRQ(ierr);

		}


	PetscFunctionReturn(0);

}
// ----------------------------------------------------------------- //
#undef __FUNCT__
#define __FUNCT__ "Set_to_zero_Vector"
PetscErrorCode Set_to_zero_Vector(JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// Set to zero the vector associated with the maximum and minimum depth of intrusion
	ierr = VecZeroEntries(jr->MEPar ->LDMin); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->MEPar->LDMax); CHKERRQ(ierr);


	// Set to zero the vector associated with the total volume extracted/injected in each cell
	ierr = VecZeroEntries(jr->MEPar->Miphase); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->MEPar->ldc); CHKERRQ(ierr);

	ierr = VecZeroEntries(jr->MEPar->dgmvvecmerge); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->MEPar->ldvecmerge); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}




