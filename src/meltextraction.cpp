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
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

/*
 * TODO:
 * - Pay attention local and global might still be inconsistent in Exchange Volume
 *
 * - So far only tested for exactly 1 extraction step
 * - What happens if the melt quantity actually decreases instead of increase? This might be not treated here yet
 * - What pressure and T should the new marker have?
 *   #The temperature is fixed and given as material parameter. The pressure must be interpolated from the surrounding.
 * - Two things are hard-coded (ctrl+F "hard-coded" to find them)
 * - How to handle the actual volume change in injection?:
 * 		1) Delete all markers from cell and inject equally distributed (and volume equal)
 *
 *      2) Reduce volume of existing markers in the cell such that volume is again correct?
 * 		3) Is it phyiscally  correct if it is just corrected for the inflated volumeas is done right now? <?>
 * 		4) Artificially massively increase the injected volume until it is correct? No
 *
 * 	POTENTIAL HELP:
 * 	- Implement a timestep criterium for the melt extraction. Only incrementally move melt up such that velocity doesn't run away.
 */

#include "LaMEM.h"
#include "cvi.h"
#include "phase.h"
#include "fdstag.h"
#include "constEq.h"
#include "advect.h"
#include "tools.h"
#include "JacRes.h"
#include "meltextraction.h"
#include "scaling.h"
#include "parsing.h"
#include "surf.h"
#include "tssolve.h"
#include "bc.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionCreate"
PetscErrorCode MeltExtractionCreate(JacRes *jr)
{ // First functions that is called
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdc)                       ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &jr->gdMoho)                    ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->Miphase)                   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &jr->gdMoho1)                   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &jr->dgmvvecmerge)              ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->DA_CELL_2D, &jr->Thickness)                 ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->fs->DA_CEN, &jr->ldc)                        ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &jr->ldvecmerge)                 ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &jr->ldMoho)                     ;     CHKERRQ(ierr);



	// UnderConstruction && This vectors are propedeutic to visualization of data that I need (Andrea).

	ierr = DMCreateLocalVector(jr->surf->DA_SURF, &jr->surf->lNewMafic)         ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF, &jr->surf->lNewContinental)   ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &jr->surf->NewMafic)         ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF, &jr->surf->NewContinental)   ;     CHKERRQ(ierr);


	ierr = DMCreateGlobalVector(jr->fs->DA_CEN,&jr->PdM)                       ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN,&jr->TdM)                       ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &jr->TC2D)                      ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->DA_CELL_2D, &jr->PC2D)                      ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->TtotMaf)          ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->TtotCon)          ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->PtotMaf)          ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->PtotCon)          ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->VcurrentMaf)      ;     CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->surf->DA_SURF,&jr->surf->VcurrentCon)      ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF,&jr->surf->TProductionCon)    ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF,&jr->surf->TProductionMaf)    ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF,&jr->surf->PProductionCon)    ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF,&jr->surf->PProductionMaf)    ;     CHKERRQ(ierr);
	ierr = DMCreateLocalVector(jr->surf->DA_SURF,&jr->surf->R_Cont)    ;     CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionDestroy"
PetscErrorCode MeltExtractionDestroy(JacRes *jr)
{// Last
	PetscErrorCode ierr;
	PetscFunctionBegin;


	ierr = VecDestroy(&jr->ldMoho)                           ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdMoho)                           ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdMoho1)                          ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->Miphase)                          ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->dgmvvecmerge)                     ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldvecmerge)                       ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->Thickness)                        ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldc)                              ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdc)                              ;     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->surf->NewMafic)                   ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->NewContinental)             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->lNewMafic)                  ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->lNewContinental)            ;     CHKERRQ(ierr);

	ierr = VecDestroy(&jr->PdM)                              ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->TdM)                              ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->TC2D)                             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->PC2D)                             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->TtotMaf)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->TtotCon)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->PtotMaf)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->PtotCon)                    ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->VcurrentMaf)                ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->VcurrentCon)                ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->TProductionCon)             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->TProductionMaf)             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->PProductionCon)             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->PProductionMaf)             ;     CHKERRQ(ierr);
	ierr = VecDestroy(&jr->surf->R_Cont)             ;     CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionSave"
PetscErrorCode MeltExtractionSave(JacRes *jr,AdvCtx *actx)
{    /* 2:Create the structure data that must be used to inject new particles, and compute the associated sink and source term.Copy dMF (see Consteq.cpp)
        into a grid whose coordinate are based on the center of cell, and which save the variable from the bulk variables */

	PetscErrorCode ierr;
	PetscFunctionBegin;


	SolVarBulk     *svBulk    ;
	SolVarCell     *svCell    ;
	FDSTAG         *fs    ;
	PData          *pd    ;
	Controls       *ctrl    ;
	PetscScalar    ***Mipbuff    ;
	PetscInt       i, j, k, nx, ny, nz, sx, sy, sz, iter, iphase    ;
	PetscInt       numPhases    ;
	PetscScalar    *phRat    ;
	PetscScalar    ***p,pc,pShift    ;
	PetscScalar    ***T,Tc    ;
	PetscScalar    mfeff,dx,dy,dz,dM,mf_temp    ;
	Material_t     *mat    ;
	Material_t     *phases;
	PetscInt       update    ;
	fs = jr->fs;
	numPhases = jr->dbm->numPhases;
	pd  = jr->Pd;
	phases = jr->dbm->phases;
	ctrl=&jr->ctrl;
	pShift    =  ctrl->pShift;
	if(ctrl->initGuess) PetscFunctionReturn(0);


	update = 0;

	// Calling Moho_Tracking (working on)
	ierr=Moho_Tracking(jr); CHKERRQ(ierr);
	ierr=Compute_Thickness(jr); CHKERRQ(ierr);

	// Initialize & get array

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp, &p); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &T); CHKERRQ(ierr);

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx);
	GET_CELL_RANGE(ny, sy, fs->dsy);
	GET_CELL_RANGE(nz, sz, fs->dsz);
	START_STD_LOOP{
		svCell = &jr->svCell[iter]; // take the central node based properties
		svBulk = &svCell->svBulk;
		svBulk->dMF=0.0;
		svBulk->Mass=0.0;
		svBulk->mf = 0.0;
		svBulk->dMass = 0.0;
		iter++;
	}END_STD_LOOP

	for(iphase=0;iphase<numPhases;iphase++)
	{
		mat=&phases[iphase];
		if(mat->MeltE>0)
		{
			//Create the buffer
			ierr = VecZeroEntries(jr->Miphase); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_CEN, jr->Miphase, &Mipbuff); CHKERRQ(ierr);
			iter=0;
			START_STD_LOOP
			{
				svCell = &jr->svCell[iter]; // take the central node based properties
				svBulk = &svCell->svBulk ; // take the bulk solution variables
				phRat = jr->svCell[iter++].phRat; // take phase ratio on the central node
				mfeff=0.0;
				mf_temp=0.0;
				if(phRat[iphase])
				{
					// Temperature&Pressure
					// access current pressure
					pc = p[k][j][i]-pShift;
					// current temperature
					Tc = T[k][j][i];
					ierr = SetDataPhaseDiagram(pd, pc, Tc, 0, mat->pdn); CHKERRQ(ierr);
					mfeff = pd->mf-svBulk->mfextot;// historical variables
					if (mfeff<0.0) mfeff=0.0;//Correction
					if( mfeff>=phases[iphase].Mtrs)
					{
						// Compute the Mass Escaping from the local volume
						dx = SIZE_CELL(i,sx,fs->dsx);
						dy = SIZE_CELL(j,sy,fs->dsy);
						dz = SIZE_CELL(k,sz,fs->dsz);
						dM=phases[iphase].Mleft;
						mf_temp=mfeff-dM;
						if(mf_temp<0.0)
						{
							dM=mfeff-0.0;
							mf_temp=0.0;
						}

						Mipbuff[k][j][i] = -phRat[iphase] * dM*dx*dy*dz;//*pd->rho_f;
						svBulk->mf +=phRat[iphase]*mf_temp;
						svBulk->dMF+=phRat[iphase]*dM;
					}
					else
					{
						Mipbuff[k][j][i] = 0.0;
						svBulk->mf +=phRat[iphase]*mfeff;
						svBulk->dMF+=phRat[iphase]*0.0;
					}
				}
			}END_STD_LOOP

			ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase, &Mipbuff); CHKERRQ(ierr);
			// Send the data to Melt Extraction Exchange volume & compute the injection
			// Update Miphase
			ierr = MeltExtractionExchangeVolume(jr,iphase,update,actx); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_CEN,jr->Miphase,&Mipbuff); CHKERRQ(ierr);
			iter=0;
			// Here we compute the right-side term to plug in the residuum (JacRes). Mass_in is the initial mass, while dMass is representing the term that has to be
			// add to svBulk->Mass. This allow us to avoid to dirty JacRes
			START_STD_LOOP{
				svCell = &jr->svCell[iter] ;// take the central node based properties
				svBulk = &svCell->svBulk ;
				svBulk->dMass +=Mipbuff[k][j][i];
				iter++;
			}END_STD_LOOP
			ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase, &Mipbuff) ; CHKERRQ(ierr);\
		}
	}
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,&p); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,&T); CHKERRQ(ierr);
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
	SolVarBulk *svBulk; // pointer 2 the solution variable defined as "bulk cell properties"
	SolVarCell *svCell; // pointer 2 the Solution variable defined in the cell
	FDSTAG *fs; // pointer 2 the grid structure and variables
	PData *pd; // pointer 2 the phase diagram structure
	Controls *ctrl;
	PetscScalar ***Mipbuff; //,***Tcon,***Pcon; // 3D structure storing the melt that has been extracted
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz, iter, iphase; // Iteration
	PetscInt numPhases ; // number of phases
	PetscScalar *phRat; // Phase Ratio
	PetscScalar ***p,***Pcon,pc,pShift; // pressure
	PetscScalar ***T,***Tcon,Tc; // temperature
	PetscScalar mfeff,dx,dy,dz,dM,mf_temp; // Effective melt extracted
	Material_t *mat; // Material properties structure
	Material_t *phases; // Phases
	PetscInt update;

	fs = jr->fs; // take the structured grid data from jr. The out put is a pointer structure
	numPhases = jr->dbm->numPhases; // take the number of phases from dbm structures
	pd  = jr->Pd; // take the structure associated to the phase diagram
	phases = jr->dbm->phases; // take the phases
	ctrl=&jr->ctrl;
	pShift    =  ctrl->pShift;       // pressure shift

	// Initialization
	update = 1;
	fs = jr->fs; // take the structured grid data from jr. The out put is a pointer structure
	numPhases = jr->dbm->numPhases ; // take the number of phases from dbm structures
	pd = jr->Pd ; // take the structure associated to the phase diagram
	phases = jr->dbm->phases; // take the phases
	ctrl=&jr->ctrl;
	pShift    =  ctrl->pShift;       // pressure shift

	// Getting the array
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,&p); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,&T); CHKERRQ(ierr);

	ierr = VecZeroEntries(jr->surf->TProductionMaf);     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->TProductionCon);     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->TtotMaf)       ;     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->TtotCon)       ;     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->PProductionMaf);     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->PProductionCon);     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->PtotMaf)       ;     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->PtotCon)       ;     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->VcurrentMaf)   ;     CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->surf->VcurrentCon)   ;     CHKERRQ(ierr);



	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx);
	GET_CELL_RANGE(ny, sy, fs->dsy);
	GET_CELL_RANGE(nz, sz, fs->dsz);

	for(iphase=0;iphase<numPhases;iphase++)
	{
		mat=&phases[iphase];
		if(mat->MeltE==1)
		{
			//Create the buffer & getting the buffer
			ierr = VecZeroEntries(jr->Miphase); CHKERRQ(ierr);
			ierr = VecZeroEntries(jr->TdM); CHKERRQ(ierr);
			ierr = VecZeroEntries(jr->PdM); CHKERRQ(ierr);


			ierr = DMDAVecGetArray(fs->DA_CEN, jr->Miphase, &Mipbuff); CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_CEN, jr->TdM, &Tcon) ; CHKERRQ(ierr);
			ierr = DMDAVecGetArray(fs->DA_CEN, jr->PdM, &Pcon) ; CHKERRQ(ierr);


			iter=0;
			START_STD_LOOP
			{
				svCell = &jr->svCell[iter] ;// take the central node based properties
				svBulk = &svCell->svBulk ;// take the bulk solution variables
				phRat = jr->svCell[iter++].phRat ; // take phase ratio on the central node
				mfeff=0.0;
				mf_temp=0.0;
				if(phRat[iphase])
				{
					// Temperature&Pressure
					// access current pressure
					pc = p[k][j][i] - pShift;
					//if(iphase==1)PetscPrintf(PETSC_COMM_SELF, "P is %6f\n",pc);

					// current temperature
					Tc = T[k][j][i];
					ierr = SetDataPhaseDiagram(pd, pc, Tc, 0, mat->pdn); CHKERRQ(ierr);
					mfeff = pd->mf-svBulk->mfextot;// historical variables
					if (mfeff<0.0) mfeff=0.0;//Correction
					if( mfeff>=phases[iphase].Mtrs)
					{
						// Compute the Mass Escaping from the local volume
						dx = SIZE_CELL(i,sx,fs->dsx);
						dy = SIZE_CELL(j,sy,fs->dsy);
						dz = SIZE_CELL(k,sz,fs->dsz);
						dM=phases[iphase].Mleft;
						mf_temp=mfeff-dM;
						if(mf_temp<0.0)
						{
							dM=mfeff-0.0;
							mf_temp=0.0;
						}
						Mipbuff[k][j][i] = -phRat[iphase]*dM*dx*dy*dz;
						Pcon[k][j][i] = -Mipbuff[k][j][i]*pc;
						Tcon[k][j][i] = -Mipbuff[k][j][i]*Tc;



					}
					else
					{
						Mipbuff[k][j][i] = 0.0;
					}

				}
			}END_STD_LOOP

			ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase, &Mipbuff) ; CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->TdM, &Tcon) ; CHKERRQ(ierr);
			ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->PdM, &Pcon) ; CHKERRQ(ierr);

			ierr = MeltExtractionExchangeVolume(jr, iphase, update, actx) ; CHKERRQ(ierr);

			ierr = MeltExtractionInterpMarker(actx, iphase) ; CHKERRQ(ierr);


		}
	}
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,&p); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,&T); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//-------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionExchangeVolume"
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr, PetscInt iphase,PetscInt update,AdvCtx *actx)
{
	FDSTAG       *fs;
	FreeSurf     *surf;
	Discret1D    *dsz;
	Material_t   *phases	; // Phases
	PetscInt     i, j, k, K, sx, sy, sz, nx, ny, nz, L, condition,K1,K2,in;
	Vec          dgmvvec ;
	PetscScalar  bz, ez;
	PetscScalar  IR, dx, dy, dz;
	PetscScalar  ***Thickness, D, D1, ***MohoG,DMin,DMax,vol;
	PetscScalar  *vdgmvvec, *vdgmvvecmerge, ***vdgmvvecmerge2, ***vdgmvvec2, ***Mipbuff;

	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	fs = jr->fs;
	dsz = &fs->dsz;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	phases = jr->dbm->phases;
	IR = phases[iphase].RelInt; // Amount of intrusion that has to be injected within the crust, at level
	surf = jr->surf;

	// Default value

	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
	// Initialize the necessary vector
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &dgmvvec); CHKERRQ(ierr);
	ierr = VecZeroEntries(dgmvvec) ; CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->dgmvvecmerge); CHKERRQ(ierr);


	// Retrieve the vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->Miphase, &Mipbuff); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dgmvvec, &vdgmvvec2); CHKERRQ(ierr);

	// scan all local cells
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		vdgmvvec2[L][j][i] += Mipbuff[k][j][i];
	}
	END_STD_LOOP

	// Restore the vector
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase , &Mipbuff)     ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dgmvvec, &vdgmvvec2)     ; CHKERRQ(ierr);
	// Perform the integral
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(dgmvvec, &vdgmvvec); CHKERRQ(ierr);
		ierr = VecGetArray(jr->dgmvvecmerge, &vdgmvvecmerge); CHKERRQ(ierr);
		ierr = MPI_Allreduce(vdgmvvec, vdgmvvecmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(dgmvvec, &vdgmvvec); CHKERRQ(ierr);
		ierr = VecRestoreArray(jr->dgmvvecmerge, &vdgmvvecmerge); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(dgmvvec,jr->dgmvvecmerge);  CHKERRQ(ierr);
	}


	if(update==1)
	{
		ierr=AverageP_T(jr); CHKERRQ(ierr);
		ierr=OutPutVolume(surf,iphase); CHKERRQ(ierr);
	}

	// Access to Vector


	// Access to the vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->Miphase, &Mipbuff)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, jr->dgmvvecmerge, &vdgmvvecmerge2)   ; CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->Thickness,&Thickness); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->gdMoho,&MohoG); CHKERRQ(ierr);
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);
	// Scan all the cells
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_PLANE_LOOP
	{
		condition = 0;
		D1 =Thickness[L][j][i];
		D = MohoG[L][j][i] + D1*phases[iphase].DInt;
		DMin=D-D1*(phases[iphase].DTol);
		DMax=D+D1*(phases[iphase].DTol);


		if(vdgmvvecmerge2[L][j][i] < 0.0)
		{
			// The grid is assumed to BE uniform.
			if(phases[iphase].DTol>0)
			{
				vol=-(phases[iphase].VolCor*IR*vdgmvvecmerge2[L][j][i])/(DMax-DMin);
				//PetscPrintf(PETSC_COMM_SELF, "IntrusionPhase =%d && iphase =%d,vol/TK %6f = && Tk_INT = %6f\n",phases[iphase].PhExt,iphase,vol*jr->scal->length*jr->scal->length,(DMax-DMin)*jr->scal->length);
				if(D1> (jr->ctrl.MinTk))
				{
				// Are the extreme of the interval in the same processor?
					if(DMin >= bz && DMin < ez && DMax >= bz && DMax < ez)
					{
						// Index Min
						K1 = FindPointInCell(dsz->ncoor, 0, dsz->ncels, DMin);
						K2 = FindPointInCell(dsz->ncoor,0,dsz->ncels,DMax);
						// Are the two point the same point?
						if(K1==K2)
						{
							Mipbuff[sz+K1][j][i] = -(phases[iphase].VolCor*IR*vdgmvvecmerge2[L][j][i]);
						}
						else
						{
							for(in=sz+K1;in<sz+K2;in++)
								{
								dz = SIZE_CELL(in,sz,fs->dsz);
								Mipbuff[in][j][i]=vol*dz;
								}
						}
					}
				//If the first statement are not true, it means that DMax and DMin MAY be in two different processor
					else if(DMin >= bz && DMin < ez )
					{
						K1 = FindPointInCell(dsz->ncoor, 0, dsz->ncels, DMin);
						for(in=sz+K1;in<sz+nz;in++)
						{
							dz = SIZE_CELL(in,sz,fs->dsz);
							Mipbuff[in][j][i]=vol*dz;
						}
					}
					else if(DMax >= bz && DMax < ez)
					{
						K2 = FindPointInCell(dsz->ncoor, 0, dsz->ncels, DMax);
						for(in=sz;in<sz+K2;in++)
						{
							dz = SIZE_CELL(in,sz,fs->dsz);
							Mipbuff[in][j][i]=vol*dz;
						}
					}
					else if (bz >=DMin && bz < DMax &&ez >=DMin && ez < DMax )
					{
						for(in=sz;in<sz+nz;in++)
						{
						dz = SIZE_CELL(in,sz,fs->dsz);
						Mipbuff[in][j][i]=vol*dz;
						}
					}

				}
				else
				{
					condition = 1 ;
				}
			}
			else
			{
			// check whether point belongs to domain
				if(D >= bz && D < ez)
				{
				// find containing cell
					K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, D);
					if(D1> (jr->ctrl.MinTk))
					{
						Mipbuff[sz+K][j][i] = -IR*phases[iphase].VolCor*vdgmvvecmerge2[L][j][i];
					}
					else
					{
						condition = 1 ;
			//			Mipbuff[sz+K][j][i] = 0.0;
					}
				}
			}
			if(update>0)
			{
				// The first time that Exchange volume is called, has as input the mass. The second time
				// it takes into account the volume. In order to retrive the thickness of the "melt extracted"
				// you need to divide for the area.
				dx = SIZE_CELL(i,sx,fs->dsx);
				dy = SIZE_CELL(j,sy,fs->dsy);
				if(condition == 0)
				{
					vdgmvvecmerge2[L][j][i]= -((1-IR)*phases[iphase].VolCor*vdgmvvecmerge2[L][j][i])/(dx*dy);
				}
				else
				{
					vdgmvvecmerge2[L][j][i]= -phases[iphase].VolCor*vdgmvvecmerge2[L][j][i]/(dx*dy);
				}
			}

		}
	}
	END_PLANE_LOOP
	// Restore Vectors
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, jr->dgmvvecmerge, &vdgmvvecmerge2); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase , &Mipbuff) ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,jr->gdMoho,&MohoG) ; CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,jr->Thickness,&Thickness) ; CHKERRQ(ierr);

	// Restore global & not useful vector
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&dgmvvec); CHKERRQ(ierr);

	// Create the local vector
	ierr = VecZeroEntries(jr->ldvecmerge); CHKERRQ(ierr);
	GLOBAL_TO_LOCAL(jr->DA_CELL_2D, jr->dgmvvecmerge, jr->ldvecmerge);

	if(update==1 && IR<1)
	{
		ierr=Extrusion_melt(surf,iphase,actx); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------------//

#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarker"
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx, PetscInt iphase)
{
	//=======================================================================
	// interpolate increments of the history field of melt extraction to markers
	//=======================================================================

	FDSTAG *fs;
	JacRes *jr;
	Marker *P;
	Material_t *phases;
	SolVarCell *svCell;// pointer 2 the Solution variable defined in the cell
	AdvVelCtx vi;
	PetscScalar UP ;
	PetscInt nx, ny, sx, sy, sz;
	PetscInt jj, ID, I, J, K,ninj;
	PetscInt mrk_inj,n_in,newM;
	PetscScalar ***vldc, Dx, Dy, Dz ;


	PetscErrorCode ierr;
	PetscFunctionBegin;
	fs = actx->fs;
	jr = actx->jr;
	phases = jr->dbm->phases;
	mrk_inj = (PetscInt)ceil(actx->nmax+actx->nmin)/2; // Take the average of max and min.
	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;
	// access 1D layouts of global vectors
	ierr = ADVelCreate(actx, &vi);  CHKERRQ(ierr);
	// communicate boundary values
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->Miphase, jr->ldc);
	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);
	ninj = 0;
	for(ID = 0; ID < fs->nCells; ID++)
	{
		svCell	= &jr->svCell[ID];
		GET_CELL_IJK(ID, I, J, K, nx, ny)
		UP = vldc[sz+K][sy+J][sx+I];

		Dx = SIZE_CELL(sx+I,sx,fs->dsx);
		Dy = SIZE_CELL(sy+J,sy,fs->dsy);
		Dz = SIZE_CELL(sz+K,sz,fs->dsz);

		if(UP>0.0)
		{
			UP=UP/(Dx*Dy*Dz);
			ninj += (PetscInt)ceil((UP*mrk_inj)); // Compute the amount of marker that has to be created

		}
	}
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);


	if(ninj>0)
	{ 	// Allocate the memory
		actx->nrecv = ninj;
		ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);
		ierr = PetscMemzero(actx->recvbuf, (size_t)actx->nrecv*sizeof(Marker)); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

		ID=0;
		n_in=0;

		for(ID=0;ID<fs->nCells;ID++)
		{
			GET_CELL_IJK(ID, I, J, K, nx, ny)
			UP = vldc[sz+K][sy+J][sx+I];

			Dx = SIZE_CELL(sx+I,sx,fs->dsx);
			Dy = SIZE_CELL(sy+J,sy,fs->dsy);
			Dz = SIZE_CELL(sz+K,sz,fs->dsz);
			if(UP>0)
			{
				UP=UP/(Dx*Dy*Dz);
				newM = (PetscInt)ceil((UP*mrk_inj));
				ierr = MeltExtractionInject(jr, actx, ID, I, J, K, UP,iphase,newM,n_in);  CHKERRQ(ierr);
				n_in += newM;
			}
		}
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

		ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

		// compute host cells for all the markers
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);
		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);
		// clear
		ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	}

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];
		// get consecutive index of the host cell
		ID = actx->cellnum[jj];
		svCell	= &jr->svCell[ID];
		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		Dx = SIZE_CELL(sx+I,sx,fs->dsx);
		Dy = SIZE_CELL(sy+J,sy,fs->dsy);
		Dz = SIZE_CELL(sz+K,sz,fs->dsz);

		// access buffer
		UP = vldc[sz+K][sy+J][sx+I];
		if(UP<0.0)
		{
			if(P->phase==iphase)
			{
				UP/=(Dx*Dy*Dz*svCell->phRat[iphase]);
				P->Mtot += -UP;  // has to increase by the amount of melt change
				P->Mvol +=  UP;  // has to decrease
				if(P->Mtot > phases[iphase].Mmax)
				{
					P->phase=phases[iphase].PhNext;
					// It is necessary to set to zero the total melt extracted from this particles. Otherwise the melt production is discontinous
					P->Mtot = 0.0;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

	ierr = ADVelDestroy(&vi); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInject"
PetscErrorCode MeltExtractionInject(JacRes *jr,AdvCtx *actx, PetscInt ID, PetscInt I, PetscInt J, PetscInt K, PetscScalar UP, PetscInt iphase, PetscInt newM, PetscInt n_in)//,PetscInt sx,PetscInt sy,PetscInt sz
{ /* 3_b
   */

	PetscInt ipn, n, ninj, PhInject,ii,sind=0;
	PetscScalar xs[3], xe[3], xp[3],*X,dx[3];
	PetscScalar x,sumind=0.0;
	PetscScalar cf_rand;
	PetscRandom rctx;
	FDSTAG *fs;
	Material_t *phases;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = actx->fs;
	phases = jr->dbm->phases; // take the phases


	PhInject = phases[iphase].PhInt;

	// get markers in cell
	n = actx->markstart[ID+1] - actx->markstart[ID];

	// We have not found a marker of the correct phase or there is still melt to be injected
	if( UP > 0)
	{
		// initialize the random number generator
		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);

		// get cell coordinates
		xs[0] = fs->dsx.ncoor[I]; xe[0] = fs->dsx.ncoor[I+1];
		xs[1] = fs->dsy.ncoor[J]; xe[1] = fs->dsy.ncoor[J+1];
		xs[2] = fs->dsz.ncoor[K];xe[2] = fs->dsz.ncoor[K+1];
		for(ipn =n_in; ipn<n_in+newM; ipn++)
		{
			// create random coordinate within this cell
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[0] = (xe[0] - xs[0]) * cf_rand + xs[0];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[1] = (xe[1] - xs[1]) * cf_rand + xs[1];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[2] = (xe[2] - xs[2])* cf_rand + xs[2];
			// calculate the closest (parent marker)
			for (ii = 0; ii < n; ii++)
			{
				X  = actx->markers[ii].X;
				dx[0] = X[0] - xp[0];
				dx[1] = X[1] - xp[1];
				dx[2] = X[2] - xp[2];
				x  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

				if (ii == 0)    { sumind = x; sind = ii; }
				if (x < sumind) { sumind = x; sind = ii; }
			}

			// create new marker
			// actx->recvbuf[ipn]      = markers[sind];

			// hard-coded new marker properties for debugging
			actx->recvbuf[ipn].phase = PhInject;
			actx->recvbuf[ipn].p = actx->markers[sind].p;
			actx->recvbuf[ipn].T = phases[iphase].TInt;
			actx->recvbuf[ipn].APS = 5;
			actx->recvbuf[ipn].Mtot = 0;

			actx->recvbuf[ipn].S.xx = 0;
			actx->recvbuf[ipn].S.xy = 0;
			actx->recvbuf[ipn].S.xz = 0;
			actx->recvbuf[ipn].S.yy = 0;
			actx->recvbuf[ipn].S.yz = 0;
			actx->recvbuf[ipn].S.zz = 0;
			actx->recvbuf[ipn].Mvol = 1;
			actx->recvbuf[ipn].U[0] = 0;
			actx->recvbuf[ipn].U[1] = 0;
			actx->recvbuf[ipn].U[2] = 0;
			actx->recvbuf[ipn].X[0] = xp[0];
			actx->recvbuf[ipn].X[1] = xp[1];
			actx->recvbuf[ipn].X[2] = xp[2];
		}
		// destroy random context
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

		// store new markers
		actx->ndel = 0;
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Moho_Tracking"
PetscErrorCode Moho_Tracking(JacRes *jr)
{
	FDSTAG *fs;
	Discret1D *dsz;
	PetscInt i, j, k,numPhases,ii;
	PetscInt sx, sy, sz, nx, ny, nz, iter,L;
	PetscScalar bz, ez, bottom;
	PetscScalar ***Mohovec2,*Mohovec22,*MMerge ,MantP;
	Material_t *phases; // Phases
	PetscScalar *phRat;

	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	fs = jr->fs;
	dsz = &fs->dsz;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	phases = jr->dbm->phases;
	numPhases = jr->dbm->numPhases; // take the number of phases from dbm structures

	// Find the bottom of the domain
	ierr = FDSTAGGetGlobalBox(fs, 0, 0, &bottom, 0, 0, 0); CHKERRQ(ierr);
	// get local coordinate bounds
	ierr = VecSet   (jr->gdMoho1,bottom);	 CHKERRQ(ierr);
	ierr = VecZeroEntries   (jr->gdMoho);	 CHKERRQ(ierr);
	ierr = DMDAVecGetArray  (jr->DA_CELL_2D,jr->gdMoho1, &Mohovec2);	CHKERRQ(ierr);
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

	ierr =DMDAVecRestoreArray(jr->DA_CELL_2D,jr->gdMoho1,&Mohovec2); CHKERRQ(ierr);

	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(jr->gdMoho1, &Mohovec22)           ; CHKERRQ(ierr);
		ierr = VecGetArray(jr->gdMoho, &MMerge) ; CHKERRQ(ierr);
		ierr = MPI_Allreduce(Mohovec22, MMerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_MAX, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(jr->gdMoho1, &Mohovec22); CHKERRQ(ierr);
		ierr = VecRestoreArray(jr->gdMoho, &MMerge); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(jr->gdMoho1,jr->gdMoho);  CHKERRQ(ierr);
	}
	// Find the values at the top
	GLOBAL_TO_LOCAL(jr->DA_CELL_2D, jr->gdMoho, jr->ldMoho);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Extrusion_melt"
PetscErrorCode Extrusion_melt(FreeSurf *surf,PetscInt iphase, AdvCtx *actx)
{
	JacRes *jr;
	FDSTAG *fs;
	PetscScalar ***topo,***lmelt;
	PetscScalar zbot, ztop, z, Melt[4], Layer;
	PetscInt L, cnt, gcnt;
	PetscInt i, j, nx, ny, sx, sy, I1, I2, J1, J2, mx, my;
	PetscErrorCode ierr;
	Material_t *phases;
	PetscFunctionBegin;

	// free surface cases only
	if(!surf->UseFreeSurf) PetscFunctionReturn(0);

	// access context
	jr = surf->jr;
	fs = jr->fs;
	L  = (PetscInt)fs->dsz.rank;
	phases = jr->dbm->phases;
	surf->MeltExtraction = 1;

	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;

	// get z-coordinates of the top and bottom boundaries
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, &ztop); CHKERRQ(ierr);

	// store the phase that is being sedimented
	surf->phaseEx = phases[iphase].PhExt;

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, jr->ldvecmerge,  &lmelt);  CHKERRQ(ierr);
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

		Melt[0] = lmelt[L][J1][I1]; if(Melt[0] < 0.0) Melt[0] = 0.0;
		Melt[1] = lmelt[L][J1][I2]; if(Melt[1] < 0.0) Melt[1] = 0.0;
		Melt[2] = lmelt[L][J2][I1]; if(Melt[2] < 0.0) Melt[2] = 0.0;
		Melt[3] = lmelt[L][J2][I2]; if(Melt[3] < 0.0) Melt[3] = 0.0;

		Layer = (Melt[0] + Melt[1] + Melt[2] + Melt[3])/4;
		// get topography
		z = topo[L][j][i];
		// uniformly advect
		z += Layer;
		// check if internal free surface goes outside the model domain
		if(z > ztop) z = ztop;
		if(z < zbot) z = zbot;

		// store advected topography
		topo[L][j][i] = z;

	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, jr->ldvecmerge,  &lmelt);  CHKERRQ(ierr);
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

	// update arrays for marker-cell interaction
	ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);
	ierr = FreeSurfGetAirPhaseRatio(surf); CHKERRQ(ierr);

	// print info

	surf->MeltExtraction = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ExchangeMassME"
PetscErrorCode ExchangeMassME(SolVarBulk *svBulk,PetscScalar dx,PetscScalar dy,PetscScalar dz, PetscScalar dt)
{
	PetscFunctionBegin;
	if(svBulk->dMass!=0.0)
	{
	svBulk->Mass=1/dt*(1-(dx*dy*dz)/(dx*dy*dz+svBulk->dMass));
	}

	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Compute_Thickness"
PetscErrorCode Compute_Thickness(JacRes *jr)
{
	FreeSurf *surf;
	FDSTAG *fs;
	PetscScalar cz[4],***h,***ntopo,z,***mg;
	PetscInt i, j, sx, sy, nx, ny, L;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	surf=jr->surf;
	fs=jr->fs;
	L = (PetscInt)fs->dsz.rank; // rank of the processor
	// Initialize the global vector
	ierr = VecZeroEntries(jr->Thickness); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &ntopo); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->Thickness,&h); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->gdMoho,&mg); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,jr->Thickness,&h); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D,jr->gdMoho,&mg); CHKERRQ(ierr);

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
	if(ctrl->initGuess) PetscFunctionReturn(0);
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
		svCell->svBulk.mfextot = 0.0;
		svCell->svBulk.mfVol  = 0.0;
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
		svCell->svBulk.mfextot += w*P->Mtot;
		svCell->svBulk.mfVol   += w*P->Mvol;

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		// normalize history variables
		svCell->svBulk.mfextot /= w;
		svCell->svBulk.mfVol   /= w;
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

		UPXY = P->Mtot;  UPXZ = P->Mtot;  UPYZ = P->Mtot;

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
	for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].svDev.mfextot = gxy[jj]/jr->svXYEdge[jj].ws;
	for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].svDev.mfextot = gxz[jj]/jr->svXZEdge[jj].ws;
	for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].svDev.mfextot = gyz[jj]/jr->svYZEdge[jj].ws;

	// restore access
	ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutPutVolume"
PetscErrorCode OutPutVolume(FreeSurf *surf,PetscInt iphase)
{
	JacRes *jr;
	FDSTAG *fs;
	Vec Crust_produced,MeanT,MeanP,BufferTemperature,BufferPressure;
	PetscScalar ***lmelt,***t2d,***p2d,***t,***p,***bt,***bp,***vcm,***vcc,***totalTMafcts,***totalTConcts,***totalPMafcts,***totalPConcts;
	PetscScalar zbot, ztop, Melt[4], T[4],P[4],***Layer,***MC,***CC;
	PetscScalar wt[4],wp[4];
	PetscInt L;
	PetscInt i, j, nx, ny, sx, sy, I1, I2, J1, J2, mx, my;
	Material_t *phases;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	// free surface cases only

	// access context
	jr = surf->jr;
	fs = jr->fs;
	L  = (PetscInt)fs->dsz.rank;
	phases = jr->dbm->phases;
	// Construct the local vector, being sure that represents the actual volume
	ierr = VecZeroEntries(jr->ldvecmerge); CHKERRQ(ierr);
	GLOBAL_TO_LOCAL(jr->DA_CELL_2D, jr->dgmvvecmerge, jr->ldvecmerge);

	// Create buffer Vector and set it equal to zero
	ierr = DMGetGlobalVector(surf->DA_SURF, &Crust_produced); CHKERRQ(ierr);
	ierr = VecZeroEntries(Crust_produced);CHKERRQ(ierr);

	ierr = DMGetGlobalVector(surf->DA_SURF, &MeanT); CHKERRQ(ierr);
	ierr = VecZeroEntries(MeanT);CHKERRQ(ierr);

	ierr = DMGetGlobalVector(surf->DA_SURF, &MeanP); CHKERRQ(ierr);
	ierr = VecZeroEntries(MeanP);CHKERRQ(ierr);

	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;

	// get z-coordinates of the top and bottom boundaries
	ierr = FDSTAGGetGlobalBox(fs, NULL, NULL, &zbot, NULL, NULL, &ztop); CHKERRQ(ierr);

	// store the phase that is being sedimented

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, jr->ldvecmerge,  &lmelt);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, Crust_produced,  &Layer);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, MeanT,  &t);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, MeanP,  &p);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->TC2D ,&t2d);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D,jr->PC2D ,&p2d);  CHKERRQ(ierr);




	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	// scan all free surface local points

	START_PLANE_LOOP
	{
		// temperature weight
		wt[0] = 1.0  ;  wt[1] = 1.0  ;  wt[2] = 1.0  ;  wt[3] = 1.0  ;
		// Pressure weight
		wp[0] = 1.0  ;  wp[1] = 1.0  ;  wp[2] = 1.0  ;  wp[3] = 1.0  ;

		I1 = i;
		I2 = i-1;
		J1 = j;
		J2 = j-1;
		// check index bounds if ghost points are undefined
		if(I1 == mx) I1--;
		if(I2 == -1) I2++;
		if(J1 == my) J1--;
		if(J2 == -1) J2++;

		// Crust Produced
		Melt[0] = -lmelt[L][J1][I1]; if(Melt[0] < 0.0) Melt[0] = 0.0;
		Melt[1] = -lmelt[L][J1][I2]; if(Melt[1] < 0.0) Melt[1] = 0.0;
		Melt[2] = -lmelt[L][J2][I1]; if(Melt[2] < 0.0) Melt[2] = 0.0;
		Melt[3] = -lmelt[L][J2][I2]; if(Melt[3] < 0.0) Melt[3] = 0.0;

		Layer[L][j][i] = (Melt[0] + Melt[1] + Melt[2] + Melt[3])/4;



		// Temperature

		T[0] = t2d[L][J1][I1]  ;
		if(T[0]<=0.0)
		{
			T[0] = 0.0  ;
			wt[0] = 0.0  ;
		}
		T[1] = t2d[L][J1][I2]  ;
		if(T[1]<=0.0)
		{
			T[1] = 0.0  ;
			wt[1] = 0.0  ;
		}
		T[2] = t2d[L][J2][I1]  ;
		if(T[2]<=0.0)
		{
			T[2] = 0.0  ;
			wt[2] = 0.0  ;
		}
		T[3] = t2d[L][J2][I2]  ;
		if(T[3]<=0.0)
		{
			T[3] = 0.0  ;
			wt[3] = 0.0  ;
		}

		if(wt[0] == 0.0 && wt[1] == 0.0 && wt[2] == 0.0 && wt[3] == 0.0)
		{
		t[L][j][i] = 0.0 ;
		}
		else
		{
		t[L][j][i] = (wt[0]*T[0]+wt[1]*T[1]+wt[2]*T[2]+wt[3]*T[3])/(wt[0]+wt[1]+wt[2]+wt[3]);
		}


		P[0] = p2d[L][J1][I1]  ;
		if(P[0]<=0.0)
		{
			P[0] = 0.0  ;
			wp[0] = 0.0  ;
		}
		P[1] = p2d[L][J1][I2]  ;
		if(P[1]<=0.0)
		{
			P[1] = 0.0  ;
			wp[1] = 0.0  ;
		}
		P[2] = p2d[L][J2][I1]  ;
		if(P[2]<=0.0)
		{
			P[2] = 0.0  ;
			wp[2] = 0.0  ;
		}
		P[3] = p2d[L][J2][I2]  ;
		if(P[3]<=0.0)
		{
			P[3] = 0.0  ;
			wp[3] = 0.0  ;
		}


		if(wp[0] == 0.0 && wp[1] == 0.0 && wp[2] == 0.0 && wp[3] == 0.0)
		{
		p[L][j][i] = 0.0 ;
		}
		else
		{
		p[L][j][i] = (wp[0]*P[0]+wp[1]*P[1]+wp[2]*P[2]+wp[3]*P[3])/(wp[0]+wp[1]+wp[2]+wp[3]);

		}

	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, jr->ldvecmerge,  &lmelt);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, Crust_produced,  &Layer);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF,MeanT,&t); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF,MeanP,&p); CHKERRQ(ierr);

	// Create Temporary Buffer for pressure and temperature
	ierr = DMGetGlobalVector(surf->DA_SURF, &BufferTemperature); CHKERRQ(ierr);
	ierr = VecZeroEntries(BufferTemperature);CHKERRQ(ierr);

	ierr = DMGetGlobalVector(surf->DA_SURF, &BufferPressure); CHKERRQ(ierr);
	ierr = VecZeroEntries(BufferPressure);CHKERRQ(ierr);

	// Save 2 Surf structure for later visualization
	if(phases[iphase].pMant==1)
	{

		ierr = DMDAVecGetArray(surf->DA_SURF, Crust_produced,  &Layer)     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, MeanT,           &t    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, MeanP,           &p    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->NewMafic,  &MC)        ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, BufferTemperature,&bt)       ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, BufferPressure,&bp)          ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->VcurrentMaf,&vcm)      ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->TtotMaf,&totalTMafcts) ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->PtotMaf,&totalPMafcts) ;   CHKERRQ(ierr);

		ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
		START_PLANE_LOOP
		{
		MC[L][j][i]  += Layer[L][j][i];
		vcm[L][j][i] += Layer[L][j][i];
		totalTMafcts[L][j][i] += Layer[L][j][i]*t[L][j][i];
		totalPMafcts[L][j][i] += Layer[L][j][i]*p[L][j][i];

		if(vcm[L][j][i]==0.0)
		{
			bt[L][j][i] = 0.0;
			bp[L][j][i] = 0.0;
		}
		else
		{
			bt[L][j][i] = totalTMafcts[L][j][i]/vcm[L][j][i];
			bp[L][j][i] = totalPMafcts[L][j][i]/vcm[L][j][i];

		}
		}
		END_PLANE_LOOP
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->NewMafic, surf->lNewMafic);
		ierr = VecZeroEntries(surf->TProductionMaf);     CHKERRQ(ierr);
		GLOBAL_TO_LOCAL(surf->DA_SURF, BufferTemperature, surf->TProductionMaf);
		ierr = VecZeroEntries(surf->PProductionMaf);     CHKERRQ(ierr);
		GLOBAL_TO_LOCAL(surf->DA_SURF, BufferPressure, surf->PProductionMaf);


		ierr = DMDAVecRestoreArray(surf->DA_SURF, Crust_produced,  &Layer)     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, MeanT,           &t    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, MeanP,           &p    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->NewMafic,  &MC)        ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, BufferTemperature,&bt)       ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, BufferPressure,   &bp)       ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->VcurrentMaf,&vcm)      ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->TtotMaf,&totalTMafcts) ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->PtotMaf,&totalPMafcts) ;   CHKERRQ(ierr);



	}
	else if(phases[iphase].pMc==1)
	{

		ierr = DMDAVecGetArray(surf->DA_SURF, Crust_produced,  &Layer)     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, MeanT,           &t    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, MeanP,           &p    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->NewContinental,  &CC)        ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, BufferTemperature,&bt)       ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, BufferPressure,&bp)       ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->VcurrentCon,&vcc)      ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->TtotCon,&totalTConcts) ;   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->PtotCon,&totalPConcts) ;   CHKERRQ(ierr);


		ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
		START_PLANE_LOOP
		{
			CC[L][j][i]  += Layer[L][j][i];
			vcc[L][j][i] += Layer[L][j][i];
			totalTConcts[L][j][i] += Layer[L][j][i]*t[L][j][i];
			totalPConcts[L][j][i] += Layer[L][j][i]*p[L][j][i];



			if(vcc[L][j][i]==0.0)
			{
				bt[L][j][i] = 0.0;
				bp[L][j][i] = 0.0;
			}
			else
			{
				bt[L][j][i] = totalTConcts[L][j][i]/vcc[L][j][i];
				bp[L][j][i] = totalPConcts[L][j][i]/vcc[L][j][i];

			}
		}
		END_PLANE_LOOP
		GLOBAL_TO_LOCAL(surf->DA_SURF, surf->NewContinental, surf->lNewContinental);

		ierr = VecZeroEntries(surf->TProductionCon);     CHKERRQ(ierr);
		GLOBAL_TO_LOCAL(surf->DA_SURF, BufferTemperature, surf->TProductionCon);
		ierr = VecZeroEntries(surf->PProductionCon);     CHKERRQ(ierr);
		GLOBAL_TO_LOCAL(surf->DA_SURF, BufferPressure, surf->PProductionCon);

		ierr = DMDAVecRestoreArray(surf->DA_SURF, Crust_produced,  &Layer)     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, MeanT,           &t    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, MeanP,           &p    )     ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->NewContinental,  &CC)        ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, BufferTemperature,&bt)       ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, BufferPressure,   &bp)       ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->VcurrentCon,&vcc)      ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->TtotCon,&totalTConcts) ;   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->PtotCon,&totalPConcts) ;   CHKERRQ(ierr);
	}

	ierr = DMRestoreGlobalVector(surf->DA_SURF,&Crust_produced); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(surf->DA_SURF,&MeanT); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(surf->DA_SURF,&MeanP); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(surf->DA_SURF,&BufferTemperature); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(surf->DA_SURF,&BufferPressure); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "AverageP_T"
PetscErrorCode AverageP_T(JacRes *jr)
{
	FDSTAG       *fs;
	Discret1D    *dsz;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, L;
	Vec          TemperaturePart, TemperatureTot ;
	Vec          PressurePart,PressureTot;
	PetscScalar  bz, ez;
	PetscScalar ***tc,***T3D,*tcm,*t,***vol;
	PetscScalar ***pc,***P3D,*pcm,*p;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	fs = jr->fs;
	dsz = &fs->dsz;
	L = (PetscInt)fs->dsz.rank; // rank of the processor

	// Default value

	ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);
	// Temperature Extraction Condition

	// Initialize the necessary vector
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &TemperaturePart); CHKERRQ(ierr);
	ierr = VecZeroEntries(TemperaturePart) ; CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &TemperatureTot); CHKERRQ(ierr);
	ierr = VecZeroEntries(TemperaturePart) ; CHKERRQ(ierr);

	// Retrieve the vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->TdM,&T3D); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, TemperaturePart, &tc); CHKERRQ(ierr);

	// scan all local cells
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		tc[L][j][i] += T3D[k][j][i];


	}
	END_STD_LOOP

	// Restore the vector
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->TdM,&T3D); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, TemperaturePart, &tc); CHKERRQ(ierr);

	// Perform the integral
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(TemperaturePart, &tcm); CHKERRQ(ierr);
		ierr = VecGetArray(TemperatureTot, &t); CHKERRQ(ierr);
		ierr = MPI_Allreduce(tcm, t, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(TemperaturePart, &tcm); CHKERRQ(ierr);
		ierr = VecRestoreArray(TemperatureTot, &t); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(TemperaturePart,TemperatureTot);  CHKERRQ(ierr);
	}

	// Pressure

		// Initialize the necessary vector
		ierr = DMGetGlobalVector(jr->DA_CELL_2D, &PressurePart); CHKERRQ(ierr);
		ierr = VecZeroEntries(TemperaturePart) ; CHKERRQ(ierr);
		ierr = DMGetGlobalVector(jr->DA_CELL_2D, &PressureTot); CHKERRQ(ierr);
		ierr = VecZeroEntries(TemperaturePart) ; CHKERRQ(ierr);

		// Retrieve the vectors
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->PdM,&P3D); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(jr->DA_CELL_2D, PressurePart, &pc); CHKERRQ(ierr);

		// scan all local cells
		ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			pc[L][j][i] += P3D[k][j][i];

		}
		END_STD_LOOP

		// Restore the vector
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->PdM,&P3D); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, TemperaturePart, &pc); CHKERRQ(ierr);

		// Perform the integral
		if(dsz->nproc != 1 )
		{
			ierr = VecGetArray(PressurePart, &pcm); CHKERRQ(ierr);
			ierr = VecGetArray(PressureTot, &p); CHKERRQ(ierr);
			ierr = MPI_Allreduce(pcm, p, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);
			ierr = VecRestoreArray(PressurePart,   &pcm); CHKERRQ(ierr);
			ierr = VecRestoreArray(PressureTot, &p  ); CHKERRQ(ierr);
		}
		else
		{
			ierr = VecCopy(PressurePart,PressureTot);  CHKERRQ(ierr);
		}

	// Computing the actual volume weightened average of P/T condition of extraction


	ierr = DMDAVecGetArray(jr->DA_CELL_2D, TemperatureTot,&tc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, PressureTot,&pc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, jr->dgmvvecmerge, &vol); CHKERRQ(ierr);

	// scan all local cells
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_PLANE_LOOP
	{
		if(vol[L][j][i]<0)
		{
			tc[L][j][i]    /= -vol[L][j][i];
			pc[L][j][i]    /= -vol[L][j][i];


		}
		else
		{
			tc[L][j][i] = 0.0;
			pc[L][j][i] = 0.0;

		}

	}
	END_PLANE_LOOP

	GLOBAL_TO_LOCAL(jr->DA_CELL_2D,TemperatureTot,jr->TC2D)
	GLOBAL_TO_LOCAL(jr->DA_CELL_2D,TemperatureTot,jr->PC2D)

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, TemperatureTot,&tc); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, PressureTot,&pc); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, jr->dgmvvecmerge, &vol); CHKERRQ(ierr);

	// Get Rid of the Buffer
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&TemperatureTot); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&TemperaturePart); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&PressureTot); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&PressurePart); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//======================================================================================================

#undef __FUNCT__
#define __FUNCT__ "Mean_Continental_Crust"
PetscErrorCode Mean_Continental_Crust(JacRes *jr)
{
	Material_t   *phases    ;
	FDSTAG       *fs;
	Discret1D    *dsz;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, L,iter,ii;
	Vec          Rel_Con,Rel_Con_Tot,Total_Ratio,buff;
	PetscScalar  bz, ez;
	PetscScalar ***cc,***cct,***cct1,*c,*ct,***tk,dx,dy,dz,***Buff;
	PetscScalar CC[4];
	PetscInt    I1, I2, J1, J2, mx, my;
	PetscScalar *phRat;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	// access context
	fs = jr->fs;
	dsz = &fs->dsz;
	phases = jr->dbm->phases;
	L = (PetscInt)fs->dsz.rank; // rank of the processor

	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &Rel_Con_Tot); CHKERRQ(ierr);
	ierr = VecZeroEntries(Rel_Con_Tot) ; CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &Rel_Con); CHKERRQ(ierr);
	ierr = VecZeroEntries(Rel_Con) ; CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->surf->DA_SURF, &Total_Ratio); CHKERRQ(ierr);
	ierr = VecZeroEntries(Total_Ratio) ; CHKERRQ(ierr);
	ierr = DMGetLocalVector(jr->DA_CELL_2D, &buff); CHKERRQ(ierr);
	ierr = VecZeroEntries(buff) ; CHKERRQ(ierr);

	ierr=Moho_Tracking(jr); CHKERRQ(ierr);
	ierr=Compute_Thickness(jr); CHKERRQ(ierr);

	// scan all local cells
	ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Rel_Con,&cc); CHKERRQ(ierr);


	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	iter=0;
	START_STD_LOOP
	{
		phRat = jr->svCell[iter++].phRat;

		for(ii=0;ii<jr->dbm->numPhases;ii++)
		{
			if(phRat[ii]>0 && phases[ii].pCc==1)
			{
				dx = SIZE_CELL(i,sx,fs->dsx);
				dy = SIZE_CELL(j,sy,fs->dsy);
				dz = SIZE_CELL(k,sz,fs->dsz);
				cc[L][j][i] += phRat[ii]*dx*dy*dz;
			}
		}

	}
	END_STD_LOOP

	// Restore the vector
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Rel_Con,&cc); CHKERRQ(ierr);

	// Perform the integral
	if(dsz->nproc != 1 )
	{
		ierr = VecGetArray(Rel_Con, &c); CHKERRQ(ierr);
		ierr = VecGetArray(Rel_Con_Tot, &ct); CHKERRQ(ierr);
		ierr = MPI_Allreduce(c, ct, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);
		ierr = VecRestoreArray(Rel_Con, &c); CHKERRQ(ierr);
		ierr = VecRestoreArray(Rel_Con_Tot, &ct); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCopy(Rel_Con,Rel_Con_Tot);  CHKERRQ(ierr);
	}
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&Rel_Con); CHKERRQ(ierr);




	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, Rel_Con_Tot,&cct); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, jr->Thickness,&tk); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		dx = SIZE_CELL(i,sx,fs->dsx);
		dy = SIZE_CELL(j,sy,fs->dsy);
		if(tk[L][j][i]==0.0)
		{
			cct[L][j][i]=0.0;
		}
		else
		{
			cct[L][j][i]/=(tk[L][j][i]*dx*dy);
		}

	}
	END_PLANE_LOOP
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, Rel_Con_Tot,&cct); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, jr->Thickness,&tk); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(jr->DA_CELL_2D,Rel_Con_Tot,buff);


	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, buff,&Buff); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->surf->DA_SURF,Total_Ratio ,&cct1); CHKERRQ(ierr);

	START_PLANE_LOOP
	{
		I1 = i;
		I2 = i-1;
		J1 = j;
		J2 = j-1;
		// check index bounds if ghost points are undefined
		if(I1 == mx) I1--;
		if(I2 == -1) I2++;
		if(J1 == my) J1--;
		if(J2 == -1) J2++;
		// Crust Produced
		CC[0] = Buff[L][J1][I1]; if(CC[0] < 0.0) CC[0] = 0.0;
		CC[1] = Buff[L][J1][I2]; if(CC[1] < 0.0) CC[1] = 0.0;
		CC[2] = Buff[L][J2][I1]; if(CC[2] < 0.0) CC[2] = 0.0;
		CC[3] = Buff[L][J2][I2]; if(CC[3] < 0.0) CC[3] = 0.0;

		cct1[L][j][i] = (CC[0] + CC[1] + CC[2] + CC[3])/4;
	}
	END_PLANE_LOOP

	ierr = VecZeroEntries(jr->surf->R_Cont) ; CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(jr->surf->DA_SURF,Total_Ratio,jr->surf->R_Cont)

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, buff,&Buff); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->surf->DA_SURF,Total_Ratio ,&cct1); CHKERRQ(ierr);

	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&Rel_Con_Tot); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(jr->DA_CELL_2D,&buff); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D,&Total_Ratio); CHKERRQ(ierr);



	PetscFunctionReturn(0);

}
