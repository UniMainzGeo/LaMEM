
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
 **    filename:   AVD.c
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
//........  .......
//---------------------------------------------------------------------------

// Routine that changes the phase as a function of the condition set by the user.
 // These phase transitions are extremely simplified wtr to the phase transition listed in
 // phase diagram. And allow to consider lower mantle, ast->lithosphere, basalt->eclogitic reaction
 //

#include "LaMEM.h"
#include "AVD.h"
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "fdstag.h"
#include "bc.h"
#include "tools.h"
#include "phase_transition.h"
#include "phase.h"
//---------------------------------------------------------------------------
// Simple function that assign the primordial phase to the marker
#undef __FUNCT__
#define __FUNCT__ "PhTr_assign_primPh"
PetscErrorCode PhTr_assign_primPh(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction
	PetscInt     i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over all local particles
	for(i = 0; i < actx->nummark; i++)
	{
		actx->markers[i].primph = actx->markers[i].phase;

	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Simple function that assign the primordial phase to the marker
#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction
	DBMat      *dbm;
	Material_t *mat;
	Ph_trans_t *PhaseTrans;
	Marker *P;
	JacRes *jr;
	PetscInt     i, ph, counter,itr,nPtr,id,newPh;

	jr = actx->jr;
	dbm = jr->dbm;
	PhaseTrans = jr->dbm->matPhtr;
	mat = dbm->phases;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over all local particles 		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);

	for(i = 0; i < actx->nummark; i++)
	{
		if(P->phase==0)
		{
			PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
		}
		P=&actx->markers[i];

		if(mat[P->primph].nPTr)
		{

			ph=P->primph;

			counter = mat[P->primph].nPTr;

			while (counter!=0)
			{
				if(mat[ph].nPTr!=0)
				{

					nPtr=mat[ph].nPTr;

					for(itr=0;itr<nPtr;itr++)
					{

						id=mat[ph].Ph_tr[itr];

						newPh = Transition(PhaseTrans, P, id, ph);

						if(newPh != ph)
						{

							ph=newPh;
							counter=mat[ph].nPTr;

							break;
						}
						else
						{
						counter--;
						}
					}
				}
				else
				{
				counter = 0;
				}
			}
			P->phase=ph;

		}


	}
	ierr = ADVProjHistMarkToGrid(actx);



	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt id,PetscInt PH)
{
	PetscInt    neq, ip,NP;
	PetscScalar Pres[2];
	Ph_trans_t *PTr;

	PTr=PhaseTrans+id;
	if(PTr->Type==1)
	{
		if(PTr->Parameter==1)
		{
			if(PTr->value[0]>0)
			{
				if(P->T>PTr->value[1])
				{
					 NP = PTr->Ph2Change;
				}
				else
				{
					 NP= PH;
				}
			}
			else
			{
				if(P->T<PTr->value[1])
				{
					NP = PTr->Ph2Change;
				}
				else
				{
					 NP = PH;
				}
			}
		}


		if(PTr->Parameter==2)
		{
			if(PTr->value[0]>0)
			{
				if(P->p>PTr->value[1])
				{
					 NP = PTr->Ph2Change;
				}
				else
				{
					 NP= PH;
				}
			}
			else
			{
				if(P->p<PTr->value[1])
				{
					 NP = PTr->Ph2Change;
				}
				else
				{
					 NP = PH;
				}
			}
		}


		if(PTr->Parameter==3)
		{
			if(PTr->value[0]>0)
			{
				if(P->X[2]>PTr->value[1])
				{
					NP = PhaseTrans->Ph2Change;
				}
				else
				{
					 NP= PH;
				}
			}
			else
			{
				if(P->X[2]<PTr->value[1])
				{
					 NP = PTr->Ph2Change;
				}
				else
				{
					 NP = PH;
				}
			}
		}
	}
	else if(PTr->Type==2)
	{
		neq = PTr->neq;
		for (ip=0;ip<neq;ip++)
		{
			Pres[ip]=(P->T-PTr->T0[ip])*PTr->gamma[ip]+PTr->P0[ip];
		}
		if(neq==1)
		{
			if(P->p>Pres[0])
			{
				NP=PTr->Ph2Change;
			}
			else
			{
				 NP = PH;
			}
		}
		else
		{
			if(P->p>Pres[0] && P->p>Pres[1])
			{
				 NP=PTr->Ph2Change;
			}
			else
			{
				 NP = PH;
			}
		}
	}
	return NP;
}



/*
 * / Simple function that assign the primordial phase to the marker
#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction
	DBMat      *dbm;
	Material_t *mat;
	Ph_trans_t *PhaseTrans;
	Marker *P;
	JacRes *jr;
	PetscInt     i, ph, counter,itr,nPtr,id,newPh;

	jr = actx->jr;
	dbm = jr->dbm;
	PhaseTrans = jr->dbm->matPhtr;
	mat = dbm->phases;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over all local particles 		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);

	for(i = 0; i < actx->nummark; i++)
	{
		if(P->phase==0)
		{
			PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
		}
		P=&actx->markers[i];

		if(mat[P->primph].nPTr)
		{

			ph=P->primph;

			counter = mat[P->primph].nPTr;

			while (counter!=0)
			{
				if(mat[ph].nPTr!=0)
				{

					nPtr=mat[ph].nPTr;

					for(itr=0;itr<nPtr;itr++)
					{

						id=mat[ph].Ph_tr[itr];

						newPh = Transition(PhaseTrans, P, id, ph);

						if(newPh != ph)
						{

							ph=newPh;
							counter=mat[ph].nPTr;

							break;
						}
						else
						{
						counter--;
						}
					}
				}
				else
				{
				counter = 0;
				}
			}

		}

		P->phase=ph;

	}
	ierr = ADVProjHistMarkToGrid(actx);



	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(AdvCtx *actx, Marker *P)
{
	// creates arrays to optimize marker-cell interaction
	DBMat      *dbm;
	Material_t *mat;
	Ph_trans_t *PhaseTrans;
	JacRes *jr;
	PetscInt     i, ph, counter,itr,nPtr,id,newPh,primPh;

	jr = actx->jr;
	dbm = jr->dbm;
	PhaseTrans = jr->dbm->matPhtr;
	mat = dbm->phases;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over all local particles 			PetscPrintf(PETSC_COMM_WORLD,"QU1 \n");	PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
		primPh = P->primph;

		nPtr = mat[primPh].nPTr;
		if(nPtr>0)
		{

			ph=primPh;

			counter = nPtr;

			while (counter!=0)
			{
				if(mat[ph].nPTr!=0)
				{

					nPtr=mat[ph].nPTr;

					for(itr=0;itr<nPtr;itr++)
					{

						id=mat[ph].Ph_tr[itr];

						newPh = Transition(PhaseTrans, P, id, ph);

						if(newPh != ph)
						{

							ph=newPh;
							counter=mat[ph].nPTr;

							break;
						}
						else
						{
						counter--;
						}
					}
				}
				else
				{
				counter = 0;
				}
			}
			P->phase=ph;
		}



	//ierr = ADVProjHistMarkToGrid(actx);



	PetscFunctionReturn(0);
}

 */

