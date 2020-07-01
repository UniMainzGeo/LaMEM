
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
#include "constEq.h"

#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(ConstEqCtx  *ctx,PetscInt ph)
{
	Material_t *mat;
	Ph_trans_t *PhaseTrans;
	PetscInt    itr,nPtr,id;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if(ctx->ctrl->initGuess==1 || ctx->ctrl->Phasetrans==0 ) PetscFunctionReturn(0);
	PhaseTrans=ctx->PhaseTrans;
	mat = &ctx->phases[ph];
	nPtr=mat->nPTr;

	itr = 0;
	for(itr=0;itr<nPtr;itr++)
		{
			id=mat->Ph_tr[itr];


			ierr = Transition(PhaseTrans, ctx, id); CHKERRQ(ierr);

		}
	PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "Transition"
PetscErrorCode Transition(Ph_trans_t *PhaseTrans, ConstEqCtx  *ctx, PetscInt id)
{
	PetscInt    neq, ip;
	PetscScalar Pres[2];
	Ph_trans_t *PTr;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	PTr=PhaseTrans+id;
	if(PTr->visual==1) 	PetscFunctionReturn(0);

	if(!strcmp(PTr->Type,"Constant"))
	{
		if(!strcmp(PTr->Parameter,"T"))
		{
			if(PTr->value[0]>0)
			{
				if(ctx->T>PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}

			}
			else
			{
				if(ctx->T<PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}

			}
		}


		if(!strcmp(PTr->Parameter,"p"))
		{
			if(PTr->value[0]>0)
			{
				if(ctx->p_lith >PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}
			}
			else
			{
				if(ctx->p_lith<PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}
			}
		}


		if(!strcmp(PTr->Parameter,"Depth"))
		{
			if(PTr->value[0]>0)
			{
				if(ctx->depth>PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}
			}
			else
			{
				if(ctx->depth<PTr->value[1])
				{
					ctx->visc_inc=PTr->visc_inc;
					ctx->rho_inc=PTr->rho_inc;
				}
			}
		}
	}
	else if(!strcmp(PTr->Type,"Clapeyron"))
	{
		neq = PTr->neq;
		for (ip=0;ip<neq;ip++)
		{
			Pres[ip]=(ctx->T-PTr->T0[ip])*PTr->gamma[ip]+PTr->P0[ip];
		}
		if(neq==1)
		{
			if(ctx->p_lith>Pres[0])
			{
				ctx->visc_inc=PTr->visc_inc;
				ctx->rho_inc=PTr->rho_inc;
			}
		}
		else
		{
			if(ctx->p_lith>Pres[0] && ctx->p>Pres[1])
			{
				ctx->visc_inc=PTr->visc_inc;
				ctx->rho_inc=PTr->rho_inc;
			}
		}
	}
	PetscFunctionReturn(0);
}
