
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
#include "parsing.h"
#include "objFunct.h"
//-----------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "DBMatReadPhaseTr"
PetscErrorCode DBMatReadPhaseTr(DBMat *dbm, FB *fb)
{
	// read softening law from file

	Ph_trans_t   *ph;
	Scaling      *scal;
	PetscInt  ID;
	PetscInt  it;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// softening law ID
	ierr = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbm->numPhtr-1); CHKERRQ(ierr);

	// get pointer to specified softening law
	ph = dbm->matPhtr + ID;
	// check ID
	if(ph->ID != -1)
	{
		 SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate phase transition law!");
	}

	// set ID
	ph->ID = ID;

	ierr = getStringParam(fb, _REQUIRED_, "Type", ph->Type,0);  CHKERRQ(ierr);
	if (!ph->Type)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type (Constant) (Clapeyron) ", (LLD)ID);
	}

	if(!strcmp(ph->Type,"Constant"))
	{
		ierr= Set_Constant_Phase_Transition(ph, dbm, fb,ID); CHKERRQ(ierr);
	}
	if(!strcmp(ph->Type,"Clapeyron"))
	{
		ierr= Set_Clapeyron_Phase_Transition(ph, dbm, fb,ID); CHKERRQ(ierr);
	}

	if(!strcmp(ph->Type,"Box_type"))
	{
		ierr= Set_Box_Within_Transition(ph, dbm, fb,ID); CHKERRQ(ierr);
	}

	ierr = getIntParam(fb, _OPTIONAL_, "number_phases", &ph->number_phases,1 , _max_num_tr_); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "PhaseBelow", ph->PhaseBelow,ph->number_phases , _max_num_phases_); CHKERRQ(ierr);
	ierr = getIntParam(fb, _OPTIONAL_, "PhaseAbove", ph->PhaseAbove,ph->number_phases , _max_num_phases_); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "DensityBelow", ph->DensityBelow,ph->number_phases , 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "DensityAbove", ph->DensityAbove,ph->number_phases, 1.0); CHKERRQ(ierr);


	if (!ph->PhaseAbove || !ph->PhaseBelow)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type (Constant) (Clapeyron) ", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"   Phase Above %d \n", (LLD)(ph->ID), ph->PhaseAbove);
	PetscPrintf(PETSC_COMM_WORLD,"   Phase Below %d \n", (LLD)(ph->ID), ph->PhaseBelow);


	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Constant_Phase_Transition"
PetscErrorCode  Set_Constant_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb,PetscInt ID)
{
	Scaling      *scal;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm -> scal;

	ierr = getStringParam(fb, _REQUIRED_, "Parameter_transition", ph->Parameter_transition, "none");  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "ConstantValue",&ph->ConstantValue, 1,1.0); CHKERRQ(ierr);
	if((!ph->Parameter_transition || !ph->ConstantValue))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using a constant parameter phase transition you need to specify the parameter: p = Pressure, T = temperature, z = depth", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] : Type = Constant , Parameter = %s . The phase is changed at %s  = %f\n", (LLD)(ph->ID), ph->Parameter_transition, ph->ConstantValue);

	if(!strcmp(ph->Parameter_transition,"T")) //Temperature
	{
		ph->ConstantValue=(ph->ConstantValue + scal->Tshift)/scal->temperature;
	}
	if(!strcmp(ph->Parameter_transition,"p"))//Pressure
	{
		ph->ConstantValue /= scal->stress_si;
	}
	if(!strcmp(ph->Parameter_transition,"Depth"))//Depth
	{
		ph->ConstantValue= scal->length;
	}
	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Constant_Phase_Transition"
PetscErrorCode  Set_Clapeyron_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb, PetscInt ID)
{
	Scaling      *scal;
	PetscInt     it=0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm -> scal;

	ierr = getStringParam(fb, _OPTIONAL_, "Name_clapeyron", ph->Name_clapeyron, "none");  CHKERRQ(ierr);
	if (ph->Name_clapeyron)
	{
		ierr = SetClapeyron_Eq(ph); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] : Type = Clapeyron, Transition law %s\n", (LLD)(ph->ID), ph->Name_clapeyron);


	}

	ierr = getIntParam(fb, _OPTIONAL_, "numberofequation", &ph->neq, 1, 2.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "clapeyron_slope", ph->clapeyron_slope, ph->neq,1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "clapeyron_slope", ph->P0_clapeyron, ph->neq,1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "T0_clapeyron", ph->T0_clapeyron, ph->neq,1.0); CHKERRQ(ierr);

	if((!ph->clapeyron_slope || !ph->T0_clapeyron || !ph->clapeyron_slope||!ph->Name_clapeyron))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using any clapeyron phase transition avaiable you need to specify P0, T0, gamma and the number of equation (P=(T-T0)*gamma+(P0)).", (LLD)ID);
	}


	if(!strcmp(ph->Type,"Clapeyron"))
	{
	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] : Type = Clapeyron, gamma = %2f [MPa], P0 = %f [Pa],T0 = %f [deg C] \n", (LLD)(ph->ID), ph->Type, ph->clapeyron_slope, ph->P0_clapeyron,ph->T0_clapeyron);
	}
	for(it==0;it< ph->neq;it++)
	{
		ph->clapeyron_slope[it] *= 1e6*(scal->temperature/scal->stress_si);
		ph->P0_clapeyron[it]/= (scal->stress_si);
		ph->T0_clapeyron[it]=(ph->T0_clapeyron[it]+scal->Tshift)/ (scal->temperature);
	}
	PetscFunctionReturn(0);

}
//=========================================================================================//
#undef __FUNCT__
#define __FUNCT__ "Set_Box_Within_Transition"
PetscErrorCode  Set_Box_Within_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb,PetscInt ID)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	Scaling      *scal;
	PetscInt      it;

	scal = dbm -> scal;

	ierr = getScalarParam(fb, _REQUIRED_, "Geometrical_box", ph->Geometric_box,6, 1.0); CHKERRQ(ierr);




	for(it=0;it<6;it++)
	{
		ph->Geometric_box[it]/=scal->length;
	}

	ierr = getScalarParam(fb, _REQUIRED_, "DeltaT_within", &ph->dT_within,1, 1.0); CHKERRQ(ierr);


	ph->dT_within=(ph->dT_within+ scal->Tshift)+scal->temperature;

	PetscFunctionReturn(0);

}
// ---------------------------------------------------------------------------------------------------------- //
#undef __FUNCT__
#define __FUNCT__ "Overwrite_Density"
PetscErrorCode  Overwrite_density(DBMat *dbm)
{
	Scaling      *scal;
	Ph_trans_t    *ph;
	Material_t   *mat;
	PetscInt    numPhTrn, nPtr,num_phases,iter, jj1, jj2;
	PetscScalar rho_above,rho_below,rho_scal;

	PetscErrorCode ierr;
	scal = dbm->scal;
	rho_scal = scal->density;
	mat = dbm->phases;
	numPhTrn= dbm->numPhtr;


	for(nPtr=0;nPtr<numPhTrn;nPtr++)
	{
		ph = dbm->matPhtr+nPtr;
		num_phases = ph->number_phases;
		for(iter=0;iter<num_phases;iter++)
		{
			rho_above = ph->DensityAbove[iter];
			rho_below = ph->DensityBelow[iter];
			if(rho_above > 0 && rho_below >0)
			{
				jj1 =ph->PhaseBelow[iter];
				mat[jj1].rho= rho_below/rho_scal;
				PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  density is  %2f\n",jj1,rho_below);
				jj2 =ph->PhaseAbove[iter];
				mat[jj2].rho= rho_above/rho_scal;
				PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  density is  %2f\n",jj2,rho_above);

			}

		}
	}

	PetscFunctionBegin;

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "SetClapeyron_Eq"
PetscErrorCode SetClapeyron_Eq(Ph_trans_t *ph)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if (!strcmp(ph->Name_clapeyron,"Eclogite"))
	{
		// Source Faccenda and Dal Zilio et al 2017 [Hacker et al. (2003) Morb+H2O]
		ph->neq          = 2;
		ph->P0_clapeyron[0] = 2e9;
		ph->P0_clapeyron[1] = 2e9;
		ph->T0_clapeyron[0] = 800;
		ph->T0_clapeyron[1] = 700;
		ph->clapeyron_slope[0] = 1.5;
		ph->clapeyron_slope[1] = -30;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_410km"))
	{
		// Source Faccenda and Dal Zilio et al 2017 [Olivine-Wadseylite Phase transition anhydrous Pyrolite+H2O Litasov and Ohtani (2003)]
		ph->neq = 1;
		ph->P0_clapeyron[0] = 13.5e9;
		ph->T0_clapeyron[0] = 1537;
		ph->clapeyron_slope[0]= 5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_510km"))
	{
		// Source Faccenda and Dal Zilio et al 2017 [WadsleyiteRingwooditePhase transition anhydrous Fo100(theoretical) Hernandez et al 2015]
		ph->neq = 1;
		ph->P0_clapeyron[0] = 18e9;
		ph->T0_clapeyron[0] = 1597;
		ph->clapeyron_slope[0]=3.5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_660km"))
	{
		// Source Faccenda and Dal Zilio et al 2017 [Post Spinel Phase transition anhydrous Pyrolite Ye et al 2014]
		ph->neq = 1;
		ph->P0_clapeyron[0] = 23e9;
		ph->T0_clapeyron[0] = 1667;
		ph->clapeyron_slope[0]=-2.5;
	}

	PetscFunctionReturn(0);
}
//===========================================================================================================//
#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction
	DBMat      *dbm;
	Ph_trans_t *PhaseTrans;
	Marker *P;
	JacRes *jr;
	PetscInt     i, ph,nPtr, numPhTrn,below,above,num_phas,outside;
	PetscInt     PH1,PH2;

	jr = actx->jr;
	dbm = jr->dbm;
	numPhTrn= dbm->numPhtr;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// loop over all local particles 		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
	nPtr=0;
	for(nPtr=0;nPtr<numPhTrn;nPtr++)
	{
		PhaseTrans = jr->dbm->matPhtr+nPtr;
		for(i = 0; i < actx->nummark; i++)
		{

			P=&actx->markers[i];
			if(!strcmp(PhaseTrans->Type,"Box_type"))
			{

				num_phas=PhaseTrans->number_phases;

				outside  = Check_Phase_above_below(PhaseTrans->PhaseBelow,P,num_phas);


				if(outside>=0)
				{

					PH1 = PhaseTrans->PhaseAbove[outside];
					ph = Transition(PhaseTrans, P, PH1, PH2,nPtr);
					P->phase=ph;
				}
			}
			else
			{
				num_phas=PhaseTrans->number_phases;
				below  = Check_Phase_above_below(PhaseTrans->PhaseBelow,P,num_phas);
				above  = Check_Phase_above_below(PhaseTrans->PhaseAbove,P,num_phas);

				if(below >= 0 || above >= 0)
				{
					if(below>=0)
					{
						PH1 = PhaseTrans->PhaseBelow[below];
						PH2 = PhaseTrans->PhaseAbove[below];
					}
					else if (above >=0)
					{
						PH1 = PhaseTrans->PhaseBelow[above];
						PH2 = PhaseTrans->PhaseAbove[above];
					}
					ph = Transition(PhaseTrans, P, PH1, PH2,nPtr);
					P->phase=ph;

				}
			}
		}
		ierr = ADVInterpMarkToCell(actx);
		}


	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1,PetscInt PH2,PetscInt ID)
{
	PetscInt ph;

	if(!strcmp(PhaseTrans->Type,"Constant"))
	{
		ph = Check_Constant_Phase_Transition(PhaseTrans,P,PH1,PH2,ID);
	}
	else if(!strcmp(PhaseTrans->Type,"Clapeyron"))
	{
		ph = Check_Clapeyron_Phase_Transition(PhaseTrans,P,PH1,PH2,ID);
	}
	else
	{
		ph = Check_Constant_Box_Transition(PhaseTrans,P,PH1, ID);
	}


	return ph;
}

//===========================================================================================================//
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,PetscInt ID) // softening parameter
{
	PetscInt ph;


	if(!strcmp(PhaseTrans->Parameter_transition,"T"))
		{
			if(P->T>=PhaseTrans->ConstantValue)
			{
				ph = PH2;
			}
			else
			{
				ph = PH1;
			}
		}

	if(!strcmp(PhaseTrans->Parameter_transition,"p"))
		{
			if(P->p>=PhaseTrans->ConstantValue)
			{
				ph = PH2;
			}
			else
			{
				ph = PH1;
			}
		}

	if(!strcmp(PhaseTrans->Parameter_transition,"Depth"))
		{
			if(P->X[2]>=PhaseTrans->ConstantValue)
			{
				ph = PH2;
			}
			else
			{
				ph = PH1;
			}
		}

	return ph;
}
//==================================================================================================================//
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,PetscInt ID) // softening parameter
{
	PetscInt ph,ip,neq;
	PetscScalar Pres[2];



	neq = PhaseTrans->neq;
	for (ip=0;ip<neq;ip++)
	{
		Pres[ip]=(P->T-PhaseTrans->T0_clapeyron[ip])*PhaseTrans->clapeyron_slope[ip]+PhaseTrans->P0_clapeyron[ip];
	}
	if(neq==1)
	{
		if(P->p>=Pres[0])
		{
			ph=PH2;
		}
		else
		{
			 ph=PH1;
		}
	}
	else
	{
		if(P->p>=Pres[0] && P->p>=Pres[1])
		{
			ph=PH2;
		}
		else
		{
			ph=PH1;
		}
	}

	return ph;
}
//=========================================================================================================== //

PetscInt Check_Constant_Box_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1,PetscInt ID) // softening parameter
{
	PetscInt ph;
	PetscScalar X[3];


		X[0]=P->X[0];
		X[1]=P->X[1];
		X[2]=P->X[2];


		if(X[0]>=PhaseTrans->Geometric_box[0] && X[0]<=PhaseTrans->Geometric_box[1] && X[1]>=PhaseTrans->Geometric_box[2] && X[1]<=PhaseTrans->Geometric_box[3] && X[2]>=PhaseTrans->Geometric_box[4] && X[2]<=PhaseTrans->Geometric_box[5])
		{

			ph = PH1;
			P->T = PhaseTrans->dT_within;

		}
		else
		{
			ph = P->phase;
		}




	return ph;
}

// ========================================================================================================== //
PetscInt Check_Phase_above_below(PetscInt *phase_array, Marker *P,PetscInt num_phas)
{
	PetscInt n,it,size;
	size = num_phas;
	it=0;
	for(it=0;it<size;it++)
	{
		n=-1;
		if(P->phase==phase_array[it])
		{
			n=it;
			break;
		}
	}

	return n;
}
