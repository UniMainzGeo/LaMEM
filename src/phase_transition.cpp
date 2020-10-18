
/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2020, JGU Mainz, Anton Popov, Boris Kaus
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
 **			Andrea Piccolo 
 ** 		Jianfeng Yang
 **         Tobias Baumann
 **         Adina Pusok
 **		
 **		Main responsible persons for this routine:
 **			Andrea Piccolo
 **			Jianfeng Yang
 **			Boris Kaus
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

/*
	The routines in this file allow changing the phase of a marker depending on conditions
	set by the user.
	They thus allow adding phase transitions to a setup in a rahther simple manner.
	Moreover, a number of phase transitions have been predefined as profiles, such as
	the basalt-eclogite reaction.
*/

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
	// read phase transitions from file
	PetscFunctionBegin;

	Ph_trans_t      *ph;
	PetscInt        ID, i;
	PetscErrorCode  ierr;
    char            str_direction[_str_len_];

	// Phase transition law ID
	ierr    =   getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbm->numPhtr-1); CHKERRQ(ierr);

	// get pointer to specified softening law
	ph      =   dbm->matPhtr + ID;

	// check ID
	if(ph->ID != -1)
	{
		 SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate phase transition law!");
	}

	// set ID
	ph->ID  =   ID;

	ierr    =   getStringParam(fb, _REQUIRED_, "Type", ph->Type,0);  CHKERRQ(ierr);
	if (!ph->Type)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type [Constant; Clapeyron; Box_type] ", (LLD)ID);
	}

	if(!strcmp(ph->Type,"Constant"))
	{
		ierr    =   Set_Constant_Phase_Transition(ph, dbm, fb,ID);    CHKERRQ(ierr);
	}
	if(!strcmp(ph->Type,"Clapeyron"))
	{
		ierr    =   Set_Clapeyron_Phase_Transition(ph, dbm, fb,ID);   CHKERRQ(ierr);
	}

	if(!strcmp(ph->Type,"Box_type"))
	{
		ierr    =   Set_Box_Within_Transition(ph, dbm, fb,ID);        CHKERRQ(ierr);
	}

	ierr = getIntParam(fb,      _OPTIONAL_, "number_phases", &ph->number_phases,1 ,                     _max_num_tr_);      CHKERRQ(ierr);
	ierr = getIntParam(fb,      _OPTIONAL_, "PhaseBelow",       ph->PhaseBelow,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
	ierr = getIntParam(fb, 	    _OPTIONAL_, "PhaseAbove",       ph->PhaseAbove,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
	ierr = getScalarParam(fb,   _OPTIONAL_, "DensityBelow",     ph->DensityBelow,   ph->number_phases , 1.0);               CHKERRQ(ierr);
	ierr = getScalarParam(fb,   _OPTIONAL_, "DensityAbove",     ph->DensityAbove,   ph->number_phases,  1.0);               CHKERRQ(ierr);

    ierr = getStringParam(fb, _OPTIONAL_, "PhaseDirection", 	str_direction, "BothWays"); 					            CHKERRQ(ierr);  
	if     	(!strcmp(str_direction, "BelowToAbove"))    ph->PhaseDirection  = 1;
	else if (!strcmp(str_direction, "AboveToBelow"))    ph->PhaseDirection  = 2;
    else if (!strcmp(str_direction, "BothWays"    ))    ph->PhaseDirection  = 0;
    else{      SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown Phase direction %s \n", str_direction);  }
        

	if (!ph->PhaseAbove || !ph->PhaseBelow)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type (Constant) (Clapeyron) ", (LLD)ID);
	}
    
    PetscPrintf(PETSC_COMM_WORLD,"     Phase Above        :  ");
    for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseAbove[i])); }
    PetscPrintf(PETSC_COMM_WORLD," \n");

    PetscPrintf(PETSC_COMM_WORLD,"     Phase Below        :  ");
    for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseBelow[i])); }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    PetscPrintf(PETSC_COMM_WORLD,"     Direction          :   %s \n", str_direction);
    
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

	ierr = getStringParam(fb, _REQUIRED_, "Parameter_transition",   ph->Parameter_transition, "none");  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "ConstantValue",          &ph->ConstantValue,        1,1.0);  CHKERRQ(ierr);
	if((!ph->Parameter_transition || !ph->ConstantValue))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using a [constant] phase transition you need to specify the parameter: P = Pressure, T = temperature, APS = PlasticStrain", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] :   Constant \n", (LLD)(ph->ID));
    PetscPrintf(PETSC_COMM_WORLD,"     Parameter          :   %s \n",    ph->Parameter_transition);
    PetscPrintf(PETSC_COMM_WORLD,"     Transition Value   :   %1.3f \n", ph->ConstantValue);

	if(!strcmp(ph->Parameter_transition,"T"))       //  Temperature [Celcius]
	{
		ph->ConstantValue   =   (ph->ConstantValue + scal->Tshift)/scal->temperature;
	}
	else if(!strcmp(ph->Parameter_transition,"P"))       //  Pressure [Pa]
	{
		ph->ConstantValue   /= scal->stress_si;
	}
	else if(!strcmp(ph->Parameter_transition,"Depth"))   //  Depth [km if geo units]
	{
		ph->ConstantValue   /= scal->length;
	}
	else if(!strcmp(ph->Parameter_transition,"APS"))   //  accumulated plastic strain
	{
		ph->ConstantValue   = ph->ConstantValue;    // is already in nd units
	}
    else{
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for [Constant] Phase transition");
    }


	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Constant_Phase_Transition"
PetscErrorCode  Set_Clapeyron_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb, PetscInt ID)
{
	PetscFunctionBegin;

	PetscErrorCode  ierr;
	Scaling         *scal;
	PetscInt        it=0;

	scal    = dbm -> scal;
	ierr    = getStringParam(fb, _OPTIONAL_, "Name_Clapeyron", ph->Name_clapeyron, "none");  CHKERRQ(ierr);
	if (ph->Name_clapeyron)
	{
		ierr = SetClapeyron_Eq(ph); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] :   Clapeyron \n", (LLD)(ph->ID));
        PetscPrintf(PETSC_COMM_WORLD,"     Transition law     :   %s\n", ph->Name_clapeyron);
    }

	ierr = getIntParam   (fb, _OPTIONAL_, "numberofequation",   &ph->neq,           1,          2.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "clapeyron_slope",    ph->clapeyron_slope,ph->neq,    1.0); CHKERRQ(ierr);    // units??
	ierr = getScalarParam(fb, _OPTIONAL_, "P0_clapeyron",       ph->P0_clapeyron,   ph->neq,    1.0); CHKERRQ(ierr);    // units??
	ierr = getScalarParam(fb, _OPTIONAL_, "T0_clapeyron",       ph->T0_clapeyron,   ph->neq,    1.0); CHKERRQ(ierr);    // units??        

	if((!ph->clapeyron_slope || !ph->T0_clapeyron || !ph->clapeyron_slope   ||  !ph->Name_clapeyron))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using any Clapeyron phase transition avaiable you need to specify P0, T0, gamma and the number of equations ( P=(T-T0)*gamma +(P0) ).", (LLD)ID);
	}

    PetscPrintf(PETSC_COMM_WORLD,"       # Equations      :   %i    [ P = P0 + gamma*(T-T0) ] \n", ph->neq);
        
	for(it=0; it<ph->neq; it++)
	{
        PetscPrintf(PETSC_COMM_WORLD,"       eq[%i]            :   gamma = %- 4.2e [MPa/C], P0 = %4.2e [Pa],  T0 = %2.1f [deg C] \n", it, ph->clapeyron_slope[it], ph->P0_clapeyron[it],ph->T0_clapeyron[it]);

		ph->clapeyron_slope[it]     *=  1e6*(scal->temperature/scal->stress_si);                    // [C/MPa]??
		ph->P0_clapeyron[it]        /=  (scal->stress_si);                                          // [Pa]
		ph->T0_clapeyron[it]        =   (ph->T0_clapeyron[it]+scal->Tshift)/ (scal->temperature);   // [Celcius]
	}
	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Box_Within_Transition"
PetscErrorCode  Set_Box_Within_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb,PetscInt ID)
{
	
	PetscFunctionBegin;

	PetscErrorCode  ierr;
    Scaling         *scal;
	PetscInt        it;

	scal = dbm -> scal;

	ierr = getScalarParam(fb, _REQUIRED_, "Geometrical_box", ph->Geometric_box,6, 1.0); CHKERRQ(ierr);

	for(it=0; it<6; it++)
	{
		ph->Geometric_box[it]   /=  scal->length;
	}

	ierr = getScalarParam(fb, _REQUIRED_, "DeltaT_within", &ph->dT_within,1, 1.0);      CHKERRQ(ierr);


	ph->dT_within   =   (ph->dT_within+ scal->Tshift)+scal->temperature;

	PetscFunctionReturn(0);

}
// ---------------------------------------------------------------------------------------------------------- //
#undef __FUNCT__
#define __FUNCT__ "Overwrite_Density"
PetscErrorCode  Overwrite_density(DBMat *dbm)
{
    /* 
        This changes the density values of the phases involved with a pre-defined phase transition such 
        that they are compatible (e.g., tthis ensures that we do not use bogus values for the basalt->eclogite transition,
        for example) 
    */

	PetscFunctionBegin;
	Scaling         *scal;
	Ph_trans_t      *ph;
	Material_t      *mat;
	PetscInt        numPhTrn,   nPtr,       num_phases, iter, jj1, jj2;
	PetscScalar     rho_above,  rho_below,  rho_scal;
    
    // Retrieve parameters:
    scal        = dbm->scal;
	rho_scal    = scal->density;
	mat         = dbm->phases;
	numPhTrn    = dbm->numPhtr;

    PetscPrintf(PETSC_COMM_WORLD,"\n   Adjusting density values due to phase transitions: \n");
    for(nPtr=0; nPtr<numPhTrn;  nPtr++)
	{
		ph          = dbm->matPhtr+nPtr;
		num_phases  = ph->number_phases;

		for(iter=0; iter<num_phases;    iter++)
		{
			rho_above = ph->DensityAbove[iter];
			rho_below = ph->DensityBelow[iter];
			if(rho_above > 0 && rho_below >0)
			{
				jj1             =   ph->PhaseBelow[iter];
				mat[jj1].rho    =   rho_below/rho_scal;
				PetscPrintf(PETSC_COMM_WORLD,"     Phase              : %4d, rho = %4.1f %s \n",jj1,rho_below, scal->lbl_density);

				jj2             =   ph->PhaseAbove[iter];
				mat[jj2].rho    =   rho_above/rho_scal;
				PetscPrintf(PETSC_COMM_WORLD,"     Phase              : %4d, rho = %4.1f %s \n",jj2,rho_above, scal->lbl_density);
			}

		}
    }
	PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "SetClapeyron_Eq"
PetscErrorCode SetClapeyron_Eq(Ph_trans_t *ph)
{
	/*
        This has predefined phase transitions that are specified as linear lines, where
        the phase transition is given by
        
        P = P0 + clapeyron_slope*(T-T0)

        Where 
            P0              :   Zero interection [Pa]
            clapeyron_slope :   Slope of the line [MPa/C]??
            T0              :   Shift in temperature [Celcius]
    */

	PetscFunctionBegin;
	if (!strcmp(ph->Name_clapeyron,"Eclogite"))
	{
		// Source Faccenda and Dal Zilio et al 2017 [Hacker et al. (2003) Morb+H2O]
		ph->neq                 =   2;
		ph->P0_clapeyron[0]     =   2e9;        // Pa
		ph->T0_clapeyron[0]     =   800;        // C
		ph->clapeyron_slope[0]  =   1.5;
		
        ph->P0_clapeyron[1]     =   2e9;        // Pa    
		ph->T0_clapeyron[1]     =   700;
		ph->clapeyron_slope[1]  =   -30;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_410km"))
	{
		// Source: Faccenda and Dal Zilio et al 2017 [Olivine-Wadseylite Phase transition anhydrous Pyrolite+H2O Litasov and Ohtani (2003)]
		ph->neq                 =   1;
		ph->P0_clapeyron[0]     =   13.5e9;
		ph->T0_clapeyron[0]     =   1537;
		ph->clapeyron_slope[0]  =   5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_510km"))
	{
		// Source: Faccenda and Dal Zilio et al 2017 [WadsleyiteRingwooditePhase transition anhydrous Fo100(theoretical) Hernandez et al 2015]
		ph->neq                 =   1;
		ph->P0_clapeyron[0]     =   18e9;
		ph->T0_clapeyron[0]     =   1597;
		ph->clapeyron_slope[0]  =   3.5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_660km"))
	{
		// Source:  Faccenda and Dal Zilio et al 2017 [Post Spinel Phase transition anhydrous Pyrolite Ye et al 2014] 
		ph->neq                 =   1;
		ph->P0_clapeyron[0]     =   23e9;
		ph->T0_clapeyron[0]     =   1667;
		ph->clapeyron_slope[0]  =   -2.5;
	}

	PetscFunctionReturn(0);
}
//===========================================================================================================//
#undef __FUNCT__
#define __FUNCT__ "Phase_Transition"
PetscErrorCode Phase_Transition(AdvCtx *actx)
{
	// creates arrays to optimize marker-cell interaction
	PetscFunctionBegin;

    PetscErrorCode  ierr;
    DBMat           *dbm;
	Ph_trans_t      *PhaseTrans;
	Marker          *P;
	JacRes          *jr;
	PetscInt        i, ph,nPtr, numPhTrn,below,above,num_phas,outside;
	PetscInt        PH1,PH2;

    // Retrieve parameters
	jr          =   actx->jr;
	dbm         =   jr->dbm;
	numPhTrn    =   dbm->numPhtr;

	// loop over all local particles 		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
	nPtr        =   0;
	for(nPtr=0; nPtr<numPhTrn; nPtr++)
	{
		PhaseTrans = jr->dbm->matPhtr+nPtr;
		for(i = 0; i < actx->nummark; i++)      // loop over all (local) particles
		{

			P   =   &actx->markers[i];      // retrieve marker

			if (!strcmp(PhaseTrans->Type,"Box_type")) // NOTE: doing a string comparison @ every particle can be time-consuming; if needed this can be changed into an integer check
			{

				num_phas    =   PhaseTrans->number_phases;

				outside     =   Check_Phase_above_below(PhaseTrans->PhaseBelow,P,num_phas);
				if(outside>=0)
				{
					PH1         =   PhaseTrans->PhaseAbove[outside];
					ph          =   Transition(PhaseTrans, P, PH1, PH2,nPtr);
					P->phase    =   ph;
				}
			}
			else
			{
				num_phas    =   PhaseTrans->number_phases;
				below       =   Check_Phase_above_below(PhaseTrans->PhaseBelow, P, num_phas);
				above       =   Check_Phase_above_below(PhaseTrans->PhaseAbove, P, num_phas);

				if  ( (below >= 0) || (above >= 0) )
				{
                    // the current phase is indeed involved in a phase transition    
					if      (   below>=0    )
					{
						PH1 = PhaseTrans->PhaseBelow[below];
						PH2 = PhaseTrans->PhaseAbove[below];
					}
					else if (   above >=0   )
					{
						PH1 = PhaseTrans->PhaseBelow[above];
						PH2 = PhaseTrans->PhaseAbove[above];
					}
					ph          =   Transition(PhaseTrans, P, PH1, PH2,nPtr);

                    if (PhaseTrans->PhaseDirection==0){
                        P->phase    =   ph;
                    }
                    else if (PhaseTrans->PhaseDirection==1 & below>=0 ){
                        P->phase    =   ph;
                    }
                    else if (PhaseTrans->PhaseDirection==2 & above>=0 ){
                        P->phase    =   ph;
                    }
                    
    			}
			}
		}

        // Interpolate markers to grid (BORIS: does this have to be done after every PhaseTransition of can it be done @ the end of the routine?):
		ierr = ADVInterpMarkToCell(actx);   CHKERRQ(ierr);

	}

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------------------
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1,PetscInt PH2,PetscInt ID)
{
	PetscInt ph;

	if(!strcmp(PhaseTrans->Type,"Constant"))    // NOTE: string comparisons can be slow; we can change this to integers if needed
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

/*------------------------------------------------------------------------------------------------------------
    Sets the values for a phase transition that occurs @ a constant value
*/
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,PetscInt ID) 
{
    
    PetscInt ph;


	if(!strcmp(PhaseTrans->Parameter_transition,"T"))   // NOTE: string comparisons can be slow; optimization possibility
		{
            // Temperature transition
            if ( P->T >= PhaseTrans->ConstantValue)     {   ph = PH2; }
			else                                        {   ph = PH1; }
		}

	if(!strcmp(PhaseTrans->Parameter_transition,"P"))
		{
            if  ( P->p >= PhaseTrans->ConstantValue)    {   ph = PH2;   }
		    else                                        {   ph = PH1;   }
          
		}

	if(!strcmp(PhaseTrans->Parameter_transition,"Depth"))
		{
          if ( P->X[2] >= PhaseTrans->ConstantValue)  {   ph = PH2;   }
          else                                        {   ph = PH1;   }
        }

	if(!strcmp(PhaseTrans->Parameter_transition,"APS")) // accumulated plastic strain
		{
            if ( P->APS >= PhaseTrans->ConstantValue)  {   ph = PH2;        }
            else                                       {   ph = PH1;        }
        }

	return ph;
}

//------------------------------------------------------------------------------------------------------------//
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,PetscInt ID) 
{
	PetscInt ph,ip,neq;
	PetscScalar Pres[2];


	neq = PhaseTrans->neq;
	for (ip=0; ip<neq; ip++)
	{
		Pres[ip]    =   (P->T - PhaseTrans->T0_clapeyron[ip]) * PhaseTrans->clapeyron_slope[ip] + PhaseTrans->P0_clapeyron[ip];
	}
	if (neq==1)
	{   
        if  ( P->p >= Pres[0]) {   ph  =   PH2;    }
        else                   {   ph  =   PH1;    }
	}
	else
	{
        // in case we have two equations to describe the phase transition:
        if  ( (P->p >= Pres[0]) && (P->p >= Pres[1]) )      {   ph  =   PH2;    }
        else                                                {   ph  =   PH1;    }
      
	}

	return ph;
}
//------------------------------------------------------------------------------------------------------------//

PetscInt Check_Constant_Box_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1,PetscInt ID) 
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

//------------------------------------------------------------------------------------------------------------//
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
