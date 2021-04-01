
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
 **		
 **		Main responsible persons for this routine:
 **			Andrea Piccolo
 **			Jianfeng Yang
 **			Boris Kaus
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/*	Bibliography reference for the phase transition
 * All the phase transition listed are coming from [1] (Tab.1).
 * [1] Manuele Faccenda, Luca Dal Zilio, The role of solid–solid phase transitions in mantle convection, Lithos,
        Volumes 268–271, 2017,
 * [2]B.R. Hacker, G.A. Abers, S.M. Peacock Subduction factory 1. Theoretical mineralogy, densities, seismic wave speeds, and H2O contents
        Journal of Geophysical Research, 108 (2003), 10.1029/2001JB001127
 * [3] Hydrous solidus of CMAS-pyrolite and melting of mantle plumes at the bottom of the upper mantle
		Geophysical Research Letters, 30 (2003), 10.1029/2003GL018318
 * [4] E.R. Hernandez, J. Brodholt, D. Alfè Structural, vibrational and thermodynamic properties of Mg2SiO4 and MgSiO3 minerals from first-principles simulations
	   Physics of the Earth and Planetary Interiors, 240 (2015), pp. 1-24
 * [5] The postspinel boundary in pyrolitic compositions determined in the laser-heated diamond anvil cell
	    Geophysical Research Letters, 41 (2014), pp. 3833-3841
 */
/*
	The routines in this file allow changing the phase of a marker depending on conditions
	set by the user.
	They thus allow adding phase transitions to a setup in a rahther simple manner.
	Moreover, a number of phase transitions have been predefined as profiles, such as
	the basalt-eclogite reaction.

	Routines that deal with phase transitions on particles. We really deal with
	two kinds of PT:
		1) change the phase number of a particle once a certain condition is reached
		   like depth of a particle of P/T conditions. The 'new' phase should
		   be defined accordingly in the input file

		2) We can take a thermodynamic phase diagram into account while computing 
			the density of a phase in the code. The same diagram can be used to
			compute melt fraction and the effective density of the partially molten
			material (taking into account that the actual melt content may be less
			than what the phase diagram indicates as a result of melt extraction)   
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
#include "surf.h"
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
    char            str_direction[_str_len_],Type_[_str_len_];

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
	ierr    =   getStringParam(fb, _REQUIRED_, "Type",Type_,NULL);  CHKERRQ(ierr);

	if(!strcmp(Type_,"Constant"))
	{
		ph->Type = _Constant_;
		ierr    =   Set_Constant_Phase_Transition(ph, dbm, fb);    	CHKERRQ(ierr);
	}
	else if(!strcmp(Type_,"Clapeyron"))
	{
		ph->Type = _Clapeyron_;
		ierr    =   Set_Clapeyron_Phase_Transition(ph, dbm, fb);   	CHKERRQ(ierr);
	}
	else if(!strcmp(Type_,"Box"))
	{
		ph->Type = _Box_;
		ierr    =   Set_Box_Phase_Transition(ph, dbm, fb);   		CHKERRQ(ierr);
	}
	
	ierr = getIntParam(fb,      _OPTIONAL_, "number_phases", &ph->number_phases,1 ,                     _max_num_tr_);      CHKERRQ(ierr);
	
	if ( ph->Type == _Box_ ){
		ierr = getIntParam(fb,      _OPTIONAL_, "PhaseOutside",     ph->PhaseOutside,	ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		ierr = getIntParam(fb, 	    _OPTIONAL_, "PhaseInside",    	ph->PhaseInside, 	ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
	}
	else{
		ierr = getIntParam(fb,      _OPTIONAL_, "PhaseBelow",       ph->PhaseBelow,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		ierr = getIntParam(fb, 	    _OPTIONAL_, "PhaseAbove",       ph->PhaseAbove,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		ierr = getScalarParam(fb,   _OPTIONAL_, "DensityBelow",     ph->DensityBelow,   ph->number_phases , 1.0);               CHKERRQ(ierr);
		ierr = getScalarParam(fb,   _OPTIONAL_, "DensityAbove",     ph->DensityAbove,   ph->number_phases,  1.0);               CHKERRQ(ierr);
	}
	

    ierr = getStringParam(fb, _OPTIONAL_, "PhaseDirection", 	str_direction, "BothWays"); 					            CHKERRQ(ierr);  
	if     	(!strcmp(str_direction, "BelowToAbove"))    ph->PhaseDirection  = 1;
	else if (!strcmp(str_direction, "OutsideToInside")) ph->PhaseDirection  = 1;
    else if (!strcmp(str_direction, "AboveToBelow"))    ph->PhaseDirection  = 2;
	else if (!strcmp(str_direction, "InsideToOutside")) ph->PhaseDirection  = 2;
	else if (!strcmp(str_direction, "BothWays"    ))    ph->PhaseDirection  = 0;
    else{      SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown Phase direction %s \n", str_direction);  }
	
	

	if (!ph->PhaseAbove || !ph->PhaseBelow)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type (Constant) (Clapeyron) ", (LLD)ID);
	}
    
	if (ph->Type == _Box_ ){
		
		if (ph->number_phases>0){
			PetscPrintf(PETSC_COMM_WORLD,"     Phase Outside      :   ");
			for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseOutside[i])); }
			PetscPrintf(PETSC_COMM_WORLD," \n");

			PetscPrintf(PETSC_COMM_WORLD,"     Phase Inside       :  ");
			for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseInside[i])); }
			PetscPrintf(PETSC_COMM_WORLD," \n");
			PetscPrintf(PETSC_COMM_WORLD,"     Direction          :   %s \n", str_direction);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD,"     No phase change    @   \n");
		}
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"     Phase Above        :  ");
		for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseAbove[i])); }
		PetscPrintf(PETSC_COMM_WORLD," \n");

		PetscPrintf(PETSC_COMM_WORLD,"     Phase Below        :  ");
		for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseBelow[i])); }
		PetscPrintf(PETSC_COMM_WORLD," \n");
		PetscPrintf(PETSC_COMM_WORLD,"     Direction          :   %s \n", str_direction);
	}
	
	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Constant_Phase_Transition"
PetscErrorCode  Set_Constant_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb)
{
	Scaling      *scal;
	char         Parameter[_str_len_];
	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm -> scal;

	ierr = getStringParam(fb, _REQUIRED_, "Parameter_transition",   Parameter, "none");  CHKERRQ(ierr);
	if(!strcmp(Parameter, "T"))
	{
		ph->Parameter_transition = _T_;
	}
	else if(!strcmp(Parameter, "P"))
	{
		ph->Parameter_transition = _Pressure_;
	}
	else if(!strcmp(Parameter, "Depth"))
	{
		ph->Parameter_transition = _Depth_;
	}
	else if(!strcmp(Parameter, "APS"))
	{
		ph->Parameter_transition = _PlasticStrain_;
	}
	else if(!strcmp(Parameter, "MeltFraction"))
	{
		ph->Parameter_transition = _MeltFraction_;
	}
	

	ierr = getScalarParam(fb, _REQUIRED_, "ConstantValue",          &ph->ConstantValue,        1,1.0);  CHKERRQ(ierr);


	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] :   Constant \n", (LLD)(ph->ID));
    PetscPrintf(PETSC_COMM_WORLD,"     Parameter          :   %s \n",    Parameter);
    PetscPrintf(PETSC_COMM_WORLD,"     Transition Value   :   %1.3f \n", ph->ConstantValue);

	if(ph->Parameter_transition==_T_)                   //  Temperature [Celcius]
	{
		ph->ConstantValue   =   (ph->ConstantValue + scal->Tshift)/scal->temperature;
	}
	else if(ph->Parameter_transition==_Pressure_)       //  Pressure [Pa]
	{
		ph->ConstantValue   /= scal->stress_si;
	}
	else if(ph->Parameter_transition==_Depth_)          //  Depth [km if geo units]
	{
		ph->ConstantValue   /= scal->length;
	}
	else if(ph->Parameter_transition==_PlasticStrain_)  //  accumulated plastic strain
	{
		ph->ConstantValue   = ph->ConstantValue;        // is already in nd units
	}
	else if(ph->Parameter_transition==_MeltFraction_)   //  melt fraction
	{
		ph->ConstantValue   = ph->ConstantValue;        // is already in nd units
	}
	else{
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for [Constant] Phase transition");
    }


	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Box_Phase_Transition"
PetscErrorCode  Set_Box_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb)
{
	Scaling      *scal;
	char         Parameter[_str_len_];
	PetscScalar  Box[6];
	PetscInt 	 i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm -> scal;

	ierr = getScalarParam(fb, _REQUIRED_, "PTBox_Bounds",   	ph->bounds,  		6, scal->length);    			CHKERRQ(ierr);

	ph->PhaseInside[0] = -1; 	// default

	for (i=0; i<6; i++){ Box[i] = ph->bounds[i]*scal->length; }		// dimensional units
	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] :   Box \n", (LLD)(ph->ID));
    PetscPrintf(PETSC_COMM_WORLD,"     Box Bounds         :   [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0],Box[1],Box[2],Box[3],Box[4],Box[5], scal->lbl_length);

	if (ph->PhaseInside[0] < 0) PetscPrintf(PETSC_COMM_WORLD,"     Don't set phase    @   \n");

	ierr = getStringParam(fb, _OPTIONAL_, "PTBox_TempType",   Parameter,  "none");    			CHKERRQ(ierr);
	if(!strcmp(Parameter, "none"))
	{
		ph->TempType = 0;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"     Don't set T inside @   \n");		CHKERRQ(ierr);
	}
	else if(!strcmp(Parameter, "constant"))
	{
		ph->TempType = 1;
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_cstTemp",   	&ph->cstTemp,  1, 1);    									CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Constant T inside  :   %1.1f %s \n", ph->cstTemp, scal->lbl_temperature);		CHKERRQ(ierr);
	 	ph->cstTemp = (ph->cstTemp + scal->Tshift)/scal->temperature;
	}
	else if(!strcmp(Parameter, "linear"))
	{
		ph->TempType = 2;
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_topTemp",   	&ph->topTemp,  1, 1);    									CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_botTemp",   	&ph->botTemp,  1, 1);    									CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Linear Temp; bot T :   %1.1f %s \n", ph->botTemp, scal->lbl_temperature);		CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Linear Temp; top T :   %1.1f %s \n", ph->topTemp, scal->lbl_temperature);		CHKERRQ(ierr);
		
		ph->topTemp = (ph->topTemp + scal->Tshift)/scal->temperature;
		ph->botTemp = (ph->botTemp + scal->Tshift)/scal->temperature;

	}
	else if(!strcmp(Parameter, "halfspace"))
	{
		ph->TempType = 3;
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_topTemp",   	&ph->topTemp,  		1, 1);    									CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_botTemp",   	&ph->botTemp,  		1, 1);    									CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_thermalAge", 	&ph->thermalAge,  	1, scal->time);    							CHKERRQ(ierr);
		
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Halfspace; top T   :   %1.1f %s \n", ph->topTemp, 				scal->lbl_temperature);		CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Halfspace; bot T   :   %1.1f %s \n", ph->botTemp, 				scal->lbl_temperature);		CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Halfspace; Age     :   %1.1f %s \n", ph->thermalAge*scal->time, 	scal->lbl_time);			CHKERRQ(ierr);
		
		ph->topTemp = (ph->topTemp + scal->Tshift)/scal->temperature;
		ph->botTemp = (ph->botTemp + scal->Tshift)/scal->temperature;

	}
	
	else{
		  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for PTBox_TempType [none; constant; linear; halfspace]");
	}
	
	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_Clapeyron_Phase_Transition"
PetscErrorCode  Set_Clapeyron_Phase_Transition(Ph_trans_t   *ph, DBMat *dbm, FB *fb)
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
	if(ph->neq>2 || ph->neq == 0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using any Clapeyron phase transition you cannot have a number of equation higher than 2, or equal to zero");

	}

	ierr = getScalarParam(fb, _OPTIONAL_, "clapeyron_slope",    ph->clapeyron_slope,ph->neq,    1.0); CHKERRQ(ierr);    // units??
	ierr = getScalarParam(fb, _OPTIONAL_, "P0_clapeyron",       ph->P0_clapeyron,   ph->neq,    1.0); CHKERRQ(ierr);    // units??
	ierr = getScalarParam(fb, _OPTIONAL_, "T0_clapeyron",       ph->T0_clapeyron,   ph->neq,    1.0); CHKERRQ(ierr);    // units??        

	if((!ph->clapeyron_slope || !ph->T0_clapeyron || !ph->clapeyron_slope   ||  !ph->Name_clapeyron))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "If you are using any Clapeyron phase transition avaiable you need to specify P0, T0, gamma and the number of equations ( P=(T-T0)*gamma +(P0) ).");
	}

    PetscPrintf(PETSC_COMM_WORLD,"       # Equations      :   %i    [ P = P0 + gamma*(T-T0) ] \n", ph->neq);
        
	for(it=0; it<ph->neq; it++)
	{
        PetscPrintf(PETSC_COMM_WORLD,"       eq[%i]            :   gamma = %- 4.2e [MPa/C], P0 = %4.2e [Pa],  T0 = %2.1f [deg C] \n", it, ph->clapeyron_slope[it], ph->P0_clapeyron[it],ph->T0_clapeyron[it]);

		ph->clapeyron_slope[it]     *=  1e6*(scal->temperature/scal->stress_si);                    // [K/MPa]
		ph->P0_clapeyron[it]        /=  (scal->stress_si);                                          // [Pa]
		ph->T0_clapeyron[it]        =   (ph->T0_clapeyron[it]+scal->Tshift)/ (scal->temperature);   // [Celcius]
	}
	PetscFunctionReturn(0);

}

/*  ierr = getStringParam(fb, _REQUIRED_, "Parameter_transition",   Parameter, "none");  CHKERRQ(ierr);
        if(!strcmp(Parameter, "T"))
        {
                ph->Parameter_transition = _T_;
		}*/    // WAS THIS HERE OR DID I ACCIDENTALLY COPY THIS HERE?

//------------------------------------------------------------------------------------------------------------//                                                                     
/* #undef __FUNCT__
#define __FUNCT__ "Check_AirPhaseRatio_Box_Transition()"
PetscErrorCode Check_AirPhaseRatio_Box_Transition(Material_t *mat, Controls ctrl, SolVarCell *svCell) //those ones only if new box feature: Ph_trans_t *ph, DBMat *dbm, FB *fb
{
        PetscInt     k, j, i, nx, ny, nz, sx, sy, sz, iter;
        PetscInt     AirPhase;
	PetscScalar  *phRat;
//     	  FDSTAG      *fs;  // likely not needed
        JacRes       *jr;
	Controls     ctrl;   // how to declare
	SolVarCell   *svCell;   // how to declare
	
        // access context                                                                                                                                                           
	AirPhase  = jr->surf->AirPhase;
	phRat     = ctx->phRat;
//        ctrl      = ctx->ctrl;

	// initialize
	AirPhase = 0;
	iter      = 0;
	
        PetscFunctionBegin;

	if(ctrl->actDike) {

	  // loop on all cells: access air phase ratio ( FROM FreeSurfGetAirPhaseRatio() )

	  START_STD_LOOP
	    {
	      // access phase ratio array                                                                                                                      
	      phRat = jr->svCell[iter++].phRat;
      
	      if(phRat[AirPhase] > 0.5) { 
		// turn dike off: 

		// get reference to material parameters table 
		//mat = &phases[i];

		mat->Mf = 0; 
		mat->Mb = 0;

	      }

	    }
	  END_STD_LOOP

	}

	} */

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
		//[1][2]
		ph->neq                 =   2;
		ph->P0_clapeyron[0]     =   2e9;        // Pa
		ph->T0_clapeyron[0]     =   800;        // C
		ph->clapeyron_slope[0]  =   1.5;
		
        ph->P0_clapeyron[1]     =   2e9;        // Pa    
		ph->T0_clapeyron[1]     =   700;
		ph->clapeyron_slope[1]  =   -30;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_WadsleyiteRingwoodite_wet"))
	{
		// [1][3]
		ph->neq                 =   1;
		ph->P0_clapeyron[0]     =   13.5e9;
		ph->T0_clapeyron[0]     =   1537;
		ph->clapeyron_slope[0]  =   5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_WadsleyiteRingwoodite_dry"))
	{
		//[1][4]
		ph->neq                 =   1;
		ph->P0_clapeyron[0]     =   18e9;
		ph->T0_clapeyron[0]     =   1597;
		ph->clapeyron_slope[0]  =   3.5;
	}
	else if(!strcmp(ph->Name_clapeyron,"Mantle_Transition_660km"))
	{
		//[1][5]
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
	PetscInt        i, ph,nPtr, numPhTrn,below,above,num_phas;
	PetscInt        PH1,PH2, ID;
	PetscScalar		T;
    PetscLogDouble  t;
	SolVarCell  	*svCell;
	Scaling      	*scal;


		
    // Retrieve parameters
	jr          =   actx->jr;
	dbm         =   jr->dbm;
	numPhTrn    =   dbm->numPhtr;
	scal 		=	dbm->scal;

	if (!numPhTrn) 	PetscFunctionReturn(0);		// only execute this function if we have phase transitions

    PrintStart(&t, "Phase_Transition", NULL);

	// loop over all local particles 		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
	nPtr        =   0;
	for(nPtr=0; nPtr<numPhTrn; nPtr++)
	{
		PhaseTrans = jr->dbm->matPhtr+nPtr;
				
		for(i = 0; i < actx->nummark; i++)      // loop over all (local) particles
		{
			// access marker
			P   =   &actx->markers[i];      

			// get consecutive index of the host cell of marker
			ID = 	actx->cellnum[i];

			// access host cell solution variables
			svCell = &jr->svCell[ID];
			
			num_phas    =   PhaseTrans->number_phases;
			if ( PhaseTrans->Type == _Box_ ){
			  below       =   Check_Phase_above_below(PhaseTrans->PhaseInside,   P, num_phas);
			  // phase inside check, if the phase doesn't belong there it returns n=below=-1 and goes to the next if condition loop
			  above       =   Check_Phase_above_below(PhaseTrans->PhaseOutside,  P, num_phas);   // phase outside check
			}
			else {
				below       =   Check_Phase_above_below(PhaseTrans->PhaseBelow,   P, num_phas);
				above       =   Check_Phase_above_below(PhaseTrans->PhaseAbove,   P, num_phas);
			}

			if  ( (below >= 0) || (above >= 0) )
			{
				PH2 = P->phase;  
				PH1 = P->phase;
                 // the current phase is indeed involved in a phase transition
				if      (   below>=0    )
				{
					if ( PhaseTrans->Type == _Box_ ){
					  PH1 = PhaseTrans->PhaseInside[below];    // PH1 is always inside the box
						PH2 = PhaseTrans->PhaseOutside[below];
					}
					else{
						PH1 = PhaseTrans->PhaseBelow[below];
						PH2 = PhaseTrans->PhaseAbove[below];
					}
				}
				else if (   above >=0   )
				{
					if ( PhaseTrans->Type == _Box_ ){
						PH1 = PhaseTrans->PhaseInside[above];    
						PH2 = PhaseTrans->PhaseOutside[above];
					}
					else{
						PH1 = PhaseTrans->PhaseBelow[above];
						PH2 = PhaseTrans->PhaseAbove[above];
					}
				}

				ph = P->phase;
				Transition(PhaseTrans, P, PH1, PH2, jr->ctrl, scal, svCell, &ph, &T, jr);

				if ( (PhaseTrans->Type == _Box_) ){
					if (PhaseTrans->PhaseInside[0]<0) ph = P->phase;				// do not change the phase
				}
			
                if (PhaseTrans->PhaseDirection==0){
                     P->phase 	=   ph;
                }
                else if ( (PhaseTrans->PhaseDirection==1) & (below>=0) ){
                    P->phase    =   ph;
                }
                else if ( (PhaseTrans->PhaseDirection==2) & (above>=0) ){
                    P->phase    =   ph;
                }
				P->T = T;	// set T

			}
		}

	}
	ierr = ADVInterpMarkToCell(actx);   CHKERRQ(ierr);

	
    PrintDone(t);

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------------------
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1, PetscInt PH2, Controls ctrl, Scaling *scal, SolVarCell *svCell, PetscInt *ph_out, PetscScalar *T_out, JacRes *jr)
{

        PetscInt    ph;
	PetscScalar T;
	//      JacRes *jr;	
	//	ConstEqCtx  *ctx;
	// access context   // NEW for the dike box option where air phase is not turned into dike phase
	//	ctrl      = ctx->ctrl;  when to use "Controls *ctrl" and when "Controls ctrl" ??
	// When to use ctrl->actDike and when ctrl.actDike?

	ph = P->phase;
	T  = P->T;

	if (PhaseTrans->Type==_Box_ && ctrl.actDike)
	{
	  Check_DikeBox_Phase_Transition(PhaseTrans,P,PH1,PH2, scal, &ph, &T, jr);            // compute phase & T within Box of dike, ignore airphase particles        
        }
	else if(PhaseTrans->Type==_Constant_)    // NOTE: string comparisons can be slow; we can change this to integers if needed
	{
		ph = Check_Constant_Phase_Transition(PhaseTrans,P,PH1,PH2, ctrl, svCell);
	}
	else if(PhaseTrans->Type==_Clapeyron_)
	{
	  ph = Check_Clapeyron_Phase_Transition(PhaseTrans,P,PH1,PH2, ctrl);
	}
	else if(PhaseTrans->Type==_Box_)                                      
	{
                Check_Box_Phase_Transition(PhaseTrans,P,PH1,PH2, scal, &ph, &T);                // compute phase & T within Box
	}
	
	*ph_out = ph;
	*T_out  = T;

	PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------------------------------------
    Sets the values for a phase transition that occurs @ a constant value
*/
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl, SolVarCell *svCell) 
{
    
    PetscInt 	ph;
	PetscScalar pShift;
	
	if (ctrl.pShift){
		pShift = ctrl.pShift;
	}
	else{
		pShift = 0.0;
	}

	ph = 0;
	if((PhaseTrans->Parameter_transition==_T_))   // NOTE: string comparisons can be slow; optimization possibility
		{
            // Temperature transition
            if ( P->T >= PhaseTrans->ConstantValue)     {   ph = PH2; }
			else                                        {   ph = PH1; }
		}

	if(PhaseTrans->Parameter_transition==_Pressure_)
		{
            if  ( (P->p+pShift) >= PhaseTrans->ConstantValue)   {   ph = PH2;   }
		    else                                       	 		{   ph = PH1;   }
          
		}

	if(PhaseTrans->Parameter_transition==_Depth_)
		{
          if ( P->X[2] >= PhaseTrans->ConstantValue)  {   ph = PH2;   }
          else                                        {   ph = PH1;   }
        }

	if(PhaseTrans->Parameter_transition==_PlasticStrain_) // accumulated plastic strain
		{
            if ( P->APS >= PhaseTrans->ConstantValue)  {   ph = PH2;        }
            else                                       {   ph = PH1;        }
        }

	if(PhaseTrans->Parameter_transition==_MeltFraction_) // melt fraction in cell
		{
            /* Different than the conditions above, the melt fraction is NOT stored
			 	on in the marker structure. Instead, it is stored on the grid. 
			 	We thus have to determine in which cell the current marker is & retrieve 
				the properties of that cell.
			*/ 
			PetscScalar mf;	// melt fraction of current cell

			mf = svCell->svBulk.mf;
			
			if (mf >= PhaseTrans->ConstantValue)  {   ph = PH2;        }
            else                                  {   ph = PH1;        }
        }

	return ph;
}

//------------------------------------------------------------------------------------------------------------//
PetscInt Check_Box_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,
			Scaling *scal, PetscInt *ph_out, PetscScalar *T_out)
{
	PetscInt 	ph;
	PetscScalar T;
	
	ph = P->phase;
	T  = P->T;
	if ( (P->X[0] >= PhaseTrans->bounds[0]) & (P->X[0] <= PhaseTrans->bounds[1]) &
		 (P->X[1] >= PhaseTrans->bounds[2]) & (P->X[1] <= PhaseTrans->bounds[3]) &
	     (P->X[2] >= PhaseTrans->bounds[4]) & (P->X[2] <= PhaseTrans->bounds[5])    ){  
	  
		// We are within the box
		ph = PH1;

		// Set the temperature structure
		if 		(PhaseTrans->TempType == 0){
			// do nothing
		}
		else if	(PhaseTrans->TempType == 1){
			// constant T inside
			T = PhaseTrans->cstTemp;
		}
		else if	(PhaseTrans->TempType == 2){
			// linear temperature profile
			PetscScalar zTop, zBot, topTemp, botTemp, d;

			zTop 	=	PhaseTrans->bounds[5];	// top 	  of domain 
			zBot 	=	PhaseTrans->bounds[4];	// bottom of domain 
			topTemp =	PhaseTrans->topTemp;	// T @ top
			botTemp =	PhaseTrans->botTemp;	// T @ bottom
			
			d 		=	(P->X[2]-zTop)/(zTop-zBot);		// normalized distance to top
			T 		=	d*(topTemp-botTemp) + topTemp;	// temperature profile	
		}
		else if	(PhaseTrans->TempType == 3){
			// halfspace cooling T inside
			PetscScalar zTop, topTemp, botTemp, T_age, d, kappa;

			zTop 	=	PhaseTrans->bounds[5];	// top 	  of domain 
			topTemp =	PhaseTrans->topTemp;	// T @ top
			botTemp =	PhaseTrans->botTemp;	// T @ bottom
			T_age 	=	PhaseTrans->thermalAge;	// thermal age

			kappa 	=	1e-6/( (scal->length_si)*(scal->length_si)/(scal->time_si));
			
			d 		=	zTop - P->X[2];
			T 		= 	(botTemp-topTemp)*erf(d/2.0/sqrt(kappa*T_age)) + topTemp;
		}
	}
	else{
		// Outside; keep T
		ph = PH2;
	}

	// return
	*ph_out = 	ph;
	*T_out 	=	T;

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------//                                                          
PetscInt Check_DikeBox_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, PetscInt *ph_out, PetscScalar *T_out, JacRes *jr)
{
        PetscInt     ph, AirPhase;
        PetscScalar  T;
	//	JacRes      *jr;   // how to initialize jr ?

	AirPhase=0.0;
	//	 PetscPrintf(PETSC_COMM_WORLD, "initialized airphase: %g \n", AirPhase);
	AirPhase  = jr->surf->AirPhase;
	//	PetscPrintf(PETSC_COMM_WORLD, "pointed airphase: %g \n", AirPhase);
        ph = P->phase;
        T  = P->T;
	
        if ( (P->X[0] >= PhaseTrans->bounds[0]) & (P->X[0] <= PhaseTrans->bounds[1]) &
                 (P->X[1] >= PhaseTrans->bounds[2]) & (P->X[1] <= PhaseTrans->bounds[3]) &
             (P->X[2] >= PhaseTrans->bounds[4]) & (P->X[2] <= PhaseTrans->bounds[5]) && ph != AirPhase    ){           // NEW condition: Dike is activated

                // We are within the box
	        ph = PH1;

                // Set the temperature structure                                                                 
                if(PhaseTrans->TempType == 0){
                        // do nothing                                                                                                                                    
                }
                else if (PhaseTrans->TempType == 1){
                        // constant T inside                                                                      
                        T = PhaseTrans->cstTemp;
                }
                else if (PhaseTrans->TempType == 2){
                        // linear temperature profile
                        PetscScalar zTop, zBot, topTemp, botTemp, d;
                        zTop    =       PhaseTrans->bounds[5];  // top    of domain                                                                                      
                        zBot    =       PhaseTrans->bounds[4];  // bottom of domain                                                                                      
                        topTemp =       PhaseTrans->topTemp;    // T @ top                                                                                               
                        botTemp =       PhaseTrans->botTemp;    // T @ bottom                                                                                           

                        d               =       (P->X[2]-zTop)/(zTop-zBot);             // normalized distance to top                                                    
                        T               =       d*(topTemp-botTemp) + topTemp;  // temperature profile                                                                   
                }
		else if (PhaseTrans->TempType == 3){
                        // halfspace cooling T inside                                                                                                                     
                        PetscScalar zTop, topTemp, botTemp, T_age, d, kappa;
                        zTop    =       PhaseTrans->bounds[5];  // top    of domain                                                                                      
                        topTemp =       PhaseTrans->topTemp;    // T @ top
			botTemp =       PhaseTrans->botTemp;    // T @ bottom                                                                                           
                        T_age   =       PhaseTrans->thermalAge; // thermal age                                                        

                        kappa   =       1e-6/( (scal->length_si)*(scal->length_si)/(scal->time_si));

                        d               =       zTop - P->X[2];
                        T               =       (botTemp-topTemp)*erf(d/2.0/sqrt(kappa*T_age)) + topTemp;
                }
        }
        else{
	  
	  // Outside; keep T
	  ph = PH2;
	}

	// return
	*ph_out =       ph;
	*T_out  =       T;

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------//
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, Controls ctrl)
{
	PetscInt ph,ip,neq;
	PetscScalar Pres[2], pShift;

  	if (ctrl.pShift){
		pShift = ctrl.pShift;
	}
	else{
		pShift = 0.0;
	}

	neq = PhaseTrans->neq;
	for (ip=0; ip<neq; ip++)
	{
		Pres[ip]    =   (P->T - PhaseTrans->T0_clapeyron[ip]) * PhaseTrans->clapeyron_slope[ip] + PhaseTrans->P0_clapeyron[ip];
	}
	if (neq==1)
	{
        if  ( (P->p+pShift) >= Pres[0]) {   ph  =   PH2;    }
        else                   			{   ph  =   PH1;    }
	}
	else
	{
        // in case we have two equations to describe the phase transition:
        if  ( ((P->p+pShift) >= Pres[0]) && ( (P->p+pShift) >= Pres[1]) )	{   ph  =   PH2;    }
        else                                                				{   ph  =   PH1;    }

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
