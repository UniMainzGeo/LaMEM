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
 * [6] Akaogi, M., Hashimoto, S., amp; Kojitani, H. (2018). Thermodynamic properties of ZrSiO 4 zircon and
 *  	reidite and of cotunnite-type ZrO 2 with application to high-pressure high-temperature phase relations
 *  	 in ZrSiO 4. Physics of the Earth and Planetary Interiors, 281,
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
#include "tssolve.h"
#include "dike.h"
//-----------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "DBMatReadPhaseTr"
PetscErrorCode DBMatReadPhaseTr(DBMat *dbm, FDSTAG *fs, FB *fb)
{
	// read phase transitions from file
	PetscFunctionBegin;

	Ph_trans_t      *ph;
	PetscInt        ID, i;
	Scaling         *scal;   // for the moving box
	PetscErrorCode  ierr;
	char            str_direction[_str_len_], Type_[_str_len_], Parameter[_str_len_];

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
		ierr    =   Set_Box_Phase_Transition(ph, dbm, fb);   	CHKERRQ(ierr);
	}
	else if(!strcmp(Type_,"NotInAirBox"))
	{
		ph->Type = _NotInAirBox_;
		ierr    =   Set_NotInAirBox_Phase_Transition(ph, dbm, fs, fb);		CHKERRQ(ierr);
	}
	
	ierr = getIntParam(fb,      _OPTIONAL_, "number_phases", &ph->number_phases,1 ,                     _max_tr_);      CHKERRQ(ierr);
	if ( ph->Type == _Box_ || ph->Type == _NotInAirBox_){
		ph->PhaseInside[0] = -1;	// default
		ierr = getIntParam(fb, 	    _OPTIONAL_, "PhaseInside",    	ph->PhaseInside, 	ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		
		ph->PhaseOutside[0] = -1;	// default
		ierr = getIntParam(fb,      _OPTIONAL_, "PhaseOutside",     ph->PhaseOutside,	ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		
	}
	else{
		ierr = getIntParam(fb,      _OPTIONAL_, "PhaseBelow",       ph->PhaseBelow,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		ierr = getIntParam(fb, 	    _OPTIONAL_, "PhaseAbove",       ph->PhaseAbove,     ph->number_phases , _max_num_phases_);  CHKERRQ(ierr);
		ierr = getScalarParam(fb,   _OPTIONAL_, "DensityBelow",     ph->DensityBelow,   ph->number_phases , 1.0);               CHKERRQ(ierr);
		ierr = getScalarParam(fb,   _OPTIONAL_, "DensityAbove",     ph->DensityAbove,   ph->number_phases,  1.0);               CHKERRQ(ierr);
	}
	
	// for moving NotInAirBox: read-in box-velocity
	if (ph->Type == _NotInAirBox_)
	  {
	    scal    =  dbm->scal;
	    ierr = getScalarParam(fb, _OPTIONAL_, "v_box", &ph->v_box, 1,  1.0); CHKERRQ(ierr);
	    ierr = getScalarParam(fb, _OPTIONAL_, "t0_box", &ph->t0_box, 1,  1.0); CHKERRQ(ierr);
	    ierr = getScalarParam(fb, _OPTIONAL_, "t1_box", &ph->t1_box, 1,  1.0); CHKERRQ(ierr);
	    ph->v_box  /= scal->velocity;
	    ph->t0_box /= scal->time;
	    ph->t1_box /= scal->time;
	  }
	
	ierr = getStringParam(fb, _OPTIONAL_, "PhaseDirection",     str_direction, "BothWays");                                          CHKERRQ(ierr);
	
	if     	(!strcmp(str_direction, "BelowToAbove"))    ph->PhaseDirection  = 1;
	else if (!strcmp(str_direction, "InsideToOutside")) ph->PhaseDirection  = 1; 
	else if (!strcmp(str_direction, "AboveToBelow"))    ph->PhaseDirection  = 2;
	else if (!strcmp(str_direction, "OutsideToInside")) ph->PhaseDirection  = 2;
	else if (!strcmp(str_direction, "BothWays"    ))    ph->PhaseDirection  = 0;
	else{      SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown Phase direction %s \n", str_direction);  }
	
	if (!ph->PhaseAbove || !ph->PhaseBelow)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You have not specify the correct phase transition type (Constant) (Clapeyron) ", (LLD)ID);
	}
    
	if (ph->Type == _Box_ || ph->Type == _NotInAirBox_){
		
		if (ph->number_phases>0){
			
			if (ph->PhaseOutside[0]>=0){
				PetscPrintf(PETSC_COMM_WORLD,"     Phase Outside      :   ");
				for (i=0; i<ph->number_phases; i++){    PetscPrintf(PETSC_COMM_WORLD," %d ", (LLD)(ph->PhaseOutside[i])); }
				PetscPrintf(PETSC_COMM_WORLD," \n");
			}
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
	
	// [Optional] parameter to reset on the tracers
	ierr = getStringParam(fb, _OPTIONAL_, "ResetParam",   Parameter,  "none");    			CHKERRQ(ierr);
	if(!strcmp(Parameter, "none"))
	{
		ph->Reset = 0;
	}
	else if(!strcmp(Parameter, "APS"))
	{
		ph->Reset = 1;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Reset Parameter    :   APS \n");		CHKERRQ(ierr);
 	}
	else{
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for PTBox_Reset [none; APS;]");
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
	else if(!strcmp(Parameter, "X"))
	{
		ph->Parameter_transition = _X_;
	}
	else if(!strcmp(Parameter, "Y"))
	{
		ph->Parameter_transition = _Y_;
	}
	else if(!strcmp(Parameter, "APS"))
	{
		ph->Parameter_transition = _PlasticStrain_;
	}
	else if(!strcmp(Parameter, "MeltFraction"))
	{
		ph->Parameter_transition = _MeltFraction_;
	}
	else if(!strcmp(Parameter, "t"))
	{
		ph->Parameter_transition = _Time_;
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
	else if(ph->Parameter_transition==_X_)          	//  X-coordinate [km if geo units]
	{
		ph->ConstantValue   /= scal->length;
	}
	else if(ph->Parameter_transition==_Y_)          	//  Y-coordinate [km if geo units]
	{
		ph->ConstantValue   /= scal->length;
	}
	else if(ph->Parameter_transition==_PlasticStrain_)  //  accumulated plastic strain
	{
		ph->ConstantValue   = ph->ConstantValue;        // 	is already in nd units
	}
	else if(ph->Parameter_transition==_MeltFraction_)   //  melt fraction
	{
		ph->ConstantValue   = ph->ConstantValue;        // is already in nd units
	}
	else if(ph->Parameter_transition==_Time_)       //  Time [s]
	{
		ph->ConstantValue   /= scal->time;
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

	// ph->PhaseInside[0] = -1; 	// default

	for (i=0; i<6; i++){ Box[i] = ph->bounds[i]*scal->length; }		// dimensional units
	PetscPrintf(PETSC_COMM_WORLD,"   Phase Transition [%lld] :   Box \n", (LLD)(ph->ID));
	PetscPrintf(PETSC_COMM_WORLD,"     Box Bounds         :   [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0],Box[1],Box[2],Box[3],Box[4],Box[5], scal->lbl_length);

	
	ierr = getIntParam(fb, _OPTIONAL_, "BoxVicinity",   &ph->BoxVicinity,  1, 1);

	if (ph->BoxVicinity==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Box Vicinity       :   Only check particles within vicinity of box (twice width) to determine inside/outside \n");	CHKERRQ(ierr);
	} 
	else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Box Vicinity       :   Use all particles to check inside/outside \n");	CHKERRQ(ierr);
	}
	
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
		  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for PTBox_TempType %s [none; constant; linear; halfspace]", Parameter);
	}

	PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "Set_NotInAirBox_Phase_Transition"
PetscErrorCode  Set_NotInAirBox_Phase_Transition(Ph_trans_t *ph, DBMat *dbm, FDSTAG *fs, FB *fb)
{
	Scaling      *scal;
	Discret1D	*dsy;
	char         Parameter[_str_len_];
	PetscInt 	 j,kk, found;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = dbm -> scal;
	dsy = &fs->dsy;

	ierr = getIntParam   (fb, _OPTIONAL_, "nsegs",  &ph->nsegs,  1,           _max_NotInAir_segs_);  CHKERRQ(ierr);
	
	//*
	if (ph->nsegs==0){
		ierr = getScalarParam(fb, _REQUIRED_, "PTBox_Bounds", ph->bounds, 6, scal->length);  CHKERRQ(ierr);
		ph->xbounds[0]=ph->bounds[0];
		ph->xbounds[1]=ph->bounds[1];
		ph->ybounds[0]=ph->bounds[2];
		ph->ybounds[1]=ph->bounds[3];
		ph->zbounds[0]=ph->bounds[4];
		ph->zbounds[1]=ph->bounds[5];
		ph->nsegs=1;		
	}
	else{
		ierr = getScalarParam(fb, _REQUIRED_, "xbounds", ph->xbounds, 2*ph->nsegs, scal->length);  CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "ybounds", ph->ybounds, 2*ph->nsegs, scal->length);  CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "zbounds", ph->zbounds, 2*ph->nsegs, scal->length);  CHKERRQ(ierr);
	}

	for (kk = 0; kk < ph->nsegs; kk++)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Phase Transition, NotInAirbox [%lld]: seg = %i, xbounds=[%g, %g], ybounds=[%g, %g], zbounds=[%g, %g] \n", \
			(LLD)(ph->ID), kk, \
		ph->xbounds[2*kk]* scal->length, ph->xbounds[2*kk+1]*scal->length,\
		ph->ybounds[2*kk]* scal->length, ph->ybounds[2*kk+1]*scal->length,\
		ph->zbounds[2*kk]* scal->length, ph->zbounds[2*kk+1]*scal->length);
	}



  	//create 1D array of xbound1 and xbound2, which define xbounds interpolated at each y-coord of cell
  	ierr = makeScalArray(&ph->cbuffL, 0, dsy->ncels+2); CHKERRQ(ierr);
  	ph->celly_xboundL = ph->cbuffL + 1;
  	ierr = makeScalArray(&ph->cbuffR, 0, dsy->ncels+2); CHKERRQ(ierr);
  	ph->celly_xboundR = ph->cbuffR + 1;


  	for(j = -1; j < dsy->ncels+1; j++)
  	{
  	   found=0;
	   for (kk = 0; kk < ph->nsegs; kk++)
	   {
  		if (dsy->ccoor[j] < ph->ybounds[0])
  		{
  			ph->celly_xboundL[j] = 1.0e12;
			ph->celly_xboundR[j] = -1.0e12;
			found=1;
			break;
		}
  		else if(ph->ybounds[2*kk] <= dsy->ccoor[j] && dsy->ccoor[j] <= ph->ybounds[2*kk+1])
  		{
			ph->celly_xboundL[j] = ph->xbounds[2*kk];
			ph->celly_xboundR[j] = ph->xbounds[2*kk+1];
			found=1;
			break;
		}
		else if (dsy->ccoor[j]>ph->ybounds[(ph->nsegs-1)*2+1])
		{
  			ph->celly_xboundL[j] = 1.0e12;
			ph->celly_xboundR[j] = -1.0e12;
			found=1;
			break;
		}
	   }

	   if (found==0) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_USER, " Cannot find NotInAirBox seg j=%i, dsy->ccoor=%g\n", \
	   		j, dsy->ccoor[j]*scal->length);
	   PetscPrintf(PETSC_COMM_WORLD, "DEBUGGING: NotInAirBox: %g, %g\n", ph->celly_xboundL[j], ph->celly_xboundR[j]); //DEBUGGING

	}


       ph->phtr_link_left = -1; 
	ierr = getIntParam(fb, _OPTIONAL_, "PhaseTransLinkLeft",   &ph->phtr_link_left,  1, dbm->numPhtr-1);
	if (ph->phtr_link_left>=0) {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"PhaseTransLinkLeft = %i\n", ph->phtr_link_left);	CHKERRQ(ierr);
	}

	ph->phtr_link_right = -1; 
	ierr = getIntParam(fb, _OPTIONAL_, "PhaseTransLinkRight",   &ph->phtr_link_right,  1, dbm->numPhtr-1);
	if (ph->phtr_link_right>=0) {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"PhaseTransLinkRight = %i\n", ph->phtr_link_right);	CHKERRQ(ierr);
	}


	
	ierr = getIntParam(fb, _OPTIONAL_, "BoxVicinity",   &ph->BoxVicinity,  1, 1);

	if (ph->BoxVicinity==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Box Vicinity       :   Only check particles within vicinity of box (twice width) to determine inside/outside \n");	CHKERRQ(ierr);
	} 
	else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"     Box Vicinity       :   Use all particles to check inside/outside \n");	CHKERRQ(ierr);
	}
	
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
		  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER, "Unknown parameter for PTBox_TempType %s [none; constant; linear; halfspace]", Parameter);
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
	else if(!strcmp(ph->Name_clapeyron,"Zircon_Reidite"))
	{
			//[6]
			ph->neq                 =   1;
			ph->P0_clapeyron[0]     =   8e9;
			ph->T0_clapeyron[0]     =   25;
			ph->clapeyron_slope[0]  =   1.9;
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
	TSSol           *ts;
	Ph_trans_t      *PhaseTrans;
	Marker          *P;
	JacRes          *jr;
	PetscInt        i, ph,nPtr, numPhTrn,below,above,num_phas;
	PetscInt        PH1,PH2, ID, InsideAbove,nphc; // nphc nophasechange condition
	PetscScalar     T, time, factor, dxBox, dyBox, dzBox;
	PetscLogDouble  t;
	SolVarCell      *svCell;
	Scaling         *scal;

    // Retrieve parameters
	jr          =   actx->jr;
	dbm         =   jr->dbm;
	ts          =   jr->ts;
	numPhTrn    =   dbm->numPhtr;
	scal 	    =	dbm->scal;
	time        =   jr->bc->ts->time;
	
	if (!numPhTrn) 	PetscFunctionReturn(0);		// only execute this function if we have phase transitions

	PrintStart(&t, "Phase_Transition", NULL);

	//For dynamic diking
	ierr = Locate_Dike_Zones(actx->jr); CHKERRQ(ierr);
	
	// loop over all phase transition laws		PetscPrintf(PETSC_COMM_WORLD,"PHASE = %d  i = %d, counter = %d\n",P->phase,i,counter);
	nPtr        =   0;
	for(nPtr=0; nPtr<numPhTrn; nPtr++)
	  {
	    PhaseTrans = jr->dbm->matPhtr+nPtr;

	    // Is the phase transition changing the phase, or other properites?
	    if((PhaseTrans->PhaseInside[0]>0 && PhaseTrans->PhaseOutside[0]>0) || (PhaseTrans->PhaseAbove[0]>0 && PhaseTrans->PhaseBelow[0]>0))
	    {
              nphc = 1;
	    }
	    else
	    {
              nphc = 0;
	    }
	    // calling the moving dike function

	    if ( PhaseTrans->Type == _NotInAirBox_ )
	    {
              if (PhaseTrans->v_box) {
                ierr = MovingBox(PhaseTrans, ts, jr); CHKERRQ(ierr);
              }

              ierr = LinkNotInAirBoxes(PhaseTrans, jr); CHKERRQ(ierr);

	    }
	    
		for(i = 0; i < actx->nummark; i++)      // loop over all (local) particles
		{
			// access marker
			P   =   &actx->markers[i];      

			// get consecutive index of the host cell of marker
			ID = 	actx->cellnum[i];

			// access host cell solution variables
			svCell = &jr->svCell[ID];
			
			num_phas    =   PhaseTrans->number_phases;

			if ( PhaseTrans->Type == _Box_ || PhaseTrans->Type == _NotInAirBox_ ){
				below       =   Check_Phase_above_below(PhaseTrans->PhaseInside,   P, num_phas);
				above       =   Check_Phase_above_below(PhaseTrans->PhaseOutside,  P, num_phas);
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
				if      (   (below>=0) && (nphc ==1))
				{
					if ( PhaseTrans->Type == _Box_ || PhaseTrans->Type == _NotInAirBox_){
						PH1 = PhaseTrans->PhaseInside[below];
						PH2 = PhaseTrans->PhaseOutside[below];
					}
					else{
						PH1 = PhaseTrans->PhaseBelow[below];
						PH2 = PhaseTrans->PhaseAbove[below];
					}
				}
				else if (   (above >=0) && (nphc==1))
				{
					if ( PhaseTrans->Type == _Box_ || PhaseTrans->Type == _NotInAirBox_){
						PH1 = PhaseTrans->PhaseInside[above];
						PH2 = PhaseTrans->PhaseOutside[above];
					}
					else{
						PH1 = PhaseTrans->PhaseBelow[above];
						PH2 = PhaseTrans->PhaseAbove[above];
					}
				}
				else
				{
					PH1 = P->phase;
					PH2 = PH1;
				}

				ph 			= P->phase;
				InsideAbove = 0;
				printf("DEBUGGING: Before Transition\n"); //DEBUGGING

				Transition(PhaseTrans, P, PH1, PH2, jr->ctrl, scal, svCell, &ph, &T, &InsideAbove, time, jr, ID);

				if ( (PhaseTrans->Type == _Box_ || PhaseTrans->Type == _NotInAirBox_ ) )
				{
					if (PhaseTrans->PhaseInside[0]<0){ 
						ph = P->phase;				// do not change the phase
					}

				
					if (PhaseTrans->BoxVicinity==1){
						factor = 1.0;
						dxBox  = (PhaseTrans->bounds[1]-PhaseTrans->bounds[0])*factor;
						dyBox  = (PhaseTrans->bounds[3]-PhaseTrans->bounds[2])*factor;
						dzBox  = (PhaseTrans->bounds[3]-PhaseTrans->bounds[2])*factor;
						
						if ( (P->X[0] < (PhaseTrans->bounds[0]-dxBox)) | (P->X[0] > (PhaseTrans->bounds[1]+dxBox)) |
							 (P->X[1] < (PhaseTrans->bounds[2]-dyBox)) | (P->X[1] > (PhaseTrans->bounds[3]+dyBox)) |
							 (P->X[2] < (PhaseTrans->bounds[4]-dzBox)) | (P->X[2] > (PhaseTrans->bounds[5]+dzBox))  )
						{
							ph = P->phase;				// do not change the phase
						}
					}

				
				}
				if (PhaseTrans->PhaseDirection==0){
					P->phase    =   ph;
				}
				else if ( (PhaseTrans->PhaseDirection==1) & (below>=0) ){
					P->phase    =   ph;
				}
				else if ( (PhaseTrans->PhaseDirection==2) & (above>=0) ){
					P->phase    =   ph;
				}
				P->T = T;	// set T

				// Reset other parameters on particles if requested
				if (PhaseTrans->PhaseDirection< 2){

					// Both ways or below2above	
					if (InsideAbove==1){
						if (PhaseTrans->Reset==1){
							P->APS = 0.0;
						}
					}
				}
				else{
					// Above to below
					if (InsideAbove==0){
						if (PhaseTrans->Reset==1){
							P->APS = 0.0;
						}
					}


				}
			}
			else{
				// allow cases in which we only reset T
				ph 			= P->phase;
				InsideAbove = 0;
				Transition(PhaseTrans, P, PH1, PH2, jr->ctrl, scal, svCell, &ph, &T, &InsideAbove, time, jr, ID);

				if ( (PhaseTrans->Type == _Box_ || PhaseTrans->Type == _NotInAirBox_ ) ){
					if (PhaseTrans->PhaseInside[0]<0){ 
						ph 		= P->phase;				// do not change the phase
					}
					
					if ((PhaseTrans->PhaseOutside[0]<0) & (PhaseTrans->PhaseDirection==2) & (InsideAbove==1)){ 	
						// PhaseOutside is set to -1 and OutsideToInside is selected, in which case we 
						// set everything inside the box to a constant phase (specified in PhaseInside)
						ph = PhaseTrans->PhaseInside[0];
						P->phase = ph;
					}
					
					P->T 	= T;	// set T
				}
			}


		}

	}
	ierr = ADVInterpMarkToCell(actx);   CHKERRQ(ierr);

    PrintDone(t);

	PetscFunctionReturn(0);
}

//----------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MovingBox"
PetscErrorCode MovingBox(Ph_trans_t *PhaseTrans, TSSol *ts, JacRes *jr)
{
  
  PetscScalar  t0_box, t1_box, v_box;
  PetscScalar  t_c, dt;
  PetscInt     j, ny;             
  FDSTAG 	*fs;  


  
  PetscFunctionBegin;
  
  dt  = ts->dt;       // time step
  t_c = ts->time;     // current time stamp, computed at the end of last time step round

  fs = jr->fs;
  ny = fs->dsy.ncels;
  
  // access the starting and end times of certain phase transition and the velocity of the phase transition-box
  t0_box = PhaseTrans->t0_box;
  t1_box = PhaseTrans->t1_box;
  v_box  = PhaseTrans->v_box;
  
  // check if the current time step is equal to the starting time of when the box is supposed to move
  if(t_c >= t0_box && t_c <= t1_box)
    {
      for(j = -1; j < ny+1; j++)
      {
         PhaseTrans->celly_xboundL[j] = PhaseTrans->celly_xboundL[j] + v_box * dt;
         PhaseTrans->celly_xboundR[j] = PhaseTrans->celly_xboundR[j] + v_box * dt;
      }
    }
  
  PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "LinkNotInAirBoxes"
PetscErrorCode LinkNotInAirBoxes(Ph_trans_t *PhaseTrans, JacRes *jr)
{
  
  Ph_trans_t   *PhaseTransLinkLeft;
  Ph_trans_t   *PhaseTransLinkRight;
  PetscScalar  Phase_Width;
  PetscInt     j, ny;             
  FDSTAG 	*fs;  


  
  PetscFunctionBegin;
  
  fs = jr->fs;
  ny = fs->dsy.ncels;
  
  if (PhaseTrans->phtr_link_left>=0)
  {
     PhaseTransLinkLeft = jr->dbm->matPhtr+PhaseTrans->phtr_link_left;
     for(j = -1; j < ny+1; j++)
      {
      	  Phase_Width = PhaseTrans->celly_xboundR[j]-PhaseTrans->celly_xboundL[j];
         PhaseTrans->celly_xboundL[j] = PhaseTransLinkLeft->celly_xboundR[j];
         PhaseTrans->celly_xboundR[j] = PhaseTrans->celly_xboundL[j]+Phase_Width;
      }
  }

  if (PhaseTrans->phtr_link_right>=0)
  {
     PhaseTransLinkRight = jr->dbm->matPhtr+PhaseTrans->phtr_link_right;
     for(j = -1; j < ny+1; j++)
      {
      	  Phase_Width = PhaseTrans->celly_xboundR[j]-PhaseTrans->celly_xboundL[j];
         PhaseTrans->celly_xboundR[j] = PhaseTransLinkRight->celly_xboundL[j];
         PhaseTrans->celly_xboundL[j] = PhaseTrans->celly_xboundR[j]-Phase_Width;
      }
  }

  PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------
PetscInt Transition(Ph_trans_t *PhaseTrans, Marker *P, PetscInt PH1, PetscInt PH2, Controls ctrl, Scaling *scal, 
		    SolVarCell *svCell, PetscInt *ph_out, PetscScalar *T_out, PetscInt *InsideAbove, PetscScalar time, JacRes *jr, PetscInt cellID)
{
	PetscInt    ph, InAbove;
	PetscScalar T;

	ph = P->phase;
	T  = P->T;
	InAbove = 0;
	
	if (PhaseTrans->Type==_NotInAirBox_ )
	{
		Check_NotInAirBox_Phase_Transition(PhaseTrans,P,PH1,PH2, scal, &ph, &T, jr, cellID);    // adjust phase according to T within Box but ignore airphase particles
	}
	else if(PhaseTrans->Type==_Constant_)    // NOTE: string comparisons can be slow; we can change this to integers if needed
	{
		Check_Constant_Phase_Transition(PhaseTrans,P,PH1,PH2, ctrl, svCell, &ph, &InAbove, time);
	}
	else if(PhaseTrans->Type==_Clapeyron_)
	{
		Check_Clapeyron_Phase_Transition(PhaseTrans,P,PH1,PH2, ctrl, &ph, &InAbove);
	}
	else if(PhaseTrans->Type==_Box_)
	{
		Check_Box_Phase_Transition(PhaseTrans,P,PH1,PH2, scal, &ph, &T, &InAbove);		// compute phase & T within Box
	}
	
	// Prepare output
	*ph_out 		= ph;
	*T_out  		= T;
	*InsideAbove 	= InAbove;

	PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------------------------------------
    Sets the values for a phase transition that occurs @ a constant value
*/
PetscInt Check_Constant_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, 
						Controls ctrl, SolVarCell *svCell, PetscInt *ph_out, PetscInt *InAbove, PetscScalar time) 
{
    
    PetscInt 	ph, InAb;
	PetscScalar pShift;
	
	if (ctrl.pShift){
		pShift = ctrl.pShift;
	}
	else{
		pShift = 0.0;
	}

	ph 		= 0;
	InAb 	= 0;
	if((PhaseTrans->Parameter_transition==_T_))   // NOTE: string comparisons can be slow; optimization possibility
		{
            // Temperature transition
            if ( P->T >= PhaseTrans->ConstantValue)     {   ph = PH2; InAb=1;  	}
			else                                        {   ph = PH1; 			}
		}

	if(PhaseTrans->Parameter_transition==_Pressure_)
		{
            if  ( (P->p+pShift) >= PhaseTrans->ConstantValue)   {   ph = PH2; InAb=1; 	}
		    else                                       	 		{   ph = PH1;  			}
          
		}
	if(PhaseTrans->Parameter_transition==_X_)
		{
          if ( P->X[0] >= PhaseTrans->ConstantValue)  {   ph = PH2; 	InAb=1;   }
          else                                        {   ph = PH1;   }
        }

	if(PhaseTrans->Parameter_transition==_Y_)
		{
          if ( P->X[1] >= PhaseTrans->ConstantValue)  {   ph = PH2; 	InAb=1;   }
          else                                        {   ph = PH1;   }
        }

	if(PhaseTrans->Parameter_transition==_Depth_)
		{
          if ( P->X[2] >= PhaseTrans->ConstantValue)  {   ph = PH2; InAb=1;	}
          else                                        {   ph = PH1; 		}
        }

	if(PhaseTrans->Parameter_transition==_PlasticStrain_) // accumulated plastic strain
		{
            if ( P->APS >= PhaseTrans->ConstantValue)  {   ph = PH2; InAb=1;}
            else                                       {   ph = PH1; 	}
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
			
			if (mf >= PhaseTrans->ConstantValue)  {   ph = PH2;   InAb=1;  }
            else                                  {   ph = PH1; 			}
        }
	if(PhaseTrans->Parameter_transition==_Time_)
		{
            if  ( time >= PhaseTrans->ConstantValue)   {   ph = PH2; InAb=1; 	}
		    else                                       	 		{   ph = PH1;  			}
          
		}

	// return
	*ph_out 	= ph;
	*InAbove 	= InAb;
	
	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------//
PetscInt Check_Box_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2,
			Scaling *scal, PetscInt *ph_out, PetscScalar *T_out, PetscInt *InAbove)
{
	PetscInt 	ph, InAb;
	PetscScalar T;
	
	ph = P->phase;
	T  = P->T;
	if ( (P->X[0] >= PhaseTrans->bounds[0]) & (P->X[0] <= PhaseTrans->bounds[1]) &
		 (P->X[1] >= PhaseTrans->bounds[2]) & (P->X[1] <= PhaseTrans->bounds[3]) &
	     (P->X[2] >= PhaseTrans->bounds[4]) & (P->X[2] <= PhaseTrans->bounds[5])  )
    {

        // We are within the box
		ph 		= PH1;
		InAb 	= 1;

		// Set the temperature structure
		if 		(PhaseTrans->TempType == 0){
		  //PetscPrintf(PETSC_COMM_WORLD,"inside box in loop in none temp \n");  // PRINT BOX TEST
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
		ph 		= PH2;
		InAb 	= 0;
	}

	// return
	*ph_out 	= 	ph;
	*T_out 		=	T;
	*InAbove 	= 	InAb;

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------//                                                          
PetscInt Check_NotInAirBox_Phase_Transition(Ph_trans_t *PhaseTrans, Marker *P,PetscInt PH1, PetscInt PH2, Scaling *scal, 
					PetscInt *ph_out, PetscScalar *T_out, JacRes *jr, PetscInt cellID)
{

	PetscInt     ph, AirPhase, I, J, K, nx, ny;             
	PetscScalar  T, xboundL, xboundR;
	FDSTAG 	*fs;
	Discret1D	*dsy;   
  
	PetscFunctionBegin;

	AirPhase  = jr->surf->AirPhase;
	ph = P->phase;
	T  = P->T;

	fs = jr->fs;
	dsy = &fs->dsy;  //dsy points to the address of jr->fs->dsy
	nx = fs->dsx.ncels;
	ny = fs->dsy.ncels;

	GET_CELL_IJK(cellID, I, J, K, nx, ny) //need to know J for celly_xboundL/R

	/*
       //particle backward of the cell center and adjacent cell is within phase trans box
       if (P->X[1] <= dsy->ccoor[J] && PhaseTrans->celly_xboundL[J-1] < PhaseTrans->celly_xboundR[J-1])  
	{
 	 	xboundL = PhaseTrans->celly_xboundL[J-1] + 
 	 	(PhaseTrans->celly_xboundL[J]-PhaseTrans->celly_xboundL[J-1])/(dsy->ccoor[J]-dsy->ccoor[J-1])*(P->X[1]-dsy->ccoor[J-1]);
 	 	xboundR = PhaseTrans->celly_xboundR[J-1] + 
 	 	(PhaseTrans->celly_xboundR[J]-PhaseTrans->celly_xboundR[J-1])/(dsy->ccoor[J]-dsy->ccoor[J-1])*(P->X[1]-dsy->ccoor[J-1]);
  	}
  	//particle is forward of the cell center and next cell is within phase trans box
  	else if (P->X[1] > dsy->ccoor[J] && PhaseTrans->celly_xboundL[J+1] < PhaseTrans->celly_xboundR[J+1]) 
  	{
  		xboundL = PhaseTrans->celly_xboundL[J] + 
  		(PhaseTrans->celly_xboundL[J+1]-PhaseTrans->celly_xboundL[J])/(dsy->ccoor[J+1]-dsy->ccoor[J])*(P->X[1]-dsy->ccoor[J]);
  		xboundR = PhaseTrans->celly_xboundR[J] + 
  		(PhaseTrans->celly_xboundR[J+1]-PhaseTrans->celly_xboundR[J])/(dsy->ccoor[J+1]-dsy->ccoor[J])*(P->X[1]-dsy->ccoor[J]);
  	}
  	//particle outside of phase transition box
  	else 
  	{
  		xboundL = PhaseTrans->celly_xboundL[J];
       	xboundR = PhaseTrans->celly_xboundR[J];
       }
       */


        	xboundL = PhaseTrans->xbounds[0];  //DEBUGGING
       	xboundR = PhaseTrans->xbounds[1]; //DEBUGGING


  	if 	( (xboundL <= P->X[0]) & (P->X[0] <= xboundR) &
       	(PhaseTrans->zbounds[0] <= P->X[2]) & (P->X[2] <= PhaseTrans->zbounds[1]) & (ph != AirPhase) )
    	{
        	// We are within the box
		ph = PH1;

		// Set the temperature structure                                                                 
		if      (PhaseTrans->TempType == 0){

		  //PetscPrintf(PETSC_COMM_WORLD,"inside NIAB in loop in none-temp \n");  // PRINT FOR NOTINAIRBOX TEST 

			// do nothing                                                                                                                                    
		}
		else if (PhaseTrans->TempType == 1){
			// constant T inside                                                                      
			T = PhaseTrans->cstTemp;
		}
		else if (PhaseTrans->TempType == 2){
			// linear temperature profile
			PetscScalar zTop, zBot, topTemp, botTemp, d;
			zTop    =       PhaseTrans->zbounds[1];  // top    of domain                                                                                      
			zBot    =       PhaseTrans->zbounds[0];  // bottom of domain                                                                                      
			topTemp =       PhaseTrans->topTemp;    // T @ top                                                                                               
			botTemp =       PhaseTrans->botTemp;    // T @ bottom                                                                                           

			d		=       (P->X[2]-zTop)/(zTop-zBot);             // normalized distance to top                                                    
			T		=       d*(topTemp-botTemp) + topTemp;  // temperature profile                                                                   
		}
		else if (PhaseTrans->TempType == 3){
			// halfspace cooling T inside                                                                                                                     
			PetscScalar zTop, topTemp, botTemp, T_age, d, kappa;
			zTop    =       PhaseTrans->zbounds[1];  // top    of domain                                                                                      
			topTemp =       PhaseTrans->topTemp;    // T @ top
			botTemp =       PhaseTrans->botTemp;    // T @ bottom                                                                                           
			T_age   =       PhaseTrans->thermalAge; // thermal age                                                        

			kappa   =       1e-6/( (scal->length_si)*(scal->length_si)/(scal->time_si));

			d		=       zTop - P->X[2];
			T		=       (botTemp-topTemp)*erf(d/2.0/sqrt(kappa*T_age)) + topTemp;
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
PetscInt Check_Clapeyron_Phase_Transition(Ph_trans_t *PhaseTrans,Marker *P,PetscInt PH1, PetscInt PH2, 
		Controls ctrl, PetscInt *ph_out, PetscInt *InAbove)
{
	PetscInt 		ph,ip,neq, InAb;
	PetscScalar 	Pres[2], pShift;

  	if (ctrl.pShift){
		pShift = ctrl.pShift;
	}
	else{
		pShift = 0.0;
	}

	neq 	= 	PhaseTrans->neq;
	InAb	=	0;
	for (ip=0; ip<neq; ip++)
	{
		Pres[ip]    =   (P->T - PhaseTrans->T0_clapeyron[ip]) * PhaseTrans->clapeyron_slope[ip] + PhaseTrans->P0_clapeyron[ip];
	}
	if (neq==1)
	{
        if  ( (P->p+pShift) >= Pres[0]) {   ph  =   PH2;  InAb=1;  }
        else                   			{   ph  =   PH1;    }
	}
	else
	{
        // in case we have two equations to describe the phase transition:
        if  ( ((P->p+pShift) >= Pres[0]) && ( (P->p+pShift) >= Pres[1]) )	{   ph  =   PH2;  InAb=1;  	}
        else                                                				{   ph  =   PH1;    		}

	}

	// return
	*ph_out 	= 	ph;
	*InAbove 	= 	InAb;

	
	PetscFunctionReturn(0);
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
//------------------------------------------------------------------------------------------------------------//
#undef __FUNCT__
#define __FUNCT__ "DynamicPhTrDestroy"
PetscErrorCode DynamicPhTrDestroy(DBMat *dbm)
{

	Ph_trans_t      *PhaseTrans;
	PetscInt   nPtr, numPhTrn;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	printf("DEBUGGING:  DynamicPhTrDestroy 1 \n");

	numPhTrn    =   dbm->numPhtr;
	nPtr        =   0;

	for(nPtr=0; nPtr<numPhTrn; nPtr++)
	{
	  printf("DEBUGGING:  DynamicPhTrDestroy 2 \n");

	  PhaseTrans = dbm->matPhtr+nPtr;
	  printf("DEBUGGING:  DynamicPhTrDestroy 3 \n");

	  ierr = PetscFree(PhaseTrans->cbuffL);        CHKERRQ(ierr);
	  ierr = PetscFree(PhaseTrans->cbuffR);       CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}