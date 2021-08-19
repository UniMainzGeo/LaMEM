/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 **    filename:   dike.cpp
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
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Jana Schierjott
 **         Garrett Ito
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/*

    This file defines properties for the dike which is defined as an additional 
    source term on the RHS of the continutiy equation

*/
//---------------------------------------------------------------------------
//.................. DIKE PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "phase.h"
#include "parsing.h"
#include "JacRes.h"
#include "dike.h"
#include "constEq.h"
#include "bc.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBDikeCreate"
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, PetscBool PrintOutput)   
{

        // read all dike parameter blocks from file
  
        PetscInt jj;

        PetscErrorCode ierr;
        PetscFunctionBegin;

        //===============                                                                                                                                               
        // DIKE PARAMETER                                                                                                               
        //===============                                                                                                                                               

        // setup block access mode                                                                                                                                      
        ierr = FBFindBlocks(fb, _OPTIONAL_, "<DikeStart>", "<DikeEnd>"); CHKERRQ(ierr);

        if(fb->nblocks)
        {
                // print overview of dike blocks from file                                                                                                           
            if (PrintOutput)
            {
		      PetscPrintf(PETSC_COMM_WORLD,"Dike blocks : \n");
            }
                // initialize ID for consistency checks                                                                                                                 
            for(jj = 0; jj < _max_num_dike_ ; jj++) dbdike->matDike[jj].ID = -1;
		// error checking
                if(fb->nblocks > _max_num_dike_)
                {
                        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many dikes specified! Max allowed: %lld", (LLD)_max_num_dike_);
                }

                // store actual number of softening laws 
                dbdike->numDike = fb->nblocks;

                if (PrintOutput){
                        PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
                }
                // read each individual softening law                                                                                                                   
                for(jj = 0; jj < fb->nblocks; jj++)
                {
                    ierr = DBReadDike(dbdike, dbm, fb, PrintOutput); CHKERRQ(ierr);
                    fb->blockID++;
                }
        }

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBReadDike"
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, PetscBool PrintOutput)
{
        // read dike parameter from file 
        Dike     *dike;
        PetscInt  ID;
	
        PetscErrorCode ierr;
        PetscFunctionBegin;

        // Dike ID                                                                                                                                                         
        ierr    = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbdike->numDike-1); CHKERRQ(ierr);
        fb->ID  = ID;

        // get pointer to specified dike parameters
        dike = dbdike->matDike + ID;

        // check ID
        if(dike->ID != -1)
        {
                 SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate of Dike option!");
        }

        // set ID 
        dike->ID = ID;

	// read and store dike  parameters. 
        ierr = getScalarParam(fb, _REQUIRED_, "Mf", &dike->Mf,    1, 1.0); CHKERRQ(ierr);
        ierr = getScalarParam(fb, _REQUIRED_, "Mb", &dike->Mb, 1, 1.0); CHKERRQ(ierr);
	ierr = getIntParam(fb, _REQUIRED_, "PhaseID", &dike->PhaseID, 1, dbm->numPhases-1); CHKERRQ(ierr);  

        if (PrintOutput)
	    {
	    PetscPrintf(PETSC_COMM_WORLD,"   Dike parameters ID[%lld] : Mf = %g, Mb = %g\n", (LLD)(dike->ID), dike->Mf, dike->Mb);
	    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
        }

        PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Dike_k_heatsource"
PetscErrorCode Dike_k_heatsource(JacRes *jr,
                                Material_t *phases,
                                PetscScalar &Tc,
                                PetscScalar *phRat,          // phase ratios in the control volume
                                PetscScalar &k,
                                PetscScalar &rho_A)
{
        BCCtx       *bc;
        Dike        *dike;
        Ph_trans_t  *PhaseTrans;
        Material_t  *M;
        PetscInt     i, j, numDike;
        PetscScalar  v_spread, left, right, kfac, dikeRHS, tempdikeRHS;

        numDike    = jr->dbdike->numDike;// number of dikes

        bc         =  jr->bc;
        PhaseTrans =  jr->dbm->matPhtr;   // phase transition

        kfac = 0.0;
        dikeRHS = 0; 

	// alternatively: call GetDikeContr() here to not re-compute the dikeRHS
	
        // loop through all dikes
        for(j = 0; j < numDike; j++)
        {
            //access the material parameters of each dike block
            dike=jr->dbdike->matDike+j;

            // access the phase ID of the dike block
            i = dike->PhaseID;

             // check if the phase ratio of a dike phase is greater than 0 in the current cell
            if(phRat[i] > 0)
            {

                if(dike->Mb == dike->Mf)
                {
                    // constant M
                    v_spread = PetscAbs(bc->velin);
                    left = PhaseTrans->bounds[0];
                    right = PhaseTrans->bounds[1];
                    tempdikeRHS = dike->Mf * 2 * v_spread / PetscAbs(left-right);
                }

            //code for along-axis variation in M goes here

                else
                {
                    tempdikeRHS = 0.0;
                } 
                // end if (dike->Mb == dike-Mf)

                dikeRHS += phRat[i]*tempdikeRHS;
	
                M = &phases[i];

                //adjust k and heat source according to Behn & Ito [2005]
                if (Tc < M->T_liq && Tc > M->T_sol)
		  {
		    kfac  += phRat[i] / ( 1 + ( M->Latent_hx/ (M->Cp*(M->T_liq-M->T_sol))) );
		    rho_A += phRat[i]*(M->rho*M->Cp)*(M->T_liq-Tc)*dikeRHS;  // Cp not used in the paper,( M->rho*M->Cp) added to conserve units
		  }
		else if (Tc <= M->T_sol)
		  {
		    rho_A += phRat[i]*( M->rho*M->Cp)*( (M->T_liq-Tc) + M->Latent_hx/M->Cp )*dikeRHS;  // this term is here to conserve units: ( M->rho*M->Cp)
		    kfac += phRat[i];
		  }
		// end adjust k and heat source according to Behn & Ito [2005]
		
	    } //end check phaseRat>0
	    
        } //end for j=0 to numDike

        k=kfac*k;  //doesn't this need to be inside the dike-loop? , smarter to pass kfac instead of k?
               	// kfac is weighted average multiplier, k is already phase-dependent, hence already weighted by phase ratio from inside JAcResGetTempParam

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetDikeContr"
PetscErrorCode GetDikeContr(ConstEqCtx *ctx,
                                  PetscScalar *phRat,          // phase ratios in the control volume
                                  PetscScalar &dikeRHS)
{

        BCCtx       *bc;
        Dike        *dike;
        Ph_trans_t  *PhaseTrans;
        PetscInt     i, j, numDike;
        PetscScalar  v_spread, M, left, right, tempDikeRHS;
	//	PetscInt     k, sx, dsx, sy, dsy;
	//	PetscScalar  front, back;

        numDike    = ctx->numDike;
        bc         = ctx->bc;
	//	fs         = bc->fs;
        PhaseTrans = ctx->PhaseTrans;

	// loop through all dike blocks
        for(j = 0; j < numDike; j++)
	  {
            // access parameters of each dike block
            dike=ctx->matDike+j;
	    
            // access the correct phase ID of the dike parameters of each dike
            i = dike->PhaseID;
	    
	    // check if the phase ratio of a dike phase is greater than 0 in the current cell
            if(phRat[i]>0)
	      {
		if(dike->Mb == dike->Mf)
		  {
		    // constant M
		    M = dike->Mf;
		    v_spread = PetscAbs(bc->velin);
		    left = PhaseTrans->bounds[0];
		    right = PhaseTrans->bounds[1];
		    tempDikeRHS = M * 2 * v_spread / PetscAbs(left-right);
		  }
		
		/*else // Mb an Mf are different
		  {
		  
		  // access context
		  y = COORD_CELL(k,sy,fs->dsy);
		  
		  front = PhaseTrans->bounds[2];
		  back  = PhaseTrans->bounds[3];
		  if(front == back)
		  {
		  // linear interpolation between different M values, Mf is in front, Mb is in back
		  M = dike->Mf + (dike->Mb - dike->Mf) * (y/(PetscAbs(front+back))); 
		  tempDikeRHS = M * 2 * v_spread / PetscAbs(left+right);
		  }
		  else
		    {
                    // linear interpolation if the dike phase is oblique
                    M = dike->Mf + (dike->Mb - dike->Mf) * (y/(PetscAbs(front+back)));
                    tempDikeRHS = M * 2 * v_spread / PetscAbs(left+right);
		    }
		    }*/
		else
		 {
		   tempDikeRHS = 0.0;
		 }
	       
		dikeRHS += phRat[i]*tempDikeRHS;

	      } // close phase ratio loop
	    
	  } // close dike block loop
	
	PetscFunctionReturn(0);
}
