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
#include "scaling.h"
#include "objFunct.h"
#include "JacRes.h"
#include "phase_transition.h"
#include "dike.h"
#include "constEq.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBDikeCreate"
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, FB *fb, PetscBool PrintOutput)   
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
                for(jj = 0; jj < _max_num_dike_; jj++) dbdike->matDike[jj].ID = -1;

                // store actual number of softening laws 
                dbdike->numDike = fb->nblocks;

                if (PrintOutput){
                        PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
                }
                // read each individual softening law                                                                                                                   
                for(jj = 0; jj < fb->nblocks; jj++)
                {
                        ierr = DBReadDike(dbdike, fb, PrintOutput); CHKERRQ(ierr);

                        fb->blockID++;
                }
        }

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBReadDike"
PetscErrorCode DBReadDike(DBPropDike *dbdike, FB *fb, PetscBool PrintOutput)
{
        // read softening law from file 
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
	ierr = getIntParam(fb, _REQUIRED_, "Phase", &dike->Phase, 1, 1); CHKERRQ(ierr);  

	
        if (PrintOutput)
	  {
	    PetscPrintf(PETSC_COMM_WORLD,"   Dike parameters ID[%lld] : Mf = %g, Mb = %g\n", (LLD)(dike->ID), dike->Mf, dike->Mb);
	    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
        }

        PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetDikeContr"
PetscErrorCode GetDikeContr(ConstEqCtx  *ctx,
                                  PetscScalar *phRat,          // phase ratios in the control volume
                                  PetscScalar &dikeRHS)
{

  BCCtx       *bc;
        Dike        *dike;
        Ph_trans_t  *PhaseTrans;
        PetscInt     i, j, numDike;
        PetscScalar  v_spread, M, left, right;

	
	numDike    = ctx->dbdike->numDike;

	PetscPrintf(PETSC_COMM_WORLD,"numdike \n");
	
        dike       = ctx->dike;

	PetscPrintf(PETSC_COMM_WORLD,"begin \n");
        bc         = ctx->bc;
	PetscPrintf(PETSC_COMM_WORLD,"begin \n");
        PhaseTrans = ctx->PhaseTrans;

	PetscPrintf(PETSC_COMM_WORLD,"begin \n");

	
          // loop through all dikes
          for(j = 0; j < numDike; j++)
            {
	      PetscPrintf(PETSC_COMM_WORLD,"2a \n");
	      // access the phase ID of the dike parameters
              i = dike->Phase;
	      PetscPrintf(PETSC_COMM_WORLD,"3a \n");
             // check if the phase ratio of a dike phase is greater than 0 in the current cell
            if(phRat[i]>0)
              {
		PetscPrintf(PETSC_COMM_WORLD,"4a \n");

               if(dike->Mb == dike->Mf)
                 {
		   PetscPrintf(PETSC_COMM_WORLD,"5a \n");
                  // constant M
                  M = dike->Mf;
                  v_spread = PetscAbs(bc->velin);
                  left = PhaseTrans->bounds[0];
                  right = PhaseTrans->bounds[1];
                  dike->dikeRHS = M * 2 * v_spread / PetscAbs(left-right);  // necessary to write dike->dikeRHS?
                 }

	  /*else                                                                                                                                                          
            {
	    // Mb an Mf are different
                // FDSTAG *fs;

                // access context
                // fs = bc->fs;
		// bdx = SIZE_NODE(i, sx, fs->dsx); // distance between two neighbouring cell centers in x-direction 
                //  cdx = SIZE_CELL(i, sx, fs->dsx); // distance between two neigbouring nodes in x-direction       
                                             
                if(front == back)
                {
                    // linear interpolation between different M values, Mf is M in front, Mb is M in back
                    M = dike.Mf + (dike.Mb - dike.Mf) * (y/(PetscAbs(front+back))); 
                    dikeRHS = M * 2 * v_spread / PetscAbs(left+right);  // [1/s] SCALE THIS TERM, now it is in km
		}
                else
                {
                    // linear interpolation if the ridge/dike phase is oblique
                    y = COORD_CELL(j,sy,fs->dsy);
                    M = Mf + (Mb - Mf) * (y/(PetscAbs(front+back)));
                    dikeRHS = M * 2 * v_spread / PetscAbs(left+right);  // [1/s] SCALE THIS TERM, now it is in km 
                }
            }*/
            else
            {
              dike->dikeRHS = 0.0;   // necessary dike->dikeRHS ?? not really right? it is always passed as a variable
            }

             dikeRHS += phRat[i]*dike->dikeRHS;   // is it correct to just use dikeRHS? still necessary to save as dike->dikeRHS before because used in cellconsteq?

	      }
        }
    PetscFunctionReturn(0);

}
