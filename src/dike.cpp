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
//#include "scaling.h"
//#include "objFunct.h"  // necessary in case of adjoint
#include "JacRes.h"
//#include "phase_transition.h"   // why not necessary but phase.h yes?, I am not using phases but phase transition boundaries
#include "dike.h"
#include "constEq.h"
#include "bc.h"                 // why necessary and phase_transition.h not? both associated structures are in consteq
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
#define __FUNCT__ "GetDikeContr"
PetscErrorCode GetDikeContr(ConstEqCtx *ctx,
                                  PetscScalar *phRat,          // phase ratios in the control volume
                                  PetscScalar &dikeRHS)
{

        BCCtx       *bc;
        Dike        *matDike;
        Ph_trans_t  *PhaseTrans;
        PetscInt     i, j, numDike;
        PetscScalar  v_spread, M, left, right;

	numDike    = ctx->numDike;
        matDike    = ctx->matDike;

        bc         = ctx->bc;
        PhaseTrans = ctx->PhaseTrans;

          // loop through all dikes
          for(j = 0; j < numDike; j++)
            {
	      // access the phase ID of the dike parameters
              i = matDike->PhaseID;

             // check if the phase ratio of a dike phase is greater than 0 in the current cell
            if(phRat[i]>0)
              {
               if(matDike->Mb == matDike->Mf)
                 {
                  // constant M
                  M = matDike->Mf;
                  v_spread = PetscAbs(bc->velin);
                  left = PhaseTrans->bounds[0];
                  right = PhaseTrans->bounds[1];
                  matDike->dikeRHS = M * 2 * v_spread / PetscAbs(left-right);  // necessary to write dike->dikeRHS?
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
              matDike->dikeRHS = 0.0;   // necessary dike->dikeRHS ?? not really right? it is always passed as a variable
            }

             dikeRHS += phRat[i]*matDike->dikeRHS;   // is it correct to just use dikeRHS? still necessary to save as dike->dikeRHS before because used in cellconsteq?

	      }
        }
    PetscFunctionReturn(0);

}
