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
#include "tssolve.h"
#include "scaling.h"
#include "fdstag.h"
#include "tools.h"
#include "surf.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBDikeCreate"
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)   
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
            if(fb->nblocks >_max_num_dike_)
            {
                        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many dikes specified! Max allowed: %lld", (LLD)_max_num_dike_ );
            }

            // store actual number of dike blocks 
            dbdike->numDike = fb->nblocks;

            if (PrintOutput){
                        PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
            }
		
            // read each individual dike block                                                                                                                   
            for(jj = 0; jj < fb->nblocks; jj++)
            {
                ierr = DBReadDike(dbdike, dbm, fb, jr, PrintOutput); CHKERRQ(ierr);
                fb->blockID++;
            }
        }
 
	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBReadDike"
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
  // read dike parameter from file 
  Dike     *dike;
  PetscInt  ID;
  Scaling  *scal;
	
  PetscErrorCode ierr;
  PetscFunctionBegin;

	// access context           
  scal = dbm->scal;

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

	// set default value for Mc in case no Mc is provided
	dike->Mc = -1.0;
	// set default value for y_Mc in case it is not used (it does not matter since it is not accessed and checked for anywhere but might be better than not setting it)
	dike->y_Mc = 0.0;

	// read and store dike  parameters. 
  ierr = getScalarParam(fb, _REQUIRED_, "Mf",      &dike->Mf,      1, 1.0);              CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Mc",      &dike->Mc,      1, 1.0);              CHKERRQ(ierr);
  ierr = getScalarParam(fb, _REQUIRED_, "Mb",      &dike->Mb,      1, 1.0);              CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "y_Mc",    &dike->y_Mc,    1, 1.0);              CHKERRQ(ierr);
	ierr = getIntParam(   fb, _REQUIRED_, "PhaseID", &dike->PhaseID, 1, dbm->numPhases-1); CHKERRQ(ierr);  
	ierr = getIntParam(   fb, _REQUIRED_, "PhaseTransID", &dike->PhaseTransID, 1, dbm->numPhtr-1); CHKERRQ(ierr);
  ierr = getIntParam(   fb, _OPTIONAL_, "dyndike",      &dike->dyndike, 1, dbm->numPhtr-1); CHKERRQ(ierr);

  ierr = getScalarParam(fb, _OPTIONAL_, "Tsol",         &dike->Tsol,    1, 1.0);              CHKERRQ(ierr);
  ierr = getScalarParam(fb, _OPTIONAL_, "zmax_magma",   &dike->zmax_magma,    1, 1.0);              CHKERRQ(ierr);
  ierr = getScalarParam(fb, _OPTIONAL_, "filtx",   &dike->filtx,    1, 1.0);              CHKERRQ(ierr);
  ierr = getScalarParam(fb, _OPTIONAL_, "drhomagma",   &dike->drhomagma,    1, 1.0);              CHKERRQ(ierr);


	// scale the location of Mc y_Mc properly:
	dike->y_Mc /= scal->length;

  if (dike->dyndike)
  {
      ierr = DMCreateLocalVector( jr->DA_CELL_2D, &dike->sxx_eff_ave);  CHKERRQ(ierr);
      ierr = DMCreateLocalVector (jr->DA_CELL_2D, &dike->dPm);  CHKERRQ(ierr);
      ierr = DMCreateLocalVector (jr->DA_CELL_2D, &dike->lthickness);  CHKERRQ(ierr);
  }

  
  if (PrintOutput)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Dike parameters ID[%lld]: PhaseTransID=%i PhaseID=%i Mf=%g, Mb=%g, Mc=%g, y_Mc=%g \n", 
      (LLD)(dike->ID), dike->PhaseTransID, dike->PhaseID, dike->Mf, dike->Mb, dike->Mc, dike->y_Mc, dike->dyndike);
    PetscPrintf(PETSC_COMM_WORLD,"                         : dyndike=%i, Tsol=%g, zmax_magma=%g, filtx=%g, drhomagma=%g \n", 
      dike->dyndike, dike->Tsol, dike->zmax_magma, dike->filtx, dike->drhomagma);
    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
  }

  PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetDikeContr"
PetscErrorCode GetDikeContr(ConstEqCtx *ctx,                                                                                                                                
                            PetscScalar *phRat,          // phase ratios in the control volume   
                            PetscInt &AirPhase,                                                                           
                            PetscScalar &dikeRHS,
                            PetscScalar &y_c,
                            PetscInt J) 
                                              
{
  
  BCCtx       *bc;
  Dike        *dike;
  Ph_trans_t  *CurrPhTr;
  PetscInt     i, nD, nPtr, numDike, numPhtr, nsegs;
  PetscScalar  v_spread, M, left, right, front, back;
  PetscScalar  y_distance, tempdikeRHS;

  PetscFunctionBegin;
  
  numDike    = ctx->numDike;
  bc         = ctx->bc;
  numPhtr    = ctx->numPhtr;


  nPtr = 0;
  nD = 0;

  
  for(nPtr=0; nPtr<numPhtr; nPtr++)   // loop over all phase transitions blocks
    {
      // access the parameters of the phasetranstion block
      CurrPhTr = ctx->PhaseTrans+nPtr;
      
      for(nD = 0; nD < numDike; nD++) // loop through all dike blocks
      {
          // access the parameters of the dike depending on the dike block
          dike = ctx->matDike+nD;
	  
	        // access the phase ID of the dike parameters of each dike
          i = dike->PhaseID;
	  
          if(CurrPhTr->ID == dike->PhaseTransID)  // compare the phaseTransID associated with the dike with the actual ID of the phase transition in this cell           
          {
	           // check if the phase ratio of a dike phase is greater than 0 in the current cell
	           if(phRat[i]>0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
		         {
                nsegs=CurrPhTr->nsegs;
		            if(dike->Mb == dike->Mf && dike->Mc < 0.0)       // constant M
		            {
		               M = dike->Mf;
		               v_spread = PetscAbs(bc->velin);
		               left = CurrPhTr->celly_xboundL[J];
		               right = CurrPhTr->celly_xboundR[J];
		               tempdikeRHS = M * 2 * v_spread / PetscAbs(left-right);
		            }
		            else if(dike->Mc >= 0.0)   // Mf, Mc and Mb
		            {
		               left = CurrPhTr->celly_xboundL[J];
		               right = CurrPhTr->celly_xboundR[J];
		               front = CurrPhTr->ybounds[0];
		               back = CurrPhTr->ybounds[2*nsegs-1];
		               v_spread = PetscAbs(bc->velin);

		               if(y_c >= dike->y_Mc)
			             {
			                 // linear interpolation between different M values, Mc is M in the middle, acts as M in front, Mb is M in back 
			                 y_distance = y_c - dike->y_Mc;
			                 M = dike->Mc + (dike->Mb - dike->Mc) * (y_distance / (back - dike->y_Mc));
			                 tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
			             }
		               else
			             {
			                 // linear interpolation between different M values, Mf is M in front, Mc acts as M in back  
			                 y_distance = y_c - front;
			                 M = dike->Mf + (dike->Mc - dike->Mf) * (y_distance / (dike->y_Mc - front));
			                 tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
			             }
		            }
		            else if(dike->Mb != dike->Mf && dike->Mc < 0.0)   // only Mf and Mb, they are different
		            {
		               left = CurrPhTr->celly_xboundL[J];
		               right = CurrPhTr->celly_xboundR[J];
		               front = CurrPhTr->ybounds[0];
                   back = CurrPhTr->ybounds[2*nsegs-1];

		               v_spread = PetscAbs(bc->velin);
		      
		               // linear interpolation between different M values, Mf is M in front, Mb is M in back
		               y_distance = y_c - front;
		               M = dike->Mf + (dike->Mb - dike->Mf) * (y_distance / (back - front));
		               tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
		            }
		            else                                         // Mb and Mf don't exist (which should not occurr)
		            {
		               tempdikeRHS = 0.0;
		            }
		  
		            dikeRHS += (phRat[i]+phRat[AirPhase])*tempdikeRHS;  // Give full divergence if cell is part dike part air

		        }  //close if phRat and xboundR>xboundL  
	        }  // close phase transition and dike phase ID comparison 
	    }  // close dike block loop
    }  // close phase transition block loop
  
  PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Dike_k_heatsource"
PetscErrorCode Dike_k_heatsource(JacRes *jr,
                                 Material_t *phases,
                                 PetscScalar &Tc,
                                 PetscScalar *phRat,          // phase ratios in the control volume 
                                 PetscScalar &k,
                                 PetscScalar &rho_A,
                                 PetscScalar &y_c,
                                 PetscInt J) 

{
  BCCtx       *bc;
  Dike        *dike;
  Ph_trans_t  *CurrPhTr;
  Material_t  *mat;
  PetscInt     i, numDike, nD, nPtr, numPhtr, nsegs;
  PetscScalar  v_spread, left, right, front, back, M, kfac, tempdikeRHS;
  PetscScalar  y_distance;
  
  PetscFunctionBegin;

  numDike    = jr->dbdike->numDike; // number of dikes
  numPhtr    = jr->dbm->numPhtr;
  bc         = jr->bc;
  
  nPtr = 0;
  nD   = 0;
  kfac = 0;
  
  for(nPtr=0; nPtr<numPhtr; nPtr++)   // loop over all phase transitions blocks                        
    {
      // access the parameters of the phasetranstion block                                             
      CurrPhTr = jr->dbm->matPhtr+nPtr;

      for(nD = 0; nD < numDike; nD++) // loop through all dike blocks                                    
        {
          // access the parameters of the dike depending on the dike block                                                    
          dike = jr->dbdike->matDike+nD;

          // access the phase ID of the dike parameters of each dike                                                       
          i = dike->PhaseID;
	  
          if(CurrPhTr->ID == dike->PhaseTransID)  // compare the phaseTransID associated with the dike with the actual ID of the phase transition in this cell
            {

              // if in the dike zone                   
              if(phRat[i]>0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
                {
                  nsegs=CurrPhTr->nsegs;
                  if(dike->Mb == dike->Mf && dike->Mc < 0.0)       // constant M                                  
                    {
                      M = dike->Mf;
                      v_spread = PetscAbs(bc->velin);
                      left = CurrPhTr->celly_xboundL[J];
                      right = CurrPhTr->celly_xboundR[J];
                      tempdikeRHS = M * 2 * v_spread / PetscAbs(left-right);
		                }
		              else if(dike->Mc >= 0.0)   // Mf, Mc and Mb            
                    {
                      left = CurrPhTr->celly_xboundL[J];
                      right = CurrPhTr->celly_xboundR[J];
                      front = CurrPhTr->ybounds[0];
                      back = CurrPhTr->ybounds[2*nsegs-1];
                      v_spread = PetscAbs(bc->velin);

                      if(y_c >= dike->y_Mc)
                        {
                          // linear interpolation between different M values, Mc is M in the middle, acts as M in front, Mb is M in back 
                          y_distance = y_c - dike->y_Mc;
                          M = dike->Mc + (dike->Mb - dike->Mc) * (y_distance / (back - dike->y_Mc));
                          tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
                        }
                      else
                        {
                          // linear interpolation between different M values, Mf is M in front, Mc acts as M in back
                          y_distance = y_c - front;
                          M = dike->Mf + (dike->Mc - dike->Mf) * (y_distance / (dike->y_Mc - front));
                          tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
                        }
                    }
                  else if(dike->Mb != dike->Mf && dike->Mc < 0.0)   // only Mf and Mb, they are different     
                    {
                      left = CurrPhTr->celly_xboundL[J];
                      right = CurrPhTr->celly_xboundR[J];
                      front = CurrPhTr->ybounds[0];
                      back = CurrPhTr->ybounds[2*nsegs-1];
                      v_spread = PetscAbs(bc->velin);

                      // linear interpolation between different M values, Mf is M in front, Mb is M in back       
                      y_distance = y_c - front;
                      M = dike->Mf + (dike->Mb - dike->Mf) * (y_distance / (back - front));
                      tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
                    }
		              else
		               {
		                  tempdikeRHS = 0.0;
		               } 
		  
		              mat = &phases[i];
		  
		              //adjust k and heat source according to Behn & Ito [2008]
		              if (Tc < mat->T_liq && Tc > mat->T_sol)
		               {
		                 kfac  += phRat[i] / ( 1 + ( mat->Latent_hx/ (mat->Cp*(mat->T_liq-mat->T_sol))) );
		                 rho_A += phRat[i]*(mat->rho*mat->Cp)*(mat->T_liq-Tc)*tempdikeRHS;  // Cp*rho not used in the paper, added to conserve units of rho_A
		               }
		              else if (Tc <= mat->T_sol)
		               {
		                 kfac  += phRat[i];
		                 rho_A += phRat[i]*( mat->rho*mat->Cp)*( (mat->T_liq-Tc) + mat->Latent_hx/mat->Cp )*tempdikeRHS;
		               }
		              else if (Tc >= mat->T_liq)
		               {
		                 kfac += phRat[i];
		               }
		              // end adjust k and heat source according to Behn & Ito [2008]
		  
		              k=kfac*k;
		  
		            }   // end if phRat and xboundR>xboundL
	          } // close phase transition and phase ID comparison	  
	      }   // end dike block loop      
    }  // close phase transition block loop
  
  PetscFunctionReturn(0);
}


//------------------------------------------------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "Locate_Dike_Zones"
PetscErrorCode Locate_Dike_Zones(JacRes *jr)
{

  Controls    *ctrl;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ctrl = &jr->ctrl;

  if (!ctrl->actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active

  ierr = Compute_sxx_eff(jr);  //compute mean effective sxx across the lithosphere

  PetscPrintf(PETSC_COMM_WORLD,"Out of Compute_sxx_eff \n");

  PetscFunctionReturn(0); 
  ierr = Smooth_sxx_eff(jr);   //smooth mean effective sxx
  
  
  PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "Compute_sxx_eff"
PetscErrorCode Compute_sxx_eff(JacRes *jr)
{
  MPI_Request srequest, rrequest;
  Vec         vbuff, vbuff2, vbuff3;
  PetscScalar ***gsxx_eff_ave;
  PetscScalar ***ibuff,***ibuff2, ***ibuff3;
  PetscScalar  *lbuff, *lbuff2, *lbuff3;
  PetscScalar dz, cumk, cumk2, ***lT, Tc, *grav, Tsol;
  PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, nD, L, ID, AirPhase, numDike;
  PetscMPIInt    rank;
  PetscScalar ***glthick, ***dPm; //for debugging only

  FDSTAG      *fs;
  Discret1D   *dsz;
  SolVarCell  *svCell;
  Dike        *dike;
  Controls    *ctrl;


/* dPm is magma pressure in excess of dynamic pressure assumed to be the magma-static head 
   at the depth of the solidus if magma pooling begins at a depth of z = zmax_magma
   sxx_eff_ave = total effective horizontal normal stress averaged across the lithosphere (ie., above Tsol).
              sxx_mean=mean(dev. sxx - effective pressure), effective pressure P'=P-Pmagma = -dPm
*/

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ctrl = &jr->ctrl;
  grav = ctrl->grav;


  numDike    = jr->dbdike->numDike; // number of dikes
  fs  =  jr->fs;
  dsz = &fs->dsz;
  L   =  (PetscInt)dsz->rank;
  AirPhase  = jr->surf->AirPhase;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for(nD = 0; nD < numDike; nD++)
  {
    dike = jr->dbdike->matDike+nD;

    if (!dike->dyndike)
    {
      PetscFunctionReturn(0);   // only execute this function if dikes is dynamic
    }
    else
    {
      //printf("ENTERING Compute_sxx_eff \n");
      // much machinery taken from JacResGetLithoStaticPressure
      // get local grid sizes
      ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

      // get integration/communication buffer (Gets a PETSc vector, vbuff, that may be used with the DM global routines)
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vbuff); CHKERRQ(ierr);
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vbuff2); CHKERRQ(ierr);
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vbuff3); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"Got vbuffs\n");

      ierr = VecZeroEntries(vbuff); CHKERRQ(ierr);
      ierr = VecZeroEntries(vbuff2); CHKERRQ(ierr);
      ierr = VecZeroEntries(vbuff3); CHKERRQ(ierr);

      // open index buffer for computation (ibuff the array that shares data with vector vbuff and is indexed with global dimensions<<G.Ito)
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vbuff, &ibuff); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vbuff2, &ibuff2); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vbuff3, &ibuff3); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"Got ibuffs\n");

      // open linear buffer for send/receive  (returns the point, lbuff, that contains this processor portion of vector data, vbuff<<G.Ito)
      ierr = VecGetArray(vbuff, &lbuff); CHKERRQ(ierr);
      ierr = VecGetArray(vbuff2, &lbuff2); CHKERRQ(ierr);
      ierr = VecGetArray(vbuff3, &lbuff3); CHKERRQ(ierr);

      //Access temperatures
      ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"Got lT\n");
      // receive from top domain (next)  dsz->grnext is the next proc up (in increasing z). Top to bottom doesn't matter here, its this way
      // because the code is patterned after GetLithoStaticPressure
      if(dsz->nproc != 1 && dsz->grnext != -1)
      {
        ierr = MPI_Irecv(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lbuff2, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lbuff3, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      }
      
      for(k = sz + nz - 1; k >= sz; k--)
      {
        dz  = SIZE_CELL(k, sz, (*dsz));
        START_PLANE_LOOP
        {
          GET_CELL_ID(ID, i, j, k, nx, ny);
          //svCell = &jr->svCell[ID]; 
          //Tc=lT[k][j][i];
          
          if ((Tc<=dike->Tsol) & (svCell->phRat[AirPhase] < 1.0))
          {
            dz  = SIZE_CELL(k, sz, (*dsz));

            //ibuff[L][j][i]+=svCell->hxx*dz;  //integrating weighted stresses
            ibuff2[L][j][i]+=dz;             //integrating thickeness
          }
          

          //interpolate depth to the solidus
          
          /*if ((k > sz) & (Tc <= dike->Tsol))
          {

            if (dike->Tsol <= lT[k-1][j][i])
            {
               cumk=dsz->ccoor[k-sz]+(dsz->ccoor[k-sz-1]-dsz->ccoor[k-sz])/(lT[k-1][j][i]-Tc)*(Tsol-Tc);
               ibuff3[L][j][i]=cumk;
            }
          } */

        }
        END_PLANE_LOOP
      }
      

      //After integrating and averaging, send it down to the next proc. 
      if(dsz->nproc != 1 && dsz->grprev != -1)
      {
        ierr = MPI_Isend(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lbuff2, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lbuff3, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      }
      
      //Now receive the answer from successive previous (underlying) procs so all procs have the answers
      if(dsz->nproc != 1 && dsz->grprev != -1)
      {
        ierr = MPI_Irecv(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lbuff2, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lbuff3, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      }

      if(dsz->nproc != 1 && dsz->grnext != -1)
      {
        ierr = MPI_Isend(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lbuff2, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lbuff3, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      }

      // (gdev is the array that shares data with devxx_mean and is indexed with global dimensions)
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->dPm, &dPm); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->lthickness, &glthick); CHKERRQ(ierr);

      //now all cores have the same solution so give that to the stress array
      START_PLANE_LOOP
        {
          glthick[L][j][i]=ibuff2[L][j][i];  //Dont need this, but using it to check solution for debugging below
          dPm[L][j][i]=(ibuff3[L][j][i]-dike->zmax_magma)*(dike->drhomagma)*grav[2];  //magmastatic pressure at solidus, note z is negative
          gsxx_eff_ave[L][j][i]=ibuff[L][j][i]/ibuff2[L][j][i]+dPm[L][j][i];  //Depth weighted mean stress + excess magma press.
         }
      END_PLANE_LOOP


      // restore buffer and mean stress vectors
      ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->lthickness, &glthick); CHKERRQ(ierr); //debugging
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->dPm, &dPm); CHKERRQ(ierr); //debugging
      

      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff, &ibuff); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff2, &ibuff2); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff3, &ibuff3); CHKERRQ(ierr);

      ierr = VecRestoreArray(vbuff, &lbuff); CHKERRQ(ierr);
      ierr = VecRestoreArray(vbuff2, &lbuff2); CHKERRQ(ierr);
      ierr = VecRestoreArray(vbuff2, &lbuff3); CHKERRQ(ierr);

      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff); CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff2); CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff3); CHKERRQ(ierr);

      //fill ghost points
      /*
      LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->sxx_eff_ave);
      LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->lthickness);
      LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->dPm);
      */

      PetscPrintf(PETSC_COMM_WORLD,"last plane looop \n");
    //Compute excess magma pressure

    //Locate xbounds
    } //end if dike->dyndike
  } //End loop over dikes

  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "Smooth_sxx_eff"
PetscErrorCode Smooth_sxx_eff(JacRes *jr)
{

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsx, *dsz;
  PetscScalar ***gsxx_eff_ave;
  PetscScalar lxmin, lxmax, filtx;
  PetscInt    i, j, sx, sy, sz, nx, ny, nz, nD, L, numDike;

  PetscScalar ***glthick, ***dPm, dum1, dum2, dum3, dum4; //for debugging only

  PetscErrorCode ierr;
  PetscFunctionBegin;

  fs  =  jr->fs;
  dsx = &fs->dsx;
  dsz = &fs->dsz;
  L   =  (PetscInt)dsz->rank;

  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

  lxmin = dsx->ccoor[0];
  lxmax = dsx->ccoor[nx];
  numDike    = jr->dbdike->numDike; // number of dikes

  for(nD = 0; nD < numDike; nD++) // loop through all dike blocks
  {
       // access the parameters of the dike depending on the dike block
     dike = jr->dbdike->matDike+nD;
     if (!dike->dyndike)
     {
        PetscFunctionReturn(0);   // only execute this function if dikes is dynamic
     }
     else
     {
      //printf("ENTERING Smooth_sxx_eff \n");
       filtx = dike->filtx;

//debugging Just checking to see that we have the values in the arrays for later use
       ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
       ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->lthickness, &glthick); CHKERRQ(ierr);
       ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->dPm, &dPm); CHKERRQ(ierr);
       START_PLANE_LOOP
       {

         if ((j==sy) && (i < sx+5))
         {
            dum1=glthick[L][j][i];;
            dum2=dPm[L][j][i];
            dum3=gsxx_eff_ave[L][j][i]-dPm[L][j][i];
            dum4=gsxx_eff_ave[L][j][i];
            printf("ranks=%i,%i,%i: i,j=%i,%i; lthick=%g, dPm=%g, devxx=%g, sxx_eff_ave=%g \n", fs->dsx.rank,fs->dsy.rank, fs->dsz.rank, i,j,dum1, dum2, dum3, dum4);  //debugging
 
         }
       }
       END_PLANE_LOOP
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->lthickness, &glthick); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->dPm, &dPm); CHKERRQ(ierr);
     }  //end else dyndike
  } //end for loop over numdike

  PetscFunctionReturn(0);  
}  
