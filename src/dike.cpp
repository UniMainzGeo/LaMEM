
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
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)   
{

        // read all dike parameter blocks from file

  Dike     *dike;
  FDSTAG   *fs;
  PetscScalar ***gsxx_eff_ave_hist;
  PetscInt jj, nD, numDike, numdyndike, istep_nave;
  PetscInt i, j, istep_count, sx, sy, sisc, nx, ny;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  if (!jr->ctrl.actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active
 
  fs = jr->fs;
  //===============                                                                                                                                               
  // DIKE PARAMETER                                                                                                               
  //===============                                                                                                                                               

  // setup block access mode                                                                                                                                      
  ierr = FBFindBlocks(fb, _OPTIONAL_, "<DikeStart>", "<DikeEnd>"); CHKERRQ(ierr);

  if(fb->nblocks)
  {
      // print overview of dike blocks from file                                                                                                           
      if (PrintOutput)
        PetscPrintf(PETSC_COMM_WORLD,"Dike blocks : \n");

      // initialize ID for consistency checks                                                                                                                 

      for(jj = 0; jj < _max_num_dike_ ; jj++) dbdike->matDike[jj].ID = -1;

      // error checking
      if(fb->nblocks >_max_num_dike_)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many dikes specified! Max allowed: %lld", (LLD)_max_num_dike_ );

      // store actual number of dike blocks 
      dbdike->numDike = fb->nblocks;

      if (PrintOutput)
        PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
		
      // read each individual dike block                                                                                                                   
      for(jj = 0; jj < fb->nblocks; jj++)
      {
        ierr = DBReadDike(dbdike, dbm, fb, jr, PrintOutput); CHKERRQ(ierr);
        fb->blockID++;
      }
  }
 
	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

  numdyndike=0;

  numDike=dbdike->numDike;
  for(nD = 0; nD < numDike; nD++) // loop through all dike blocks
  {
      dike = dbdike->matDike+nD;
      if (dike->dyndike_start)
      {
        numdyndike++;
        if (numdyndike ==1)
        {
            // DM for 1D cell center vector  (take this out of this loop because it will be repeated with >1 dynamic dike)
            ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
            fs->dsx.tcels, fs->dsy.nproc, fs->dsz.nproc, 
            fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc, 1, 1,
            0, 0, 0, &jr->DA_CELL_1D); CHKERRQ(ierr);


            //DM for 2D cell center vector, with istep_nave planes for time averaging
            ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
            fs->dsx.tcels, fs->dsy.tcels, fs->dsz.nproc*dike->istep_nave, 
            fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc, 1, 1,
            0, 0, 0, &jr->DA_CELL_2D_tave); CHKERRQ(ierr);
        }

        //creating local vectors and inializing the history vector
        ierr = DMCreateLocalVector( jr->DA_CELL_2D, &dike->sxx_eff_ave);  CHKERRQ(ierr);
        ierr = DMCreateLocalVector( jr->DA_CELL_2D_tave, &dike->sxx_eff_ave_hist);  CHKERRQ(ierr);
        ierr = DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);
        ierr = DMDAGetCorners(jr->DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, &istep_nave); CHKERRQ(ierr);    

        for(j = sy; j < sy+ny; j++) 
        {
          for(i = sx; i < sx+nx; i++)
          {
            for (istep_count=sisc; istep_count<sisc+istep_nave; istep_count++)
            {
              gsxx_eff_ave_hist[istep_count][j][i]=0.0;
            }
          }
        }

        ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);
      }  //End if dyndike->start
  }  //End loop through dikes


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
  // read dike parameter from file 
  Dike     *dike;
  FDSTAG   *fs;
  PetscInt  ID;
  Scaling  *scal;
  	
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

	// access context           
  scal = dbm->scal;
  fs = jr->fs;

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
  ierr = getIntParam(   fb, _OPTIONAL_, "dyndike_start",&dike->dyndike_start, 1, -1); CHKERRQ(ierr);

  if (dike->dyndike_start)
  {
    dike->Tsol = 1000;
    dike->zmax_magma = -15.0;
    dike->drhomagma = 500;
    dike->npseg = 4.0; 
    dike->filtx = 2.0;
    dike->istep_nave = 2;

    ierr = getScalarParam(fb, _OPTIONAL_, "Tsol",         &dike->Tsol,         1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "zmax_magma",   &dike->zmax_magma,   1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "filtx",        &dike->filtx,        1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "drhomagma",    &dike->drhomagma,    1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "npseg",        &dike->npseg,        1, 1.0); CHKERRQ(ierr);
    ierr = getIntParam(fb, _OPTIONAL_, "istep_nave",     &dike->istep_nave,        1, 50); CHKERRQ(ierr);

    dike->istep_count=dike->istep_nave;   //initialize so that when istep=dike_start, it is set to 0

    dike->npseg0=fs->dsy.tcels - floor((PetscScalar)fs->dsy.tcels/dike->npseg)*dike->npseg;
    if (dike->npseg0==0) dike->npseg0=dike->npseg;

  }


	// scale the location of Mc y_Mc properly:
	dike->y_Mc /= scal->length;

  
  if (PrintOutput)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Dike parameters ID[%lld]: PhaseTransID=%i PhaseID=%i Mf=%g, Mb=%g, Mc=%g, y_Mc=%g \n", 
      (LLD)(dike->ID), dike->PhaseTransID, dike->PhaseID, dike->Mf, dike->Mb, dike->Mc, dike->y_Mc);
    if (dike->dyndike_start)
    {
      PetscPrintf(PETSC_COMM_WORLD,"                         dyndike_start=%i, Tsol=%g, zmax_magma=%g, filtx=%g, drhomagma=%g \n", 
        dike->dyndike_start, dike->Tsol, dike->zmax_magma, dike->filtx, dike->drhomagma);
      PetscPrintf(PETSC_COMM_WORLD,"                         npseg=%g, npseg0=%g, istep_nave=%i, istep_count=%i \n",
        dike->npseg, dike->npseg0, dike->istep_nave, dike->istep_count);

    }
    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");    
  }

   if (dike->dyndike_start)
  {
      dike->Tsol = (dike->Tsol +  jr->scal->Tshift)/jr->scal->temperature;
      dike->filtx /= jr->scal->length;
      dike->drhomagma /= jr->scal->density;
      dike->zmax_magma /= jr->scal->length;
  }

  PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------------
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

  PetscFunctionBeginUser;
  
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
  
  PetscFunctionBeginUser;

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
PetscErrorCode Locate_Dike_Zones(JacRes *jr)
{
//Ugh:  history arrays are from before the last timestep and so stresses used here are a step behind

  Controls    *ctrl;
  Dike        *dike;
  PetscInt   nD, numDike, iwrite_counter1, iwrite_counter2, icounter;
  PetscErrorCode ierr; 

  PetscFunctionBeginUser;

  ctrl = &jr->ctrl;
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  if (!ctrl->actDike || jr->ts->istep+1 == 0) PetscFunctionReturn(0);   // only execute this function if dikes are active

  numDike    = jr->dbdike->numDike; // number of dikes
  iwrite_counter1=0;
  iwrite_counter2=0;
  icounter=0;

  for(nD = 0; nD < numDike; nD++)
  {
    dike = jr->dbdike->matDike+nD;

    if (dike->dyndike_start && jr->ts->istep+1 >= dike->dyndike_start)
    {
       // compute lithostatic pressure
       if (icounter==0) ierr = JacResGetLithoStaticPressure(jr); CHKERRQ(ierr);
       icounter++;
       
       ierr = Compute_sxx_eff(jr,nD, iwrite_counter1); CHKERRQ(ierr);  //compute mean effective sxx across the lithosphere

       ierr = Smooth_sxx_eff(jr,nD, iwrite_counter2); CHKERRQ(ierr);  //smooth mean effective sxx  //debugging

       ierr = Set_dike_zones(jr,nD); CHKERRQ(ierr); //centered on peak sxx_eff_ave
    }
    else
    {
       PetscFunctionReturn(0); 
    }

  }
  PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------------

PetscErrorCode Compute_sxx_eff(JacRes *jr, PetscInt nD, PetscInt iwrite_counter)
{
  MPI_Request srequest, rrequest;
  Vec         vsxx, vliththick, vzsol;
  PetscScalar ***gsxx_eff_ave, ***p_lith;
  PetscScalar ***sxx,***liththick, ***zsol;
  PetscScalar  *lsxx, *lliththick, *lzsol;
  PetscScalar dz, ***lT, Tc, *grav, Tsol, dPmag;
  PetscScalar dbug1, dbug2, dbug3, xcell, ycell;
  PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, L, ID, AirPhase;
  PetscMPIInt    rank;

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz;
  SolVarCell  *svCell;
  Controls    *ctrl;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;


  dbug1=(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out);  //debugging
  dbug2=floor(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out); //debugging
  dbug3=dbug1-dbug2; //debugging

  ctrl = &jr->ctrl;
  grav = ctrl->grav;


  fs  =  jr->fs;
  dsz = &fs->dsz;
  //dsx = &fs->dsx;  //debugging
  //scal = fs->scal; //debugging
  L   =  (PetscInt)dsz->rank;
  AirPhase  = jr->surf->AirPhase;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);


  dike = jr->dbdike->matDike+nD;

       // get local grid sizes
      ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

      // get integration/communication buffer (Gets a PETSc vector, vsxx, that may be used with the DM global routines)
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vliththick); CHKERRQ(ierr);
      ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vzsol); CHKERRQ(ierr);

      ierr = VecZeroEntries(vsxx); CHKERRQ(ierr);
      ierr = VecZeroEntries(vliththick); CHKERRQ(ierr);
      ierr = VecZeroEntries(vzsol); CHKERRQ(ierr);

      // open index buffer for computation (sxx the array that shares data with vector vsxx and is indexed with global dimensions<<G.Ito)
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);



      // open linear buffer for send/receive  (returns the point, lsxx, that contains this processor portion of vector data, vsxx<<G.Ito)
      ierr = VecGetArray(vsxx, &lsxx); CHKERRQ(ierr);
      ierr = VecGetArray(vliththick, &lliththick); CHKERRQ(ierr);
      ierr = VecGetArray(vzsol, &lzsol); CHKERRQ(ierr);

      //Access temperatures
      ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);

      // receive from top domain (next)  dsz->grnext is the next proc up (in increasing z). Top to bottom doesn't matter here, its this way
      // because the code is patterned after GetLithoStaticPressure
      if(dsz->nproc != 1 && dsz->grnext != -1)
      {
        ierr = MPI_Irecv(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      }
      Tsol=dike->Tsol;

      for(k = sz + nz - 1; k >= sz; k--)
      {
        dz  = SIZE_CELL(k, sz, (*dsz));
        START_PLANE_LOOP
          GET_CELL_ID(ID, i-sx, j-sy, k-sz, nx, ny);  //GET_CELL_ID needs local indices
          svCell = &jr->svCell[ID]; 
          Tc=lT[k][j][i];
 
          
          if ((Tc<=Tsol) & (svCell->phRat[AirPhase] < 1.0))
          {
            dz  = SIZE_CELL(k, sz, (*dsz));

            sxx[L][j][i]+=(svCell->hxx - svCell->svBulk.pn + p_lith[k][j][i])*dz;  //integrating dz-weighted non-lithostatic stress

            liththick[L][j][i]+=dz;             //integrating thickeness
          }
          
          //interpolate depth to the solidus
          if ((Tc <= Tsol) & (Tsol < lT[k-1][j][i]))
          {
            zsol[L][j][i]=dsz->ccoor[k-sz]+(dsz->ccoor[k-sz-1]-dsz->ccoor[k-sz])/(lT[k-1][j][i]-Tc)*(Tsol-Tc); 
               //zsol[L][j][i]=dsz->ccoor[k-sz]+(dsz->ccoor[k-sz-1]-dsz->ccoor[k-sz]);
          }
          //if (j==1 && i==61) printf("k=%i, Tc=%g, Tsol=%g, lT=%g, zsol=%g \n",k,Tc,Tsol,lT[k-1][j][i],zsol[L][j][i]);
        END_PLANE_LOOP
      }
      

      //After integrating and averaging, send it down to the next proc. 
      if(dsz->nproc != 1 && dsz->grprev != -1)
      {
        ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      }
      
      //Now receive the answer from successive previous (underlying) procs so all procs have the answers
      if(dsz->nproc != 1 && dsz->grprev != -1)
      {
        ierr = MPI_Irecv(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Irecv(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      }

      if(dsz->nproc != 1 && dsz->grnext != -1)
      {
        ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

        ierr = MPI_Isend(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
        ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      }

      // (gdev is the array that shares data with devxx_mean and is indexed with global dimensions)
      ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);

      //now all cores have the same solution so give that to the stress array
      START_PLANE_LOOP

        dPmag=-(zsol[L][j][i]-dike->zmax_magma)*(dike->drhomagma)*grav[2];  //excess magmastatic pressure at solidus, note z AND grav[2] <0
        if (dPmag<0) dPmag=dPmag*1e3;                                  //Keep dike over the magma. But caution with Smooth_sxx_eff    
        gsxx_eff_ave[L][j][i]=sxx[L][j][i]/liththick[L][j][i]+dPmag;  //Depth weighted mean effective stress: (total stress)+(magma press)
                                                                      // (magma press)=lp_lith+dPmagma
      END_PLANE_LOOP  

      if (L==0 && dbug3 < 0.05/jr->ts->nstep_out && iwrite_counter==0)  //debugging
      {
        iwrite_counter++;
        START_PLANE_LOOP
          xcell=COORD_CELL(i, sx, fs->dsx);
          ycell=COORD_CELL(j, sy, fs->dsy);
          dPmag=-(zsol[L][j][i]-dike->zmax_magma)*(dike->drhomagma)*grav[2];  //effective pressure is P-Pmagma= negative of magmastatic pressure at solidus, note z AND grav[2] <0
          if (dPmag<0) dPmag=dPmag*1e3;                                  //Keep dike over the magma. But caution with Smooth_sxx_eff    
          PetscSynchronizedPrintf(PETSC_COMM_WORLD,"101010.1010 %i %g %g %g %g\n", jr->ts->istep+1,xcell, ycell, gsxx_eff_ave[L][j][i], dPmag);   //debugging    
        END_PLANE_LOOP  
      }
      PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);   //debugging    
 

      // restore buffer and mean stress vectors
      ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);      

      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);

      ierr = VecRestoreArray(vsxx, &lsxx); CHKERRQ(ierr);
      ierr = VecRestoreArray(vliththick, &lliththick); CHKERRQ(ierr);
      ierr = VecRestoreArray(vzsol, &lzsol); CHKERRQ(ierr);

      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vliththick); CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vzsol); CHKERRQ(ierr);

      //fill ghost points
      
      LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->sxx_eff_ave);

  ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------
// Apply a box filter to the depth-averaged effective stress.  
// **NOTE** There is NO message passing between adjacent procs in x, so this wont work well if the zone of high 
// stress is spit by a processor boundary. This will work completely if cpu_x = 1

PetscErrorCode Smooth_sxx_eff(JacRes *jr, PetscInt nD, PetscInt iwrite_counter)
{

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz, *dsy;
  PetscScalar ***gsxx_eff_ave, ***gsxx_eff_ave_hist;
  PetscScalar x, sum_sxx, sum_dx, dx, xx;
  PetscInt    i, ii, j, sx, sy, sz, nx, ny, nz;
  PetscInt    L, M, rank, nseg, npseg, npseg0, jback, jj;
  PetscInt    sisc, istep_count, istep_nave;
  Vec         vvec1d;
  PetscScalar ***vec1d, *lvec1d;
  PetscScalar xcell, ycell, dbug1, dbug2, dbug3;

  //Scaling     *scal;  //debugging

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Request srequest, rrequest;

  fs  =  jr->fs;
  dsz = &fs->dsz;
  dsy = &fs->dsy;
  L   =  (PetscInt)dsz->rank;
  M   =  (PetscInt)dsy->rank;
  //scal = fs->scal;  //debugging
  dbug1=(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out);  //debugging
  dbug2=floor(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out); //debugging
  dbug3=dbug1-dbug2; //debugging

  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

  dike = jr->dbdike->matDike+nD;

//-------------------------------------------------------------------------------------------------------
// Set up a temporary (1D) array for
// (1) 1st: storing data while smoothing across x and..
// (2) 2nd: storing data to pass between proc for averaging in y
// get communication buffer (Gets a PETSc vector, vvec1d, that may be used with the DM global routines)
       ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vvec1d); CHKERRQ(ierr);
       ierr = VecZeroEntries(vvec1d); CHKERRQ(ierr);
// open index buffer for computation (vec1d is the array that shares data with vector vvec1d & indexed with global dimensions)
       ierr = DMDAVecGetArray(jr->DA_CELL_1D, vvec1d, &vec1d); CHKERRQ(ierr);
// open linear buffer for send/receive  (returns the pointer, lsxx..., that contains this processor portion of vector data, vvec1d)
       ierr = VecGetArray(vvec1d, &lvec1d); CHKERRQ(ierr);
//-------------------------------------------------------------------------------------------------------

//access depth-averaged effective sxx array
       ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
//---------------------------------------------------------------------------------------------
//  moving box filter in x
//---------------------------------------------------------------------------------------------
       for(j = sy; j < sy+ny; j++) 
       {
          for(i = sx; i < sx+nx; i++)
          {
            x = COORD_CELL(i, sx, fs->dsx);
            sum_sxx=0.0;
            sum_dx=0.0;

            for(ii = sx; ii < sx+nx; ii++) 
            {
               xx = COORD_CELL(ii, sx, fs->dsx);
               if (fabs(xx-x) <= 2*dike->filtx)
               {                 
                  dx  = SIZE_CELL(ii, sx, fs->dsx);
                  sum_sxx+=gsxx_eff_ave[L][j][ii]*exp(-0.5*pow(((xx-x)/dike->filtx),2))*dx;
                  sum_dx+=exp(-0.5*pow(((xx-x)/dike->filtx),2))*dx;
                }     
            }
            vec1d[L][M][i]=sum_sxx/sum_dx;

            /*if (L==0 && dbug3 < 0.05/jr->ts->nstep_out)  //debugging
            {  
               ydebugging = COORD_CELL(j, sy, fs->dsy);  //debugging           
               x = COORD_CELL(i, sx, fs->dsx);
               printf("10101010 %i %g %g %g %g %i %g \n", jr->ts->istep+1,x, ydebugging, gsxx_eff_ave[L][j][i], sum_sxx);   //debugging 
            } */
          } //end 1st loop over x to build filtered row

          for(i = sx; i < sx+nx; i++)  //set the smoothed values into permanent array here
          {
            gsxx_eff_ave[L][j][i]=vec1d[L][M][i];
            vec1d[L][M][i]=0.0;  //clean this up for use below
          }
       } //end 1st loop over j for smoothing in x
//---------------------------------------------------------------------------------------------
//  averaging in y-increments (i.e., "segments") along axis
//---------------------------------------------------------------------------------------------
       npseg0=(PetscInt)dike->npseg0;

       for (j=sy; j<sy+ny; j++)
       {
         nseg=floor(((PetscScalar)j-dike->npseg0)/dike->npseg)+2;  //the current moveable dike segment
         if (nseg==1)  //the first global segment could have a different number of cells
         {
           jj=j+1;                         //cell number within segment 
           npseg=npseg0;
         }
         else
         {
           npseg=(PetscInt)dike->npseg;
           jj=j-((nseg-2)*npseg+npseg0)+1;  //cell number within segment 

         }

         //if starting within segmented from previous proc, get array from previous proc
         if (j==sy && jj>1 && dsy->nproc > 1 && dsy->grprev != -1)  
         {
           ierr = MPI_Irecv(lvec1d, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
           ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
         }

         for (i=sx; i<sx+nx; i++) //for each value of i, over j for jj cells 
         {
           if (jj==1) vec1d[L][M][i]=0;

           vec1d[L][M][i]=vec1d[L][M][i]+gsxx_eff_ave[L][j][i]/(PetscScalar)npseg;
 
           if  (jj==npseg && j-sy>=npseg-1) //if at end of segment and segment is contained in current proc
           {
              for (jback=0; jback<jj; jback++)  //fill in the mean values for all cells in segment
              {
                gsxx_eff_ave[L][j-jback][i]=vec1d[L][M][i];
              }//end loop through jback
           }
         } //end loop over i

         // if at end of segment that started in previous core then pass the array back
         if ((jj==npseg) && (j-sy < npseg-1) && (dsy->nproc != 1) &&  (dsy->grprev != -1))
         {
           for (jback=0; jback<=j-sy; jback++)  //fill in the mean values for all cells in segment
           {
             for (i=sx; i<=sx+nx; i++)
             {
               gsxx_eff_ave[L][j-jback][i]=vec1d[L][M][i];
             }
           }

           ierr = MPI_Isend(lvec1d, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
           ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
         }

         //if at end of y domain and not end of segment send to next and receive from next
         if ((j-sy==ny-1) && (jj<npseg) && (dsy->nproc != 1) &&  (dsy->grnext != -1))  
         {
           ierr = MPI_Isend(lvec1d, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
           ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
 
           ierr = MPI_Irecv(lvec1d, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
           ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
           for (jback=0; jback<jj; jback++) //fill in mean values in segment received from next core
           {
             for (i=sx; i<=sx+nx; i++)
             {
               gsxx_eff_ave[L][j-jback][i]=vec1d[L][M][i];           
             }
           }
         }

       } //end 2nd loop over j for averaging in y-increments

       ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vvec1d, &vec1d); CHKERRQ(ierr);
       ierr = VecRestoreArray(vvec1d, &lvec1d); CHKERRQ(ierr);
       ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vvec1d); CHKERRQ(ierr);

//---------------------------------------------------------------------------------------------
//  averaging over time
//---------------------------------------------------------------------------------------------
       ierr = DMDAGetCorners(jr->DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, &istep_nave); CHKERRQ(ierr);

       if(istep_nave!=dike->istep_nave) 
       {
          SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Problems: istep_nave=%i, dike->istep_nave=%i\n", istep_nave, dike->istep_nave);
       }

       ierr = DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);

       dike->istep_count++;

       if (dike->istep_count+1 > istep_nave) 
       {
         dike->istep_count=0;
       }  
    
       START_PLANE_LOOP
         sum_sxx=0;
         gsxx_eff_ave_hist[sisc+dike->istep_count][j][i]=gsxx_eff_ave[L][j][i];  //array for current step

         for (istep_count=sisc; istep_count<sisc+istep_nave; istep_count++)
         {
           sum_sxx+=gsxx_eff_ave_hist[istep_count][j][i];
         }

         gsxx_eff_ave[L][j][i]=sum_sxx/((PetscScalar)istep_nave);
 
       END_PLANE_LOOP

       if (L==0 && dbug3 < 0.05/jr->ts->nstep_out && iwrite_counter==0)  //debugging
       { 
          iwrite_counter++;
          START_PLANE_LOOP
            xcell=COORD_CELL(i, sx, fs->dsx);
            ycell=COORD_CELL(j, sy, fs->dsy);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"202020.2020 %i %g %g %g\n", jr->ts->istep+1,xcell, ycell, gsxx_eff_ave[L][j][i]);   //debugging    
          END_PLANE_LOOP  
       }            
       PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);   //debugging    


       //restore arrays
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);

  PetscFunctionReturn(0);  
}  

//----------------------------------------------------------------------------------------------------
// Set bounds of NotInAir box based on peak sxx_eff_ave
// NOTE that NOW, this only works if cpu_x =1
//

PetscErrorCode Set_dike_zones(JacRes *jr, PetscInt nD)
{

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsx, *dsz;
  Ph_trans_t  *CurrPhTr;
  PetscScalar ***gsxx_eff_ave;
  PetscScalar xcenter, sxx_max, dike_width, mindist, xshift, xcell, dx;
  PetscInt    i, lj, j, sx, sy, sz, nx, ny, nz, nPtr, numPhtr, L, Lx, numDike, ixcenter;
  PetscScalar dbug1, dbug2, dbug3, ydebugging;
  PetscScalar sxxm, sxxp, dx12, dsdx1, dsdx2, x_maxsxx;   
  PetscInt    ixmax;
 
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  fs  =  jr->fs;
  dsz = &fs->dsz;
  dsx = &fs->dsx;
  L   =  dsz->rank;
  Lx  =  dsx->rank;
  numDike    = jr->dbdike->numDike; // number of dikes
  numPhtr    = jr->dbm->numPhtr;

  dbug1=(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out);  //debugging
  dbug2=floor(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out); //debugging
  dbug3=dbug1-dbug2; //debugging

  dike = jr->dbdike->matDike+nD;

       if (Lx>0)
       {
         PetscPrintf(PETSC_COMM_WORLD,"Set_dike_zones requires cpu_x = 1 Lx = %i \n", Lx);
         SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Set_dike_zones requires cpu_x = 1 Lx = %i \n", Lx);
       }
       ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
       ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

       for(nPtr=0; nPtr<numPhtr; nPtr++)   // loop over all phase transitions blocks                        
       {                                               
         CurrPhTr = jr->dbm->matPhtr+nPtr;
         if(CurrPhTr->ID == dike->PhaseTransID)  // compare the phaseTransID associated with the dike with the actual ID of the phase transition in this cell
         {
            for(lj = 0; lj < ny; lj++)  //local index
            {
              sxx_max=-1e12;
              mindist=1e12;
              ixcenter = 0;
              xshift=0;
              ixmax=sx+1;

              j=sy+lj;  //global index
              dike_width=CurrPhTr->celly_xboundR[lj]-CurrPhTr->celly_xboundL[lj];
              xcenter=(CurrPhTr->celly_xboundR[lj] + CurrPhTr->celly_xboundL[lj])/2;
              ydebugging=COORD_CELL(j, sy, fs->dsy);  //debugging

              for(i=sx+1; i < sx+nx-1; i++) 
              {
                xcell=COORD_CELL(i, sx, fs->dsx);
                if (fabs(xcell-xcenter) <= mindist) //find indice of dike zone center (xcenter)
                {
                  ixcenter=i;
                  mindist=fabs(xcell-xcenter);
                }
              } //end loop to find ixcenter
 
              for(i=ixcenter-2; i <= ixcenter+2; i++) //find max gsxx_eff at each value of y
              {
                if ((gsxx_eff_ave[L][j][i] > sxx_max))
                {
                  sxx_max=gsxx_eff_ave[L][j][i];
                  //xshift=COORD_CELL(i, sx, fs->dsx)-xcenter;
                  ixmax=i;
                }
   
              } 
              //finding where slope of dsxx/dx=0
              sxxm =  gsxx_eff_ave[L][j][ixmax-1];  //left of maximum point
              sxxp =  gsxx_eff_ave[L][j][ixmax+1]; ;  //right of max. point
 
              dsdx1=2*(sxx_max-sxxm)/(SIZE_CELL(ixmax-1, sx, fs->dsx)+SIZE_CELL(ixmax, sx, fs->dsx));  //slope left of max
              dsdx2=2*(sxxp-sxx_max)/(SIZE_CELL(ixmax+1, sx, fs->dsx)+SIZE_CELL(ixmax, sx, fs->dsx));  //slope right of max
              dx12=(COORD_CELL(ixmax+1, sx, fs->dsx)-COORD_CELL(ixmax-1, sx, fs->dsx))/2;

              if ((dsdx1>0) & (dsdx2<0))  //if local maximum, interpolate to find where dsdx=0;
              {
                x_maxsxx=(COORD_CELL(ixmax-1, sx, fs->dsx)+COORD_CELL(ixmax, sx, fs->dsx))/2-dsdx1/(dsdx2-dsdx1)*dx12;
              }
              else  //just higher on either side of dike
              {
                x_maxsxx=COORD_CELL(ixmax,sx,fs->dsx);
              }

              xshift=x_maxsxx-xcenter;


              dx=SIZE_CELL(ixcenter,sx, fs->dsx);  
              if (xshift>0 && fabs(xshift) > 0.5*SIZE_CELL(ixcenter, sx, fs->dsx)) //ensure new center is within width of cell to right of center
              {
                xshift=0.5*SIZE_CELL(ixcenter, sx, fs->dsx);
              }
              else if (xshift<0 && fabs(xshift) > 0.5*SIZE_CELL(ixcenter-1, sx, fs->dsx)) //ensure its within the width of cell left of center
              {

                xshift=-0.5*SIZE_CELL(ixcenter-1, sx, fs->dsx);
              }

              CurrPhTr->celly_xboundL[lj]=xcenter+xshift-dike_width/2;  
              CurrPhTr->celly_xboundR[lj]=xcenter+xshift+dike_width/2;  

              if (L==0 && dbug3 < 0.05/jr->ts->nstep_out)   //debugging
              {
                ydebugging = COORD_CELL(j, sy, fs->dsy);  //debugging
                xcell=(COORD_CELL(ixmax-1, sx, fs->dsx)+COORD_CELL(ixmax, sx, fs->dsx))/2;
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"303030.3030 %i %g %g %g %g %g %g %g\n", jr->ts->istep+1, ydebugging, xcenter+xshift, 
                CurrPhTr->celly_xboundL[lj], CurrPhTr->celly_xboundR[lj], x_maxsxx, xcell, COORD_CELL(ixmax, sx, fs->dsx));  //debugging
              }

            }//end loop over j cell row
            PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
          }
       }  //end loop over nPtr
       ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
     
  PetscFunctionReturn(0);  
}

//---------------------------------------------------------------------------

PetscErrorCode DynamicDike_ReadRestart(DBPropDike *dbdike, DBMat *dbm, JacRes *jr, FB *fb, FILE *fp, PetscBool PrintOutput)  
{
  Dike        *dike;
  PetscInt   nD, numDike;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;


  numDike    = dbdike->numDike; // number of dikes

  // create dike database
  ierr = DBDikeCreate(dbdike, dbm, fb, jr, PETSC_TRUE);   CHKERRQ(ierr);

  for(nD = 0; nD < numDike; nD++)
  {
    dike = jr->dbdike->matDike+nD;
    if (dike->dyndike_start)
     {
      // read mean stress history, 2D array (local vector created with DA_CELL_2D_tave in DBReadDike)
      ierr = VecReadRestart(dike->sxx_eff_ave_hist, fp); CHKERRQ(ierr);
     }
  }
  PetscFunctionReturn(0);
}
//------------------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------

PetscErrorCode DynamicDike_WriteRestart(JacRes *jr, FILE *fp)
{

  Dike        *dike;
  PetscInt   nD, numDike;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  numDike    = jr->dbdike->numDike; // number of dikes

  for(nD = 0; nD < numDike; nD++)
  {
    dike = jr->dbdike->matDike+nD;
    if (dike->dyndike_start)
    {
      // WRITE mean stress history 2D array (local vector created with DA_CELL_2D_tave in DBReadDike)
       ierr = VecWriteRestart(dike->sxx_eff_ave_hist, fp); CHKERRQ(ierr);
    }
  }


  PetscFunctionReturn(0);
}
  
//---------------------------------------------------------------------------

PetscErrorCode DynamicDike_Destroy(JacRes *jr)
{
  
  Dike        *dike;
  PetscInt   nD, numDike, dyndike_on;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  numDike    = jr->dbdike->numDike; // number of dikes
  dyndike_on = 0;

  for(nD = 0; nD < numDike; nD++)
  {
     dike = jr->dbdike->matDike+nD;
     ierr = VecDestroy(&dike->sxx_eff_ave_hist); CHKERRQ(ierr);
     dyndike_on=1;

  }

  if (dyndike_on==1)
  {
    ierr = DMDestroy(&jr->DA_CELL_2D_tave); CHKERRQ(ierr);
    ierr = DMDestroy(&jr->DA_CELL_1D); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}