/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
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
#include "advect.h"

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
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many dikes specified! Max allowed: %lld", (LLD)_max_num_dike_ );

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
			if (numdyndike ==1) //(take this out of this loop because it will be repeated with >1 dynamic dike)
			{
				// DM for 1D cell center vector. vector is ny+1 long, but for whatever reason that must appear in the first,
				// not second entry on line 2 below 
				ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
				fs->dsy.tnods, fs->dsy.nproc, fs->dsz.nproc, 
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
	PetscInt  ID;
	Scaling  *scal;
	
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
	ierr = getIntParam(   fb, _OPTIONAL_, "dyndike_start",&dike->dyndike_start, 1, -1); CHKERRQ(ierr);

  if (dike->dyndike_start)
  {
    dike->Tsol = 1000;
    dike->zmax_magma = -15.0;
    dike->drhomagma = 500;
    dike->filtx = 1.0;
    dike->filty = 4.0;
    dike->istep_nave = 2;

    ierr = getScalarParam(fb, _OPTIONAL_, "Tsol",         &dike->Tsol,         1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "zmax_magma",   &dike->zmax_magma,   1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "filtx",        &dike->filtx,        1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "filty",        &dike->filty,        1, 1.0); CHKERRQ(ierr);
    ierr = getScalarParam(fb, _OPTIONAL_, "drhomagma",    &dike->drhomagma,    1, 1.0); CHKERRQ(ierr);
    ierr = getIntParam(fb, _OPTIONAL_, "istep_nave",     &dike->istep_nave,        1, 50); CHKERRQ(ierr);

    dike->istep_count=dike->istep_nave;   //initialize so that when istep=dike_start, it is set to 0
  }


	// scale the location of Mc y_Mc properly:
	dike->y_Mc /= scal->length;

  
  if (PrintOutput)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Dike parameters ID[%lld]: PhaseTransID=%i PhaseID=%i Mf=%g, Mb=%g, Mc=%g, y_Mc=%g \n", 
      (LLD)(dike->ID), dike->PhaseTransID, dike->PhaseID, dike->Mf, dike->Mb, dike->Mc, dike->y_Mc);
    if (dike->dyndike_start)
    {
      PetscPrintf(PETSC_COMM_WORLD,"                         dyndike_start=%i, Tsol=%g, zmax_magma=%g, drhomagma=%g \n", 
        dike->dyndike_start, dike->Tsol, dike->zmax_magma, dike->drhomagma);
      PetscPrintf(PETSC_COMM_WORLD,"                         filtx=%g, filty=%g, istep_nave=%i, istep_count=%i \n",
        dike->filtx, dike->filty, dike->istep_nave, dike->istep_count);

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
PetscErrorCode Locate_Dike_Zones(AdvCtx *actx)
{
//Ugh:  history arrays are from before the last timestep and so stresses used here are a step behind

  Controls    *ctrl;
  JacRes      *jr;
  Dike        *dike;
  Ph_trans_t  *CurrPhTr;
  FDSTAG      *fs;
  PetscInt   nD, numDike, numPhtr, nPtr, n, icounter;
  PetscInt   j, j1, j2, sx, sy, sz, ny, nx, nz;
  PetscErrorCode ierr; 

  PetscFunctionBeginUser;

  jr = actx->jr;
  ctrl = &jr->ctrl;

  if (!ctrl->actDike || jr->ts->istep+1 == 0) PetscFunctionReturn(0);   // only execute this function if dikes are active
  fs = jr->fs;

  PetscPrintf(PETSC_COMM_WORLD, "\n");

  numDike    = jr->dbdike->numDike; // number of dikes
  numPhtr    = jr->dbm->numPhtr;

  icounter=0;
  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

  for(nD = 0; nD < numDike; nD++)
  {
    dike = jr->dbdike->matDike+nD;

    if (dike->dyndike_start && jr->ts->istep+1 >= dike->dyndike_start)
    {
       // compute lithostatic pressure
       if (icounter==0) 
       {
         ierr = JacResGetLithoStaticPressure(jr); CHKERRQ(ierr);
         ierr = ADVInterpMarkToCell(actx);   CHKERRQ(ierr);
       }
       icounter++;
       //---------------------------------------------------------------------------------------------
       //  Find dike phase transition
       //---------------------------------------------------------------------------------------------
       nPtr=-1;
       for(n=0; n<numPhtr; n++)                          
       {                                               
          CurrPhTr = jr->dbm->matPhtr+n;
          if(CurrPhTr->ID == dike->PhaseTransID)  
          {
            nPtr=n;
          }
       }//end loop over Phtr
       if (nPtr==-1) 
       SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "PhaseTransID problems with dike %i, nPtr=%i\n", nD, nPtr);
       CurrPhTr = jr->dbm->matPhtr+nPtr;
       //---------------------------------------------------------------------------------------------
       //  Find y-bounds of current dynamic dike
       //---------------------------------------------------------------------------------------------
       j1=ny-1;
       j2=0;
       for(j = 0; j < ny; j++)
       {
         if (CurrPhTr->celly_xboundR[j] > CurrPhTr->celly_xboundL[j])
         {
           j1=min(j1,j);
           j2=max(j2,j);
         }
       }

      ierr = Compute_sxx_eff(jr,nD); CHKERRQ(ierr);  //compute mean effective sxx across the lithosphere

      ierr = Smooth_sxx_eff(jr,nD, j1, j2); CHKERRQ(ierr);  //smooth mean effective sxx

      ierr = Set_dike_zones(jr, nD, nPtr,j1, j2); CHKERRQ(ierr); //centered on peak sxx_eff_ave  //commented out for debugging
    }

  }
  PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------------

PetscErrorCode Compute_sxx_eff(JacRes *jr, PetscInt nD)
{
  MPI_Request srequest, rrequest;
  Vec         vsxx, vliththick, vzsol;
  PetscScalar ***gsxx_eff_ave, ***p_lith;
  PetscScalar ***sxx,***liththick, ***zsol;
  PetscScalar  *lsxx, *lliththick, *lzsol;
  PetscScalar dz, ***lT, Tc, *grav, Tsol, dPmag;
  PetscScalar dbug1, dbug2, dbug3, xcell, ycell;
  PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, L, ID, AirPhase;
//  PetscInt    M;
  PetscMPIInt    rank;
  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz;
  //Discret1D   *dsy;
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
  //dsy = &fs->dsy;
  //M   =  (PetscInt)dsy->rank;
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
  // DMDAVecGetArray(DM da,Vec vec,void *array) Returns a multiple dimension array that shares data with the underlying vector 
  // and is indexed using the global dimension
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);

  // open linear buffer for send/receive  (returns the point, lsxx, that contains this processor portion of vector data, vsxx<<G.Ito)
  //Returns a pointer to a contiguous array that contains this processors portion of the vector data.
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
    //if (dPmag<0) dPmag=dPmag*1e3;                                  //Keep dike over the magma. But caution with Smooth_sxx_eff    
    gsxx_eff_ave[L][j][i]=sxx[L][j][i]/liththick[L][j][i]+dPmag;  //Depth weighted mean effective stress: (total stress)+(magma press)                                                                      // (magma press)=lp_lith+dPmagma
  END_PLANE_LOOP

  if (L==0 && dbug3 < 0.05/jr->ts->nstep_out)  //debugging
  {
     START_PLANE_LOOP
        xcell=COORD_CELL(i, sx, fs->dsx);
        ycell=COORD_CELL(j, sy, fs->dsy);
        dPmag=-(zsol[L][j][i]-dike->zmax_magma)*(dike->drhomagma)*grav[2];  //effective pressure is P-Pmagma= negative of magmastatic pressure at solidus, note z AND grav[2] <0
        //if (dPmag<0) dPmag=dPmag*1e3;                                  //Keep dike over the magma. But caution with Smooth_sxx_eff    
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"101010.1010 %i %g %g %g %g %i\n", jr->ts->istep+1, xcell, ycell, gsxx_eff_ave[L][j][i], dPmag, nD);   //debugging    
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

PetscErrorCode Smooth_sxx_eff(JacRes *jr, PetscInt nD, PetscInt  j1, PetscInt j2)
{

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz, *dsy;
  //Discret1D   *dsx;

  //Ph_trans_t  *CurrPhTr;

  PetscScalar ***gsxx_eff_ave, ***gsxx_eff_ave_hist;
  PetscScalar xc, yc, xx, yy, dx, dy, sum_sxx, sum_w;
  PetscInt    j, jj, j1prev, j2prev, j1next, j2next, jj1, jj2; 
  PetscInt    i,ii, ii1, ii2;
  PetscInt    sx, sy, sz, nx, ny, nz;
  PetscInt    L, M;
  PetscMPIInt rank;
  PetscInt    sisc, istep_count, istep_nave;
  Vec         vycoors, vycoors_prev, vycoors_next;
  Vec         vybound, vybound_prev, vybound_next;
  Vec         vsxx, vsxx_prev, vsxx_next;
  PetscScalar ***ycoors, *lycoors, ***ycoors_prev, *lycoors_prev, ***ycoors_next, *lycoors_next;
  PetscScalar ***ybound, *lybound, ***ybound_prev, *lybound_prev, ***ybound_next, *lybound_next;
  PetscScalar ***sxx,***sxx_prev, ***sxx_next;
  PetscScalar  *lsxx, *lsxx_prev, *lsxx_next;
  PetscScalar filtx, filty, w, dfac;
  PetscScalar dbug1, dbug2, dbug3;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Request srequest, rrequest, srequest2, rrequest2, srequest3, rrequest3;

  fs  =  jr->fs;
  dsz = &fs->dsz;
  dsy = &fs->dsy;
  L   =  (PetscInt)dsz->rank;
  M   =  (PetscInt)dsy->rank;
  dbug1=(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out);  //debugging
  dbug2=floor(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out); //debugging
  dbug3=dbug1-dbug2; //debugging

  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

  dike = jr->dbdike->matDike+nD;
  filtx=dike->filtx;
  filty=dike->filty;
  dfac=2.0; //maximum distance for Gaussian weights is dfac*filtx and dfac*filty
//-------------------------------------------------------------------------------------------------------
// get communication buffer (Gets a PETSc vector, vycoors, that may be used with the DM global routines)
//y node coords
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vycoors); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vycoors_prev); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vycoors_next); CHKERRQ(ierr);

  ierr = VecZeroEntries(vycoors); CHKERRQ(ierr);
  ierr = VecZeroEntries(vycoors_prev); CHKERRQ(ierr);
  ierr = VecZeroEntries(vycoors_next); CHKERRQ(ierr);

//celly_xbound info
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vybound); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vybound_prev); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vybound_next); CHKERRQ(ierr);

  ierr = VecZeroEntries(vybound); CHKERRQ(ierr);
  ierr = VecZeroEntries(vybound_prev); CHKERRQ(ierr);
  ierr = VecZeroEntries(vybound_next); CHKERRQ(ierr);

//sxx_ave info
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vsxx_prev); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vsxx_next); CHKERRQ(ierr);

  ierr = VecZeroEntries(vsxx); CHKERRQ(ierr);
  ierr = VecZeroEntries(vsxx_prev); CHKERRQ(ierr);
  ierr = VecZeroEntries(vsxx_next); CHKERRQ(ierr); 


// open index buffer for computation (ycoors is the array that shares data with vector vycoors & indexed with global dimensions)
//y node coords
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vycoors, &ycoors); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vycoors_prev, &ycoors_prev); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vycoors_next, &ycoors_next); CHKERRQ(ierr);
//celly_xbound info
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vybound, &ybound); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vybound_prev, &ybound_prev); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_1D, vybound_next, &ybound_next); CHKERRQ(ierr);
//sxx_ave info
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx_prev, &sxx_prev); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx_next, &sxx_next); CHKERRQ(ierr);


// open linear buffer for send/receive  (returns the pointer, lsxx..., that contains this processor portion of vector data, vycoors)
//y node coords
  ierr = VecGetArray(vycoors, &lycoors); CHKERRQ(ierr);
  ierr = VecGetArray(vycoors_prev, &lycoors_prev); CHKERRQ(ierr);
  ierr = VecGetArray(vycoors_next, &lycoors_next); CHKERRQ(ierr);
//celly_xbound info
  ierr = VecGetArray(vybound, &lybound); CHKERRQ(ierr);
  ierr = VecGetArray(vybound_prev, &lybound_prev); CHKERRQ(ierr);
  ierr = VecGetArray(vybound_next, &lybound_next); CHKERRQ(ierr);
//sxx_ave info
  ierr = VecGetArray(vsxx, &lsxx); CHKERRQ(ierr);
  ierr = VecGetArray(vsxx_prev, &lsxx_prev); CHKERRQ(ierr);
  ierr = VecGetArray(vsxx_next, &lsxx_next); CHKERRQ(ierr);

//-------------------------------------------------------------------------------------------------------
//access depth-averaged effective sxx array
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
  
  for(j = sy; j < sy+ny; j++) 
  { for(i = sx; i < sx+nx; i++) 
    {
      sxx[L][j][i]=gsxx_eff_ave[L][j][i];
    }
  }
  
//---------------------------------------------------------------------------------------------
//  Save info on y bounds on current dike
//---------------------------------------------------------------------------------------------
  for(j = j1; j <=j2; j++)
  {       
    ybound[L][M][j] = j+10;
  }
//---------------------------------------------------------------------------------------------
// passing arrays between previous and next y proc
//---------------------------------------------------------------------------------------------

  for(j = 0; j <= ny; j++)
  {
      ycoors[L][M][j]=COORD_NODE(j+sy,sy,fs->dsy);  //can put j in last entry because ny<nx
  } 

  if (dsy->nproc > 1 && dsy->grprev != -1)  //not the first proc
  {
      ierr = MPI_Irecv(lycoors_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Irecv(lybound_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Irecv(lsxx_prev, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest3); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      ierr = MPI_Isend(lycoors, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Isend(lybound, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest3); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
  }

  if ((dsy->nproc != 1) &&  (dsy->grnext != -1))  //the first proc in y
  {
      ierr = MPI_Isend(lycoors, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Isend(lybound, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest3); CHKERRQ(ierr);
      ierr = MPI_Wait(&srequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

      ierr = MPI_Irecv(lycoors_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Irecv(lybound_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
      ierr = MPI_Irecv(lsxx_next, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest3); CHKERRQ(ierr);
      ierr = MPI_Wait(&rrequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
  }

  //---------------------------------------------------------------------------------
  // Gaussian filter
  //---------------------------------------------------------------------------------

  for(j = j1+sy; j <= j2+sy; j++) //loop over ybounds of current dike
  {
    yc = COORD_CELL(j, sy, fs->dsy);

    
    //identify global y index of cells within dfac*filty of yc
    j1prev=ny+sy-1; j2prev=sy;
    j1next=ny+sy-1; j2next=sy;
    jj1=sy+ny-1; jj2=sy;
    for(jj = sy; jj < sy+ny; jj++)
    {
      
      yy=(ycoors_prev[L][M][jj-sy+1]+ycoors_prev[L][M][jj-sy])/2;
      if ( dsy->grprev != -1 && fabs(yc-yy) <= dfac*filty && ybound_prev[L][M][jj-sy]==jj-sy+10) 
      {
        j1prev=min(j1prev,jj);   
        j2prev=max(j2prev,jj);
      }

      yy=(ycoors_next[L][M][jj-sy+1]+ycoors_next[L][M][jj-sy])/2;
      if (dsy->grnext != -1 && fabs(yy-yc) <= dfac*filty && ybound_next[L][M][jj-sy]==jj-sy+10)
      {
        j1next=min(j1next,jj);   
        j2next=max(j2next,jj);
      }
      
      yy=COORD_CELL(jj, sy, fs->dsy);
      if (fabs(yy-yc) <= dfac*filty && ybound[L][M][jj-sy]==jj-sy+10)
      {
        jj1=min(jj1,jj);
        jj2=max(jj2,jj);
      }
    }


    for (i = sx; i < sx+nx; i++)  //loop to assign filtered value in cell j,i
    {
      sum_sxx=0.0;
      sum_w=0.0;

      xc =  COORD_CELL(i, sx, fs->dsx);
      
      //identify x cells within dfac*filtx of xc
      ii1=sx+nx-1; ii2=sx;
      for (ii = sx; ii < sx+nx; ii++)
      {
        xx = COORD_CELL(ii, sx, fs->dsx);
        if (fabs(xx-xc) <= dfac*filtx)
        {
          ii1=min(ii1,ii);
          ii2=max(ii2,ii);
        }
      }
      //weighted mean of values from previous node
      for (jj = j1prev; jj <= j2prev; jj++) 
      {
         dy=ycoors_prev[L][M][jj+1-sy]-ycoors_prev[L][M][jj-sy];
         yy = (ycoors_prev[L][M][jj+1-sy] + ycoors_prev[L][M][jj-sy])/2;
         for (ii = ii1; ii <= ii2; ii++)
         {
           dx = SIZE_CELL(ii, sx, fs->dsx);
           xx = COORD_CELL(ii, sx, fs->dsx);
           w=exp(-0.5*(pow(((xx-xc)/filtx),2) + pow(((yy-yc)/filty),2)))*dx*dy;
           sum_sxx+=sxx_prev[L][jj][ii]*w;
           sum_w+=w;
         }
      }

      //weighted mean of values on current node
      for (jj = jj1; jj <= jj2; jj++)
      {
         dy=SIZE_CELL(jj,sy,fs->dsy);
         yy = COORD_CELL(jj,sy,fs->dsy);
 
         for (ii = ii1; ii <= ii2; ii++)
         {
           dx = SIZE_CELL(ii,sx,fs->dsx);
           xx = COORD_CELL(ii, sx, fs->dsx);
           w=exp(-0.5*(pow(((xx-xc)/filtx),2) + pow(((yy-yc)/filty),2)) )*dx*dy;
           sum_sxx+=sxx[L][jj][ii]*w;
           sum_w+=w;
         }
      }

      //weighted mean of values from next node
      for (jj = j1next; jj <= j2next; jj++)
      {
         dy=ycoors_next[L][M][jj+1-sy]-ycoors_next[L][M][jj-sy];
         yy = (ycoors_next[L][M][jj+1-sy]-ycoors_next[L][M][jj-sy])/2;

         for (ii = ii1; ii <= ii2; ii++)
         {
           dx = SIZE_CELL(ii,sx,fs->dsx);
           xx = COORD_CELL(ii, sx, fs->dsx);
           w=exp(-0.5*(pow(((xx-xc)/filtx),2) + pow(((yy-yc)/filty),2)))*dx*dy;
           sum_sxx+=sxx_next[L][jj][ii]*w;
           sum_w+=w;
         }
      }

      sum_w=max(sum_w,0.0);
      gsxx_eff_ave[L][j][i]=sum_sxx/sum_w;

    }//End loop over i
  }// End loop over j
  

  //restore ycoors arrays
  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vycoors_prev, &ycoors_prev); CHKERRQ(ierr);
  ierr = VecRestoreArray(vycoors_prev, &lycoors_prev); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vycoors_prev); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vycoors, &ycoors); CHKERRQ(ierr);
  ierr = VecRestoreArray(vycoors, &lycoors); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vycoors); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vycoors_next, &ycoors_next); CHKERRQ(ierr);
  ierr = VecRestoreArray(vycoors_next, &lycoors_next); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vycoors_next); CHKERRQ(ierr);

  //restore ybound arrays
  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vybound_prev, &ybound_prev); CHKERRQ(ierr);
  ierr = VecRestoreArray(vybound_prev, &lybound_prev); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vybound_prev); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vybound, &ybound); CHKERRQ(ierr);
  ierr = VecRestoreArray(vybound, &lybound); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vybound); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vybound_next, &ybound_next); CHKERRQ(ierr);
  ierr = VecRestoreArray(vybound_next, &lybound_next); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vybound_next); CHKERRQ(ierr);

  //restore sxx arrays
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx_prev, &sxx_prev); CHKERRQ(ierr);
  ierr = VecRestoreArray(vsxx_prev, &lsxx_prev); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx_prev); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = VecRestoreArray(vsxx, &lsxx); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx_next, &sxx_next); CHKERRQ(ierr);
  ierr = VecRestoreArray(vsxx_next, &lsxx_next); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx_next); CHKERRQ(ierr);

//---------------------------------------------------------------------------------------------
//  averaging over time
//---------------------------------------------------------------------------------------------
  if (dike->istep_nave>1)
  {
   ierr = DMDAGetCorners(jr->DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, &istep_nave); CHKERRQ(ierr);

   if(istep_nave!=dike->istep_nave) 
   {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Problems: istep_nave=%i, dike->istep_nave=%i\n", istep_nave, dike->istep_nave);
   }

   ierr = DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);

   dike->istep_count++;

   if (dike->istep_count+1 > istep_nave) 
   {
      dike->istep_count=0;
   }  
    
   for (j = j1+sy; j <= j2+sy; j++ )  //Global coordinates
   { 
      for (i=sx; i < sx+ny; i++)
      {
         sum_sxx=0;
         gsxx_eff_ave_hist[sisc+dike->istep_count][j][i]=gsxx_eff_ave[L][j][i];  //array for current step

         for (istep_count=sisc; istep_count<sisc+istep_nave; istep_count++)
         {
            sum_sxx+=gsxx_eff_ave_hist[istep_count][j][i];
         }

         gsxx_eff_ave[L][j][i]=sum_sxx/((PetscScalar)istep_nave);
      }
    }

   ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);
  }// end if nstep_ave>1

  if (L==0 && dbug3 < 0.05/jr->ts->nstep_out)  //debugging
  { 
    START_PLANE_LOOP
      xc=COORD_CELL(i, sx, fs->dsx);
      yc=COORD_CELL(j, sy, fs->dsy);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"202020.2020 %i %g %g %g %i\n", jr->ts->istep+1,xc, yc, gsxx_eff_ave[L][j][i], nD);   //debugging    
    END_PLANE_LOOP  
  }            
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);   //debugging    

  //restore arrays
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);

  PetscFunctionReturn(0);  
}  

//----------------------------------------------------------------------------------------------------
// Set bounds of NotInAir box based on peak sxx_eff_ave
// NOTE that NOW, this only works if cpu_x =1
//

PetscErrorCode Set_dike_zones(JacRes *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2)
{

  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsx, *dsz;
  Ph_trans_t  *CurrPhTr;
  PetscScalar ***gsxx_eff_ave;
  PetscScalar xcenter, sxx_max, dike_width, mindist, xshift, xcell;
  //PetscScalar dx;
  PetscInt    i, lj, j, sx, sy, sz, nx, ny, nz, L, Lx, ixcenter;
  PetscScalar dbug1, dbug2, dbug3, ycell;
  PetscScalar sxxm, sxxp, dx12, dsdx1, dsdx2, x_maxsxx;   
  PetscInt    ixmax;
 
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  fs  =  jr->fs;
  dsz = &fs->dsz;
  dsx = &fs->dsx;
  L   =  dsz->rank;
  Lx  =  dsx->rank;

  dbug1=(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out);  //debugging
  dbug2=floor(((PetscScalar)jr->ts->istep+1)/jr->ts->nstep_out); //debugging
  dbug3=dbug1-dbug2; //debugging

  dike = jr->dbdike->matDike+nD;
  CurrPhTr = jr->dbm->matPhtr+nPtr;

  if (Lx>0)
  {
     PetscPrintf(PETSC_COMM_WORLD,"Set_dike_zones requires cpu_x = 1 Lx = %i \n", Lx);
     SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Set_dike_zones requires cpu_x = 1 Lx = %i \n", Lx);
  }
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
                                        
  for(lj = j1; lj <= j2; lj++)  //local index
  {
     sxx_max=-1e12;
     mindist=1e12;
     ixcenter = 0;
     xshift=0;
     ixmax=sx+1;

     j=sy+lj;  //global index
     dike_width=CurrPhTr->celly_xboundR[lj]-CurrPhTr->celly_xboundL[lj];
     xcenter=(CurrPhTr->celly_xboundR[lj] + CurrPhTr->celly_xboundL[lj])/2;

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


     //dx=SIZE_CELL(ixcenter,sx, fs->dsx);
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
        ycell = COORD_CELL(j, sy, fs->dsy);  //debugging
        xcell=(COORD_CELL(ixmax-1, sx, fs->dsx)+COORD_CELL(ixmax, sx, fs->dsx))/2;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"303030.3030 %i %g %g %g %g %g %g %g\n", jr->ts->istep+1, ycell, xcenter+xshift, 
        CurrPhTr->celly_xboundL[lj], CurrPhTr->celly_xboundR[lj], x_maxsxx, xcell, COORD_CELL(ixmax, sx, fs->dsx),nD);  //debugging
     }

  }//end loop over j cell row

  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
     
  PetscFunctionReturn(0);  
}

//---------------------------------------------------------------------------
PetscErrorCode DynamicDike_ReadRestart(DBPropDike *dbdike, DBMat *dbm, JacRes *jr, FB *fb, FILE *fp)
//PetscErrorCode DynamicDike_ReadRestart(DBPropDike *dbdike, DBMat *dbm, JacRes *jr, FB *fb, FILE *fp, PetscBool PrintOutput)
{
  Controls    *ctrl;
  Dike        *dike;
  PetscInt   nD, numDike;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ctrl = &jr->ctrl;
  if (!ctrl->actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active

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

  Controls    *ctrl;
  Dike        *dike;
  PetscInt   nD, numDike;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ctrl = &jr->ctrl;
  if (!ctrl->actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active

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
  Controls    *ctrl;
  PetscInt   nD, numDike, dyndike_on;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ctrl = &jr->ctrl;

  if (!ctrl->actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active


  numDike    = jr->dbdike->numDike; // number of dikes
  dyndike_on=0;
  for(nD = 0; nD < numDike; nD++)
  {
     dike = jr->dbdike->matDike+nD;
     if (dike->dyndike_start)
     {
       ierr = VecDestroy(&dike->sxx_eff_ave_hist); CHKERRQ(ierr);
       dyndike_on=1;
     }
  }

  if (dyndike_on==1)
  {
    ierr = DMDestroy(&jr->DA_CELL_2D_tave); CHKERRQ(ierr);
    ierr = DMDestroy(&jr->DA_CELL_1D); CHKERRQ(ierr);    
  }

  PetscFunctionReturn(0);
}
