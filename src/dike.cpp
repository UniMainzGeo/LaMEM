
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
			ierr = DMCreateLocalVector( jr->DA_CELL_2D, &dike->magPressure);  CHKERRQ(ierr);
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
		dike->filtx = 1.5;
		dike->filty = 1.5;
		dike->istep_nave = 2;
		dike->nstep_locate = 1;
		dike->out_stress = 0;
		dike->out_dikeloc = 0;
		dike->magPfac=1.0;
		dike->magPwidth=1e+30; 
		//dike->ymindyn=-1e+30;
		//dike->ymaxdyn=1e+30; 

		ierr = getScalarParam(fb, _OPTIONAL_, "Tsol",		&dike->Tsol,		1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "zmax_magma",	&dike->zmax_magma,	1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "filtx",		&dike->filtx,		1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "filty",		&dike->filty,		1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "drhomagma",	&dike->drhomagma,	1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "magPfac",	&dike->magPfac,		1, 1.0); CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "magPwidth",	&dike->magPwidth,	1, 1.0); CHKERRQ(ierr);
		//ierr = getScalarParam(fb, _OPTIONAL_, "ymindyn",	&dike->ymindyn,		1, 1.0); CHKERRQ(ierr);
		//ierr = getScalarParam(fb, _OPTIONAL_, "ymaxdyn",	&dike->ymaxdyn,		1, 1.0); CHKERRQ(ierr);

		ierr = getIntParam(fb, _OPTIONAL_, "istep_nave",	&dike->istep_nave,	1, 50); CHKERRQ(ierr);
		ierr = getIntParam(fb, _OPTIONAL_, "nstep_locate",	&dike->nstep_locate,1, 1000); CHKERRQ(ierr);
		ierr = getIntParam(fb, _OPTIONAL_, "out_stress",	&dike->out_stress,	1, 50); CHKERRQ(ierr);
		ierr = getIntParam(fb, _OPTIONAL_, "out_dikeloc",	&dike->out_dikeloc,	1, 50); CHKERRQ(ierr);

		dike->istep_count=dike->istep_nave;   //initialize so that when istep=dike_start, it is set to 0
	}


	// scale the location of Mc y_Mc properly:
	dike->y_Mc /= scal->length;

  
	if (PrintOutput)
	{
		PetscPrintf(PETSC_COMM_WORLD,"  Dike parameters ID[%lld]: PhaseTransID=%lld PhaseID=%lld Mf=%g, Mb=%g, Mc=%g, y_Mc=%g \n", 
		(LLD)(dike->ID), (LLD)(dike->PhaseTransID), (LLD)(dike->PhaseID), dike->Mf, dike->Mb, dike->Mc, dike->y_Mc);
	if (dike->dyndike_start)
	{
		PetscPrintf(PETSC_COMM_WORLD,"       dyndike_start=%lld, Tsol=%g, zmax_magma=%g,drhomagma=%g, magPfac=%g, magPwidth=%g\n", 
			(LLD)(dike->dyndike_start), dike->Tsol, dike->zmax_magma, dike->drhomagma, dike->magPfac, dike->magPwidth);
		PetscPrintf(PETSC_COMM_WORLD,"       filtx=%g, filty=%g, istep_nave=%lld, istep_count=%lld \n",
			dike->filtx, dike->filty, (LLD)(dike->istep_nave), (LLD)(dike->istep_count));
		PetscPrintf(PETSC_COMM_WORLD,"       nstep_locate=%lld, out_stress=%lld, out_dikeloc=%lld\n",
			(LLD)(dike->nstep_locate), (LLD)(dike->out_stress), (LLD)(dike->out_dikeloc));
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
  fs  =  jr->fs;
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

	if (dike->dyndike_start && (jr->ts->istep+1 >= dike->dyndike_start) && ((jr->ts->istep+1) % dike->nstep_locate) == 0) 
	//if (dike->dyndike_start && (jr->ts->istep+1 >= dike->dyndike_start)) //debugging
    {
	   PetscPrintf(PETSC_COMM_WORLD, "Locating Dike zone: istep=%lld dike # %lld\n", (LLD)(jr->ts->istep + 1),(LLD)(nD));
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
       SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "PhaseTransID problems with dike %lld, nPtr=%lld\n", (LLD)(nD), (LLD)(nPtr));
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

      ierr = Compute_sxx_magP(jr, nD); CHKERRQ(ierr);  //compute mean effective sxx across the lithosphere

      ierr = Smooth_sxx_eff(jr,nD, nPtr, j1, j2); CHKERRQ(ierr);  //smooth mean effective sxx

      ierr = Set_dike_zones(jr, nD, nPtr,j1, j2); CHKERRQ(ierr); //centered on peak sxx_eff_ave 
    }

  }
  PetscFunctionReturn(0);

}
//------------------------------------------------------------------------------------------------------------------

PetscErrorCode Compute_sxx_magP(JacRes *jr, PetscInt nD)
{
  MPI_Request srequest, rrequest;
  Vec         vsxx, vPmag, vliththick, vzsol;
  PetscScalar ***gsxx_eff_ave, ***p_lith, ***magPressure;
  PetscScalar ***sxx, ***Pmag, ***liththick, ***zsol;
  PetscScalar  *lsxx, *lPmag, *lliththick, *lzsol;
  PetscScalar dz, ***lT, Tc, *grav, Tsol, dPmag, magma_presence;
  PetscScalar xcell, ycell;
  PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, L, ID, AirPhase;
  PetscInt 	  istep, nstep_out;
  PetscMPIInt    rank;


  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz;
  SolVarCell  *svCell;
  Controls    *ctrl;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  istep=jr->ts->istep+1; 
  nstep_out=jr->ts->nstep_out;

  ctrl = &jr->ctrl;
  grav = ctrl->grav;

  fs  =  jr->fs;
  dsz = &fs->dsz;
  L   =  (PetscInt)dsz->rank;
  AirPhase  = jr->surf->AirPhase;


  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);

  dike = jr->dbdike->matDike+nD;

  // get local grid sizes
  ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

  // get integration/communication buffer (Gets a PETSc vector, vsxx, that may be used with the DM global routines)
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vPmag); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vliththick); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vzsol); CHKERRQ(ierr);

  ierr = VecZeroEntries(vsxx); CHKERRQ(ierr);
  ierr = VecZeroEntries(vPmag); CHKERRQ(ierr);
  ierr = VecZeroEntries(vliththick); CHKERRQ(ierr);
  ierr = VecZeroEntries(vzsol); CHKERRQ(ierr);

  // open index buffer for computation (sxx the array that shares data with vector vsxx and is indexed with global dimensions<<G.Ito)
  // DMDAVecGetArray(DM da,Vec vec,void *array) Returns a multiple dimension array that shares data with the underlying vector 
  // and is indexed using the global dimension
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vPmag, &Pmag); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);

  // open linear buffer for send/receive  (returns the point, lsxx, that contains this processor portion of vector data, vsxx<<G.Ito)
  //Returns a pointer to a contiguous array that contains this processors portion of the vector data.
  ierr = VecGetArray(vsxx, &lsxx); CHKERRQ(ierr);
  ierr = VecGetArray(vPmag, &lPmag); CHKERRQ(ierr);
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

     ierr = MPI_Irecv(lPmag, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
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
            sxx[L][j][i]+=(svCell->hxx - svCell->svBulk.pn)*dz;  //integrating dz-weighted total stress
            Pmag[L][j][i]+=p_lith[k][j][i]*dz;

            liththick[L][j][i]+=dz;             //integrating thickeness
          }
          
          //interpolate depth to the solidus
          if ((Tc <= Tsol) & (Tsol < lT[k-1][j][i]))
          {
            zsol[L][j][i]=dsz->ccoor[k-sz]+(dsz->ccoor[k-sz-1]-dsz->ccoor[k-sz])/(lT[k-1][j][i]-Tc)*(Tsol-Tc); 
          }
      END_PLANE_LOOP
  } 

      //After integrating and averaging, send it down to the next proc. 
  if(dsz->nproc != 1 && dsz->grprev != -1)
  {
     ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Isend(lPmag, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Isend(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Isend(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
  }
      
	//Now receive/send the answer from successive previous (underlying) procs so all procs have the answers
  if(dsz->nproc != 1 && dsz->grprev != -1)
  {
     ierr = MPI_Irecv(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Irecv(lPmag, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
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

     ierr = MPI_Isend(lPmag, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Isend(lliththick, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

     ierr = MPI_Isend(lzsol, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

  }

  // (gdev is the array that shares data with devxx_mean and is indexed with global dimensions)
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);

  //now all cores in z have the same solution so give that to the stress array

/*	START_PLANE_LOOP  //testing this option
		zsol_max=max(zsol[L][j][i],zsol_max);     //Need to find global max accross procs via MPI_Reduce                                                    
	END_PLANE_LOOP
*/
	magma_presence=1.0;  //testing
	
	START_PLANE_LOOP
		dPmag=(dike->zmax_magma-zsol[L][j][i])*(dike->drhomagma)*grav[2];  //excess static magma pressure at solidus, z & grav[2] <0 so this is positive
		//magma_presence=dike->magPfac*(zsol[L][j][i]-dike->zmax_magma)/(zsol_max-dike->zmax_magma);  //this feature undergoing testing
		magPressure[L][j][i] = (Pmag[L][j][i]/liththick[L][j][i]+dPmag)*magma_presence;
		gsxx_eff_ave[L][j][i]= sxx[L][j][i]/liththick[L][j][i];  //Depth weighted mean total stress                                                                
	END_PLANE_LOOP

	//output means stress array to std on timesteps of other output
	if (((istep % nstep_out)==0) && (dike->out_stress>0))  
	{
		if (L==0)
		{
			START_PLANE_LOOP
				xcell=COORD_CELL(i, sx, fs->dsx);
				ycell=COORD_CELL(j, sy, fs->dsy);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"101010.1010 %lld %g %g %g %g %lld %g %g\n", (LLD)(jr->ts->istep+1), 
				xcell, ycell, gsxx_eff_ave[L][j][i], magPressure[L][j][i], (LLD)(nD), zsol[L][j][i],magma_presence);     
			END_PLANE_LOOP 
		} 
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT); //All procs must run this
	}

  // restore buffer and mean stress vectors
  ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);      
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);


  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx, &Pmag); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);

  ierr = VecRestoreArray(vsxx, &lsxx); CHKERRQ(ierr);
  ierr = VecRestoreArray(vsxx, &lPmag); CHKERRQ(ierr);
  ierr = VecRestoreArray(vliththick, &lliththick); CHKERRQ(ierr);
  ierr = VecRestoreArray(vzsol, &lzsol); CHKERRQ(ierr);

  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vPmag); CHKERRQ(ierr);  
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vliththick); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vzsol); CHKERRQ(ierr);

  //fill ghost points
      
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->sxx_eff_ave);

  ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &p_lith); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------
// Apply elliptical Gaussian smoothing of depth-averaged effective stress.  
// **NOTE** There is NO message passing between adjacent procs in x

PetscErrorCode Smooth_sxx_eff(JacRes *jr, PetscInt nD, PetscInt nPtr, PetscInt  j1, PetscInt j2)
{

	FDSTAG      *fs;
	Dike        *dike;
	Discret1D   *dsz, *dsy;
	Ph_trans_t  *CurrPhTr;

	PetscScalar ***gsxx_eff_ave, ***gsxx_eff_ave_hist, ***magPressure;
	PetscScalar ***ycoors, *lycoors, ***ycoors_prev, *lycoors_prev, ***ycoors_next, *lycoors_next;
	PetscScalar ***ybound, *lybound, ***ybound_prev, *lybound_prev, ***ybound_next, *lybound_next;
	PetscScalar ***sxx, *lsxx, ***sxx_prev, *lsxx_prev, ***sxx_next, *lsxx_next;
	PetscScalar ***magP, *lmagP, ***magP_prev, *lmagP_prev, ***magP_next, *lmagP_next;
	PetscScalar xc, yc, xx, yy, dx, dy, sum_sxx, sum_magP, sum_w;
	PetscScalar filtx, filty, filtxy, w, dfac, magPfac, magPwidth,xcenter;
	PetscScalar xcenter_north, xcenter_south, ycenter_north, ycenter_south, xcenter_search, ycenter_search;
	PetscScalar azim, dalong, dxazim, dyazim, radbound, sumslope, sumadd;

	Vec         vycoors, vycoors_prev, vycoors_next;
	Vec         vybound, vybound_prev, vybound_next;
	Vec         vsxx, vsxx_prev, vsxx_next;
	Vec         vmagP, vmagP_prev, vmagP_next;

	PetscInt    j, jj, j1prev, j2prev, j1next, j2next, jj1, jj2, jc; 
	PetscInt    i,ii, ii1, ii2;
	PetscInt    sx, sy, sz, nx, ny, nz;
	PetscInt    L, M;
	PetscMPIInt rank;
	PetscInt    sisc, istep_count, istep_nave, istep, nstep_out;

  
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Request srequest, rrequest, srequest2, rrequest2, srequest3, rrequest3, srequest4, rrequest4;

	fs  =  jr->fs;
	dsz = &fs->dsz;
	dsy = &fs->dsy;
	L   =  (PetscInt)dsz->rank;
	M   =  (PetscInt)dsy->rank;

	istep=jr->ts->istep+1; 
	nstep_out=jr->ts->nstep_out;

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	//CurrPhTr = jr->dbm->matPhtr+nPtr;
	dike = jr->dbdike->matDike+nD;
	filtx=dike->filtx;
	filty=dike->filty;
	filtxy=max(filtx,filty);
	dfac=2.0; //maximum distance for Gaussian weights is dfac*filtx and dfac*filty

	magPfac=dike->magPfac;
	magPwidth=dike->magPwidth;
	CurrPhTr = jr->dbm->matPhtr+nPtr;

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

//magP info
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vmagP); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vmagP_prev); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_2D, &vmagP_next); CHKERRQ(ierr);

	ierr = VecZeroEntries(vmagP); CHKERRQ(ierr);
	ierr = VecZeroEntries(vmagP_prev); CHKERRQ(ierr);
	ierr = VecZeroEntries(vmagP_next); CHKERRQ(ierr); 

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

//magP info
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, vmagP, &magP); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, vmagP_prev, &magP_prev); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, vmagP_next, &magP_next); CHKERRQ(ierr);


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

//magP info
	ierr = VecGetArray(vmagP, &lmagP); CHKERRQ(ierr);
	ierr = VecGetArray(vmagP_prev, &lmagP_prev); CHKERRQ(ierr);
	ierr = VecGetArray(vmagP_next, &lmagP_next); CHKERRQ(ierr);

//access depth-averaged arrays on current proc
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);

  
	START_PLANE_LOOP

		sxx[L][j][i]=gsxx_eff_ave[L][j][i];
		magP[L][j][i]=magPressure[L][j][i];

	END_PLANE_LOOP
  
//  Save info on y bounds of current dike and y coords of nodes
	for(j = j1; j <=j2; j++)
	{       
		ybound[L][M][j] = (PetscScalar)(j+10); //non-zero values are along the dike. Zero otherwise. 
	}

	for(j = 0; j <= ny; j++)
	{
		ycoors[L][M][j]=COORD_NODE(j+sy,sy,fs->dsy);  //can put j in last entry because ny<nx
	} 

//--------------------------------------------------
// passing arrays between previous and next y proc
//--------------------------------------------------
	if (dsy->nproc > 1 && dsy->grprev != -1)  //Exchange arrays from previous proc, if not the first proc
	{
		ierr = MPI_Irecv(lycoors_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lybound_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lsxx_prev, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lmagP_prev, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);


		ierr = MPI_Isend(lycoors, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lybound, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lmagP, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	}

	if ((dsy->nproc != 1) &&  (dsy->grnext != -1))  //Exhange arrays with next proc, if not the last proc
	{
		ierr = MPI_Isend(lycoors, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lybound, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lmagP, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Irecv(lycoors_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lybound_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lsxx_next, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lmagP_next, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

	}

//--------------------------------------------------------------------------------------
// Gaussian filter with one elliptical axis oriented with local azimuth of dike zone
//--------------------------------------------------------------------------------------
	//loop over ybounds of current dike
	for(j = j1+sy; j <= j2+sy; j++) 
	{

		//Local azimuth of dike as the mean of all dike points within distance 
		//of filty of current dike point (xcenter, yc)
		xcenter=(CurrPhTr->celly_xboundR[j-sy] + CurrPhTr->celly_xboundL[j-sy])/2;
		yc = COORD_CELL(j, sy, fs->dsy);
		sumslope=0;
		sumadd = 0; 
		for (jc=j1+sy; jc<=j2+sy; jc++)
		{
			xcenter_search=(CurrPhTr->celly_xboundR[jc-sy] + CurrPhTr->celly_xboundL[jc-sy])/2;
			ycenter_search=COORD_CELL(jc, sy, fs->dsy);
			dalong=sqrt(pow((xcenter-xcenter_search),2)+pow((yc-ycenter_search),2));

			if (jc<j && dalong<=filty) //if south of current point
			{
				xcenter_north=(CurrPhTr->celly_xboundR[jc-sy+1] + CurrPhTr->celly_xboundL[jc-sy+1])/2;
				xcenter_south=xcenter_search;
				ycenter_north=COORD_CELL(jc+1, sy, fs->dsy);
				ycenter_south=ycenter_search;
				sumslope += (xcenter_north-xcenter_south)/(ycenter_north-ycenter_south);
				sumadd += 1;
			}
			else if (jc >j && dalong<=filty)  //if north of current point
			{
				xcenter_north=xcenter_search;
				xcenter_south=(CurrPhTr->celly_xboundR[jc-sy-1] + CurrPhTr->celly_xboundL[jc-sy-1])/2;
				ycenter_north=ycenter_search;
				ycenter_south=COORD_CELL(jc-1, sy, fs->dsy);
				sumslope += (xcenter_north-xcenter_south)/(ycenter_north-ycenter_south);
				sumadd += 1;
			}
		}
		azim=atan(sumslope/sumadd);

		//identify global y index of cells within dfac*filtxy of yc on local and adjacent processors
		j1prev=ny+sy-1; j2prev=sy;
		j1next=ny+sy-1; j2next=sy;
		jj1=sy+ny-1; jj2=sy;
		for(jj = sy; jj < sy+ny; jj++)
		{
 			yy=(ycoors_prev[L][M][jj-sy+1]+ycoors_prev[L][M][jj-sy])/2;
			if ( dsy->grprev != -1 && fabs(yc-yy) <= dfac*filtxy && ybound_prev[L][M][jj-sy]==(PetscScalar)(jj-sy+10)) 
			{
				j1prev=min(j1prev,jj);   
				j2prev=max(j2prev,jj);
			}

			yy=(ycoors_next[L][M][jj-sy+1]+ycoors_next[L][M][jj-sy])/2;
			if (dsy->grnext != -1 && fabs(yy-yc) <= dfac*filtxy && ybound_next[L][M][jj-sy]==(PetscScalar)(jj-sy+10))
			{
				j1next=min(j1next,jj);   
				j2next=max(j2next,jj);
			}
      
			yy=COORD_CELL(jj, sy, fs->dsy);
			if (fabs(yy-yc) <= dfac*filtxy && ybound[L][M][jj-sy]==(PetscScalar)(jj-sy+10))
			{
				jj1=min(jj1,jj);
				jj2=max(jj2,jj);
			}
		}

		//Loop over i to assign filtered value in cell j,i (again, one proc across all x dimension)
		for (i = sx; i < sx+nx; i++)  
		{
			sum_sxx=0.0;
			sum_magP=0.0;
			sum_w=0.0;

			xc =  COORD_CELL(i, sx, fs->dsx);
      
			//identify x cells within dfac*filtxy of xc
			ii1=sx+nx-1; ii2=sx;
			for (ii = sx; ii < sx+nx; ii++)
			{
				xx = COORD_CELL(ii, sx, fs->dsx);
				if (fabs(xx-xc) <= dfac*filtxy)
				{
					ii1=min(ii1,ii);
					ii2=max(ii2,ii);
				}
			}

			//weighted mean of values from previous proc
			for (jj = j1prev; jj <= j2prev; jj++) 
			{
				dy=ycoors_prev[L][M][jj+1-sy]-ycoors_prev[L][M][jj-sy];
				yy = (ycoors_prev[L][M][jj+1-sy] + ycoors_prev[L][M][jj-sy])/2;
				for (ii = ii1; ii <= ii2; ii++)
				{
					dx = SIZE_CELL(ii, sx, fs->dsx);
					xx = COORD_CELL(ii, sx, fs->dsx);

					dxazim=cos(azim)*(xx-xc)-sin(azim)*(yy-yc);
					dyazim=sin(azim)*(xx-xc)+cos(azim)*(yy-yc);

					radbound=(pow((dxazim/(dfac*filtx)),2) + pow((dyazim/(dfac*filty)),2));					
					if (radbound<=1) //limit area of summing to within radbound of cell
					{
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/filty),2)))*dx*dy;
						sum_sxx += sxx_prev[L][jj][ii]*w;
						sum_magP += magP_prev[L][jj][ii]*w;
						sum_w+=w;
					}
				}
			}//end loop over cells on previous proc

			//weighted mean of values on current proc
			for (jj = jj1; jj <= jj2; jj++)
			{
				dy=SIZE_CELL(jj,sy,fs->dsy);
				yy = COORD_CELL(jj,sy,fs->dsy);
 
				for (ii = ii1; ii <= ii2; ii++)
				{
					dx = SIZE_CELL(ii,sx,fs->dsx);
					xx = COORD_CELL(ii, sx, fs->dsx);

					dxazim=cos(azim)*(xx-xc)-sin(azim)*(yy-yc);
					dyazim=sin(azim)*(xx-xc)+cos(azim)*(yy-yc);

					radbound=(pow((dxazim/(dfac*filtx)),2) + pow((dyazim/(dfac*filty)),2));					
					if (radbound<=1)  //limit area of summing to within radbound of cell
					{
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/filty),2)))*dx*dy;
						sum_sxx += sxx[L][jj][ii]*w;
						sum_magP += magP[L][jj][ii]*w;
						sum_w+=w;
					}				
				}
			}//end loop over cells on current proc

			//weighted mean of values from next proc
			for (jj = j1next; jj <= j2next; jj++)
			{
				dy=ycoors_next[L][M][jj+1-sy]-ycoors_next[L][M][jj-sy];
				yy = (ycoors_next[L][M][jj+1-sy] + ycoors_next[L][M][jj-sy])/2;

				for (ii = ii1; ii <= ii2; ii++)
				{
					dx = SIZE_CELL(ii,sx,fs->dsx);
					xx = COORD_CELL(ii, sx, fs->dsx);

					dxazim=cos(azim)*(xx-xc)-sin(azim)*(yy-yc);
					dyazim=sin(azim)*(xx-xc)+cos(azim)*(yy-yc);

					radbound=(pow((dxazim/(dfac*filtx)),2) + pow((dyazim/(dfac*filty)),2));					
					if (radbound<=1)  //limit area of summing to within radbound of cell
					{
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/filty),2)))*dx*dy;
						sum_sxx += sxx_next[L][jj][ii]*w;
						sum_magP += magP_next[L][jj][ii]*w;
						sum_w+=w;
					}
				}
			} //end loop over cells from next proc

			sum_w=max(sum_w,0.0);
			magPressure[L][j][i]=(sum_magP/sum_w);
			gsxx_eff_ave[L][j][i]=(sum_sxx/sum_w) + magPressure[L][j][i]*magPfac*exp(-0.5*(pow((cos(azim)*(xcenter-xc)/magPwidth),2)));
			//gsxx_eff_ave[L][j][i]=(sum_sxx/sum_w) + magPressure[L][j][i];

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

	//restore sxx arrays
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP_prev, &magP_prev); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP_prev, &lmagP_prev); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP_prev); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP, &magP); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP, &lmagP); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP_next, &magP_next); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP_next, &lmagP_next); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP_next); CHKERRQ(ierr);

//--------------------------------------------------
//  TIME Averaging
//--------------------------------------------------
	if (dike->istep_nave>1)
	{
		ierr = DMDAGetCorners(jr->DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, &istep_nave); CHKERRQ(ierr);

		if(istep_nave!=dike->istep_nave) 
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Problems: istep_nave=%lld, dike->istep_nave=%lld\n", 
			(LLD)(istep_nave), (LLD)(dike->istep_nave));
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

	if (((istep % nstep_out)==0) && (dike->out_stress>0))  
	{
		if (L==0)
		{ 
			START_PLANE_LOOP
				xc=COORD_CELL(i, sx, fs->dsx);
				yc=COORD_CELL(j, sy, fs->dsy);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"202020.2020 %lld %g %g %g %g %lld %lld\n", (LLD)(jr->ts->istep+1),xc, yc, gsxx_eff_ave[L][j][i], 
				magPressure[L][j][i],(LLD)(nD), (LLD)(dike->istep_count));       
			END_PLANE_LOOP  
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT); //All procs must run this
	}            

	//restore arrays
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);

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
	Discret1D   *dsx, *dsy, *dsz;
	Ph_trans_t  *CurrPhTr;
	PetscScalar ***gsxx_eff_ave;
	PetscScalar xcenter, sxx_max, dike_width, mindist, xshift, xcell;
	PetscScalar ***xboundL_pass, *lxboundL_pass, ***xboundR_pass, *lxboundR_pass;
	Vec         vxboundL_pass, vxboundR_pass;
	PetscInt    i, lj, j, sx, sy, sz, nx, ny, nz, L, Lx, M, ixcenter;
	PetscScalar sxxm, sxxp, dx12, dsdx1, dsdx2, x_maxsxx, ycell, dtime;   
	PetscInt    ixmax, istep, nstep_out;
 	MPI_Request srequest, rrequest;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  =  jr->fs;
	dsz = &fs->dsz;
	dsy = &fs->dsy;
	dsx = &fs->dsx;
	L   =  (PetscInt)dsz->rank;
	M   =  (PetscInt)dsy->rank;
	Lx  =  (PetscInt)dsx->rank;

	istep=jr->ts->istep+1; 
	nstep_out=jr->ts->nstep_out;

	dike = jr->dbdike->matDike+nD;
	CurrPhTr = jr->dbm->matPhtr+nPtr;
	dtime=jr->scal->time*jr->ts->time;


	if (Lx>0)
	{
		PetscPrintf(PETSC_COMM_WORLD,"Set_dike_zones requires cpu_x = 1 Lx = %lld \n", (LLD)(Lx));
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Set_dike_zones requires cpu_x = 1 Lx = %lld \n", (LLD)(Lx));
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
 
		dsdx1=(sxx_max-sxxm)/(COORD_CELL(ixmax, sx, fs->dsx)-COORD_CELL(ixmax-1, sx, fs->dsx));  //slope left of max
		dsdx2=(sxxp-sxx_max)/(COORD_CELL(ixmax+1, sx, fs->dsx)-COORD_CELL(ixmax, sx, fs->dsx));  //slope right of max
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

		if (xshift>0 && fabs(xshift) > 0.5*SIZE_CELL(ixcenter, sx, fs->dsx)) //ensure new center is within width of cell to right of center
		{
        	xshift=0.5*SIZE_CELL(ixcenter, sx, fs->dsx);
		}
		else if (xshift<0 && fabs(xshift) > 0.5*SIZE_CELL(ixcenter-1, sx, fs->dsx)) //ensure its within the width of cell left of center
		{
        	xshift=-0.5*SIZE_CELL(ixcenter-1, sx, fs->dsx);
		}

		//relocating dike bounds here
		CurrPhTr->celly_xboundL[lj]=xcenter+xshift-dike_width/2; 
		CurrPhTr->celly_xboundR[lj]=xcenter+xshift+dike_width/2; 

		if (L==0 &&  ((istep % nstep_out)==0) && (dike->out_dikeloc > 0)) 
		{
        	ycell = COORD_CELL(j, sy, fs->dsy);  
        	xcell=(COORD_CELL(ixmax-1, sx, fs->dsx)+COORD_CELL(ixmax, sx, fs->dsx))/2;
        	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"303030.3030 %lld %g %g %g %g %g %g %g %lld %g \n", 
				(LLD)(jr->ts->istep+1), ycell, xcenter, xshift, 
        		x_maxsxx, COORD_CELL(ixmax, sx, fs->dsx), 
				CurrPhTr->celly_xboundL[lj], CurrPhTr->celly_xboundR[lj], (LLD)(nD), dtime);  
		}

	}//end loop over j cell row

	if (((istep % nstep_out)==0) && (dike->out_dikeloc > 0))  
	{
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	}
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);

//-----------------------------------------------------------------------------------
// Set locations of ghost nodes
//-----------------------------------------------------------------------------------
	ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vxboundL_pass); CHKERRQ(ierr);
	ierr = VecZeroEntries(vxboundL_pass); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_1D, vxboundL_pass, &xboundL_pass); CHKERRQ(ierr);
	ierr = VecGetArray(vxboundL_pass, &lxboundL_pass); CHKERRQ(ierr);

	ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vxboundR_pass); CHKERRQ(ierr);
	ierr = VecZeroEntries(vxboundR_pass); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_1D, vxboundR_pass, &xboundR_pass); CHKERRQ(ierr);
	ierr = VecGetArray(vxboundR_pass, &lxboundR_pass); CHKERRQ(ierr);

	//Northernmost (top) ghost coord of northernmost (top) proc
	if (dsy->grnext == -1) 
	{
		CurrPhTr->celly_xboundL[ny] = CurrPhTr->celly_xboundL[ny-1];
		CurrPhTr->celly_xboundR[ny] = CurrPhTr->celly_xboundR[ny-1];
	}
	//Southernmost ghost coord of southernmost (bottom) proc
	if (dsy->grprev == -1)  
	{
		CurrPhTr->celly_xboundL[-1] = CurrPhTr->celly_xboundL[0];
		CurrPhTr->celly_xboundR[-1] = CurrPhTr->celly_xboundR[0];
	}

	//Receive from 2nd northernmost (top in y) proc to southernmost (bottom in y) and set bottom ghost node
	if (dsy->nproc > 1 && dsy->grnext != -1) //MPI_Wait will make this run in sequence from top to bottom
	{
		ierr = MPI_Irecv(lxboundL_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Irecv(lxboundR_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		CurrPhTr->celly_xboundL[ny] = xboundL_pass[L][M][0];
		CurrPhTr->celly_xboundR[ny] = xboundR_pass[L][M][0];
	}

	//Send down from northernmost (top) to southmost (bottom)
	if(dsy->nproc != 1 && dsy->grprev != -1)
  	{
		for(lj = 0; lj < ny; lj++)
		{       
			xboundL_pass[L][M][lj] = CurrPhTr->celly_xboundL[lj];  
			xboundR_pass[L][M][lj] = CurrPhTr->celly_xboundR[lj];  
		}
		ierr = MPI_Isend(lxboundL_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Isend(lxboundR_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
  	}

	if(dsy->nproc != 1 && dsy->grprev != -1)  //Receive coordinates from previous node & set BOTTOM ghost node
	{
		ierr = MPI_Irecv(lxboundL_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Irecv(lxboundR_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
	
		CurrPhTr->celly_xboundL[-1] = xboundL_pass[L][M][ny-1];
		CurrPhTr->celly_xboundR[-1] = xboundR_pass[L][M][ny-1];
  	}

	if(dsy->nproc != 1 && dsy->grnext != -1)
  	{
		for(lj = 0; lj < ny; lj++)
		{       
			xboundL_pass[L][M][lj] = CurrPhTr->celly_xboundL[lj];  
			xboundR_pass[L][M][lj] = CurrPhTr->celly_xboundR[lj];  
		}
     	ierr = MPI_Isend(lxboundL_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     	ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Isend(lxboundR_pass, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
     	ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
  	}

	ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vxboundL_pass, &xboundL_pass); CHKERRQ(ierr);
	ierr = VecRestoreArray(vxboundL_pass, &lxboundL_pass); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vxboundL_pass); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vxboundR_pass, &xboundR_pass); CHKERRQ(ierr);
	ierr = VecRestoreArray(vxboundR_pass, &lxboundR_pass); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vxboundR_pass); CHKERRQ(ierr);

  PetscFunctionReturn(0);  
}

//---------------------------------------------------------------------------
PetscErrorCode DynamicDike_ReadRestart(DBPropDike *dbdike,  DBMat *dbm, JacRes *jr, TSSol *ts, FILE *fp)
{
	FB              *fb;
	Controls    *ctrl;
	Dike        *dike;
	PetscInt   nD, numDike;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ctrl = &jr->ctrl;
	if (!ctrl->actDike) PetscFunctionReturn(0);   // only execute this function if dikes are active

	ierr = FBLoad(&fb, PETSC_TRUE);
	ierr = TSSolCreate(ts, fb); 				CHKERRQ(ierr);

	numDike    = dbdike->numDike; // number of dikes

// create dike database
	ierr = DBDikeCreate(dbdike, dbm, fb, jr, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = FBDestroy(&fb); CHKERRQ(ierr);

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
