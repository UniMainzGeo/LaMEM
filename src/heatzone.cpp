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

	This file defines properties for a heatzone intended to be used for heating
	within a zone without relying on phase ratio (i.e. set PhaseInside = -1
	in phase transition law). However, this does still rely on density and
	specific heat values from the associated phase material properties.

*/
//---------------------------------------------------------------------------
//.................. HEATZONE PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "phase.h"
#include "parsing.h"
#include "JacRes.h"
#include "heatzone.h"
#include "dike.h"
#include "constEq.h"
#include "bc.h"
#include "tssolve.h"
#include "scaling.h"
#include "fdstag.h"
#include "tools.h"
#include "surf.h"
#include "advect.h"

/* // debugging tools *djking
#include <iostream>
		// Pause for debugging (doesn't keep variable values)
		std::cout << "Press Enter to continue..."; std::cin.get();  return 0; */

//---------------------------------------------------------------------------
PetscErrorCode DBHeatZoneCreate(DBPropHeatZone *dbheatzone, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
	// read all heatzone parameter blocks from file

	PetscInt jj; // numHeatZone;

	PetscFunctionBeginUser;

	if (!jr->ctrl.actHeatZone)
		PetscFunctionReturn(0); // only execute this function if heatzones are active

	//====================
	// HeatZone Parameters
	//====================

	// setup block access mode
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<HeatZoneStart>", "<HeatZoneEnd>"));

	if (fb->nblocks)
	{
		// print overview of heatzone blocks from file
		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "HeatZone blocks : \n");

		// initialize ID for consistency checks
		for (jj = 0; jj < _max_num_heatzone_; jj++)
			dbheatzone->matHeatZone[jj].ID = -1;

		// error checking
		if (fb->nblocks > _max_num_heatzone_)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many heatzones specified! Max allowed: %lld", (LLD)_max_num_heatzone_);

		// store actual number of heatzone blocks
		dbheatzone->numHeatZone = fb->nblocks;

		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		// read each individual heatzone block
		for (jj = 0; jj < fb->nblocks; jj++)
		{
			PetscCall(DBReadHeatZone(dbheatzone, dbm, fb, jr, PrintOutput));
			fb->blockID++;
		}
	}

	PetscCall(FBFreeBlocks(fb));

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode DBReadHeatZone(DBPropHeatZone *dbheatzone, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
	// read heatzone parameter from file
	BCCtx *bc;
	HeatZone *heatzone;
	PetscInt ID;
	Scaling *scal;
	char Parameter[_str_len_];
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	scal = dbm->scal;
	bc = jr->bc;

	// HeatZone ID
	ierr = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbheatzone->numHeatZone - 1);
	CHKERRQ(ierr);
	fb->ID = ID;

	// get pointer to specified heatzone parameters
	heatzone = dbheatzone->matHeatZone + ID;

	// check ID
	if (heatzone->ID != -1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate of HeatZone option!");
	}

	// set ID
	heatzone->ID = ID;

	// read and store heatzone  parameters
	ierr = getIntParam(fb, _REQUIRED_, "PhaseID", &heatzone->PhaseID, 1, dbm->numPhases - 1);
	CHKERRQ(ierr);
	ierr = getIntParam(fb, _REQUIRED_, "PhaseTransID", &heatzone->PhaseTransID, 1, dbm->numPhtr - 1);
	CHKERRQ(ierr);
	ierr = getStringParam(fb, _REQUIRED_, "HeatFunction", Parameter, "q_hotspot");

	if (!strcmp(Parameter, "q_hotspot"))
	{
		heatzone->HeatFunction = 0;
		heatzone->heatRate = 1e-10;
		ierr = getScalarParam(fb, _REQUIRED_, "AsthenoTemp", &heatzone->asthenoTemp, 1, 1);
		CHKERRQ(ierr);
		ierr = getScalarParam(fb, _REQUIRED_, "HeatRate", &heatzone->heatRate, 1, 1);
		CHKERRQ(ierr);

		heatzone->asthenoTemp = (heatzone->asthenoTemp + scal->Tshift) / scal->temperature; // scaling temperature input
	}
	else if (!strcmp(Parameter, "q_ridge"))
	{
		heatzone->HeatFunction = 1;
		ierr = getScalarParam(fb, _REQUIRED_, "AsthenoTemp", &heatzone->asthenoTemp, 1, 1);
		CHKERRQ(ierr);
		ierr = getScalarParam(fb, _OPTIONAL_, "SpreadingRate", &heatzone->spreadingRate, 1, scal->velocity);
		CHKERRQ(ierr); // R

		heatzone->asthenoTemp = (heatzone->asthenoTemp + scal->Tshift) / scal->temperature; // scaling temperature input

		if (heatzone->spreadingRate) // needs dependence on local spreading rate to be useful
		{
			heatzone->spreadingRate = 2 * heatzone->spreadingRate; // convert to full spreading rate
		}
		else
		{
			heatzone->spreadingRate = 2 * PetscAbs(bc->velin); // get bc input and convert to full spreading rate
		}
	}
	else
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown parameter for HeatFunction %s [q_hotspot; q_ridge]", Parameter);
	}

	// print HeatZone info to terminal
	if (PrintOutput)
	{
		if (heatzone->HeatFunction == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "  HeatZone parameters ID[%lld]: PhaseTransID=%lld PhaseID=%lld AsthenoTemp=%g, HeatRate=%g \n",
						(LLD)(heatzone->ID), (LLD)(heatzone->PhaseTransID), (LLD)(heatzone->PhaseID), heatzone->asthenoTemp * scal->temperature - scal->Tshift, heatzone->heatRate);
		}
		else if (heatzone->HeatFunction == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "  HeatZone parameters ID[%lld]: PhaseTransID=%lld PhaseID=%lld AsthenoTemp=%g, SpreadingRate=%g \n",
						(LLD)(heatzone->ID), (LLD)(heatzone->PhaseTransID), (LLD)(heatzone->PhaseID), heatzone->asthenoTemp * scal->temperature - scal->Tshift, heatzone->spreadingRate * scal->velocity);
		}
		PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------

PetscErrorCode GetHeatZoneSource(JacRes *jr,
								 Material_t *phases,
								 PetscScalar &Tc,
								 PetscScalar *phRat, // phase ratios in the control volume
								 PetscScalar &k,
								 PetscScalar &rho_A,
								 PetscScalar &y_c,
								 PetscInt J)

{
	HeatZone *heatzone;
	Dike *dike;
	Ph_trans_t *CurrPhTr;
	Material_t *mat;
	PetscInt numHeatZone, numPhtr, numPhases, numDike;
	PetscInt nHZ, nPtr, nPhase, nD;
	PetscScalar heatFunct, asthenoTemp, heatRate, spreadingRate;
	PetscScalar rho, Cp;
	PetscScalar hz_left, hz_right, hz_width, hz_x_cent;
//	PetscScalar hz_front, hz_back, hz_length, hz_y_cent; // *3D
	PetscScalar pi, st_dev, F_x, delta_x; // *revisit
	PetscScalar x_cellc, y_cellc, z_cellc; // *revisit

//	PetscScalar v_spread, left, right, front, back, x_distance; // SubtractDikeHeatSource *revisit

PetscInt debug;

	PetscFunctionBeginUser;

debug = 0;

	numPhtr = jr->dbm->numPhtr;				   // number of phase transitions
	numHeatZone = jr->dbheatzone->numHeatZone; // number of heatzones
	numDike = jr->dbdike->numDike;			   // number of dikes
	numPhases = jr->dbm->numPhases;

	x_cellc = *(jr->fs->dsx.ccoor);
	y_cellc = *(jr->fs->dsy.ccoor);
	z_cellc = *(jr->fs->dsz.ccoor);

if (debug == 1) 
{
	PetscPrintf(PETSC_COMM_WORLD, "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	PetscPrintf(PETSC_COMM_WORLD, "------------------ Compute Heating ----------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n \n");
	PetscPrintf(PETSC_COMM_WORLD, "... checking for phase transition laws of heating zones\n");
}
	for (nPtr = 0; nPtr < numPhtr; nPtr++) // loop over all phase transitions blocks
	{
		// access the parameters of the phasetranstion block
		CurrPhTr = jr->dbm->matPhtr + nPtr; // important params: ID, type (NotInAirBox), and bounds (left = bounds[0], right = bounds [1])
		hz_left = CurrPhTr->bounds[0];
		hz_right = CurrPhTr->bounds[1];
//		hz_front = CurrPhTr->bounds[2]; *revisit
//		hz_back = CurrPhTr->bounds[3]; *revisit

if (debug == 1) {PetscPrintf(PETSC_COMM_WORLD, "... checking phase transition %lld\n", (LLD)(nPtr));}

		for (nHZ = 0; nHZ < numHeatZone; nHZ++) // loop through all heatzone blocks
		{
			// access the necessary parameters of the heatzone block
			heatzone = jr->dbheatzone->matHeatZone + nHZ; // getting hz params
			heatFunct = heatzone->HeatFunction;
			asthenoTemp = heatzone->asthenoTemp;
			heatRate = heatzone->heatRate;
			spreadingRate = heatzone->spreadingRate;


			// if in the HeatZone 
			if (CurrPhTr->ID == heatzone->PhaseTransID) // its a HeatZone
			{

if (debug == 1) PetscPrintf(PETSC_COMM_WORLD, "Phase Transition %lld is HeatZone %lld!!!!!\n", (LLD)(nPtr), (LLD)(nHZ));
				
/* 				// if this cell is within HeatZone
				if (x_cellc >)
				{
					// code
				}
				 */
       
				// loop through phases
				for (nPhase = 0; nPhase < numPhases; nPhase++)
				{
if (debug == 1) PetscPrintf(PETSC_COMM_WORLD, "... checking phase ID %lld for heating contributions\n", (LLD)(nPhase));
					
					// if the current phase exists within HeatZone
					if(phRat[nPhase]>0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J]) // what is J? *revisit
					{
						// compute heating
						mat = &phases[nPhase]; // getting phase params
						rho = mat->rho;
						Cp = mat->Cp;
						
						// compute environmental params
						pi = 3.14159265358979323846;
						hz_width = sqrt(pow(hz_right - hz_left, 2));
						hz_x_cent = hz_right - hz_left/2;
						delta_x = sqrt(pow(hz_x_cent - x_cellc, 2)); // *revisit
						st_dev = hz_width/(2*sqrt(2*log(2)));
						F_x = (hz_width / (st_dev * sqrt(2 * pi))) * exp(-pow(delta_x, 2) / (2 * pow(st_dev, 2)));
/*						hz_width = sqrt(pow(hz_right - hz_left, 2) + pow(hz_back - hz_front, 2)); // *3D
						hz_length = sqrt(pow(hz_back - hz_front, 2)); // *3D
						hz_y_cent = hz_back - hz_front/2; // *3D */


if (debug == 1) PetscPrintf(PETSC_COMM_WORLD, "!!!!!!!!!!!!!!!! COMPUTE !!!!!!!!!!!!!!!\n");
						
						// calculate phaseid heating contribution from phRat
						if (heatFunct == 0) // q_hotspot = rho_A;
						{						
							rho_A += phRat[nPhase]*rho*Cp*heatRate*F_x*(asthenoTemp-Tc);
						}
						else if (heatFunct == 1) // q_ridge *revisit
						{
							rho_A += phRat[nPhase]*rho*Cp*F_x*(asthenoTemp-Tc)*spreadingRate/hz_width; // *revisit
						}
						else // should not be possible
						{
							SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "What did you do?!");
						}
						

						// is dike heating used in the simulation?
						if (jr->ctrl.actDike && jr->ctrl.dikeHeat)
						{
							// loop through dike blocks to find PhaseID
							for (nD = 0; nD < numDike; nD++)
							{
								dike = jr->dbdike->matDike + nD;

								// if phase is dike zone phase
								if (dike->PhaseID == nPhase)
								{
									// Subtract heat added via diking

if (debug == 1) PetscPrintf(PETSC_COMM_WORLD, "###### Run SubtractDikeHeatSource ######\n");
								}
							}
						}
					}
				}
			}
		}
	}

if (debug == 1) 
{
	PetscPrintf(PETSC_COMM_WORLD, "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	PetscPrintf(PETSC_COMM_WORLD, "------------ Heating Calculation Complete ---------\n");
	PetscPrintf(PETSC_COMM_WORLD, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n \n");
}
	/* 	for (nPtr = 0; nPtr < numPhtr; nPtr++) // loop over all phase transitions blocks
		{

			// access the parameters of the phasetranstion block
			CurrPhTr = jr->dbm->matPhtr + nPtr;

			for (nHZ = 0; nHZ < numHeatZone; nHZ++) // loop through all heatzone blocks
			{
				// access the parameters of the heatzone depending on the heatzone block
				heatzone = jr->dbheatzone->matHeatZone + nHZ; */

	/* 			// access the phase ID of the heatzone parameters of each heatzone
				i = heatzone->PhaseID; */
	// not needed if using phase transitions in the cell... *revisit

	/* 			if (CurrPhTr->ID == heatzone->PhaseTransID) // compare the phaseTransID associated with the heatzone with the actual ID of the phase transition in this cell
	//			if (CurrPhTr->ID is within Bounds of HeatZone) // compare the phaseTransID associated with the heatzone with the actual ID of the phase transition in this cell
				{

					// if in the heatzone
					if (CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
					{
						nsegs = CurrPhTr->nsegs;
						if (heatzone->Mb == heatzone->Mf && heatzone->Mc < 0.0) // constant M
						{
							M = heatzone->Mf;
							v_spread = PetscAbs(bc->velin);
							left = CurrPhTr->celly_xboundL[J];
							right = CurrPhTr->celly_xboundR[J];
							tempheatzoneRHS = M * 2 * v_spread / PetscAbs(left - right);
						}
						else if (heatzone->Mc >= 0.0) // Mf, Mc and Mb
						{
							left = CurrPhTr->celly_xboundL[J];
							right = CurrPhTr->celly_xboundR[J];
							front = CurrPhTr->ybounds[0];
							back = CurrPhTr->ybounds[2 * nsegs - 1];
							v_spread = PetscAbs(bc->velin);

							if (y_c >= heatzone->y_Mc)
							{
								// linear interpolation between different M values, Mc is M in the middle, acts as M in front, Mb is M in back
								y_distance = y_c - heatzone->y_Mc;
								M = heatzone->Mc + (heatzone->Mb - heatzone->Mc) * (y_distance / (back - heatzone->y_Mc));
								tempheatzoneRHS = M * 2 * v_spread / PetscAbs(left - right);
							}
							else
							{
								// linear interpolation between different M values, Mf is M in front, Mc acts as M in back
								y_distance = y_c - front;
								M = heatzone->Mf + (heatzone->Mc - heatzone->Mf) * (y_distance / (heatzone->y_Mc - front));
								tempheatzoneRHS = M * 2 * v_spread / PetscAbs(left - right);
							}
						}
						else if (heatzone->Mb != heatzone->Mf && heatzone->Mc < 0.0) // only Mf and Mb, they are different
						{
							left = CurrPhTr->celly_xboundL[J];
							right = CurrPhTr->celly_xboundR[J];
							front = CurrPhTr->ybounds[0];
							back = CurrPhTr->ybounds[2 * nsegs - 1];
							v_spread = PetscAbs(bc->velin);

							// linear interpolation between different M values, Mf is M in front, Mb is M in back
							y_distance = y_c - front;
							M = heatzone->Mf + (heatzone->Mb - heatzone->Mf) * (y_distance / (back - front));
							tempheatzoneRHS = M * 2 * v_spread / PetscAbs(left - right);
						}
						else
						{
							tempheatzoneRHS = 0.0;
						}

						mat = &phases[i];

						// adjust k and heat source according to Behn & Ito [2008]
						if (Tc < mat->T_liq && Tc > mat->T_sol)
						{
							kfac += phRat[i] / (1 + (mat->Latent_hx / (mat->Cp * (mat->T_liq - mat->T_sol))));
							rho_A += phRat[i] * (mat->rho * mat->Cp) * (mat->T_liq - Tc) * tempheatzoneRHS; // Cp*rho not used in the paper, added to conserve units of rho_A
						}
						else if (Tc <= mat->T_sol)
						{
							kfac += phRat[i];
							rho_A += phRat[i] * (mat->rho * mat->Cp) * ((mat->T_liq - Tc) + mat->Latent_hx / mat->Cp) * tempheatzoneRHS;
						}
						else if (Tc >= mat->T_liq)
						{
							kfac += phRat[i];
						}
						// end adjust k and heat source according to Behn & Ito [2008]

						k = kfac * k;

					} // end if phRat and xboundR>xboundL
				}	  // close phase transition and phase ID comparison
			}		  // end heatzone block loop
		}			  // close phase transition block loop */

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------

/* PetscErrorCode SubtractDikeHeatSource(JacRes *jr,
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
  PetscInt     i, numDike, nHZ, nPtr, numPhtr, nsegs;
  PetscScalar  v_spread, left, right, front, back, M, kfac, tempdikeRHS;
  PetscScalar  y_distance;

  PetscFunctionBeginUser;

  numDike    = jr->dbdike->numDike; // number of dikes
  numPhtr    = jr->dbm->numPhtr;
  bc         = jr->bc;

  nPtr = 0;
  nHZ   = 0;
  kfac = 0;

  for(nPtr=0; nPtr<numPhtr; nPtr++)   // loop over all phase transitions blocks
	{

	  // access the parameters of the phasetranstion block
	  CurrPhTr = jr->dbm->matPhtr+nPtr;

	  for(nHZ = 0; nHZ < numDike; nHZ++) // loop through all dike blocks
		{
		  // access the parameters of the dike depending on the dike block
		  dike = jr->dbdike->matDike+nHZ;

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
} */