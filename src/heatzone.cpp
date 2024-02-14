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

		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	PetscCall(FBFreeBlocks(fb));

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode DBReadHeatZone(DBPropHeatZone *dbheatzone, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
	// read heatzone parameter from file
	BCCtx *bc;
	Scaling *scal;
	HeatZone *heatzone;
	PetscInt i, ID;
	PetscScalar  invt, Box[6];
	char Parameter[_str_len_];
	
	PetscFunctionBeginUser;

	// access context
	scal = dbm->scal;
	bc = jr->bc;

	// HeatZone ID
	PetscCall(getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbheatzone->numHeatZone - 1));
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
	invt = 1/jr->scal->time;

	// read and store heatzone  parameters
	PetscCall(getStringParam(fb, _REQUIRED_, "HeatFunction", Parameter, "q_hotspot"));
	PetscCall(getScalarParam(fb, _REQUIRED_, "HZ_Bounds", heatzone->bounds, 6, scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "AsthenoTemp", &heatzone->asthenoTemp, 1, 1));
	PetscCall(getScalarParam(fb, _REQUIRED_, "rho", &heatzone->rho, 1, scal->density));
	PetscCall(getScalarParam(fb, _REQUIRED_, "Cp", &heatzone->Cp, 1, scal->cpecific_heat));

	if (!strcmp(Parameter, "q_hotspot"))
	{
		heatzone->HeatFunction = 0;
		PetscCall(getScalarParam(fb, _REQUIRED_, "HeatRate", &heatzone->heatRate, 1, scal->strain_rate));

	}
	else if (!strcmp(Parameter, "q_ridge"))
	{
		heatzone->HeatFunction = 1;
		PetscCall(getScalarParam(fb, _OPTIONAL_, "SpreadingRate", &heatzone->spreadingRate, 1, scal->velocity)); // R

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

	// scaling
	heatzone->asthenoTemp = (heatzone->asthenoTemp + scal->Tshift) / scal->temperature;

	// print HeatZone info to terminal
	if (PrintOutput)
	{
		// get heatzone bounds
		for (i=0; i<6; i++){ Box[i] = heatzone->bounds[i]*scal->length; }
		
		// print block
		if (heatzone->HeatFunction == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   HeatZone [%lld]: hotspot heating\n", (LLD)(heatzone->ID));
			PetscPrintf(PETSC_COMM_WORLD, "     Bounds     : [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0],Box[1],Box[2],Box[3],Box[4],Box[5], scal->lbl_length);
			PetscPrintf(PETSC_COMM_WORLD, "     Parameters : AsthenoTemp = %1.0f %s, HeatRate = %g %s\n", heatzone->asthenoTemp * scal->temperature - scal->Tshift, scal->lbl_temperature, heatzone->heatRate, scal->lbl_strain_rate);
			PetscPrintf(PETSC_COMM_WORLD, "                  rho = %1.0f %s, Cp = %g %s\n", heatzone->rho * scal->density, scal->lbl_density, heatzone->Cp * scal->cpecific_heat, scal->lbl_cpecific_heat);
		}
		else if (heatzone->HeatFunction == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   HeatZone [%lld]: ridge heating\n", (LLD)(heatzone->ID));
			PetscPrintf(PETSC_COMM_WORLD, "     Bounds     : [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0],Box[1],Box[2],Box[3],Box[4],Box[5], scal->lbl_length);
			PetscPrintf(PETSC_COMM_WORLD, "     Parameters : AsthenoTemp = %1.0f %s, SpreadingRate = %1.1f %s\n", heatzone->asthenoTemp * scal->temperature - scal->Tshift, scal->lbl_temperature, heatzone->spreadingRate * scal->velocity, scal->lbl_velocity);
			PetscPrintf(PETSC_COMM_WORLD, "                  rho = %1.0f %s, Cp = %g %s\n", heatzone->rho * scal->density, scal->lbl_density, heatzone->Cp * scal->cpecific_heat, scal->lbl_cpecific_heat);
		}
	}

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------

PetscErrorCode GetHeatZoneSource(JacRes *jr,
								 Material_t *phases,
								 PetscScalar &Tc,
								 PetscScalar *phRat,
								 PetscScalar &rho_A,
								 PetscScalar &y_c,
								 PetscScalar &x_c,
								 PetscScalar &z_c,
								 PetscInt J)

{
	HeatZone *heatzone;
	PetscInt nHZ, numHeatZone, AirPhase;
	PetscScalar asthenoTemp, heatRate, spreadingRate;
	PetscScalar rho, Cp, hzRat, st_dev, F_x, delta_x;
	PetscScalar hz_left, hz_right, hz_width, hz_x_cent;
	PetscScalar hz_front, hz_back, hz_bottom, hz_top;
//	PetscScalar hz_length, hz_y_cent; // *3D

//	Dike *dike;
//	Ph_trans_t *CurrPhTr;
//	PetscInt nPtr, numPhtr, nD, numDike, nPhase, numPhases;
//	Material_t *mat;
//	PetscInt nHZ, nPtr, nPhase, nD;
//	PetscScalar v_spread, left, right, front, back, x_distance;

	PetscFunctionBeginUser;

//	numPhtr = jr->dbm->numPhtr;				   // number of phase transitions
//	numDike = jr->dbdike->numDike;			   // number of dikes
//	numPhases = jr->dbm->numPhases;

	numHeatZone = jr->dbheatzone->numHeatZone; // number of heatzones

	for (nHZ = 0; nHZ < numHeatZone; nHZ++) // loop through all heatzone blocks
	{
		// access the necessary parameters of the heatzone block
		heatzone = jr->dbheatzone->matHeatZone + nHZ; // getting hz params

		// material parameters
		asthenoTemp = heatzone->asthenoTemp;
		heatRate = heatzone->heatRate;
		spreadingRate = heatzone->spreadingRate;
		rho = heatzone->rho;
		Cp = heatzone->Cp;

		// geometric parameters
		hz_left = heatzone->bounds[0];
		hz_right = heatzone->bounds[1];
		hz_front = heatzone->bounds[2];
		hz_back = heatzone->bounds[3];
		hz_bottom = heatzone->bounds[4];
		hz_top = heatzone->bounds[5];
		hz_width = PetscSqrtScalar(pow(hz_right - hz_left, 2));
		st_dev = hz_width / (2 * PetscSqrtScalar(2 * log(2)));
		hz_x_cent = (hz_right + hz_left) / 2;

		// if we are close to the heatzone
		if (x_c > (hz_x_cent - hz_width) && x_c < (hz_x_cent + hz_width) && y_c > hz_front && y_c < hz_back && z_c > hz_bottom && z_c < hz_top)
		{
			// compute environmental parameters
			delta_x = PetscSqrtScalar(pow(hz_x_cent - x_c, 2)); // *revisit
			F_x = (hz_width / (st_dev * PetscSqrtScalar(2 * PETSC_PI))) * exp(-pow(delta_x, 2) / (2 * pow(st_dev, 2)));
			/* 		hz_width = sqrt(pow(hz_right - hz_left, 2) + pow(hz_back - hz_front, 2)); // *3D
					hz_length = sqrt(pow(hz_back - hz_front, 2)); // *3D
					hz_y_cent = hz_back - hz_front/2; // *3D */

			// calculate not in air phase ratio
			hzRat = 1;
			AirPhase = jr->surf->AirPhase;
			if (AirPhase != -1)
			{
				hzRat -= phRat[AirPhase];
			}

			// calculate heating contribution
			if (heatzone->HeatFunction == 0) // q_hotspot
			{
				rho_A += hzRat * rho * Cp * heatRate * F_x * (asthenoTemp - Tc) * jr->scal->dissipation_rate; // * (jr->scal->stress_si / jr->scal->time) -> at 1 yr  // jr->ts->dt
			}
			else if (heatzone->HeatFunction == 1) // q_ridge
			{
				rho_A += hzRat * rho * Cp * F_x * (asthenoTemp - Tc) * spreadingRate / hz_width;
			}
			else // should not be possible
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "What did you do?!");
			}

			/* 		// is dike heating used in the simulation?
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
								PetscPrintf(PETSC_COMM_WORLD, "###### Run SubtractDikeHeatSource ######\n");
							}
						}
					}

					for (nPtr = 0; nPtr < numPhtr; nPtr++) // loop over all phase transitions blocks
					{
						// access the parameters of the phasetranstion block
						CurrPhTr = jr->dbm->matPhtr + nPtr; // important params: ID, type (NotInAirBox), and bounds (left = bounds[0], right = bounds [1])

						// Get HeatZone[nHZ] parameters
						if (CurrPhTr->ID == heatzone->PhaseTransID)
						{
							{
								for (nPtr = 0; nPtr < numPhtr; nPtr++) // loop over all phase transitions blocks again
								{
									// access the parameters of the phasetranstion block
									CurrPhTr = jr->dbm->matPhtr + nPtr; // now we're looking for phase trasition of current cell

									// loop through phases
									for (nPhase = 0; nPhase < numPhases; nPhase++)
									{
										// if the current phase exists within HeatZone
										if (phRat[nPhase] > 0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
										{
										}
									}
								}
							}
						}
					} */
		}
	}
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