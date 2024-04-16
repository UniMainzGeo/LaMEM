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
	PetscScalar Box[6];
	char Funct[_str_len_], Dim[_str_len_];

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

	// read and store heatzone  parameters
	PetscCall(getStringParam(fb, _REQUIRED_, "HeatFunction", Funct, NULL));
	PetscCall(getStringParam(fb, _REQUIRED_, "FunctType", Dim, NULL));
	PetscCall(getScalarParam(fb, _REQUIRED_, "HZ_Bounds", heatzone->bounds, 6, scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "AsthenoTemp", &heatzone->asthenoTemp, 1, 1));
	PetscCall(getScalarParam(fb, _REQUIRED_, "rho", &heatzone->rho, 1, scal->density));
	PetscCall(getScalarParam(fb, _REQUIRED_, "Cp", &heatzone->Cp, 1, scal->cpecific_heat));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "TimeStart", &heatzone->timeStart, 1, scal->time));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "TempStart", &heatzone->tempStart, 1, 1));

	// error checking bounds
	if ((heatzone->bounds[1] < heatzone->bounds[0]) | (heatzone->bounds[3] < heatzone->bounds[2]) | (heatzone->bounds[5] < heatzone->bounds[4]))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Negative width detected in x, y, or z-direction: Check HZ_Bounds");
	}

	// check heat function inputs
	if (!strcmp(Funct, "q_hotspot"))
	{
		heatzone->HeatFunction = 0;
		PetscCall(getScalarParam(fb, _REQUIRED_, "HeatRate", &heatzone->heatRate, 1, scal->strain_rate));
	}
	else if (!strcmp(Funct, "q_ridge"))
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
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown parameter for HeatFunction %s [q_hotspot; q_ridge]", Funct);
	}

	// check heat zone dimensionality
	if (!strcmp(Dim, "1d_x-gauss"))
	{
		heatzone->FunctType = 0; // default
	}
	else if (!strcmp(Dim, "2d_xy-gauss"))
	{
		heatzone->FunctType = 1;
	}
	else if (!strcmp(Dim, "3d_xyz-gauss"))
	{
		heatzone->FunctType = 2;
	}
	else
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown parameter for HZ_Type %s [1d_gauss; 2d_gauss]", Dim); // 1d_ydir, 2d_elliptical
	}

	// temperature scaling
	heatzone->asthenoTemp = (heatzone->asthenoTemp + scal->Tshift) / scal->temperature;
	heatzone->tempStart = (heatzone->tempStart + scal->Tshift) / scal->temperature;

	// print HeatZone info to terminal
	if (PrintOutput)
	{
		// get heatzone bounds
		for (i = 0; i < 6; i++)
		{
			Box[i] = heatzone->bounds[i] * scal->length;
		}

		// print block
		if (heatzone->HeatFunction == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   HeatZone [%lld] : hotspot heating, %s\n", (LLD)(heatzone->ID), Dim);
			PetscPrintf(PETSC_COMM_WORLD, "     Bounds     : [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0], Box[1], Box[2], Box[3], Box[4], Box[5], scal->lbl_length);
			PetscPrintf(PETSC_COMM_WORLD, "     Parameters : AsthenoTemp = %1.0f %s, HeatRate = %g %s\n", heatzone->asthenoTemp * scal->temperature - scal->Tshift, scal->lbl_temperature, heatzone->heatRate, scal->lbl_strain_rate);
			PetscPrintf(PETSC_COMM_WORLD, "                  rho = %1.0f %s, Cp = %g %s\n", heatzone->rho * scal->density, scal->lbl_density, heatzone->Cp * scal->cpecific_heat, scal->lbl_cpecific_heat);
			PetscPrintf(PETSC_COMM_WORLD, "                  startTime = %g %s, startTemp = %1.0f %s\n", heatzone->timeStart * scal->time, scal->lbl_time, heatzone->tempStart * scal->temperature - scal->Tshift, scal->lbl_temperature);
		}
		else if (heatzone->HeatFunction == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   HeatZone [%lld] : ridge heating, %s\n", (LLD)(heatzone->ID), Dim);
			PetscPrintf(PETSC_COMM_WORLD, "     Bounds     : [%1.1f; %1.1f; %1.1f; %1.1f; %1.1f; %1.1f] %s \n", Box[0], Box[1], Box[2], Box[3], Box[4], Box[5], scal->lbl_length);
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
								 PetscInt J,
								 PetscScalar sxx_eff_ave_cell)

{
	HeatZone *heatzone;
	PetscInt nHZ, numHeatZone, AirPhase;
	PetscScalar rho, Cp, asthenoTemp, heatRate, spreadingRate;
	PetscScalar hzRat, st_dev, F_x, delta_hz_cent, hz_ind;
	PetscScalar hz_left, hz_right, hz_width, hz_x_cent;
	PetscScalar hz_front, hz_back, hz_y_cent;
	PetscScalar hz_bottom, hz_top, hz_z_cent;
	PetscScalar hz_contr, timeRat;

	PetscFunctionBeginUser;

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

		// heatzone geometry
		hz_left = heatzone->bounds[0];	 // left
		hz_right = heatzone->bounds[1];	 // right
		hz_front = heatzone->bounds[2];	 // front
		hz_back = heatzone->bounds[3];	 // back
		hz_bottom = heatzone->bounds[4]; // top
		hz_top = heatzone->bounds[5];	 // bottom

		hz_width = hz_right - hz_left; // all gaussian dependent on x-dir width!
		st_dev = hz_width / (2 * PetscSqrtScalar(2 * log(2)));
		hz_x_cent = (hz_right + hz_left) / 2;
		hz_y_cent = (hz_back + hz_front) / 2;
		hz_z_cent = (hz_top + hz_bottom) / 2;

		// is timestep past TimeStart?
		timeRat = 1; // determines the amount of heating applied over timestep; default of 1 implies heating over entire timestep
		if (heatzone->timeStart > 0)
		{
			if (jr->ts->time >= heatzone->timeStart && (jr->ts->time - heatzone->timeStart) < jr->ts->dt && (jr->ts->time - heatzone->timeStart) != 0)
			{
				timeRat = (jr->ts->time - heatzone->timeStart) / jr->ts->dt;
			}
			else if (jr->ts->time < heatzone->timeStart)
			{
				timeRat = 0;
			}
		}
		
		// check if within heatzone and if cell temperature is above TempStart but less than AsthenoTemp
		hz_ind = 0;					  // heatzone boolean
		if (heatzone->FunctType == 0) // 1d_x-gauss
		{
			if (x_c > (hz_x_cent - hz_width) && x_c < (hz_x_cent + hz_width) && y_c > hz_front && y_c < hz_back && z_c > hz_bottom && z_c < hz_top && Tc >= heatzone->tempStart && Tc <= heatzone->asthenoTemp)
			{
				hz_ind = 1;
				delta_hz_cent = abs(hz_x_cent - x_c); // distance from the center of hz
			}
		}
		else if (heatzone->FunctType == 1) // 2d_xy-gauss
		{
			if (x_c > (hz_x_cent - hz_width) && x_c < (hz_x_cent + hz_width) && y_c > (hz_y_cent - hz_width) && y_c < (hz_y_cent + hz_width) && z_c > hz_bottom && z_c < hz_top && Tc >= heatzone->tempStart && Tc <= heatzone->asthenoTemp)
			{
				hz_ind = 1;
				delta_hz_cent = PetscSqrtScalar(pow(hz_x_cent - x_c, 2) + pow(hz_y_cent - y_c, 2)); // distance from the center of hz
			}
		}
		else if (heatzone->FunctType == 2) // 3d_xyz-gauss "little ball of heat"
		{
			if (x_c > (hz_x_cent - hz_width) && x_c < (hz_x_cent + hz_width) && y_c > (hz_y_cent - hz_width) && y_c < (hz_y_cent + hz_width) && z_c > (hz_z_cent - hz_width) && z_c < (hz_z_cent + hz_width) && Tc >= heatzone->tempStart && Tc <= heatzone->asthenoTemp)
			{
				hz_ind = 1;
				delta_hz_cent = PetscSqrtScalar(pow(hz_x_cent - x_c, 2) + pow(hz_y_cent - y_c, 2) + pow(hz_z_cent - z_c, 2)); // distance from the center of hz
			}
		}

		// if we are close to the heatzone bounds
		if (hz_ind == 1)
		{
			// compute environmental parameters
			F_x = (hz_width / (st_dev * PetscSqrtScalar(2 * PETSC_PI))) * exp(-pow(delta_hz_cent, 2) / (2 * pow(st_dev, 2)));

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
				hz_contr = timeRat * hzRat * rho * Cp * heatRate * F_x * (asthenoTemp - Tc); // * invt; // jr->scal->dissipation_rate; // * (jr->scal->stress_si / jr->scal->time) -> at 1 yr  // jr->ts->dt
			}
			else if (heatzone->HeatFunction == 1) // q_ridge
			{
				hz_contr = timeRat * hzRat * rho * Cp * F_x * (asthenoTemp - Tc) * spreadingRate / hz_width;
			}
			else // should not be possible
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "What did you do?!");
			}
			// is dike heating used in the simulation?
			if (jr->ctrl.actDike && jr->ctrl.dikeHeat)
			{
				// Subtract heat added via diking
				PetscCall(SubtractDikeHeatSource(jr, phases, Tc, phRat, hz_contr, y_c, J, sxx_eff_ave_cell));
			}

			rho_A += hz_contr; // add heating to energy equation as source term
		}
	}
	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------

PetscErrorCode SubtractDikeHeatSource(JacRes *jr,
									  Material_t *phases,
									  PetscScalar &Tc,
									  PetscScalar *phRat,
									  PetscScalar &hz_contr,
									  PetscScalar &y_c,
									  PetscInt J,
									  PetscScalar sxx_eff_ave_cell)

{
	// parameters to determine dilation term
	BCCtx *bc;
	Dike *dike;
	Ph_trans_t *CurrPhTr;
	PetscInt i, nD, nPtr, numDike, numPhtr, nsegs;
	PetscScalar v_spread, M, left, right, front, back;
	PetscScalar y_distance, tempdikeRHS;
	PetscScalar P_comp, div_max, M_rat, zeta;
	PetscScalar dike_contr;

	// heating parameters
	Material_t *mat;

	PetscFunctionBeginUser;

	numDike = jr->dbdike->numDike; // number of dikes
	numPhtr = jr->dbm->numPhtr;
	bc = jr->bc;

	nPtr = 0;
	nD = 0;

	for (nPtr = 0; nPtr < numPhtr; nPtr++) // loop over all phase transitions blocks
	{

		// access the parameters of the phasetranstion block
		CurrPhTr = jr->dbm->matPhtr + nPtr;

		for (nD = 0; nD < numDike; nD++) // loop through all dike blocks
		{
			// access the parameters of the dike depending on the dike block
			dike = jr->dbdike->matDike + nD;

			// access the phase ID of the dike parameters of each dike
			i = dike->PhaseID;

			if (CurrPhTr->ID == dike->PhaseTransID) // compare the phaseTransID associated with the dike with the actual ID of the phase transition in this cell
			{
				// if in the dike zone
				if (phRat[i] > 0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
				{
					nsegs = CurrPhTr->nsegs;

					if (dike->Mb == dike->Mf && dike->Mc < 0.0) // spatially constant M
					{
						M = dike->Mf;
						v_spread = PetscAbs(bc->velin);
						left = CurrPhTr->celly_xboundL[J];
						right = CurrPhTr->celly_xboundR[J];

						if (jr->ctrl.var_M)
						{
							P_comp = sxx_eff_ave_cell - dike->Ts;
							M_rat = M; // M ratio *revisit
							div_max = M_rat * 2 * (v_spread / (right - left));

							if (P_comp > 0) // diking occurs
							{
								if (dike->Mf == 0)
								{
									zeta = dike->A * (dike->zeta_0 / P_comp);
								}
								else
								{
									zeta = dike->A * (dike->zeta_0 / P_comp) + P_comp / div_max;
								}

								tempdikeRHS = P_comp / zeta;
							}
							else // diking DOES NOT occur
							{
								tempdikeRHS = 0.0;
							}
						}
						else // not using var_M
						{
							tempdikeRHS = M * 2 * v_spread / (right - left);
						}
					}
					else if (dike->Mc >= 0.0) // Mf, Mc and Mb are all user defined
					{
						if (jr->ctrl.var_M) // check varaible M option isn't used
						{
							SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid option: var_M option requires uniform M");
						}

						left = CurrPhTr->celly_xboundL[J];
						right = CurrPhTr->celly_xboundR[J];
						front = CurrPhTr->ybounds[0];
						back = CurrPhTr->ybounds[2 * nsegs - 1];
						v_spread = PetscAbs(bc->velin);

						if (y_c >= dike->y_Mc)
						{
							// linear interpolation between different M values, Mc is M in the middle, acts as M in front, Mb is M in back
							y_distance = y_c - dike->y_Mc;
							M = dike->Mc + (dike->Mb - dike->Mc) * (y_distance / (back - dike->y_Mc));
							tempdikeRHS = M * 2 * v_spread / (right - left);
						}
						else
						{
							// linear interpolation between different M values, Mf is M in front, Mc acts as M in back
							y_distance = y_c - front;
							M = dike->Mf + (dike->Mc - dike->Mf) * (y_distance / (dike->y_Mc - front));
							tempdikeRHS = M * 2 * v_spread / (right - left);
						}
					}
					else if (dike->Mb != dike->Mf && dike->Mc < 0.0) // only Mf and Mb, they are different
					{
						if (jr->ctrl.var_M) // check varaible M option isn't used
						{
							SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid option: var_M option requires uniform M");
						}

						left = CurrPhTr->celly_xboundL[J];
						right = CurrPhTr->celly_xboundR[J];
						front = CurrPhTr->ybounds[0];
						back = CurrPhTr->ybounds[2 * nsegs - 1];
						v_spread = PetscAbs(bc->velin);

						// linear interpolation between different M values, Mf is M in front, Mb is M in back
						y_distance = y_c - front;
						M = dike->Mf + (dike->Mb - dike->Mf) * (y_distance / (back - front));
						tempdikeRHS = M * 2 * v_spread / (right - left);
					}
					else // Mb and Mf don't exist (which should not occurr)
					{
						SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No values for Mb and Mf. Dike option invalid!");
					}

					mat = &phases[i];

					// find dike heat source used (Behn and Ito, 2008)
					dike_contr = 0;
					if (Tc < mat->T_liq && Tc > mat->T_sol)
					{
						dike_contr = phRat[i] * (mat->rho * mat->Cp) * (mat->T_liq - Tc) * tempdikeRHS;
					}
					else if (Tc <= mat->T_sol)
					{
						dike_contr = phRat[i] * (mat->rho * mat->Cp) * ((mat->T_liq - Tc) + mat->Latent_hx / mat->Cp) * tempdikeRHS;
					}

					// compare heat sources to limit heatzone contribution where diking
					if (dike_contr >= hz_contr)
					{
						hz_contr = 0;
					}
					else
					{
						hz_contr -= dike_contr;
					}

				} // end if phRat and xboundR>xboundL
			}	  // close phase transition and phase ID comparison
		}		  // end dike block loop
	}			  // close phase transition block loop

	PetscFunctionReturn(0);
}