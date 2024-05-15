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

// added for file output *djking
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

//---------------------------------------------------------------------------
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{

	// read all dike parameter blocks from file
	Dike *dike;
	FDSTAG *fs;
	PetscScalar ***gsxx_eff_ave_hist, ***raw_gsxx_ave_hist, ***smooth_gsxx_ave_hist;
	PetscInt jj, nD, numDike, numdyndike, istep_nave;
	PetscInt i, j, istep_count, sx, sy, sisc, nx, ny;

	PetscFunctionBeginUser;

	if (!jr->ctrl.actDike)
		PetscFunctionReturn(0); // only execute this function if dikes are active

	fs = jr->fs;
	//===============
	// DIKE PARAMETER
	//===============

	// setup block access mode
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<DikeStart>", "<DikeEnd>"));

	if (fb->nblocks)
	{
		// print overview of dike blocks from file
		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "Dike blocks : \n");

		// initialize ID for consistency checks

		for (jj = 0; jj < _max_num_dike_; jj++)
			dbdike->matDike[jj].ID = -1;

		// error checking
		if (fb->nblocks > _max_num_dike_)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many dikes specified! Max allowed: %lld", (LLD)_max_num_dike_);

		// store actual number of dike blocks
		dbdike->numDike = fb->nblocks;

		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

		// read each individual dike block
		for (jj = 0; jj < fb->nblocks; jj++)
		{
			PetscCall(DBReadDike(dbdike, dbm, fb, jr, PrintOutput));
			fb->blockID++;
		}

		if (PrintOutput)
			PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	}

	PetscCall(FBFreeBlocks(fb));

	numdyndike = 0;

	numDike = dbdike->numDike;
	for (nD = 0; nD < numDike; nD++) // loop through all dike blocks
	{
		dike = dbdike->matDike + nD;
		if (dike->dyndike_start)
		{
			numdyndike++;
			if (numdyndike == 1) //(take this out of this loop because it will be repeated with >1 dynamic dike)
			{
				// DM for 1D cell center vector. vector is ny+1 long, but for whatever reason that must appear in the first,
				// not second entry on line 2 below
				PetscCall(DMDACreate3dSetUp(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
											fs->dsy.tnods, fs->dsy.nproc, fs->dsz.nproc,
											fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc, 1, 1,
											0, 0, 0, &jr->DA_CELL_1D));

				// DM for 2D cell center vector, with istep_nave planes for time averaging
				PetscCall(DMDACreate3dSetUp(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
											fs->dsx.tcels, fs->dsy.tcels, fs->dsz.nproc * dike->istep_nave,
											fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc, 1, 1,
											0, 0, 0, &jr->DA_CELL_2D_tave));
			}

			// creating local vectors and inializing the history vector
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->magPressure));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->focused_magPressure));

			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->sxx_eff_ave));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->raw_sxx));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->raw_sxx_ave));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->smooth_sxx));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->smooth_sxx_ave));

			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->solidus)); // *djking
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D, &dike->magPresence)); // *djking
			
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D_tave, &dike->sxx_eff_ave_hist));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D_tave, &dike->raw_sxx_ave_hist));
			PetscCall(DMCreateLocalVector(jr->DA_CELL_2D_tave, &dike->smooth_sxx_ave_hist));
			
			PetscCall(DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist));
			PetscCall(DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->raw_sxx_ave_hist, &raw_gsxx_ave_hist));
			PetscCall(DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->smooth_sxx_ave_hist, &smooth_gsxx_ave_hist));

			PetscCall(DMDAGetCorners(jr->DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, &istep_nave));

			for (j = sy; j < sy + ny; j++)
			{
				for (i = sx; i < sx + nx; i++)
				{
					for (istep_count = sisc; istep_count < sisc + istep_nave; istep_count++)
					{
						gsxx_eff_ave_hist[istep_count][j][i] = 0.0;
						raw_gsxx_ave_hist[istep_count][j][i] = 0.0;
						smooth_gsxx_ave_hist[istep_count][j][i] = 0.0;
					}
				}
			}

			PetscCall(DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist));
			PetscCall(DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->raw_sxx_ave_hist, &raw_gsxx_ave_hist));
			PetscCall(DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->smooth_sxx_ave_hist, &smooth_gsxx_ave_hist));
		} // End if dyndike->start
	}	  // End loop through dikes

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput)
{
	// read dike parameter from file
	Dike *dike;
	PetscInt ID;
	Scaling *scal;

	PetscFunctionBeginUser;

	// access context
	scal = dbm->scal;

	// Dike ID
	PetscCall(getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbdike->numDike - 1));
	fb->ID = ID;

	// get pointer to specified dike parameters
	dike = dbdike->matDike + ID;

	// check ID
	if (dike->ID != -1)
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
	PetscCall(getScalarParam(fb, _REQUIRED_, "Mf", &dike->Mf, 1, 1.0));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "Mc", &dike->Mc, 1, 1.0));
	PetscCall(getScalarParam(fb, _REQUIRED_, "Mb", &dike->Mb, 1, 1.0));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "y_Mc", &dike->y_Mc, 1, scal->length));
	PetscCall(getIntParam(fb, _REQUIRED_, "PhaseID", &dike->PhaseID, 1, dbm->numPhases - 1));
	PetscCall(getIntParam(fb, _REQUIRED_, "PhaseTransID", &dike->PhaseTransID, 1, dbm->numPhtr - 1));
	PetscCall(getIntParam(fb, _OPTIONAL_, "dyndike_start", &dike->dyndike_start, 1, -1));

	// variable M parameters
	if (jr->ctrl.var_M)
	{

		dike->A = 1e3;		 // default smoothing
		dike->Ts = 10e6;	 // default tensile strength
		dike->zeta_0 = 1e22; // default reference bulk viscosity

		PetscCall(getScalarParam(fb, _OPTIONAL_, "A", &dike->A, 1, scal->stress_si));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "Ts", &dike->Ts, 1, scal->stress_si));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "zeta_0", &dike->zeta_0, 1, scal->viscosity));
	}

	// parameters for average lithospheric stress calculations (includes magma pressure)
	if (dike->dyndike_start || jr->ctrl.var_M || jr->ctrl.sol_track)
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
		dike->magPfac = 1.0;
		dike->magPwidth = 1e+30;
		// dike->ymindyn=-1e+30;
		// dike->ymaxdyn=1e+30;

		PetscCall(getScalarParam(fb, _OPTIONAL_, "Tsol", &dike->Tsol, 1, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "filtx", &dike->filtx, 1, scal->length));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "filty", &dike->filty, 1, scal->length));

		PetscCall(getScalarParam(fb, _OPTIONAL_, "zmax_magma", &dike->zmax_magma, 1, scal->length));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "drhomagma", &dike->drhomagma, 1, scal->density));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "magPfac", &dike->magPfac, 1, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "magPwidth", &dike->magPwidth, 1, scal->length));
		// PetscCall(getScalarParam(fb, _OPTIONAL_, "ymindyn",	&dike->ymindyn,		1, 1.0));
		// PetscCall(getScalarParam(fb, _OPTIONAL_, "ymaxdyn",	&dike->ymaxdyn,		1, 1.0));

		PetscCall(getIntParam(fb, _OPTIONAL_, "istep_nave", &dike->istep_nave, 1, 50));
		PetscCall(getIntParam(fb, _OPTIONAL_, "nstep_locate", &dike->nstep_locate, 1, 1000));
		PetscCall(getIntParam(fb, _OPTIONAL_, "out_stress", &dike->out_stress, 1, 50));
		PetscCall(getIntParam(fb, _OPTIONAL_, "out_dikeloc", &dike->out_dikeloc, 1, 50));

		// scaling
		dike->Tsol = (dike->Tsol + jr->scal->Tshift) / jr->scal->temperature;

		// initialize so that when istep=dike_start, it is set to 0
		dike->istep_count = dike->istep_nave;
	}

	// print diking info to terminal
	if (PrintOutput)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Dike [%lld]: ", (LLD)(dike->ID));
		if (dike->dyndike_start)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Dynamic starting at timestep = %lld\n", (LLD)(dike->dyndike_start));
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Static \n");
		}
		PetscPrintf(PETSC_COMM_WORLD, "   PhaseTransID=%lld PhaseID=%lld Mf=%g, Mb=%g, Mc=%g, y_Mc=%g \n",
					(LLD)(dike->PhaseTransID), (LLD)(dike->PhaseID), dike->Mf, dike->Mb, dike->Mc, dike->y_Mc);
		if (jr->ctrl.var_M)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   Variable M option used:\n");
			PetscPrintf(PETSC_COMM_WORLD, "     A = %lld %s, Ts = %lld %s, zeta_0 = %g %s\n",
						(LLD)(dike->A * scal->stress), scal->lbl_stress, (LLD)(dike->Ts * scal->stress), scal->lbl_stress, dike->zeta_0 * scal->viscosity, scal->lbl_viscosity);
		}
		if (dike->dyndike_start || jr->ctrl.var_M)
		{
			PetscPrintf(PETSC_COMM_WORLD, "   gsxx_eff_ave smoothing parameters:\n");
			PetscPrintf(PETSC_COMM_WORLD, "     Tsol = %1.0f %s, filtx = %1.2f %s, filty = %1.2f %s\n",
						dike->Tsol * scal->temperature - scal->Tshift, scal->lbl_temperature, dike->filtx * scal->length, scal->lbl_length, dike->filty * scal->length, scal->lbl_length);
			PetscPrintf(PETSC_COMM_WORLD, "     nstep_locate = %lld, istep_nave = %lld, istep_count = %lld\n",
						(LLD)(dike->nstep_locate), (LLD)(dike->istep_nave), (LLD)(dike->istep_count));
			PetscPrintf(PETSC_COMM_WORLD, "   excess magma pressure options:\n");
			PetscPrintf(PETSC_COMM_WORLD, "     zmax_magma = %1.0f %s, magPwidth = %1.0f %s\n",
						dike->zmax_magma*scal->length, scal->lbl_length, dike->magPwidth*scal->length, scal->lbl_length);
			PetscPrintf(PETSC_COMM_WORLD, "     drhomagma = %1.0f %s, magPfac = %1.1f\n",
						dike->drhomagma*scal->density, scal->lbl_density, dike->magPfac);
		}
	}

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------
PetscErrorCode GetDikeContr(JacRes *jr,
							PetscScalar *phRat, // phase ratios in the control volume
							PetscInt &AirPhase,
							PetscScalar &dikeRHS,
							PetscScalar &y_c,
							PetscScalar &z_c,
							PetscInt J,
							PetscScalar sxx_eff_ave_cell,
							PetscScalar zsolidus) // *revisit (add PetscInt I if more than 1 dike?)

{

	BCCtx *bc;
	Dike *dike;
	Ph_trans_t *CurrPhTr;
	PetscInt i, nD, nPtr, numDike, numPhtr, nsegs;
	PetscScalar v_spread, M, left, right, front, back;
	PetscScalar y_distance, tempdikeRHS;
	PetscScalar P_comp, div_max, M_rat, zeta;

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
				// if (phRat[i] > 0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J] && z_c >= zsolidus) // solidus *djking
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
								zeta = dike->A * (dike->zeta_0 / P_comp) + P_comp / div_max;
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

					// Divergence
					dikeRHS += (phRat[i] + phRat[AirPhase]) * tempdikeRHS; // Give full divergence if cell is part dike part air (*revisit Why??)

				} // close if phRat and xboundR>xboundL
			}	  // close phase transition and dike phase ID comparison
		}		  // close dike block loop
	}			  // close phase transition block loop
	PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------------------------------------------
PetscErrorCode Dike_k_heatsource(JacRes *jr,
								 Material_t *phases,
								 PetscScalar &Tc,
								 PetscScalar *phRat,
								 PetscScalar &k,
								 PetscScalar &rho_A,
								 PetscScalar &y_c,
								 PetscScalar &z_c,
								 PetscInt J,
								 PetscScalar sxx_eff_ave_cell,
								 PetscScalar zsolidus)

{
	// parameters to determine dilation term
	BCCtx *bc;
	Dike *dike;
	Ph_trans_t *CurrPhTr;
	PetscInt i, nD, nPtr, numDike, numPhtr, nsegs;
	PetscScalar v_spread, M, left, right, front, back;
	PetscScalar y_distance, tempdikeRHS;
	PetscScalar P_comp, div_max, M_rat, zeta;

	// heating parameters
	Material_t *mat;
	PetscScalar kfac;

	PetscFunctionBeginUser;

	numDike = jr->dbdike->numDike; // number of dikes
	numPhtr = jr->dbm->numPhtr;
	bc = jr->bc;

	nPtr = 0;
	nD = 0;
	kfac = 0;

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
				if (phRat[i] > 0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J]) // solidus *djking
				// if (phRat[i] > 0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J] && z_c >= zsolidus) // solidus *djking
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
								zeta = dike->A * (dike->zeta_0 / P_comp) + P_comp / div_max;
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
							SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid option: var_M option requires single global M (i.e. Mf = Mb)");
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
							SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid option: var_M option requires single global M (i.e. Mf = Mb)");
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

					// adjust k and heat source according to Behn & Ito [2008]
					if (Tc < mat->T_liq && Tc > mat->T_sol)
					{
						kfac += phRat[i] / (1 + (mat->Latent_hx / (mat->Cp * (mat->T_liq - mat->T_sol))));
						rho_A += phRat[i] * (mat->rho * mat->Cp) * (mat->T_liq - Tc) * tempdikeRHS; // Cp*rho not used in the paper, added to conserve units of rho_A
					}
					else if (Tc <= mat->T_sol)
					{
						kfac += phRat[i];
						rho_A += phRat[i] * (mat->rho * mat->Cp) * ((mat->T_liq - Tc) + mat->Latent_hx / mat->Cp) * tempdikeRHS;
					}
					else if (Tc >= mat->T_liq)
					{
						kfac += phRat[i];
					}
					// end adjust k and heat source according to Behn & Ito [2008]

					k = kfac * k;

				} // end if phRat and xboundR>xboundL
			}	  // close phase transition and phase ID comparison
		}		  // end dike block loop
	}			  // close phase transition block loop

	PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------------------------------------------
PetscErrorCode Locate_Dike_Zones(AdvCtx *actx)
{

	Controls *ctrl;
	JacRes *jr;
	Dike *dike;
	Ph_trans_t *CurrPhTr;
	FDSTAG *fs;
	PetscInt nD, numDike, numPhtr, nPtr, n, icounter;
	PetscInt i, j, j1, j2, sx, sy, sz, ny, nx, nz;
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	jr = actx->jr;
	fs = jr->fs;
	ctrl = &jr->ctrl;

	if (!ctrl->actDike || jr->ts->istep + 1 == 0) PetscFunctionReturn(0); // only execute this function if dikes are active
	fs = jr->fs;

	PetscPrintf(PETSC_COMM_WORLD, "\n");
	numDike = jr->dbdike->numDike; // number of dikes
	numPhtr = jr->dbm->numPhtr;

	icounter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	if (ctrl->actDike)
	{
		for (nD = 0; nD < numDike; nD++)
		{
			dike = jr->dbdike->matDike + nD; // sets dike to point to the nD-th element of matDike array

			if ((dike->dyndike_start && (jr->ts->istep + 1 >= dike->dyndike_start)) || (jr->ctrl.var_M))
			{
				PetscPrintf(PETSC_COMM_WORLD, "Locating Dike zone: istep=%lld dike # %lld\n", (LLD)(jr->ts->istep + 1), (LLD)(nD));
				// compute lithostatic pressure
				if (icounter == 0)
				{
					ierr = JacResGetLithoStaticPressure(jr); CHKERRQ(ierr);
					ierr = ADVInterpMarkToCell(actx); CHKERRQ(ierr);
				}
				icounter++;

				//---------------------------------------------------------------------------------------------
				//  Find dike phase transition
				//---------------------------------------------------------------------------------------------
				nPtr = -1;
				for (n = 0; n < numPhtr; n++)
				{
					CurrPhTr = jr->dbm->matPhtr + n;
					if (CurrPhTr->ID == dike->PhaseTransID)
					{
						nPtr = n;
					}
				} // end loop over Phtr

				if (nPtr == -1)
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "PhaseTransID problems with dike %lld, nPtr=%lld\n", (LLD)(nD), (LLD)(nPtr));

				CurrPhTr = jr->dbm->matPhtr + nPtr;

				//---------------------------------------------------------------------------------------------
				//  Find y-bounds of current dynamic dike
				//---------------------------------------------------------------------------------------------
				j1 = ny - 1;
				j2 = 0;
				for (j = 0; j < ny; j++)
				{
					if (CurrPhTr->celly_xboundR[j] > CurrPhTr->celly_xboundL[j])
					{
						j1 = (PetscInt)min(j1, j);
						j2 = (PetscInt)max(j2, j);
					}
				}

				ierr = Compute_sxx_magP(jr, nD); CHKERRQ(ierr); // compute mean effective sxx across the lithosphere

				ierr = Smooth_sxx_eff(jr, nD, nPtr, j1, j2); CHKERRQ(ierr); // smooth mean effective sxx

				// Only relocate dike zone if dynamic diking is on and if on an nstep_locate timestep
				if (dike->dyndike_start && (jr->ts->istep + 1 >= dike->dyndike_start) && ((jr->ts->istep + 1) % dike->nstep_locate) == 0)
				{
					ierr = Set_dike_zones(jr, nD, nPtr, j1, j2); CHKERRQ(ierr); // centered on peak sxx_eff_ave
				}

				// change z boundary of dike box if trackSolidus is set
				if (jr->ctrl.sol_track) // *djking
				{
					ierr = Set_dike_base(jr, nD, nPtr, j1, j2); CHKERRQ(ierr); // use solidus values to set zbound of dike
				}
			}
			PetscFunctionReturn(0);
		}
	}
	else if (jr->ctrl.sol_track) // gets solidus array (as well as average sxx, etc...)
	{
		nD = 0;
		ierr = Compute_sxx_magP(jr, nD);
		CHKERRQ(ierr); // compute mean effective sxx across the lithosphere
	}
}
//------------------------------------------------------------------------------------------------------------------

PetscErrorCode Compute_sxx_magP(JacRes *jr, PetscInt nD)
{
  MPI_Request srequest, rrequest;
  Vec         vsxx, vPmag, vliththick, vzsol; 
  PetscScalar ***magPressure, ***focused_magPressure;
  PetscScalar ***gsxx_eff_ave, ***p_lith;
  PetscScalar ***raw_gsxx, ***smooth_gsxx;
  PetscScalar ***raw_gsxx_ave, ***smooth_gsxx_ave;
  PetscScalar ***sxx, ***Pmag, ***liththick, ***zsol;
  PetscScalar ***solidus, ***magPresence; // *djking
  PetscScalar zsol_max_local = -PETSC_MAX_REAL,  zsol_max_global; // *djking
  PetscScalar  *lsxx, *lPmag, *lliththick, *lzsol;
  PetscScalar dz, ***lT, Tc, *grav, Tsol, dPmag, magma_presence;
  PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, L, ID, AirPhase;
  PetscMPIInt rank;


  FDSTAG      *fs;
  Dike        *dike;
  Discret1D   *dsz;
  SolVarCell  *svCell;
  Controls    *ctrl;

  PetscErrorCode ierr;
  PetscFunctionBeginUser;

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
  
  ///// INTEGRATE LITHSPHERIC SXX, LITHOSTATIC PRESSURE (Pmag), AND THICKNESS /////
  Tsol=dike->Tsol;
  for(k = sz + nz - 1; k >= sz; k--)
  {
     dz  = SIZE_CELL(k, sz, (*dsz));
     START_PLANE_LOOP
          GET_CELL_ID(ID, i-sx, j-sy, k-sz, nx, ny);  //GET_CELL_ID needs local indices
          svCell = &jr->svCell[ID]; 
          Tc=lT[k][j][i];
 
          if ((Tc<=Tsol) && (svCell->phRat[AirPhase] < 1.0))
          {
            dz  = SIZE_CELL(k, sz, (*dsz));
            sxx[L][j][i]+=(svCell->hxx - svCell->svBulk.pn)*dz;  //integrating dz-weighted total stress (*revisit deviatoric??)
            Pmag[L][j][i]+=p_lith[k][j][i]*dz; // removed from magPressure due to varM dike initiation criteria *djking

            liththick[L][j][i]+=dz;             //integrating thickness
          }
          
          //interpolate depth to the solidus
          if ((Tc <= Tsol) && (Tsol < lT[k-1][j][i]))
          {
            zsol[L][j][i]=dsz->ccoor[k-sz]+(dsz->ccoor[k-sz-1]-dsz->ccoor[k-sz])/(lT[k-1][j][i]-Tc)*(Tsol-Tc); 
          }
      END_PLANE_LOOP
  } 

      //After integrating thickness and dz-weighted total stress, send it down to the next proc. 
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

  //now all cores in z have the same solution so give that to the stress array
  
  ///// CALCULATE AVERAGE LITHOSPHERIC VALUES /////

  // (gdev is the array that shares data with devxx_mean and is indexed with global dimensions)
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->focused_magPressure, &focused_magPressure); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->raw_sxx, &raw_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->raw_sxx_ave, &raw_gsxx_ave); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->smooth_sxx, &smooth_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->smooth_sxx_ave, &smooth_gsxx_ave); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr); // *djking
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPresence, &magPresence); CHKERRQ(ierr); // *djking

	// store solidus in dike structure and find solidus max *djking
	START_PLANE_LOOP
	solidus[L][j][i] = zsol[L][j][i]; // store zsol in the solidus array outside function
	zsol_max_local = PetscMax(zsol[L][j][i], zsol_max_local); // finding local max solidus (thinnest lithosphere)
	END_PLANE_LOOP
	MPI_Allreduce(&zsol_max_local, &zsol_max_global, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); // find solidus global max

/* 	// to find 2 peaks if we switch to using focused_magPressure *djking testing
	PetscScalar zsol_max_local[2] = {-PETSC_MAX_REAL, -PETSC_MAX_REAL}; // Array to store local top two maxima
    PetscScalar zsol_max_global[2];
	    START_PLANE_LOOP
        solidus[L][j][i] = zsol[L][j][i];
        if (zsol[L][j][i] > zsol_max_local[0]) {
            zsol_max_local[1] = zsol_max_local[0];
            zsol_max_local[0] = zsol[L][j][i];
        } else if (zsol[L][j][i] > zsol_max_local[1]) {
            zsol_max_local[1] = zsol[L][j][i];
        }
    END_PLANE_LOOP
	MPI_Allreduce(&zsol_max_local, &zsol_max_global, 2, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD); */

	// calculate depth average stress (sxx) and excess magma pressure
	START_PLANE_LOOP

		dPmag = 0; // set dPmag to zero
		magma_presence= 0;  //testing
		if (dike->zmax_magma - zsol[L][j][i] < 0) // if negative, then postive magma pressure at solidus exists
		{
			dPmag = (dike->zmax_magma - zsol[L][j][i]) * (dike->drhomagma) * grav[2]; // excess static magma pressure at solidus
			magma_presence=dike->magPfac*(zsol[L][j][i]-dike->zmax_magma)/(zsol_max_global-dike->zmax_magma);  //undergoing testing
		}
		//magPressure[L][j][i] = (Pmag[L][j][i]/liththick[L][j][i]+dPmag)*magma_presence;
		magPressure[L][j][i] = dPmag;	// removed average lithostatic pressure *djking
		focused_magPressure[L][j][i] = dPmag * magma_presence;	// revisit and turn on magma_presence *djking
		magPresence[L][j][i] = magma_presence; // testing *djking
	
		gsxx_eff_ave[L][j][i] = sxx[L][j][i] / liththick[L][j][i]; // Depth weighted mean total stress
		raw_gsxx[L][j][i] = sxx[L][j][i] / liththick[L][j][i]; // Depth weighted mean total stress
		raw_gsxx_ave[L][j][i] = sxx[L][j][i] / liththick[L][j][i]; // Depth weighted mean total stress
		smooth_gsxx[L][j][i] = sxx[L][j][i] / liththick[L][j][i]; // Depth weighted mean total stress
		smooth_gsxx_ave[L][j][i] = sxx[L][j][i] / liththick[L][j][i]; // Depth weighted mean total stress
		
	END_PLANE_LOOP

/* 	// output mean stress array to .txt file on timesteps of standard output *djking
	if (((istep % nstep_out) == 0 || istep == 1) && (dike->out_stress > 0))
	{
		if (L == 0)
		{
			// Form the filename based on jr->ts->istep+1
			std::ostringstream oss;
			oss << "gsxx_Timestep_" << (jr->ts->istep + 1) << ".txt";

			std::string filename = oss.str();

			// Open a file with the formed filename
			std::ofstream outFile(filename);
			if (outFile)
			{
				START_PLANE_LOOP
				xcell = COORD_CELL(i, sx, fs->dsx);
				ycell = COORD_CELL(j, sy, fs->dsy);

				// Writing space delimited data
				outFile
					<< " " << xcell << " " << ycell
					<< " " << gsxx_eff_ave[L][j][i]
					<< " " << magPressure[L][j][i]
					<< " " << nD << " " << zsol[L][j][i]
					<< " " << magma_presence << "\n";

				END_PLANE_LOOP
			}
			else
			{
				std::cerr << "Error creating file: " << filename << std::endl;
			}
		}
	} */

  // restore buffer and mean stress vectors
  ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT);  CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->focused_magPressure, &focused_magPressure); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->raw_sxx, &raw_gsxx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->raw_sxx_ave, &raw_gsxx_ave); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->smooth_sxx, &smooth_gsxx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->smooth_sxx_ave, &smooth_gsxx_ave); CHKERRQ(ierr); 

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr); // *djking
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPresence, &magPresence); CHKERRQ(ierr); // *djking

  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vsxx, &sxx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vPmag, &Pmag); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vliththick, &liththick); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vzsol, &zsol); CHKERRQ(ierr);

  ierr = VecRestoreArray(vsxx, &lsxx); CHKERRQ(ierr);
  ierr = VecRestoreArray(vPmag, &lPmag); CHKERRQ(ierr);
  ierr = VecRestoreArray(vliththick, &lliththick); CHKERRQ(ierr);
  ierr = VecRestoreArray(vzsol, &lzsol); CHKERRQ(ierr);

  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vsxx); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vPmag); CHKERRQ(ierr);  
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vliththick); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vzsol); CHKERRQ(ierr);


  //fill ghost points

  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->magPressure);
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->focused_magPressure);
      
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->sxx_eff_ave);
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->raw_sxx);
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->raw_sxx_ave);
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->smooth_sxx);
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->smooth_sxx_ave);

  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->solidus); // *djking
  LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->magPresence); // *djking
//  ierr = VecGhostUpdateBegin(dike->solidus, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(dike->solidus, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

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

	PetscScalar ***magPressure, ***focused_magPressure;
	PetscScalar ***solidus, ***magPresence; // *djking
	PetscScalar ***gsxx_eff_ave, ***gsxx_eff_ave_hist;
	PetscScalar ***raw_gsxx, ***smooth_gsxx;
	PetscScalar ***raw_gsxx_ave, ***raw_gsxx_ave_hist;
	PetscScalar ***smooth_gsxx_ave, ***smooth_gsxx_ave_hist;
	PetscScalar ***ycoors, *lycoors, ***ycoors_prev, *lycoors_prev, ***ycoors_next, *lycoors_next;
	PetscScalar ***xcenter, *lxcenter, ***xcenter_prev, *lxcenter_prev, ***xcenter_next, *lxcenter_next;
	PetscScalar ***sxx, *lsxx, ***sxx_prev, *lsxx_prev, ***sxx_next, *lsxx_next;
	PetscScalar ***magP, *lmagP, ***magP_prev, *lmagP_prev, ***magP_next, *lmagP_next;
	PetscScalar xc, yc, xx, yy, dx, dy, sum_sxx, sum_magP, sum_w;
	PetscScalar filtx, filty, w, dfac, magPfac, magPwidth;
	PetscScalar xcent, xcent_north, xcent_south, ycent_north, ycent_south, xcent_search, ycent_search;
	PetscScalar azim, dalong, dxazim, dyazim, radbound, sumslope, sumadd;
	PetscScalar dx_tot, dy_tot, str_y;
	//PetscScalar dyazmin, dyazmax, dyaz;

	Vec         vycoors, vycoors_prev, vycoors_next;
	Vec         vxcenter, vxcenter_prev, vxcenter_next;
	Vec         vsxx, vsxx_prev, vsxx_next;
	Vec         vmagP, vmagP_prev, vmagP_next;

	PetscInt    j, jj, j1prev, j2prev, j1next, j2next, jj1, jj2; 
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

	dike = jr->dbdike->matDike+nD;
	filtx=dike->filtx;
	filty=dike->filty;
	dfac=1.0; //maximum distance for Gaussian weights is dfac*filtx and dfac*filty

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

//for dike center
	ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vxcenter); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vxcenter_prev); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(jr->DA_CELL_1D, &vxcenter_next); CHKERRQ(ierr);

	ierr = VecZeroEntries(vxcenter); CHKERRQ(ierr);
	ierr = VecZeroEntries(vxcenter_prev); CHKERRQ(ierr);
	ierr = VecZeroEntries(vxcenter_next); CHKERRQ(ierr);

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
//dike center info
	ierr = DMDAVecGetArray(jr->DA_CELL_1D, vxcenter, &xcenter); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_1D, vxcenter_prev, &xcenter_prev); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_1D, vxcenter_next, &xcenter_next); CHKERRQ(ierr);
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
	ierr = VecGetArray(vxcenter, &lxcenter); CHKERRQ(ierr);
	ierr = VecGetArray(vxcenter_prev, &lxcenter_prev); CHKERRQ(ierr);
	ierr = VecGetArray(vxcenter_next, &lxcenter_next); CHKERRQ(ierr);
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
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->raw_sxx, &raw_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->raw_sxx_ave, &raw_gsxx_ave); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->smooth_sxx, &smooth_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->smooth_sxx_ave, &smooth_gsxx_ave); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->focused_magPressure, &focused_magPressure); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr); // *djking
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->magPresence, &magPresence); CHKERRQ(ierr); // *djking
  
	START_PLANE_LOOP

		magP[L][j][i]=magPressure[L][j][i];
		sxx[L][j][i]=gsxx_eff_ave[L][j][i];

	END_PLANE_LOOP
  
//  Set up y-node coord and dike center arrays for passing between procs
	for(j = 0; j <= ny; j++)
	{
		xcenter[L][M][j]=1e+12;
		ycoors[L][M][j]=COORD_NODE(j+sy,sy,fs->dsy);  //can put j in last entry because ny<nx
	} 
//Dike center is given only on the current dike, i.e., j=j1 to j2
	for(j = j1; j <=j2; j++)
	{
		xcenter[L][M][j]=(CurrPhTr->celly_xboundR[j] + CurrPhTr->celly_xboundL[j])/2;    
	}

//--------------------------------------------------
// passing arrays between previous and next y proc
//--------------------------------------------------
	if (dsy->nproc > 1 && dsy->grprev != -1)  //Exchange arrays from previous proc, if not the first proc
	{
		ierr = MPI_Irecv(lycoors_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lxcenter_prev, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lsxx_prev, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lmagP_prev, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &rrequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);


		ierr = MPI_Isend(lycoors, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lxcenter, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grprev, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
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
		ierr = MPI_Isend(lxcenter, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest2); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest2, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lsxx, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest3); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest3, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Isend(lmagP, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &srequest4); CHKERRQ(ierr);
		ierr = MPI_Wait(&srequest4, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);

		ierr = MPI_Irecv(lycoors_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest); CHKERRQ(ierr);
		ierr = MPI_Wait(&rrequest, MPI_STATUSES_IGNORE);  CHKERRQ(ierr);
		ierr = MPI_Irecv(lxcenter_next, (PetscMPIInt)(ny+1), MPIU_SCALAR, dsy->grnext, 0, PETSC_COMM_WORLD, &rrequest2); CHKERRQ(ierr);
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
		//of filty of current dike point (xcent, yc)
		xcent=xcenter[L][M][j-sy];
		yc = COORD_CELL(j, sy, fs->dsy);
		sumslope=0;
		sumadd = 0; 
		//loop through full y domain to find all points of dike near (xcent,yc)
		for(jj = sy; jj < sy+ny; jj++)
		{
			//Current proc
			xcent_search=xcenter[L][M][jj-sy];  //beyond dike end this will be 1e12 so dalong>filty
			ycent_search=COORD_CELL(jj, sy, fs->dsy);
			dalong=sqrt(pow((xcent-xcent_search),2)+pow((yc-ycent_search),2));
			if (jj<j && dalong<=filty) 			//if south of current point
			{
				xcent_north=xcenter[L][M][jj-sy+1];
				ycent_north=COORD_CELL(jj+1, sy, fs->dsy);
				xcent_south=xcent_search;
				ycent_south=ycent_search;
				sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
				sumadd += 1;
			}
			else if (jj >j && dalong<=filty)  //if north of current point
			{
				xcent_north=xcent_search;
				ycent_north=ycent_search;
				xcent_south=xcenter[L][M][jj-sy-1];
				ycent_south=COORD_CELL(jj-1, sy, fs->dsy);
				sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
				sumadd += 1;
			}
			
			//NEXT proc
			if ( dsy->grnext != -1)
			{
				xcent_search=xcenter_next[L][M][jj-sy];  //if beyond dike end this will be 1e12 so dalong>filty
 				ycent_search=(ycoors_next[L][M][jj-sy+1]+ycoors_next[L][M][jj-sy])/2;
				dalong=sqrt(pow((xcent-xcent_search),2)+pow((yc-ycent_search),2));
				if (jj==sy && dalong<=filty) 			//if at southernmost cell of next proc
				{
					xcent_north=xcent_search;   
					ycent_north=ycent_search;
					xcent_south=xcenter[L][M][ny-1];  	//northernmost point of current proc (local index)
					ycent_south=COORD_CELL(ny+sy-1, sy, fs->dsy);  //uses global indexing
					sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
					sumadd += 1;
				}
				if (jj > sy && dalong<=filty)   		 //if north of southernmost cell of next proc
				{
					xcent_north=xcent_search;
					ycent_north=ycent_search;
					xcent_south=xcenter_next[L][M][jj-sy-1];
					ycent_south=(ycoors_next[L][M][jj-sy]+ycoors_next[L][M][jj-sy-1])/2;
					sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
					sumadd += 1;
				}
			}

			//Previous proc
			if ( dsy->grprev != -1)
			{
				xcent_search=xcenter_prev[L][M][jj-sy];  //if beyond dike end this will be 1e12 and dalong>filty
 				ycent_search=(ycoors_prev[L][M][jj-sy+1]+ycoors_prev[L][M][jj-sy])/2;
				dalong=sqrt(pow((xcent-xcent_search),2)+pow((yc-ycent_search),2));

				if (jj==sy+ny-1 && dalong<=filty) 		//if at northern most cell of prev proc
				{
					xcent_north=xcenter[L][M][0];   	//southernmost cell of current proc (local indexing)
					ycent_north=COORD_CELL(sy, sy, fs->dsy);  //uses global indexing
					xcent_south=xcent_search;  
					ycent_south=ycent_search;
					sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
					sumadd += 1;
				}
				if (jj < sy+ny-1 && dalong<=filty)    	//if south of northernmost cell of prev proc
				{
					xcent_north=xcenter_prev[L][M][jj-sy+1];  
					ycent_north=(ycoors_prev[L][M][jj-sy+2]+ycoors_prev[L][M][jj-sy+1])/2;
					xcent_south=xcent_search;  //northern most cell of prev proc
					ycent_south=ycent_search;
					sumslope += (xcent_north-xcent_south)/(ycent_north-ycent_south);
					sumadd += 1;
				}
			}
			
		} //done with loop over j to get mean azimuth
		azim=atan(sumslope/sumadd);

		//identify global y index of cells within dy_tot of yc on local and adjacent processors
		j1prev=ny+sy-1; j2prev=sy;
		j1next=ny+sy-1; j2next=sy;
		jj1=sy+ny-1; jj2=sy;
		//projected from slanted axis coords to get x & y grid distances needed to encompass dfac*filtx and dfac*filty
		dx_tot=(fabs(dfac*filtx*cos(azim))+fabs(dfac*filty*sin(azim)));
		dy_tot=(fabs(dfac*filtx*sin(azim))+fabs(dfac*filty*cos(azim)));  


		//dyazmin=1e6; dyazmax=-1e6;  //for detecting if near dike zone end
		//Loop over y to define area of Gaussian smoothing patch
		for(jj = sy; jj < sy+ny; jj++)
		{
			//Previous proc
 			yy=(ycoors_prev[L][M][jj-sy+1]+ycoors_prev[L][M][jj-sy])/2;
			if ( dsy->grprev != -1 && fabs(yc-yy) <= dy_tot && xcenter_prev[L][M][jj-sy] < 1.0e+12) 
			{
				j1prev=(PetscInt)min(j1prev,jj);   
				j2prev=(PetscInt)max(j2prev,jj);
			}
			/* dyaz=(yy-yc)/cos(azim);  //for stretching: if distance oriented with "azim" is within filty 
			if ( dsy->grprev != -1 && fabs(dyaz) <= filty && xcenter_prev[L][M][jj-sy] < 1.0e+12) 
			{
				dyazmin=(PetscScalar)min(dyaz,dyazmin);
			}*/

			//Next proc
			yy=(ycoors_next[L][M][jj-sy+1]+ycoors_next[L][M][jj-sy])/2;
			if (dsy->grnext != -1 && fabs(yy-yc) <= dy_tot && xcenter_next[L][M][jj-sy] < 1.0e+12)
			{
				j1next=(PetscInt)min(j1next,jj);   
				j2next=(PetscInt)max(j2next,jj);
			}
			/*dyaz=(yy-yc)/cos(azim);  //for stretching: if distance oriented with "azim" is within filty
			if (dsy->grnext != -1 && fabs(dyaz)<=filty && xcenter_next[L][M][jj-sy] < 1.0e+12)
			{
				dyazmax=(PetscScalar)max(dyaz,dyazmax);
			}*/

			//Current proc
			yy=COORD_CELL(jj, sy, fs->dsy);
			if (fabs(yy-yc) <= dy_tot && xcenter[L][M][jj-sy] < 1.0e+12)
			{
				jj1=(PetscInt)min(jj1,jj);
				jj2=(PetscInt)max(jj2,jj);
			}
			/*dyaz=(yy-yc)/cos(azim);  //for stretching: if distance oriented with "azim" is within filty 
			if (fabs(dyaz) <= filty && xcenter[L][M][jj-sy] < 1.0e+12)
			{
				dyazmin=(PetscScalar)min(dyaz,dyazmin);
				dyazmax=(PetscScalar)max(dyaz,dyazmax);
			}*/
		}  //end y loop for defining area of Gaussian smoothing patch

		str_y=1;
		//if ((dyazmax-dyazmin)<2*filty)       //if dike zone end limits the distance to < dfac*filty north or south
		//	str_y=2*filty/(dyazmax-dyazmin);  //then stretch filty smoothing extends a total distance 2*dfac*filty
		//str_y=(PetscScalar)min(str_y,2.0);
		/*if (L==0)  //debugging
		{ 
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"212121.2121 %lld M=%i: %i %g prev: %g to %g, curr: %g to %g next: %g to %g dytot=%g \n", 
			(LLD)(jr->ts->istep+1), M, j, yc, 
			(ycoors_prev[L][M][j1prev-sy+1]+ycoors_prev[L][M][j1prev-sy])/2, (ycoors_prev[L][M][j2prev-sy+1]+ycoors_prev[L][M][j2prev-sy])/2,
												COORD_CELL(jj1, sy, fs->dsy), COORD_CELL(jj2, sy, fs->dsy),
			(ycoors_next[L][M][j1next-sy+1]+ycoors_next[L][M][j1next-sy])/2,(ycoors_next[L][M][j2next-sy+1]+ycoors_next[L][M][j2next-sy])/2, dy_tot);

			//PetscSynchronizedPrintf(PETSC_COMM_WORLD,"2121.2121 %lld M=%i: j=%i, %g, %g \n", (LLD)(jr->ts->istep+1), M, j, ycoors_prev[L][M][j-sy], ycoors_prev[L][M][j-sy+1]); 
		}          //debugging
		*/

		//Loop over i to assign filtered value in cell j,i (again, one proc across all x dimension)
		for (i = sx; i < sx+nx; i++)  
		{
			sum_sxx=0.0;
			sum_magP=0.0;
			sum_w=0.0;

			xc =  COORD_CELL(i, sx, fs->dsx);
      
			//identify x cells within dfac*filtx of xc
			ii1=sx+nx-1; ii2=sx;
			for (ii = sx; ii < sx+nx; ii++)
			{
				xx = COORD_CELL(ii, sx, fs->dsx);
				if (fabs(xx-xc) <= dx_tot)
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
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/(str_y*filty)),2)))*dx*dy;
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
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/(str_y*filty)),2)))*dx*dy;
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
						w=exp(-0.5*(pow((dxazim/filtx),2) + pow((dyazim/(str_y*filty)),2)))*dx*dy;
						sum_sxx += sxx_next[L][jj][ii]*w;
						sum_magP += magP_next[L][jj][ii]*w;
						sum_w+=w;
					}
				}
			} //end loop over cells from next proc

			//sum_w=max(sum_w,0.0);  //why would sum_w be <0???!
			magPressure[L][j][i]=(sum_magP/sum_w);
			focused_magPressure[L][j][i]=(sum_magP/sum_w)*magPfac*exp(-0.5*(pow((cos(azim)*(xcent-xc)/magPwidth),2)));

			smooth_gsxx[L][j][i]=(sum_sxx/sum_w);
			smooth_gsxx_ave[L][j][i]=(sum_sxx/sum_w);
			gsxx_eff_ave[L][j][i]=(sum_sxx/sum_w) + magPressure[L][j][i];
			//gsxx_eff_ave[L][j][i]=(sum_sxx/sum_w) + focused_magPressure[L][j][i]; // testing

		}//End loop over i
	}// End loop over j

  	//PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT); // debugging All procs must run this

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

	//restore xcenter arrays
	ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vxcenter_prev, &xcenter_prev); CHKERRQ(ierr);
	ierr = VecRestoreArray(vxcenter_prev, &lxcenter_prev); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vxcenter_prev); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vxcenter, &xcenter); CHKERRQ(ierr);
	ierr = VecRestoreArray(vxcenter, &lxcenter); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vxcenter); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_1D, vxcenter_next, &xcenter_next); CHKERRQ(ierr);
	ierr = VecRestoreArray(vxcenter_next, &lxcenter_next); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_1D, &vxcenter_next); CHKERRQ(ierr);

	//restore magP arrays
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP_prev, &magP_prev); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP_prev, &lmagP_prev); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP_prev); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP, &magP); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP, &lmagP); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, vmagP_next, &magP_next); CHKERRQ(ierr);
	ierr = VecRestoreArray(vmagP_next, &lmagP_next); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(jr->DA_CELL_2D, &vmagP_next); CHKERRQ(ierr);

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
		ierr = DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->raw_sxx_ave_hist, &raw_gsxx_ave_hist); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(jr->DA_CELL_2D_tave, dike->smooth_sxx_ave_hist, &smooth_gsxx_ave_hist); CHKERRQ(ierr);

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
				raw_gsxx_ave_hist[sisc+dike->istep_count][j][i]=raw_gsxx_ave[L][j][i];  //array for current step
				smooth_gsxx_ave_hist[sisc+dike->istep_count][j][i]=smooth_gsxx_ave[L][j][i];  //array for current step

				for (istep_count=sisc; istep_count<sisc+istep_nave; istep_count++)
				{
 					sum_sxx+=gsxx_eff_ave_hist[istep_count][j][i];
 					sum_sxx+=raw_gsxx_ave_hist[istep_count][j][i];
 					sum_sxx+=smooth_gsxx_ave_hist[istep_count][j][i];
				}

				gsxx_eff_ave[L][j][i]=sum_sxx/((PetscScalar)istep_nave);
				raw_gsxx_ave[L][j][i]=sum_sxx/((PetscScalar)istep_nave);
				smooth_gsxx_ave[L][j][i]=sum_sxx/((PetscScalar)istep_nave);
			}
		}

		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->sxx_eff_ave_hist, &gsxx_eff_ave_hist); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->raw_sxx_ave_hist, &raw_gsxx_ave_hist); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D_tave, dike->smooth_sxx_ave_hist, &smooth_gsxx_ave_hist); CHKERRQ(ierr);
	}// end if nstep_ave>1

  // output smoothed stress array to .txt file on timesteps of other output
  if (((istep % nstep_out) == 0 || istep == 1) && (dike->out_stress > 0)) 
  {
    if (L == 0)
    {
      // Form the filename based on jr->ts->istep+1
      std::ostringstream oss;
      oss << "sxx_outputs_Timestep_" << std::setfill('0') << std::setw(8) << (jr->ts->istep+1) << ".txt";
      std::string filename = oss.str();

      // Open a file with the formed filename
      std::ofstream outFile(filename);
      if (outFile)
      {
        START_PLANE_LOOP
        xc = COORD_CELL(i, sx, fs->dsx);
        yc = COORD_CELL(j, sy, fs->dsy);

        // Writing space delimited data
        outFile
          << " " << xc << " " << yc 
          << " " << magPressure[L][j][i] 
          << " " << focused_magPressure[L][j][i] 
          << " " << raw_gsxx[L][j][i] 
          << " " << raw_gsxx_ave[L][j][i] 
          << " " << smooth_gsxx[L][j][i] 
          << " " << smooth_gsxx_ave[L][j][i] 
          << " " << gsxx_eff_ave[L][j][i] 
          << " " << jr->ts->istep+1 << " " << jr->ts->time * jr->scal->time    
          << " " << solidus[L][j][i] << " " << magPresence[L][j][i] << "\n";    

        END_PLANE_LOOP
      }
      else
      {
        std::cerr << "Error creating file: " << filename << std::endl;
      }
    }
  }  

	//restore arrays
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPressure, &magPressure); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->focused_magPressure, &focused_magPressure); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->raw_sxx, &raw_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->raw_sxx_ave, &raw_gsxx_ave); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->smooth_sxx, &smooth_gsxx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->smooth_sxx_ave, &smooth_gsxx_ave); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr); // *djking
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->magPresence, &magPresence); CHKERRQ(ierr); // *djking

	LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->sxx_eff_ave);
	LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->solidus); // *djking
	LOCAL_TO_LOCAL(jr->DA_CELL_2D, dike->magPresence); // *djking

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

 // dike location to .txt file on timesteps of other output
    if (L==0 &&  ((istep % nstep_out) == 0 || istep == 1) && (dike->out_dikeloc > 0))
    {
      // Form the filename based on jr->ts->istep
      std::ostringstream oss;
      oss << "dikeloc_Timestep_" << (jr->ts->istep+1) << ".txt";

      std::string filename = oss.str();

      // Create/open file
      std::ofstream outFile(filename);
      if (outFile)
      {
        ycell = COORD_CELL(j, sy, fs->dsy);
        xcell=(COORD_CELL(ixmax-1, sx, fs->dsx)+COORD_CELL(ixmax, sx, fs->dsx))/2;

        // Writing space delimited data
        outFile
          << " " << xcell << " " << ycell 
          << " " << xcenter << " " << xshift << " " << x_maxsxx 
          << " " << COORD_CELL(ixmax, sx, fs->dsx) 
          << " " << CurrPhTr->celly_xboundL[lj] 
          << " " << CurrPhTr->celly_xboundR[lj] 
          << " " << nD << " " << dtime << "\n"; 
      }
      else
      {
        std::cerr << "Error creating file: " << filename << std::endl;
      }
    }
  } //end loop over j cell row

	if (((istep % nstep_out)==0 || istep == 1) && (dike->out_dikeloc > 0))  
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
//----------------------------------------------------------------------------------------------------
// Set bottom bounds of NotInAir box based on solidus
// NOTE: this only tested with cpu_x = 1
//

PetscErrorCode Set_dike_base(JacRes *jr, PetscInt nD, PetscInt nPtr, PetscInt j1, PetscInt j2)
{

	FDSTAG      *fs;
	Dike        *dike;
	Discret1D   *dsz;
	Ph_trans_t  *CurrPhTr;
	PetscScalar ldikeSolidus = PETSC_MAX_REAL;
	PetscScalar gdikeSolidus, ***solidus, xc;
	PetscInt    sx, nx, sy, ny, sz, nz, i, j, L;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  =  jr->fs;
	dsz = &fs->dsz;
	L   =  (PetscInt)dsz->rank;

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	dike = jr->dbdike->matDike+nD;
	CurrPhTr = jr->dbm->matPhtr+nPtr;

	// get solidus array
	ierr = DMDAVecGetArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr);

	// find local minimum solidus in dike zone nD
	for (j = j1; j <= j2; j++)
	{
		for (i = sx; i < sx + nx; i++)
		{
			xc =  COORD_CELL(i, sx, fs->dsx);
			if (xc >= CurrPhTr->celly_xboundL[j] && xc <= CurrPhTr->celly_xboundR[j])
			{
				ldikeSolidus = PetscMin(ldikeSolidus, solidus[L][j][i]);
			}
		}
	}
	
	// find the global minimum across all processors (should be the same across all z-procs anyways)
	ierr = MPI_Allreduce(&ldikeSolidus, &gdikeSolidus, 1, MPIU_SCALAR, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);

	// restore solidus array
	ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dike->solidus, &solidus); CHKERRQ(ierr);

	// set zbounds[0] so that divergence only occurs in brittle lithosphere
	CurrPhTr->zbounds[0] = gdikeSolidus;

	// output *djking
	if (L == 0)
	{
		std::ostringstream oss;
		oss << "dikeSolidus_Timestep_" << std::setfill('0') << std::setw(8) << (jr->ts->istep + 1) << ".txt";
		std::string filename = oss.str();

		// open output file
		::ofstream outFile(filename, std::ios_base::app); // append if already created
		if (outFile)
		{
			// write data
			outFile
				<< "dike " << nD << " " << gdikeSolidus << " " << CurrPhTr->zbounds[0] << "\n";
		}
		else
		{
			std::cerr << "Error opening file: " << filename << std::endl;
		}
	}

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
