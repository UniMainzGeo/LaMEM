/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 **    filename:   matProps.c
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
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "solVar.h"
#include "fdstag.h"
#include "bc.h"
#include "JacRes.h"
#include "parsing.h"
#include "matProps.h"
#include "tools.h"
//---------------------------------------------------------------------------
// read material parameter limits
PetscErrorCode MatParLimRead(
		FB        *fb,
		Scaling   *scal,
		MatParLim *matLim)
{
/*
// scale gas constant with characteristic temperature
	matLim->Rugc        *= scal->temperature;
	matLim->rho_fluid   /= scal->density;
	matLim->rho_lithos  /= scal->density;
	matLim->theta_north /= scal->angle;

 	// viscosity limits
	PetscScalar eta_min;
	PetscScalar eta_max;
	// reference viscosity (initial guess)
	PetscScalar eta_ref;
	// reference temperature
	PetscScalar TRef;
	// universal gas constant
	PetscScalar Rugc;
	// viscosity & strain-rate tolerances
	PetscScalar eta_atol; // viscosity absolute tolerance
	PetscScalar eta_rtol; // viscosity relative tolerance
	PetscScalar DII_atol; // strain rate absolute tolerance
	PetscScalar DII_rtol; // strain rate relative tolerance
	// background (reference) strain-rate
	PetscScalar DII_ref;
	// plasticity parameters limits
	PetscScalar minCh;  // minimum cohesion
	PetscScalar minFr;  // maximum friction
	PetscScalar tauUlt; // ultimate yield stress
	// thermo-mechanical coupling controls
	PetscScalar shearHeatEff; // shear heating efficiency parameter [0 - 1]
	// rheology controls
	PetscBool   quasiHarmAvg; // quasi-harmonic averaging regularization flag (plasticity)
	PetscScalar cf_eta_min;   // visco-plastic regularization parameter (plasticity)
	PetscScalar n_pw;         // power-law regularization parameter (plasticity)
	PetscBool   initGuessFlg; // initial guess computation flag
	PetscBool   presLimFlg;   // pressure limit flag for plasticity
	PetscBool   presLimAct;   // activate pressure limit flag
	// fluid density for depth-dependent density model
	PetscScalar  rho_fluid;
	// rock density if we want to use lithostatic pressure in viscosit calculations
	PetscScalar  rho_lithos;
	// direction to the North for stress orientation
	// counter-clockwise positive measured from x-axis
	PetscScalar  theta_north;
	// print warning messages
	PetscBool    warn;
	// matrix-free closed-form jacobian
	PetscBool   jac_mat_free;
// ===
	// initialize material parameter limits
	PetscBool flg;
	PetscInt  cnt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

//	matLim->eta_min      = usr->LowerViscosityCutoff;
//	matLim->eta_max      = usr->UpperViscosityCutoff;
//	matLim->eta_ref      = usr->InitViscosity;

	matLim->TRef         = 0.0;
	matLim->Rugc         = 8.3144621;
	matLim->eta_atol     = 0.0;
	matLim->eta_rtol     = 1e-8;
	matLim->DII_atol     = 0.0;
	matLim->DII_rtol     = 1e-8;
	matLim->minCh        = 0.0;
	matLim->minFr        = 0.0;
	matLim->tauUlt       = DBL_MAX;
	matLim->shearHeatEff = 1.0;

	matLim->quasiHarmAvg = PETSC_FALSE;
	matLim->cf_eta_min   = 0.0;
	matLim->n_pw         = 0.0;

	matLim->initGuessFlg = PETSC_TRUE;
	matLim->rho_fluid    = 0.0;
	matLim->rho_lithos 	 = 0.0;	 // lithostatic density
	matLim->theta_north  = 90.0; // by default y-axis
	matLim->warn         = PETSC_TRUE;
	matLim->jac_mat_free = PETSC_FALSE;

//	if(usr->DII_ref) matLim->DII_ref = usr->DII_ref;
//	else
//	{
//		matLim->DII_ref = 1.0;
//		PetscPrintf(PETSC_COMM_WORLD," WARNING: Reference strain rate DII_ref is not defined. Use a non-dimensional reference value of DII_ref =%f \n", matLim->DII_ref);
//	}

	cnt = 0;

	// plasticity stabilization parameters
	ierr = PetscOptionsHasName(NULL, NULL, "-quasi_harmonic", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) { matLim->quasiHarmAvg = PETSC_TRUE; cnt++; }

	ierr = PetscOptionsGetScalar(NULL, NULL, "-cf_eta_min",  &matLim->cf_eta_min, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) { cnt++; }

	ierr = PetscOptionsGetScalar(NULL, NULL, "-n_pw",  &matLim->n_pw, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) { cnt++; }

	if(cnt > 1)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot combine plasticity stabilization methods (-quasi_harmonic -cf_eta_min -n_pw) \n");
	}

	// set Jacobian flag
	ierr = PetscOptionsHasName(NULL, NULL, "-jac_mat_free", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) matLim->jac_mat_free = PETSC_TRUE;

	if(cnt &&  matLim->jac_mat_free == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Analytical Jacobian is not available for plasticity stabilizations (-jac_mat_free -quasi_harmonic -cf_eta_min -n_pw) \n");
	}

	ierr = PetscOptionsGetScalar(NULL, NULL, "-rho_fluid",  &matLim->rho_fluid, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-rho_lithos", &matLim->rho_lithos, NULL); CHKERRQ(ierr);		// specify lithostatic density on commandline (if not set, we don't use this)

	ierr = PetscOptionsGetScalar(NULL, NULL, "-theta_north", &matLim->theta_north, NULL); CHKERRQ(ierr);

	ierr = PetscOptionsHasName(NULL, NULL, "-stop_warnings", &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE) matLim->warn = PETSC_FALSE;

	ierr = PetscOptionsGetScalar(NULL, NULL, "-shearHeatEff", &matLim->shearHeatEff, NULL); CHKERRQ(ierr);

	if(matLim->shearHeatEff > 1.0) matLim->shearHeatEff = 1.0;
	if(matLim->shearHeatEff < 0.0) matLim->shearHeatEff = 0.0;

*/

	/*

	// initialize
	matLim->eta_min      = 1e-8;
	matLim->eta_max      = 1e28;
	matLim->Rugc         = 8.3144621;
	matLim->eta_rtol     = 1e-8;
	matLim->DII_rtol     = 1e-8;
	matLim->minCh        = 1e-8;
	matLim->tauUlt       = 1e12;
	matLim->shearHeatEff = 1.0;
	matLim->init_guess   = PETSC_TRUE;
	matLim->plast_reg    = PETSC_FALSE;

	// read from options
	ierr = getScalarParam(fb, _OPTIONAL_, "eta_min",      &matLim->eta_min,      1, 1e-8,  1e28,  scal->viscosity  ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "eta_max",      &matLim->eta_max,      1, 1e-8,  1e28,  scal->viscosity  ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "TRef",         &matLim->TRef,         1, 0.0,   1e3,   scal->temperature); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "eta_atol",     &matLim->eta_atol,     1, 0.0,   1e5,   scal->viscosity  ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "eta_rtol",     &matLim->eta_rtol,     1, 1e-16, 1e-5,  1.0              ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "DII_atol",     &matLim->DII_atol,     1, 0.0,   1e-5,  scal->strain_rate); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "DII_rtol",     &matLim->DII_rtol,     1, 1e-16, 1e-5,  1.0              ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "DII_ref",      &matLim->DII_ref,      1, 1e-18, 1e-5,  scal->strain_rate); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "minCh",        &matLim->minCh,        1, 0.0,   1e12,  scal->stress     ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "minFr",        &matLim->minFr,        1, 0.0,   60.0,  scal->angle      ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "tauUlt",       &matLim->tauUlt,       1, 1e3,   1e12,  scal->stress     ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "shearHeatEff", &matLim->shearHeatEff, 1, 0.0,   1.0,   1.0              ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_fluid",    &matLim->rho_fluid,    1, 1e2,   1e4,   scal->density    ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "eta_init",     &matLim->eta_init,     1, 1e-8,  1e28,  scal->viscosity  ); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "plast_ratio",  &matLim->plast_ratio,  1, 5.0,   1e2,   1.0              ); CHKERRQ(ierr);

	// scale gas constant with characteristic temperature
	matLim->Rugc *= scal->temperature;

	// set viscosity regularization flag
	if(matLim->plast_ratio) matLim->plast_reg = PETSC_TRUE;

	 */
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropsReadAll"
PetscErrorCode MatPropsReadAll(
		FB         *fb,
		Scaling    *scal,
		PetscInt   *numPhases,
		Material_t *phases,
		PetscInt   *numSoft,
		Soft_t     *matSoft)
{
	// read all material phases and softening laws from file

	PetscInt jj;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//===============
	// SOFTENING LAWS
	//===============

	// print overview of softening laws from file
	PetscPrintf(PETSC_COMM_WORLD,"Reading softening laws: \n\n");

	// setup block access mode
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<SofteningStart>", "<SofteningEnd>"); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(jj = 0; jj < max_num_soft; jj++) matSoft[jj].ID = -1;

	// error checking
	if(fb->nblocks > max_num_soft)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many softening laws specified! Max allowed: %lld", (LLD)max_num_soft);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// read each individual softening law
	for(jj = 0; jj < fb->nblocks; jj++)
	{
		ierr = MatSoftRead(fb, fb->nblocks, matSoft); CHKERRQ(ierr);

		fb->blockID++;
	}

	// store actual number of softening laws
	(*numSoft) = fb->nblocks;

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	//================
	// MATERIAL PHASES
	//================

	// print overview of material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Reading material parameters: \n\n");

	// setup block access mode
	ierr = FBFindBlocks(fb, _REQUIRED_, "<MaterialStart>", "<MaterialEnd>"); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(jj = 0; jj < max_num_phases; jj++) phases[jj].ID = -1;

	// error checking
	if(fb->nblocks > max_num_phases)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)max_num_phases);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// read each individual phase
	for(jj = 0; jj < fb->nblocks; jj++)
	{
		ierr = MatPhaseRead(fb, scal, fb->nblocks, phases, (*numSoft), matSoft); CHKERRQ(ierr);

		fb->blockID++;

		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}

	// store actual number of phases
	(*numPhases) = fb->nblocks;

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPhaseRead"
PetscErrorCode MatPhaseRead(
		FB         *fb,
		Scaling    *scal,
		PetscInt    numPhases,
		Material_t *phases,
		PetscInt    numSoft,
		Soft_t     *matSoft)
{
	// read material properties from file with error checking

	Material_t *m;
	PetscScalar eta, eta0, e0;
	PetscInt    ID = -1, chSoftID, frSoftID, MSN;
	char        ndiff[MAX_NAME_LEN], ndisl[MAX_NAME_LEN], npeir[MAX_NAME_LEN];

	// output labels
	char        lbl_rho  [_lbl_sz_];
	char        lbl_eta  [_lbl_sz_];
	char        lbl_Bd   [_lbl_sz_];
	char        lbl_E    [_lbl_sz_];
	char        lbl_V    [_lbl_sz_];
	char        lbl_Bn   [_lbl_sz_];
	char        lbl_Bp   [_lbl_sz_];
	char        lbl_tau  [_lbl_sz_];
	char        lbl_fr   [_lbl_sz_];
	char        lbl_alpha[_lbl_sz_];
	char        lbl_cp   [_lbl_sz_];
	char        lbl_k    [_lbl_sz_];
	char        lbl_A    [_lbl_sz_];
    char        lbl_beta [_lbl_sz_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize additional parameters
	eta      =  0.0;
	eta0     =  0.0;
	e0       =  0.0;
	chSoftID = -1;
	frSoftID = -1;
	MSN      =  numSoft - 1;

	// phase ID
	ierr = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, numPhases-1); CHKERRQ(ierr);

	// get pointer to specified phase
	m = phases + ID;

	// check ID
	if(m->ID != -1)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Duplicate phase definition!");
	}

	// set ID
	m->ID = ID;

	//=================================================================================
	// density
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "rho0",     &m->rho,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_n",    &m->rho_n, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_c",    &m->rho_c, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "beta",     &m->beta,  1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// creep profiles
	//=================================================================================
	// set predefined diffusion creep profile
	ierr = GetProfileName(fb, scal, ndiff, "diff_profile");               CHKERRQ(ierr);
	ierr = SetDiffProfile(m, ndiff);                                      CHKERRQ(ierr);
	// set predefined dislocation creep profile
	ierr = GetProfileName(fb, scal, ndisl, "disl_profile");               CHKERRQ(ierr);
	ierr = SetDislProfile(m, ndisl);                                      CHKERRQ(ierr);
	// set predefined Peierls creep profile
	ierr = GetProfileName(fb, scal, npeir, "peir_profile");               CHKERRQ(ierr);
	ierr = SetPeirProfile(m, npeir);                                      CHKERRQ(ierr);
	//=================================================================================
	// Newtonian linear diffusion creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "eta",      &eta,      1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Bd",       &m->Bd,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Ed",       &m->Ed,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Vd",       &m->Vd,    1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// power-law (dislocation) creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "eta0",     &eta0,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "e0",       &e0,       1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Bn",       &m->Bn,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "n",        &m->n,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "En",       &m->En,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Vn",       &m->Vn,    1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// Peierls creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "Bp",       &m->Bp,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "taup",     &m->taup,  1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "gamma",    &m->gamma, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "q",        &m->q,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Ep",       &m->Ep,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Vp",       &m->Vp,    1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// elasticity
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "shear",    &m->G,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "bulk",     &m->K,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Kp",       &m->Kp,    1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// plasticity (Drucker-Prager)
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "cohesion", &m->ch,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "friction", &m->fr,    1, 1.0); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "chSoftID", &chSoftID, 1, MSN); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "frSoftID", &frSoftID, 1, MSN); CHKERRQ(ierr);
	//=================================================================================
	// energy
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "alpha",    &m->alpha, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "cp",       &m->Cp,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "k",        &m->k,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "A",        &m->A,     1, 1.0); CHKERRQ(ierr);
	//=================================================================================

	// check depth-dependent density parameters
	if((!m->rho_n && m->rho_c) || (m->rho_n && !m->rho_c))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "rho_n & rho_c must be specified simultaneously for phase %lld", (LLD)ID);
	}

	if(m->fr && !m->ch)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Nonzero cohesion must be specified for phase %lld", (LLD)ID);
	}

	// set pointers to softening laws
	if(chSoftID != -1) m->chSoft = matSoft + chSoftID;
	if(frSoftID != -1) m->frSoft = matSoft + frSoftID;

	// check strain-rate dependent creep
	if((!eta0 && e0) || (eta0 && !e0))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)ID);
	}

	// check power-law exponent
	if(!m->n && ((eta0 && e0) || m->Bn))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)ID);
	}

	// check Peierls creep
	if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)ID);
	}

	// recompute creep parameters
	if(eta)        m->Bd = 1.0/(2.0*eta);
	if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

	// check that at least one essential deformation mechanism is specified
	if(!m->Bd && !m->Bn && !m->G)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)ID);
	}

	// print
	if(scal->utype == _NONE_)
	{
		sprintf(lbl_rho,   "[ ]"         );
		sprintf(lbl_eta,   "[ ]"         );
		sprintf(lbl_Bd,    "[ ]"         );
		sprintf(lbl_E,     "[ ]"         );
		sprintf(lbl_V,     "[ ]"         );
		sprintf(lbl_Bn,    "[ ]"         );
		sprintf(lbl_Bp,    "[ ]"         );
		sprintf(lbl_tau,   "[ ]"         );
		sprintf(lbl_fr,    "[ ]"         );
		sprintf(lbl_alpha, "[ ]"         );
		sprintf(lbl_beta,  "[ ]"         );
		sprintf(lbl_cp,    "[ ]"         );
		sprintf(lbl_k,     "[ ]"         );
		sprintf(lbl_A,     "[ ]"         );
	}
	else
	{
		sprintf(lbl_rho,   "[kg/m3]"     );
		sprintf(lbl_eta,   "[Pa*s]"      );
		sprintf(lbl_Bd,    "[1/(Pa*s)]"  );
		sprintf(lbl_E,     "[J/mol]"     );
		sprintf(lbl_V,     "[m3/mol]"    );
		sprintf(lbl_Bn,    "[1/(Pa^n*s)]");
		sprintf(lbl_Bp,    "[1/s]"       );
		sprintf(lbl_tau,   "[Pa]"        );
		sprintf(lbl_fr,    "[deg]"       );
		sprintf(lbl_alpha, "[1/K]"       );
		sprintf(lbl_beta,  "[1/Pa]"      );
		sprintf(lbl_cp,    "[J/kg/K]"    );
		sprintf(lbl_k,     "[W/m/K]"     );
		sprintf(lbl_A,     "[W/m3]"      );
	}

	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: rho = %g %s, eta = %g %s, beta = %g %s\n", (LLD)(m->ID), m->rho, lbl_rho,  eta, lbl_eta, m->beta, lbl_beta);
	if (strlen(ndiff)) PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) diffusion creep profile: %s \n",(LLD)(m->ID), ndiff);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) Bd = %g %s, Ed = %g %s, Vd = %g %s \n", (LLD)(m->ID), m->Bd, lbl_Bd, m->Ed, lbl_E, m->Vd, lbl_V);
	if (strlen(ndisl)) PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) dislocation creep profile: %s \n",(LLD)(m->ID), ndisl);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) Bn = %g %s, En = %g %s, Vn = %g %s, n = %g [ ] \n", (LLD)(m->ID), m->Bn, lbl_Bn, m->En , lbl_E, m->Vn, lbl_V, m->n);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (peirl) Bp = %g %s, Ep = %g %s, Vp = %g %s, taup = %g %s, gamma = %g [ ], q = %g [ ] \n", (LLD)(m->ID), m->Bp, lbl_Bp, m->Ep, lbl_E, m->Vp, lbl_V, m->taup, lbl_tau, m->gamma, m->q);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (elast) G = %g %s, K = %g %s, Kp = %g [ ] \n", (LLD)(m->ID), m->G, lbl_tau, m->K, lbl_tau, m->Kp);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g %s, friction angle = %g %s \n", (LLD)(m->ID),m->ch, lbl_tau, m->fr, lbl_fr);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [ ], friction SoftLaw = %lld [ ] \n", (LLD)(m->ID),(LLD)chSoftID, (LLD)frSoftID);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g %s, cp = %g %s, k = %g %s, A = %g %s \n", (LLD)(m->ID),m->alpha, lbl_alpha, m->Cp, lbl_cp,m->k, lbl_k, m->A, lbl_A);
	PetscPrintf(PETSC_COMM_WORLD,"    \n");

	// scale
	// NOTE: [1] activation energy is not scaled
	//       [2] activation volume is multiplied with characteristic stress in SI units

	m->rho     /= scal->density;
	m->rho_c   *= scal->length_si;
	m->beta    *= scal->stress_si; // [1/Pa]

	// diffusion creep
	m->Bd      *= scal->viscosity;
	m->Vd      *= scal->stress_si;

	// dislocation creep (power-law)
	m->Bn      *= pow(scal->stress_si, m->n)*scal->time_si;
	m->Vn      *= scal->stress_si;

	// Peierls creep
	m->Bp      /=  scal->strain_rate;
	m->Vp      *=  scal->stress_si;
	m->taup    /=  scal->stress_si;

	// elasticity
	m->G       /= scal->stress_si;
	m->K       /= scal->stress_si;

	// plasticity
	m->ch      /= scal->stress_si;
	m->fr      /= scal->angle;

	// temperature
	m->alpha   /= scal->expansivity;
	m->Cp      /= scal->cpecific_heat;
	m->k       /= scal->conductivity;
	m->A       /= scal->heat_production;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatSoftRead"
PetscErrorCode MatSoftRead(
		FB       *fb,
		PetscInt  numSoft,
		Soft_t   *matSoft)
{
	// read softening law from file

	Soft_t   *s;
	PetscInt  ID;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// softening law ID
	ierr = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, numSoft-1); CHKERRQ(ierr);

	// get pointer to specified softening law
	s = matSoft + ID;

	// check ID
	if(s->ID != -1)
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Duplicate softening law!");
	}

	// set ID
	s->ID = ID;

	// read and store softening law parameters
	ierr = getScalarParam(fb, _REQUIRED_, "A",    &s->A,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "APS1", &s->APS1, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "APS2", &s->APS2, 1, 1.0); CHKERRQ(ierr);

	if(!s->A || !s->APS1 || !s->APS2)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "All parameters must be nonzero for softening law %lld", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"SoftLaw [%lld]: A = %g, APS1 = %g, APS2 = %g \n", (LLD)(s->ID), s->A, s->APS1, s->APS2);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//............ PREDEFINED RHEOLOGICAL PROFILES (from literature) ............
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetProfileName"
PetscErrorCode GetProfileName(FB *fb, Scaling *scal, char name[], const char key[])
{
	// read profile name from file

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscMemzero(name, sizeof(char)*MAX_NAME_LEN); CHKERRQ(ierr);

	ierr = getStringParam(fb, _OPTIONAL_, key, name, MAX_NAME_LEN);  CHKERRQ(ierr);

	if(strlen(name) && scal->utype == _NONE_)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Predefined creep profile is not supported for non-dimensional setup");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetDiffProfile"
PetscErrorCode SetDiffProfile(Material_t *m, char name[])
{
	// set diffusion creep profiles from literature

	TensorCorrection tensorCorrection;
	PetscInt         MPa;
	PetscScalar      d0, p;
	PetscScalar      C_OH_0, r;

	// We assume that the creep law has the form:
	// Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
	// Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
	//
	// In addition, we take into account that the creep-laws are typically measured under uniaxial or simple shear,
	// whereas we need them in tensorial format (tensorCorrection and F2) as defined in T. Gerya book.
	//
	// The resulting expressions for effective viscosity:
	// Diffusion:   inv_eta_diff = 2 * [Bd * exp(-(Ed + P*Vd)/(R*T))]
	// Dislocation: inv_eta_disl = 2 * [Bn * exp(-(En + P*Vn)/(R*T))]^(1/n) * eII^(1-1/n)
	//
	// In LaMEM we include the effect of grain size, H2O and tensor correction in the pre-factor (Bd,Bn) such that:
	// Diffusion:   Bd  = (2*F2)^(-1) * Ad [Pa] * d^-p * C_OH^r
	// Dislocation: Bn  = (2*F2)^(-n) * An [Pa]        * C_OH^r
	//
	//   eII     -   strain rate             [1/s]
	//   Tau     -   stress                  [Pa]
	//   P       -   pressure                [Pa]
	//   R       -   gas constant
	//   Ad, An  -   prefactor (Bn before taking into account grain size and water fugacity) [Pa^(-n)s^(-1)]
	//   Bd, Bn  -   prefactor               [Pa^(-n)s^(-1)]
	//   n       -   power-law exponent (n=1 for diffusion creep)
	//   Ed, En  -   activation Energy       [J/MPA/mol]
	//   Vd, Vn  -   activation volume       [m^3/mol]
	//   d       -   grain size              [in micro-meters (1e-6 meter)]
	//   p       -   exponent of grain size
	//   C_OH    -   water fugacity in H/10^6 Si  (see Hirth & Kohlstedt 2003 for a description)
	//   r       -   power-law exponent of C_OH term
	//   MPa     -   transform units: 0 - units in Pa; 1 - units in MPa

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check for empty string
	if(!strlen(name)) PetscFunctionReturn(0);

	if (!strcmp(name,"Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		m->Bd            =   1.5e9;
		m->Ed            =   375e3;
		m->Vd            =   5e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		d0               =   10e3;
		p                =   3;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003_constant_C_OH"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		m->Bd            =   1.0e6;
		m->Ed            =   335e3;
		m->Vd            =   4e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		d0               =   10e3;
		p                =   3;
		C_OH_0           =   1000;
		r                =   1;
	}

	else if (!strcmp(name,"Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		m->Bd            =   2.5e7;
		m->Ed            =   375e3;
		m->Vd            =   10e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;

		d0               =   10e3;
		p                =   3;
		C_OH_0           =   1000;
		r                =   0.8;
	}

	else
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "No such diffusion creep profile: %s! ",name);
	}

	// make tensor correction and transform units from MPa if necessary
	ierr = SetProfileCorrection(&m->Bd,1,tensorCorrection,MPa); CHKERRQ(ierr);

	// take into account grain size and water content
	m->Bd *= pow(d0,-p)*pow(C_OH_0,r);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetDislProfile"
PetscErrorCode SetDislProfile(Material_t *m, char name[])
{
	// set dislocation creep profiles from literature

	TensorCorrection tensorCorrection;
	PetscInt         MPa;
	PetscScalar      C_OH_0, r;

	// We assume that the creep law has the form:
	// Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
	// Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
	//
	// In addition, we take into account that the creep-laws are typically measured under uniaxial or simple shear,
	// whereas we need them in tensorial format (tensorCorrection and F2) as defined in T. Gerya book.
	//
	// The resulting expressions for effective viscosity:
	// Diffusion:   inv_eta_diff = 2 * [Bd * exp(-(Ed + P*Vd)/(R*T))]
	// Dislocation: inv_eta_disl = 2 * [Bn * exp(-(En + P*Vn)/(R*T))]^(1/n) * eII^(1-1/n)
	//
	// In LaMEM we include the effect of grain size (d,p), H2O fugacity (C_OH_0,r) and tensor correction (F2) in the pre-factor (Bd,Bn) such that:
	// Diffusion:   Bd  = (2*F2)^(-1) * Ad [Pa] * d^-p * C_OH^r
	// Dislocation: Bn  = (2*F2)^(-n) * An [Pa]        * C_OH^r
	//
	//   eII     -   strain rate             [1/s]
	//   Tau     -   stress                  [Pa]
	//   P       -   pressure                [Pa]
	//   R       -   gas constant
	//   Ad, An  -   prefactor (Bn before taking into account grain size and water fugacity) [Pa^(-n)s^(-1)]
	//   Bd, Bn  -   prefactor               [Pa^(-n)s^(-1)]
	//   n       -   power-law exponent (n=1 for diffusion creep)
	//   Ed, En  -   activation Energy       [J/MPA/mol]
	//   Vd, Vn  -   activation volume       [m^3/mol]
	//   d       -   grain size              [in micro-meters (1e-6 meter)]
	//   p       -   exponent of grain size
	//   C_OH    -   water fugacity in H/10^6 Si  (see Hirth & Kohlstedt 2003 for a description)
	//   r       -   power-law exponent of C_OH term
	//   MPa     -   transform units: 0 - units in Pa; 1 - units in MPa

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check for empty string
	if(!strlen(name)) PetscFunctionReturn(0);

	if (!strcmp(name,"Dry_Olivine-Ranalli_1995"))
	{
		// after Ranalli 1995
		m->Bn            =   2.5e4;
		m->n             =   3.5;
		m->En            =   532e3;
		m->Vn            =   17e-6;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Olivine-Ranalli_1995"))
	{
		// after Ranalli 1995
		m->Bn            =   2.0e3;
		m->n             =   4.0;
		m->En            =   471e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Quartz_Diorite-Hansen_Carter_1982"))
	{
		// taken from Carter and Tsenn (1986). Flow properties of continental lithosphere - page 18.
		m->Bn            =   pow(10,-1.5);
		m->n             =   2.4;
		m->En            =   212e3;
		m->Vn            =   0;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Diabase-Caristan_1982"))
	{
		// Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
		m->Bn            =   6e-2;
		m->n             =   3.05;
		m->En            =   276e3;
		m->Vn            =   1;
		tensorCorrection =   _UniAxial_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Tumut_Pond_Serpentinite-Raleigh_Paterson_1965"))
	{
		// Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
		m->Bn            =   6.3e-7;
		m->n             =   2.8;
		m->En            =   66e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Quarzite-Ranalli_1995"))
	{
		// used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
		// Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
		// Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
		// in Ranalli 1995 (page 334 Table 10.3)
		m->Bn            =   3.2e-4;
		m->n             =   2.3;
		m->En            =   154e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Quarzite-Ranalli_1995"))
	{
		// used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
		// Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
		// Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
		// in Ranalli 1995 (page 334 Table 10.3)
		m->Bn            =   6.7e-6;
		m->n             =   2.4;
		m->En            =   156e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Mafic_Granulite-Ranalli_1995"))
	{
		// used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
		// Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
		// Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
		// in Ranalli 1995 (page 334 Table 10.3)
		m->Bn            =   1.4e4;
		m->n             =   4.2;
		m->En            =   445e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Maryland_strong_diabase-Mackwell_et_al_1998"))
	{
		// Mackwell, Zimmerman & Kohlstedt (1998). High-temperature deformation
		// of dry diabase with application to tectonics on Venus. JGR 103. B1. 975-984. page 980
		m->Bn            =   8;
		m->n             =   4.7;
		m->En            =   485e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Quarzite-Ueda_et_al_2008"))
	{
		// Parameters used in Ueda et al (PEPI 2008)
		m->Bn            =   pow(10,-3.5);
		m->n             =   2.3;
		m->En            =   154e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Diabase-Huismans_et_al_2001"))
	{
		// parameters used in Huismans et al 2001
		m->Bn            =   3.2e-20;
		m->n             =   3.05;
		m->En            =   276e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Granite-Huismans_et_al_2001"))
	{
		// parameters used in Huismans et al 2001
		m->Bn            =   3.16e-26;
		m->n             =   3.3;
		m->En            =   186.5e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Dry_Upper_Crust-Schmalholz_Kaus_Burg_2009"))
	{
		// granite - Burg And Podladchikov (1999)
		m->Bn            =   3.16e-26;
		m->n             =   3.3;
		m->En            =   190e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Weak_Lower_Crust-Schmalholz_Kaus_Burg_2009"))
	{
		// diabase - Burg And Podladchikov (1999)
		m->Bn            =   3.2e-20;
		m->n             =   3.0;
		m->En            =   276e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Plagioclase_An75-Ranalli_1995"))
	{
		m->Bn            =   3.3e-4;
		m->n             =   3.2;
		m->En            =   238e3;
		m->Vn            =   0;
		tensorCorrection =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		// Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
		m->Bn            =   1600;
		m->n             =   3.5;
		m->En            =   520e3;
		m->Vn            =   22e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		C_OH_0           =   1000;
		r                =   1.2;
	}

	else if (!strcmp(name,"Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003_constant_C_OH"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		// Inside the subduction Factory 83?105. Table 1, "wet dislocation (constant C_OH)" parameters
		m->Bn            =   90;
		m->n             =   3.5;
		m->En            =   480e3;
		m->Vn            =   11e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		C_OH_0           =   1000;
		r                =   1.2;
	}

	else if (!strcmp(name,"Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		// Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
		m->Bn            =   1.1e5;
		m->n             =   3.5;
		m->En            =   530e3;
		m->Vn            =   15e-6;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Olivine-Burg_Podladchikov_1999"))
	{
		// after Burg and Podladchikov 1999
		m->Bn            =   7.1e-14;
		m->n             =   3.0;
		m->En            =   510e3;
		m->Vn            =   0;
		tensorCorrection =   _SimpleShear_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Upper_Mantle-Burg_Schmalholz_2008"))
	{
		// used in  SchmalholzKausBurg(2009), Geology (wet olivine)
		m->Bn            =   2e-21;
		m->n             =   4.0;
		m->En            =   471e3;
		m->Vn            =   0;
		tensorCorrection =   _SimpleShear_;
		MPa              =   0;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Granite-Tirel_et_al_2008"))
	{
		// used in  SchmalholzKausBurg(2009), Geology
		m->Bn            =   1.25e-9;
		m->n             =   3.2;
		m->En            =   123e3;
		m->Vn            =   0;
		tensorCorrection =   _SimpleShear_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}
	
	else if (!strcmp(name,"Ara_rocksalt-Urai_et_al.(2008)"))
    {
        // Ara rocksalt as published in Urai et al.(2008)
        m->Bn            =   1.82e-9;
        m->n             =   5;
        m->En            =   32.4e3;
        m->Vn            =   0;
        tensorCorrection =   _UniAxial_;
        MPa              =   1;
        C_OH_0           =   1;
        r                =   0;
    }

    else if (!strcmp(name,"Polycrystalline_Anhydrite-Mueller_and_Briegel(1978)"))
    {
        //
        m->Bn            =   3.16228e1;
        m->n             =   2;
        m->En            =   152.3e3;
        m->Vn            =   0;
        tensorCorrection =   _UniAxial_;
        MPa              =   1;
        C_OH_0           =   1;
        r                =   0;
    }
    
	else
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "No such dislocation creep profile: %s! ",name);
	}

	// make tensor correction and transform units from MPa if necessary
	ierr = SetProfileCorrection(&m->Bn,m->n,tensorCorrection,MPa); CHKERRQ(ierr);

	// take into account grain size and water content
	m->Bn *= pow(C_OH_0,r);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetPeirProfile"
PetscErrorCode SetPeirProfile(Material_t *m, char name[])
{
	// set Peierls creep profiles from literature

	// We assume that the creep law has the form:
	// Peierls:   eII = Bp * exp( - (EP + P*VP)/(R*T)*(1-gamma)^q) * (Tau/gamma/taup)^s
	//            s   = (Ep+p*Vp)/(R*T)*(1-gamma)^(q-1)*q*gamma
	//
	// where:
	// Bp         - pre-exponential constant for the Peierls mechanism [1/s]
	// Ep         - activation energy [J/mol K]
	// Vp         - activation volume [m3/mol ]
	// taup       - Peierl stress [Pa]
	// gamma      - adjustable constant [-]
	// q          - stress dependence for Peierls creep [-]
	// s          - Peierls creep exponent (typical values between 7-11) [-]

	PetscFunctionBegin;

	// check for empty string
	if(!strlen(name)) PetscFunctionReturn(0);

	if (!strcmp(name,"Olivine_Peierls-Kameyama_1999"))
	{
		// used in Kameyama et al 1999 (EPSL), vol 168., pp. 159-172
		// original source: Guyot and Dorn (1967) and Poirier (1985)
		m->Bp            = 5.7e11;
		m->Ep            = 5.4e5;
		m->Vp            = 0.0;
		m->taup          = 8.5e9;
		m->gamma         = 0.1;
		m->q             = 2;
	}

	else
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "No such Peierls creep profile: %s! ",name);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetProfileCorrection"
PetscErrorCode SetProfileCorrection(PetscScalar *B, PetscScalar n, TensorCorrection tensorCorrection, PetscInt MPa)
{
	// set tensor and units correction for rheological profiles

	PetscScalar F2, Bi;
	// Lab. experiments are typically done under simple shear or uni-axial
	// compression, which requires a correction in order to use them in tensorial format.
	// An explanation is given in the textbook of Taras Gerya, chapter 6, p. 71-78.

	PetscFunctionBegin;

	Bi = *B;

	// Tensor correction
	// In LaMEM this is added to the pre-factor and not to the effective viscosity as in T. Gerya
	if      (tensorCorrection == _UniAxial_)    F2 = pow(0.5,(n-1)/n) / pow(3,(n+1)/(2*n)); //  F2 = 1/2^((n-1)/n)/3^((n+1)/2/n);
	else if (tensorCorrection == _SimpleShear_) F2 = pow(0.5,(2*n-1)/n);                    //  F2 = 1/2^((2*n-1)/n);
	else if (tensorCorrection == _None_)        F2 = 0.5;
	else
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Unknown tensor correction in creep mechanism profile!");
	}

	// Units correction from [MPa^(-n)s^(-1)] to [Pa^(-n)s^(-1)] if required
	if (MPa) Bi = pow(2*F2,-n) * pow(1e6*pow(Bi,-1/n),-n);
	else     Bi = pow(2*F2,-n) * Bi;

	(*B) = Bi;

	PetscFunctionReturn(0);
}
/*
//---------------------------------------------------------------------------
// This needs to be updated for the use in the inversion routines
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropSetFromLibCall"
PetscErrorCode MatPropSetFromLibCall(JacRes *jr, ModParam *mod)
{
	// overwrite MATERIAL PARAMETERS with model parameters provided by a calling function

	PetscInt 	id,im;
	PetscScalar eta, eta0, e0;
	Material_t  *m;

	PetscFunctionBegin;

	if(mod == NULL) PetscFunctionReturn(0);

	// does a calling function provide model parameters?
	if(mod->use == 0) PetscFunctionReturn(0);

	// set material properties
	if(mod->use == 1) {
		PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"# Material properties set from calling function: \n");


		for(im=0;im<mod->mdN;im++)
		{

			id = mod->phs[im];
			// get pointer to specified phase
			m = jr->phases + id;

			// linear viscosity
			if(mod->typ[im] == _ETA_)
			{

				// initialize additional parameters
				eta      =  0.0;
				eta0     =  0.0;
				e0       =  0.0;
				eta = mod->val[im];

				// check strain-rate dependent creep
				if((!eta0 && e0) || (eta0 && !e0))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)id);
				}

				// check power-law exponent
				if(!m->n && ((eta0 && e0) || m->Bn))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)id);
				}

				// check Peierls creep
				if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)id);
				}

				// recompute creep parameters
				if(eta)        m->Bd = 1.0/(2.0*eta);
				if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

				// check that at least one essential deformation mechanism is specified
				if(!m->Bd && !m->Bn && !m->G)
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)id);
				}

				PetscPrintf(PETSC_COMM_WORLD,"#    eta[%lld] = %g \n",(LLD)id,eta);
			}

			// constant density
			else if(mod->typ[im] == _RHO0_)
			{
				m->rho = mod->val[im];
				PetscPrintf(PETSC_COMM_WORLD,"#    rho0[%lld] = %5.5f \n",(LLD)id,m->rho);
			}

			else
			{
				PetscPrintf(PETSC_COMM_WORLD,"WARNING: inversion parameter type is not implemented \n");
			}

		}
		PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropSetFromCL"
PetscErrorCode MatPropSetFromCL(JacRes *jr)
{
	// overwrite MATERIAL PARAMETERS with command line options

	PetscErrorCode 	ierr;
	PetscBool		flg,get_options;
	PetscInt 		id;
	char 			matprop_opt[MAX_PATH_LEN];
	PetscScalar eta, eta0, e0;
	Material_t *m;

	PetscFunctionBegin;

	flg = PETSC_FALSE;
	get_options = PETSC_FALSE;

	ierr = PetscOptionsGetBool(NULL, NULL, "-SetMaterialProperties", &get_options, NULL ); 					CHKERRQ(ierr);

	if(get_options) {
		PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"# Material properties set from command line: \n");

		for(id=0;id<jr->numPhases;id++){

			// get pointer to specified phase
			m = jr->phases + id;

			// initialize additional parameters
			eta      =  0.0;
			eta0     =  0.0;
			e0       =  0.0;

			// linear viscosity
			sprintf(matprop_opt,"-eta_%lld",(LLD)id);
			ierr = PetscOptionsGetReal(NULL, NULL ,matprop_opt,&eta	, &flg); 				CHKERRQ(ierr);

				// check strain-rate dependent creep
				if((!eta0 && e0) || (eta0 && !e0))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)id);
				}

				// check power-law exponent
				if(!m->n && ((eta0 && e0) || m->Bn))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)id);
				}

				// check Peierls creep
				if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)id);
				}

				// recompute creep parameters
				if(eta)        m->Bd = 1.0/(2.0*eta);
				if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

				// check that at least one essential deformation mechanism is specified
				if(!m->Bd && !m->Bn && !m->G)
				{
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)id);
				}

			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    eta[%lld]	= %g \n",(LLD)id,eta);

			// constant density
			sprintf(matprop_opt,"-rho0_%lld",(LLD)id);
			ierr = PetscOptionsGetReal(NULL, NULL ,matprop_opt,&m->rho	, &flg);			CHKERRQ(ierr);
			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    rho0[%lld]	= %5.5f \n",(LLD)id,m->rho);

		}
		PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
