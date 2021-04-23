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
 **    filename:   phase.c
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
 **			Andrea Piccolo
 **         Georg Reuber
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/*
	This file defined various material properties for the phases

*/

//---------------------------------------------------------------------------
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "phase.h"
#include "parsing.h"
#include "scaling.h"
#include "objFunct.h"
#include "JacRes.h"
#include "phase_transition.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBMatCreate"
PetscErrorCode DBMatCreate(DBMat *dbm, FB *fb, PetscBool PrintOutput)
{
	// read all material phases and softening laws from file

	PetscInt jj;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//===============
	// SOFTENING LAWS
	//===============

	// setup block access mode
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<SofteningStart>", "<SofteningEnd>"); CHKERRQ(ierr);

	if(fb->nblocks)
	{
		// print overview of softening laws from file
		if (PrintOutput){
			PetscPrintf(PETSC_COMM_WORLD,"Softening laws: \n");
		}
		// initialize ID for consistency checks
		for(jj = 0; jj < _max_num_soft_; jj++) dbm->matSoft[jj].ID = -1;

		// error checking
		if(fb->nblocks > _max_num_soft_)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many softening laws specified! Max allowed: %lld", (LLD)_max_num_soft_);
		}

		// store actual number of softening laws
		dbm->numSoft = fb->nblocks;

		if (PrintOutput){
			PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
		}
		// read each individual softening law
		for(jj = 0; jj < fb->nblocks; jj++)
		{
			ierr = DBMatReadSoft(dbm, fb, PrintOutput); CHKERRQ(ierr);

			fb->blockID++;
		}
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);



	//================
	// MATERIAL PHASES
	//================
	if (PrintOutput){
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

		// print overview of material parameters read from file
		PetscPrintf(PETSC_COMM_WORLD,"Material parameters: \n");

		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}

	// setup block access mode
	ierr = FBFindBlocks(fb, _REQUIRED_, "<MaterialStart>", "<MaterialEnd>"); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(jj = 0; jj < _max_num_phases_; jj++) dbm->phases[jj].ID = -1;

	// error checking
	if(fb->nblocks > _max_num_phases_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)_max_num_phases_);
	}

	// store actual number of phases
	dbm->numPhases = fb->nblocks;

	// read each individual phase
	for(jj = 0; jj < fb->nblocks; jj++)
	{
		ierr = DBMatReadPhase(dbm, fb, PrintOutput); CHKERRQ(ierr);

		fb->blockID++;

	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);
	if (PrintOutput){
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}

	// setup block access mode
		ierr = FBFindBlocks(fb, _OPTIONAL_, "<PhaseTransitionStart>", "<PhaseTransitionEnd>"); CHKERRQ(ierr);

		if(fb->nblocks)
		{
			// print overview of softening laws from file
			PetscPrintf(PETSC_COMM_WORLD,"Phase Transition laws: \n");

			// initialize ID for consistency checks
			for(jj = 0; jj < _max_num_soft_; jj++) dbm->matPhtr[jj].ID = -1;

			// error checking
			if(fb->nblocks > _max_num_tr_)
			{
				SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many phase_transition specified! Max allowed: %lld", (LLD)_max_num_tr_);
			}

			// store actual number of Phase Transition laws
			dbm->numPhtr = fb->nblocks;

			PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

			// read each individual softening law
			for(jj = 0; jj < fb->nblocks; jj++)
			{
				ierr = DBMatReadPhaseTr(dbm, fb); CHKERRQ(ierr);

				fb->blockID++;
			}

			// adjust density if needed
			ierr = Overwrite_density(dbm);CHKERRQ(ierr);
		
		}
		ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	

    //=================================================
	// OVERWRITE MATERIAL PARAMETERS WITH GLOBAL VARIABLES
	//=================================================
    ierr = DBMatOverwriteWithGlobalVariables(dbm, fb); CHKERRQ(ierr);

	if (PrintOutput){
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBMatReadSoft"
PetscErrorCode DBMatReadSoft(DBMat *dbm, FB *fb, PetscBool PrintOutput)
{
	// read softening law from file
	Scaling  *scal;
	Soft_t   *s;
	PetscInt  ID;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	scal = dbm->scal;

	// softening law ID
	ierr 	= getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbm->numSoft-1); CHKERRQ(ierr);
	fb->ID  = ID;

	// get pointer to specified softening law
	s = dbm->matSoft + ID;

	// check ID
	if(s->ID != -1)
	{
		 SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duplicate softening law!");
	}

	// set ID
	s->ID = ID;

	// read and store softening law parameters
	ierr = getScalarParam(fb, _REQUIRED_, "A",    &s->A,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "APS1", &s->APS1, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _REQUIRED_, "APS2", &s->APS2, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Lm",   &s->Lm,   1, 1.0); CHKERRQ(ierr);

	if(!s->A || !s->APS1 || !s->APS2)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "A, APS1, APS2 parameters must be nonzero for softening law %lld", (LLD)ID);
	}

	if (PrintOutput){
		if(s->Lm)
		{
			PetscPrintf(PETSC_COMM_WORLD,"   SoftLaw [%lld] : A = %g, APS1 = %g, APS2 = %g, Lm = %g\n", (LLD)(s->ID), s->A, s->APS1, s->APS2, s->Lm);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD,"   SoftLaw [%lld] : A = %g, APS1 = %g, APS2 = %g\n", (LLD)(s->ID), s->A, s->APS1, s->APS2);
		}
	}

	// SCALE

	s->Lm /= scal->length;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBMatReadPhase"
PetscErrorCode DBMatReadPhase(DBMat *dbm, FB *fb, PetscBool PrintOutput)
{
	// read material properties from file with error checking
	Scaling    *scal;
	Material_t *m;
	PetscInt    ID = -1, visID = -1, chSoftID, frSoftID, MSN, print_title;
	size_t 	    StringLength;
	PetscScalar eta, eta0, e0, Kb, G, E, nu, Vp, Vs, eta_st;
	PetscScalar healTau;   // NEW FOR HEALING
	char        ndiff[_str_len_], ndisl[_str_len_], npeir[_str_len_], title[_str_len_];
	char        PhaseDiagram[_str_len_], PhaseDiagram_Dir[_str_len_], Name[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	scal = dbm->scal;

	// initialize additional parameters
	eta      =  0.0;
	eta0     =  0.0;
	e0       =  0.0;
	//K        =  0.0;	// note: will be deprecated and renamed to Kb; we spit put an error message for now if we still find it in the input file
	Kb    	 =  0.0;	// bulk modulus		
	G        =  0.0;
	E        =  0.0;
	nu       =  0.0;
	Vp       =  0.0;
	Vs       =  0.0;
	eta_st   =  0.0;
	healTau = 1e30;   // NEW FOR HEALING, default value so we don't need an if-loop, healTau is always set
	chSoftID = -1;
	frSoftID = -1;
	MSN      =  dbm->numSoft - 1;

	// phase ID
	ierr 	 = getIntParam(fb, _REQUIRED_, "ID", &ID, 1, dbm->numPhases-1); CHKERRQ(ierr);
	fb->ID	 = ID;
	
	// get pointer to specified phase
	m = dbm->phases + ID;

	// check ID
	if(m->ID != -1)
	{
		 SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, "Duplicate phase definition!");
	}

	// set ID
	m->ID = ID;

	// visualization ID
	ierr = getIntParam(fb, _OPTIONAL_, "visID", &visID, 1, 0); CHKERRQ(ierr);

	if(visID != -1) m->visID = visID;
	else            m->visID = ID;

	// Name of the phase (mostly for visualization purposes)
	ierr = getStringParam(fb, _OPTIONAL_, "Name", Name, "none"); CHKERRQ(ierr);
	if(strcmp(Name, "none"))
	{
		strcpy(m->Name, Name);
	}

	//============================================================
	// density & phase diagram info
	//============================================================
	// Get the name of the phase diagram
	ierr = getStringParam(fb, _OPTIONAL_, "rho_ph",   PhaseDiagram, "none");          CHKERRQ(ierr);
	if (strcmp(PhaseDiagram, "none"))
	{
		// Note: the maximum length of the string PhaseDiagram is _str_len_
		// internally, however, a smaller string length is employed to save spac
		StringLength = strlen(PhaseDiagram)+3;		// 3, because we will add ".in" to the filename 

		// implies we are loading a phase diagram file from disk
		m->pdAct = 1;
		
		// Get the directory of the phase diagram if specified
		ierr = getStringParam(fb, _OPTIONAL_, "rho_ph_file", PhaseDiagram_Dir, "none"); CHKERRQ(ierr);
		if(strcmp(PhaseDiagram_Dir, "none"))
		{
			StringLength = StringLength + strlen(PhaseDiagram_Dir);	
			strcpy(m->pdf, PhaseDiagram_Dir);
		}
		// check that the length of the directory and the length of the file name does not exceed 	
		if (StringLength>_pd_name_sz_){
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"The length of the Phase Diagram Name and directory exceeds the maximum allowed length of %i /n", _pd_name_sz_);
		}

		// copy string
		strcat(m->pdf, PhaseDiagram);
		strcat(m->pdf, ".in");		// add the file ending

		strcpy(m->pdn, PhaseDiagram);

	    // Take into account only the melt, and not the density from a phase diagram
		ierr = getIntParam(fb, _OPTIONAL_, "Phase_Melt", &m->Phase_Diagram_melt, 1, 1); CHKERRQ(ierr);

	}
	else
	{
		m->pdAct = 0;	// no phase diagram is used
	}
	
	// Default Melt_Parametrization value


	//============================================================
	// Creep profiles
	//============================================================
	// set predefined diffusion creep profile
	ierr = GetProfileName(fb, scal, ndiff, "diff_prof"); CHKERRQ(ierr);
	ierr = SetDiffProfile(m, ndiff);                     CHKERRQ(ierr);
	// set predefined dislocation creep profile
	ierr = GetProfileName(fb, scal, ndisl, "disl_prof"); CHKERRQ(ierr);
	ierr = SetDislProfile(m, ndisl);                     CHKERRQ(ierr);
	// set predefined Peierls creep profile
	ierr = GetProfileName(fb, scal, npeir, "peir_prof"); CHKERRQ(ierr);
	ierr = SetPeirProfile(m, npeir);                     CHKERRQ(ierr);
	//=================================================================================
	// density
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "rho",      &m->rho,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_n",    &m->rho_n, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_c",    &m->rho_c, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "beta",     &m->beta,  1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// elasticity
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "G",        &G,        1, 1.0); CHKERRQ(ierr);
//	ierr = getScalarParam(fb, _OPTIONAL_, "K",        &K,        1, 1.0); CHKERRQ(ierr); // note-> will be removed (avoid confusion with k)
	ierr = getScalarParam(fb, _OPTIONAL_, "Kb",       &Kb,       1, 1.0); CHKERRQ(ierr); // note-> new nomenclature of bulk modulus (avoid confusion with k)
	ierr = getScalarParam(fb, _OPTIONAL_, "E",        &E,        1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "nu",       &nu,       1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Kp",       &m->Kp,    1, 1.0); CHKERRQ(ierr);
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
	ierr = getScalarParam(fb, _OPTIONAL_, "En",       &m->En,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Vn",       &m->Vn,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "n",        &m->n,     1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// Peierls creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "Bp",       &m->Bp,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Ep",       &m->Ep,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Vp",       &m->Vp,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "taup",     &m->taup,  1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "gamma",    &m->gamma, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "q",        &m->q,     1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// dc-creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "Bdc",      &m->Bdc,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Edc",      &m->Edc,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Rdc",      &m->Rdc,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "mu",       &m->mu,    1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// ps-creep
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "Bps",      &m->Bps,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Eps",      &m->Eps,   1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "d",        &m->d,     1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// plasticity (Drucker-Prager)
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "ch",       &m->ch,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "fr",       &m->fr,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "eta_st",   &eta_st,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rp",       &m->rp,     1, 1.0); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "chSoftID", &chSoftID,  1, MSN); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "frSoftID", &frSoftID,  1, MSN); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "healTau",  &m->healTau, 1, 1.0); CHKERRQ(ierr); // NEW FOR HEALING
	//=================================================================================
	// energy
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "alpha",    &m->alpha, 1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "Cp",       &m->Cp,    1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "k",        &m->k,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "A",        &m->A,     1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "T",        &m->T,     1, 1.0); CHKERRQ(ierr);
	//=================================================================================
	// melt fraction viscosity parametrization
	//=================================================================================
	ierr = getScalarParam(fb, _OPTIONAL_, "mfc",      &m->mfc,    1, 1.0);  CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "rho_melt", &m->rho_melt,1, 1.0);  CHKERRQ(ierr);


	// DEPTH-DEPENDENT

	// check depth-dependent density parameters
	if((!m->rho_n && m->rho_c) || (m->rho_n && !m->rho_c))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Depth-dependent density parameters must be specified simultaneously for phase %lld (rho_n + rho_c)", (LLD)ID);
	}

	if(m->rp < 0.0 || m->rp > 1.0)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Pore pressure ratio must be between 0 and 1 for phase %lld (rp)", (LLD)ID);
	}

	// PLASTICITY

	if(m->fr && !m->ch)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cohesion must be specified for phase %lld (fr + ch)", (LLD)ID);
	}

	if(!m->fr && frSoftID != -1)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Friction angle must be specified for phase %lld (frSoftID + fr)", (LLD)ID);
	}

	if(!m->ch && chSoftID != -1)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cohesion must be specified for phase %lld (chSoftID + ch)", (LLD)ID);
	}
	
	if((!m->rho_melt && m->Phase_Diagram_melt))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "You need to specify the density of the melting phase for phase %lld", (LLD)ID);
	}




    m->eta_st   = eta_st;
   
	// set softening law IDs
	m->chSoftID = chSoftID;
	m->frSoftID = frSoftID;

	// DIFFUSION

	if(!(( eta && !m->Bd)   // eta
	||   (!eta &&  m->Bd)   // Bd
	||   (!eta && !m->Bd))) // nothing
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Diffusion creep parameters are not unique for phase %lld (eta, Bd)\n", (LLD)ID);
	}

	// compute diffusion creep constant
	if(eta) m->Bd = 1.0/(2.0*eta);

	// DISLOCATION

	if(!(( eta0 &&  e0 &&  m->n && !m->Bn)   // eta0, e0, n
	||   (!eta0 && !e0 &&  m->n &&  m->Bn)   // Bn, n
	||   (!eta0 && !e0 && !m->n && !m->Bn))) // nothing
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Dislocation creep parameters are not unique for phase %lld (eta0 + e0 + n, Bn + n)\n", (LLD)ID);
	}

	// compute dislocation creep constant
	if(eta0 && e0 && m->n) m->Bn = pow(2.0*eta0, -m->n)*pow(e0, 1 - m->n);

	// PEIERLS

	if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Peierls creep parameters are incomplete for phase %lld (Bp + taup + gamma + q + Ep)", (LLD)ID);
	}

	if(m->Bp && !m->Bn)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Peierls creep requires dislocation creep for phase %lld (Bp, Bn)", (LLD)ID);
	}

	// DC

	if(m->Bdc && (!m->Edc || !m->Rdc || !m->mu))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "dc-creep parameters are incomplete for phase %lld (Bdc + Edc + Rdc + mu)", (LLD)ID);
	}

	if(m->Bdc && m->Bn)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot combine dc-creep with dislocation creep for phase %lld (Bdc + Bn)", (LLD)ID);
	}

	// PS

	if(m->Bps && (!m->Eps || !m->d))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "ps-creep parameters are incomplete for phase %lld (Bps + Eps + d)", (LLD)ID);
	}

	if(m->Bps && !m->Bdc && !m->Bn)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "ps-creep requires either dc-creep or dislocation creep for phase %lld (Bps + Bdc, Bps + Bn)", (LLD)ID);
	}

	if(m->Bps && m->Bd)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot combine ps-creep with diffusion creep for phase %lld (Bps + Bd)", (LLD)ID);
	}

	// ELASTICITY
	// I'm taking out this warning message as it interfers with Adjoint and defining k (conductivity) from the command-line
	//if (K){
	//	SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "The bulk modulus parameter is now called 'Kb' and no longer 'K'; change your ParamFile accordingly");
	//}

	if(!(( G && !Kb && !E && !nu)   // G
	||   (!G &&  Kb && !E && !nu)   // Kb
	||   ( G &&  Kb && !E && !nu)   // G & Kb
	||   ( G && !Kb && !E &&  nu)   // G & nu
	||   (!G &&  Kb && !E &&  nu)   // Kb & nu
	||   (!G && !Kb &&  E &&  nu)   // E & nu
	||   (!G && !Kb && !E && !nu))) // nothing
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unsupported or nonunique combination of elasticity parameters for phase %lld (G, Kb, G + Kb, G + nu, Kb + nu, E + nu)\n", (LLD)ID);
	}

	if(m->Kp && !Kb)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Bulk modulus (Kb) must be specified for phase %lld (Kb + Kp)", (LLD)ID);
	}

	if(m->beta && Kb)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Density pressure dependence parameters are not unique for phase %lld (beta, Kb)", (LLD)ID);
	}

	// compute elastic parameters
	if( G  && nu)          Kb  = 2*G*(1 + nu)/(3*(1 - 2*nu));
	if( Kb && nu)          G   = (3*Kb*(1 - 2*nu))/(2*(1 + nu));
	if( E  && nu)        { Kb  = E/(3*(1 - 2*nu)); G = E/(2*(1 + nu)); }
	if(!E  && Kb && G)      E  = 9*Kb*G/(3*Kb + G);
	if(!nu && Kb && G)      nu = (3*Kb - 2*G)/(2*(3*Kb + G));
	if( Kb  && G && m->rho) Vp = sqrt((Kb + 4.0*G/3.0)/m->rho);
	if( G  && m->rho)       Vs = sqrt((G/m->rho));

	// store elastic moduli
	m->G  = G;
	m->Kb = Kb;


	// check that at least one essential deformation mechanism is specified
	if(!m->Bd && !m->Bn && !m->G && !m->Bdc)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), Bdc, G must be specified for phase %lld", (LLD)ID);
	}

	// PRINT (optional)
	if (PrintOutput){
		if (strlen(m->Name)>0){
			PetscPrintf(PETSC_COMM_WORLD,"   Phase ID : %lld     --   %s ",(LLD)(m->ID), m->Name);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD,"   Phase ID : %lld  %s ",(LLD)(m->ID), m->Name);
			
		}

		if(strlen(ndiff)) PetscPrintf(PETSC_COMM_WORLD,"\n   diffusion creep profile  : %s", ndiff);
		if(strlen(ndisl)) PetscPrintf(PETSC_COMM_WORLD,"\n   dislocation creep profile: %s", ndisl);

		sprintf(title, "   (dens)   : "); print_title = 1;
		MatPrintScalParam(m->rho,   "rho",   "[kg/m^3]", scal, title, &print_title);
		if(m->pdAct == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD,"- Employing phase diagram: %s", PhaseDiagram);
		}
		else
		MatPrintScalParam(m->rho_n, "rho_n", "[ ]",      scal, title, &print_title);
		MatPrintScalParam(m->rho_c, "rho_c", "[1/m]",    scal, title, &print_title);
		MatPrintScalParam(m->beta,  "beta",  "[1/Pa]",   scal, title, &print_title);
		MatPrintScalParam(m->rho_melt, "rho",     "[kg/m^3]",      scal, title, &print_title);


		sprintf(title, "   (elast)  : "); print_title = 1;
		MatPrintScalParam(G,     "G",  "[Pa]",  scal, title, &print_title);
		MatPrintScalParam(Kb,    "Kb",  "[Pa]",  scal, title, &print_title);
		MatPrintScalParam(E,     "E",  "[Pa]",  scal, title, &print_title);
		MatPrintScalParam(nu,    "nu", "[ ]",   scal, title, &print_title);
		MatPrintScalParam(m->Kp, "Kp", "[ ]",   scal, title, &print_title);
		MatPrintScalParam(Vp,    "Vp", "[m/s]", scal, title, &print_title);
		MatPrintScalParam(Vs,    "Vs", "[m/s]", scal, title, &print_title);

		sprintf(title, "   (diff)   : "); print_title = 1;
		MatPrintScalParam(eta,   "eta", "[Pa*s]",    scal, title, &print_title);
		MatPrintScalParam(m->Bd, "Bd",  "[1/Pa/s]",  scal, title, &print_title);
		MatPrintScalParam(m->Ed, "Ed",  "[J/mol]",   scal, title, &print_title);
		MatPrintScalParam(m->Vd, "Vd",  "[m^3/mol]", scal, title, &print_title);

		sprintf(title, "   (disl)   : "); print_title = 1;
		MatPrintScalParam(eta0,  "eta0", "[Pa*s]",     scal, title, &print_title);
		MatPrintScalParam(e0,    "e0",   "[1/s]",      scal, title, &print_title);
		MatPrintScalParam(m->Bn, "Bn",   "[1/Pa^n/s]", scal, title, &print_title);
		MatPrintScalParam(m->En, "En",   "[J/mol]",    scal, title, &print_title);
		MatPrintScalParam(m->Vn, "Vn",   "[m^3/mol]",  scal, title, &print_title);
		MatPrintScalParam(m->n,  "n",    "[ ]",        scal, title, &print_title);

		sprintf(title, "   (peirl)  : "); print_title = 1;
		MatPrintScalParam(m->Bp,    "Bp",    "[1/s]",     scal, title, &print_title);
		MatPrintScalParam(m->Ep,    "Ep",    "[J/mol]",   scal, title, &print_title);
		MatPrintScalParam(m->Vp,    "Vp",    "[m^3/mol]", scal, title, &print_title);
		MatPrintScalParam(m->taup,  "taup",  "[Pa]",      scal, title, &print_title);
		MatPrintScalParam(m->gamma, "gamma", "[ ]",       scal, title, &print_title);
		MatPrintScalParam(m->q,     "q",     "[ ]",       scal, title, &print_title);

		sprintf(title, "   (dc)     : "); print_title = 1;
		MatPrintScalParam(m->Bdc,   "Bdc",  "[1/s]",   scal, title, &print_title);
		MatPrintScalParam(m->Edc,   "Edc",  "[J/mol]", scal, title, &print_title);
		MatPrintScalParam(m->Rdc,   "Rdc",  "[ ]",     scal, title, &print_title);
		MatPrintScalParam(m->mu,    "mu",   "[Pa]",    scal, title, &print_title);

		sprintf(title, "   (ps)     : "); print_title = 1;
		MatPrintScalParam(m->Bps,   "Bps",  "[K*m^3/Pa/s]", scal, title, &print_title);
		MatPrintScalParam(m->Eps,   "Eps",  "[J/mol]",      scal, title, &print_title);
		MatPrintScalParam(m->d,     "d",    "[m]",          scal, title, &print_title);

		sprintf(title, "   (plast)  : "); print_title = 1;
		MatPrintScalParam(m->ch,     "ch",     "[Pa]",   scal, title, &print_title);
		MatPrintScalParam(m->fr,     "fr",     "[deg]",  scal, title, &print_title);
		MatPrintScalParam(m->eta_st, "eta_st", "[Pa*s]", scal, title, &print_title);
		MatPrintScalParam(m->rp,     "rp",     "[ ]",    scal, title, &print_title);
		MatPrintScalParam(m->healTau,"healTau","[Myr]",    scal, title, &print_title);   // NEW FOR HEALING
		if(frSoftID != -1) PetscPrintf(PETSC_COMM_WORLD, "frSoftID = %lld ", (LLD)frSoftID);
		if(chSoftID != -1) PetscPrintf(PETSC_COMM_WORLD, "chSoftID = %lld ", (LLD)chSoftID);

		sprintf(title, "   (temp)   : "); print_title = 1;
		MatPrintScalParam(m->alpha, "alpha", "[1/K]",    scal, title, &print_title);
		MatPrintScalParam(m->Cp,    "Cp",    "[J/kg/K]", scal, title, &print_title);
		MatPrintScalParam(m->k,     "k",     "[W/m/k]",  scal, title, &print_title);
		MatPrintScalParam(m->A,     "A",     "[W/kg]",   scal, title, &print_title);
		MatPrintScalParam(m->T,     "T",     "[C]",      scal, title, &print_title);
		PetscPrintf(PETSC_COMM_WORLD,"\n\n");
	}

	// SCALE

	// NOTE: [1] activation energy is not scaled, gas constant is scaled with temperature (RT)
	//       [2] activation volume is scaled with characteristic stress in SI units       (pV)

	m->rho    /= scal->density;
	m->rho_c  *= scal->length_si;
	m->beta   *= scal->stress_si; // [1/Pa]
	m->rho_melt /= scal->density;

	// diffusion creep
	m->Bd     *= scal->viscosity;
	m->Vd     *= scal->stress_si;

	// dislocation creep (power-law)
	m->Bn     *= pow(scal->stress_si, m->n)*scal->time_si;
	m->Vn     *= scal->stress_si;

	// Peierls creep
	m->Bp     /=  scal->strain_rate;
	m->Vp     *=  scal->stress_si;
	m->taup   /=  scal->stress_si;

	// dc-creep
	m->Bdc    /=  scal->strain_rate;
	m->mu     /=  scal->stress_si;

	// ps-creep
	m->Bps   *= scal->viscosity/scal->volume_si/scal->temperature;
	m->d     /= scal->length_si;

	// elasticity
	m->G      /= scal->stress_si;
	m->Kb     /= scal->stress_si;

	// plasticity
	m->ch     /= scal->stress_si;
	m->fr     /= scal->angle;
	m->healTau /= scal->time;            // NEW FOR HEALING

    m->eta_st /= scal->viscosity;
    
	// temperature
	m->alpha  /= scal->expansivity;
	m->Cp     /= scal->cpecific_heat;
	m->k      /= scal->conductivity;
	m->A      /= scal->heat_production;

	// phase-temperature
	if(m->T) m->T = (m->T + scal->Tshift)/scal->temperature;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void MatPrintScalParam(
		PetscScalar par,  const char key[],   const char label[],
		Scaling    *scal, const char title[], PetscInt   *print_title)
{
	// monitor parameter value

	if(par == 0.0) return;

	if((*print_title))
	{
		PetscPrintf(PETSC_COMM_WORLD, "\n%s", title);

		(*print_title) = 0;
	}

	if(scal->utype == _NONE_)
	{
		PetscPrintf(PETSC_COMM_WORLD, "%s = %g [ ]  ", key, par);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "%s = %g %s  ", key, par, label);
	}
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

	ierr = getStringParam(fb, _OPTIONAL_, key, name, NULL);  CHKERRQ(ierr);

	if(strlen(name) && scal->utype == _NONE_)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Predefined creep profile is not supported for non-dimensional setup");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetDiffProfile"
PetscErrorCode SetDiffProfile(Material_t *m, char name[])
{
	// set diffusion creep profiles from literature

	// We assume that the creep law has the form:
	// Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
	// Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
	//
	// In addition, we take into account that the creep-laws are typically measured under uniaxial or simple shear,
	// whereas we need them in tensorial format (for correction see e.g. T. Gerya book).
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

	ExpType          type;
	PetscInt         MPa;
	PetscScalar      d0, p;
	PetscScalar      C_OH_0, r;

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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
		MPa              =   1;
		d0               =   10e3;
		p                =   3;
		C_OH_0           =   1000;
		r                =   0.8;
	}
	else if (!strcmp(name,"Dry_Plagioclase_RybackiDresen_2000"))
	{
		// after Rybacki and Dresen, 2000, JGR 
		m->Bd            =   1.2589e12;		
		m->Ed            =   460e3;
		m->Vd            =   24e-6;
		type             =   _UniAxial_;
		MPa              =   1;
		d0               =   100;					// in microns in their paper
		p                =   3;
		C_OH_0           =   1;
		r                =   0.0;
	}
	else if (!strcmp(name,"Wet_Plagioclase_RybackiDresen_2000"))
	{
		// after Rybacki and Dresen, 2000, JGR 
		m->Bd            =   0.1995;		
		m->Ed            =   159e3;
		m->Vd            =   38e-6;
		type             =   _UniAxial_;
		MPa              =   1;
		d0               =   100;					// in microns in their paper
		p                =   3;
		C_OH_0           =   158.4893;				// fugacity water 10^2.2, not sure whether that is implemented correctly
		r                =   1.0;					// 1.0 in Rybacki paper, but to be checked that units match
	}
	else
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No such diffusion creep profile: %s! ",name);
	}

	// correct experimental creep prefactor to tensor units
	ierr = CorrExpPreFactor(m->Bd, 1, type, MPa); CHKERRQ(ierr);

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

	// We assume that the creep law has the form:
	// Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
	// Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
	//
	// In addition, we take into account that the creep-laws are typically measured under uniaxial or simple shear,
	// whereas we need them in tensorial format (for correction see e.g. T. Gerya book).
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

	ExpType          type;
	PetscInt         MPa;
	PetscScalar      C_OH_0, r;

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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _SimpleShear_;
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
		type             =   _UniAxial_;
		MPa              =   1;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
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
		type             =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Dry_Plagioclase_RybackiDresen_2000"))
	{
		m->Bn            =   5.0119e12;		
		m->n             =   3.0;
		m->En            =   641e3;
		m->Vn            =   24e-6;
		type             =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   1;
		r                =   0;
	}

	else if (!strcmp(name,"Wet_Plagioclase_RybackiDresen_2000"))
	{
		m->Bn            =   1.5849;		
		m->n             =   3.0;
		m->En            =   345e3;
		m->Vn            =   38e-6;
		type             =   _UniAxial_;
		MPa              =   1;
		C_OH_0           =   158.4893;				// fugacity water 10^2.2, not sure whether that is implemented correctly
		r                =   1.0;					// 1.0 in Rybacki paper, but to be checked that units match

	}

	else if (!strcmp(name,"Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003"))
	{
		// after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
		// Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
		m->Bn            =   1600;
		m->n             =   3.5;
		m->En            =   520e3;
		m->Vn            =   22e-6;
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
		type             =   _SimpleShear_;
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
        type             =   _UniAxial_;
        MPa              =   1;
        C_OH_0           =   1;
        r                =   0;
    }

	else if (!strcmp(name,"RockSaltReference_BGRa_class3-Braeumer_et_al_2011"))
    {
        // Taken from Bräuer et al. (2011) Description of the Gorleben site (PART 4): Geotechnical exploration of the Gorleben salt dome - page 126
        m->Bn            =   5.2083e-7; // 2^(class==3)/32 * 0.18 * 1/24/3600 * (1MPa^(-1))^5
        m->n             =   5;
        m->En            =   54e3;
        m->Vn            =   0;
        type             =   _UniAxial_;
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
        type             =   _UniAxial_;
        MPa              =   1;
        C_OH_0           =   1;
        r                =   0;
    }
    
	else
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No such dislocation creep profile: %s! ",name);
	}

	// correct experimental creep prefactor to tensor units
	ierr = CorrExpPreFactor(m->Bn, m->n, type, MPa); CHKERRQ(ierr);

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
	// taup       - Peierls stress [Pa]
	// gamma      - adjustable constant [-]
	// q          - stress dependence for Peierls creep [-]
	// s          - Peierls creep exponent (typical values between 7-11) [-]

	ExpType  type;
	PetscInt MPa;

	PetscErrorCode ierr;
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
        type             = _None_; // what to use for indenter-type tests? (set to none for now)
        MPa              = 0;
	}

	else
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No such Peierls creep profile: %s! ",name);
	}

	// correct Peierls prefactor & stress to tensor units
	ierr = CorrExpStressStrainRate(m->Bp, m->taup, type, MPa); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CorrExpPreFactor"
PetscErrorCode CorrExpPreFactor(PetscScalar &B, PetscScalar n, ExpType type, PetscInt MPa)
{
	// correct experimental creep prefactor to tensor units
	// correction factor depends on experiment type (see e.g. Gerya, 2010, chapter 6, p. 71-78)

	// B - creep prefactor
	// n - power law exponent

	PetscFunctionBegin;

	// apply experimental to tensor correction
	if      (type == _UniAxial_)    B *= pow(3.0, (n + 1.0)/2.0)/2.0; //  F = (3^(n+1)/2)/2
	else if (type == _SimpleShear_) B *= pow(2.0,  n - 1.0);          //  F =  2^(n-1)
	else if (type != _None_)
	{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown rheology experiment type!");
	}

	// apply correction from [MPa^(-n) s^(-1)] to [Pa^(-n) s^(-1)] if required
	if(MPa) B *= pow(10, -6.0*n);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CorrExpStressStrainRate"
PetscErrorCode CorrExpStressStrainRate(PetscScalar &D, PetscScalar &S, ExpType type, PetscInt MPa)
{
	// correct experimental stress and strain rate parameters to tensor units
	// correction factor depends on experiment type (see e.g. Gerya, 2010, chapter 6, p. 71-78)

	// D - creep strain rate parameter (e.g. Peierls prefactor)
	// S - creep stress parameter (e.g. Peierls stress)

	PetscFunctionBegin;

	// apply experimental to tensor correction
	if      (type == _UniAxial_)    { D *= sqrt(3.0)/2.0; S /= sqrt(3.0); }
	else if (type == _SimpleShear_) { D /= 2.0;           S /= 2.0;       }
	else if (type != _None_)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown rheology experiment type!");
	}

	// apply correction from [MPa] to [Pa] if required
	if(MPa) S *= 1e6;

	PetscFunctionReturn(0);
}
//------------------------------------------------------------------

//--------------------------------------------------------------------------

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
		PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD," Material properties set from calling function: \n");


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
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)id);
				}

				// check power-law exponent
				if(!m->n && ((eta0 && e0) || m->Bn))
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)id);
				}

				// check Peierls creep
				if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)id);
				}

				// recompute creep parameters
				if(eta)        m->Bd = 1.0/(2.0*eta);
				if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

				// check that at least one essential deformation mechanism is specified
				if(!m->Bd && !m->Bn && !m->G)
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)id);
				}

				PetscPrintf(PETSC_COMM_WORLD,"    eta[%lld] = %g \n",(LLD)id,eta);
			}

			// constant density
			else if(mod->typ[im] == _RHO0_)
			{
				m->rho = mod->val[im];
				PetscPrintf(PETSC_COMM_WORLD,"    rho0[%lld] = %5.5f \n",(LLD)id,m->rho);
			}

			else
			{
				PetscPrintf(PETSC_COMM_WORLD,"WARNING: inversion parameter type is not implemented \n");
			}

		}
		PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");
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
		PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD," Material properties set from command line: \n");

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
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)id);
				}

				// check power-law exponent
				if(!m->n && ((eta0 && e0) || m->Bn))
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)id);
				}

				// check Peierls creep
				if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep))
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)id);
				}

				// recompute creep parameters
				if(eta)        m->Bd = 1.0/(2.0*eta);
				if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

				// check that at least one essential deformation mechanism is specified
				if(!m->Bd && !m->Bn && !m->G)
				{
					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)id);
				}

			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"    eta[%lld]	= %g \n",(LLD)id,eta);

			// constant density
			sprintf(matprop_opt,"-rho0_%lld",(LLD)id);
			ierr = PetscOptionsGetReal(NULL, NULL ,matprop_opt,&m->rho	, &flg);			CHKERRQ(ierr);
			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"    rho0[%lld]	= %5.5f \n",(LLD)id,m->rho);

		}
		PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PrintMatProp"
PetscErrorCode PrintMatProp(Material_t *MatProp)
{
	// Prints an overview of the material properties specified for a certain phase (for debugging)
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,">>> Material properties for phase %i with visId=%i : \n",MatProp->ID, MatProp->visID);
	PetscPrintf(PETSC_COMM_WORLD,">>> Density:          rho   = %1.7e,  rho_n = %1.7e,    rho_c = %1.7e,   beta = %1.7e \n",  MatProp->rho,MatProp->rho_c, MatProp->rho_c, MatProp->beta);
	PetscPrintf(PETSC_COMM_WORLD,">>> Elasticity:       Kb    = %1.7e,  Kp    = %1.7e,    G     = %1.7e \n",                MatProp->Kb, MatProp->Kp, MatProp->G);
	PetscPrintf(PETSC_COMM_WORLD,">>> Diffusion Cr.:    Bd    = %1.7e,  Ed    = %1.7e,    Vd    = %1.7e \n",                MatProp->Bd, MatProp->Ed, MatProp->Vd);
	PetscPrintf(PETSC_COMM_WORLD,">>> Dislocation Cr.:  Bn    = %1.7e,  n     = %1.7e,    Vn    = %1.7e,    En  = %1.7e \n", MatProp->Bn, MatProp->n, MatProp->Vn, MatProp->En);
	PetscPrintf(PETSC_COMM_WORLD,">>> Peierls Cr.:      Bp    = %1.7e,  Ep    = %1.7e,    Vp    = %1.7e,    taup= %1.7e,    gamma  = %1.7e,     q        = %1.7e \n", MatProp->Bp, MatProp->Ep, MatProp->Vp, MatProp->taup, MatProp->gamma, MatProp->q);
	PetscPrintf(PETSC_COMM_WORLD,">>> dc Cr.:           Bdc   = %1.7e,  Edc   = %1.7e,    Rdc   = %1.7e,    mu  = %1.7e \n", MatProp->Bdc, MatProp->Edc, MatProp->Rdc, MatProp->mu);
    PetscPrintf(PETSC_COMM_WORLD,">>> ps Cr.:           Bps   = %1.7e,  Eps   = %1.7e,    d     = %1.7e,    \n", MatProp->Bps, MatProp->Eps, MatProp->d);
    
    PetscPrintf(PETSC_COMM_WORLD,">>> Plasticity:       fr    = %1.7e,  ch    = %1.7e,    eta_st= %1.7e,    rp= %1.7e,    frSoftID = %i,                chSoftID = %i \n", MatProp->fr, MatProp->ch, MatProp->eta_st, MatProp->rp, MatProp->frSoftID, MatProp->chSoftID);
	PetscPrintf(PETSC_COMM_WORLD,">>> Thermal:          alpha = %1.7e,  Cp    = %1.7e,    k     = %1.7e,    A = %1.7e,    T        = %1.7e \n", MatProp->alpha, MatProp->Cp, MatProp->k, MatProp->A, MatProp->T);

	PetscPrintf(PETSC_COMM_WORLD," \n");

  
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DBMatOverwriteWithGlobalVariables"
PetscErrorCode DBMatOverwriteWithGlobalVariables(DBMat *dbm, FB *fb)
{
    PetscFunctionBegin;

    PetscErrorCode  ierr;
    PetscScalar     eta_min;
    PetscInt        ID;
    Material_t      *m;
    Scaling         *scal;
	PetscFunctionBegin;

	// access context
	scal    = dbm->scal;

    eta_min = 0;
    ierr    = getScalarParam(fb, _OPTIONAL_, "eta_min",         &eta_min,        1, 1.0); CHKERRQ(ierr);

	for(ID = 0; ID < dbm->numPhases; ID++)
	{   
    	// get pointer to specified phase
        m = &dbm->phases[ID];
      
        // Plasticity stabilization viscosity: set to eta_min, if not defined 
        if(!m->eta_st) m->eta_st = eta_min/scal->viscosity; // set default stabilization viscosity if not defined    
        
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
