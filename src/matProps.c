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
#include "solVar.h"
#include "fdstag.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "parsing.h"
#include "matProps.h"
#include "tools.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropInit"
PetscErrorCode MatPropInit(JacRes *jr, FILE *fp)
{
	// initialize MATERIAL PARAMETERS from file

	PetscInt  *ls, *le;
	PetscInt   i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Reading material parameters: \n\n");

	// clear memory
	ierr = PetscMemzero(jr->phases, sizeof(Material_t)*(size_t)max_num_phases); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(i = 0; i < max_num_phases; i++)
	{
		jr->phases[i].ID = -1;
	}

	// allocate memory for arrays to store line info
	ierr = makeIntArray(&ls, NULL, max_num_phases); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, max_num_phases); CHKERRQ(ierr);

	// read number of entries
	getLineStruct(fp, ls, le, max_num_phases, &count_starts, &count_ends, "<MaterialStart>","<MaterialEnd>");

	// error checking
	if(count_starts > max_num_phases || count_ends > max_num_phases)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)max_num_phases);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incomplete material structures! <MaterialStart> & <MaterialEnd> don't match");
	}

	// store actual number of materials
	jr->numPhases = count_starts;

	// read each individual phase
	for(i = 0; i < jr->numPhases; i++)
	{
		ierr = MatPropGetStruct(fp,
				jr->numPhases, jr->phases,
				jr->numSoft, jr->matSoft,
				ls[i], le[i], jr->scal.utype); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropGetStruct"
PetscErrorCode MatPropGetStruct(FILE *fp,
		PetscInt numPhases, Material_t *phases,
		PetscInt numSoft,   Soft_t     *matSoft,
		PetscInt ils, PetscInt ile, UnitsType utype)
{
	// read material properties from file with error checking
	// WARNING! This function assumes correctly defined softening parameters

	Material_t *m;
	PetscScalar eta, eta0, e0;
	PetscInt    ID = -1, chSoftID, frSoftID, found;
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

	// phase ID
	getMatPropInt(fp, ils, ile, "ID", &ID, &found);

	// error checking
	if(!found)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "No phase ID specified! ");
	}
	if(ID > numPhases - 1)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Incorrect phase numbering!");
	}

	// initialize additional parameters
	eta      =  0.0;
	eta0     =  0.0;
	e0       =  0.0;
	chSoftID = -1;
	frSoftID = -1;

	// get pointer to specified phase
	m = phases + ID;

	// check ID
	if(m->ID != -1)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Incorrect phase numbering!");
	}

	// set ID
	m->ID = ID;

	//============================================================
	// density
	//============================================================
	getMatPropScalar(fp, ils, ile, "rho0",      &m->rho,   NULL);
	getMatPropScalar(fp, ils, ile, "rho_n",     &m->rho_n, NULL);
	getMatPropScalar(fp, ils, ile, "rho_c",     &m->rho_c, NULL);
	getMatPropScalar(fp, ils, ile, "beta",      &m->beta,  NULL); // pressure-dependence of density

	//============================================================
	// Creep profiles
	//============================================================
	// set predefined diffusion creep profile
	getMatPropString(fp, ils, ile, "diff_profile", ndiff, MAX_NAME_LEN, &found);
	if(found) { ierr = SetDiffProfile(m, ndiff); CHKERRQ(ierr); }

	// set predefined dislocation creep profile
	getMatPropString(fp, ils, ile, "disl_profile", ndisl, MAX_NAME_LEN, &found);
	if(found) { ierr = SetDislProfile(m, ndisl); CHKERRQ(ierr); }

	// set predefined Peierls creep profile
	getMatPropString(fp, ils, ile, "peir_profile", npeir, MAX_NAME_LEN, &found);
	if(found) { ierr = SetPeirProfile(m, npeir); CHKERRQ(ierr); }

	//============================================================
	// Newtonian linear diffusion creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "eta",       &eta,      NULL);
	getMatPropScalar(fp, ils, ile, "Bd",        &m->Bd,    NULL);
	getMatPropScalar(fp, ils, ile, "Ed",        &m->Ed,    NULL);
	getMatPropScalar(fp, ils, ile, "Vd",        &m->Vd,    NULL);

	//============================================================
	// power-law (dislocation) creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "eta0",      &eta0,     NULL);
	getMatPropScalar(fp, ils, ile, "e0",        &e0,       NULL);
	getMatPropScalar(fp, ils, ile, "Bn",        &m->Bn,    NULL);
	getMatPropScalar(fp, ils, ile, "n",         &m->n,     NULL);
	getMatPropScalar(fp, ils, ile, "En",        &m->En,    NULL);
	getMatPropScalar(fp, ils, ile, "Vn",        &m->Vn,    NULL);

	//============================================================
	// Peierls creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "Bp",        &m->Bp,    NULL);
	getMatPropScalar(fp, ils, ile, "taup",      &m->taup,  NULL);
	getMatPropScalar(fp, ils, ile, "gamma",     &m->gamma, NULL);
	getMatPropScalar(fp, ils, ile, "q",         &m->q,     NULL);
	getMatPropScalar(fp, ils, ile, "Ep",        &m->Ep,    NULL);
	getMatPropScalar(fp, ils, ile, "Vp",        &m->Vp,    NULL);
	//============================================================
	// elasticity
	//============================================================
	getMatPropScalar(fp, ils, ile, "shear",     &m->G,     NULL);
	getMatPropScalar(fp, ils, ile, "bulk",      &m->K,     NULL);
	getMatPropScalar(fp, ils, ile, "Kp",        &m->Kp,    NULL);
	//============================================================
	// plasticity (Drucker-Prager)
	//============================================================
	getMatPropScalar(fp, ils, ile, "cohesion",  &m->ch,    NULL);
	getMatPropScalar(fp, ils, ile, "friction",  &m->fr,    NULL);
	getMatPropScalar(fp, ils, ile, "lambda",    &m->rp,    NULL);
	getMatPropInt   (fp, ils, ile, "chSoftID",  &chSoftID, NULL);
	getMatPropInt   (fp, ils, ile, "frSoftID",  &frSoftID, NULL);
	//============================================================
	// energy
	//============================================================
	getMatPropScalar(fp, ils, ile, "alpha",     &m->alpha, NULL);
	getMatPropScalar(fp, ils, ile, "cp",        &m->Cp,    NULL);
	getMatPropScalar(fp, ils, ile, "k",         &m->k,     NULL);
	getMatPropScalar(fp, ils, ile, "A",         &m->A,     NULL);
	//============================================================

	// check depth-dependent density parameters
	if((!m->rho_n && m->rho_c) || (m->rho_n && !m->rho_c))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "rho_n & rho_c must be specified simultaneously for phase %lld", (LLD)ID);
	}

	// check softening laws
	if(chSoftID > numSoft-1)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect cohesion softening law specified for phase %lld", (LLD)ID);
	}
	if(frSoftID > numSoft-1)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect friction softening law specified for phase %lld", (LLD)ID);
	}

	if(m->fr && !m->ch)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Nonzero cohesion must be specified for phase %lld", (LLD)ID);
	}

	if((m->rp>1) || (m->rp<0))
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "pore pressure ratio must be between 0 and 1 for phase %lld", (LLD)ID);
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

	// check units for predefined profile
	if((strlen(ndiff) || strlen(ndisl)) && utype == _NONE_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot have a predefined creep profile (phase %lld), in a non-dimensional setup!", (LLD)ID);
	}

	// print
	if(utype == _NONE_)
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
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g %s, friction angle = %g %s, pore pressure ratio = %g [ ] \n", (LLD)(m->ID),m->ch, lbl_tau, m->fr, lbl_fr, m->rp);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [ ], friction SoftLaw = %lld [ ] \n", (LLD)(m->ID),(LLD)chSoftID, (LLD)frSoftID);
	PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g %s, cp = %g %s, k = %g %s, A = %g %s \n", (LLD)(m->ID),m->alpha, lbl_alpha, m->Cp, lbl_cp,m->k, lbl_k, m->A, lbl_A);
	PetscPrintf(PETSC_COMM_WORLD,"    \n");

	PetscFunctionReturn(0);
}

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
#undef __FUNCT__
#define __FUNCT__ "MatSoftInit"
PetscErrorCode MatSoftInit(JacRes *jr, FILE *fp)
{
	// initialize SOFTENING LAWS from file

	PetscInt    *ls,*le;
	PetscInt     i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of softening laws from file
	PetscPrintf(PETSC_COMM_WORLD,"Reading softening laws: \n\n");

	// clear memory
	ierr = PetscMemzero(jr->matSoft, sizeof(Soft_t)*(size_t)max_num_soft); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(i = 0; i < max_num_soft; i++)
	{
		jr->matSoft[i].ID = -1;
	}

	// allocate memory for arrays to store line info
	ierr = makeIntArray(&ls, NULL, max_num_soft); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, max_num_soft); CHKERRQ(ierr);

	// read number of entries
	getLineStruct(fp, ls, le, max_num_soft, &count_starts, &count_ends, "<SofteningStart>", "<SofteningEnd>");

	// error checking
	if(count_starts > max_num_soft || count_ends > max_num_soft)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many softening laws specified! Max allowed: %lld", (LLD)max_num_soft);
	}
	if(count_starts != count_ends)
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incomplete material structures! <SofteningStart> & <SofteningEnd> don't match");
	}

	// store actual number of softening laws
	jr->numSoft = count_starts;

	// read each individual softening law
	for(i = 0; i < jr->numSoft; i++)
	{
		ierr = MatSoftGetStruct(fp, jr->numSoft, jr->matSoft, ls[i], le[i]); CHKERRQ(ierr);
	}

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatSoftGetStruct"
PetscErrorCode MatSoftGetStruct(FILE *fp,
		PetscInt numSoft, Soft_t *matSoft,
		PetscInt ils, PetscInt ile)
{
	// read softening law from file

	Soft_t   *s;
	PetscInt  found, ID;

	PetscFunctionBegin;

	// softening law ID
	getMatPropInt(fp, ils, ile, "softID", &ID, &found);

	// error checking
	if(!found)
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "No softening law ID specified! ");
	}
	if(ID > numSoft - 1)
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect softening law numbering!");
	}

	// get pointer to specified softening law
	s = matSoft + ID;

	// check ID
	if(s->ID != -1)
	{
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect softening law numbering!");
	}

	// set ID
	s->ID = ID;

	// read and store softening law parameters
	getMatPropScalar(fp, ils, ile, "A",    &s->A,    &found);
	getMatPropScalar(fp, ils, ile, "APS1", &s->APS1, &found);
	getMatPropScalar(fp, ils, ile, "APS2", &s->APS2, &found);

	if(!s->A || !s->APS1 || !s->APS2)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "All parameters must be specified simultaneously for softening law %lld", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"SoftLaw [%lld]: A = %g, APS1 = %g, APS2 = %g \n", (LLD)(s->ID), s->A, s->APS1, s->APS2);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// set diffusion creep profiles from literature
#undef __FUNCT__
#define __FUNCT__ "SetDiffProfile"
PetscErrorCode SetDiffProfile(Material_t *m, char name[])
{
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
// set dislocation creep profiles from literature
#undef __FUNCT__
#define __FUNCT__ "SetDislProfile"
PetscErrorCode SetDislProfile(Material_t *m, char name[])
{
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
// set Peierls creep profiles from literature
#undef __FUNCT__
#define __FUNCT__ "SetPeirProfile"
PetscErrorCode SetPeirProfile(Material_t *m, char name[])
{
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
// set tensor and units correction for rheological profiles
#undef __FUNCT__
#define __FUNCT__ "SetProfileCorrection"
PetscErrorCode SetProfileCorrection(PetscScalar *B, PetscScalar n, TensorCorrection tensorCorrection, PetscInt MPa)
{
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
//---------------------------------------------------------------------------
void getLineStruct(
		FILE *fp, PetscInt *ls, PetscInt *le, PetscInt max_num,
		PetscInt *count_starts, PetscInt *count_ends,
		const char key[], const char key_end[])
{
	// get the positions in the file for the material structures

	char line[MAX_LINE_LEN];
	PetscInt comment, start_match, end_match;
	PetscInt i = 0, ii = 0;

	// reset to start of file
	rewind(fp);

	while( !feof(fp) )
	{
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// if line is blank
		if( strlen(line) == 0 ) { continue; }

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) { continue; }

		// find start of material structure
		start_match = material_key_matches( key, line );
		if( start_match == _TRUE )
		{
			// check bounds
			if(i+1 > max_num) { i++; return; }

			// mark file position and scan until end_match is found
			ls[i++] = (PetscInt)ftell( fp );
		}

		// find end of material structure
		end_match = material_key_matches( key_end, line );
		if( end_match == _TRUE )
		{
			// check bounds
			if(ii+1 > max_num) { ii++; return; }

			// mark file position
			le[ii++] = (PetscInt)ftell( fp );
		}

		// if no match continue
		if( (start_match == _FALSE) && (end_match == _FALSE) ) { continue; }
	}

	(*count_starts)  = i;
	(*count_ends)    = ii;
}
//---------------------------------------------------------------------------
void getMatPropInt(FILE *fp, PetscInt ils, PetscInt ile,
		const char key[], PetscInt *value, PetscInt *found)
{
	// get integer within specified positions of the file

	char line[MAX_LINE_LEN];
	PetscInt comment, pos;
	PetscInt match, int_val;

	// init flag
	if(found) (*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			int_val = (PetscInt)strtol( line, NULL, 0 );

			if(found)
				(*found) = _TRUE;
				(*value) = int_val;

			return;
		}
	}
}
//---------------------------------------------------------------------------
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile,
		const char key[], PetscScalar *value, PetscInt *found)
{
	// get scalar within specified positions of the file

	char          line[MAX_LINE_LEN];
	PetscInt      comment, pos;
	PetscInt      match;
	PetscScalar   double_val;

	// init flag
	if(found) (*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			double_val = (PetscScalar)strtod( line, NULL );

			if(found)
				(*found) = _TRUE;
				(*value) = double_val;

			return;
		}
	}
}
//---------------------------------------------------------------------------
void getMatPropString( FILE *fp, PetscInt ils, PetscInt ile, const char key[], char value[], PetscInt max_L, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment, pos;
	PetscInt match, line_L;
	char LINE[MAX_LINE_LEN];

	// init flag
	if(found) (*found) = _FALSE;

	memset( value, 0, sizeof(char)*(size_t)max_L );

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);
			strip_all_whitespace(line, LINE);

			trim_past_comment(LINE);
			line_L = (PetscInt)strlen( LINE );

			strncpy( value, LINE, (size_t)line_L );
			if( line_L > max_L ) {
				printf("parse_GetString: Error, input string is not large enough to hold result \n");
				return;
			}

			if(found) (*found) = _TRUE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
void getMatPropIntArray(FILE *fp, PetscInt ils, PetscInt ile,const char key[],
		PetscInt *nvalues, PetscInt values[], PetscInt *found)
{
	// get scalar within specified positions of the file

	char          line[MAX_LINE_LEN];
	PetscInt      comment, pos, count;
	PetscInt      match;
	PetscInt      int_val;
	char 		  *_line;

	// init flag
	if(found) (*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			count = 0;
			_line = line;
			for(;;)
			{
				char *endp;
				int_val = (PetscInt)strtod(_line, &endp);
				values[count] = int_val;

				if(endp == _line) break;

				_line = endp;
				count++;
			}

			*nvalues = count;
			*found 	 = _TRUE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
void getMatPropScalArray(FILE *fp, PetscInt ils, PetscInt ile,const char key[],
		PetscInt *nvalues, PetscScalar values[], PetscInt *found)
{
	// get scalar within specified positions of the file

	char          line[MAX_LINE_LEN];
	PetscInt      comment, pos, count;
	PetscInt      match;
	PetscScalar   double_val;
	char 		  *_line;

	// init flag
	if(found) (*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			count = 0;
			_line = line;
			for(;;)
			{
				char *endp;
				double_val = strtod(_line, &endp);
				values[count] = double_val;

				if(endp == _line) break;

				_line = endp;
				count++;
			}

			*nvalues = count;
			*found 	 = _TRUE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
