//---------------------------------------------------------------------------
//..............   LaMEM - FDSTAG CANONICAL INTERFACE ROUTINES   ............
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interface.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InitMaterialProps"
PetscErrorCode InitMaterialProps(JacRes *jr, UserContext *usr)
{
	// initialize material properties in the FDSTAG data structures

	PetscInt    i, viscLaw;
	PetscScalar eta, D, n;
	Material_t *phases;
	PhaseProps *PhaseProperties;
	MatParLim  *matLim;

	PetscFunctionBegin;

	// access phase properties in user context variables
	phases          = jr->phases;
	PhaseProperties = &usr->PhaseProperties;
	matLim          = &jr->matLim;

	// clear phase properties array
	memset(phases, 0, sizeof(Material_t)*(size_t)jr->numPhases);

	// copy phase parameters
	for(i = 0; i < jr->numPhases; i++)
	{
		// density and power law parameters
		phases[i].ID     = i;
		phases[i].rho    = PhaseProperties->rho[i];
		phases[i].frSoft = NULL;
		phases[i].chSoft = NULL;

		// get viscosity law
		viscLaw = PhaseProperties->ViscosityLaw[i];

		if(viscLaw == 1)
		{
			// set diffusion creep constant
			phases[i].Bd = 1.0/(2.0*PhaseProperties->mu[i]);
		}
		else if(viscLaw == 2)
		{
			// power-law creep parameters
			eta = PhaseProperties->mu[i];
			D   = PhaseProperties->Powerlaw_e0[i];
			n   = PhaseProperties->n_exponent[i];

			// convert power-law creep parameters to dislocation creep parameters
			phases[i].Bn = pow(2.0*eta, -n)*pow(D, 1-n);
			phases[i].n  = n;

			// store reference strain rate
			matLim->DII_ref = D;
		}
		else
		{	// unsupported viscosity law
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Unsupported viscosity law used");
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
