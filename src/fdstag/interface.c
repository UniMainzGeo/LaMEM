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

	PetscInt     i, viscLaw, plastLaw, numPhases, numSoft;
	PetscScalar  eta, D, n;
	Material_t  *phases;
	PhaseProps  *PhaseProperties;
	MatParLim   *matLim;
	Soft_t      *matSoft;
	Scaling     *scal;
	PetscBool   quasi_harmonic;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsHasName(PETSC_NULL, "-use_quasi_harmonic_viscosity", &quasi_harmonic); CHKERRQ(ierr);

	// access phase properties in user context variables
	numPhases       = usr->num_phases;
	PhaseProperties = &usr->PhaseProperties;
	matLim          = &jr->matLim;
	scal            = &jr->scal;

	// compute number of material softening laws
	for(i = 0, numSoft = 0; i < numPhases; i++)
	{
		if(PhaseProperties->PlasticityLaw[i])
		{
			if(PhaseProperties->CohesionAfterWeakening[i])      numSoft++;
			if(PhaseProperties->FrictionAngleAfterWeakening[i]) numSoft++;
		}
	}

	// allocate material parameters & softening laws
	ierr = PetscMalloc(sizeof(Material_t)*(size_t)numPhases, &phases); CHKERRQ(ierr);
	ierr = PetscMemzero(phases, sizeof(Material_t)*(size_t)numPhases); CHKERRQ(ierr);

	ierr = PetscMalloc(sizeof(Soft_t)*(size_t)numSoft, &matSoft);      CHKERRQ(ierr);
	ierr = PetscMemzero(matSoft, sizeof(Soft_t)*(size_t)numSoft);      CHKERRQ(ierr);

	// read phase parameters
	for(i = 0, numSoft = 0; i < numPhases; i++)
	{
		// phase identifier
		phases[i].ID = i;

		// density parameters
		phases[i].rho = PhaseProperties->rho[i];

		// viscosity parameters
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
		else if(viscLaw == 4)
		{
			// temperature-dependent power-law creep parameters
			eta = PhaseProperties->A[i];
			D   = PhaseProperties->Powerlaw_e0[i];
			n   = PhaseProperties->n_exponent[i];

			// convert power-law creep parameters to dislocation creep parameters
			phases[i].Bn = pow(2.0*eta, -n)*pow(D, 1-n);
			phases[i].n  = n;

			// store activation energy
			phases[i].En = PhaseProperties->E[i];

			// store reference strain rate
			matLim->DII_ref = D;
		}
		else
		{	// unsupported viscosity law
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unsupported viscosity law used");
		}

		// elasticity parameters
		phases[i].K = PhaseProperties->ElasticBulkModule[i];  // bulk modulus
		phases[i].G = PhaseProperties->ElasticShearModule[i]; // shear modulus

		// plasticity parameters
		plastLaw = PhaseProperties->PlasticityLaw[i];

		if(plastLaw)
		{
			// cohesion & friction angle
			phases[i].ch = PhaseProperties->Cohesion[i];
			phases[i].fr = PhaseProperties->FrictionAngle[i]/scal->angle;

			// softening laws
			if(PhaseProperties->CohesionAfterWeakening[i])
			{
				// set cohesion softening law
				matSoft[numSoft].APS1 = PhaseProperties->Weakening_PlasticStrain_Begin[i];
				matSoft[numSoft].APS2 = PhaseProperties->Weakening_PlasticStrain_End  [i];
				matSoft[numSoft].A    = 1.0 - PhaseProperties->CohesionAfterWeakening[i]/phases[i].ch;

				// store cohesion softening law
				phases[i].chSoft = &matSoft[numSoft++];
			}
			if(PhaseProperties->FrictionAngleAfterWeakening[i])
			{
				// set friction softening law
				matSoft[numSoft].APS1 = PhaseProperties->Weakening_PlasticStrain_Begin[i];
				matSoft[numSoft].APS2 = PhaseProperties->Weakening_PlasticStrain_End  [i];
				matSoft[numSoft].A    = 1.0 - (PhaseProperties->FrictionAngleAfterWeakening[i]/scal->angle)/phases[i].fr;

				// store friction softening law
				phases[i].frSoft = &matSoft[numSoft++];
			}

			if(quasi_harmonic == PETSC_TRUE)
			{
				phases[i].quasi_harmonic = 1;
			}
		}
	}

	// store material parameters & softening laws
	jr->numPhases = numPhases;
	jr->phases    = phases;
	jr->numSoft   = numSoft;
	jr->matSoft   = matSoft;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
