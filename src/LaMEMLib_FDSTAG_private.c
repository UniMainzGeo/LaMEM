//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Assembly_FDSTAG.h"
#include "Elements.h"
#include "GetGravityField.h"
#include "GetSurfaceVelocity.h"
#include "LaMEMLib_FDSTAG_private.h"
//---------------------------------------------------------------------------
// Tobias, what should I do with this? This code was part of solution routine.
// Let's think of incorporating your gravity stuff in canonical version.
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CalculateMisfitValues"
PetscErrorCode CalculateMisfitValues(
	UserContext        *user,
	LaMEMVelPressureDA  C,
	PetscInt            itime,
	PetscScalar        *LaMEM_OutputParameters)
{
	PetscScalar OptMfit_L1, OptMfit_L2, OptMfit_ChiSquare;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Calculate gravity field + misfit

	if(user->GravityField.GetIt==1)
	{
		PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- START evaluating gravity field --------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- Check your results for reliability ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- This is a development version ---------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");

		ierr = GetGravityField(user, C->ngp_vel, itime); CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- END evaluating gravity field ----------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
	}

	//==========================================================================================
	// Extract surface velocity + calculate misfit

	if(user->SurfVelField.GetIt == 1)
	{
		PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- START extracting surface velocity -----\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- This is a development version ---------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- Check your results for reliability ----\n");
		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");

		if (user->ErosionParameters.UseInternalFreeSurface == 1)
		{
			ierr = GetSurfaceVelocity_InternalFreeSurface(user, itime);CHKERRQ(ierr);
		}
		else
		{
			ierr = GetSurfaceVelocity(user, user->DA_Vel, user->sol_advect);CHKERRQ(ierr);
		}

		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- END extracting surface velocity -------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
	}

	if(user->Optimisation.GetIt == 1)
	{
		if(user->SurfVelField.GetIt == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "==> Misfit SurfVelField (L1)/N = %g \n", user->Optimisation.SumAbsVel/user->Optimisation.NSurfVel);
			PetscPrintf(PETSC_COMM_WORLD, "==> Misfit SurfVelField (L2)/sqrt(N) = %g \n", sqrt(user->Optimisation.SumSqrsVel/user->Optimisation.NSurfVel));
		}
		else
		{
			user->Optimisation.SumAbsVel  = 0.0;
			user->Optimisation.SumSqrsVel = 0.0;
			user->Optimisation.NSurfVel   = 0.0;
		}
		if(user->GravityField.GetIt == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD,"==> Misfit GravityField (L1)/N = %g \n",user->Optimisation.SumAbsGrav/user->Optimisation.NGrav);
			PetscPrintf(PETSC_COMM_WORLD,"==> Misfit GravityField (L2)/sqrt(N) = %g \n",sqrt(user->Optimisation.SumSqrsGrav/user->Optimisation.NGrav));
		}
		else
		{
			user->Optimisation.SumAbsGrav  = 0.0;
			user->Optimisation.SumSqrsGrav = 0.0;
			user->Optimisation.NGrav       = 0.0;
		}

		// Misfit based on L1 norm (arithmetic mean)
		OptMfit_L1 = (user->Optimisation.SumAbsGrav + user->Optimisation.SumAbsVel)/(user->Optimisation.NGrav + user->Optimisation.NSurfVel);

		// Misfit based on L2 norm (rms mean)
		OptMfit_L2 = sqrt((user->Optimisation.SumSqrsGrav + user->Optimisation.SumSqrsVel)/(user->Optimisation.NGrav + user->Optimisation.NSurfVel));

		// Chi-square misfit based
		OptMfit_ChiSquare = (user->Optimisation.SumSqrsGrav + user->Optimisation.SumSqrsVel)/(user->Optimisation.NGrav + user->Optimisation.NSurfVel);

		PetscPrintf(PETSC_COMM_WORLD,"==> Misfit (L1)/N = %g \n",OptMfit_L1);
		PetscPrintf(PETSC_COMM_WORLD,"==> Misfit (L2)/sqrt(N) = %g \n",OptMfit_L2);
		PetscPrintf(PETSC_COMM_WORLD,"==> Misfit chi-square = %g \n",OptMfit_ChiSquare);

		*LaMEM_OutputParameters = OptMfit_ChiSquare;

		PetscPrintf(PETSC_COMM_WORLD,"#******************************************\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
