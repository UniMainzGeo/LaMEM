//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Assembly_FDSTAG.h"
#include "Elements.h"
#include "GetGravityField.h"
#include "GetSurfaceVelocity.h"
#include "LaMEMLib_FDSTAG_private.h"
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "ViewLinearSolveResidual"
PetscErrorCode ViewLinearSolveResidual(UserContext *user)
{
	PetscInt    loc;
	PetscScalar val;
	Vec 	    Div;
	Vec         mRes, prod;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// compute divergence residual
	ierr = VecDuplicate(user->Pressure, &Div);  CHKERRQ(ierr);
	ierr = MatMult(user->PV_MAT,user->sol,Div); CHKERRQ(ierr);
	ierr = VecAXPY(Div, -1.0, user->rhs_p);     CHKERRQ(ierr);
	ierr = VecAbs(Div);                         CHKERRQ(ierr);

	// compute momentum residual
	ierr = VecDuplicate(user->rhs, &mRes);              CHKERRQ(ierr);
	ierr = VecCopy(user->rhs, mRes);                    CHKERRQ(ierr);
	ierr = VecDuplicate(user->rhs, &prod);              CHKERRQ(ierr);
	ierr = MatMult(user->VV_MAT, user->sol, prod);      CHKERRQ(ierr);
	ierr = VecAXPY(mRes, -1.0, prod);                   CHKERRQ(ierr);
	ierr = MatMult(user->VP_MAT, user->Pressure, prod); CHKERRQ(ierr);
	ierr = VecAXPY(mRes, -1.0, prod);                   CHKERRQ(ierr);
	ierr = VecAbs(mRes);                                CHKERRQ(ierr);

	PetscPrintf( PETSC_COMM_WORLD, "------------------------------------------\n" );
	PetscPrintf( PETSC_COMM_WORLD, "Stokes solver summary: \n" );

	PetscPrintf( PETSC_COMM_WORLD, "  divergence: \n" );
	ierr = VecMin( Div, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    Div_min = %12.12e \n", val );
	ierr = VecMax( Div, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    Div_max = %12.12e \n", val );
	ierr = VecNorm( Div, NORM_2, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |Div|_2 = %12.12e \n", val );
	ierr = VecNorm( Div, NORM_1, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |Div|_1 = %12.12e \n", val );

	PetscPrintf( PETSC_COMM_WORLD, "  momentum: \n" );
	ierr = VecMin( mRes, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    mRes_min = %12.12e \n", val );
	ierr = VecMax( mRes, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    mRes_max = %12.12e \n", val );
	ierr = VecNorm( mRes, NORM_2, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |mRes|_2 = %12.12e \n", val );
	ierr = VecNorm( mRes, NORM_1, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |mRes|_1 = %12.12e \n", val );

	PetscPrintf( PETSC_COMM_WORLD, "  velocity: \n" );
	ierr = VecMin(user->sol, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    u_min   = %12.12e \n", val );
	ierr = VecMax(user->sol, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    u_max   = %12.12e \n", val );
	ierr = VecNorm(user->sol, NORM_2, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |u|_2   = %12.12e \n", val );
	ierr = VecNorm(user->sol, NORM_1, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |u|_1   = %12.12e \n", val );

	PetscPrintf( PETSC_COMM_WORLD, "  pressure: \n" );
	ierr = VecMin(user->Pressure, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    p_min   = %12.12e \n", val );
	ierr = VecMax(user->Pressure, &loc, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    p_max   = %12.12e \n", val );
	ierr = VecNorm(user->Pressure, NORM_2, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |p|_2   = %12.12e \n", val );
	ierr = VecNorm(user->Pressure, NORM_1, &val ); CHKERRQ(ierr);
	PetscPrintf( PETSC_COMM_WORLD, "    |p|_1   = %12.12e \n", val );
	PetscPrintf( PETSC_COMM_WORLD, "------------------------------------------\n" );

	ierr = VecDestroy(&Div);  CHKERRQ(ierr);
	ierr = VecDestroy(&mRes); CHKERRQ(ierr);
	ierr = VecDestroy(&prod); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
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
#undef __FUNCT__
#define __FUNCT__ "CalculateTimeStep"
PetscErrorCode CalculateTimeStep(UserContext *user, PetscInt itime)
{

	PetscInt    mx, my, mz;
	PetscScalar MaxVel, MinVel, dx, dy, dz, dmin, fac_time;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//-------------------------------------
	// compute length of the next time step
	//-------------------------------------

	ierr     = VecMax(user->sol_advect, PETSC_NULL, &MaxVel); CHKERRQ(ierr);
	ierr     = VecMin(user->sol_advect, PETSC_NULL, &MinVel); CHKERRQ(ierr);
	MaxVel   = PetscMax(MaxVel, PetscAbsScalar(MinVel));
	ierr     = DMDAGetInfo(user->DA_Vel, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	dx       = user->W/((double) mx);
	dy       = user->L/((double) my);
	dz       = user->H/((double) mz);
	dmin     = PetscMin(dx,dy);
	dmin     = PetscMin(dmin,dz);
	fac_time = 1.0;
	if(itime < 10)
	{
		fac_time = ((PetscScalar)(itime+1))/((PetscScalar)10);
	}

	// compute Courant time step
	user->dt = user->CFL*fac_time*dmin/MaxVel;

	// limit time step
	if(user->dt > user->dt_max) user->dt = user->dt_max;

	//--------------------------------
	// print velocity & time step info
	//--------------------------------

	if(user->Characteristic.Length > 1.0)
	{
		PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user->Characteristic.Velocity*user->Characteristic.cmYear);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [cm/year]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD," Maximum velocity  %g  ",MaxVel*user->Characteristic.Velocity);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [m/s]  \n");	} else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}

	if(user->Characteristic.Length > 1.0)
	{
		// Most likely a setup that runs in natural length scales
		PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ",user->time*user->Characteristic.Time/user->Characteristic.SecYear, user->dt*user->Characteristic.Time/user->Characteristic.SecYear);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [Years]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}
	else
	{
		// Lab time scale or non-dimensional units
		PetscPrintf(PETSC_COMM_WORLD," Time = %g, dt=%g ", user->time*user->Characteristic.Time, user->dt*user->Characteristic.Time);
		if (user->DimensionalUnits == 1) { PetscPrintf(PETSC_COMM_WORLD," [s]  \n"); } else { PetscPrintf(PETSC_COMM_WORLD,"  \n"); }
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CheckVelocityError"
PetscErrorCode CheckVelocityError(UserContext *user)
{
	PetscScalar MaxVel, MinVel;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Error detection
	ierr = VecMax(user->sol_advect, PETSC_NULL, &MaxVel);	CHKERRQ(ierr);
	ierr = VecMin(user->sol_advect, PETSC_NULL, &MinVel); CHKERRQ(ierr);
	MaxVel = PetscMax(MaxVel, PetscAbsScalar(MinVel));

	if(isnan(MaxVel))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "  *** Emergency stop! Maximum velocity is NaN ***  \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
