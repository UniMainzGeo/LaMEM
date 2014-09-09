//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Assembly_FDSTAG.h"
#include "Elements.h"
#include "GetGravityField.h"
#include "GetSurfaceVelocity.h"
#include "LaMEMLib_FDSTAG_private.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "CreateSolutionVectors"
PetscErrorCode CreateSolutionVectors(UserContext *user)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Generate and initialize global Vectors for Stokes
	ierr = DMCreateGlobalVector(user->DA_Vel, &user->sol);               CHKERRQ(ierr);

	// copy of velocity solution vector used to advect properties
	ierr = DMCreateGlobalVector(user->DA_Vel,  &user->sol_advect);       CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(user->DA_Pres, &user->Pressure);         CHKERRQ(ierr);

	// Create a vector that holds the inverse of viscosity (for scaling)
	ierr = DMCreateGlobalVector(user->DA_Pres, &user->ViscosityScaling); CHKERRQ(ierr);

	// initialize
	ierr = VecSet(user->sol,              0.0); CHKERRQ(ierr);
	ierr = VecSet(user->sol_advect,       0.0); CHKERRQ(ierr);
	ierr = VecSet(user->Pressure,         0.0); CHKERRQ(ierr);
	ierr = VecSet(user->ViscosityScaling, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DestroySolutionObjects"
PetscErrorCode DestroySolutionObjects(UserContext *user, LaMEMVelPressureDA *C)
{
	PetscInt    i;
	PetscMPIInt rank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// DM's
	ierr = DMDestroy(&user->DA_Vel);        CHKERRQ(ierr);
	ierr = DMDestroy(&user->DA_Pres);       CHKERRQ(ierr);
	ierr = DMDestroy(&user->DA_Materials);  CHKERRQ(ierr);
	ierr = DMDestroy(&user->DA_Processors); CHKERRQ(ierr);

	if (user->VTKOutputFiles == 1)
	{
		ierr = DMDestroy( &user->DA_Quadrature ); CHKERRQ(ierr);
	}
	ierr = DMDestroy(&user->DA_SurfaceTopography); CHKERRQ(ierr);
	ierr = DMDestroy(&user->DA_BottomTopography);  CHKERRQ(ierr);

	// Vectors
	ierr = VecDestroy(&user->Materials);            CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography);    CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vx); CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vy); CHKERRQ(ierr);
	ierr = VecDestroy(&user->SurfaceTopography_Vz); CHKERRQ(ierr);
	ierr = VecDestroy(&user->BottomTopography);     CHKERRQ(ierr);

	// Cleanup FD erosion code
	if ((user->ErosionParameters.ErosionModel==2) && (rank==0))
	{
		// not good! implement erosion-code-data-structure create and destroy functions
		ierr = DMDestroy (&user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode);CHKERRQ(ierr);
		ierr = VecDestroy(&user->ErosionParameters.FE_ErosionCode.ErosionSurface);   CHKERRQ(ierr);
	}

	for (i=0; i<user->num_phases; i++)
	{
		ierr = VecDestroy(&user->FDSTAG.Center_PhaseProportions[i]);		 CHKERRQ(ierr);
		ierr = VecDestroy(&user->FDSTAG.Corner_PhaseProportions[i]);		 CHKERRQ(ierr);
		ierr = VecDestroy(&user->FDSTAG.Corner_PhaseProportions_local[i]);CHKERRQ(ierr);
		ierr = VecDestroy(&user->FDSTAG.XYPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
		ierr = VecDestroy(&user->FDSTAG.XZPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
		ierr = VecDestroy(&user->FDSTAG.YZPoints_PhaseProportions[i]);	 CHKERRQ(ierr);
	}
	ierr = VecDestroy(&user->FDSTAG.Center_Temperature);                  CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_Pressure);                     CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_Strain);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_PlasticStrain);                CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_T2nd);                         CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_E2nd);                         CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_EffectiveViscosity);           CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_Density);                      CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Center_NumParticles);                 CHKERRQ(ierr);

	ierr = VecDestroy(&user->FDSTAG.Corner_Temperature);                  CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Corner_Pressure);                     CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.Corner_Density);                      CHKERRQ(ierr);

	ierr = VecDestroy(&user->FDSTAG.XYPoints_Temperature);                CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_Pressure);                   CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_Strain);                     CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_PlasticStrain);              CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_T2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_E2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_EffectiveViscosity);         CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XYPoints_Density);                    CHKERRQ(ierr);

	ierr = VecDestroy(&user->FDSTAG.XZPoints_Temperature);                CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_Pressure);                   CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_Strain);                     CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_PlasticStrain);              CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_T2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_E2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_EffectiveViscosity);         CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.XZPoints_Density);                    CHKERRQ(ierr);

	ierr = VecDestroy(&user->FDSTAG.YZPoints_Temperature);                CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_Pressure);                   CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_Strain);                     CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_PlasticStrain);              CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_T2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_E2nd);                       CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_EffectiveViscosity);         CHKERRQ(ierr);
	ierr = VecDestroy(&user->FDSTAG.YZPoints_Density);                    CHKERRQ(ierr);

	ierr = DMDestroy(&user->FDSTAG.DA_CENTER);                            CHKERRQ(ierr);
	ierr = DMDestroy(&user->FDSTAG.DA_CORNER);                            CHKERRQ(ierr);
	ierr = DMDestroy(&user->FDSTAG.DA_XY_POINTS);                         CHKERRQ(ierr);
	ierr = DMDestroy(&user->FDSTAG.DA_XZ_POINTS);                         CHKERRQ(ierr);
	ierr = DMDestroy(&user->FDSTAG.DA_YZ_POINTS);                         CHKERRQ(ierr);

	// particles
	ierr = PetscFree(user->ParticlesLocal); CHKERRQ(ierr);

	ierr = VecDestroy(&user->Pressure);        CHKERRQ(ierr);

	ierr = VecDestroy(&user->ViscosityScaling); CHKERRQ(ierr);
	ierr = VecDestroy(&user->sol);              CHKERRQ(ierr);
	ierr = VecDestroy(&user->sol_advect);       CHKERRQ(ierr);

	// Matrices
	ierr = MatDestroy(&user->VV_MAT); CHKERRQ(ierr);
	ierr = MatDestroy(&user->PP_MAT); CHKERRQ(ierr);

	ierr = MatDestroy(&user->VP_MAT); CHKERRQ(ierr);
	ierr = MatDestroy(&user->PV_MAT); CHKERRQ(ierr);

	ierr = MatDestroy(&user->approx_S); CHKERRQ(ierr);

	ierr = LaMEMVelPressureDADestroy(C); CHKERRQ(ierr);

	ierr = PetscFree(user->ParticlesLocal);    CHKERRQ(ierr);
	ierr = PetscFree(user->TimeDependentData); CHKERRQ(ierr);

	if(user->InputParamFile)
	{
		MaterialDestroy(user->PhaseMaterialProperties);
	}

	PetscFunctionReturn(0);
}
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
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCompPrecond"
PetscErrorCode FDSTAGCompPrecond(
		Mat VV_MAT, Mat VP_MAT, Mat PV_MAT, Mat PP_MAT, Mat approx_S,
		Vec ViscosityScaling, PetscScalar *nrmVV, UserContext *user)
{
	PetscBool		eliminate_stickyair_from_system, subtract_volumetric_strainrates;
	PetscInt        nel_x, nel_y, nel_z, nnode_x, nnode_y, nnode_z, number_coefficients;
	PetscInt        i,j,k,xm,ym,zm,xs,ys,zs;
	PetscInt        iel_x, iel_y, iel_z;
	PetscInt        xmp,ymp,zmp,xsp,ysp,zsp,zs_FreeSurface;
	PetscInt 		*rowid_DirichletBC_array, *rowid_addPoints_array, numDirichlet, numAddpoints;
	PetscInt        AirPhase;
	const PetscInt *P_globalindices, *Vel_globalindices;
	PetscInt        FreeSurfaceCells[7];
	PetscScalar     dx_P[2], dy_P[2], dz_P[2], dx_vec[7], dy_vec[7], dz_vec[7], Eta_Center[7], Eta_XY[4], Eta_XZ[4], Eta_YZ[4], dRho_dxdydz[3],z_FreeSurface;
	PetscScalar     ***viscosity_center, ***viscosity_XY,  ***viscosity_XZ,  ***viscosity_YZ, ***density_center, ***PhaseProportionsAir_Center, ***LocalSurfaceTopography;
	PetscScalar 	z_center, FreeSurface_Fraction;
	PetscScalar 	ScalingParameterDiagonal;
	DM              cda, cda_SurfaceTopo;
	Vec             gc, Viscosity_Center_local, Viscosity_XY_local, Viscosity_XZ_local, Viscosity_YZ_local, Density_Center_local, PhaseProportionsAir_Center_Vec;
	Vec 			LocalSurfaceTopography_vec, gc_SurfaceTopo;
	DMDACoor3d      ***coords;
	MatStencil      row_pp[1];            // indices for pp matrix
	PetscScalar     v_pp[1];              // values  for pp matrix
	MatStencil      col_pv[6];            // indices for pv matrix
	PetscScalar     v_pv[6];              // values  for pv matrix
	MatStencil      col_vp[2];            // indices for vp matrix
	PetscScalar     v_vp[2];              // values  for vp matrix
	MatStencil      row_vv[1],col_vv[32]; // indices for vv matrix
	PetscScalar     v_vv[32];             // values  for vv matrix
	PetscInt        rowidx, vgidx[32],pgidx[2];
	DMDACoor3d		***coors_SurfaceTopo;
	DM              da, da_pres;
	ISLocalToGlobalMapping P_ltogm, Vel_ltogm;

	PetscErrorCode  ierr;
	PetscFunctionBegin;

	// set DM objects
	da      = user->DA_Vel;
	da_pres = user->DA_Pres;

	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);

	/* print some info */
	ierr = DMDAGetInfo(user->DA_Processors, 0, &nel_x, &nel_y, &nel_z,	0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr); 	// # of elements in all directions
	ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y, &nnode_z,	0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  // # of nodes 	 in all directions
	PetscPrintf(PETSC_COMM_WORLD,"#  Forming FD stiffness matrix mx,my,mz=(%lld, %lld, %lld) ...  \n",(LLD)nnode_x,(LLD)nnode_y,(LLD)nnode_z);

	/* Corners */
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

	/* Obtain mapping from current proc to global numbering for both velocity and pressure */
//	ierr = DMDAGetGlobalIndices(da_pres, &np,   &P_globalindices);   CHKERRQ(ierr);
//	ierr = DMDAGetGlobalIndices(da,      &nvel, &Vel_globalindices); CHKERRQ(ierr);

	ierr = DMGetLocalToGlobalMapping(da_pres, &P_ltogm);   CHKERRQ(ierr);
	ierr = DMGetLocalToGlobalMapping(da,      &Vel_ltogm); CHKERRQ(ierr);

	ierr = ISLocalToGlobalMappingGetIndices(P_ltogm,   &P_globalindices); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(Vel_ltogm, &Vel_globalindices); CHKERRQ(ierr);

	/* Get viscosity @ center points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER,&Viscosity_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_EffectiveViscosity,INSERT_VALUES,Viscosity_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_EffectiveViscosity,INSERT_VALUES,Viscosity_Center_local); CHKERRQ(ierr);

	/* Get density @ center points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER,&Density_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_Center_local); CHKERRQ(ierr);

	/* Get viscosity @ corner points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_XY_POINTS,&Viscosity_XY_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_XY_POINTS,user->FDSTAG.XYPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XY_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_XY_POINTS,user->FDSTAG.XYPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XY_local); CHKERRQ(ierr);

	ierr = DMGetLocalVector     (user->FDSTAG.DA_XZ_POINTS,&Viscosity_XZ_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_XZ_POINTS,user->FDSTAG.XZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XZ_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_XZ_POINTS,user->FDSTAG.XZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XZ_local); CHKERRQ(ierr);

	ierr = DMGetLocalVector     (user->FDSTAG.DA_YZ_POINTS,&Viscosity_YZ_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_YZ_POINTS,user->FDSTAG.YZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_YZ_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_YZ_POINTS,user->FDSTAG.YZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_YZ_local); CHKERRQ(ierr);

	/* Get arrays for viscosity @ center, XY, XZ & YZ points */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,		Viscosity_Center_local, &viscosity_center ); CHKERRQ(ierr); // Viscosity @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,		Density_Center_local,   &density_center ); CHKERRQ(ierr); 	// Density @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XY_POINTS, 	Viscosity_XY_local,     &viscosity_XY     ); CHKERRQ(ierr); // Viscosity @ Sxy points
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XZ_POINTS, 	Viscosity_XZ_local,     &viscosity_XZ     ); CHKERRQ(ierr); // Viscosity @ Sxz points
	ierr = DMDAVecGetArray(user->FDSTAG.DA_YZ_POINTS, 	Viscosity_YZ_local,     &viscosity_YZ     ); CHKERRQ(ierr); // Viscosity @ Syz points

	/* In case we use an internal free surface, we need this: */
	if (user->ErosionParameters.UseInternalFreeSurface==1){
		AirPhase = user->ErosionParameters.StickyAirPhase;
	}
	else{
		AirPhase = 0;
	}
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER, &PhaseProportionsAir_Center_Vec);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER, user->FDSTAG.Center_PhaseProportions[AirPhase],INSERT_VALUES,PhaseProportionsAir_Center_Vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER, user->FDSTAG.Center_PhaseProportions[AirPhase],INSERT_VALUES,PhaseProportionsAir_Center_Vec); CHKERRQ(ierr);
	ierr = DMDAVecGetArray		(user->FDSTAG.DA_CENTER, PhaseProportionsAir_Center_Vec, 	&PhaseProportionsAir_Center	 ); CHKERRQ(ierr);

	/* Get free surface information [to correct stencil for internal free surface] */
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,	 &cda_SurfaceTopo);   CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography, &gc_SurfaceTopo);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,	 &coors_SurfaceTopo); CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_SurfaceTopography,PETSC_NULL,PETSC_NULL,&zs_FreeSurface,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);

	/* Do we eliminate sticky air or not? */
	eliminate_stickyair_from_system = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL, "-eliminate_stickyair_from_system", &eliminate_stickyair_from_system,PETSC_NULL ); CHKERRQ(ierr);

	// Do we extract the volumetric component of strain rates or not (Exx_dev = Exx-1/3*(Exx+Eyy+Ezz), or Exx=Exx)
	// YES! ALWAYS!
	subtract_volumetric_strainrates = PETSC_TRUE;
//	ierr = PetscOptionsGetBool( PETSC_NULL, "-subtract_volumetric_strainrates", &subtract_volumetric_strainrates,PETSC_NULL ); CHKERRQ(ierr);

	if (subtract_volumetric_strainrates){
		number_coefficients = 32;			// we have more equations as we subtract the volumetric components
	}
	else{
		number_coefficients = 20;
	}

	/* Make a loop over all local elements, construct the element stiffness matrix */
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++) {
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++) {
			for (iel_x=xsp; iel_x<xsp+xmp; iel_x++) {
				i = iel_x; j = iel_y; k = iel_z;		// Initialize variables

				/* Extract coordinates and spacing of the local element & neighboring elements
				 * (required in order to form the FD discretization)
				 */
				ierr = FDSTAG_ComputeSpacing(coords, iel_x,iel_y,iel_z, user, dx_vec, dy_vec, dz_vec, dx_P, dy_P, dz_P, &z_center); CHKERRQ(ierr);

				/* Extract material parameters */
				ierr = FDSTAG_ExtractMaterialParameters(viscosity_center, viscosity_XY, viscosity_YZ, viscosity_XZ, density_center, PhaseProportionsAir_Center, LocalSurfaceTopography,
														iel_x,iel_y,iel_z,zs_FreeSurface,dx_P, dy_P, dz_P, user, Eta_Center, Eta_XY, Eta_YZ, Eta_XZ, dRho_dxdydz, FreeSurfaceCells, &z_FreeSurface); CHKERRQ(ierr);

				/* If we have an internal free surface, compute the height of that relatve to the center pressure node */
				FreeSurface_Fraction=1.0;
				if ((FreeSurfaceCells[1]==0) && (FreeSurfaceCells[6]==1) && (eliminate_stickyair_from_system)){
					// Current cell is below FS & cell above is above FS
					FreeSurface_Fraction = (z_FreeSurface-z_center)/dz_P[1];
				}

				/* The nomenclature for the global stiffness matrixes is:
				 *
				 * | VV_MAT VP_MAT | | V |    |Rhs |
				 * |        	   | |	 |  = |    |
				 * | PV_MAT PP_MAT | | P |    | 0  |
				 *
				 */

				/* (1) Create FD stencil for the incompressibility equation
				 *
				 * 		(Vx(i+1)-Vx(i))/dx + (Vy(j+1)-Vy(j))/dy + (Vz(k+1)-Vz(k))/dz = 1/gamma*P(i,j,k)
				 *
				 *	BC's are not imposed on this equation, so we don't have to worry about it.
				 */
				if ( (i < (nnode_x-1)) && (j < (nnode_y-1)) && (k < (nnode_z-1)) ) {
					ierr = FDSTAG_FDStencil_Incompressibility(i,j,k,dx_vec, dy_vec, dz_vec, row_pp,v_pp,col_pv, v_pv); CHKERRQ(ierr);

					/* Find global indices */
					ierr = DAGetGlobalIndex(da_pres,row_pp[0].i,row_pp[0].j,row_pp[0].k,row_pp[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,6,col_pv,0,PETSC_NULL,vgidx,pgidx); CHKERRQ(ierr);

					/* Add pressure BC's for the free surface */
					if ((eliminate_stickyair_from_system) && (FreeSurfaceCells[1]==1)){
						v_pv[0] = v_pv[1] = v_pv[2] = v_pv[3] = v_pv[4] = v_pv[5] = 0;		// add a constant pressure BC here
					}

					ierr = MatSetValues(PV_MAT,1,&rowidx,6,vgidx,v_pv,ADD_VALUES); CHKERRQ(ierr);

					/* Create the approximate Schur matrix, with 1/eta on the diagonal */
					v_pp[0] = -1.0/Eta_Center[1];
					ierr = MatSetValue(approx_S,rowidx,rowidx,v_pp[0],INSERT_VALUES); CHKERRQ(ierr);

					/* Store vector with viscosity scaling */
					ierr = VecSetValue(ViscosityScaling,rowidx,1/Eta_Center[1],INSERT_VALUES); CHKERRQ(ierr);

					/* Create the PP matrix, with 0 on the diagonal */
					if (user->StokesSolver==1){
						v_pp[0] = -1;		// powell-hesteness scales it later
					}
					else {
						v_pp[0] = 0;		// as the system is incompressible
					}

					/* Add pressure BC's for the free surface */
					if ((eliminate_stickyair_from_system) && (FreeSurfaceCells[1]==1)){
						v_pp[0] = 1.0;		// add a constant pressure BC for sticky-air cells
					}
					ierr = MatSetValue(PP_MAT,rowidx,rowidx,v_pp[0],INSERT_VALUES); CHKERRQ(ierr);
				}

				/* (2)  Create FD stencil for the first force balance equation:
				 * -(P(i+1)-P(i))/dx + (Txx(i+1)-Txx(i))/dx + (Txy(j+1)-Txy(j))/dy + (Txz(k+1)-Txz(k))/dz = 0
				 *
				 *	This equation is located @ Vx nodes.
				 */
				if ((i > 0) && (i < (nnode_x-1)) ) {
					/* Compute FD stencil (w/out BCs) */
					ierr = FDSTAG_FDStencil_ForceBalance1(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XY, Eta_XZ, row_vv, col_vv, v_vv, col_vp,  v_vp); CHKERRQ(ierr);

					/* Correct stencil for BCs (front/back/lower/upper) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance1(j,k,nel_y,nel_z, v_vv, user); CHKERRQ(ierr);

					/* Correct stencil for internal free surface */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);
				}

				/* (3)  Create FD stencil for the second force balance equation:
				 * -(P(j+1)-P(j))/dy + (Txy(i+1)-Txy(i))/dx + (Tyy(j+1)-Tyy(j))/dy + (Tyz(k+1)-Tyz(k))/dz = 0
				 *
				 * This equation is located @ Vy nodes.
				 */
				if ((j > 0) && (j < (nnode_y-1))) {
					FDSTAG_FDStencil_ForceBalance2(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XY, Eta_YZ,	row_vv, col_vv, v_vv, col_vp,  v_vp);

					/* Correct stencil for BCs (left/right/lower/upper) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance2(i,k,nel_x,nel_z, v_vv, user); CHKERRQ(ierr);

					/* Global indices */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

				}

				/* (4)  Create FD stencil for the third force balance equation:
				 * -(P(k)-P(k-1))/dz + (Txz(i+1)-Txz(i))/dx + (Syz(j+1)-Syz(j))/dy + (Tzz(k+1)-Tzz(k))/dz = 0
				 *
				 * This equation is located @ Vz nodes
				 */
				if ( (k > 0) && (k < (nnode_z-1)) ) {

					FreeSurface_Fraction=1.0;
					FDSTAG_FDStencil_ForceBalance3(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XZ, Eta_YZ, row_vv, col_vv, v_vv, col_vp,  v_vp, FreeSurface_Fraction);

					/* Correct stencil for BCs (left/right/front/back) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance3(i,j,nel_x,nel_y, v_vv, v_vp, FreeSurfaceCells, user, eliminate_stickyair_from_system); CHKERRQ(ierr);

					/* Find global indices */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

				}

			}
		}
	}

	ScalingParameterDiagonal = 1.0; // The value set on the diagonal; it should be on

	// Allocate array that will contain indices of the boundary conditions
	PetscMalloc((size_t)(2*xmp*ymp + 2*ymp*zmp + 2*xmp*zmp)*sizeof(PetscInt),&rowid_DirichletBC_array);

	/* Set Dirichlet BCs for velocity normal to the side boundaries
	 *
	 * Note:  we do set "dummy" values here as PETSC requires the full matrix to be assembled in order to symmetrize it.
	 * Yet, these values will be overwritten at a later stage once we symmetrize the matrix
	 * */
	numDirichlet=0;
	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			for (i=xs; i<xs+xm; i++) {

				if ( ((i == 0) || (i == nnode_x-1)) && (j < nnode_y-1) && (k < nnode_z-1) ) {
					/* Left or right boundary [Vx=constant] */
					row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
					v_vv[0] = ScalingParameterDiagonal;
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

					rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
					numDirichlet 							=	numDirichlet + 1;

				}

				if ( ((j == 0) || (j == nnode_y-1)) && (i < nnode_x-1) && (k < nnode_z-1)) {
					/* Front or back boundary [Vy=constant]*/
					row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
					v_vv[0] = ScalingParameterDiagonal;
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

					rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
					numDirichlet 							=	numDirichlet + 1;

				}

				if ( (k == nnode_z-1) && (i < nnode_x-1) && (j < nnode_y-1) ) {
					/* Upper boundary */

					if (user->BC.UpperBound == 1) {
						//  Free slip: Vz=0
						row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
						v_vv[0] = ScalingParameterDiagonal;
						ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
						ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

						rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
						numDirichlet 							=	numDirichlet + 1;

					} else if (user->BC.UpperBound ==0 ) {
						// stress-free; Szz=0

						if ( (j == 0) || (j == nnode_y-1) || (i == 0) || (i == nnode_x-1) ) {
							// nodes at the side boundaries (not sure whether this is necessary)
							row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
							v_vv[0] = ScalingParameterDiagonal;
							//					ierr = MatSetValuesStencil(VV_MAT,1,row_vv,1,row_vv,v_vv,ADD_VALUES);		CHKERRQ(ierr);
							ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
							ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

							rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
							numDirichlet 							=	numDirichlet + 1;

						} else {
							// center nodes

							FreeSurface_Fraction=1.0;
							FDSTAG_FDStencil_ForceBalance3(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XZ, Eta_YZ, row_vv, col_vv, v_vv, col_vp,  v_vp, FreeSurface_Fraction);

							/* Correct stencil for Stress-free upper BC */
							// delete coeffs by hand
							v_vp[1] 	= 0;				//	P(iz+1/2)=0
							v_vv[0] 	= 0;				// 	tau_zz(iz+1/2)=0
							v_vv[1] 	= 0;
							col_vv[0] = col_vv[1];		// 	Put 'fake' coefficients as the routines below otherwise complain that it is outside the domain
							col_vp[1] = col_vp[0];

							/* Global indices */
							ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
							ierr = FDSTAG_GetStencilValues(da,da_pres,20,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

							ierr = MatSetValues(VV_MAT,1,&rowidx,20,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
							ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

							rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
							numDirichlet 							=	numDirichlet + 1;

						}


					}
				} else if ( (k == 0) && (i < nnode_x-1) && (j < nnode_y-1) ) {
					/* Lower or upper boundary [Vz=constant]*/
					row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
					v_vv[0] 								= 	ScalingParameterDiagonal;
					ierr 									= 	DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr 									=	MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
					rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
					numDirichlet 							=	numDirichlet + 1;

				}

			}
		}
	}


	/* The FDSTAG formulation, in combination with the DA, results in some velocity points that are not used in the discretization.
	 * Yet, the matrix should not contain any rows that are zero.
	 * Therefore, we have to set the diagonal of the matrix to 1 for these velocity points */
	ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  // # of nodes 	 in all directions
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

	// Allocate array that will contain indices of the boundary conditions
	PetscMalloc((size_t)(2*xm*ym + 2*ym*zm + 2*xm*zm)*sizeof(PetscInt),&rowid_addPoints_array);

	numAddpoints = 0;

	/* Set additional Vx points */
	for (j=ys; j<ys+ym; j++) {
		for	(i=xs; i<xs+xm; i++) {
			if (zs+zm == nnode_z) { // the local processor borders this boundary
				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;

			}
		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (i=xs; i<xs+xm; i++) {
			if (ys+ym == nnode_y) {
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	/* Set additional Vy points */
	for (j=ys; j<ys+ym; j++) {
		for (i=xs; i<xs+xm; i++) {
			if (zs+zm == nnode_z) {
				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
				v_vv[0] = 1.0;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			if (xs+xm == nnode_x) {
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	/* Set additional Vz points */
	for (k=zs; k<zs+zm; k++) {
		for (i=xs; i<xs+xm; i++) {
			if (ys+ym == nnode_y) {
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}

		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			if (xs+xm == nnode_x) {
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx);CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}
	// end of setting Dirichlet BC

	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,      Viscosity_Center_local, &viscosity_center); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,      Density_Center_local,   &density_center); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XY_POINTS,   Viscosity_XY_local,     &viscosity_XY);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XZ_POINTS,   Viscosity_XZ_local,     &viscosity_XZ);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_YZ_POINTS,   Viscosity_YZ_local,     &viscosity_YZ);     CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,  	PhaseProportionsAir_Center_Vec, &PhaseProportionsAir_Center); CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&Viscosity_Center_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&Density_Center_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_XY_POINTS	,&Viscosity_XY_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_XZ_POINTS,	&Viscosity_XZ_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_YZ_POINTS,	&Viscosity_YZ_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&PhaseProportionsAir_Center_Vec);	CHKERRQ(ierr);

	/* Finalize the matrices */
	ierr = MatAssemblyBegin(VV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global velocity stiffness matrix
	ierr = MatAssemblyEnd  (VV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global velocity stiffness matrix

	ierr = MatAssemblyBegin(PP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global pressure mass matrix
	ierr = MatAssemblyEnd  (PP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global pressure mass matrix

	ierr = MatAssemblyBegin(VP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global VP matrix
	ierr = MatAssemblyEnd  (VP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global VP matrix

	ierr = MatAssemblyBegin(PV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global PV matrix
	ierr = MatAssemblyEnd  (PV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global PV matrix

	ierr = MatAssemblyBegin(approx_S,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble schur preconditioning matrix
	ierr = MatAssemblyEnd  (approx_S, MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr); // Assemble schur preconditioning matrix

	ierr = VecAssemblyBegin(ViscosityScaling); CHKERRQ(ierr);	// Assemble ViscosityScaling
	ierr = VecAssemblyEnd  (ViscosityScaling); CHKERRQ(ierr);  // Assemble ViscosityScaling

	//coordinates
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);

	/* Cleaning up */
	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);

	/* 				Make the VV and VP matrix symmetric
	 *
	 * 	Doing this requires us to know on which rows we have set BC's; we have stored that above for
	 * 		rowid_DirichletBC_array - Dirichlet BC's  (essentially normal velocities) to the boundary
	 * 		rowid_addPoints_array 	- Additional points added since we use a DMDA with 3 dof
	 *
	 * */

	IS is_dir, is_add;

	// create index sets
	ierr = ISCreateGeneral(PETSC_COMM_WORLD, numDirichlet, rowid_DirichletBC_array, PETSC_COPY_VALUES, &is_dir); CHKERRQ(ierr);
	ierr = ISCreateGeneral(PETSC_COMM_WORLD, numAddpoints, rowid_addPoints_array,   PETSC_COPY_VALUES, &is_add); CHKERRQ(ierr);

	// correct VP_MAT
	ierr = MatZeroRowsIS(VP_MAT, is_dir, 0.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
	ierr = MatZeroRowsIS(VP_MAT, is_add, 0.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

	// correct VV_MAT & and put zero on diagonal to compute unbiased matrix norm
	ierr = MatZeroRowsColumnsIS(VV_MAT, is_dir, 0.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
	ierr = MatZeroRowsColumnsIS(VV_MAT, is_add, 0.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

	// get unbiased norm
	ierr = MatNorm(VV_MAT, NORM_INFINITY, nrmVV); CHKERRQ(ierr);

	// correct VV_MAT again & put units on diagonal
	ierr = MatZeroRowsColumnsIS(VV_MAT, is_dir, 1.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
	ierr = MatZeroRowsColumnsIS(VV_MAT, is_add, 1.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

	// Clean up
	ierr = PetscFree(rowid_DirichletBC_array); CHKERRQ(ierr);
	ierr = PetscFree(rowid_addPoints_array);   CHKERRQ(ierr);
	ierr = ISDestroy(&is_dir);                 CHKERRQ(ierr);
	ierr = ISDestroy(&is_add);                 CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
