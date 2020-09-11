
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
 **    filename:   AVD.c
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

// Routine that changes the phase as a function of the condition set by the user.
 // These phase transitions are extremely simplified wtr to the phase transition listed in
 // phase diagram. And allow to consider lower mantle, ast->lithosphere, basalt->eclogitic reaction
 //

#include "LaMEM.h"
#include "AVD.h"
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "fdstag.h"
#include "bc.h"
#include "tools.h"
#include "phase_transition.h"
#include "phase.h"
#include "constEq.h"
#include "parsing.h"
#include "objFunct.h"
//-----------------------------------------------------------------//
PetscScalar  Compute_Melt_Fraction(PetscScalar   P, PetscScalar T ,Material_t *phase, ConstEqCtx *ctx)
{

	Scaling *scal;
	PetscScalar T_C,P_GPa;
	PetscScalar mf=0.0;

	if(ctx->ctrl->initGuess==1 || !strcmp(phase->Melt_Parametrization,"none")) return mf=0.0;

	scal  = ctx->scal;
	T_C   = T*scal->temperature-scal->Tshift;
	P_GPa = (P*scal->stress_si)/1e9;


	if(!strcmp(phase->Melt_Parametrization,"Katz_2003"))
	{
		PetscScalar A[3],B[3],C[3],BETA[2],CHI[2],KAPPA,GAMMA,r[2],LAMBDA;

		A[0]    =   1085.7  ;
		A[1]    =   132.9   ;
		A[2]    =   -5.1    ;
		B[0]    =   1475.0  ;
		B[1]    =   80.0    ;
		B[2]    =   -3.2    ;
		C[0]    =   1780.0  ;
		C[1]    =   45.0    ;
		C[2]    =   -2.0    ;
		r[0]    =   0.5     ;
		r[1]    =   0.08    ;
		BETA[0] =   1.5     ;
		BETA[1] =   1.5     ;
		KAPPA   =   43.0    ;
		GAMMA   =   0.75    ;
		CHI[0]  =   0.12    ;
		CHI[1]  =   0.01    ;
		LAMBDA  =   0.6     ;

		PetscScalar F;
		PetscScalar rho_s = 3300, rho_l=2900;
		PetscScalar Fcpx = phase->Cpx_mode/(r[0] + r[1]*P_GPa);
		PetscScalar Tsol, Tlhz, Tcpx, Tliq, dT, Cf_SAT,Cf;

		Cf = phase->Melt_Water_Par;

		Cf_SAT = CHI[0]*pow(P_GPa,LAMBDA) + CHI[1]*P_GPa;
		if(Cf_SAT<Cf) Cf=Cf_SAT;


		dT = KAPPA*pow(100.0*Cf,GAMMA);
		Tsol = A[0] + A[1]*P_GPa + A[2]*P_GPa*P_GPa;
		Tlhz = B[0] + B[1]*P_GPa + B[2]*P_GPa*P_GPa;
		Tcpx = pow(Fcpx,1.0/BETA[0])*(Tlhz - Tsol) + Tsol;
		Tliq = C[0] + C[1]*P_GPa + C[2]*P_GPa*P_GPa;

		if(T_C<Tsol-dT)
		{
		   F = 0.0;
		}
		else if(T_C<Tcpx-dT)
		{
		   F = pow(((T_C-(Tsol-dT))/(Tlhz-Tsol)),BETA[0]);
		 }
		else if(T_C<Tliq-dT)
		{
		   F = Fcpx + (1-Fcpx)*pow(((T_C-(Tcpx-dT))/(Tliq-Tcpx)),BETA[1]);
		}
		else
		{
		   F = 1.0;
		}

		mf = F/((1-F)*(rho_l/rho_s)+F);
	}
// Convert weight fraction to volumetric fraction

	return mf;



}
