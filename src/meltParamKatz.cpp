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
 **    filename:   meltParamKatz.cpp
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
//...................   PARAMETERIZED MELT FRACTION   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "meltParamKatz.h"

// The original code for Katz et al., 2003, G3. A new parameterization of hydrous mantle melting.
// Modified to solve the problem of solidus line curl back at large pressure.
// here I use a linear relationship between the solidus and pressure when pressure larger than
// a critical value Pc instead of the quadratic polynomials

// global value: pressure cutoff
PetscScalar Pc = 10.0;

//---------------------------------------------------------------------------
void setMeltParamsToDefault_Katz(meltPar_Katz *mp)
{
	// Initialize a meltPar_Katz struct to the default values

	mp->A1=DEFAULT_A1; mp->A2=DEFAULT_A2; mp->A3=DEFAULT_A3;
	mp->B1=DEFAULT_B1; mp->B2=DEFAULT_B2; mp->B3=DEFAULT_B3;
	mp->C1=DEFAULT_C1; mp->C2=DEFAULT_C2; mp->C3=DEFAULT_C3;
	mp->r1=DEFAULT_R1;           mp->r2=DEFAULT_R2;
	mp->beta1=DEFAULT_BETA1;     mp->beta2=DEFAULT_BETA2;
	mp->K=DEFAULT_K;             mp->gamma=DEFAULT_GAMMA;
	mp->D_water=DEFAULT_D_WATER; mp->lambda=DEFAULT_LAMBDA;
	mp->chi1=DEFAULT_CHI1;       mp->chi2=DEFAULT_CHI2;
	mp->Cp=DEFAULT_CP;           mp->DS=DEFAULT_DS;
}
//---------------------------------------------------------------------------
PetscScalar MPgetFReactive(PetscScalar P,PetscScalar T,PetscScalar Cf,PetscScalar M,meltPar_Katz *mp)
{
	// Returns the equilibrium wt fraction degree of melting given pressure
	// (GPa), temperature (Centrigrade), bulk water (weight fraction) and
	// modal cpx (weight fraction).

	PetscScalar Fcpx = M/(mp->r1 + mp->r2*P);
	PetscScalar Tsol, Tlhz, Tcpx, Tliq, dT, Cf_SAT, F, T0, dTdP;

	Cf_SAT = mp->chi1*pow(P,mp->lambda) + mp->chi2*P;
	Cf = (Cf<Cf_SAT) ? Cf:Cf_SAT;
	dT = mp->K*pow(100.0*Cf,mp->gamma);
	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;
	if(P>Pc){
		T0   = mp->A1 + mp->A2*Pc + mp->A3*Pc*Pc;
		dTdP = mp->A2 + 2.0*mp->A3*Pc;
		Tsol = T0 + dTdP*(P-Pc);
	}
	Tlhz = mp->B1 + mp->B2*P + mp->B3*P*P;
	if(P>Pc){
		T0   = mp->B1 + mp->B2*Pc + mp->B3*Pc*Pc;
		dTdP = mp->B2 + 2.0*mp->B3*Pc;
		Tlhz = T0 + dTdP*(P-Pc);
	}
	Tcpx = pow(Fcpx,1.0/mp->beta1)*(Tlhz - Tsol) + Tsol;
	Tliq = mp->C1 + mp->C2*P + mp->C3*P*P;
	if(P>Pc){
		T0   = mp->C1 + mp->C2*Pc + mp->C3*Pc*Pc;
		dTdP = mp->C2 + 2.0*mp->C3*Pc;
		Tliq = T0 + dTdP*(P-Pc);
	}

	if (T<Tsol-dT) {
		F = 0.0;
	} else if (T<Tcpx-dT) {
		F = pow(((T-(Tsol-dT))/(Tlhz-Tsol)),mp->beta1);
	} else if (T<Tliq-dT) {
		F = Fcpx + (1-Fcpx)*pow(((T-(Tcpx-dT))/(Tliq-Tcpx)),mp->beta2);
	} else {
		F = 1.0;
	}

	return F;
}
//---------------------------------------------------------------------------
PetscScalar MPgetFEquilib(PetscScalar P,PetscScalar T,PetscScalar X,PetscScalar M,meltPar_Katz *mp)
{
	// Returns the equilibrium wt fraction degree of melting given pressure
	// (GPa), temperature (Centrigrade), bulk water (weight fraction) and
	// modal cpx (weight fraction).

	PetscScalar Fcpx = M/(mp->r1 + mp->r2*P);
	PetscScalar Tsol, Tlhz, Tcpx, Tliq, dT[3], T0, dTdP;

	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;
	if(P>Pc){
		T0   = mp->A1 + mp->A2*Pc + mp->A3*Pc*Pc;
		dTdP = mp->A2 + 2.0*mp->A3*Pc;
		Tsol = T0 + dTdP*(P-Pc);
	}
	Tlhz = mp->B1 + mp->B2*P + mp->B3*P*P;
	if(P>Pc){
		T0   = mp->B1 + mp->B2*Pc + mp->B3*Pc*Pc;
		dTdP = mp->B2 + 2.0*mp->B3*Pc;
		Tlhz = T0 + dTdP*(P-Pc);
	}
	Tcpx = pow(Fcpx,1.0/mp->beta1)*(Tlhz - Tsol) + Tsol;
	Tliq = mp->C1 + mp->C2*P + mp->C3*P*P;
	if(P>Pc){
		T0   = mp->C1 + mp->C2*Pc + mp->C3*Pc*Pc;
		dTdP = mp->C2 + 2.0*mp->C3*Pc;
		Tliq = T0 + dTdP*(P-Pc);
	}

	dT[0] = calcDT(P,X,0.0,mp);  /* compute dT for F=0.0 */
	dT[1] = calcDT(P,X,Fcpx,mp); /* compute dT for F=Fcpx_out */
	dT[2] = calcDT(P,X,1.0,mp);  /* compute dT for F=1.0 */

	if (T<=(Tsol-dT[0])) {
		return 0.0;
	} else if (T<=(Tcpx-dT[1])) {
		return FX_bal(0.0,Fcpx,T,P,X,Fcpx,mp);
	} else if (T<=(Tliq-dT[2])) {
		return FX_bal(Fcpx,1.0,T,P,X,Fcpx,mp);
	} else
		return 1.0;
}
//---------------------------------------------------------------------------
PetscScalar MPgetTSolidus(PetscScalar P,PetscScalar X,meltPar_Katz *mp)
{
	// Returns the solidus temperature (Centigrade) given pressure, bulk water content.

	PetscScalar Tsol, dT, T0, dTdP;

	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;
	if(P>Pc){
		T0   = mp->A1 + mp->A2*Pc + mp->A3*Pc*Pc;
		dTdP = mp->A2 + 2.0*mp->A3*Pc;
		Tsol = T0 + dTdP*(P-Pc);
	}
	dT   = calcDT(P,X,0.0,mp);  /* compute dT for F=0.0 */
	return Tsol - dT;
}
//---------------------------------------------------------------------------
PetscScalar MPgetTEquilib(PetscScalar P,PetscScalar F,PetscScalar X,PetscScalar M,meltPar_Katz *mp)
{
	// Returns the temperature (Centigrade) given pressure, degree of melting,
	// bulk water content and modal cpx.

	PetscScalar Fcpx = M/(mp->r1 + mp->r2*P);
	PetscScalar Tsol, Tlhz, Tcpx, Tliq, T0, dTdP;

	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;
	if(P>Pc){
		T0   = mp->A1 + mp->A2*Pc + mp->A3*Pc*Pc;
		dTdP = mp->A2 + 2.0*mp->A3*Pc;
		Tsol = T0 + dTdP*(P-Pc);
	}
	Tlhz = mp->B1 + mp->B2*P + mp->B3*P*P;
	if(P>Pc){
		T0   = mp->B1 + mp->B2*Pc + mp->B3*Pc*Pc;
		dTdP = mp->B2 + 2.0*mp->B3*Pc;
		Tlhz = T0 + dTdP*(P-Pc);
	}
	Tcpx = pow(Fcpx,1.0/mp->beta1)*(Tlhz - Tsol) + Tsol;
	Tliq = mp->C1 + mp->C2*P + mp->C3*P*P;
	if(P>Pc){
		T0   = mp->C1 + mp->C2*Pc + mp->C3*Pc*Pc;
		dTdP = mp->C2 + 2.0*mp->C3*Pc;
		Tliq = T0 + dTdP*(P-Pc);
	}


	if (F<=0.0) {
		return Tsol - calcDT(P,X,0.0,mp);
	} else if (F<=Fcpx) {
		return pow(F,1.0/mp->beta1)*(Tlhz-Tsol) + Tsol - calcDT(P,X,F,mp);
	} else if (F<1.0) {
		return pow((F-Fcpx)/(1.0-Fcpx),1.0/mp->beta2)*(Tliq-Tcpx) + Tcpx - calcDT(P,X,F,mp);
	} else
		return Tliq - calcDT(P,X,1.0,mp);
}  
//---------------------------------------------------------------------------
PetscScalar MPgetFconsH(PetscScalar P,PetscScalar Ti,PetscScalar X,PetscScalar M,PetscScalar *Tf,meltPar_Katz *mp)
{
	// The function whose root is the melt fraction that conserves enthalpy.

	PetscScalar Tsol,F;

	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;

	if (Ti<(Tsol-calcDT(P,X,0.0,mp))) {
		*Tf = Ti;
		return 0.0;
	} else {
		F  = FT_bal(0.0,1.0,Ti,P,X,M,mp);
		*Tf = MPgetTEquilib(P,F,X,M,mp);
		return F;
	}
}
//---------------------------------------------------------------------------
PetscScalar calcDT(PetscScalar P,PetscScalar X,PetscScalar F,meltPar_Katz *mp)
{
	// Compute the solidus lowering effect of bulk water content GIVEN a degree of melting.

	PetscScalar Xsat, Xmelt;

	Xsat  = mp->chi1*pow(P,mp->lambda) + mp->chi2*P; /* saturation content */
	Xmelt = X/( mp->D_water + F*(1.0-mp->D_water) ); /* melt water content */
	Xmelt = (Xmelt<Xsat) ? Xmelt:Xsat;               /* min of the two */
	return mp->K*pow(100.0*Xmelt,mp->gamma);         /* delta T due to water */
}
//---------------------------------------------------------------------------
PetscScalar calcF(PetscScalar T,PetscScalar dT,PetscScalar P,PetscScalar Fcpx,meltPar_Katz *mp)
{
	// The equations that give F from T, dT, P, Fcpx.

	PetscScalar Tsol, Tlhz, Tcpx, Tliq;

	Tsol = mp->A1 + mp->A2*P + mp->A3*P*P;
	Tlhz = mp->B1 + mp->B2*P + mp->B3*P*P;
	Tcpx = pow(Fcpx,1.0/mp->beta1)*(Tlhz - Tsol) + Tsol;
	Tliq = mp->C1 + mp->C2*P + mp->C3*P*P;

	if (T<=(Tsol-dT)) {
		return 0.0;
	} else if (T<=(Tcpx-dT)) {
		return pow((T - (Tsol-dT))/(Tlhz-Tsol),mp->beta1);
	} else if (T<=(Tliq-dT)) {
		return Fcpx + (1.0-Fcpx)*pow((T - (Tcpx-dT))/(Tliq-Tcpx),mp->beta2);
	} else
		return 1.0;
}
//---------------------------------------------------------------------------
PetscScalar FZero(PetscScalar F, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar Fcpx,meltPar_Katz *mp)
{
	// The function whose zero is the correct melt fraction for a given set of state variables.
	return calcF(T,calcDT(P,X,F,mp),P,Fcpx,mp) - F;
}
//---------------------------------------------------------------------------
PetscScalar HZero(PetscScalar F,PetscScalar T,PetscScalar P,PetscScalar X,PetscScalar M,meltPar_Katz *mp)
{
	// The function whose root is the melt fraction that conserves enthalpy.
	return (MPgetTEquilib(P,F,X,M,mp)+273.0)*(mp->Cp+mp->DS*F) - (T+273.0)*mp->Cp;
}
//---------------------------------------------------------------------------
PetscScalar FX_bal(PetscScalar x1,PetscScalar x2,PetscScalar T,PetscScalar P,PetscScalar X,PetscScalar Fcpx,meltPar_Katz *mp)
{
	// Adapted from Ridder's Method as described by Numerical Recipes in C.

	PetscInt    j;
	PetscScalar fh,fl,fm,fnew,s,xh,xl,xm,xnew,ans;

	fl = FZero(x1,T,P,X,Fcpx,mp);
	fh = FZero(x2,T,P,X,Fcpx,mp);

	if ( (fl>0.0 && fh<0.0) || (fl<0.0 && fh>0.0) ) {
		xl=x1; xh=x2; ans=UNUSED;

		for (j=1;j<=MAXITS;j++) {
			xm=0.5*(xl+xh);
			fm=FZero(xm,T,P,X,Fcpx,mp);
			s=sqrt(fm*fm-fl*fh);
			if (s==0.0) return ans;
			xnew=xm+(xm-xl)*((fl>=fh ? 1.0 : -1.0)*fm/s);
			if(fabs(xnew-ans) <= X_ACC) return ans;
			ans=xnew;
			fnew=FZero(ans,T,P,X,Fcpx,mp);
			if (fnew==0.0) return ans;
			if (SIG(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIG(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIG(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else
				PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: never get here (1)\n");

			if (fabs(xh-xl) <= X_ACC) return ans;
		}
		PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: exceed max iterations\n");

	} else {
		if (fl==0.0) return x1;
		if (fh==0.0) return x2;
		PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: never get here (2)\n");
	}
	return 0.0; /* never used */
}
//---------------------------------------------------------------------------
PetscScalar FT_bal(PetscScalar x1,PetscScalar x2,PetscScalar T,PetscScalar P,PetscScalar X,PetscScalar M,meltPar_Katz *mp)
{
	// Adapted from Ridder's Method as described by Numerical Recipes in C.

	PetscInt    j;
	PetscScalar fh,fl,fm,fnew,s,xh,xl,xm,xnew,ans;

	fl = HZero(x1,T,P,X,M,mp);
	fh = HZero(x2,T,P,X,M,mp);

	if ( (fl>0.0 && fh<0.0) || (fl<0.0 && fh>0.0) ) {
		xl=x1; xh=x2; ans=UNUSED;

		for (j=1;j<=MAXITS;j++) {
			xm=0.5*(xl+xh);
			fm=HZero(xm,T,P,X,M,mp);
			s=sqrt(fm*fm-fl*fh);
			if (s==0.0) return ans;
			xnew=xm+(xm-xl)*((fl>=fh ? 1.0 : -1.0)*fm/s);
			if(fabs(xnew-ans) <= X_ACC) return ans;
			ans=xnew;
			fnew=HZero(ans,T,P,X,M,mp);
			if (fnew==0.0) return ans;
			if (SIG(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIG(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIG(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else
				PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: never get here (1)\n");

			if (fabs(xh-xl) <= X_ACC) return ans;
		}
		PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: exceed max iterations\n");

	} else {
		if (fl==0.0) return x1;
		if (fh==0.0) return x2;
		PetscPrintf(PETSC_COMM_WORLD,"FX_bal error: never get here (2)\n");
	}
	return 0.0; /* never used */
}
//---------------------------------------------------------------------------
