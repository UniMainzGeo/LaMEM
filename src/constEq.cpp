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
 **    filename:   constEq.c
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
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "constEq.h"
#include "phase.h"
#include "JacRes.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ConstEqCtxSetup"
PetscErrorCode ConstEqCtxSetup(
		ConstEqCtx  *ctx,      // evaluation context
		Material_t  *mat,      // phase parameters
		Soft_t      *soft,     // material softening laws
		Controls    *ctrl,     // parameters and controls
		PetscScalar  DII,      // effective strain-rate
		PetscScalar  APS,      // accumulated plastic strain
		PetscScalar  dt,       // time step
		PetscScalar  p,        // pressure
		PetscScalar  p_lithos, // lithostatic pressure
		PetscScalar  p_pore,   // pore pressure
		PetscScalar  T)        // temperature
{
	// setup nonlinear constitutive equation evaluation context
	// evaluate dependence on constant parameters (pressure, temperature)

	PetscInt    pd;
	PetscScalar Q, RT, ch, fr, p_visc, p_upper, p_lower, dP, p_total;

	PetscFunctionBegin;

	// set RT
	RT         = ctrl->Rugc*T;
	if(!RT) RT = -1.0;

	// use reference strain-rate instead of zero
	if(DII == 0.0) DII = ctrl->DII_ref;

	// initialize
	ctx->DII   = DII; // effective strain-rate
	ctx->A_els = 0.0; // elasticity constant
	ctx->A_dif = 0.0; // diffusion constant
	ctx->A_dis = 0.0; // dislocation constant
	ctx->N_dis = 1.0; // dislocation exponent
	ctx->A_prl = 0.0; // Peierls constant
	ctx->N_prl = 1.0; // Peierls exponent
	ctx->taupl = 0.0; // plastic yield stress
	ctx->fr    = 0.0; // effective friction coefficient
	pd         = 0;   // pressure-dependence flag

	//===============
	// TOTAL PRESSURE
	//===============

	if(ctrl->gwType == _GW_NONE_) p_pore = 0.0;

	p_total = p + ctrl->biot*p_pore;

	//=================
	// VISCO-ELASTICITY
	//=================

	// assign pressure for viscous laws
	if(ctrl->pLithoVisc)  p_visc = p_lithos;
	else                  p_visc = p_total;

	// ELASTICITY
	if(mat->G)
	{
		// Elasticity correction can only DECREASE the viscosity.
		// eta/G << dt (viscous regime)  eta*(dt/(dt + eta/G)) -> eta
		// eta/G >> dt (elastic regime)  eta*(dt/(dt + eta/G)) -> G*dt < eta
		// Elasticity doesn't normally interact with the bottom viscosity limit,
		// instead it rather acts as a smooth limiter for maximum viscosity.

		ctx->A_els = 0.5/(mat->G*dt);
	}

	// LINEAR DIFFUSION CREEP (NEWTONIAN)
	if(mat->Bd)
	{
		Q          = (mat->Ed + p_visc*mat->Vd)/RT;
		ctx->A_dif =  mat->Bd*exp(-Q);
	}

	// DISLOCATION CREEP (POWER LAW)
	if(mat->Bn)
	{
		Q          = (mat->En + p_visc*mat->Vn)/RT;
		ctx->N_dis =  mat->n;
		ctx->A_dis =  mat->Bn*exp(-Q);
	}

	// PEIERLS CREEP (LOW TEMPERATURE RATE-DEPENDENT PLASTICITY, POWER-LAW APPROXIMATION)
	// ONLY EVALUATE FOR TEMPERATURE-DEPENDENT CASES
	if(mat->Bp && T)
	{
		Q          = (mat->Ep + p_visc*mat->Vp)/RT;
		ctx->N_prl =  Q*pow(1.0-mat->gamma, mat->q-1.0)*mat->q*mat->gamma;
		ctx->A_prl =  mat->Bp/pow(mat->gamma*mat->taup, ctx->N_prl)*exp(-Q*pow(1.0-mat->gamma, mat->q));
	}

	//===========
	// PLASTICITY
	//===========

	// apply strain softening to friction and cohesion
	ch = ApplyStrainSoft(soft, mat->chSoftID, APS, mat->ch);
	fr = ApplyStrainSoft(soft, mat->frSoftID, APS, mat->fr);

	// fit to limits
	if(ch < ctrl->minCh) ch = ctrl->minCh;
	if(fr < ctrl->minFr) fr = ctrl->minFr;

	// override total pressure with lithostatic if requested
	if(ctrl->pLithoPlast)
	{
		// Use lithostatic, rather than dynamic pressure to evaluate yield stress
		// This converges better, but does not result in localization of deformation & shear banding,
		// so only apply it for large-scale simulations where plasticity does not matter much

		p_total = p_lithos;
	}
	else if(ctrl->pLimPlast)
	{
		// apply pressure limits

		// yielding surface: (S1-S3)/2 = (S1+S3)/2*sin(phi) + C*cos(phi)
		// pressure can be written as: P = (S1+S2+S3)/3 and P~=S2,then P=(S1+S3)/2
		// so the yield surface can be rewritten as:
		// P-S3=P*sin(phi) + C*cos(phi)   --> compression
		// S1-P=P*sin(phi) + C*cos(phi)   --> extension
		// under pure shear compression, S3=P_Lithos and S1=P_Lithos when extension
		// P = -( S3+C*cos(phi))/(sin(phi)-1)  --> compression
		// P = -(-S1+C*cos(phi))/(sin(phi)+1)  --> extension

		p_upper = -( p_lithos + ch * cos(fr))/(sin(fr) - 1.0); // compression
		p_lower = -(-p_lithos + ch * cos(fr))/(sin(fr) + 1.0); // extension

		if(p_total > p_upper) p_total = p_upper;
		if(p_total < p_lower) p_total = p_lower;
	}

	// compute cohesion and friction coefficient
	ch = cos(fr)*ch;
	fr = sin(fr);

	// compute effective mean stress
	dP = (p_total - p_pore);

	// compute yield stress
	if(dP < 0.0) { ctx->taupl =         ch; pd = 0; } // Von-Mises model for extension
	else         { ctx->taupl = dP*fr + ch; pd = 1; } // Drucker-Prager model for compression

	// correct for ultimate yield stress (if defined)
	if(ctrl->tauUlt) { if(ctx->taupl > ctrl->tauUlt) { ctx->taupl = ctrl->tauUlt; pd = 0; } }

	// store friction coefficient for a pressure-dependent plasticity model
	if(pd) ctx->fr = fr;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar GetConsEqRes(PetscScalar eta, void *pctx)
{
	// Compute residual of the nonlinear visco-elastic constitutive equation:

	PetscScalar tauII, DII_els, DII_dif, DII_dis, DII_prl;

	// access context
	ConstEqCtx *ctx = (ConstEqCtx*)pctx;

	// compute stress
	tauII = 2.0*eta*ctx->DII;

	// creep strain rates
	DII_els = ctx->A_els*tauII;                  // elasticity
	DII_dif = ctx->A_dif*tauII;                  // diffusion
	DII_dis = ctx->A_dis*pow(tauII, ctx->N_dis); // dislocation
	DII_prl = ctx->A_prl*pow(tauII, ctx->N_prl); // Peierls

	// residual function (r)
	// r < 0 if eta > solution (negative on overshoot)
	// r > 0 if eta < solution (positive on undershoot)

	return ctx->DII - (DII_els + DII_dif + DII_dis + DII_prl);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetEffVisc"
PetscErrorCode GetEffVisc(
		ConstEqCtx  *ctx,
		Controls    *ctrl,        // parameters and controls
		PetscScalar *eta_total,
		PetscScalar *eta_creep,
		PetscScalar *eta_vp,
		PetscScalar *DIIpl,
		PetscScalar *dEta,
		PetscScalar *fr)
{
	// stabilization parameters
	PetscScalar eta_ve, eta_pl, eta_pw, eta_vp_reg, eta_st, eta_dis, eta_prl, cf;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_dis, inv_eta_prl, inv_eta_max;

	PetscFunctionBegin;

	// initialize
	(*DIIpl) = 0.0;
	(*dEta)  = 0.0;
	(*fr)    = 0.0;

	inv_eta_max = ctrl->inv_eta_max;

	//==============
	// INITIAL GUESS
	//==============

	// set reference viscosity as initial guess
	if(ctrl->initGuess)
	{
		(*eta_total) = ctrl->eta_ref;
		(*eta_creep) = ctrl->eta_ref;
		(*eta_vp)    = ctrl->eta_ref;
		
		PetscFunctionReturn(0);
	}

	//=====================
	// ISOLATED VISCOSITIES
	//=====================

	inv_eta_els = 0.0;
	inv_eta_dif = 0.0;
	inv_eta_dis = 0.0;
	inv_eta_prl = 0.0;

	// elasticity
	if(ctx->A_els) inv_eta_els = 2.0*ctx->A_els;
	// diffusion
	if(ctx->A_dif) inv_eta_dif = 2.0*ctx->A_dif;
	// dislocation
	if(ctx->A_dis) inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(ctx->DII, 1.0 - 1.0/ctx->N_dis);
	// Peierls
	if(ctx->A_prl) inv_eta_prl = 2.0*pow(ctx->A_prl, 1.0/ctx->N_prl)*pow(ctx->DII, 1.0 - 1.0/ctx->N_prl);

	// error handling
	if(PetscIsInfOrNanScalar(inv_eta_dif)) inv_eta_dif = 0.0;
	if(PetscIsInfOrNanScalar(inv_eta_dis)) inv_eta_dis = 0.0;
	if(PetscIsInfOrNanScalar(inv_eta_prl)) inv_eta_prl = 0.0;

	//================
	// CREEP VISCOSITY
	//================

	// store creep & viscoplastic viscosity for output
	(*eta_creep) = ctrl->eta_min + 1.0/(inv_eta_dif + inv_eta_dis + inv_eta_prl + inv_eta_max);
	(*eta_vp)    = (*eta_creep);

	//========================
	// VISCO-ELASTIC VISCOSITY
	//========================

	// compute visco-elastic viscosity
	eta_ve = ctrl->eta_min + 1.0/(inv_eta_els + inv_eta_dif + inv_eta_dis + inv_eta_prl + inv_eta_max);

	// set visco-elastic prediction
	(*eta_total) = eta_ve;

	// get viscosity derivatives
	if(inv_eta_dis)
	{
		eta_dis  = 1.0/inv_eta_dis;
		cf       = (eta_ve - ctrl->eta_min)/eta_dis;
		(*dEta) += cf*cf*(1.0/ctx->N_dis - 1.0)*eta_dis;
	}
	if(inv_eta_prl)
	{
		eta_prl  = 1.0/inv_eta_prl;
		cf       = (eta_ve - ctrl->eta_min)/eta_prl;
		(*dEta) += cf*cf*(1.0/ctx->N_prl - 1.0)*eta_prl;
	}

	//===========
	// PLASTICITY
	//===========

	if(ctx->taupl && !ctrl->initGuess)
	{
		// compute true plastic viscosity
		eta_pl = ctx->taupl/(2.0*ctx->DII);

		// ultimate viscosity cutoff for plasticity
		if(eta_pl < ctrl->eta_min) eta_pl = ctrl->eta_min;

		//==============================================
		// compute total viscosity
		// minimum viscosity (true) model is the default
		//==============================================

		if(ctrl->quasiHarmAvg)
		{
			// quasi-harmonic mean
			(*eta_total) = 1.0/(1.0/eta_pl + 1.0/eta_ve);
		}
		else if(ctrl->n_pw)
		{
			// rate-dependent power-law stabilization
			eta_pw = (ctx->taupl/2.0)*pow(ctx->DII, 1/ctrl->n_pw - 1.0);

			if(eta_pw < eta_ve) (*eta_total) = eta_pw;
		}
		else if(ctrl->cf_eta_min)
		{
			// rate-dependent visco-plastic regularization
			eta_st = eta_ve/ctrl->cf_eta_min;

			eta_vp_reg = (eta_st + eta_pl)/(1.0 + eta_st/eta_ve);

			if(eta_vp_reg < eta_ve) (*eta_total) = eta_vp_reg;
		}
		else if(eta_pl < eta_ve)
		{
			// minimum viscosity model (unstable) (default)
			(*eta_total) = eta_pl;
		}

		if(eta_pl < eta_ve)
		{
			// store plastic strain rate, viscosity derivative & effective friction
			(*eta_vp) =  eta_pl;
			(*DIIpl)  =  ctx->DII*(1.0 - (*eta_total)/eta_ve);
			(*dEta)   = -eta_pl;
			(*fr)     =  ctx->fr;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar ApplyStrainSoft(Soft_t *soft, PetscInt ID, PetscScalar APS, PetscScalar par)
{
	// Apply strain softening to a parameter (friction, cohesion)
	PetscScalar  k;
	Soft_t      *s;
	// check whether softening is defined
	if(ID == -1) return par;
	// access parameters
	s = soft + ID;
	// compute scaling ratio
	if(APS <= s->APS1)
		k = 1.0;
	if(APS > s->APS1 && APS < s->APS2)
		k = 1.0 - s->A*((APS - s->APS1)/(s->APS2 - s->APS1));
	if(APS >= s->APS2)
		k = 1.0 - s->A;
	// apply strain softening
	return par*k;
}
//---------------------------------------------------------------------------
// compute inverse deviatoric elastic parameter
PetscScalar GetI2Gdt(
	PetscInt     numPhases,
	Material_t  *phases,
	PetscScalar *phRat,
	PetscScalar  dt)
{
	PetscInt    i;
	PetscScalar I2Gdt, Gavg;

	Gavg  = 0.0;
	I2Gdt = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// average elastic materials only
		if(phases[i].G)
		{
			Gavg += phRat[i]*phases[i].G;
		}
	}

	if(Gavg) I2Gdt = 1.0/Gavg/dt/2.0;

	return I2Gdt;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DevConstEq"
PetscErrorCode DevConstEq(
		SolVarDev   *svDev,     // solution variables
		PetscScalar *eta_creep, // creep viscosity (for output)
		PetscScalar *eta_vp,    // viscoplastic viscosity (for output)
		PetscInt     numPhases, // number phases
		Material_t  *phases,    // phase parameters
		Soft_t      *soft,      // material softening laws
		PetscScalar *phRat,     // phase ratios
		Controls    *ctrl,       // parameters and controls
		PetscScalar  p_lithos,  // lithostatic pressure
		PetscScalar  p_pore,    // pore pressure
		PetscScalar  dt,        // time step
		PetscScalar  p,         // pressure
		PetscScalar  T)         // temperature
{
	// Evaluate deviatoric constitutive equations in control volume

	PetscInt     i;
	ConstEqCtx   ctx;
	Material_t  *mat;
	PetscScalar  DII, APS, eta_total, eta_creep_phase, eta_viscoplastic_phase, DIIpl, dEta, fr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize total effective strain rate & APS
	DII = svDev->DII;
	APS = svDev->APS;

	// initialize effective viscosity & plastic strain-rate
	svDev->eta   = 0.0;
	svDev->DIIpl = 0.0;
	(*eta_creep) = 0.0;
	(*eta_vp)    = 0.0;

	svDev->dEta  = 0.0;
	svDev->fr    = 0.0;
	svDev->yield = 0.0;
	dEta         = 0.0;
	fr           = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// get reference to material parameters table
			mat = &phases[i];

			// setup nonlinear constitutive equation evaluation context
			ierr = ConstEqCtxSetup(&ctx, mat, soft, ctrl, DII, APS, dt, p, p_lithos, p_pore, T); CHKERRQ(ierr);

			// solve effective viscosity & plastic strain rate
			ierr = GetEffVisc(&ctx, ctrl, &eta_total, &eta_creep_phase, &eta_viscoplastic_phase, &DIIpl, &dEta, &fr); CHKERRQ(ierr);

			// average parameters
			svDev->eta   += phRat[i]*eta_total;
			svDev->DIIpl += phRat[i]*DIIpl;
			(*eta_creep) += phRat[i]*eta_creep_phase;
			(*eta_vp)    += phRat[i]*eta_viscoplastic_phase;

			svDev->dEta  += phRat[i]*dEta;
			svDev->fr    += phRat[i]*fr;
			svDev->yield += phRat[i]*ctx.taupl;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VolConstEq"
PetscErrorCode VolConstEq(
		SolVarBulk  *svBulk,    // solution variables
		PetscInt     numPhases, // number phases
		Material_t  *phases,    // phase parameters
		PetscScalar *phRat,     // phase ratios
		Controls    *ctrl,      // parameters and controls
		PetscScalar  depth,     // depth for depth-dependent density model
		PetscScalar  dt,        // time step
		PetscScalar  p,         // pressure
		PetscScalar  T)         // temperature
{
	// Evaluate volumetric constitutive equations in control volume
	PetscInt     i;
	Material_t  *mat;
	PetscScalar  cf_comp, cf_therm, Kavg, rho;

	PetscFunctionBegin;

	// initialize effective density, thermal expansion & inverse bulk elastic parameter
	svBulk->rho   = 0.0;
	svBulk->alpha = 0.0;
	svBulk->IKdt  = 0.0;
	Kavg          = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// get reference to material parameters table
			mat = &phases[i];

			// initialize
			cf_comp  = 1.0;
			cf_therm = 1.0;

			// elastic compressiblility correction (Murnaghan's equation)
			// ro/ro_0 = (1 + K'*P/K)^(1/K')
			if(mat->K)
			{
				Kavg += phRat[i]*mat->K;

				if(mat->Kp) cf_comp = pow(1.0 + mat->Kp*(p/mat->K), 1.0/mat->Kp);
				else        cf_comp = 1.0 + p/mat->K;
			}

			// ro/ro_0 = (1 + beta*P)
			if(mat->beta)
			{
				// negative sign as compressive pressures (increasing depth) is negative in LaMEM
				cf_comp = 1.0 + p*mat->beta;
			}

			// thermal expansion correction
			// ro/ro_0 = 1 - alpha*(T - TRef)
			if(mat->alpha)
			{
				cf_therm  = 1.0 - mat->alpha*(T - ctrl->TRef);
			}

			// get density
			if(mat->rho_n)
			{
				// depth-dependent density (ad-hoc)
				rho = mat->rho - (mat->rho - ctrl->rho_fluid)*mat->rho_n*exp(-mat->rho_c*depth);
			}
			else
			{
				// temperature & pressure-dependent density
				rho = mat->rho*cf_comp*cf_therm;
			}

			// update density, thermal expansion & inverse bulk elastic parameter
			svBulk->rho   += phRat[i]*rho;
			svBulk->alpha += phRat[i]*mat->alpha;
		}
	}

	if(Kavg) svBulk->IKdt = 1.0/Kavg/dt;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetStressCell"
PetscErrorCode GetStressCell(
		SolVarCell  *svCell, // solution variables
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    //
		PetscScalar  dzz)   //
{
	// compute stress, plastic strain-rate and shear heating term on cell

	SolVarDev   *svDev;
	PetscScalar  DII, cfpl, txx, tyy, tzz;

	PetscFunctionBegin;

	// access deviatoric variables
	svDev = &svCell->svDev;

	// compute deviatoric stresses
	svCell->sxx = 2.0*svDev->eta*dxx;
	svCell->syy = 2.0*svDev->eta*dyy;
	svCell->szz = 2.0*svDev->eta*dzz;

	// get strain-rate invariant
	DII = svDev->DII;

	// compute plastic scaling coefficient
	if(DII) cfpl = svDev->DIIpl/DII;
	else    cfpl = 0.0;

	// compute plastic strain-rate components
	txx = cfpl*dxx;
	tyy = cfpl*dyy;
	tzz = cfpl*dzz;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = 0.5*(txx*txx + tyy*tyy + tzz*tzz);

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	txx = svCell->dxx - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tyy = svCell->dyy - svDev->I2Gdt*(svCell->syy - svCell->hyy);
	tzz = svCell->dzz - svDev->I2Gdt*(svCell->szz - svCell->hzz);

	// compute shear heating term contribution
	svDev->Hr = (txx*svCell->sxx + tyy*svCell->syy + tzz*svCell->szz);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// compute stress, plastic strain-rate and shear heating term on edge
#undef __FUNCT__
#define __FUNCT__ "GetStressEdge"
PetscErrorCode GetStressEdge(
		SolVarEdge  *svEdge, // solution variables
		PetscScalar  d)      // effective shear strain rate component
{

	SolVarDev   *svDev;
	PetscScalar  DII, cfpl, t;

	PetscFunctionBegin;

	// access deviatoric variables
	svDev = &svEdge->svDev;

	// compute shear stress
	svEdge->s = 2.0*svDev->eta*d;

	// get strain-rate invariant
	DII = svDev->DII;

	// compute plastic scaling coefficient
	if(DII) cfpl = svDev->DIIpl/DII;
	else    cfpl = 0.0;

	// compute plastic strain-rate components
	t = cfpl*d;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = t*t;

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	t = svEdge->d - svDev->I2Gdt*(svEdge->s - svEdge->h);

	// compute shear heating term contribution
	svDev->Hr = 2.0*t*svEdge->s;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
