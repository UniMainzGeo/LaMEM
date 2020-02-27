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
#include "tools.h"
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "volConstEq"
PetscErrorCode volConstEq(SolVarBulk *svBulk, ConstEqCtx *ctx)
{
	// evaluate volumetric constitutive equations in control volume

	Controls    *ctrl;
	PData       *pd;
	Material_t  *mat, *phases;
	PetscInt     i, numPhases;
	PetscScalar *phRat, dt, p, depth, T, cf_comp, cf_therm, Kavg, rho;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ctrl      = ctx->ctrl;
	numPhases = ctx->numPhases;
	phases    = ctx->phases;
	phRat     = ctx->phRat;
	pd        = ctx->pd;
	depth     = ctx->depth;
	dt        = ctx->dt;
	p         = ctx->p;
	T         = ctx->T;

	// initialize effective density, thermal expansion & inverse bulk elastic parameter
	svBulk->rho    = 0.0;
	svBulk->alpha  = 0.0;
	svBulk->IKdt   = 0.0;
	Kavg           = 0.0;
	svBulk->mf     = 0.0;
	svBulk->rho_pf = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// get reference to material parameters table
			mat = &phases[i];

			if(mat->pdAct == 1)
			{
				// compute melt fraction from phase diagram
				ierr = setDataPhaseDiagram(pd, p, T, mat->pdn); CHKERRQ(ierr);

				svBulk->mf     += phRat[i]*pd->mf;
				svBulk->rho_pf += phRat[i]*pd->rho_f;
			}

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
			// phase diagram
			else if(mat->pdAct == 1)
			{
				rho = (pd->mf * pd->rho_f) + ((1-pd->mf) * pd->rho);
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
#define __FUNCT__ "devConstEq"
PetscErrorCode devConstEq(SolVarDev *svDev, ConstEqCtx *ctx)
{
	// evaluate deviatoric constitutive equations in control volume

	Controls    *ctrl;
	PetscScalar *phRat;
	PetscInt     i, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ctrl      = ctx->ctrl;
	numPhases = ctx->numPhases;
	phRat     = ctx->phRat;

	// zero out results
	ctx->eta    = 0.0; // effective viscosity
	ctx->eta_cr = 0.0; // creep viscosity
	ctx->eta_vp = 0.0; // visco-plastic viscosity
	ctx->DIIdif = 0.0; // diffusion creep strain rate
	ctx->DIIdis = 0.0; // dislocation creep strain rate
	ctx->DIIprl = 0.0; // Peierls creep strain rate
	ctx->DIIpl  = 0.0; // plastic strain rate
	ctx->yield  = 0.0; // yield stress

	// viscous initial guess
	if(ctrl->initGuess)
	{
		svDev->eta    = ctrl->eta_ref;
		ctx  ->eta    = ctrl->eta_ref;
		ctx  ->eta_cr = ctrl->eta_ref;
		ctx  ->eta_vp = ctrl->eta_ref;
		ctx  ->DIIdif = 1.0;

		PetscFunctionReturn(0);
	}

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// setup phase parameters
			ierr = setUpPhase(ctx, i); CHKERRQ(ierr);

			// compute phase viscosities and strain rate partitioning
			ierr = getPhaseVisc(ctx, i); CHKERRQ(ierr);

		}
	}

	// normalize strain rates
	if(ctx->DII)
	{
		ctx->DIIdif /= ctx->DII;
		ctx->DIIdis /= ctx->DII;
		ctx->DIIprl /= ctx->DII;
		ctx->DIIpl  /= ctx->DII;
	}

	// store viscosity
	svDev->eta = ctx->eta + ctrl->eta_min;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getPhaseVisc"
PetscErrorCode getPhaseVisc(ConstEqCtx *ctx, PetscInt ID)
{
	Controls    *ctrl;
	PetscInt    it, fail;
	PetscScalar eta_min, eta_mean, eta, eta_cr, eta_vp, tauII, taupl, DII;
	PetscScalar DIIdif, DIImax, DIIdis, DIIprl, DIIpl, DIIvs, DIIvp, phRat;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_max, inv_eta_dis, inv_eta_prl, inv_eta_min;

	PetscFunctionBegin;

	// access context
	ctrl  = ctx->ctrl;      // global controls
	phRat = ctx->phRat[ID]; // phase ratio
	taupl = ctx->taupl;     // plastic yield stress
	DII   = ctx->DII;       // effective strain rate

	// initialize
	it     = 1;
	fail   = 0;
	DIIpl  = 0.0;
	eta_cr = 0.0;
	eta_vp = 0.0;

	//===========
	// PLASTICITY
	//===========
	if(taupl && DII)
	{
		// compute plastic stress and viscosity
		tauII = taupl;
		eta   = tauII/(2.0*DII);

		// compute plastic strain rate
		DIIpl = getConsEqRes(eta, ctx);

		// reset if plasticity is not active
		if(DIIpl < 0.0)
		{
			DIIpl = 0.0;
		}
	}

	//=================
	// VISCO-ELASTICITY
	//=================
	if(!DIIpl)
	{
		// get isolated viscosities
		inv_eta_els = 0.0;
		inv_eta_dif = 0.0;
		inv_eta_max = 0.0;
		inv_eta_dis = 0.0;
		inv_eta_prl = 0.0;

		// elasticity
		if(ctx->A_els) inv_eta_els = 2.0*ctx->A_els;
		// diffusion
		if(ctx->A_dif) inv_eta_dif = 2.0*ctx->A_dif;
		// upper bound
		if(ctx->A_max) inv_eta_max = 2.0*ctx->A_max;
		// dislocation
		if(ctx->A_dis) inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(DII, 1.0 - 1.0/ctx->N_dis);
		// Peierls
		if(ctx->A_prl) inv_eta_prl = 2.0*pow(ctx->A_prl, 1.0/ctx->N_prl)*pow(DII, 1.0 - 1.0/ctx->N_prl);

		// get minimum viscosity (upper bound)
		inv_eta_min                               = inv_eta_els;
		if(inv_eta_dif > inv_eta_min) inv_eta_min = inv_eta_dif;
		if(inv_eta_max > inv_eta_min) inv_eta_min = inv_eta_max;
		if(inv_eta_dis > inv_eta_min) inv_eta_min = inv_eta_dis;
		if(inv_eta_prl > inv_eta_min) inv_eta_min = inv_eta_prl;
		eta_min = 1.0/inv_eta_min;

		// get quasi-harmonic mean (lower bound)
		eta_mean = 1.0/(inv_eta_els + inv_eta_dif + inv_eta_max + inv_eta_dis + inv_eta_prl);

		// NOTE: if closed-form solution exists, it is equal to lower bound
		// If only one mechanism is active, then both bounds are coincident
		// Bisection will return immediately if closed-form solution exists

		if(DII)
		{
			// apply bisection algorithm to nonlinear scalar equation
			fail = solveBisect(eta_mean, eta_min, ctrl->lrtol, ctrl->lmaxit, eta, it, getConsEqRes, ctx);
		}
		else
		{
			eta = eta_mean;
		}

		// compute stress
		tauII = 2.0*eta*DII;
	}

	// update iteration statistics
	ctx->stats[0] += 1.0;               // starts
	ctx->stats[1] += (PetscScalar)fail; // fails
	ctx->stats[2] += (PetscScalar)it;   // iterations

	// compute strain rates
	DIIdif = ctx->A_dif*tauII;                  // diffusion
	DIImax = ctx->A_max*tauII;                  // upper bound
	DIIdis = ctx->A_dis*pow(tauII, ctx->N_dis); // dislocation
	DIIprl = ctx->A_prl*pow(tauII, ctx->N_prl); // Peierls
	DIIvs  = DIIdif + DIImax + DIIdis + DIIprl; // viscous (total)
	DIIvp  = DIIvs + DIIpl;                     // visco-plastic

	// compute viscosities
	if(DIIvs) eta_cr = tauII/DIIvs/2.0;
	if(DIIvp) eta_vp = tauII/DIIvp/2.0;

	// update results
	ctx->eta    += phRat*eta;    // effective viscosity
	ctx->eta_cr += phRat*eta_cr; // creep viscosity
	ctx->eta_vp += phRat*eta_vp; // visco-plastic viscosity
	ctx->DIIdif += phRat*DIIdif; // diffusion creep strain rate
	ctx->DIIdis += phRat*DIIdis; // dislocation creep strain rate
	ctx->DIIprl += phRat*DIIprl; // Peierls creep strain rate
	ctx->DIIpl  += phRat*DIIpl;  // plastic strain rate
	ctx->yield  += phRat*taupl;  // plastic yield rate


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "setUpPhase"
PetscErrorCode setUpPhase(ConstEqCtx *ctx, PetscInt ID)
{
	// setup nonlinear constitutive equation evaluation context
	// evaluate dependence on constant parameters (pressure, temperature)
	Material_t  *mat;
	Soft_t      *soft;
	Controls    *ctrl;
	PData       *pd;
	PetscScalar  APS, dt, p, p_lith, p_pore, T, mf, mfd, mfn;
	PetscScalar  Q, RT, ch, fr, p_visc, p_upper, p_lower, dP, p_total;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	mat    = ctx->phases + ID;
	soft   = ctx->soft;
	ctrl   = ctx->ctrl;
	pd     = ctx->pd;
	APS    = ctx->APS;
	dt     = ctx->dt;
	p      = ctx->p;
	p_lith = ctx->p_lith;
	p_pore = ctx->p_pore;
	T      = ctx->T;
	mf     = 0.0;

	if(mat->pdAct == 1)
	{
		// compute melt fraction from phase diagram
		ierr = setDataPhaseDiagram(pd, p, T, mat->pdn); CHKERRQ(ierr);

		// store melt fraction
		mf = pd->mf;
	}

	// set RT
	RT         =  ctrl->Rugc*T;
	if(!RT) RT = -1.0;

	// initialize phase parameters
	ctx->A_els = 0.0; // elasticity constant
	ctx->A_dif = 0.0; // diffusion constant
	ctx->A_max = 0.0; // upper bound constant
	ctx->A_dis = 0.0; // dislocation constant
	ctx->N_dis = 1.0; // dislocation exponent
	ctx->A_prl = 0.0; // Peierls constant
	ctx->N_prl = 1.0; // Peierls exponent
	ctx->taupl = 0.0; // plastic yield stress

	// MELT FRACTION
	mfd = 1.0;
	mfn = 1.0;

	if(mf)
	{
		// limit melt fraction
		if(mf > ctrl->mfmax) mf = ctrl->mfmax;

		// compute corrections factors for diffusion & dislocation creep
		mfd = exp(mat->mfc*mf);
		mfn = exp(mat->mfc*mf*mat->n);
	}

	// PRESSURE

	// pore pressure
	if(ctrl->gwType == _GW_NONE_) p_pore = 0.0;

	// total pressure
	p_total = p + ctrl->biot*p_pore;

	// assign pressure for viscous laws
	if(ctrl->pLithoVisc)  p_visc = p_lith;
	else                  p_visc = p_total;

	// ELASTICITY
	if(mat->G)
	{
		// Elasticity correction can only DECREASE the viscosity.
		// eta/G << dt (viscous regime)  eta*(dt/(dt + eta/G)) -> eta
		// eta/G >> dt (elastic regime)  eta*(dt/(dt + eta/G)) -> G*dt < eta
		// Elasticity doesn't normally interact with the bottom viscosity limit,
		// instead it rather acts as a smooth limiter for maximum viscosity.

		ctx->A_els = 1.0/(mat->G*dt)/2.0;
	}

	// LINEAR DIFFUSION CREEP (NEWTONIAN)
	if(mat->Bd)
	{
		Q          = (mat->Ed + p_visc*mat->Vd)/RT;
		ctx->A_dif =  mat->Bd*exp(-Q)*mfd;
	}

	// PS-CREEP
	else if(mat->Bps && T)
	{
		Q          = mat->Eps/RT;
		ctx->A_dif = mat->Bps*exp(-Q)/T/pow(mat->d, 3.0);
	}

	// UPPER BOUND CREEP
	if(ctrl->eta_max)
	{
		ctx->A_max = 1.0/(ctrl->eta_max)/2.0;
	}

	// DISLOCATION CREEP (POWER LAW)
	if(mat->Bn)
	{
		Q          = (mat->En + p_visc*mat->Vn)/RT;
		ctx->N_dis =  mat->n;
		ctx->A_dis =  mat->Bn*exp(-Q)*mfn;
	}

	// DC-CREEP
	else if(mat->Bdc && T)
	{
		Q          = mat->Edc/RT;
		ctx->N_dis = Q;
		ctx->A_dis = mat->Bdc*exp(-Q*log(mat->Rdc))*pow(mat->mu, -Q);
	}

	// PEIERLS CREEP (LOW TEMPERATURE RATE-DEPENDENT PLASTICITY, POWER-LAW APPROXIMATION)
	if(mat->Bp && T)
	{
		Q           = (mat->Ep + p_visc*mat->Vp)/RT;
		ctx->N_prl =  Q*pow(1.0-mat->gamma, mat->q-1.0)*mat->q*mat->gamma;
		ctx->A_prl =  mat->Bp/pow(mat->gamma*mat->taup, ctx->N_prl)*exp(-Q*pow(1.0-mat->gamma, mat->q));
	}

	if(PetscIsInfOrNanScalar(ctx->A_dif)) ctx->A_dif = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_dis)) ctx->A_dis = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_prl)) ctx->A_prl = 0.0;

	// PLASTICITY
	if(!mat->ch && !mat->fr)
	{
		PetscFunctionReturn(0); // return if no plasticity is set
	}

	// apply strain softening to friction and cohesion
	ch = applyStrainSoft(soft, mat->chSoftID, APS, mat->ch);
	fr = applyStrainSoft(soft, mat->frSoftID, APS, mat->fr);

	// fit to limits
	if(ch < ctrl->minCh) ch = ctrl->minCh;
	if(fr < ctrl->minFr) fr = ctrl->minFr;

	// override total pressure with lithostatic if requested
	if(ctrl->pLithoPlast)
	{
		// Use lithostatic, rather than dynamic pressure to evaluate yield stress
		// This converges better, but does not result in localization of deformation & shear banding,
		// so only apply it for large-scale simulations where plasticity does not matter much

		p_total = p_lith;
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

		p_upper = -( p_lith + ch * cos(fr))/(sin(fr) - 1.0); // compression
		p_lower = -(-p_lith + ch * cos(fr))/(sin(fr) + 1.0); // extension

		if(p_total > p_upper) p_total = p_upper;
		if(p_total < p_lower) p_total = p_lower;
	}

	// compute cohesion and friction coefficient
	ch = cos(fr)*ch;
	fr = sin(fr);

	// compute effective mean stress
	dP = (p_total - p_pore);

	// compute yield stress
	if(dP < 0.0) ctx->taupl =         ch; // Von-Mises model for extension
	else         ctx->taupl = dP*fr + ch; // Drucker-Prager model for compression

	// correct for ultimate yield stress (if defined)
	if(ctrl->tauUlt) { if(ctx->taupl > ctrl->tauUlt) ctx->taupl = ctrl->tauUlt; }

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getStressCell"
PetscErrorCode getStressCell(
		SolVarCell  *svCell, // solution variables
		ConstEqCtx  *ctx,    // evaluation context
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz)    // ...
{
	// compute stress, plastic strain-rate and shear heating term on cell

	SolVarDev   *svDev;
	PetscScalar  txx, tyy, tzz;

	PetscFunctionBegin;

	// access deviatoric variables
	svDev = &svCell->svDev;

	// compute deviatoric stresses
	svCell->sxx = 2.0*svDev->eta*dxx;
	svCell->syy = 2.0*svDev->eta*dyy;
	svCell->szz = 2.0*svDev->eta*dzz;

	// compute plastic strain-rate components
//	txx = cfpl*dxx;
//	tyy = cfpl*dyy;
//	tzz = cfpl*dzz;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = 0.5*(txx*txx + tyy*tyy + tzz*tzz);

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	txx = svCell->dxx - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tyy = svCell->dyy - svDev->I2Gdt*(svCell->syy - svCell->hyy);
	tzz = svCell->dzz - svDev->I2Gdt*(svCell->szz - svCell->hzz);

	// compute shear heating term contribution
	svDev->Hr = (txx*svCell->sxx + tyy*svCell->syy + tzz*svCell->szz);


/*
	svDev->eta = ctx->eta + ctrl->eta_min;

	// update results
	ctx->eta    += phRat*eta;    // effective viscosity
	ctx->eta_cr += phRat*eta_cr; // creep viscosity
	ctx->eta_vp += phRat*eta_vp; // visco-plastic viscosity
	ctx->DIIdif += phRat*DIIdif; // diffusion creep strain rate
	ctx->DIIdis += phRat*DIIdis; // dislocation creep strain rate
	ctx->DIIprl += phRat*DIIprl; // Peierls creep strain rate
	ctx->DIIpl  += phRat*DIIpl;  // plastic strain rate
	ctx->yield  += phRat*taupl;  // plastic yield rate

*/



	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// compute stress, plastic strain-rate and shear heating term on edge
#undef __FUNCT__
#define __FUNCT__ "getStressEdge"
PetscErrorCode getStressEdge(
		SolVarEdge  *svEdge, // solution variables
		PetscScalar  d,      // effective shear strain rate component
		PetscScalar  cfpl)   // plastic scaling coefficient
{

	SolVarDev   *svDev;
	PetscScalar  t;

	PetscFunctionBegin;

	// access deviatoric variables
	svDev = &svEdge->svDev;

	// compute shear stress
	svEdge->s = 2.0*svDev->eta*d;

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
PetscScalar getConsEqRes(PetscScalar eta, void *pctx)
{
	// Compute residual of the nonlinear visco-elastic constitutive equation:

	PetscScalar tauII, DIIels, DIIdif, DIImax, DIIdis, DIIprl;

	// access context
	ConstEqCtx *ctx = (ConstEqCtx*)pctx;

	// compute stress
	tauII = 2.0*eta*ctx->DII;

	// creep strain rates
	DIIels = ctx->A_els*tauII;                  // elasticity
	DIIdif = ctx->A_dif*tauII;                  // diffusion
	DIImax = ctx->A_max*tauII;                  // upper bound
	DIIdis = ctx->A_dis*pow(tauII, ctx->N_dis); // dislocation
	DIIprl = ctx->A_prl*pow(tauII, ctx->N_prl); // Peierls

	// residual function (r)
	// r < 0 if eta > solution (negative on overshoot)
	// r > 0 if eta < solution (positive on undershoot)

	return 1.0 - (DIIels + DIIdif + DIImax + DIIdis + DIIprl)/ctx->DII;
}
//---------------------------------------------------------------------------
PetscScalar applyStrainSoft(
		Soft_t      *soft, // material softening laws
		PetscInt     ID,   // softening law ID
		PetscScalar  APS,  // accumulated plastic strain
		PetscScalar  par) // softening parameter
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
PetscScalar getI2Gdt(
		PetscInt     numPhases, // number phases
		Material_t  *phases,    // phase parameters
		PetscScalar *phRat,     // phase ratios in the control volume
		PetscScalar  dt)        // time step
{
	// compute inverse deviatoric elastic parameter

	PetscInt    i;
	PetscScalar I2Gdt, Gavg;

	Gavg  = 0.0;
	I2Gdt = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		Gavg += phRat[i]*phases[i].G;
	}

	if(Gavg) I2Gdt = 1.0/Gavg/dt/2.0;

	return I2Gdt;
}
//---------------------------------------------------------------------------
//.............................. PHASE DIAGRAM  .............................
//---------------------------------------------------------------------------
// get the density from a phase diagram
#undef __FUNCT__
#define __FUNCT__ "setDataPhaseDiagram"
PetscErrorCode setDataPhaseDiagram(
		PData       *pd,
		PetscScalar  p,
		PetscScalar  T,
		char         pdn[])
{
    PetscInt       	i,j,i_pd,indT[2],indP[2],ind[4],found;
    PetscScalar    	fx0,fx1,weight[4];
	PetscScalar 	minP, dP, minT, dT;

	PetscFunctionBegin;

	// Get the correct phase diagram
	for(i=0; i<_max_num_pd_; i++)
	{
		i_pd  = -1;
		found = 1;
		if(!pd->rho_pdns[0][i])
		{
			// We found an empty phase diagram spot
		}
		else
		{
			for(j=0; j<_pd_name_sz_; j++)
			{
				if((pd->rho_pdns[j][i] != pdn[j]))
				{
					found = 0;
					break;
				}
			}
			if(found==1)
			{
				i_pd = i;  // Store the column of the buffer
				break;
			}
		}
	}

	if(i_pd<0)
	{
		pd->rho = 0;
		PetscFunctionReturn(0);
	}

	// Take absolute value of pressure
	if(p<0)
	{
		// p = -1*p;
		p = 0;
	}

	// copy in temporary variables for code readability (and speed in fact as less reading from memory is triggered)
	minP 	=	pd->minP[i_pd];
	minT 	=	pd->minT[i_pd];
	dP 		=	pd->dP[i_pd];
	dT 		=	pd->dT[i_pd];

	indT[0] = (PetscInt)floor((T-minT)/dT);
	indT[1] = indT[0] + 1;

	indP[0] = (PetscInt)floor((p-(minP))/dP);
	indP[1] = indP[0] + 1;

	weight[0] = (    (dP*((PetscScalar)indP[1]) + minP) - p) / ((dP*((PetscScalar)indP[1]) + minP) - (dP*((PetscScalar)indP[0]) +minP) );
	weight[1] = (p - (dP*((PetscScalar)indP[0]) + minP)    ) / ((dP*((PetscScalar)indP[1]) + minP) - (dP*((PetscScalar)indP[0]) +minP) );
	weight[2] = (    (dT*((PetscScalar)indT[1]) + minT) - T) / ((dT*((PetscScalar)indT[1]) + minT) - (dT*((PetscScalar)indT[0]) +minT) );
	weight[3] = (T - (dT*((PetscScalar)indT[0]) + minT)    ) / ((dT*((PetscScalar)indT[1]) + minT) - (dT*((PetscScalar)indT[0]) +minT) );

	if(indT[1]>(pd->nT[i_pd]))
	{
		indT[0] = pd->nT[i_pd]-1;
		indT[1] = pd->nT[i_pd];
		weight[2] = 1;
		weight[3] = 0;
	}
	if(indP[1]>(pd->nP[i_pd]))
	{
		indP[0] = pd->nP[i_pd]-1;
		indP[1] = pd->nP[i_pd];
		weight[0] = 1;
		weight[1] = 0;
	}
	if(indT[0]<1)
	{
		indT[0] = 0;
		indT[1] = 1;
		weight[2] = 0;
		weight[3] = 1;
	}
	if(indP[0]<1)
	{
		indP[0] = 0;
		indP[1] = 1;
		weight[0] = 0;
		weight[1] = 1;
	}
	ind[0] = pd->nT[i_pd] * (indP[0]-1) + indT[0];
	ind[1] = pd->nT[i_pd] * (indP[0]-1) + indT[1];
	ind[2] = pd->nT[i_pd] * (indP[1]-1) + indT[0];
	ind[3] = pd->nT[i_pd] * (indP[1]-1) + indT[1];
	if(ind[0]<0)
	{
		ind[0] = 0;
		ind[1] = 1;
	}
	if(ind[3]>pd->nT[i_pd]*pd->nP[i_pd])
	{
		ind[2] = pd->nT[i_pd]*pd->nP[i_pd]-1;
		ind[3] = pd->nT[i_pd]*pd->nP[i_pd];
	}

	// Interpolate density
	fx0 = weight[0] * pd->rho_v[ind[0]][i_pd] + weight[1] * pd->rho_v[ind[2]][i_pd];
	fx1 = weight[0] * pd->rho_v[ind[1]][i_pd] + weight[1] * pd->rho_v[ind[3]][i_pd];
	pd->rho = weight[2] * fx0           + weight[3] * fx1;

	// Interpolate melt fraction if present
	if(pd->numProps[i_pd] == 4 )
	{
		fx0 = weight[0] * pd->Me_v[ind[0]][i_pd] + weight[1] * pd->Me_v[ind[2]][i_pd];
		fx1 = weight[0] * pd->Me_v[ind[1]][i_pd] + weight[1] * pd->Me_v[ind[3]][i_pd];
		pd->mf  = weight[2] * fx0      + weight[3] * fx1;
	}
	// Interpolate mf + rho fluid if present
	else if(pd->numProps[i_pd] == 5)
	{
		fx0 = weight[0] * pd->Me_v[ind[0]][i_pd] + weight[1] * pd->Me_v[ind[2]][i_pd];
		fx1 = weight[0] * pd->Me_v[ind[1]][i_pd] + weight[1] * pd->Me_v[ind[3]][i_pd];
		pd->mf  = weight[2] * fx0      + weight[3] * fx1;

		fx0 = weight[0] * pd->rho_f_v[ind[0]][i_pd] + weight[1] * pd->rho_f_v[ind[2]][i_pd];
		fx1 = weight[0] * pd->rho_f_v[ind[1]][i_pd] + weight[1] * pd->rho_f_v[ind[3]][i_pd];
		pd->rho_f  = weight[2] * fx0   + weight[3] * fx1;

	}
	// No melt fraction
	else
	{
		pd->mf = 0;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
