/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "constEq.h"
#include "tssolve.h"
#include "surf.h"
#include "phase.h"
#include "JacRes.h"
#include "meltParam.h"
#include "tools.h"
#include "phase_transition.h"
#include "scaling.h"
#include "parsing.h"
#include "bc.h"
#include "dike.h"
//---------------------------------------------------------------------------
PetscErrorCode setUpConstEq(ConstEqCtx *ctx, JacRes *jr)
{
	// setup constitutive equation evaluation context parameters

	PetscFunctionBeginUser;

	ctx->bc        =  jr->bc;              // boundary conditions for inflow velocity
	ctx->numPhases =  jr->dbm->numPhases;  // number phases
	ctx->phases    =  jr->dbm->phases;     // phase parameters
	ctx->numDike   =  jr->dbdike->numDike; // number of dikes
	ctx->matDike   =  jr->dbdike->matDike; // dike properties
	ctx->soft      =  jr->dbm->matSoft;    // material softening laws
	ctx->ctrl      = &jr->ctrl;            // control parameters
	ctx->Pd        =  jr-> Pd;             // phase diagram data
	ctx->dt        =  jr->ts->dt;          // time step
	ctx->PhaseTrans = jr->dbm->matPhtr;    // phase transition
	ctx->numPhtr   =  jr->dbm->numPhtr;    // number of phase transition laws
	ctx->scal      =  jr->scal;            // scaling
	ctx->stats[0]  =  0.0;                 // total number of [starts, ...
	ctx->stats[1]  =  0.0;                 // ... successes,
	ctx->stats[2]  =  0.0;                 // ... iterations]

	// set average surface topography for depth computation
	ctx->avg_topo = DBL_MAX;

	if(jr->surf->AirPhase != -1)
	{
		ctx->avg_topo = jr->surf->avg_topo;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode setUpCtrlVol(
		ConstEqCtx  *ctx,    // context
		PetscScalar *phRat,  // phase ratios in the control volume
		SolVarDev   *svDev,  // deviatoric variables
		SolVarBulk  *svBulk, // volumetric variables
		PetscScalar  p,      // pressure
		PetscScalar  p_lith, // lithostatic pressure
		PetscScalar  p_pore, // pore pressure
		PetscScalar  T,      // temperature
		PetscScalar  DII,    // effective strain rate
		PetscScalar  z,      // z-coordinate of control volume
		PetscScalar  Le)     // characteristic element size
{
	// setup control volume parameters

	PetscFunctionBeginUser;

	ctx->phRat  = phRat;  // phase ratios in the control volume
	ctx->svDev  = svDev;  // deviatoric variables
	ctx->svBulk = svBulk; // volumetric variables
	ctx->p      = p;      // pressure
	ctx->p_lith = p_lith; // lithostatic pressure
	ctx->p_pore = p_pore; // pore pressure
	ctx->T      = T;      // temperature
	ctx->DII    = DII;    // effective strain rate
	ctx->Le     = Le;     // characteristic element size


	// compute depth below the free surface
	// WARNING! "depth" is loosely defined for large topography variations
	ctx->depth = 0.0;

	if(ctx->avg_topo != DBL_MAX && z != DBL_MAX)
	{
		ctx->depth = ctx->avg_topo - z;

		if(ctx->depth < 0.0)
		{
			ctx->depth = 0.0;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode setUpPhase(ConstEqCtx *ctx, PetscInt ID)
{
	// setup phase parameters for deviatoric constitutive equation
	// evaluate dependence on constant parameters (pressure, temperature)

	Material_t  *mat;
	Soft_t      *soft;
	Controls    *ctrl;
	PData       *Pd;
	PetscScalar  APS, Le, dt, p, p_lith, p_pore, T, mf, mfd, mfn;
	PetscScalar  Q, RT, ch, fr, p_visc, p_upper, p_lower, dP, p_total;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	mat    = ctx->phases + ID;
	soft   = ctx->soft;
	ctrl   = ctx->ctrl;
	Pd     = ctx->Pd;
	APS    = ctx->svDev->APS;
	Le     = ctx->Le;
	dt     = ctx->dt;
	p      = ctx->p;
	p_lith = ctx->p_lith;
	p_pore = ctx->p_pore;
	T      = ctx->T;
	mf     = 0.0;

	p 	   = p + ctrl->pShift;		// add pressure shift to pressure field

	if(mat->pdAct == 1)
	{
		// compute melt fraction from phase diagram
		ierr = setDataPhaseDiagram(Pd, p, T, mat->pdn);CHKERRQ(ierr);

		// store melt fraction
		mf = Pd->mf;
	}

	//if(strcmp(mat->Melt_Parametrization,"none") & ctrl->melt_feedback == 1)
//	{
	//	mf = Compute_Melt_Fraction(p, T ,mat,ctx);
	//}

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
	ctx->A_fk  = 0.0; // Frank-Kamenetzky constant
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
		ctx->A_dif = mat->Bd*exp(-Q)*mfd;
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

	// Frank-Kamenetzky Viscosity
	if(mat->gamma_fk && T)
	{
		ctx->A_fk = 1.0/(mat->eta_fk*exp(-mat->gamma_fk*(T-mat->TRef_fk)))/2.0;
	}

	if(PetscIsInfOrNanScalar(ctx->A_dif)) ctx->A_dif = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_dis)) ctx->A_dis = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_prl)) ctx->A_prl = 0.0;
	if(PetscIsInfOrNanScalar(ctx->A_fk))  ctx->A_fk  = 0.0;

	// PLASTICITY
	if(!mat->ch && !mat->fr)
	{
		PetscFunctionReturn(0); // return if no plasticity is set
	}

	// apply strain softening to friction and cohesion
	ch = applyStrainSoft(soft, mat->chSoftID, APS, Le, mat->ch);
	fr = applyStrainSoft(soft, mat->frSoftID, APS, Le, mat->fr);

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

	// store regularization viscosity
	ctx->eta_vp = mat->eta_vp;

	// correct for ultimate yield stress (if defined)
	if(ctrl->tauUlt) { if(ctx->taupl > ctrl->tauUlt) ctx->taupl = ctrl->tauUlt; }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode devConstEq(ConstEqCtx *ctx)
{
	// evaluate deviatoric constitutive equations in control volume

	Controls    *ctrl;
	PetscScalar *phRat;
	SolVarDev   *svDev;
	Material_t  *mat;
	PetscInt     i, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ctrl      = ctx->ctrl;
	numPhases = ctx->numPhases;
	phRat     = ctx->phRat;
	svDev     = ctx->svDev;

	// zero out results
	ctx->eta    = 0.0; // effective viscosity
	ctx->eta_cr = 0.0; // creep viscosity
	ctx->DIIdif = 0.0; // diffusion creep strain rate
	ctx->DIIdis = 0.0; // dislocation creep strain rate
	ctx->DIIprl = 0.0; // Peierls creep strain rate
	ctx->DIIfk  = 0.0; // Frank-Kamenetzky strain rate
	ctx->DIIpl  = 0.0; // plastic strain rate
	ctx->yield  = 0.0; // yield stress

	// zero out stabilization and viscoplastic viscosity
	svDev->eta_st = 0.0;

	// viscous initial guess
	if(ctrl->initGuess)
	{
		ctx->eta    = ctrl->eta_ref;
		ctx->eta_cr = ctrl->eta_ref;
		ctx->DIIdif = 1.0;

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

			// update stabilization and viscoplastic viscosity
			mat            = ctx->phases + i;
			svDev->eta_st += phRat[i]*mat->eta_st;
		}
	}

	// normalize strain rates
	if(ctx->DII)
	{
		ctx->DIIdif /= ctx->DII;
		ctx->DIIdis /= ctx->DII;
		ctx->DIIprl /= ctx->DII;
		ctx->DIIfk  /= ctx->DII;
		ctx->DIIpl  /= ctx->DII;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode getPhaseVisc(ConstEqCtx *ctx, PetscInt ID)
{
	// compute phase viscosities and strain rate partitioning

	Controls    *ctrl;
	PetscInt    it, conv;
	PetscScalar eta_min, eta_mean, eta, eta_cr, tauII, taupl, DII;
	PetscScalar DIIdif, DIImax, DIIdis, DIIprl, DIIpl, DIIplc, DIIfk, DIIvs, phRat;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_max, inv_eta_dis, inv_eta_prl, inv_eta_fk, inv_eta_min;

	PetscFunctionBeginUser;

	// access context
	ctrl   = ctx->ctrl;      // global controls
	phRat  = ctx->phRat[ID]; // phase ratio
	taupl  = ctx->taupl;     // plastic yield stress
	DII    = ctx->DII;       // effective strain rate

	// initialize
	it     = 1;
	conv   = 1;
	DIIpl  = 0.0;
	eta_cr = 0.0;

	//===========
	// PLASTICITY
	//===========

	if(taupl && DII)
	{
		// get initial yield stress and viscosity
		tauII = taupl;
		eta   = tauII/(2.0*DII);

		// compute initial plastic strain rate
		DIIpl = getConsEqRes(eta, ctx);

		// reset if plasticity is not active
		if(DIIpl < 0.0)
		{
			DIIpl = 0.0;
		}
		else if(ctx->eta_vp)
		{
			//================================================================
			// solve regularized visco-plastic strain by fixed-point iteration
			//================================================================
			do
			{
				// get regularized yield stress and viscosity
				tauII = taupl + 2.0*ctx->eta_vp*DIIpl;
				eta   = tauII/(2.0*DII);

			    // store current strain strain rate
			    DIIplc = DIIpl;

				// compute updated plastic strain rate
				DIIpl = getConsEqRes(eta, ctx);

				// set convergence flag
				conv = (PetscAbsScalar((DIIpl - DIIplc)/DII) <= ctrl->lrtol);

			} while(!conv && ++it < ctrl->lmaxit);
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
		inv_eta_fk  = 0.0;

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
		// Frank-Kamenetzky
		if(ctx->A_fk)  inv_eta_fk  = 2.0*ctx->A_fk;

		// get minimum viscosity (upper bound)
		inv_eta_min                               = inv_eta_els;
		if(inv_eta_dif > inv_eta_min) inv_eta_min = inv_eta_dif;
		if(inv_eta_max > inv_eta_min) inv_eta_min = inv_eta_max;
		if(inv_eta_dis > inv_eta_min) inv_eta_min = inv_eta_dis;
		if(inv_eta_prl > inv_eta_min) inv_eta_min = inv_eta_prl;
		if(inv_eta_fk  > inv_eta_min) inv_eta_min = inv_eta_fk;
		eta_min = 1.0/inv_eta_min;

		// get quasi-harmonic mean (lower bound)
		eta_mean = 1.0/(inv_eta_els + inv_eta_dif + inv_eta_max + inv_eta_dis + inv_eta_prl + inv_eta_fk);

		// NOTE: if closed-form solution exists, it is equal to lower bound
		// If only one mechanism is active, then both bounds are coincident
		// Bisection will return immediately if closed-form solution exists

		// apply bisection algorithm to nonlinear scalar equation
		conv = solveBisect(eta_mean, eta_min, ctrl->lrtol*DII, ctrl->lmaxit, eta, it, getConsEqRes, ctx);

		// compute stress
		tauII = 2.0*eta*DII;
	}

	// update iteration statistics
	ctx->stats[0] += 1.0;               // start counter
	ctx->stats[1] += (PetscScalar)conv; // convergence counter
	ctx->stats[2] += (PetscScalar)it;   // iteration counter

	// compute strain rates
	DIIdif = ctx->A_dif*tauII;                  // diffusion
	DIImax = ctx->A_max*tauII;                  // upper bound
	DIIdis = ctx->A_dis*pow(tauII, ctx->N_dis); // dislocation
	DIIprl = ctx->A_prl*pow(tauII, ctx->N_prl); // Peierls
	DIIfk  = ctx->A_fk*tauII;                   // Frank-Kamenetzky
	DIIvs  = DIIdif + DIImax + DIIdis + DIIprl + DIIfk; // viscous (total)

	// compute creep viscosity
	if(DIIvs) eta_cr = tauII/DIIvs/2.0;

	// update results
	ctx->eta    += phRat*eta;    // effective viscosity
	ctx->eta_cr += phRat*eta_cr; // creep viscosity
	ctx->DIIdif += phRat*DIIdif; // diffusion creep strain rate
	ctx->DIIdis += phRat*DIIdis; // dislocation creep strain rate
	ctx->DIIprl += phRat*DIIprl; // Peierls creep strain rate
	ctx->DIIfk  += phRat*DIIfk;  // Frank-Kamenetzky
	ctx->DIIpl  += phRat*DIIpl;  // plastic strain rate
	ctx->yield  += phRat*taupl;  // plastic yield stress

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscScalar getConsEqRes(PetscScalar eta, void *pctx)
{
	// compute residual of the nonlinear visco-elastic constitutive equation

	PetscScalar tauII, DIIels, DIIdif, DIImax, DIIdis, DIIprl, DIIfk;

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
	DIIfk  = ctx->A_fk*tauII;                   // Frank-Kamenetzky

	// residual function (r)
	// r < 0 if eta > solution (negative on overshoot)
	// r > 0 if eta < solution (positive on undershoot)

	return ctx->DII - (DIIels + DIIdif + DIImax + DIIdis + DIIprl + DIIfk);
}
//---------------------------------------------------------------------------
PetscScalar applyStrainSoft(
		Soft_t      *soft, // material softening laws
		PetscInt     ID,   // softening law ID
		PetscScalar  APS,  // accumulated plastic strain
		PetscScalar  Le,   // characteristic element size
		PetscScalar  par)  // softening parameter
{
	// apply strain softening to a parameter (friction, cohesion)

	PetscScalar  k;  // dt
	PetscScalar  A, APS1, APS2, Lm; 
	Soft_t      *s;

	// check whether softening is defined
	if(ID == -1) return par;

	// access parameters
	s    = soft + ID;
	APS1 = s->APS1;
	APS2 = s->APS2;
	A    = s->A;
	Lm   = s->Lm;

	// Fracture Energy Regularization
	if(Lm)
	{
		APS1 *= Le/Lm;
		APS2 *= Le/Lm;
	}

	// compute scaling ratio
	if(APS <= APS1)               k = 1.0;
	if(APS >  APS1 && APS < APS2) k = 1.0 - A*((APS - APS1)/(APS2 - APS1));
	if(APS >= APS2)               k = 1.0 - A;

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
PetscErrorCode volConstEq(ConstEqCtx *ctx)
{
	// evaluate volumetric constitutive equations in control volume
	Controls    *ctrl;
	PData       *Pd;
	SolVarBulk  *svBulk;
	Material_t  *mat, *phases;
	PetscInt     i, numPhases;
	PetscScalar *phRat, dt, p, depth, T, cf_comp, cf_therm, Kavg, rho;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	ctrl      = ctx->ctrl;
	Pd        = ctx->Pd;
	svBulk    = ctx->svBulk;
	numPhases = ctx->numPhases;
	phases    = ctx->phases;
	phRat     = ctx->phRat;
	depth     = ctx->depth;
	dt        = ctx->dt;
	p         = ctx->p;
	T         = ctx->T;

	p         = p+ctrl->pShift;

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
				ierr = setDataPhaseDiagram(Pd, p, T, mat->pdn); CHKERRQ(ierr);

				svBulk->mf     += phRat[i]*Pd->mf;

				if(mat->rho_melt)
				{
					svBulk->rho_pf += phRat[i]*mat->rho_melt;
				}
				else
				{
					svBulk->rho_pf += phRat[i]*Pd->rho_f;
				}
			}

			// initialize
			cf_comp  = 1.0;
			cf_therm = 1.0;

			// elastic compressibility correction (Murnaghan's equation)
			// ro/ro_0 = (1 + Kb'*P/Kb)^(1/Kb')
			if(mat->Kb)
			{
				Kavg += phRat[i]*mat->Kb;

				if(mat->Kp) cf_comp = pow(1.0 + mat->Kp*(p/mat->Kb), 1.0/mat->Kp);
				else        cf_comp = 1.0 + p/mat->Kb;
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
			else if(mat->pdAct == 1 && !mat->Phase_Diagram_melt)
			{
				// Compute density from phase diagram, while also taking the actual melt content into account
				PetscScalar mf;
 				
				mf = Pd->mf;
				if (mf > ctrl->mfmax){ mf = ctrl->mfmax; }

				rho = (mf * Pd->rho_f) + ((1.0 - mf ) * Pd->rho);
			}
			else if(mat->pdAct == 1 && mat->Phase_Diagram_melt)
			{
				PetscScalar mf;

				mf = Pd->mf;
				if (mf > ctrl->mfmax){ mf = ctrl->mfmax; }
				rho = mat->rho*cf_comp*cf_therm;
				rho = (Pd->mf * mat->rho_melt) + ((1-Pd->mf) * rho);

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
PetscErrorCode cellConstEq(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarCell  *svCell, // solution variables
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz,    // ...
		PetscScalar &sxx,    // Cauchy stress components
		PetscScalar &syy,    // ...
		PetscScalar &szz,    // ...
		PetscScalar &gres,   // volumetric residual
		PetscScalar &rho,    // effective density
		PetscScalar &dikeRHS) // dike RHS for gres calculation
{
	// evaluate constitutive equations on the cell

	SolVarDev   *svDev;
	SolVarBulk  *svBulk;
	Controls    *ctrl;
	PetscScalar  eta_st, ptotal, txx, tyy, tzz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	svDev  = ctx->svDev;
	svBulk = ctx->svBulk;
	ctrl   = ctx->ctrl;

	// evaluate deviatoric constitutive equation
	ierr = devConstEq(ctx); CHKERRQ(ierr);

	// evaluate volumetric constitutive equation
	ierr = volConstEq(ctx); CHKERRQ(ierr);

	// get stabilization viscosity
	if(ctrl->initGuess) eta_st = 0.0;
	else                eta_st = svDev->eta_st;

	// compute stabilization stresses
	sxx = 2.0*eta_st*svCell->dxx;
	syy = 2.0*eta_st*svCell->dyy;
	szz = 2.0*eta_st*svCell->dzz;

	// compute history shear stress
	svCell->sxx = 2.0*ctx->eta*dxx;
	svCell->syy = 2.0*ctx->eta*dyy;
	svCell->szz = 2.0*ctx->eta*dzz;

	// compute plastic strain-rate components
	txx = ctx->DIIpl*dxx;
	tyy = ctx->DIIpl*dyy;
	tzz = ctx->DIIpl*dzz;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = 0.5*(txx*txx + tyy*tyy + tzz*tzz);

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	txx = svCell->dxx - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tyy = svCell->dyy - svDev->I2Gdt*(svCell->syy - svCell->hyy);
	tzz = svCell->dzz - svDev->I2Gdt*(svCell->szz - svCell->hzz);

	// compute shear heating term contribution
	svDev->Hr =
		txx*svCell->sxx + tyy*svCell->syy + tzz*svCell->szz +
		sxx*svCell->dxx + syy*svCell->dyy + szz*svCell->dzz;

	// compute total viscosity
	svDev->eta = ctx->eta + eta_st;

	// get total pressure (effective pressure + pore pressure)
	ptotal = ctx->p + ctrl->biot*ctx->p_pore;

	// compute total Cauchy stresses
	sxx += svCell->sxx - ptotal;
	syy += svCell->syy - ptotal;
	szz += svCell->szz - ptotal;

	// save output variables
	svCell->eta_cr = ctx->eta_cr; // creep viscosity
	svCell->DIIdif = ctx->DIIdif; // relative diffusion creep strain rate
	svCell->DIIdis = ctx->DIIdis; // relative dislocation creep strain rate
	svCell->DIIprl = ctx->DIIprl; // relative Peierls creep strain rate
	svCell->DIIfk  = ctx->DIIfk;  // relative Frank-Kamenetzky strain rate
	svCell->DIIpl  = ctx->DIIpl;  // relative plastic strain rate
	svCell->yield  = ctx->yield;  // average yield stress in control volume


	if(ctrl->actExp && ctrl->actDike)
    {
        gres= -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta + svBulk->alpha*(ctx->T - svBulk->Tn)/ctx->dt + dikeRHS;
    }
	else if(ctrl->actDike)
    {
        gres = -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta + dikeRHS;
    }
	else if(ctrl->actExp)
    {
        gres = -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta + svBulk->alpha*(ctx->T - svBulk->Tn)/ctx->dt;
    }
	else
    {
        gres = -svBulk->IKdt*(ctx->p - svBulk->pn) - svBulk->theta;
    }

	// store effective density
	rho = svBulk->rho;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode edgeConstEq(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarEdge  *svEdge, // solution variables
		PetscScalar  d,      // effective shear strain rate component
		PetscScalar &s)      // Cauchy stress component
{
	// evaluate constitutive equations on the edge

	SolVarDev   *svDev;
	PetscScalar  t, eta_st;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	svDev = &svEdge->svDev;

	// evaluate deviatoric constitutive equation
	ierr = devConstEq(ctx); CHKERRQ(ierr);

	// get stabilization viscosity
	if(ctx->ctrl->initGuess) eta_st = 0.0;
	else                     eta_st = svDev->eta_st;

	// compute stabilization stress
	s = 2.0*eta_st*svEdge->d;

	// compute history shear stress
	svEdge->s = 2.0*ctx->eta*d;

	// compute plastic strain-rate component
	t = ctx->DIIpl*d;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = t*t;

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	t = svEdge->d - svDev->I2Gdt*(svEdge->s - svEdge->h);

	// compute shear heating term contribution
	svDev->Hr = 2.0*t*svEdge->s + 2.0*svEdge->d*s;

	// compute total viscosity
	svDev->eta = ctx->eta + eta_st;

	// compute total stress
	s += svEdge->s;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode checkConvConstEq(ConstEqCtx *ctx)
{
	// check convergence of constitutive equations
	LLD         ndiv, nit;
	PetscScalar stats[3] = {1.0, 1.0, 1.0};

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// exchange convergence statistics
	// total number of [starts, successes, iterations]
	ierr = MPI_Reduce(stats, ctx->stats, 3, MPIU_SCALAR, MPI_SUM, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

	// compute number of diverged equations and average iteration count
	ndiv = (LLD)(stats[0] - stats[1]);
	nit  = (LLD)(stats[2] / stats[0]);

	if(ndiv)
	{
		PetscPrintf(PETSC_COMM_WORLD,"*****************************************************************************\n");
		PetscPrintf(PETSC_COMM_WORLD,"Warning! Number of diverged points : %lld \n", ndiv);
		PetscPrintf(PETSC_COMM_WORLD,"Average iteration count            : %lld \n", nit);
		PetscPrintf(PETSC_COMM_WORLD,"*****************************************************************************\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//.............................. PHASE DIAGRAM  .............................
//---------------------------------------------------------------------------
// get the density from a phase diagram
PetscErrorCode setDataPhaseDiagram(
		PData       *pd,
		PetscScalar  p,
		PetscScalar  T,
		char         pdn[])
{
    PetscInt       	i,j,i_pd,indT[2],indP[2],ind[4],found;
    PetscScalar    	fx0,fx1,weight[4];
	PetscScalar 	minP, dP, minT, dT;

	PetscFunctionBeginUser;

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
