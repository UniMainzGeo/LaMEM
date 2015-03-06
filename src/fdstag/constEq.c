//---------------------------------------------------------------------------
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "solVar.h"
#include "constEq.h"
#include "dfzero.h"
//---------------------------------------------------------------------------
// * add different viscosity averaging methods (echo info to output)
// * make sure that Peierls creep is deactivated for isothermal analysis
// ...
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "ConstEqCtxSetup"
PetscErrorCode ConstEqCtxSetup(
	ConstEqCtx  *ctx,  // evaluation context
	Material_t  *mat,  // phase parameters
	MatParLim   *lim,  // phase parameters limits
	PetscScalar  DII,  // effective strain-rate
	PetscScalar  APS,  // accumulated plastic strain
	PetscScalar  dt,   // time step
	PetscScalar  p,    // pressure
	PetscScalar  T)    // temperature
{
	// setup nonlinear constitutive equation evaluation context
	// evaluate dependence on constant parameters (pressure, temperature)

	PetscInt    ln, nl;
	PetscScalar Q, RT, ch, fr;

	PetscFunctionBegin;

	if(T) RT = lim->Rugc*T;
	else  RT = 1.0;

	ln = 0;
	nl = 0;

	// use reference strain-rate instead of zero
	if(DII == 0.0) DII = lim->DII_ref;

	// initialize
	ctx->DII   = DII;         // effective strain-rate
	ctx->A_els = 0.0;         // elasticity constant
	ctx->A_dif = 0.0;         // diffusion constant
	ctx->A_dis = 0.0;         // dislocation constant
	ctx->N_dis = 1.0;         // dislocation exponent
	ctx->A_prl = 0.0;         // Peierls constant
	ctx->N_prl = 1.0;         // Peierls exponent
	ctx->taupl = 0.0;         // plastic yield stress
	ctx->cfsol = PETSC_TRUE;  // closed-form solution flag

	// ELASTICITY
	if(mat->G)
	{
		// Elasticity correction can only DECREASE the viscosity.
		// eta/G << dt (viscous regime)  eta*(dt/(dt + eta/G)) -> eta
		// eta/G >> dt (elastic regime)  eta*(dt/(dt + eta/G)) -> G*dt < eta
		// Elasticity doesn't normally interact with the bottom viscosity limit,
		// instead it rather acts as a smooth limiter for maximum viscosity.

		ctx->A_els = 0.5/(mat->G*dt);
		ln++;
	}

	// DIFFUSION CREEP (NEWTONIAN)
	if(mat->Bd)
	{
		Q          = (mat->Ed + p*mat->Vd)/RT;
		ctx->A_dif =  mat->Bd*exp(-Q);
		ln++;
	}

	// DISLOCATION CREEP (POWER LAW)
	if(mat->Bn)
	{
// ACHTUNG
//		if(lim->initGuessFlg == PETSC_TRUE)
//		Q          = 0.0;
//		else
		Q          = (mat->En + p*mat->Vn)/RT;
		ctx->N_dis =  mat->n;
		ctx->A_dis =  mat->Bn*exp(-Q);
		nl++;
	}

	// PEIERLS CREEP (LOW TEMPERATURE RATE-DEPENDENT PLASTICITY, POWER-LAW APPROXIMATION)
	if(mat->Bp)
	{
		Q          = (mat->Ep + p*mat->Vp)/RT;
		ctx->N_prl =  Q*pow(1.0-mat->gamma, mat->q-1.0)*mat->q*mat->gamma;
		ctx->A_prl =  mat->Bp/pow(mat->gamma*mat->taup, ctx->N_prl)*exp(-Q*pow(1.0-mat->gamma, mat->q));
		nl++;
	}

	// apply strain softening to friction and cohesion
	ch = ApplyStrainSoft(mat->chSoft, APS, mat->ch);
	fr = ApplyStrainSoft(mat->frSoft, APS, mat->fr);

	// compute cohesion and friction coefficient
	ch = cos(fr)*ch;
	fr = sin(fr);

	// fit to limits
	if(ch < lim->minCh) ch = lim->minCh;
	if(fr < lim->minFr) fr = lim->minFr;

	// compute yield stress
	if(p < 0.0) ctx->taupl = ch;        // Von-Mises model for extension
	else        ctx->taupl = ch + p*fr; // Drucker-Prager model for compression

	// correct for ultimate yield stress
	if(ctx->taupl > lim->tauUlt) ctx->taupl = lim->tauUlt;

	// set iteration flag (linear + nonlinear, or more than one nonlinear)
	if((nl && ln) || nl > 1) ctx->cfsol = PETSC_FALSE;

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
	MatParLim   *lim,
	PetscScalar *eta,
	PetscScalar *DIIpl)
{
	// ACHTUNG! THIS WILL ONLY WORK FOR:
	//  * LINEAR VISCO-ELASTO-PLASTIC MATERIAL
	//  * POWER-LAW MATERIAL

	PetscScalar eta_ve, inv_eta_els, inv_eta_dif, inv_eta_dis, eta_pl;

	PetscFunctionBegin;

	// initialize
	(*eta)   = 0.0;
	(*DIIpl) = 0.0;

	// set reference viscosity as initial guess (is )
	if(lim->eta_ref && lim->initGuessFlg == PETSC_TRUE)
	{
		(*eta) = lim->eta_ref;

		PetscFunctionReturn(0);
	}

	if(ctx->A_dis)
	{
		// compute power-law viscosity (use reference strain-rate as initial guess)
		if(lim->initGuessFlg == PETSC_TRUE)
		{
			inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(lim->DII_ref, 1.0 - 1.0/ctx->N_dis);
		}
		else
		{
			inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(ctx->DII, 1.0 - 1.0/ctx->N_dis);
		}

		(*eta) = 1.0/inv_eta_dis;
	}
	else
	{
		// elasticity
		inv_eta_els = 2.0*ctx->A_els;

		// diffusion
		inv_eta_dif = 2.0*ctx->A_dif;

		// compute visco-elastic viscosity
		eta_ve = 1.0/(inv_eta_els + inv_eta_dif);

		// get plastic viscosity
		eta_pl = ctx->taupl/(2.0*ctx->DII);

		// check plasticity condition (not for initial guess)
		if(lim->initGuessFlg != PETSC_TRUE && eta_pl && eta_ve > eta_pl)
		{
			if(lim->quasiHarmAvg == PETSC_TRUE)
			{
				(*eta) = 1.0/(1.0/eta_pl + 1.0/eta_ve);
			}
			else
			{
				(*eta) = eta_pl;
			}

			// compute plastic strain rate
			(*DIIpl) = ctx->DII - ctx->taupl/(2.0*eta_ve);

		}
		else
		{
			(*eta) = eta_ve;
		}
	}

	// enforce constraints
	if((*eta) < lim->eta_min) (*eta) = lim->eta_min;
	if((*eta) > lim->eta_max) (*eta) = lim->eta_max;

	PetscFunctionReturn(0);

	//=============================================

/*
	// "isolated" viscosity for each creep mechanism can be computed by
	// assuming that single creep mechanism consumes entire strain rate.

	// Observations:
	// [A] viscosity is LARGER (or equal) than quasi-harmonic mean of isolated viscosities
	// [B] viscosity is SMALLER (or equal) than minimum isolated viscosity
	// [C] residual is POSITIVE for viscosities SMALLER then solution, and vice-versa

	PetscInt    IFLAG;
	PetscScalar B, C, R, RE, AE;
	PetscScalar eta_min, eta_max, inv_eta_min, inv_eta_max;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_dis, inv_eta_prl;
	PetscScalar inv_eta_top, eta_top, inv_eta_bot, eta_bot, etapl, res;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize viscosity limits
	eta_min = lim->eta_min;
	eta_max = lim->eta_max;

	// initialize plastic strain-rate
	(*DIIpl) = 0.0;

	// check for plasticity
	if(ctx->taupl)
	{
		// compute plastic viscosity
		etapl = ctx->taupl/(2.0*ctx->DII);

		// compute residual (which, if positive, has a meaning of plastic strain-rate)
		res = GetConsEqRes(etapl, ctx);

		// check for positive plastic strain-rate
		if(res > 0.0)
		{
			(*eta) = etapl; (*DIIpl) = res; PetscFunctionReturn(0);
		}
		else
		{
			// correct viscosity limits (cannot be larger than plastic viscosity)
			if(eta_min > etapl)
			{
				// minimum is larger than plastic viscosity (return minimum)
				(*eta) = eta_min; PetscFunctionReturn(0);
			}
			if(eta_max > etapl)
			{
				// maximum is larger than plastic viscosity (correcting)
				eta_max = etapl;
			}
		}
	}

	// compute inverses of viscosity limits
	inv_eta_min = 1.0/eta_min;
	inv_eta_max = 1.0/eta_max;

	// initialize inverses of isolated viscosities
	inv_eta_els = 0.0;
	inv_eta_dif = 0.0;
	inv_eta_dis = 0.0;
	inv_eta_prl = 0.0;

	// compute inverses of isolated viscosities
	// NOTE: these are well-defined functions, safe to evaluate

	// elasticity
	if(ctx->A_els) inv_eta_els = 2.0*ctx->A_els;
	// diffusion
	if(ctx->A_dif) inv_eta_dif = 2.0*ctx->A_dif;
	// dislocation
	if(ctx->A_dis) inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(ctx->DII, 1.0 - 1.0/ctx->N_dis);
	// Peierls
	if(ctx->A_prl) inv_eta_prl = 2.0*pow(ctx->A_prl, 1.0/ctx->N_prl)*pow(ctx->DII, 1.0 - 1.0/ctx->N_prl);

	//=================================
	// bottom end of viscosity interval
	//=================================

	// find bottom viscosity estimation (quasi-harmonic mean viscosity)
	// this is also a closed-form solution for simple cases

	inv_eta_bot = inv_eta_els + inv_eta_dif + inv_eta_dis + inv_eta_prl;

	// truncate bottom viscosity to maximum
	if(inv_eta_bot < inv_eta_max)
	{
		// bottom limit is larger than maximum (no iterations due to [A])
		(*eta) = eta_max; PetscFunctionReturn(0);
	}

	// truncate bottom viscosity to minimum
	if(inv_eta_bot > inv_eta_min)
	{
		// bottom limit is smaller than minimum (correcting)
		inv_eta_bot = inv_eta_min;
	}

	// initialize bottom end of solution interval
	eta_bot = 1.0/inv_eta_bot;

	// never iterate if closed-form solution is available
	if(ctx->cfsol == PETSC_TRUE)
	{
		(*eta) = eta_bot; PetscFunctionReturn(0);
	}

	// check for truncation from the bottom
	if(GetConsEqRes(eta_bot, ctx) < 0.0)
	{
		// no iterations due to [C]
		(*eta) = eta_bot; PetscFunctionReturn(0);
	}

	//==============================
	// top end of viscosity interval
	//==============================

	// initialize minimum isolated viscosity
	inv_eta_top = 0.0;

	// find top viscosity estimation (minimum isolated viscosity)
	if(inv_eta_top < inv_eta_els) inv_eta_top = inv_eta_els;
	if(inv_eta_top < inv_eta_dif) inv_eta_top = inv_eta_dif;
	if(inv_eta_top < inv_eta_dis) inv_eta_top = inv_eta_dis;
	if(inv_eta_top < inv_eta_prl) inv_eta_top = inv_eta_prl;

	// truncate top viscosity to minimum
	if(inv_eta_top > inv_eta_min)
	{
		// top limit is smaller than minimum (no iterations due to [B])
		(*eta) = eta_min; PetscFunctionReturn(0);
	}

	// truncate top viscosity to maximum
	if(inv_eta_top < inv_eta_max)
	{
		// top limit is larger than maximum (correcting)
		inv_eta_top = inv_eta_max;
	}

	// initialize top end of solution interval
	eta_top = 1.0/inv_eta_top;

	// check for truncation from the top
	if(GetConsEqRes(eta_top, ctx) > 0.0)
	{
		// no iterations due to [C]
		(*eta) = eta_top; PetscFunctionReturn(0);
	}

	//====================
	// I T E R A T I O N S
	//====================

	// set initial guess, interval & tolerances
	B  = eta_bot;
	C  = eta_top;
	R  = 0.5*(eta_bot + eta_top); // bisection rule
	RE = lim->eta_rtol;
	AE = lim->eta_atol;

	// solve constitutive equation in residual form
	DFZERO(&GetConsEqRes, ctx, &B, &C, R, RE, AE, &IFLAG);

	// check return code
	if(IFLAG != 1)

	// assign converged viscosity
	(*eta) = B;

	PetscFunctionReturn(0);
*/
}
//---------------------------------------------------------------------------
PetscScalar ApplyStrainSoft(Soft_t *sl, PetscScalar APS, PetscScalar par)
{
	// Apply strain softening to a parameter (friction, cohesion)
	PetscScalar k;
	if(!sl) return par;
	// compute scaling ratio
	if(APS <= sl->APS1)
		k = 1.0;
	if(APS > sl->APS1 && APS < sl->APS2)
		k = 1.0 - sl->A*((APS - sl->APS1)/(sl->APS2 - sl->APS1));
	if(APS >= sl->APS2)
		k = 1.0 - sl->A;
	// apply strain softening
	return par*k;
}
//---------------------------------------------------------------------------
// compute inverse deviatoric elastic viscosity
PetscScalar GetI2Gdt(
	PetscInt     numPhases,
	Material_t  *phases,
	PetscScalar *phRat,
	PetscScalar  dt)
{
	PetscInt    i;
	PetscScalar I2Gdt;

	// initialize inverse deviatoric elastic viscosity
	I2Gdt = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{	// average elastic materials only
		if(phases[i].G) I2Gdt += 0.5*phRat[i]/phases[i].G/dt;
	}
	return I2Gdt;
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DevConstEq"
PetscErrorCode DevConstEq(
	SolVarDev   *svDev,     // solution variables
	PetscInt     numPhases, // number phases
	Material_t  *phases,    // phase parameters
	PetscScalar *phRat,     // phase ratios
	MatParLim   *lim,       // phase parameters limits
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T)         // temperature
{
	// Evaluate deviatoric constitutive equations in control volume

	PetscInt     i;
	ConstEqCtx   ctx;
	Material_t  *mat;
	PetscScalar  DII, APS, eta, DIIpl;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize total effective strain rate & APS
	DII = svDev->DII;
	APS = svDev->APS;

	// initialize effective viscosity & plastic strain-rate
	svDev->eta   = 0.0;
	svDev->DIIpl = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// get reference to material parameters table
			mat = &phases[i];

			// setup nonlinear constitutive equation evaluation context
			ierr = ConstEqCtxSetup(&ctx, mat, lim, DII, APS, dt, p, T); CHKERRQ(ierr);

			// solve effective viscosity & plastic strain rate
			ierr = GetEffVisc(&ctx, lim, &eta, &DIIpl); CHKERRQ(ierr);

			//=============================
			// ADD GEOMETRIC AVERAGING HERE
			//=============================
			svDev->eta   += phRat[i]*eta;
			svDev->DIIpl += phRat[i]*DIIpl;

//			svDev->eta   += phRat[i]*log(eta);
//			svDev->DIIpl += phRat[i]*log(DIIpl);

		}
	}

//	svDev->eta   = exp(svDev->eta);
//	svDev->DIIpl = exp(svDev->DIIpl);

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
	MatParLim   *lim,       // phase parameters limits
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T)         // temperature
{
	// Evaluate volumetric constitutive equations in control volume

	PetscInt     i;
	Material_t  *mat;
	PetscScalar  cf_comp, cf_therm, IKdt;

	PetscFunctionBegin;

	// initialize effective density, thermal expansion & inverse bulk elastic viscosity
	svBulk->rho   = 0.0;
	svBulk->alpha = 0.0;
	svBulk->IKdt  = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{
		// update present phases only
		if(phRat[i])
		{
			// get reference to material parameters table
			mat = &phases[i];

			// initilaize
			IKdt     = 0.0;
			cf_comp  = 1.0;
			cf_therm = 1.0;

			// elastic compressiblility correction (Murnaghan's equation)
			// ro/ro_0 = (1 + K'*P/K)^(1/K')
			if(mat->K)
			{
				IKdt = 1.0/mat->K/dt;
				if(mat->Kp) cf_comp = pow(1.0 + mat->Kp*(p/mat->K), 1.0/mat->Kp);
				else        cf_comp = 1.0 + p/mat->K;
			}

			// thermal expansion correction
			// ro/ro_0 = 1 - alpha*(T - TRef)
			if(mat->alpha)
			{
				cf_therm  = 1.0 - mat->alpha*(T - lim->TRef);
			}

			// update density, thermal expansion & inverse bulk elastic viscosity
			svBulk->rho   += phRat[i]*mat->rho*cf_comp*cf_therm;
			svBulk->alpha += phRat[i]*mat->alpha;
			svBulk->IKdt  += phRat[i]*IKdt;

		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetStressCell"
PetscErrorCode GetStressCell(
	SolVarCell  *svCell, // solution variables
	MatParLim   *lim,    // phase parameters limits
	PetscScalar  dxx,    // effective normal strain rate components
	PetscScalar  dyy,    // ...
	PetscScalar  dzz)    // ...
{
	// compute stress, plastic strain-rate and shear heating term on cell

	SolVarDev   *svDev;
	PetscScalar  DII, cfpl, txx, tyy, tzz;

	// access deviatoric variables
	svDev = &svCell->svDev;

	// compute deviatoric stresses
	svCell->sxx = 2.0*svDev->eta*dxx;
	svCell->syy = 2.0*svDev->eta*dyy;
	svCell->szz = 2.0*svDev->eta*dzz;

	// get strain-rate invariant
	DII = svDev->DII;

	// use reference strain-rate instead of zero
	if(DII == 0.0) DII = lim->DII_ref;

	// compute plastic scaling coefficient
	cfpl = svDev->DIIpl/DII;

	// compute plastic strain-rate components
	txx = cfpl*dxx;
	tyy = cfpl*dyy;
	tzz = cfpl*dzz;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = 0.5*(txx*txx + tyy*tyy + tzz*tzz);

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	txx = svCell->dxx - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tyy = svCell->dyy - svDev->I2Gdt*(svCell->sxx - svCell->hxx);
	tzz = svCell->dzz - svDev->I2Gdt*(svCell->sxx - svCell->hxx);

	// compute shear heating term contribution
	svDev->Hr = (txx*svCell->sxx + tyy*svCell->syy + tzz*svCell->szz)*lim->shearHeatEff;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// compute stress, plastic strain-rate and shear heating term on edge
#undef __FUNCT__
#define __FUNCT__ "GetStressEdge"
PetscErrorCode GetStressEdge(
	SolVarEdge  *svEdge, // solution variables
	MatParLim   *lim,    // phase parameters limits
	PetscScalar  d)      // effective shear strain rate component
{

	SolVarDev   *svDev;
	PetscScalar  DII, cfpl, t;

	// access deviatoric variables
	svDev = &svEdge->svDev;

	// compute shear stress
	svEdge->s = 2.0*svDev->eta*d;

	// get strain-rate invariant
	DII = svDev->DII;

	// use reference strain-rate instead of zero
	if(DII == 0.0) DII = lim->DII_ref;

	// compute plastic scaling coefficient
	cfpl = svDev->DIIpl/DII;

	// compute plastic strain-rate components
	t = cfpl*d;

	// store contribution to the second invariant of plastic strain-rate
	svDev->PSR = t*t;

	// compute dissipative part of total strain rate (viscous + plastic = total - elastic)
	t = svEdge->d - svDev->I2Gdt*(svEdge->s - svEdge->h);

	// compute shear heating term contribution
	svDev->Hr = 2.0*t*svEdge->s*lim->shearHeatEff;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Elastic stress rotation functions
//---------------------------------------------------------------------------
void GetRotationMatrix(
	Tensor2RN   *R,  // rotation matrix
	PetscScalar  dt, // time step
	PetscScalar  wx, // vorticity vector components
	PetscScalar  wy, // ...
	PetscScalar  wz) // ...
{
	// compute rotation matrix from axis & angle (Euler-Rodrigues formula)
	// WARNING! Courant criterion for rotation angle should be implemented
	// angle tolerance probably should be also larger than machine epsilon

	PetscScalar w, theta, ct, st, cf;

	// get length of the instantaneous rotation axis (vorticity intensity)
	w = sqrt(wx*wx + wy*wy + wz*wz);

	// get finite rotation angle (vorticity intensity is twice the average angular velocity)
	theta = dt*(w/2.0);

	// round-off
	if(theta <  2.0*DBL_EPSILON)
	{
		R->xx = 1.0;   R->xy = 0.0;   R->xz = 0.0;
		R->yx = 0.0;   R->yy = 1.0;   R->yz = 0.0;
		R->zx = 0.0;   R->zy = 0.0;   R->zz = 1.0;

		return;
	}

	// get unit rotation axis
	wx /= w;
	wy /= w;
	wz /= w;

	// compute rotation operator using Euler-Rodrigues formula
	ct = cos(theta);
	st = sin(theta);
	cf = 1.0 - ct;

	R->xx =  ct    + cf*wx*wx;   R->xy = -st*wz + cf*wx*wy;   R->xz =  st*wy + cf*wx*wz;
	R->yx =  st*wz + cf*wy*wx;   R->yy =  ct    + cf*wy*wy;   R->yz = -st*wx + cf*wy*wz;
	R->zx = -st*wy + cf*wz*wx;   R->zy =  st*wx + cf*wz*wy;   R->zz =  ct    + cf*wz*wz;
}
//---------------------------------------------------------------------------
void RotateStress(Tensor2RN *R, Tensor2RS *S, Tensor2RS *SR)
{
	// rotate stress tensor
	// [SR] = [R] * [S] * [R]^t
	// SRij = Rik * Skl * Rjl

	PetscScalar dx, dy, dz;

	dx = R->xx*S->xx + R->xy*S->xy + R->xz*S->xz;
	dy = R->xx*S->xy + R->xy*S->yy + R->xz*S->yz;
	dz = R->xx*S->xz + R->xy*S->yz + R->xz*S->zz;

	SR->xx = dx*R->xx + dy*R->xy + dz*R->xz;
	SR->xy = dx*R->yx + dy*R->yy + dz*R->yz;
	SR->xz = dx*R->zx + dy*R->zy + dz*R->zz;

	dx = R->yx*S->xx + R->yy*S->xy + R->yz*S->xz;
	dy = R->yx*S->xy + R->yy*S->yy + R->yz*S->yz;
	dz = R->yx*S->xz + R->yy*S->yz + R->yz*S->zz;

	SR->yy = dx*R->yx + dy*R->yy + dz*R->yz;
	SR->yz = dx*R->zx + dy*R->zy + dz*R->zz;

	dx = R->zx*S->xx + R->zy*S->xy + R->zz*S->xz;
	dy = R->zx*S->xy + R->zy*S->yy + R->zz*S->yz;
	dz = R->zx*S->xz + R->zy*S->yz + R->zz*S->zz;

	SR->zz = dx*R->zx + dy*R->zy + dz*R->zz;
}
//---------------------------------------------------------------------------
void Tensor2RSCopy(Tensor2RS *A, Tensor2RS *B)
{
	// copy symmetric second order tensor B = A
	B->xx = A->xx;
	B->xy = A->xy; B->yy = A->yy;
    B->xz = A->xz; B->yz = A->yz; B->zz = A->zz;
}
//---------------------------------------------------------------------------


// ERROR HANDLING FOR EVALUATION CONTEXT ROUTINE

/*

 #include <math.h>
#if defined(math_errhandling) \
  && (math_errhandling & MATH_ERREXCEPT)
#include <fenv.h>
#endif



#if defined(math_errhandling) \
  && (math_errhandling & MATH_ERREXCEPT)
  feclearexcept(FE_ALL_EXCEPT);
#endif
errno = 0;

// call the function

double x;
double y;
double result;

if (((x == 0.f) && islessequal(y, 0)) || (isless(x, 0))) {
  // handle domain error
}

result = pow(x, y);

#if !defined(math_errhandling) \
  || (math_errhandling & MATH_ERRNO)
if (errno != 0) {
  // handle range error
}
#endif
#if defined(math_errhandling) \
  && (math_errhandling & MATH_ERREXCEPT)
if (fetestexcept(FE_INVALID
               | FE_DIVBYZERO
               | FE_OVERFLOW
               | FE_UNDERFLOW) != 0)
{
  // handle range error
}
#endif

*/

//---------------------------------------------------------------------------

// ELASTIC VISCOSITY COMPUTATION WITH RESCALING

/*
	PetscInt    i;
	PetscScalar sum = 0.0;

	// initialize inverse deviatoric elastic viscosity
	svDev->I2Gdt = 0.0;

	// scan all phases
	for(i = 0; i < numPhases; i++)
	{	// average elastic materials only
		if(phases[i].G)
		{	svDev->I2Gdt += 0.5*svDev->phRat[i]/phases[i].G/dt;
			sum += svDev->phRat[i];
		}
	}

	// re-scale for purely inelastic materials
	if(sum) svDev->I2Gdt /= sum;
*/

//---------------------------------------------------------------------------

// THERMAL CONSTITUTIVE EQUATIONS EXAMPLE

/*
//---------------------------------------------------------------------------
void ConstEq::heat(
			  Material  & mat,   // material properties
              double      T,     // updated temperature   (element average)
              double      Tn,    // converged temperature (element average)
              double    * grad,  // temperature gradient  (element average)
              double      step,  // time step
              double    * q,     // heat flux vector
              double    & Ud)    // approximation for internal energy rate
{
	// compute approximation for internal energy rate
	Ud = mat.Cp*(T - Tn)/step;
	// compute heat flux vector
	q[0] = - mat.k*grad[0];
	q[1] = - mat.k*grad[1];
	q[2] = - mat.k*grad[2];
}
//---------------------------------------------------------------------------
void ConstEq::heat(Material  & mat,   // material properties
              double    * grad,  // temperature gradient  (element average)
              double    * q)     // heat flux vector
{
	// compute heat flux vector
	q[0] = - mat.k*grad[0];
	q[1] = - mat.k*grad[1];
	q[2] = - mat.k*grad[2];
}
//---------------------------------------------------------------------------

*/
