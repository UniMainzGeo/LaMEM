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
#include "solVar.h"
#include "constEq.h"
#include "dfzero.h"
#include "tools.h"

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

	if(T) RT =  lim->Rugc*T;
	else  RT = -1.0;

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

	// LINEAR DIFFUSION CREEP (NEWTONIAN)
	if(mat->Bd)
	{
		Q          = (mat->Ed + p*mat->Vd)/RT;
		ctx->A_dif =  mat->Bd*exp(-Q);
		ln++;
	}

	// DISLOCATION CREEP (POWER LAW)
	if(mat->Bn)
	{
		Q          = (mat->En + p*mat->Vn)/RT;
		ctx->N_dis =  mat->n;
		ctx->A_dis =  mat->Bn*exp(-Q);
		nl++;
	}

	// PEIERLS CREEP (LOW TEMPERATURE RATE-DEPENDENT PLASTICITY, POWER-LAW APPROXIMATION)
	// ONLY EVALUATE FOR TEMPERATURE-DEPENDENT CASES
	if(mat->Bp && T)
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
	PetscScalar *eta_total,
	PetscScalar *eta_creep,
	PetscScalar *DIIpl)
{
	// stabilization parameters
	PetscScalar cf_eta_min = 10.0;
	PetscScalar eta_ve, DIIve, H;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_dis, inv_eta_prl, inv_eta_creep;

	PetscFunctionBegin;

	//==============
	// INITIAL GUESS
	//==============

	// zero out plastic strain rate
	(*DIIpl) = 0.0;

	// set reference viscosity as initial guess
	if(lim->eta_ref && lim->initGuessFlg == PETSC_TRUE)
	{
		(*eta_total) = lim->eta_ref;
		(*eta_creep) = lim->eta_ref;

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
	// NONLINEAR CREEP
	//================

	// get creep viscosity (old and good quasi-harmonic mean)
	inv_eta_creep = inv_eta_dif + inv_eta_dis + inv_eta_prl;

	// enforce limits
	if(inv_eta_creep < 1.0/lim->eta_max) inv_eta_creep = 1.0/lim->eta_max;
	if(inv_eta_creep > 1.0/lim->eta_min) inv_eta_creep = 1.0/lim->eta_min;

	// store creep viscosity for output
	(*eta_creep) = 1.0/inv_eta_creep;

	//===========
	// ELASTICITY
	//===========

	// compute visco-elastic viscosity
	eta_ve = 1.0/(inv_eta_els + inv_eta_creep);

	// visco-elastic prediction
	(*eta_total) = eta_ve;

	//===========
	// PLASTICITY
	//===========

	if(ctx->taupl && lim->initGuessFlg != PETSC_TRUE)
	{
		// compute visco-elastic strain rate
		DIIve = ctx->taupl/(2.0*eta_ve);

		if(lim->quasiHarmAvg == PETSC_TRUE)
		{
			//====================================
			// regularized rate-dependent approach
			//====================================

			// check for nonzero plastic strain rate
			if(DIIve < ctx->DII)
			{
				// store plastic strain rate & viscosity
				H            = eta_ve/cf_eta_min;
				(*eta_total) = 1.0/(1.0/eta_ve + 1.0/H) + (ctx->taupl/(2.0*ctx->DII))/(1.0 + H/eta_ve);
				(*DIIpl)     = ctx->DII*(1.0 - (*eta_total)/eta_ve);
			}
		}
		else
		{
			//====================================
			// classical rate-independent approach
			//====================================

			// check for nonzero plastic strain rate
			if(DIIve < ctx->DII)
			{

				// store plastic strain rate & viscosity
				(*eta_total)  = ctx->taupl/(2.0*ctx->DII);
				(*DIIpl)      = ctx->DII - DIIve;

				// lower viscosity bound
				if((*eta_total) < lim->eta_min)
				{
					(*eta_total) = lim->eta_min;
					(*DIIpl)     = ctx->DII*(1.0 - lim->eta_min/eta_ve);
				}

			}
		}
	}

	PetscFunctionReturn(0);
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
	PetscScalar *eta_creep, // creep viscosity (for output)
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
	PetscScalar  DII, APS, eta_total, eta_creep_phase, DIIpl;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize total effective strain rate & APS
	DII = svDev->DII;
	APS = svDev->APS;

	// initialize effective viscosity & plastic strain-rate
	svDev->eta   = 0.0;
	svDev->DIIpl = 0.0;
	(*eta_creep) = 0.0;

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
			ierr = GetEffVisc(&ctx, lim, &eta_total, &eta_creep_phase, &DIIpl); CHKERRQ(ierr);

			//=============================
			// ADD GEOMETRIC AVERAGING HERE
			//=============================
			svDev->eta   += phRat[i]*eta_total;
			svDev->DIIpl += phRat[i]*DIIpl;
			(*eta_creep) += phRat[i]*eta_creep_phase;

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
	PetscScalar  depth,     // depth for depth-dependent density model
	PetscScalar  dt,        // time step
	PetscScalar  p,         // pressure
	PetscScalar  T)         // temperature
{
	// Evaluate volumetric constitutive equations in control volume

	PetscInt     i;
	Material_t  *mat;
	PetscScalar  cf_comp, cf_therm, IKdt, rho;

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
				cf_therm  = 1.0 - mat->alpha*(T - lim->TRef);
			}

			// get density
			if(mat->rho_n)
			{
				// depth-dependent density (ad-hoc)
				rho = mat->rho - (mat->rho - lim->rho_fluid)*mat->rho_n*exp(-mat->rho_c*depth);
			}
			else
			{
				// temperature & pressure-dependent density
				rho = mat->rho*cf_comp*cf_therm;
			}

			// update density, thermal expansion & inverse bulk elastic viscosity
			svBulk->rho   += phRat[i]*rho;
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

	PetscFunctionBegin;

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

	PetscFunctionBegin;

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
		Tensor2RNUnit(R);

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
// Temperature parameters functions
//---------------------------------------------------------------------------
void GetTempParam(
	PetscInt     numPhases,
	Material_t  *phases,
	PetscScalar *phRat,
	PetscScalar *k_,  // conductivity
	PetscScalar *Cp_, // capacity
	PetscScalar *A_)  // radiogenic heat
{
	// compute effective energy parameters in the cell

	PetscInt    i;
    Material_t  *M;
	PetscScalar cf, k, Cp, A;

	// initialize
	k  = 0.0;
    Cp = 0.0;
	A  = 0.0;

	// average all phases
	for(i = 0; i < numPhases; i++)
	{
		M   = &phases[i];
		cf  = phRat[i];
		k  += cf*M->k;
		Cp += cf*M->Cp;
		A  += cf*M->A;
	}

	// store
	if(k_)  (*k_)  = k;
    if(Cp_) (*Cp_) = Cp;
	if(A_)  (*A_)  = A;
}
//---------------------------------------------------------------------------
// Infinite Strain Axis (ISA) calculation functions

// Kaminski et. al, 2004. D-Rex, a program for calculation of seismic
// anisotropy due to crystal lattice preferred orientation in the convective
// upper mantle, Gophys. J. Int, 158, 744-752.

//---------------------------------------------------------------------------
PetscInt Tensor2RNEigen(Tensor2RN *L, PetscScalar tol, PetscScalar eval[])
{
	//=======================================================================
	//
	// compute eigenvalues of a nonsymmetric tensor with zero trace
	//
	// WARNING! TENSOR TRACE MUST BE ZERO (NOT CHECKED HERE)
	//
	// return codes:
	//
	//    0 - three nearly zero eigenvalues (up to a tolerance)
	//    1 - three real eigenvalues (can be multiple)
	//    2 - one positive real eigenvalue & a complex conjugate pair
	//    3 - one negative real eigenvalue & a complex conjugate pair
	//
	// three real eigenvalues are sorted in descending order
	// real eigenvalue always precedes the complex conjugate pair
	//
	//=======================================================================

	PetscInt    code;
	PetscScalar I2, I3, p, q, D, theta, l1, l2, l3, cx, t, sd, r, s;

	// get invariants
	I2 = L->xx*L->yy + L->yy*L->zz + L->xx*L->zz
	-    L->xy*L->yx - L->yz*L->zy - L->xz*L->zx;

	I3 = L->xx*(L->yy*L->zz - L->yz*L->zy)
	+    L->xy*(L->yz*L->zx - L->yx*L->zz)
	+    L->xz*(L->yx*L->zy - L->yy*L->zx);

	// get discriminant
	p =  I2;
	q = -I3;
	D = (q*q)/4.0 + (p*p*p)/27.0;

	if(fabs(D) < tol)
	{
		//================================
		// three (nearly) zero eigenvalues
		//================================

		l1   = 0.0;
		l2   = 0.0;
		l3   = 0.0;
		cx   = 0.0;

		// set return code
		code = 0;
	}
	else if(D < 0.0)
	{
		//=======================
		// three real eigenvalues
		//=======================

		theta = ARCCOS((3.0*q)/(2.0*p)*sqrt(-3.0/p));

		l1   = 2.0*sqrt(-p/3.0)*cos( theta            /3.0);
		l2   = 2.0*sqrt(-p/3.0)*cos((theta - 2.0*M_PI)/3.0);
		l3   = 2.0*sqrt(-p/3.0)*cos((theta - 4.0*M_PI)/3.0);
		cx   = 0.0;

		// set return code
		code = 1;

		// sort eigenvalues
		if(l2 > l1) { t = l1; l1 = l2; l2 = t; }
		if(l3 > l1) { t = l1; l1 = l3; l3 = t; }
		if(l3 > l2) { t = l2; l2 = l3; l3 = t; }
	}
	else
	{
		//=============================================
		// one real eigenvalue & complex conjugate pair
		//=============================================

		sd = sqrt(D);

		r  = ODDROOT(-q/2.0 + sd, 1.0/3.0);
		s  = ODDROOT(-q/2.0 - sd, 1.0/3.0);

		// get real parts of eigenvalues
		l1 =  r + s;
		l2 = -l1/2.0;
		l3 =  l2;

		// get modulus of imaginary part of complex conjugate pair
		cx = fabs(r-s)*sqrt(3.0)/2.0;

		// set return code
		if(l1 > 0.0) code = 2; // positive real root
		else         code = 3; // negative real root

		// complex eigenvalues are:
		// l2 = -(r+s)/2 + cx*i
		// l3 = -(r+s)/2 - cx*i
	}

	// store result
	eval[0] = l1;
	eval[1] = l2;
	eval[2] = l3;
	eval[3] = cx;

	return code;
}
//---------------------------------------------------------------------------
PetscInt getISA(Tensor2RN *pL, PetscScalar ISA[], PetscScalar *plnrm)
{
	// compute direction of Infinite Strain Axis

	// return codes:
	//    -2 - spectral decomposition failed to converge
	//    -1 - ISA is undefined
	//     0 - ISA is defined, computed, and returned
	//     1 - simple shear case (ISA has same direction as velocity)

	Tensor2RS   Cs;
	Tensor2RN   L, L2, I, F, Ft, C;
	PetscScalar l1, l2, l3, cx, D, lnrm, eval[4], evect[9], ltol, ttol;
	PetscInt    maxit, code;

	// WARNING! set tolerances via command line using MatParLim structure
	ltol  = 1e-9;  // loose tolerance
	ttol  = 1e-13; // tight tolerance
	maxit = 30;    // maximum number of Jacobi rotations

	// initialize
	ISA[0] = 0.0;
	ISA[1] = 0.0;
	ISA[2] = 0.0;

	// copy velocity gradient, compute norm
	Tensor2RNCopy(pL, &L);
	Tensor2RNNorm(&L, &lnrm);

	// return norm of the velocity gradient if necessary
	if(plnrm) (*plnrm) = lnrm;

	//========================================
	// *** zero gradient, ISA is undefined ***
	//========================================
	if(!lnrm) return -1;

	// normalize velocity gradient, remove trace
	Tensor2RNDivide(&L, lnrm);
	Tensor2RNTrace(&L);

	// get eigenvalues of the velocity gradient
	code = Tensor2RNEigen(&L, ltol, eval);

	//==================================================
	// *** three zero eigenvalues, simple shear case ***
	//==================================================
	if(code == 0) return 1;

	//===================================================================
	// *** negative real + complex pair eigenvalues, ISA is undefined ***
	//===================================================================
	if(code == 3) return -1;

	// get denominator of Sylvester's formula
	l1 = eval[0];
	l2 = eval[1];
	l3 = eval[2];
	cx = eval[3];
	D  = (l1 - l2)*(l1 - l3) + cx*cx;

	//===============================================
	// *** multiple eigenvalues, ISA is undefined ***
	//===============================================
	if(fabs(D) < ttol) return -1;

	// three distinct real or one positive real + complex pair eigenvalues
	// ISA is defined by l2 & l3 eigenvalues
	// compute deformation gradient
	// scaling doesn't affect eigenvectors (denominator is set to unit)
	// F = (L - l2*I)*(L - l3*I)/D

	Tensor2RNUnit(&I);
	Tensor2RNProduct(&L, &L, &L2);
	Tensor2RNSum3(&L2, 1.0/D, &L, -(l2 + l3)/D, &I, (l2*l3 + cx*cx)/D, &F);

	// compute right Cauchy-Green deformation tensor C = F^t*F
	Tensor2RNTranspose(&F, &Ft);
	Tensor2RNProduct(&Ft, &F, &C);
	Tensor2RNCopySym(&C, &Cs);

	// perform spectral decomposition
	code = Tensor2RSSpectral(&Cs, eval, evect, ttol, ltol, maxit);

	//=====================================================================
	// *** spectral decomposition failed to converge, ISA is undefined  ***
	//=====================================================================
	if(code) return -2;

	// ISA is the eigenvector corresponding to the largest eigenvalue
	ISA[0] = evect[0];
	ISA[1] = evect[1];
	ISA[2] = evect[2];

	//===============================================
	// *** ISA is defined, computed, and returned ***
	//===============================================
	return 0;
}
//---------------------------------------------------------------------------
void Tensor2RNNorm(Tensor2RN *A, PetscScalar *pk)
{
	// k = |A|

	PetscScalar s, k;

	s = fabs(A->xx) + fabs(A->xy) + fabs(A->xz);           k = s;
	s = fabs(A->yx) + fabs(A->yy) + fabs(A->yz); if(s > k) k = s;
	s = fabs(A->zx) + fabs(A->zy) + fabs(A->zz); if(s > k) k = s;

	(*pk) = k;
}
//---------------------------------------------------------------------------
void Tensor2RNDivide(Tensor2RN *A, PetscScalar k)
{
	// A = A/k

	A->xx /= k; A->xy /= k; A->xz /= k;
	A->yx /= k; A->yy /= k; A->yz /= k;
	A->zx /= k; A->zy /= k; A->zz /= k;
}
//---------------------------------------------------------------------------
void Tensor2RNTrace(Tensor2RN *A)
{
	// A = A - tr[A]/3*I

	PetscScalar tr;

	tr = (A->xx + A->yy + A->zz)/3.0;

	A->xx -= tr;
	A->yy -= tr;
	A->zz -= tr;
}
//---------------------------------------------------------------------------
void Tensor2RNSym(Tensor2RN *A, Tensor2RN *B)
{
	// B = (A + A^t)/2

	B->xx =  A->xx;                B->xy = (A->yx + A->xy)/2.0;   B->xz = (A->zx + A->xz)/2.0;
	B->yx = (A->xy + A->yx)/2.0;   B->yy =  A->yy;                B->yz = (A->zy + A->yz)/2.0;
	B->zx = (A->xz + A->zx)/2.0;   B->zy = (A->yz + A->zy)/2.0;   B->zz =  A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNProduct(Tensor2RN *A, Tensor2RN *B, Tensor2RN *C)
{
	// C = A*B

	C->xx = A->xx*B->xx + A->xy*B->yx + A->xz*B->zx;
	C->xy = A->xx*B->xy + A->xy*B->yy + A->xz*B->zy;
	C->xz = A->xx*B->xz + A->xy*B->yz + A->xz*B->zz;
	C->yx = A->yx*B->xx + A->yy*B->yx + A->yz*B->zx;
	C->yy = A->yx*B->xy + A->yy*B->yy + A->yz*B->zy;
	C->yz = A->yx*B->xz + A->yy*B->yz + A->yz*B->zz;
	C->zx = A->zx*B->xx + A->zy*B->yx + A->zz*B->zx;
	C->zy = A->zx*B->xy + A->zy*B->yy + A->zz*B->zy;
	C->zz = A->zx*B->xz + A->zy*B->yz + A->zz*B->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNTranspose(Tensor2RN *A, Tensor2RN *B)
{
	// B = A^t

	B->xx = A->xx; B->xy = A->yx; B->xz = A->zx;
	B->yx = A->xy; B->yy = A->yy; B->yz = A->zy;
	B->zx = A->xz; B->zy = A->yz; B->zz = A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNCopy(Tensor2RN *A, Tensor2RN *B)
{
	// B = A

	B->xx = A->xx; B->xy = A->xy; B->xz = A->xz;
	B->yx = A->yx; B->yy = A->yy; B->yz = A->yz;
	B->zx = A->zx; B->zy = A->zy; B->zz = A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNCopySym(Tensor2RN *A, Tensor2RS *B)
{
	// B = A

	B->xx = A->xx;
	B->xy = A->xy; B->yy = A->yy;
	B->xz = A->xz; B->yz = A->yz; B->zz = A->zz;
}

//---------------------------------------------------------------------------
void Tensor2RNUnit(Tensor2RN *A)
{
	// A = I

	A->xx = 1.0; A->xy = 0.0; A->xz = 0.0;
	A->yx = 0.0; A->yy = 1.0; A->yz = 0.0;
	A->zx = 0.0; A->zy = 0.0; A->zz = 1.0;
}
//---------------------------------------------------------------------------
void Tensor2RNSum3(
	Tensor2RN *A, PetscScalar ka,
	Tensor2RN *B, PetscScalar kb,
	Tensor2RN *C, PetscScalar kc,
	Tensor2RN *R)
{
	// R = ka*A + kb*B + kc*C

	R->xx = ka*A->xx + kb*B->xx + kc*C->xx;
	R->xy = ka*A->xy + kb*B->xy + kc*C->xy;
	R->xz = ka*A->xz + kb*B->xz + kc*C->xz;
	R->yx = ka*A->yx + kb*B->yx + kc*C->yx;
	R->yy = ka*A->yy + kb*B->yy + kc*C->yy;
	R->yz = ka*A->yz + kb*B->yz + kc*C->yz;
	R->zx = ka*A->zx + kb*B->zx + kc*C->zx;
	R->zy = ka*A->zy + kb*B->zy + kc*C->zy;
	R->zz = ka*A->zz + kb*B->zz + kc*C->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNView(Tensor2RN *A, const char *msg)
{
	printf("%s: \n\n", msg);
	printf("%g %g %g \n",   A->xx, A->xy, A->xz);
	printf("%g %g %g \n",   A->yx, A->yy, A->yz);
	printf("%g %g %g \n\n", A->zx, A->zy, A->zz);
}
//---------------------------------------------------------------------------
void Tensor2RSView(Tensor2RS *A, const char *msg)
{
	printf("%s: \n\n", msg);
	printf("%g %g %g \n",   A->xx, A->xy, A->xz);
	printf("%g %g %g \n",   A->xy, A->yy, A->yz);
	printf("%g %g %g \n\n", A->xz, A->yz, A->zz);
}
//---------------------------------------------------------------------------
PetscInt Tensor2RSSpectral(
	Tensor2RS   *A,      // symmetric tensor
	PetscScalar eval[],  // eigenvalues (sorted)
	PetscScalar evect[], // eigenvectors (corresponding)
	PetscScalar ttol,    // tight tolerance (convergence condition)
	PetscScalar ltol,    // loose tolerance (divergence condition)
	PetscInt    itmax)   // maximum number rotations
{
	//=====================================================
	// Jacobi rotation algorithm for spectral decomposition
	//=====================================================

	// return codes:
	// 	 0 - converged to tight tolerance within maximum rotations
	// 	 1 - failed to converge to loose tolerance within maximum rotations

	PetscInt    iter, opt, code;
	PetscScalar atmp, ntmp[3];

	PetscScalar f, max, theta, t, c, s, tau, w, z;
	PetscScalar a1, a2, a3, a12, a13, a23, *n1, *n2, *n3;

	// macro for single Jacobi rotation
	#define JAC_ROT(pp, qq, pq, rp, rq, vp, vq) \
	{	theta = 0.5*(qq - pp)/pq; \
		t     = 1.0/(fabs(theta) + sqrt(theta*theta + 1.0)); \
		if(theta < 0.0) t = -t; \
		c = 1.0/sqrt(t*t + 1.0); s = t*c; tau = s/(1.0 + c); \
		pp -= t*pq;  qq += t*pq;   pq     = 0.0; \
		w  =  rp;     z  = rq;     rp    -= s*(z + tau*w);  rq    += s*(w - tau*z); \
		w  =  vp[0];  z  = vq[0];  vp[0] -= s*(z + tau*w);  vq[0] += s*(w - tau*z); \
		w  =  vp[1];  z  = vq[1];  vp[1] -= s*(z + tau*w);  vq[1] += s*(w - tau*z); \
		w  =  vp[2];  z  = vq[2];  vp[2] -= s*(z + tau*w);  vq[2] += s*(w - tau*z); \
	}

	// macro for swapping two principal values & principal directions
	#define SWAP_EIG_PAIR(ai, aj, ni, nj) \
	{	atmp    = ai;    ai    = aj;    aj    = atmp; \
		ntmp[0] = ni[0]; ni[0] = nj[0]; nj[0] = ntmp[0]; \
		ntmp[1] = ni[1]; ni[1] = nj[1]; nj[1] = ntmp[1]; \
		ntmp[2] = ni[2]; ni[2] = nj[2]; nj[2] = ntmp[2]; \
	}

	// set return code
	code = 0;

	// copy tensor components
	a1  = A->xx;
	a12 = A->xy; a2  = A->yy;
	a13 = A->xz; a23 = A->yz; a3 = A->zz;

	// initialize principal directions
	n1 = evect;
	n2 = evect + 3;
	n3 = evect + 6;

	n1[0] = 1.0; n2[0] = 0.0; n3[0] = 0.0;
	n1[1] = 0.0; n2[1] = 1.0; n3[1] = 0.0;
	n1[2] = 0.0; n2[2] = 0.0; n3[2] = 1.0;

	// zero out off-diagonal component by Jacobi rotations
	iter = 0;
	do
	{	// select maximum off-diagonal component
		f = fabs(a12);               max = f;  opt = 1;
		f = fabs(a13); if(f > max) { max = f;  opt = 2; }
		f = fabs(a23); if(f > max) { max = f;  opt = 3; }

		// check convergence
		if(max < ttol) break;

		// perform Jacobi rotation
		if     (opt == 1) JAC_ROT(a1, a2, a12, a13, a23, n1, n2) // a12 term
		else if(opt == 2) JAC_ROT(a1, a3, a13, a12, a23, n1, n3) // a13 term
		else              JAC_ROT(a2, a3, a23, a12, a13, n2, n3) // a23 term

	} while(++iter < itmax);

	// check divergence
	if(iter == itmax)
	{
		// select maximum off-diagonal component
		f = fabs(a12);               max = f;  opt = 1;
		f = fabs(a13); if(f > max) { max = f;  opt = 2; }
		f = fabs(a23); if(f > max) { max = f;  opt = 3; }

		// check whether algorithm failed to converge to loose
		// tolerance within prescribed number of iterations
		if(max > ltol) code = 1;
	}

	// sort principal values in descending order & permute principal directions
	if(a2 > a1) SWAP_EIG_PAIR(a1, a2, n1, n2)
	if(a3 > a1) SWAP_EIG_PAIR(a1, a3, n1, n3)
	if(a3 > a2) SWAP_EIG_PAIR(a2, a3, n2, n3)

	// store eigenvalues
	eval[0] = a1;
	eval[1] = a2;
	eval[2] = a3;

	return code;
}
//---------------------------------------------------------------------------
/*
// ERROR HANDLING FOR CONTEXT EVALUATION ROUTINE

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
/*
// ELASTIC VISCOSITY COMPUTATION WITH RESCALING
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
/*
// THERMAL CONSTITUTIVE EQUATIONS EXAMPLE
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
*/
//---------------------------------------------------------------------------
/*
{
	// ACHTUNG! THIS WILL ONLY WORK FOR:
	//  * LINEAR VISCO-ELASTO-PLASTIC MATERIAL
	//  * POWER-LAW MATERIAL

	// stabilization parameters
	PetscScalar cf_eta_min = 10.0;

	PetscScalar eta_ve, DIIve, H;
	PetscScalar inv_eta_els, inv_eta_dif, inv_eta_dis;

	PetscFunctionBegin;

	// zero out plastic strain rate
	(*DIIpl) = 0.0;

	// set reference viscosity as initial guess
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
			inv_eta_dis = 2.0*pow(ctx->A_dis, 1.0/ctx->N_dis)*pow(ctx->DII,     1.0 - 1.0/ctx->N_dis);
		}

		(*eta) = 1.0/inv_eta_dis;

	}
	else
	{
		// elasticity
		inv_eta_els = 2.0*ctx->A_els;

		// linear viscous creep
		inv_eta_dif = 2.0*ctx->A_dif;

		// compute visco-elastic viscosity
		eta_ve = 1.0/(inv_eta_els + inv_eta_dif);

		// visco-elastic prediction
		(*eta) = eta_ve;

		//===========
		// PLASTICITY
		//===========

		if(ctx->taupl && lim->initGuessFlg != PETSC_TRUE)
  		{
			// compute visco-elastic strain rate
			DIIve = ctx->taupl/(2.0*eta_ve);

			if(lim->quasiHarmAvg == PETSC_TRUE)
			{
				//====================================
				// regularized rate-dependent approach
				//====================================

				// check for nonzero plastic strain rate
				if(DIIve < ctx->DII)
				{
					// store plastic strain rate & viscosity
					H        = eta_ve/cf_eta_min;
					(*eta)   = 1.0/(1.0/eta_ve + 1.0/H) + (ctx->taupl/(2.0*ctx->DII))/(1.0 + H/eta_ve);
					(*DIIpl) = ctx->DII*(1.0 - (*eta)/eta_ve);
				}
			}
			else
			{
				//====================================
				// classical rate-independent approach
				//====================================

				// check for nonzero plastic strain rate
				if(DIIve < ctx->DII)
				{
					// store plastic strain rate & viscosity
					(*eta)   = ctx->taupl/(2.0*ctx->DII);
					(*DIIpl) = ctx->DII - DIIve;
				}
			}
		}
	}

	// enforce constraints
	if((*eta) < lim->eta_min) (*eta) = lim->eta_min;
	if((*eta) > lim->eta_max) (*eta) = lim->eta_max;

	PetscFunctionReturn(0);

	//=============================================


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
}
*/
//---------------------------------------------------------------------------
/*
	ABANDONED PLASTICITY STABILIZATIONS
	PetscScalar a          = (1.0 + sqrt(5.0))/2.0;
	PetscScalar cf_eta_max = 1000.0;
	PetscScalar cf_eta_min = 10.0;
	PetscScalar k          = log(cf_eta_max - 1.0)/(a - 1.0);

	// penalty function
	H = (eta_ve/cf_eta_min)*(1.0 + exp(k*(a - ctx->DII/DIIve)));

	// regularized viscosity
	(*eta)   = 1.0/(1.0/eta_ve + 1.0/H) + (ctx->taupl/(2.0*ctx->DII))/(1.0 + H/eta_ve);

	// plastic strain rate
	if((*eta) > eta_ve) (*DIIpl) = 0.0;
	else                (*DIIpl) = ctx->DII*(1.0 - (*eta)/eta_ve);

	// plasticity regularization parameters (WARNING: add input file options!)
	k = 7.0;
	a = (1.0 + sqrt(5.0))/2.0;
	k = 3.0;
	a = 1.0;

	// compute shifted approximate Heaviside function
	H = 1.0/(1.0 + exp(k*(a - ctx->DII/DIIve)));
	// compute regularized effective viscosity
	(*eta) = 1.0/((1.0-H)/eta_ve + H/eta_pl);
	(*eta) = 1.0/(1.0/eta_ve + 1.0/eta_pl);

 	if(ctx->taupl && lim->initGuessFlg != PETSC_TRUE)
  	{
		// compute visco-elastic strain rate
		DIIve = ctx->taupl/(2.0*eta_ve);

		// check for nonzero plastic strain rate
		if(DIIve < ctx->DII)
		{
			// store plastic strain rate
			(*DIIpl) = ctx->DII - DIIve;

			if(lim->quasiHarmAvg == PETSC_TRUE)
			{
				//====================================
				// regularized rate-dependent approach
				//====================================
				(*eta) = (ctx->taupl + 2.0*lim->eta_plast*(*DIIpl))/(2.0*ctx->DII);
			}
			else
			{
				//====================================
				// classical rate-independent approach
				//====================================

				(*eta) = ctx->taupl/(2.0*ctx->DII);
			}
		}
	}
*/
//---------------------------------------------------------------------------
/*
Tensor2RS A, SR;
Tensor2RN L, Lt, LLt, R, D;
L.xx =  1.0;     L.xy =  2.0;     L.xz =  3.0;
L.yx =  4.0;     L.yy =  5.0;     L.yz =  6.0;
L.zx =  7.0;     L.zy =  8.0;     L.zz =  9.0;
A.xx =  1.0000; A.xy = 0.5000; A.xz = 1.0/3.0;
A.yy =  1.0000; A.yz = 2.0/3.0;
A.zz = -2.0000;
Tensor2RNView(&L, "A");
Tensor2RNUnit(&R); Tensor2RNView(&R, "B = I");
Tensor2RNCopy(&L, &R); Tensor2RNView(&R, "B = A");
Tensor2RNDivide(&R, 2.0); Tensor2RNView(&R, "B = A/2");
Tensor2RNTranspose(&L, &Lt); Tensor2RNView(&Lt, "B = A^t");
Tensor2RNProduct(&L, &Lt, &LLt); Tensor2RNView(&LLt, "C = A*B");
Tensor2RNSum3(&L, 1.0, &Lt, 1.0, &LLt, 1.0, &R); Tensor2RNView(&R, "D = A + B + C");
Tensor2RSView(&A, "S");
RotateStress(&L, &A, &SR); Tensor2RSView(&SR, "S' = A*S*A^t");
Tensor2RNDeviator(&L, &D); Tensor2RNView(&D, "D = A - tr[A]/3*I");
Tensor2RSDeviator(&L, &SR); Tensor2RSView(&SR, "D = (A + A^t)/2 - tr[A]/3*I");



*/
