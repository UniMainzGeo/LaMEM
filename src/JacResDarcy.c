/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Boris Kaus, Anton Popov
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
 **    filename:   JacResDarcy.c
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
//......................   DARCY FUNCTIONS   ..........................
//---------------------------------------------------------------------------

/*
 * In this routine, we setup the functions required to solve the Darcy equation,
 * namely:
 *
 *
 * 		div(Kphi/mu grad(P_l)-rhol*g) =  Ss dP_l/dt + H
 *
 *	where:
 *		mu 		=   liquid viscosity
 *		rhol	=   liquid density
 *		Kphi 	=	permeability
 *		Ss		=	specific storage
 *		P_l 	=	liquid pressure
 *		g 		=	gravitational acceleration (acting in the z-direction )
 *		H		=   source
 *		Ts      =   Tensile Strength
 *
 */

#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "JacResDarcy.h"
#include "matrix.h"
#include "constEq.h"
#include "tools.h"
#include "interpolate.h"
#include "matProps.h"
#include "parsing.h"

//---------------------------------------------------------------------------

#define SCATTER_FIELD(da, vec, FIELD) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec)

// liquid density
#define GET_rhol \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, &rhol, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	rhol;

// liquid viscosity
#define GET_mul \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, &mul, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	mul;

// permeability
#define GET_Kphi \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, &Kphi, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Kphi;

// matrix compressibility
#define GET_betam \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, &betam, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	betam;

// liquid compressibility
#define GET_betal \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, &betal, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	betal;

// porosity
#define GET_Phi \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, NULL, &Phi, NULL, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Phi;

//// specific storage
//#define GET_Ssl \
//	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, &Ssl, NULL, NULL); CHKERRQ(ierr); \
//	buff[k][j][i] 	= 	Ssl;

// Tensile Strength
#define GET_Tsl \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, NULL, NULL, &Tsl, NULL, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Tsl;

// Undrained Poisson's ratio
#define GET_nuu \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &nuu, NULL, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	nuu;

// Undrained permeability
#define GET_Kphiu \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &Kphiu, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Kphiu;

// Undrained porosity
#define GET_Phiu \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &Phiu); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Phiu;

//---------------------------------------------------------------------------
// This extracts coefficients from named vectors attached to a da
#define ExtractCoefficientsFromDA(da, name,local, global, data)	\
	ierr = 	DMGetLocalVector(da,&local); 							CHKERRQ(ierr);	\
	ierr = 	DMGetNamedGlobalVector(da,name,&global); 				CHKERRQ(ierr);	\
	ierr = 	DMGlobalToLocalBegin(da, global,INSERT_VALUES, local); 	CHKERRQ(ierr);	\
	ierr = 	DMGlobalToLocalEnd  (da, global,INSERT_VALUES, local); 	CHKERRQ(ierr);	\
	ierr = 	DMDAVecGetArray(da, local, &data);  					CHKERRQ(ierr);	\

//---------------------------------------------------------------------------
// This restores the coefficients to a named vector in a DA
#define RestoreCoefficientsFromDA(da, name,local, global, data)	\
	ierr = 	DMDAVecRestoreArray(da, local, &data);  				CHKERRQ(ierr);	\
	ierr = 	DMRestoreNamedGlobalVector(da,name,&global); 			CHKERRQ(ierr);	\
	ierr = 	DMRestoreLocalVector(da,&local); 						CHKERRQ(ierr);	\


//---------------------------------------------------------------------------
// Darcy parameters functions
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetDarcyParam"
PetscErrorCode JacResGetDarcyParam(
		JacRes      *jr,
		PetscScalar *phRat,
		PetscScalar *rhol_,      // Liquid density
		PetscScalar *mul_,       // Liquid viscosity
		PetscScalar *Kphi_,      // Permeability
		PetscScalar *betam_,     // Matrix compressibility
		PetscScalar *betal_,     // Liquid compressibility
		PetscScalar *Phi_,       // Effective porosity
		PetscScalar *Ts_,        // Tensile strength
		PetscScalar *nuu_,       // Undrained poisson's ratio
		PetscScalar *Kphiu_,     // Undrained permeability
		PetscScalar *Phiu_)      // Undrained porosity
//PetscScalar *Ss_,      // Specific storage
{

	// compute effective darcy parameters in the cell

	PetscInt    i, numPhases;
    Material_t  *phases, *M;

	PetscScalar cf, rhol, mul, Kphi, betam, betal, Phi, Ts, nuu, Kphiu, Phiu; //Ss,

	PetscFunctionBegin;

	// initialize
	rhol    = 0.0;
	mul     = 0.0;
	Kphi    = 0.0;
	//Ss     = 0.0;
	betam   = 0.0;
	betal   = 0.0;
	Phi     = 0.0;
	Ts      = 0.0;
	nuu     = 0.0;
	Kphiu   = 0.0;
	Phiu    = 0.0;

	numPhases 	= 	jr->numPhases;
	phases    	= 	jr->phases;

	// average all phases
	for(i = 0; i < numPhases; i++)
	{
		M       = &phases[i];
		cf      =  phRat[i];

		// compute average permeability, liquid density and liquid viscosity
		rhol   += cf*M->rhol;
		mul    += cf*M->mul;
		Kphi   += cf*M->Kphi;
		//Ss   += cf*M->Ss;
		betam  += cf*M->betam;
		betal  += cf*M->betal;
		Phi    += cf*M->Phi;
		Ts     += cf*M->TS;
		nuu    += cf*M->nuu;
		Kphiu  += cf*M->Kphiu;
		Phiu   += cf*M->Phiu;
	}

	// store
	if(rhol_)   (*rhol_)  = rhol;
	if(mul_) 	(*mul_)   = mul;
	if(Kphi_) 	(*Kphi_)  = Kphi;
	//if(Ss_)   	(*Ss_)   = Ss;
	if(betam_) 	(*betam_) = betam;
	if(betal_) 	(*betal_) = betal;
	if(Phi_) 	(*Phi_)   = Phi;
	if(Ts_)     (*Ts_)    = Ts;
	if(nuu_) 	(*nuu_)   = nuu;
	if(Kphiu_) 	(*Kphiu_) = Kphiu;
	if(Phiu_) 	(*Phiu_)  = Phiu;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCheckDarcyParam"
PetscErrorCode JacResCheckDarcyParam(JacRes *jr)
{
	// check whether all Darcy material parameters are properly defined

	PetscInt    i, numPhases;
	Material_t  *phases, *M;

	PetscFunctionBegin;

	// Only for cases where Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	// initialize
	numPhases = jr->numPhases;
	phases    = jr->phases;

	// check all phases
	for(i = 0; i < numPhases; i++)
	{
		M = &phases[i];
		if(M->rhol  == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define liquid density of phase %lld\n",        (LLD)i);
		if(M->mul   == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define liquid viscosity of phase %lld\n",      (LLD)i);
		if(M->Kphi  == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define permeability of phase %lld\n",          (LLD)i);
		//if(M->Ss    == 0.0)    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define specific storage of phase %lld\n", 	(LLD)i);
		if(M->betal == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define liquid compressibility of phase %lld\n",(LLD)i);
		if(M->betam == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define matrix compressibility of phase %lld\n",(LLD)i);
		if(M->Phi   == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define porosity of phase %lld\n",              (LLD)i);
		if(M->TS    == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define tensile strength of phase %lld\n",      (LLD)i);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCreateDarcyParam"
PetscErrorCode JacResCreateDarcyParam(JacRes *jr)
{
	// setup temperature parameters

	FDSTAG *fs;
	const PetscInt *lx, *ly, *lz, dof = 1; //, dof=2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// only when Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create Darcy DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
		DMDA_STENCIL_BOX,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		dof, 1, lx, ly, lz, &jr->DA_Pl); CHKERRQ(ierr);

	// create darcy preconditioner matrix
	ierr = DMCreateMatrix(jr->DA_Pl, &jr->App); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// LiquidPressure solution vector
	ierr = DMCreateGlobalVector(jr->DA_Pl, &jr->dPl); 		CHKERRQ(ierr);

	// create local LiquidPressure vector using box-stencil central DMDA
	ierr = DMCreateLocalVector(jr->DA_Pl, &jr->lPl); 		CHKERRQ(ierr);

	// LiquidPressure residual (global)
	ierr = DMCreateGlobalVector(jr->DA_Pl, &jr->r_Pl); 		CHKERRQ(ierr);

	// create linear solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &jr->Pl_ksp); 		CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(jr->Pl_ksp,"darcy_");    	CHKERRQ(ierr);
	ierr = KSPSetFromOptions(jr->Pl_ksp);            		CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResDestroyDarcyParam"
PetscErrorCode JacResDestroyDarcyParam(JacRes *jr)
{
	// Destroy Darcy/Liquid Pressure parameters

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//ierr = VecDestroy(&jr->lPl);   	CHKERRQ(ierr);

	// only if Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	ierr = DMDestroy (&jr->DA_Pl); 			CHKERRQ(ierr);
	ierr = MatDestroy(&jr->App);   			CHKERRQ(ierr);
	ierr = VecDestroy(&jr->dPl);   			CHKERRQ(ierr);		// solution vector
	ierr = VecDestroy(&jr->lPl);			CHKERRQ(ierr);		// local solution vector
	ierr = KSPDestroy(&jr->Pl_ksp);			CHKERRQ(ierr);		// Darcy KSP context
	ierr = VecDestroy(&jr->r_Pl);   		CHKERRQ(ierr);		// residual vector

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResUpdateDarcy"
PetscErrorCode JacResUpdateDarcy(JacRes *jr)
{
	// correct LiquidPressure for diffusion (Newton update)

	PetscScalar ***sol_local, ***incr_global;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	ierr = DMDAVecGetArray(jr->DA_Pl,  jr->lPl, &sol_local ); 		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_Pl,  jr->dPl, &incr_global); 		CHKERRQ(ierr);

	ierr = DMDAGetCorners(jr->DA_Pl, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// update local solution vector
		sol_local[k][j][i] 	-=	incr_global[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(jr->DA_Pl,  jr->dPl, &incr_global); 		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_Pl,  jr->lPl, &sol_local); 		CHKERRQ(ierr);

	// apply two-point constraints
	ierr = JacResApplyDarcyBC(jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResApplyDarcyBC"
PetscErrorCode JacResApplyDarcyBC(JacRes *jr)
{
	// apply Darcy two-point constraints

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscScalar pmdof;
	PetscScalar	***sol;
	PetscScalar ***bcPl;
	PetscInt    mcx, mcy, mcz;
	PetscInt    I, J, K, fi, fj, fk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// exchange internal ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lPl)

	// access local solution & boundary constraints
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lPl,  &sol);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcPl, &bcPl); CHKERRQ(ierr);

	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = sol[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC_DARCY(bcPl, sol, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC_DARCY(bcPl, sol, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC_DARCY(bcPl, sol, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC_DARCY(bcPl, sol, k, J, i, pmdof) }
		if(k == 0)   {
			fk = 1; K = k-1;
			SET_TPC_DARCY(bcPl, sol, K, j, i, pmdof)
		}
		if(k == mcz) {
			fk = 1; K = k+1;
			SET_TPC_DARCY(bcPl, sol, K, j, i, pmdof)
		}

		if(fi*fj)    SET_EDGE_CORNER_DARCY(n, sol, k, J, I, k, j, i, pmdof)
		if(fi*fk)    SET_EDGE_CORNER_DARCY(n, sol, K, j, I, k, j, i, pmdof)
		if(fj*fk)    SET_EDGE_CORNER_DARCY(n, sol, K, J, i, k, j, i, pmdof)
		if(fi*fj*fk) SET_EDGE_CORNER_DARCY(n, sol, K, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_Pl, jr->lPl,  &sol);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_Pl, bc->bcPl, &bcPl); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetDarcyRes"
PetscErrorCode JacResGetDarcyRes(JacRes *jr)
{
	// compute Darcy residual vector

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	BCCtx      		*bc;
	PetscInt    	iter, num, *list, numPhases;
	PetscInt    	Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    	i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, l;
	//PetscInt 		sx_fine,sy_fine,sz_fine;
 	PetscScalar 	bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar 	bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar	 	bqx, fqx, bqy, fqy, bqz, fqz;
 	PetscScalar 	dt, dx, dy, dz;
 	PetscScalar 	Pl_c, Pl_h, pc, gz, Pl_n, *grav;
 	PetscScalar 	***sol, ***hydro, ***p, ***res, *e;

 	PetscScalar 	lrho, mu, Kphi, betam, betal, Phi, Ts, nuu, Kphiu, Phiu, Ss;
	PetscScalar 	***rhol, ***mul, ***lKphi, ***lbetam, ***lbetal, ***lPhi, ***Tsl, ***lnuu, ***lKphiu, ***lPhiu; //, ***Ss

	Vec 			rhol_local, mul_local, local_Kphi, betam_local, betal_local, Phi_local, Ts_local, nuu_local, Kphiu_local, Phiu_local; //Ss_local, ;
	Vec 			rhol_vec, mul_vec, Kphi_vec, betam_vec, betal_vec, Phi_vec, Ts_vec, nuu_vec, Kphiu_vec, Phiu_vec;   //, Ss_vec

	PetscScalar     H, magnitude, increment, tini, tfin;

	PetscScalar     p_total, dP, biot, Kd, yield;
	PetscScalar     min_overpressure, max_overpressure, dPmin, dPmax;
	PetscInt src;
	MatParLim  *matLim;
	Material_t *phases;
	PetscScalar scal;



	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs        = jr->fs;
	grav      = jr->grav;       // gravity acceleration
	gz 		  =	grav[2];		// gravitational acceleration in z-direction
	bc        = jr->bc;
	dt        = jr->ts.dt;      // time step
	num       = bc->Pl_NumSPC;
	list      = bc->Pl_SPCList;


	scal   = jr->scal.length;

	// access residual context variables
	numPhases =  jr->numPhases; 	// number phases
	phases    =  jr->phases;    	// phase parameters
	matLim    = &jr->matLim;    	// phase parameters limiters

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	// Extract required material parameters
	ExtractCoefficientsFromDA(jr->DA_Pl, "rhol",  rhol_local,   rhol_vec, 	rhol);		// liquid density
	ExtractCoefficientsFromDA(jr->DA_Pl, "mul",	  mul_local,    mul_vec, 	mul);		// Liquid viscosity
	ExtractCoefficientsFromDA(jr->DA_Pl,"Kphi",	  local_Kphi,  	Kphi_vec,  	lKphi);		// Permeability K
	//ExtractCoefficientsFromDA(jr->DA_Pl, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);	// specific storage
	ExtractCoefficientsFromDA(jr->DA_Pl, "betam", betam_local,  betam_vec, 	lbetam);    // matrix compressibility
	ExtractCoefficientsFromDA(jr->DA_Pl, "betal", betal_local,  betal_vec, 	lbetal);    // liquid compressibility
	ExtractCoefficientsFromDA(jr->DA_Pl, "Phi",   Phi_local,    Phi_vec, 	lPhi);      // porosity
	ExtractCoefficientsFromDA(jr->DA_Pl, "Tsl",   Ts_local,     Ts_vec, 	Tsl);       // tensile strength
	ExtractCoefficientsFromDA(jr->DA_Pl, "nuu",   nuu_local, 	nuu_vec, 	lnuu);      // undrained Poisson's ratio
	ExtractCoefficientsFromDA(jr->DA_Pl, "Kphiu", Kphiu_local, 	Kphiu_vec, 	lKphiu);    // undrained permeability
	ExtractCoefficientsFromDA(jr->DA_Pl, "Phiu",  Phiu_local, 	Phiu_vec, 	lPhiu);     // undrained porosity

	// access work vectors
	ierr = DMDAVecGetArray(jr->DA_Pl,  jr->r_Pl,          &res);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lPl,           &sol);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->hydro_lPl,     &hydro);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,            &p);      CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	iter = 0;
	START_STD_LOOP
	{

		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		Pl_n  = svBulk->Pln;  	// liquid pressure history
		Pl_c  = sol[k][j][i]; 	// current liquid pressure
		Pl_h  = hydro[k][j][i]; // hydrostatic liquid pressure
		pc    = p[k][j][i];     // current pressure

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		lrho	=   rhol[k][j][i];
		mu		=	mul[k][j][i];
		Kphi 	=	lKphi[k][j][i];
		//Ss		= 	Ssl[k][j][i];
		betam 	=	lbetam[k][j][i];
		betal 	=	lbetal[k][j][i];
		Phi 	=	lPhi[k][j][i];
		Ts 	    =	Tsl[k][j][i];
		nuu 	=	lnuu[k][j][i];
		Kphiu 	=	lKphiu[k][j][i];
		Phiu 	=	lPhiu[k][j][i];

		//pShift?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		matLim    = &jr->matLim;       // phase parameters limiters
		biot      = matLim->biot;      // Biot pressure parameter
		p_total   = pc + biot*Pl_c;
		dP        = p_total - Pl_c;

		dPmin= pc+biot*Pl_h-Pl_h;
		dPmax= Ts;

		// Liquid density
		svBulk->Rhol = lrho;

		// Permeability and porosity
		if (dP <= dPmin && dP > dPmax)
		{
			if (Kphiu && Kphi != Kphiu)
			{
				Kphi = Kphi + (Kphiu - Kphi)/(dPmax-dPmin)*(dP-dPmin); //Kphi = svBulk->Kphi * exp(-(dP/Kd));
			}
			//if (Phiu  && Phi  != Phiu)
			//{
			//	Phi  = Phi + (Phiu - Phi)/(dPmax-dPmin)*(dP-dPmin);
			//}
		}
		else if (dP > dPmin)
		{
			//!!!!
		}
		else
		{
			// Tensile failure
			Kphi = Kphiu;
			//Phi  = Phiu;
		}


		svBulk->Kphi = Kphi;
		svBulk->Phi  = Phi;

		// Specific storage
		Ss = (betam + Phi*betal);

		bkx 	= 	((Kphi + lKphi[k][j][Im1])/2.0);       fkx = ((Kphi + lKphi[k][j][Ip1])/2.0);
		bky 	= 	((Kphi + lKphi[k][Jm1][i])/2.0);       fky = ((Kphi + lKphi[k][Jp1][i])/2.0);
		bkz 	= 	((Kphi + lKphi[Km1][j][i])/2.0);       fkz = ((Kphi + lKphi[Kp1][j][i])/2.0);

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);


		// SOURCE ----------------------------------------------------------------------------------------------------------------------------
		H = 0.0;
		// read each source
		for(src = 0; src < jr->NumDarcySources; src++)
		{
			if (JacResGetTime(jr) >= jr->DarcySources[src].tini && JacResGetTime(jr) <= jr->DarcySources[src].tfin)
			{
				if (jr->DarcySources[src].i == i && jr->DarcySources[src].j == j && jr->DarcySources[src].k == k)
				{
					increment = jr->DarcySources[src].increment;
					magnitude = jr->DarcySources[src].magnitude *      (1.0 +(JacResGetTime(jr)-jr->DarcySources[src].tini)*increment); // [m^3/s]
					// If source is Volumetric flow [m^3/s] then H = magnitude / volume of the cell [1/s]
					//H = H/(dx*dy*dz);
					//H = H + magnitude/(dx*scal*dy*scal*dz*scal);
					H = H + magnitude/(dx*dy*dz);
				}
			}
		}
		// ------------------------------------------------------------------------------------------------------------------------------------

		// compute LiquidPressure fluxes
		bqx = bkx/mu*(Pl_c - sol[k][j][i-1])/bdx;   									fqx = fkx/mu*(sol[k][j][i+1] - Pl_c)/fdx;
		bqy = bky/mu*(Pl_c - sol[k][j-1][i])/bdy;   									fqy = fky/mu*(sol[k][j+1][i] - Pl_c)/fdy;
		//bqz = bkz/mu*( (   Pl_c - sol[k-1][j][i]  )/bdz ); 	fqz = fkz/mu*( (   sol[k+1][j][i] - Pl_c)/fdz );
		//bqz = bkz/mu*( (   (Pl_c-Pl_h) - (sol[k-1][j][i]-hydro[k-1][j][i])  )/bdz ); 	fqz = fkz/mu*((   (sol[k+1][j][i]-hydro[k+1][j][i]) - (Pl_c-Pl_h))/fdz );
		bqz = bkz/mu*( (   Pl_c - sol[k-1][j][i]  )/bdz - lrho*gz); 	fqz = fkz/mu*( (   sol[k+1][j][i] - Pl_c)/fdz - lrho*gz); // change if liquid density is not constant


		// compute Liquid velocity flow
		svBulk->liquidvelocity[0] = -(fqx+bqx)/2.0;
		svBulk->liquidvelocity[1] = -(fqy+bqy)/2.0;
		svBulk->liquidvelocity[2] = -(fqz+bqz)/2.0;

		// Residual
		res[k][j][i] = Ss*(Pl_c-Pl_n)/dt -   (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz - H;


	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lPl,   	 &sol);  	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_Pl,  jr->r_Pl,   	 &res);  	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->hydro_lPl,&hydro);  	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,       &p);  	    CHKERRQ(ierr);

	// restore material coefficients on DA
	RestoreCoefficientsFromDA(jr->DA_Pl, "rhol", rhol_local,    rhol_vec,   rhol);
	RestoreCoefficientsFromDA(jr->DA_Pl, "mul",	 mul_local,     mul_vec,    mul);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Kphi", local_Kphi,    Kphi_vec,   lKphi);
	//RestoreCoefficientsFromDA(jr->DA_Pl, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);
	RestoreCoefficientsFromDA(jr->DA_Pl, "betam",betam_local,   betam_vec,  lbetam);
	RestoreCoefficientsFromDA(jr->DA_Pl, "betal",betal_local,   betal_vec,  lbetal);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Phi",  Phi_local,     Phi_vec,    lPhi);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Tsl",  Ts_local,      Ts_vec,     Tsl);
	RestoreCoefficientsFromDA(jr->DA_Pl, "nuu",  nuu_local,     nuu_vec,    lnuu);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Kphiu",Kphiu_local,   Kphiu_vec,  lKphiu);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Phiu", Phiu_local,    Phiu_vec,    lPhiu);

	// impose primary pressure constraints
	ierr = VecGetArray(jr->r_Pl, &e); CHKERRQ(ierr);

	for(i = 0; i < num; i++) e[list[i]] = 0.0;

	ierr = VecRestoreArray(jr->r_Pl, &e); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetDarcyMat"
PetscErrorCode JacResGetDarcyMat(JacRes *jr)
{
	// Darcy matrix
	// Ss*Pl/dt - d/dxi (Kphi/mu (dPl/dxj))

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	BCCtx      *bc;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dt, dx, dy, dz;
 	PetscScalar 	Pl_c, Pl_h, pc, gz, Pl_n;
 	 PetscScalar 	***sol, ***hydro, ***p;
	PetscScalar v[7], cf[6], lrho, mu, Kphi, betam, betal, Phi, Ts, nuu, Kphiu, Phiu, Ss;
	MatStencil  row[1], col[7];
	PetscScalar 	***rhol, ***mul, ***lKphi, ***lbetam, ***lbetal, ***lPhi, ***Tsl, ***lnuu, ***lKphiu, ***lPhiu, ***bcPl;
	//Vec         rhol_local, mul_local, local_Kphi;
	//Vec 		rhol_vec,   mul_vec,   Kphi_vec;
	Vec 			rhol_local, mul_local, local_Kphi, betam_local, betal_local, Phi_local, Ts_local, nuu_local, Kphiu_local, Phiu_local; //Ss_local, ;
	Vec 			rhol_vec, mul_vec, Kphi_vec, betam_vec, betal_vec, Phi_vec, Ts_vec, nuu_vec, Kphiu_vec, Phiu_vec;   //, Ss_vec


	PetscScalar     p_total, dP, biot, Kd, yield;
	PetscScalar     min_overpressure, max_overpressure, dPmin, dPmax;
	PetscInt src;
	MatParLim  *matLim;
	Material_t *phases;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// access residual context variables
	fs        = jr->fs;
	bc        = jr->bc;
	num       = bc->Pl_NumSPC;
	list      = bc->Pl_SPCList;
	dt        = jr->ts.dt;     	// time step

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	// clear matrix coefficients
	ierr = MatZeroEntries(jr->App); CHKERRQ(ierr);

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcPl,  &bcPl); 	CHKERRQ(ierr);

	// Extract required material parameters
	ExtractCoefficientsFromDA(jr->DA_Pl, "rhol",  rhol_local,   rhol_vec, 	rhol);		// liquid density
	ExtractCoefficientsFromDA(jr->DA_Pl, "mul",	  mul_local,    mul_vec, 	mul);		// Liquid viscosity
	ExtractCoefficientsFromDA(jr->DA_Pl,"Kphi",	  local_Kphi,  	Kphi_vec,  	lKphi);		// Permeability K
	//ExtractCoefficientsFromDA(jr->DA_Pl, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);	// specific storage
	ExtractCoefficientsFromDA(jr->DA_Pl, "betam", betam_local,  betam_vec, 	lbetam);    // matrix compressibility
	ExtractCoefficientsFromDA(jr->DA_Pl, "betal", betal_local,  betal_vec, 	lbetal);    // liquid compressibility
	ExtractCoefficientsFromDA(jr->DA_Pl, "Phi",   Phi_local,    Phi_vec, 	lPhi);      // porosity
	ExtractCoefficientsFromDA(jr->DA_Pl, "Tsl",   Ts_local,     Ts_vec, 	Tsl);       // tensile strength
	ExtractCoefficientsFromDA(jr->DA_Pl, "nuu",   nuu_local, 	nuu_vec, 	lnuu);      // undrained Poisson's ratio
	ExtractCoefficientsFromDA(jr->DA_Pl, "Kphiu", Kphiu_local, 	Kphiu_vec, 	lKphiu);    // undrained permeability
	ExtractCoefficientsFromDA(jr->DA_Pl, "Phiu",  Phiu_local, 	Phiu_vec, 	lPhiu);     // undrained porosity

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lPl,           &sol);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->hydro_lPl,     &hydro);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,            &p);      CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		Pl_n  = svBulk->Pln;  	// liquid pressure history
		Pl_c  = sol[k][j][i]; 	// current liquid pressure
		Pl_h  = hydro[k][j][i]; // hydrostatic liquid pressure
		pc    = p[k][j][i];     // current pressure

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcPl[k][j][i-1] != DBL_MAX) cf[0] = -1.0;}
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcPl[k][j][i+1] != DBL_MAX) cf[1] = -1.0;}
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcPl[k][j-1][i] != DBL_MAX) cf[2] = -1.0;}
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcPl[k][j+1][i] != DBL_MAX) cf[3] = -1.0;}
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcPl[k-1][j][i] != DBL_MAX) cf[4] = -1.0;}
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcPl[k+1][j][i] != DBL_MAX) cf[5] = -1.0;}

		lrho	=   rhol[k][j][i];
		mu		=	mul[k][j][i];
		Kphi 	=	lKphi[k][j][i];
		//Ss		= 	Ssl[k][j][i];
		betam 	=	lbetam[k][j][i];
		betal 	=	lbetal[k][j][i];
		Phi 	=	lPhi[k][j][i];
		Ts 	    =	Tsl[k][j][i];
		nuu 	=	lnuu[k][j][i];
		Kphiu 	=	lKphiu[k][j][i];
		Phiu 	=	lPhiu[k][j][i];

		//pShift?!!
		matLim    = &jr->matLim;       // phase parameters limiters
		biot      = matLim->biot;      // Biot pressure parameter
		p_total   = pc + biot*Pl_c;
		dP        = p_total - Pl_c;

		dPmin= pc+biot*Pl_h-Pl_h;
		dPmax= Ts;

		/*//// Permeability and porosity
		if (dP <= dPmin && dP > dPmax)
		{
			if (Kphiu && Kphi != Kphiu)
			{
				Kphi = Kphi + (Kphiu - Kphi)/(dPmax-dPmin)*(dP-dPmin); //Kphi = svBulk->Kphi * exp(-(dP/Kd));
			}
			if (Phiu  && Phi  != Phiu)
			{
				Phi  = Phi + (Phiu - Phi)/(dPmax-dPmin)*(dP-dPmin);
			}
		}
		else if (dP > dPmin)
		{
			//!!!!
		}
		else
		{
			// Tensile failure
			Kphi = Kphiu;
			Phi  = Phiu;
		}*/

		// Specific storage
		Ss = (betam + Phi*betal);

		// compute average permeabilities
		bkx 	= 	(Kphi + lKphi[k][j][Im1])/2.0;       fkx = (Kphi + lKphi[k][j][Ip1])/2.0;
		bky 	= 	(Kphi + lKphi[k][Jm1][i])/2.0;       fky = (Kphi + lKphi[k][Jp1][i])/2.0;
		bkz 	= 	(Kphi + lKphi[Km1][j][i])/2.0;       fkz = (Kphi + lKphi[Kp1][j][i])/2.0;

		// liquid viscosity
		bkx /= mu;							fkx /= mu;
		bky /= mu;							fky /= mu;
		bkz /= mu;							fkz /= mu;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC multipliers
		v[0] =  -bkx/bdx/dx*cf[0];
		v[1] =  -fkx/fdx/dx*cf[1];
		v[2] =  -bky/bdy/dy*cf[2];
		v[3] =  -fky/fdy/dy*cf[3];
		v[4] =  -bkz/bdz/dz*cf[4];
		v[5] =  -fkz/fdz/dz*cf[5];
		v[6] =  Ss/dt
		+       (bkx/bdx + fkx/fdx)/dx
		+       (bky/bdy + fky/fdy)/dy
		+       (bkz/bdz + fkz/fdz)/dz;

		// set matrix coefficients
		ierr = MatSetValuesStencil(jr->App, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(jr->DA_Pl, bc->bcPl, &bcPl); CHKERRQ(ierr);

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lPl,   	 &sol);  	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->hydro_lPl,&hydro);  	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,       &p);  	    CHKERRQ(ierr);

	// restore material coefficients on DA
	RestoreCoefficientsFromDA(jr->DA_Pl, "rhol", rhol_local,    rhol_vec,   rhol);
	RestoreCoefficientsFromDA(jr->DA_Pl, "mul",	 mul_local,     mul_vec,    mul);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Kphi", local_Kphi,    Kphi_vec,   lKphi);
	//RestoreCoefficientsFromDA(jr->DA_Pl, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);
	RestoreCoefficientsFromDA(jr->DA_Pl, "betam",betam_local,   betam_vec,  lbetam);
	RestoreCoefficientsFromDA(jr->DA_Pl, "betal",betal_local,   betal_vec,  lbetal);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Phi",  Phi_local,     Phi_vec,    lPhi);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Tsl",  Ts_local,      Ts_vec, 	Tsl);
	RestoreCoefficientsFromDA(jr->DA_Pl, "nuu",  nuu_local,     nuu_vec,    lnuu);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Kphiu",Kphiu_local,   Kphiu_vec,  lKphiu);
	RestoreCoefficientsFromDA(jr->DA_Pl, "Phiu", Phiu_local,    Phiu_vec,   lPhiu);

	// assemble LiquidPressure/Darcy matrix
	ierr = MatAIJAssemble(jr->App, num, list, 1.0); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "UpdateDarcy_DA"
PetscErrorCode UpdateDarcy_DA(JacRes *jr)
{
	// This routine
	//	1) Updates the coordinates attached to the Darcy DA (as the grid might have
	//		been deformed
	//
	//	2) Updates material properties, we later need to construct the residual and Jacobian
	//		These properties will be attached to the DA using a hook, as this automatically takes
	//		care of refining them to coarser levels if we use MG

	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscInt 	   	iter, sx, sy, sz;
	PetscInt       	i, j, k, nx, ny, nz;
	Vec 			rhol_vec,   local_rhol;
	Vec 			Mul_vec,    local_Mul;
	Vec 			Kphi_vec,   local_KPhi;
	//Vec 			Ss_vec,	    local_Ss;
	Vec 			betam_vec,  local_betam;
	Vec 			betal_vec,  local_betal;
	Vec 			Phi_vec,    local_Phi;
	Vec 			Ts_vec,     local_Ts;
	Vec 			nuu_vec,	local_nuu;
	Vec 			Kphiu_vec,  local_Kphiu;
	Vec 			Phiu_vec,   local_Phiu;
	Vec 			BC_vec,     local_BC;
	PetscScalar    	***buff, rhol, mul, Kphi, betam, betal, Phi, Tsl, nuu, Kphiu, Phiu; //, Ssl

	//---------------------------------------------
	// (1) Update coordinates on 3D DMDA object
	ierr =  UpdateCoordinatesOnCellCenterDA(jr->DA_Pl, jr->fs);	CHKERRQ(ierr);
	//---------------------------------------------

	//---------------------------------------------
	// (2) Update material properties on local vector & attach to DM

	// Create temporary local vectors
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_rhol);				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Mul); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_KPhi); 				CHKERRQ(ierr);
	//ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Ss);					CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_betam);				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_betal); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Phi); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Ts);					CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_nuu); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Kphiu); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Phiu); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_BC); 					CHKERRQ(ierr);

	// Set data on local vectors
	SCATTER_FIELD(jr->DA_Pl, local_rhol, GET_rhol);            // liquid density
	SCATTER_FIELD(jr->DA_Pl, local_Mul,  GET_mul);             // compute 3D array of liquid viscosity
	SCATTER_FIELD(jr->DA_Pl, local_KPhi, GET_Kphi);            // 3D array of permeability pre-coefficient
	//SCATTER_FIELD(jr->DA_Pl, local_Ss,   GET_Ssl);            // specific storage
	SCATTER_FIELD(jr->DA_Pl, local_betam,GET_betam);
	SCATTER_FIELD(jr->DA_Pl, local_betal,GET_betal);
	SCATTER_FIELD(jr->DA_Pl, local_Phi,  GET_Phi);
	SCATTER_FIELD(jr->DA_Pl, local_Ts,   GET_Tsl);				// tensile strength
	SCATTER_FIELD(jr->DA_Pl, local_nuu,  GET_nuu);
	SCATTER_FIELD(jr->DA_Pl, local_Kphiu,GET_Kphiu);
	SCATTER_FIELD(jr->DA_Pl, local_Phiu, GET_Phiu);

	// Get/Create global vector, attached to DM
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"rhol",        &rhol_vec);         CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"mul",         &Mul_vec);          CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Kphi",        &Kphi_vec);         CHKERRQ(ierr);
	//ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Ssl",			&Ss_vec);         CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"betam",       &betam_vec);        CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"betal",       &betal_vec);        CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Phi",         &Phi_vec);          CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Tsl",         &Ts_vec);           CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"nuu",         &nuu_vec);          CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Kphiu",       &Kphiu_vec);        CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Phiu",        &Phiu_vec);         CHKERRQ(ierr);

	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"BC",          &BC_vec);           CHKERRQ(ierr);

	// Move local values to global vector

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_rhol,INSERT_VALUES ,rhol_vec);		CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_rhol,INSERT_VALUES ,rhol_vec);		CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Mul,INSERT_VALUES ,Mul_vec);		CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Mul,INSERT_VALUES ,Mul_vec);		CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_KPhi,INSERT_VALUES ,Kphi_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_KPhi,INSERT_VALUES ,Kphi_vec);	CHKERRQ(ierr);

	//ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Ss,INSERT_VALUES ,Ss_vec);	CHKERRQ(ierr);
	//ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Ss,INSERT_VALUES ,Ss_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_betam,INSERT_VALUES ,betam_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_betam,INSERT_VALUES ,betam_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_betal,INSERT_VALUES ,betal_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_betal,INSERT_VALUES ,betal_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Phi,INSERT_VALUES ,Phi_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Phi,INSERT_VALUES ,Phi_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Ts,INSERT_VALUES ,Ts_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Ts,INSERT_VALUES ,Ts_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_nuu,INSERT_VALUES ,nuu_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_nuu,INSERT_VALUES ,nuu_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Kphiu,INSERT_VALUES ,Kphiu_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Kphiu,INSERT_VALUES ,Kphiu_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Phiu,INSERT_VALUES ,Phiu_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Phiu,INSERT_VALUES ,Phiu_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,jr->bc->bcPl,INSERT_VALUES ,BC_vec);			CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,jr->bc->bcPl,INSERT_VALUES ,BC_vec);			CHKERRQ(ierr);

	// attach the named vector to the DM

	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"rhol",		&rhol_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"mul",			&Mul_vec); 				CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Kphi",		&Kphi_vec); 			CHKERRQ(ierr);
	//ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Ssl",			&Ss_vec); 				CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"betam",		&betam_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"betal",		&betal_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Phi",		    &Phi_vec); 			    CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Tsl",			&Ts_vec); 				CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"nuu",         &nuu_vec);              CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Kphiu",		&Kphiu_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Phiu",		&Phiu_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"BC",			&BC_vec); 				CHKERRQ(ierr);

	// cleaning up
	ierr = 	VecDestroy(&local_rhol);		CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Mul);			CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_KPhi);		CHKERRQ(ierr);
	//ierr = 	VecDestroy(&local_Ss);			CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_betam);		CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_betal);		CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Phi);         CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Ts);			CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_nuu);         CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Kphiu);       CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Phiu);        CHKERRQ(ierr);

	ierr = 	VecDestroy(&local_BC);			CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SolveDarcyKSP"
PetscErrorCode SolveDarcyKSP(JacRes *jr)
{
	PetscFunctionBegin;
	PetscErrorCode ierr;

	ierr = JacResGetDarcyRes(jr);   		CHKERRQ(ierr);
	ierr = JacResGetDarcyMat(jr);             				CHKERRQ(ierr);
	ierr = KSPSetOperators(jr->Pl_ksp, jr->App, jr->App); 	CHKERRQ(ierr);
	ierr = KSPSetUp(jr->Pl_ksp);                            CHKERRQ(ierr);
	ierr = KSPSolve(jr->Pl_ksp, jr->r_Pl, jr->dPl);      	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitDarcy"
PetscErrorCode JacResInitDarcy(JacRes *jr)
{
	// initialize liquid pressure from markers, or from zero

	FDSTAG     *fs;
	BCCtx       *bc;
	PetscScalar ***lPl, ***bcPl, Pl;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;
	bc = jr->bc;

	ierr = VecZeroEntries(jr->lPl); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lPl,  &lPl);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcPl, &bcPl); CHKERRQ(ierr);

	iter = 0;

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{

		Pl = bcPl[k][j][i];

		if(Pl == DBL_MAX) Pl = jr->svCell[iter].svBulk.Pln;

		lPl[k][j][i] = Pl;

		iter++;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lPl,  &lPl);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcPl, &bcPl); CHKERRQ(ierr);

	// apply two-point constraints
	ierr = JacResApplyDarcyBC(jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
// Source
#undef __FUNCT__
#define __FUNCT__ "GetCellCoordinatesDarcySources"
PetscErrorCode GetCellCoordinatesDarcySources(JacRes *jr)
{
	PetscInt src, i, j, k, nx, ny, nz, sx, sy, sz, ii, jj, kk;
	PetscScalar  xs[3], xe[3], ys[3], ye[3], zs[3], ze[3];
	FDSTAG     *fs;
	PetscScalar scal, x,y,z,dx,dy,dz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	scal   = jr->scal.length;

	// read each source
	for(src = 0; src < jr->NumDarcySources; src++)
	{
		ii = -1;
		jj = -1;
		kk = -1;

		//-------------------------------
		// central points
		//-------------------------------
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)
		START_STD_LOOP
		{
			// get coordinate of center of cell
			x  = COORD_CELL(i, sx, fs->dsx);
			y  = COORD_CELL(j, sy, fs->dsy);
			z  = COORD_CELL(k, sz, fs->dsz);

			// size of cell
			dx = SIZE_CELL(i, sx, fs->dsx);
			dy = SIZE_CELL(j, sy, fs->dsy);
			dz = SIZE_CELL(k, sz, fs->dsz);

			// get cell coordinates
			xs[0] = x-dx/2; xe[0] = x+dx/2;
			ys[1] = y-dy/2; ye[1] = y+dy/2;
			zs[2] = z-dz/2; ze[2] = z+dz/2;

			if (jr->DarcySources[src].x/scal > xs[0] && jr->DarcySources[src].x/scal <= xe[0] && jr->DarcySources[src].y/scal > ys[1] && jr->DarcySources[src].y/scal <= ye[1] && jr->DarcySources[src].z/scal > zs[2] && jr->DarcySources[src].z/scal <= ze[2])
			{
				ii=i;
				jj=j;
				kk=k;
			}
		}
		END_STD_LOOP

		if(ISParallel(PETSC_COMM_WORLD))
		{
			ierr = MPI_Allreduce(&ii, &jr->DarcySources[src].i, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
			ierr = MPI_Allreduce(&jj, &jr->DarcySources[src].j, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
			ierr = MPI_Allreduce(&kk, &jr->DarcySources[src].k, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
		}
		else
		{
			jr->DarcySources[src].i = ii;
			jr->DarcySources[src].j = jj;
			jr->DarcySources[src].k = kk;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DarcySourcePropInit"
PetscErrorCode DarcySourcePropInit(JacRes *jr, FILE *fp)
{
	// initialize source properties from file

	PetscInt  *ls, *le;
	PetscInt   i, count_starts, count_ends;
	PetscInt found, sources;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Reading liquid source parameters: \n\n");

	// clear memory
	ierr = PetscMemzero(jr->DarcySources, sizeof(DarcySourceParam)*(size_t)Max_Num_Darcy_Sources); CHKERRQ(ierr);

	parse_GetInt(fp, "NumDarcySources",&sources, &found);
	if (found==PETSC_TRUE)
	{
		jr->NumDarcySources=sources;
	}else
	{
		jr->NumDarcySources=0;
	}

	// initialize ID for consistency checks
	for(i = 0; i < Max_Num_Darcy_Sources; i++)
	{
		jr->DarcySources[i].ID = -1;
	}

	// allocate memory for arrays to store line info
	ierr = makeIntArray(&ls, NULL, Max_Num_Darcy_Sources); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, Max_Num_Darcy_Sources); CHKERRQ(ierr);

	// read number of entries
	getLineStruct(fp, ls, le, Max_Num_Darcy_Sources, &count_starts, &count_ends, "<DarcySourceStart>","<DarcySourceEnd>");

	// error checking
	if(count_starts > Max_Num_Darcy_Sources || count_ends > Max_Num_Darcy_Sources)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many source specified for Darcy simulation! Max allowed: %lld", (LLD)Max_Num_Darcy_Sources);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incomplete source structures for Darcy simulation! <DarcySourceStart> & <DarcySourceEnd> don't match");
	}

	if(count_starts < jr->NumDarcySources || count_ends < jr->NumDarcySources)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too few source parameters specified for Darcy simulation! Max allowed: %lld", (LLD)Max_Num_Darcy_Sources);
	}

	// read each individual phase
	for(i = 0; i < jr->NumDarcySources; i++)
	{
		ierr = SourcePropGetStruct(fp,
				jr->NumDarcySources, jr->DarcySources,
				ls[i], le[i], jr->scal.utype); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SourcePropGetStruct"
PetscErrorCode SourcePropGetStruct(FILE *fp,
		PetscInt numSources, DarcySourceParam *sources,
		PetscInt ils, PetscInt ile, UnitsType utype)
{
	// read source properties from file with error checking

	DarcySourceParam *m;
	PetscInt    ID = -1, found;

	// output labels
	char        lbl_x  			[_lbl_sz_];
	char        lbl_magnitude  	[_lbl_sz_];
	char        lbl_tini  		[_lbl_sz_];
	char        lbl_tfin  		[_lbl_sz_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// phase ID
	getMatPropInt(fp, ils, ile, "ID", &ID, &found);

	// error checking
	if(!found)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "No source ID specified for Darcy code! ");
	}
	if(ID > numSources - 1)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Incorrect source numbering for Darcy code!");
	}

	// get pointer to specified phase
	m = sources + ID;

	// check ID
	if(m->ID != -1)
	{
		 SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Incorrect source numbering for Darcy code!");
	}

	// set ID
	m->ID = ID;

	getMatPropScalar(fp, ils, ile, "x",         &m->x,          NULL);
	getMatPropScalar(fp, ils, ile, "y",         &m->y,          NULL);
	getMatPropScalar(fp, ils, ile, "z",         &m->z,          NULL);
	getMatPropScalar(fp, ils, ile, "magnitude", &m->magnitude,  NULL); // pressure-dependence of density
	getMatPropScalar(fp, ils, ile, "increment", &m->increment,  NULL); // pressure-dependence of density
	getMatPropScalar(fp, ils, ile, "tini",      &m->tini,       NULL); // Time to start the source
	getMatPropScalar(fp, ils, ile, "tfin",      &m->tfin,       NULL); // Time to stop the source

	// print
	if(utype == _NONE_)
	{
		sprintf(lbl_x,   		 "[ ]"     );
		sprintf(lbl_magnitude,   "[ ]"     );
		sprintf(lbl_tini,        "[ ]"     );
		sprintf(lbl_tfin,        "[ ]"     );
	}
	else
	{
		sprintf(lbl_x,   		"[m]"     );
		sprintf(lbl_magnitude,  "[1/s]"   );
		sprintf(lbl_tini,       "[s]"     );
		sprintf(lbl_tfin,       "[s]"     );
	}

	PetscPrintf(PETSC_COMM_WORLD,"    Source [%lld]: x = %g %s, y = %g %s, z = %g %s, magnitude = %g %s, tini = %g %s, tfin = %g %s", (LLD)(m->ID), m->x, lbl_x, m->y, lbl_x, m->z, lbl_x, m->magnitude, lbl_magnitude, m->tini, lbl_tini, m->tfin, lbl_tfin);
	PetscPrintf(PETSC_COMM_WORLD,"    \n");

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "UpdateFailureType"
PetscErrorCode UpdateFailureType(JacRes *jr)
{
	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	MatParLim  *matLim;
	PetscInt    iter, numPhases;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, l;
	PetscScalar XX, XX1, XX2, XX3, XX4;
	PetscScalar YY, YY1, YY2, YY3, YY4;
	PetscScalar ZZ, ZZ1, ZZ2, ZZ3, ZZ4;
	PetscScalar XY, XY1, XY2, XY3, XY4;
	PetscScalar XZ, XZ1, XZ2, XZ3, XZ4;
	PetscScalar YZ, YZ1, YZ2, YZ3, YZ4;
	PetscScalar pc, pc_pore, pc_hydro, ptotal, biot, dt, dP, TensileS, intersection, fr, ch, sensitivity;
	PetscScalar ***dxx, ***dyy, ***dzz, ***dxy, ***dxz, ***dyz, ***p, ***p_pore, ***p_hydro;
	PetscBool actDarcy;
	PetscInt src, frac, frac2;

	PetscScalar sxx,syy,szz;
	PetscScalar DevStressII;

	// To plot mohr's cercles
	PetscScalar C1, C2, C3, R1, R2, R3;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	numPhases =  jr->numPhases; 	// number phases
	phases    =  jr->phases;    	// phase parameters
	matLim    = &jr->matLim;    	// phase parameters limiters
	dt        =  jr->ts.dt;     	// time step
	biot      =  matLim->biot;      // Biot pressure parameter
	actDarcy  = jr->actDarcy;
	sensitivity = 0.0;
	src = 0;
	frac=0; frac2=0;

	if (actDarcy == PETSC_TRUE) {

		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dxx);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dyy);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldzz, &dzz);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &dxy);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &dxz);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &dyz);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,       &p);        CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp_pore,  &p_pore);   CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->hydro_lPl,  &p_hydro);   CHKERRQ(ierr);

		iter = 0;

		//-------------------------------
		// xy edge points
		//-------------------------------
		iter = 0;
		GET_NODE_RANGE(nx, sx, fs->dsx)
		GET_NODE_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			svEdge = &jr->svXYEdge[iter++];
			dxy[k][j][i] = svEdge->s;
			svEdge->svDev.failS = 0.0;
			svEdge->svDev.failT = 0.0;
			svEdge->svDev.failTS = 0.0;

			//dxy[k][j][i] = jr->svXYEdge[iter++].s;


		}
		END_STD_LOOP

		//ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &dxy); CHKERRQ(ierr);

		LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);

		//-------------------------------
		// xz edge points
		//-------------------------------
		iter = 0;
		GET_NODE_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_NODE_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			svEdge = &jr->svXZEdge[iter++];
			dxz[k][j][i] = svEdge->s;
			svEdge->svDev.failS = 0.0;
			svEdge->svDev.failT = 0.0;
			svEdge->svDev.failTS = 0.0;

			// access solution variables
			//dxz[k][j][i] = jr->svXZEdge[iter++].s;

		}
		END_STD_LOOP

		//ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &dxz); CHKERRQ(ierr);

		LOCAL_TO_LOCAL(fs->DA_XZ, jr->ldxz);

		//-------------------------------
		// yz edge points
		//-------------------------------
		iter = 0;
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_NODE_RANGE(ny, sy, fs->dsy)
		GET_NODE_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			svEdge = &jr->svYZEdge[iter++];
			dyz[k][j][i] = svEdge->s;
			svEdge->svDev.failS = 0.0;
			svEdge->svDev.failT = 0.0;
			svEdge->svDev.failTS = 0.0;

			// access solution variables
			//dyz[k][j][i] = jr->svYZEdge[iter++].s;

		}
		END_STD_LOOP

		//ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &dyz); CHKERRQ(ierr);

		LOCAL_TO_LOCAL(fs->DA_YZ, jr->ldyz);

		//-------------------------------
		// central points
		//-------------------------------
		iter = 0;
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			frac = 0;


			// access solution variables
			svCell = &jr->svCell[iter++];
			svDev  = &svCell->svDev;
			svBulk = &svCell->svBulk;

			//=================
			// SECOND INVARIANT of deviatoric stress
			//=================

			// x-y plane, i-j indices
			XY1 = dxy[k][j][i];
			XY2 = dxy[k][j+1][i];
			XY3 = dxy[k][j][i+1];
			XY4 = dxy[k][j+1][i+1];

			// x-z plane, i-k indices
			XZ1 = dxz[k][j][i];
			XZ2 = dxz[k+1][j][i];
			XZ3 = dxz[k][j][i+1];
			XZ4 = dxz[k+1][j][i+1];

			// y-z plane, j-k indices
			YZ1 = dyz[k][j][i];
			YZ2 = dyz[k+1][j][i];
			YZ3 = dyz[k][j+1][i];
			YZ4 = dyz[k+1][j+1][i];

			// compute second invariant
			DevStressII = 0.5*(svCell->sxx*svCell->sxx + svCell->syy*svCell->syy + svCell->szz*svCell->szz) +
						  0.25*(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4) +
						  0.25*(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4) +
						  0.25*(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4);

			// store square root of second invariant
			DevStressII = sqrt(DevStressII);

			// Radio of Mohr's cercles
			//R1=1/2*sqrt( (svCell->sxx-svCell->syy)*(svCell->sxx-svCell->syy)) + (XY1+syx)^2)/2.0;
			//R2=1/2*sqrt((Sxx-Szz)^2); %+(sxz+szx)^2);
			//R3=1/2*sqrt((Syy-Szz)^2); %+(sxz+szx)^2);



			// access current pressures
			pc = p[k][j][i];

			// access current pore pressure (zero if deactivated)
			pc_pore  = p_pore[k][j][i];
			pc_hydro  = p_hydro[k][j][i];

			// get total pressure (effective pressure, computed by LaMEM, plus pore pressure)
			//ptotal = pc + biot*(pc_pore-pc_hydro);
			ptotal = pc + biot*pc_pore;

			dP = ptotal - pc_pore; // effective mean stress
			//dP = ptotal - (pc_pore-pc_hydro); // effective mean stress

			//if (k==25 && j==0 && i==127) {
			//	dP = ptotal - pc_pore; // effective mean stress
			//}

			TensileS    = -svBulk->Ts;
			fr = svDev->fr;
			ch = svDev->ch;
			if (ch > 0.0) { // if ch is 0 means that is not yet calculated

				if (fr-1.0 == 0.0)  intersection = 0.0;
				else                intersection = (TensileS-ch)/(fr-1.0);

				//svDev->failT = 0.0;
				//svDev->failS = 0.0;

				// sensitivity
				sensitivity =jr->matLim.stress_sensitivity_for_failure; // Pascals

					if (dP < -(TensileS+sensitivity)) {// + 1e+6) {
						//svDev->failTS = 1.0;
						if (jr->ts.istep > jr->num_of_first_dt) svDev->failT = 1.0;
						frac = 1;
					}
					else
					{
						if (DevStressII > svDev->yield + sensitivity) {
							if (dP>=intersection) {
								// Shear
								if (jr->ts.istep > jr->num_of_first_dt) svDev->failS = 1.0;
								frac = 1;
							}else {
								// Tensile
								if (jr->ts.istep > jr->num_of_first_dt) svDev->failT = 1.0;
								frac = 1;
							}
						}
					}
			}

			// Save pore pressure in the source point, just to analyze data after the simulation
			for(src = 0; src < jr->NumDarcySources; src++)
			{
				if (jr->DarcySources[src].i == i && jr->DarcySources[src].j == j && jr->DarcySources[src].k == k)
				{
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"------------------------------------------ \n");
					sxx = svCell->sxx-ptotal; syy=svCell->syy-ptotal; szz=svCell->szz-ptotal;
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"In point source %lld, Liquid pressure Pl= %g; Effective pressure dP= %g; TAU_II = %g;\n", (LLD)src, pc_pore, dP,DevStressII);

					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"         P=%g; Sxx=%g; Syy=%g; Szz=%g; A=%g;\n", ptotal, sxx,syy,szz,sxx/szz);
					//PetscSynchronizedPrintf(PETSC_COMM_WORLD,"         sxy=%g; sxz=%g; syz=%g;\n", 0.25*(XY1*XY1 + XY2*XY2 + XY3*XY3 + XY4*XY4), 0.25*(XZ1*XZ1 + XZ2*XZ2 + XZ3*XZ3 + XZ4*XZ4), 0.25*(YZ1*YZ1 + YZ2*YZ2 + YZ3*YZ3 + YZ4*YZ4));
					if (frac==1) {
						PetscSynchronizedPrintf(PETSC_COMM_WORLD,"          Fracture !!!!!!!!!!!!!!!! \n");
					}
				}

			}
		}
		END_STD_LOOP

		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

		// restore vectors
		ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &dxy); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &dxz); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &dyz); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dxx); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dyy); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldzz, &dzz); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,   &p);   CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore,   &p_pore); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->hydro_lPl,   &p_hydro); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
