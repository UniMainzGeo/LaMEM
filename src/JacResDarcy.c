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
 * 		div(Kphi/mu grad(P_l)) =  grad(Kphi rhol*g/mu)
 *
 *	where:
 *		mu 	=       liquid viscosity
 *		rhol 	=	liquid density
 *		Kphi 	=	permeability
 *		Pl 	=	liquid pressure
 *		g 	=	gravitational acceleration (acting in the z-direction )
 *
 *		Once Pl is known, the liquid velocity (vl) can be computed from the Darcy velocity:
 *
 *			\phi*(vl - vs) = K_d/mu*(\grad(Pl) - rhol g
 *
 *		vs 	=	solid velocity
 *		vl 	=	liquid velocity
 *		\phi 	=	porosity
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

 typedef struct {
	 PetscScalar Pl;
	 // new//PetscScalar Phi;
} DarcyDOF;

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
						/*#define GET_Kphi0 \
							ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, &kc, NULL, NULL); CHKERRQ(ierr); \
							buff[k][j][i] 	= 	kc;

						// liquid density
						#define GET_rhol \
							ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, &rhol, NULL); CHKERRQ(ierr); \
							buff[k][j][i] 	= 	rhol;

						// liquid viscosity
						#define GET_mul \
							ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, &mul); CHKERRQ(ierr); \
							buff[k][j][i] 	= 	mul;

						// porosity
						#define GET_Phi \
							Phi 			=	jr->svCell[iter++].svBulk.Phi; \
							buff[k][j][i] 	= 	Phi;

						// solid viscosity:
						#define GET_Viscosity \
							eta_creep 		=	jr->svCell[iter++].eta_creep; \
							buff[k][j][i] 	= 	eta_creep;
						*/
// New
#define GET_Kphi \
	kphil = jr->svCell[iter++].svBulk.Kphi; \
	buff[k][j][i] 	= 	kphil;		// get kphi from svBulk or from Material_t??

// liquid viscosity
#define GET_mul \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, &mul, NULL); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	mul;

// liquid viscosity
#define GET_Ssl \
	ierr = JacResGetDarcyParam(jr, jr->svCell[iter++].phRat, NULL, NULL, &Ssl); CHKERRQ(ierr); \
	buff[k][j][i] 	= 	Ssl;

// Pressure at last time step
#define GET_PlOld \
	PlOld 			=	jr->svCell[iter++].svBulk.PlOld; \
	buff[k][j][i] 	= 	PlOld;
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
		PetscScalar *Kphi_,     	// Permeability prefactor
						//PetscScalar *rhol_, 		// Liquid density
		PetscScalar *mu_,			// Liquid viscosity
		PetscScalar *Ss_) 			// New: Specific storage
{
	// compute effective darcy parameters in the cell

	PetscInt    i, numPhases;
    Material_t  *phases, *M;

	PetscScalar Kphi, mu, cf, Kphi_local, Ss; 		//, Kphi_0, rhol

	PetscFunctionBegin;

	// initialize
	Kphi    	= 0.0;
						//rhol    	= 0.0;
	mu      	= 0.0;

	// New
	Ss = 0.0;

	numPhases 	= 	jr->numPhases;
	phases    	= 	jr->phases;

	// average all phases
	for(i = 0; i < numPhases; i++)
	{
		M       = &phases[i];
		cf      =  phRat[i];

		// Compute permeability law
							//Kphi_0     = M->Kphi;			// basic permeability
							//Kphi_local = Kphi_0;	        // local permeability
		// New
		Kphi_local = M->Kphi;	        // local permeability

		// compute average permeability, liquid density and liquid viscosity
							//Kphi += cf*Kphi_local;
							//rhol += cf*M->rhol;
		Kphi += cf*Kphi_local;
		mu   += cf*M->mu;

		// New
		Ss   += cf*M->Ss;
	}

	// store
	if(Kphi_) (*Kphi_) = Kphi;
							//if(rhol_) (*rhol_) = rhol;
	if(mu_)   (*mu_)   = mu;

	// New
	if(Ss_)   (*Ss_)   = Ss;

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

		// check density of the rock phases
										//if(M->rhol  == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define liquid density of phase %lld\n", 	(LLD)i);
		if(M->Kphi  == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define permeability of phase %lld\n", 		(LLD)i);
		if(M->mu    == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define liquid viscosity of phase %lld\n", 	(LLD)i);
		// New
		if(M->Ss    == 0.0)  SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define specific storage of phase %lld\n", 	(LLD)i);
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
	const PetscInt *lx, *ly, *lz; //, dof=2;

	// New
	const dof = 1;

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

	// Set field names
	ierr = DMDASetFieldName(jr->DA_Pl,0,"Pl");	CHKERRQ(ierr);		// liquid pressure
							//ierr = DMDASetFieldName(jr->DA_Pl,1,"Phi");	CHKERRQ(ierr);		// porosity

	// create darcy preconditioner matrix
	ierr = DMCreateMatrix(jr->DA_Pl, &jr->App); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(jr->App, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// LiquidPressure solution vector
	ierr = DMCreateGlobalVector(jr->DA_Pl, &jr->Pl); 		CHKERRQ(ierr);

	// create local LiquidPressure vector using box-stencil central DMDA
	ierr = DMCreateLocalVector(jr->DA_Pl, &jr->lPl); CHKERRQ(ierr);

	// LiquidPressure residual (global)
	ierr = DMCreateGlobalVector(jr->DA_Pl, &jr->r_Pl); 		CHKERRQ(ierr);


	// New
	// create local b vector using box-stencil central DMDA
	ierr = DMCreateLocalVector(jr->DA_Pl, &jr->Darcylb); CHKERRQ(ierr);

	// b (global)
	ierr = DMCreateGlobalVector(jr->DA_Pl, &jr->Darcyb); 		CHKERRQ(ierr);


	// Setup SNES context
	ierr = SNESCreate(PETSC_COMM_WORLD, &jr->Pl_snes);		CHKERRQ(ierr);								// Create SNES
	ierr = SNESSetDM(jr->Pl_snes, jr->DA_Pl);				CHKERRQ(ierr);								// attach DMDA to SNES
	ierr = SNESSetFunction(jr->Pl_snes,jr->r_Pl,		&FormFunction_DARCY, jr); 	CHKERRQ(ierr);		// attach residual funciom
	ierr = SNESSetJacobian(jr->Pl_snes,jr->App,jr->App, &FormJacobian_DARCY, jr);	CHKERRQ(ierr);		// attach picard jacobian (optional)
	ierr = SNESSetOptionsPrefix(jr->Pl_snes,"darcy_");    	CHKERRQ(ierr);
	ierr = SNESSetFromOptions(jr->Pl_snes);					CHKERRQ(ierr);

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

	ierr = VecDestroy(&jr->lPl);   	CHKERRQ(ierr);

	// only if Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	ierr = DMDestroy (&jr->DA_Pl); 		CHKERRQ(ierr);

	// Darcy parameters
	ierr = MatDestroy(&jr->App);   		CHKERRQ(ierr);

	ierr = VecDestroy(&jr->Pl);   		CHKERRQ(ierr);		// solution vector
	ierr = VecDestroy(&jr->r_Pl);   	CHKERRQ(ierr);		// residual vector
	ierr = VecDestroy(&jr->lPl);		CHKERRQ(ierr);		// local solution vector
	ierr = SNESDestroy(&jr->Pl_snes);	CHKERRQ(ierr);

	// New
	ierr = VecDestroy(&jr->Darcylb);	CHKERRQ(ierr);
	ierr = VecDestroy(&jr->Darcyb);		CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitDarcy"
PetscErrorCode JacResInitDarcy(JacRes *jr)
{
	// initialize liquid pressure from markers, or from zero

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscScalar ***lPl, ***bcPl, Pl;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;
	bc = jr->bc;

	ierr = VecZeroEntries(jr->lPl); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(jr->DA_Pl, jr->lPl,  &lPl);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_Pl, bc->bcPl, &bcPl); CHKERRQ(ierr);

	iter = 0;

	ierr = DMDAGetCorners(jr->DA_Pl, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		Pl = bcPl[k][j][i];

		// in the case of Darcy, we do not actually advect Pl on particles, which
		// is why iyt is initialized to zero here
		if(Pl == DBL_MAX) Pl = 0.0	;

		lPl[k][j][i] = Pl;

		iter++;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(jr->DA_Pl, jr->lPl,  &lPl);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_Pl, bc->bcPl, &bcPl); CHKERRQ(ierr);


	// apply two-point constraints
	ierr = JacResApplyDarcyBC(jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResUpdateDarcy"
PetscErrorCode JacResUpdateDarcy(JacRes *jr)
{
	// correct LiquidPressure for diffusion (Newton update)
	SolVarCell 	*svCell;
	SolVarBulk 	*svBulk;
	FDSTAG      *fs;
	DarcyDOF	***sol_local, ***sol_global;
	PetscScalar ***lPl, ***Pl, ***lKphi;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
						//Vec 		local_Kphi0,  	Kphi0_vec;
	// New
	Vec 		local_Kphi,  	Kphi_vec;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	ierr = DMDAVecGetArray(jr->DA_Pl, jr->lPl, &sol_local ); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_Pl,  jr->Pl, &sol_global); 	CHKERRQ(ierr);

	ierr = DMDAGetCorners(jr->DA_Pl, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

								//ExtractCoefficientsFromDA(jr->DA_Pl, "Kphi0",		local_Kphi0,  	Kphi0_vec,  	lKphi);				// Permeability coefficients
								//ExtractCoefficientsFromDA(jr->DA_Pl, "Kphi",		local_Kphi,  	Kphi_vec,  	lKphi);				// Permeability coefficients

	iter = 0;
	START_STD_LOOP
	{
								//svCell 			= 	&jr->svCell[iter++];
								//svBulk 			= 	&svCell->svBulk;

		// update local solution vector
		sol_local[k][j][i].Pl 	=	sol_global[k][j][i].Pl;
								//sol_local[k][j][i].Phi 	=	sol_global[k][j][i].Phi;

								//svBulk->->Kphi			=	lKphi[k][j][i];					// watch out -> this is only the precoefficient!
								//svBulk->Phi				=	sol_global[k][j][i].Phi;		// Update for next dt

	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(jr->DA_Pl,  jr->Pl, &sol_global); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_Pl, jr->lPl, &sol_local); 	CHKERRQ(ierr);

								//RestoreCoefficientsFromDA(jr->DA_Pl, "Kphi0",		local_Kphi0, 	Kphi0_vec, 		lKphi);		// Permeability coefficients

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
	DarcyDOF	***sol;
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
	LOCAL_TO_LOCAL(jr->DA_Pl, jr->lPl)

	// access local solution & boundary constraints
	ierr = DMDAVecGetArray(jr->DA_Pl, jr->lPl,  &sol);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_Pl, bc->bcPl, &bcPl); CHKERRQ(ierr);

	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = sol[k][j][i].Pl;

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC_DARCY(bcPl, sol, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC_DARCY(bcPl, sol, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC_DARCY(bcPl, sol, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC_DARCY(bcPl, sol, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC_DARCY(bcPl, sol, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC_DARCY(bcPl, sol, K, j, i, pmdof) }

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
#define __FUNCT__ "FormFunction_DARCY"
PetscErrorCode FormFunction_DARCY(SNES snes,Vec x,Vec f, JacRes *jr)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// compute Residual
	ierr = JacResGetDarcyRes(snes, x, f, jr); 	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetDarcyRes"
PetscErrorCode JacResGetDarcyRes(SNES snes, Vec x, Vec f, JacRes *jr)
{
	// compute Darcy residual vector

	FDSTAG     		*fs;
	BCCtx      		*bc;
	DM 				da, coordDA;
	Vec				local_sol, local_res, coordinates;
	PetscInt    	iter, num, *list;
	PetscInt    	Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    	i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscInt 		sx_fine,sy_fine,sz_fine;
 	PetscScalar 	bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar 	bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar	 	bqx, fqx, bqy, fqy, bqz, fqz;
	PetscScalar 	brz, frz;
 	PetscScalar 	dx, dy, dz;
 	PetscScalar 	Kphi, mu, Pl_c, gz, *grav; 		// kc,  rhol_param

	PetscScalar 	***mul, ***lKphi; 				//, ***lKphi0, ***rhol, ***buff, *Pl, ***Phi_old,  ***EtaCreep;
	PetscScalar		dt; 							//	K, K0, eta_bulk, phi0,n,m,R, Phi, Phi_bulk;
													//PetscScalar		K, K0, eta_bulk, phi0,n,m,R, Phi,Phi_bulk;
	DarcyDOF 		***sol, ***res;
													//Vec 			Kphi0_vec,local_Kphi0;
	Vec 			local_Kphi;
													//Vec 			PhiOld_local, 	PhiOld_vec;
													//Vec 			rhol_local, 	rhol_vec;
	Vec 			mul_local, 		mul_vec;
													//Vec 			EtaCreep_local, EtaCreep_vec;
	DMDACoor3d     ***coords;

	//New
	Vec 			Kphi_vec;
	Vec 			Ss_local,  Ss_vec;
	PetscScalar		***Ssl;											// specific storage

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Get DA for the CURRENT level (not necessarily the fine grid!)
	ierr 		=	SNESGetDM(snes,&da); 				CHKERRQ(ierr);

	// Get local vectors attached with the solution at the local DA
	ierr 		= 	DMGetLocalVector(da,&local_sol); 	CHKERRQ(ierr);
	ierr 		= 	DMGetLocalVector(da,&local_res); 	CHKERRQ(ierr);
						//ierr 		= 	DMGetLocalVector(da,&local_Kphi); 	CHKERRQ(ierr);

	// access residual context variables
	fs        = jr->fs;
	grav      = jr->grav;       // gravity acceleration
	gz 		  =	grav[2];		// gravitational acceleration in z-direction
	bc        = jr->bc;
	dt        = jr->ts.dt;      // time step
	num       = bc->Pl_NumSPC;
	list      = bc->Pl_SPCList;

	// Dimensions of the global grid at the current level
	DMDAGetInfo(da, 0, &nx, &ny,&nz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE, PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
	mx 			= 	nx-1;
	my 			= 	ny-1;
	mz 			= 	nz-1;

	// Extract required material parameters
							//ExtractCoefficientsFromDA(da, "Kphi0",		local_Kphi0,  	Kphi0_vec,  	lKphi0);		// Permeability K0
							//ExtractCoefficientsFromDA(da, "rhol",		rhol_local, 	rhol_vec, 		rhol);				// Liquid density
							//ExtractCoefficientsFromDA(da, "Phi_old",	PhiOld_local, 	PhiOld_vec, 	Phi_old);			// Porosity at last timestep
	ExtractCoefficientsFromDA(da, "mul",	mul_local, 		mul_vec, 	mul);				// Liquid viscosity
							//ExtractCoefficientsFromDA(da, "eta_creep",	EtaCreep_local, EtaCreep_vec, 	EtaCreep);			// Solid viscosity
	// New
	ExtractCoefficientsFromDA(da,"Kphi",	local_Kphi,  	Kphi_vec,  	lKphi);		// Permeability K
	//ExtractCoefficientsFromDA(da, "Pl_old",	PlOld_local, 	PlOld_vec, 	Pl_old);
	ExtractCoefficientsFromDA(da, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);

							//ierr = DMGlobalToLocalBegin(da, Kphi0_vec,INSERT_VALUES, local_Kphi); CHKERRQ(ierr);
							//ierr = DMGlobalToLocalEnd  (da, Kphi0_vec,INSERT_VALUES, local_Kphi); CHKERRQ(ierr);

	// Copy global solution vector into local  vector
	ierr = DMGlobalToLocalBegin(da, x,INSERT_VALUES, local_sol); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (da, x,INSERT_VALUES, local_sol); CHKERRQ(ierr);

	// Copy global residual into local  vector
	ierr = DMGlobalToLocalBegin(da, f,INSERT_VALUES, local_res); 	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (da, f,INSERT_VALUES, local_res); 	CHKERRQ(ierr);

	// Make coordinates available
	ierr = DMGetCoordinateDM(da, &coordDA);							CHKERRQ(ierr);		// coordinates DA
	ierr = DMGetCoordinatesLocal(da, &coordinates);					CHKERRQ(ierr);		// vector with coordinates
	ierr = DMDAVecGetArray(coordDA, coordinates, &coords);			CHKERRQ(ierr);		// array with coordinates

	// access work vectors
	ierr = DMDAVecGetArray(da, local_res,  	&res); 					CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_sol,  	&sol);  				CHKERRQ(ierr);
							//ierr = DMDAVecGetArray(da, local_Kphi,  &lKphi);  				CHKERRQ(ierr);

	// Update ghost points for the correct boundary values
	ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
	ierr = DMDAGetCorners(jr->DA_Pl, &sx_fine, &sy_fine, &sz_fine, NULL,NULL,NULL); CHKERRQ(ierr);


	// Put this in an other place (Vec x is locked read only)////////////////////////////////////////////////////////////
	/*START_STD_LOOP
	{
		if (k==0){
			sol[k-1][j][i].Pl =  2*bc->Pl_bot - sol[k][j][i].Pl;    // bottom boundary, dirichlet value constant
		}
		if (k==mz){
			sol[k+1][j][i].Pl =  2*bc->Pl_top - sol[k][j][i].Pl  ;	// top boundary, dirichlet value constant
		}

		if (j==0){
			sol[k][j-1][i].Pl =	sol[k][j][i].Pl;		// zero flux
		}
		if (j==my){
			sol[k][j+1][i].Pl = sol[k][j][i].Pl;		// zero flux
		}
		if (i==0){
			sol[k][j][i-1].Pl = sol[k][j][i].Pl;		// zero flux
		}
		if (i==mx){
			sol[k][j][i+1].Pl = sol[k][j][i].Pl;		// zero flux
		}
	}
	END_STD_LOOP*/
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

							/*//---------------
							// Update (nonlinear) material properties @ all points (including ghost points)
							//---------------
							ierr = DMDAGetGhostCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

							phi0 			= 		0.01;
							m 				=		1.0;
							n 				=		3.0;
							R 				=		1.0;	// decompaction_weakening
							PetscOptionsGetScalar(PETSC_NULL, PETSC_NULL,"-PorosityWaves_R",&R,PETSC_NULL);
							START_STD_LOOP
							{
								K0 				=		lKphi0[k][j][i];
								Phi 			=		sol[k][j][i].Phi;
								lKphi[k][j][i] 	= 		K0*PetscPowScalar((Phi/phi0),n);

							//	PetscPrintf(PETSC_COMM_WORLD,"eta_bulk=%f  phi=%e K0=%e KPhi=%e Phi_old=%e\n",eta_bulk, Phi, K0, lKphi[k][j][i],Phi_old[k][j][i] );
							//	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Stop here");
							}
							END_STD_LOOP
							*/

	//---------------
	// central points
	//---------------
	ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
	START_STD_LOOP
	{
		// access solution vector
		Pl_c  = sol[k][j][i].Pl; // current liquid pressure

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// compute average permeabilities at edges from the defined values @ centers
		Kphi 	=	 lKphi[k][j][i];
		bkx 	= 	(Kphi + lKphi[k][j][Im1])/2.0;       fkx = (Kphi + lKphi[k][j][Ip1])/2.0;
		bky 	= 	(Kphi + lKphi[k][Jm1][i])/2.0;       fky = (Kphi + lKphi[k][Jp1][i])/2.0;
		bkz 	= 	(Kphi + lKphi[Km1][j][i])/2.0;       fkz = (Kphi + lKphi[Kp1][j][i])/2.0;

							/*
							// compute liquid densities at the necessary points, from their definition at center points
							brz 	= (rhol[k][j][i] + rhol[Km1][j][i])/2.0; frz = (rhol[k][j][i] + rhol[Kp1][j][i])/2.0;

							// compute liquid densities at the necessary points, from their definition at center points
							brz 	= (rhol[k][j][i] + rhol[Km1][j][i])/2.0; frz = (rhol[k][j][i] + rhol[Kp1][j][i])/2.0;
							*/

	    // get mesh steps
		bdx 	= 	coords[k  ][j  ][i  ].x - coords[k  ][j  ][i-1].x;
		fdx 	= 	coords[k  ][j  ][i+1].x - coords[k  ][j  ][i  ].x;
		bdy 	= 	coords[k  ][j  ][i  ].y - coords[k  ][j-1][i  ].y;
		fdy 	= 	coords[k  ][j+1][i  ].y - coords[k  ][j  ][i  ].y;
		bdz 	= 	coords[k  ][j  ][i  ].z - coords[k-1][j  ][i  ].z;
		fdz 	= 	coords[k+1][j  ][i  ].z - coords[k  ][j  ][i  ].z;

		// deal with boundaries of the domain
		if (i==0 ) bdx= fdx;	if (i==mx) fdx= bdx;
		if (j==0 ) bdy= fdy;	if (j==my) fdy= bdy;
		if (k==0 ) bdz= fdz;	if (k==mz) fdz= bdz;

		//	extract liquid viscosity
		mu		=	mul[k][j][i];		// liquid viscosity

		// compute LiquidPressure fluxes
		bqx 	= bkx/mu*(Pl_c - sol[k][j][i-1].Pl)/bdx;   fqx = fkx/mu*(sol[k][j][i+1].Pl - Pl_c)/fdx;
		bqy 	= bky/mu*(Pl_c - sol[k][j-1][i].Pl)/bdy;   fqy = fky/mu*(sol[k][j+1][i].Pl - Pl_c)/fdy;
		bqz 	= bkz/mu*(Pl_c - sol[k-1][j][i].Pl)/bdz;   fqz = fkz/mu*(sol[k+1][j][i].Pl - Pl_c)/fdz;

		// get mesh steps
		dx 	= (bdx+fdx)/2.0;
		dy 	= (bdy+fdy)/2.0;
		dz 	= (bdz+fdz)/2.0;

							/*// original balance equation:
							Phi_bulk 		=	sol[k][j][i].Phi;
							if (Phi_bulk<1e-6) Phi_bulk = 1e-6;
							eta_bulk 		= 	EtaCreep[k][j][i]/(Phi_bulk);


							// correct bulk viscosity for decompaction weakening
							if 		(Pl_c< 0){
								eta_bulk = eta_bulk*R;
							}

							if ((i==0) && (j==0)){
									//			PetscPrintf(PETSC_COMM_SELF,"i,j,k=[%i,%i,%i] bkz=%e, fkz=%e\n",i,j,k,bkz,fkz);
							}
							// (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz = (frz*g - brz*g)/dz


							// to get positive diagonal in the preconditioner matrix
							// put right hand side to the left, which gives the following:

							//g_Pl[k][j][i] = (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz ;
							res[k][j][i].Pl 	=  Pl_c/eta_bulk - (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz + 1*(fkz/mu*frz*gz - bkz/mu*brz*gz)/dz;

							// Form residual for porosity
							Phi 				=	sol[k][j][i].Phi;
							res[k][j][i].Phi 	=	(Phi-Phi_old[k][j][i])/dt + (1-Phi)*Pl_c/eta_bulk;
							*/


		// new
		// Residual for Pl
		//res[k][j][i].Pl = Ss[k][j][i]*(Pl_c-Pl_old[k][j][i])/dt- (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz; // + 1*(fkz/mu*frz*gz - bkz/mu*brz*gz)/dz;
		res[k][j][i].Pl = Ssl[k][j][i]*Pl_c/dt- (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz;


	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(da, local_res, 	&res); 			CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, local_sol,   &sol);  		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);	CHKERRQ(ierr);	// array with coordinates
						//ierr = DMDAVecRestoreArray(da, local_Kphi,  &lKphi);  		CHKERRQ(ierr);


	// restore material coefficients on DA
						//RestoreCoefficientsFromDA(da, "Kphi0",	local_Kphi0, 	Kphi0_vec, 		lKphi0);	// Permeability coefficients
						//RestoreCoefficientsFromDA(da, "rhol",		rhol_local, 	rhol_vec, 		rhol);		// Liquid density
						//RestoreCoefficientsFromDA(da, "Phi_old",	PhiOld_local, 	PhiOld_vec, 	Phi_old);	// Porosity at last timestep
	RestoreCoefficientsFromDA(da, "mul",	mul_local, 		mul_vec, 	mul);		// Liquid viscosity
						//RestoreCoefficientsFromDA(da, "eta_creep",	EtaCreep_local, EtaCreep_vec, 	EtaCreep);	// Solid viscosity

	//New
	//RestoreCoefficientsFromDA(da, "Pl_old",	PlOld_local, 	PlOld_vec, 	Pl_old);	// Liquid pressure at last timestep
	RestoreCoefficientsFromDA(da, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);		// Specific storage
	RestoreCoefficientsFromDA(da, "Kphi",	local_Kphi, 	Kphi_vec, 	lKphi);		// Permeability coefficients


	// copy back local vectors to global ones
	ierr = DMLocalToGlobalBegin(da, local_res,INSERT_VALUES, f); 	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (da, local_res,INSERT_VALUES, f); 	CHKERRQ(ierr);

	// We need that??? (Vec x is locked read only)
	//ierr = DMLocalToGlobalBegin(da, local_sol,INSERT_VALUES, x); 	CHKERRQ(ierr);
	//ierr = DMLocalToGlobalEnd  (da, local_sol,INSERT_VALUES, x); 	CHKERRQ(ierr);
	//////////////////////////////////////////

	ierr= DMRestoreLocalVector (da, &local_Kphi); 					CHKERRQ(ierr);
	ierr= DMRestoreLocalVector (da,	&local_sol); 					CHKERRQ(ierr);
	ierr= DMRestoreLocalVector (da,	&local_res); 					CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormJacobian_DARCY"
PetscErrorCode FormJacobian_DARCY(SNES snes,Vec x, Mat P, Mat J, JacRes *jr)
{
	// calls routine to create the jacobian

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// compute jacobian
	ierr = JacResGetDarcyMat(snes, x, jr); 		CHKERRQ(ierr);

	// set back residual
	MatDuplicate(jr->App,MAT_COPY_VALUES, &P);
	MatDuplicate(jr->App,MAT_COPY_VALUES, &J);

	{
		PetscViewer view;
		 PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Matrixes.bin",FILE_MODE_WRITE,&view);
		 MatView(P, view);
		 PetscViewerDestroy(&view);
	}


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetDarcyMat"
PetscErrorCode JacResGetDarcyMat(SNES snes, Vec x, JacRes *jr)
{
	// assemble Darcy preconditioner matrix
	// COMPLETE SINGLE-POINT CONSTRAINT IMPLEMENTATION !!!

	FDSTAG     *fs;
	BCCtx      *bc;
	DM 			da, coordDA;
	Vec			local_sol, local_res, coordinates;
	SolVarCell *svCell;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;
	PetscScalar v[7], cf[6], mu,  Kphi, Ss;	//, kc, rhol
	MatStencil  row[1], col[7];
	PetscScalar ***lKphi, ***bcPl, ***mul;	//,***lrhol, ***buff

	PetscScalar		n,m,R; 	//K, K0, eta_bulk, phi0,Phi_bulk, Phi;
	DarcyDOF 		***sol;
									//Vec 			Kphi0_vec,		local_Kphi0;
	Vec 			local_Kphi;
	Vec 			PhiOld_local, 	PhiOld_vec;
									//Vec 			rhol_local, 	rhol_vec;
	Vec 			mul_local, 		mul_vec;
									//Vec 			EtaCreep_local, EtaCreep_vec;
	DMDACoor3d 		***coords;

	PetscScalar dt;

	//New
	Vec 			Kphi_vec;
	Vec 			Ss_local,  Ss_vec;
	PetscScalar		***Ssl;											// specific storage


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs        = jr->fs;
	bc        = jr->bc;
	num       = bc->Pl_NumSPC;
	list      = bc->Pl_SPCList;
	dt        = jr->ts.dt;     // time step

	// Get DA for the CURRENT level (not necessarily the fine grid!)
	ierr 		=	SNESGetDM(snes,&da); 				CHKERRQ(ierr);

	// Get local vectors attached with the solution at the local DA
	ierr 		= 	DMGetLocalVector(da,&local_sol); 	CHKERRQ(ierr);
	ierr 		= 	DMGetLocalVector(da,&local_res); 	CHKERRQ(ierr);
							//ierr 		= 	DMGetLocalVector(da,&local_Kphi); 	CHKERRQ(ierr);

	// Dimensions of the global grid at the current level
	DMDAGetInfo(da, 0, &nx, &ny,&nz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE, PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
	mx 			= 	nx-1;
	my 			= 	ny-1;
	mz 			= 	nz-1;

	// clear matrix coefficients
	ierr = MatZeroEntries(jr->App); CHKERRQ(ierr);

	// Extract required material parameters
							//ExtractCoefficientsFromDA(da, "Kphi0",		local_Kphi0,  	Kphi0_vec,  	lKphi0);			// Permeability K0
							//ExtractCoefficientsFromDA(da, "rhol",		rhol_local, 	rhol_vec, 		rhol);				// Liquid density
							//ExtractCoefficientsFromDA(da, "Phi_old",	PhiOld_local, 	PhiOld_vec, 	Phi_old);			// Porosity at last timestep
	ExtractCoefficientsFromDA(da, "mul"		,	mul_local, 		mul_vec, 		mul);				// Liquid viscosity
							//ExtractCoefficientsFromDA(da, "eta_creep",	EtaCreep_local, EtaCreep_vec, 	EtaCreep);			// Solid viscosity

	// New
	ExtractCoefficientsFromDA(da,"Kphi",	local_Kphi,  	Kphi_vec,  	lKphi);		// Permeability K
	ExtractCoefficientsFromDA(da, "Ssl",	Ss_local, 		Ss_vec, 	Ssl);

							//ierr = DMGlobalToLocalBegin(da, Kphi0_vec,INSERT_VALUES, local_Kphi); CHKERRQ(ierr);
							//ierr = DMGlobalToLocalEnd  (da, Kphi0_vec,INSERT_VALUES, local_Kphi); CHKERRQ(ierr);

	// Copy global solution vector into local  vector
	ierr = DMGlobalToLocalBegin(da, x,INSERT_VALUES, local_sol); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (da, x,INSERT_VALUES, local_sol); CHKERRQ(ierr);

	// Copy global residual into local  vector
//	ierr = DMGlobalToLocalBegin(da, f,INSERT_VALUES, local_res); 	CHKERRQ(ierr);
//	ierr = DMGlobalToLocalEnd  (da, f,INSERT_VALUES, local_res); 	CHKERRQ(ierr);

	// Make coordinates available
	ierr = DMGetCoordinateDM(da, &coordDA);							CHKERRQ(ierr);		// coordinates DA
	ierr = DMGetCoordinatesLocal(da, &coordinates);					CHKERRQ(ierr);		// vector with coordinates
	ierr = DMDAVecGetArray(coordDA, coordinates, &coords);			CHKERRQ(ierr);		// array with coordinates

	// access work vectors
//	ierr = DMDAVecGetArray(da, local_res,  	&res); 					CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_sol,  	&sol);  				CHKERRQ(ierr);
							//ierr = DMDAVecGetArray(da, local_Kphi,  &lKphi);  				CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_Pl, bc->bcPl, &bcPl); 	CHKERRQ(ierr);

								/*
								// Update ghost points for the correct boundary values
								ierr = DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
								//ierr = DMDAGetCorners(jr->DA_Pl, &sx_fine, &sy_fine, &sz_fine, NULL,NULL,NULL); CHKERRQ(ierr);

								//---------------
								// Update (nonlinear) material properties @ all points (including ghost points)
								//---------------
								ierr 			= 	DMDAGetGhostCorners(da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
								phi0 			= 	0.01;
								m 				=	1.0;
								n 				=	3.0;
								R 				=	1.0;	// decompaction_weakening
								PetscOptionsGetScalar(PETSC_NULL,PETSC_NULL,"-PorosityWaves_R",&R,PETSC_NULL);
								START_STD_LOOP
								{
									K0 				=		lKphi0[k][j][i];
									Phi 			=		sol[k][j][i].Phi;
									lKphi[k][j][i] 	= 		K0*PetscPowScalar((Phi/phi0),n);

									//	PetscPrintf(PETSC_COMM_WORLD,"eta_bulk=%f  phi=%e K0=%e KPhi=%e Phi_old=%e\n",eta_bulk, Phi, K0, lKphi[k][j][i],Phi_old[k][j][i] );
									//	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Stop here");
								}
								END_STD_LOOP
								*/
	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(jr->DA_Pl, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables @ fine grid
		svCell = &jr->svCell[iter++];

		// permeability, liquid density and viscosity
		// and specific storage
		ierr = JacResGetDarcyParam(jr, svCell->phRat, &Kphi, &mu, &Ss); CHKERRQ(ierr);

		// check index bounds and TPC multipliers
		/*
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcPl[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }	// need to have bc set!!
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcPl[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcPl[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcPl[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcPl[k-1][j][i] != DBL_MAX) cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcPl[k+1][j][i] != DBL_MAX) cf[5] = -1.0; }
		 */

		// hardcoded dirichlet Top/Bottom
		//Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcPl[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }	// need to have bc set!!
		//Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcPl[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		//Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcPl[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		//Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcPl[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		//Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++;  cf[4] = -1.0; }
		//Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--;  cf[5] = -1.0; }

		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++;  }	// need to have bc set!!
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--;  }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++;  }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--;  }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++;  cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--;  cf[5] = -1.0; }

		// compute average permeabilities at edges from the defined values @ centers
				Kphi 	=	 lKphi[k][j][i];
				bkx 	= 	(Kphi + lKphi[k][j][Im1])/2.0;       fkx = (Kphi + lKphi[k][j][Ip1])/2.0;
				bky 	= 	(Kphi + lKphi[k][Jm1][i])/2.0;       fky = (Kphi + lKphi[k][Jp1][i])/2.0;
				bkz 	= 	(Kphi + lKphi[Km1][j][i])/2.0;       fkz = (Kphi + lKphi[Kp1][j][i])/2.0;

				// we need K/mu in the governing equation
				mu		=	mul[k][j][i];		// liquid viscosity
				bkx /= mu;							fkx /= mu;
				bky /= mu;							fky /= mu;
				bkz /= mu;							fkz /= mu;

				// get mesh steps
				// if we had access to an accordingly coarsened 1D mesh object, it would be preferable
				bdx 	= 	coords[k  ][j  ][i  ].x - coords[k  ][j  ][i-1].x;
				fdx 	= 	coords[k  ][j  ][i+1].x - coords[k  ][j  ][i  ].x;
				bdy 	= 	coords[k  ][j  ][i  ].y - coords[k  ][j-1][i  ].y;
				fdy 	= 	coords[k  ][j+1][i  ].y - coords[k  ][j  ][i  ].y;
				bdz 	= 	coords[k  ][j  ][i  ].z - coords[k-1][j  ][i  ].z;
				fdz 	= 	coords[k+1][j  ][i  ].z - coords[k  ][j  ][i  ].z;

				// deal with boundaries of the domain
				if (i==0 ) bdx= fdx;	if (i==mx) fdx= bdx;
				if (j==0 ) bdy= fdy;	if (j==my) fdy= bdy;
				if (k==0 ) bdz= fdz;	if (k==mz) fdz= bdz;

				// get mesh steps
				dx 		= 	(bdx+fdx)/2.0;
				dy 		= 	(bdy+fdy)/2.0;
				dz 		= 	(bdz+fdz)/2.0;

												/*// original balance equation:
																	//Phi_bulk =	sol[k][j][i].Phi;
												if (Phi_bulk<1e-6) Phi_bulk = 1e-6;
																	//eta_bulk = 	EtaCreep[k][j][i]/(Phi_bulk);

												if 		(sol[k][j][i].Pl< 0){
													eta_bulk = eta_bulk*R;
												}
												*/

				//==================================================================
				// Fluid pressure
				//==================================================================
				// set row/column indices for the fluid pressure equations
				row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
				col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
				col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
				col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
				col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
				col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
				col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
				col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

				Ss		=	Ssl[k][j][i];		// specific storage

				// set values including TPC multipliers
				v[0] =  -bkx/bdx/dx*cf[0];
				v[1] =  -fkx/fdx/dx*cf[1];
				v[2] =  -bky/bdy/dy*cf[2];
				v[3] =  -fky/fdy/dy*cf[3];
				v[4] =  -bkz/bdz/dz*cf[4];
				v[5] =  -fkz/fdz/dz*cf[5];
											//v[6] =  1.0/eta_bulk
				v[6] =  Ss/dt
				+       (bkx/bdx + fkx/fdx)/dx
				+       (bky/bdy + fky/fdy)/dy
				+       (bkz/bdz + fkz/fdz)/dz;

														//// Coupling term with porosity
														//col[7].k = k;   col[7].j = j;   col[7].i = i;   col[7].c = 1;

														//// THIS MIGHT BE A MORE COMPLICATED COUPLING ACTUALLY..
														//v[7] = sol[k][j][i].Pl/eta_bulk;


				// set matrix coefficients
														//ierr = MatSetValuesStencil(jr->App, 1, row, 8, col, v, ADD_VALUES); CHKERRQ(ierr);
				ierr = MatSetValuesStencil(jr->App, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

														/*
														//==================================================================
														// Porosity
														//==================================================================
														row[0].k 	= 	k;   row[0].j = j;   row[0].i = i;   row[0].c = 1;

														col[0].k 	= 	k;   col[0].j = j;   col[0].i = i;   col[0].c = 1;
														v[0] 		=  	1/dt;												// diagonal porosity term

														// we presumably need to incorporate the coupling terms as well...
														col[1].k 	= 	k;   col[1].j = j;   col[1].i = i;   col[1].c = 0;		// Pf[k,j,i]
																				//v[1] 		=  	 (1-sol[k][j][i].Phi)/eta_bulk;											// diagonal porosity term

														// set matrix coefficients
														ierr = MatSetValuesStencil(jr->App, 1, row, 2, col, v, ADD_VALUES); CHKERRQ(ierr);
														*/

				//res[k][j][i].Phi 	=	(Phi-Phi_old[k][j][i])/dt + (1-Phi)*Pl_c/eta_bulk;

				// NOTE! since only TPC are active, no SPC modification is necessary
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(da, local_sol,   &sol);  		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);	CHKERRQ(ierr);	// array with coordinates
	ierr = DMDAVecRestoreArray(da, local_Kphi,  &lKphi);  		CHKERRQ(ierr);

	// restore material coefficients on DA
								//RestoreCoefficientsFromDA(da, "Kphi0",		local_Kphi0, 	Kphi0_vec, 		lKphi0);	// Permeability coefficients
								//RestoreCoefficientsFromDA(da, "rhol",		rhol_local, 	rhol_vec, 		rhol);		// Liquid density
								//RestoreCoefficientsFromDA(da, "Phi_old",	PhiOld_local, 	PhiOld_vec, 	Phi_old);	// Porosity at last timestep
	RestoreCoefficientsFromDA(da, "mul"		,	mul_local, 		mul_vec, 		mul);		// Liquid viscosity
								//RestoreCoefficientsFromDA(da, "eta_creep",	EtaCreep_local, EtaCreep_vec, 	EtaCreep);	// Solid viscosity


	// New
	RestoreCoefficientsFromDA(da, "Kphi",	local_Kphi, 	Kphi_vec, 		lKphi);		// Permeability coefficients
	RestoreCoefficientsFromDA(da, "Ssl" ,	Ss_local, 		Ss_vec, 		Ssl);		// specific storage

	// copy back local vectors to global ones
	ierr = DMLocalToGlobalBegin(da, local_sol,INSERT_VALUES, x); 	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (da, local_sol,INSERT_VALUES, x); 	CHKERRQ(ierr);

								//ierr= DMRestoreLocalVector (da, &local_Kphi); 					CHKERRQ(ierr);
	ierr= DMRestoreLocalVector (da,	&local_sol); 					CHKERRQ(ierr);
	ierr= DMRestoreLocalVector (da,	&local_res); 					CHKERRQ(ierr);


	// restore access
	//ierr = DMDAVecRestoreArray(jr->DA_Pl, bc->bcPl, &bcPl);   CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(jr->DA_Pl, jr->ldyy, &lrhol);  CHKERRQ(ierr);
	//ierr = DMDAVecRestoreArray(jr->DA_Pl, jr->ldxx, &lKphi);    CHKERRQ(ierr);

	// assemble LiquidPressure/Darcy matrix
	ierr = MatAIJAssemble(jr->App, num, list, 1.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
Diffusion term expansion

		bqx = bkx*(Tc - T[k][j][i-1])/bdx;   fqx = fkx*(T[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - T[k][j-1][i])/bdy;   fqy = fky*(T[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - T[k-1][j][i])/bdz;   fqz = fkz*(T[k+1][j][i] - Tc)/fdz;

[1]
		-(fqx - bqx)/dx
		-(fqy - bqy)/dy
		-(fqz - bqz)/dz
[2]
		(bqx - fqx)/dx
		(bqy - fqy)/dy
		(bqz - fqz)/dz
[3]
		(bkx*(Tc - T[k][j][i-1])/bdx - fkx*(T[k][j][i+1] - Tc)/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy - fky*(T[k][j+1][i] - Tc)/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz - fkz*(T[k+1][j][i] - Tc)/fdz)/dz
[4]
		(bkx*(Tc - T[k][j][i-1])/bdx + fkx*(Tc - T[k][j][i+1])/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy + fky*(Tc - T[k][j+1][i])/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz + fkz*(Tc - T[k+1][j][i])/fdz)/dz
[5]
		bkx/bdx/dx*(Tc - T[k][j][i-1]) + fkx/fdx/dx*(Tc - T[k][j][i+1])
		bky/bdy/dy*(Tc - T[k][j-1][i]) + fky/fdy/dy*(Tc - T[k][j+1][i])
		bkz/bdz/dz*(Tc - T[k-1][j][i]) + fkz/fdz/dz*(Tc - T[k+1][j][i])
[6]
		(bkx/bdx/dx + fkx/fdx/dx)*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy/dy + fky/fdy/dy)*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz/dz + fkz/fdz/dz)*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
[7]

		(bkx/bdx + fkx/fdx)/dx*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy + fky/fdy)/dy*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz + fkz/fdz)/dz*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
*/
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCCreateDarcy"
PetscErrorCode BCCreateDarcy(JacRes *jr, BCCtx *bc)
{

	PetscFunctionBegin;

	PetscErrorCode 	ierr;
	FDSTAG 			*fs;

	// Only for cases where Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	// create BC context for DARCY
	fs = bc->fs;

	// create boundary conditions vectors (velocity, pressure, temperature, liquid pressure)
	ierr = DMCreateLocalVector(jr->DA_Pl,  &bc->bcPl);  CHKERRQ(ierr);

	// SPC (LiquidPressure/Darcy)
	ierr = makeIntArray (&bc->Pl_SPCList, NULL, fs->dof.lnp); CHKERRQ(ierr);
	ierr = makeScalArray(&bc->Pl_SPCVals, NULL, fs->dof.lnp); CHKERRQ(ierr);


	// Mark as unconstrained:
	ierr = VecSet(bc->bcPl, DBL_MAX); CHKERRQ(ierr);

	bc->Pl_NumSPC 	= 0;

	bc->fs   = fs;

	// exchange ghost point constraints
	// AVOID THIS BY SETTING CONSTRAINTS REDUNDANTLY
	// IN MULTIGRID ONLY REPEAT BC COARSENING WHEN THINGS CHANGE
	LOCAL_TO_LOCAL(jr->DA_Pl, bc->bcPl)


	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BCApplyBound_DARCY"
PetscErrorCode BCApplyBound_DARCY(BCCtx *bc, JacRes *jr)
{
	// apply boundary conditions for DARCY
	FDSTAG      *fs;
	PetscScalar Pl_bot, Pl_top;
	PetscInt    mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcPl;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// Only for cases where Darcy is active
	if(jr->actDarcy != PETSC_TRUE) PetscFunctionReturn(0);

	// access context
	fs = bc->fs;

	// initialize maximal index in z-direction
	mcz = fs->dsz.tcels - 1;

	//-----------------------------------------------------
	// LiquidPressure points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	Pl_bot 	=	bc->Pl_bot;
	Pl_top 	=	bc->Pl_top;
	ierr = DMDAVecGetArray(jr->DA_Pl, bc->bcPl,  &bcPl);  CHKERRQ(ierr);
	if(Pl_bot >= 0.0 || Pl_top >= 0.0)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// only positive or zero fluid pressure BC's!
			// negative will set zero-flux BC automatically
			if(Pl_bot >= 0.0 && k == 0)   { bcPl[k-1][j][i] = Pl_bot; }
			if(Pl_top >= 0.0 && k == mcz) { bcPl[k+1][j][i] = Pl_top; }
		}
		END_STD_LOOP
	}
	ierr = DMDAVecRestoreArray(jr->DA_Pl, bc->bcPl,  &bcPl);  CHKERRQ(ierr);

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
	//DM             da    = jr->DA_Pl;		// get's Darcy CA
	PetscInt 	   	iter, sx, sy, sz;
	PetscInt       	i, j, k, nx, ny, nz;
						//Vec            	local_rhol, 	rhol_vec;
						//Vec 			Kphi0_vec, 		local_KPhi0;
	Vec 			Kphi_vec, 		local_KPhi;
						//Vec 			Phi_vec, 		local_Phi;								// porosity of last timestep
	Vec 			Mul_vec, 		local_Mul;								// liquid viscosity
						//Vec 			EtaCreep_vec, 	local_EtaCreep;							// solid (creep) viscosity
						//Vec 			local_sol;
	PetscScalar    	***buff, kphil, mul; // kc, Phi, rhol, eta_creep;
					//DarcyDOF		***sol;
	//New
	PetscScalar Ssl, PlOld;
	Vec 			PlOld_vec, local_PlOld;
	Vec 			Ss_vec,	local_Ss;



	//---------------------------------------------
	// (1) Update coordinates on 3D DMDA object
	ierr =  UpdateCoordinatesOnCellCenterDA(jr->DA_Pl, jr->fs);	CHKERRQ(ierr);
	//---------------------------------------------

	//---------------------------------------------
	// (2) Update material properties on local vector & attach to DM

	// Create temporary local vectors
									//ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_KPhi0); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_KPhi); 				CHKERRQ(ierr);
									//ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_rhol); 				CHKERRQ(ierr);
									//ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Phi); 				CHKERRQ(ierr);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Mul); 				CHKERRQ(ierr);
									//ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_EtaCreep); 			CHKERRQ(ierr);
	// New
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_Ss);
	ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_PlOld);

	// Set data on local vectors
									//SCATTER_FIELD(jr->DA_Pl, local_KPhi0, GET_Kphi0);			// 3D array of permeability pre-coefficient
									//SCATTER_FIELD(jr->DA_Pl, local_rhol,  GET_rhol);			// compute 3D array of liquid density
									//SCATTER_FIELD(jr->DA_Pl, local_Phi,   GET_Phi);				// compute 3D array of porosity (of last timestep)
	SCATTER_FIELD(jr->DA_Pl, local_Mul,   GET_mul);				// compute 3D array of liquid viscosity
									//SCATTER_FIELD(jr->DA_Pl, local_EtaCreep,  GET_Viscosity);	// compute 3D array of solid viscosity
	//New
	SCATTER_FIELD(jr->DA_Pl, local_KPhi, GET_Kphi);				// 3D array of permeability pre-coefficient
	SCATTER_FIELD(jr->DA_Pl, local_Ss,   GET_Ssl);
	SCATTER_FIELD(jr->DA_Pl, local_PlOld,GET_PlOld);

	// Get/Create global vector, attached to DM
									//ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Kphi0",		&Kphi0_vec); 			CHKERRQ(ierr);
									//ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"rhol",		&rhol_vec); 			CHKERRQ(ierr);
									//ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Phi_old",		&Phi_vec); 				CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"mul",			&Mul_vec); 				CHKERRQ(ierr);
									//ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"eta_creep",	&EtaCreep_vec); 		CHKERRQ(ierr);
	//New
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Kphi",		&Kphi_vec); 			CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Ssl",			&Ss_vec); 				CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(jr->DA_Pl,"Pl_old",		&PlOld_vec); 			CHKERRQ(ierr);

	// Move local values to global vector
									/*
									ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_KPhi0,INSERT_VALUES ,Kphi0_vec);	CHKERRQ(ierr);
									ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_KPhi0,INSERT_VALUES ,Kphi0_vec);	CHKERRQ(ierr);

									ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_rhol,INSERT_VALUES ,rhol_vec);		CHKERRQ(ierr);
									ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_rhol,INSERT_VALUES ,rhol_vec);		CHKERRQ(ierr);

									ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Phi,INSERT_VALUES ,Phi_vec);		CHKERRQ(ierr);
									ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Phi,INSERT_VALUES ,Phi_vec);		CHKERRQ(ierr);
									*/

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Mul,INSERT_VALUES ,Mul_vec);		CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Mul,INSERT_VALUES ,Mul_vec);		CHKERRQ(ierr);

	// New
	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_KPhi,INSERT_VALUES ,Kphi_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_KPhi,INSERT_VALUES ,Kphi_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_Ss,INSERT_VALUES ,Ss_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_Ss,INSERT_VALUES ,Ss_vec);	CHKERRQ(ierr);

	ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_PlOld,INSERT_VALUES ,PlOld_vec);	CHKERRQ(ierr);
	ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_PlOld,INSERT_VALUES ,PlOld_vec);	CHKERRQ(ierr);

									/*
									ierr = 	DMLocalToGlobalBegin(jr->DA_Pl,local_EtaCreep,INSERT_VALUES ,EtaCreep_vec);		CHKERRQ(ierr);
									ierr = 	DMLocalToGlobalEnd  (jr->DA_Pl,local_EtaCreep,INSERT_VALUES ,EtaCreep_vec);		CHKERRQ(ierr);
									*/

	// attach the named vector to the DM

									//ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Kphi0",		&Kphi0_vec); 			CHKERRQ(ierr);
									//ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"rhol",		&rhol_vec); 			CHKERRQ(ierr);
									//ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Phi_old",		&Phi_vec); 				CHKERRQ(ierr);

	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"mul",			&Mul_vec); 				CHKERRQ(ierr);
									//ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"eta_creep",	&EtaCreep_vec); 		CHKERRQ(ierr);
	// New
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Kphi",		&Kphi_vec); 			CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Ssl",			&Ss_vec); 				CHKERRQ(ierr);
	ierr = 	DMRestoreNamedGlobalVector(jr->DA_Pl,"Pl_old",		&PlOld_vec); 			CHKERRQ(ierr);

	// cleaning up
									//ierr = 	VecDestroy(&local_KPhi0);		CHKERRQ(ierr);
									//ierr = 	VecDestroy(&local_rhol);		CHKERRQ(ierr);
									//ierr = 	VecDestroy(&local_Phi);			CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Mul);			CHKERRQ(ierr);
									//ierr = 	VecDestroy(&local_EtaCreep);	CHKERRQ(ierr);
	// New
	ierr = 	VecDestroy(&local_KPhi);		CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Ss);			CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_PlOld);		CHKERRQ(ierr);

//	ierr = 	VecDestroy(&Kphi0_vec);			CHKERRQ(ierr);
//	ierr = 	VecDestroy(&rhol_vec);			CHKERRQ(ierr);
//	ierr = 	VecDestroy(&Phi_vec);			CHKERRQ(ierr);
//	ierr = 	VecDestroy(&Mul_vec);			CHKERRQ(ierr);
//	ierr = 	VecDestroy(&EtaCreep_vec);		CHKERRQ(ierr);
		//---------------------------------------------

	ierr = DMCoarsenHookAdd(jr->DA_Pl,DMCoarsenHook_DARCY,PETSC_NULL,PETSC_NULL);

	//---------------------------------------------
									/*
									// (3) Initialize solution vector for Porosity of last timestep

									// Copy global solution vector into local  vector
									ierr = 	DMCreateLocalVector(jr->DA_Pl, &local_sol); 		 		 	  CHKERRQ(ierr);
									ierr = 	DMGlobalToLocalBegin(jr->DA_Pl, jr->Pl,INSERT_VALUES, local_sol); CHKERRQ(ierr);
									ierr = 	DMGlobalToLocalEnd  (jr->DA_Pl, jr->Pl,INSERT_VALUES, local_sol); CHKERRQ(ierr);

									// access work vectors
									ierr = 	DMDAVecGetArray(jr->DA_Pl, local_sol,  	&sol);  				  CHKERRQ(ierr);

									iter = 0;
									START_STD_LOOP
									{
										Phi 				= jr->svCell[iter++].svBulk.Phi;
										sol[k][j][i].Phi  	= Phi;
									}
									END_STD_LOOP


									ierr = DMDAVecRestoreArray(jr->DA_Pl, local_sol,  	&sol);  				CHKERRQ(ierr);
									ierr = DMLocalToGlobalBegin(jr->DA_Pl, local_sol,INSERT_VALUES, jr->Pl); 	CHKERRQ(ierr);
									ierr = DMLocalToGlobalEnd  (jr->DA_Pl, local_sol,INSERT_VALUES, jr->Pl); 	CHKERRQ(ierr);

									ierr = 	VecDestroy(&local_sol);												CHKERRQ(ierr);
									*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetDarcy_b"
PetscErrorCode GetDarcy_b(JacRes *jr)
{
	// This routine
	// Create b vector for liquid pressure of last time step

	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscScalar    dt    = jr->ts.dt;       // time step

	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	Vec 			local_sol, local_b, local_Ss;
	Vec 			Ss_vec	;

	DarcyDOF		***sol, ***lb;
	PetscScalar    	***Ssl;


	// Create temporary local vectors
	ierr = 	DMCreateLocalVector(jr->DA_Pl,&local_b); 				CHKERRQ(ierr);
	// Copy global solution vector into local vector
	ierr = 	DMGlobalToLocalBegin(jr->DA_Pl, jr->Darcyb,INSERT_VALUES, local_b); CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalEnd  (jr->DA_Pl, jr->Darcyb,INSERT_VALUES, local_b); CHKERRQ(ierr);

	ierr = 	DMCreateLocalVector (jr->DA_Pl, &local_sol); 		 		 	  CHKERRQ(ierr);
	// Copy global solution vector into local vector
	ierr = 	DMGlobalToLocalBegin(jr->DA_Pl, jr->Pl,INSERT_VALUES, local_sol); CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalEnd  (jr->DA_Pl, jr->Pl,INSERT_VALUES, local_sol); CHKERRQ(ierr);

	ExtractCoefficientsFromDA(jr->DA_Pl, "Ssl"		,	local_Ss, 		Ss_vec, 		Ssl);


	// access work vectors
	ierr = 	DMDAVecGetArray(jr->DA_Pl, local_sol,  &sol);  				CHKERRQ(ierr);
	ierr = 	DMDAVecGetArray(jr->DA_Pl, local_b,  	&lb);  				CHKERRQ(ierr);

	ierr = DMDAGetCorners(jr->DA_Pl, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
	iter = 0;
	START_STD_LOOP
	{
		lb[k][j][i].Pl = Ssl[k][j][i] * sol[k][j][i].Pl /dt;
	}
	END_STD_LOOP

	RestoreCoefficientsFromDA(jr->DA_Pl, "Ssl"		,	local_Ss, 		Ss_vec, 		Ssl);


	ierr = DMDAVecRestoreArray (jr->DA_Pl, local_sol,  	&sol);  			CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(jr->DA_Pl, local_sol,INSERT_VALUES, jr->Pl); 	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (jr->DA_Pl, local_sol,INSERT_VALUES, jr->Pl); 	CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray (jr->DA_Pl, local_b,  	&lb);  			CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(jr->DA_Pl, local_b,INSERT_VALUES, jr->Darcyb); 	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (jr->DA_Pl, local_b,INSERT_VALUES, jr->Darcyb); 	CHKERRQ(ierr);

	ierr = 	VecDestroy(&local_sol);									CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_b);									CHKERRQ(ierr);
	ierr = 	VecDestroy(&local_Ss);									CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DMCoarsenHook_DARCY"
PetscErrorCode DMCoarsenHook_DARCY(DM dmf,DM dmc,void *ctx)
{
	// This interpolates material properties from a finer DM to a coarser DM
	//	We use the interpolation matrix attached with the DM
	//	And the material properties must be specified on named vectors in the DA
							//PetscInt       	rlevel,clevel;
	PetscErrorCode 	ierr;
	Vec 			Rscale;
							//Vec 			Kphi0_fine, 	Kphi0_coarse;
							//Vec 			Phi_old_fine, 	Phi_old_coarse;
							//Vec 			rhol_fine, 		rhol_coarse;
	Vec 			mul_fine, 		mul_coarse;
							//Vec 			EtaCreep_fine, 	EtaCreep_coarse;
	Mat 			A;

	// New
	Vec 			Ssl_fine, 		Ssl_coarse;
	Vec 			Kphi_fine, 		Kphi_coarse;
	Vec 			Pl_old_fine, 	Pl_old_coarse;

	PetscFunctionBegin;

	// Get interpolation matrix from fine->coarse (or vice-versa if transpose is used)
	ierr = 	DMCreateInterpolation(dmc,dmf,&A,&Rscale);					CHKERRQ(ierr);

	// Get application data vectors
							//ierr = 	DMGetNamedGlobalVector(dmf,"Kphi0",&Kphi0_fine); 				CHKERRQ(ierr);	// fine level
							//ierr = 	DMGetNamedGlobalVector(dmc,"Kphi0",&Kphi0_coarse); 				CHKERRQ(ierr);	// coarse level

							//ierr = 	DMGetNamedGlobalVector(dmf,"Phi_old",&Phi_old_fine); 			CHKERRQ(ierr);	// fine level
							//ierr = 	DMGetNamedGlobalVector(dmc,"Phi_old",&Phi_old_coarse); 			CHKERRQ(ierr);	// coarse level

							//ierr = 	DMGetNamedGlobalVector(dmf,"rhol",&rhol_fine); 					CHKERRQ(ierr);	// fine level
							//ierr = 	DMGetNamedGlobalVector(dmc,"rhol",&rhol_coarse); 				CHKERRQ(ierr);	// coarse level

	// New
	ierr = 	DMGetNamedGlobalVector(dmf,"Kphi",&Kphi_fine); 				CHKERRQ(ierr);	// fine level
	ierr = 	DMGetNamedGlobalVector(dmc,"Kphi",&Kphi_coarse); 			CHKERRQ(ierr);	// coarse level


	ierr = 	DMGetNamedGlobalVector(dmf,"mul",&mul_fine); 				CHKERRQ(ierr);	// fine level
	ierr = 	DMGetNamedGlobalVector(dmc,"mul",&mul_coarse); 				CHKERRQ(ierr);	// coarse level

	// New
	ierr = 	DMGetNamedGlobalVector(dmf,"Ssl",&Ssl_fine); 				CHKERRQ(ierr);	// fine level
	ierr = 	DMGetNamedGlobalVector(dmc,"Ssl",&Ssl_coarse); 				CHKERRQ(ierr);	// coarse level

	ierr = 	DMGetNamedGlobalVector(dmf,"Pl_old",&Pl_old_fine); 			CHKERRQ(ierr);	// fine level
	ierr = 	DMGetNamedGlobalVector(dmc,"Pl_old",&Pl_old_coarse); 		CHKERRQ(ierr);	// coarse level

							//ierr = 	DMGetNamedGlobalVector(dmf,"eta_creep",&EtaCreep_fine); 		CHKERRQ(ierr);	// fine level
							//ierr = 	DMGetNamedGlobalVector(dmc,"eta_creep",&EtaCreep_coarse); 		CHKERRQ(ierr);	// coarse level

	// Perform the interpolation:
							//ierr = MatRestrict(A,Kphi0_fine,Kphi0_coarse);							CHKERRQ(ierr);
							//ierr = VecPointwiseMult(Kphi0_coarse,Kphi0_coarse, Rscale);				CHKERRQ(ierr);

							//ierr = MatRestrict(A,Phi_old_fine,Phi_old_coarse);						CHKERRQ(ierr);
							//ierr = VecPointwiseMult(Phi_old_coarse,Phi_old_coarse, Rscale);			CHKERRQ(ierr);

							//ierr = MatRestrict(A,rhol_fine,rhol_coarse);							CHKERRQ(ierr);
							//ierr = VecPointwiseMult(rhol_coarse,rhol_coarse, Rscale);				CHKERRQ(ierr);

	ierr = MatRestrict(A,mul_fine,mul_coarse);								CHKERRQ(ierr);
	ierr = VecPointwiseMult(mul_coarse,mul_coarse, Rscale);					CHKERRQ(ierr);

	// New
	ierr = MatRestrict(A,Ssl_fine,Ssl_coarse);								CHKERRQ(ierr);
	ierr = VecPointwiseMult(Ssl_coarse,Ssl_coarse, Rscale);					CHKERRQ(ierr);

							//ierr = MatRestrict(A,EtaCreep_fine,EtaCreep_coarse);					CHKERRQ(ierr);
							//ierr = VecPointwiseMult(EtaCreep_coarse,EtaCreep_coarse, Rscale);		CHKERRQ(ierr);

	// Restore vectors
							//ierr = 	DMRestoreNamedGlobalVector(dmc,"Kphi0",&Kphi0_coarse); 			CHKERRQ(ierr);	// coarse level
							//ierr = 	DMRestoreNamedGlobalVector(dmf,"Kphi0",&Kphi0_fine); 			CHKERRQ(ierr);	// fine level

							//ierr = 	DMRestoreNamedGlobalVector(dmf,"Phi_old",&Phi_old_fine); 		CHKERRQ(ierr);	// fine level
							//ierr = 	DMRestoreNamedGlobalVector(dmc,"Phi_old",&Phi_old_coarse); 		CHKERRQ(ierr);	// coarse level

							//ierr = 	DMRestoreNamedGlobalVector(dmf,"rhol",&rhol_fine); 				CHKERRQ(ierr);	// fine level
							//ierr = 	DMRestoreNamedGlobalVector(dmc,"rhol",&rhol_coarse); 			CHKERRQ(ierr);	// coarse level

	ierr = 	DMRestoreNamedGlobalVector(dmf,"mul",&mul_fine); 				CHKERRQ(ierr);	// fine level
	ierr = 	DMRestoreNamedGlobalVector(dmc,"mul",&mul_coarse); 				CHKERRQ(ierr);	// coarse level

	// New
	ierr = 	DMRestoreNamedGlobalVector(dmf,"Ssl",&Ssl_fine); 				CHKERRQ(ierr);	// fine level
	ierr = 	DMRestoreNamedGlobalVector(dmc,"Ssl",&Ssl_coarse); 				CHKERRQ(ierr);	// coarse level

	ierr = 	DMRestoreNamedGlobalVector(dmf,"Pl_old",&Pl_old_fine); 			CHKERRQ(ierr);	// fine level
	ierr = 	DMRestoreNamedGlobalVector(dmc,"Pl_old",&Pl_old_coarse); 		CHKERRQ(ierr);	// coarse level

	ierr = 	DMRestoreNamedGlobalVector(dmc,"Kphi",&Kphi_coarse); 			CHKERRQ(ierr);	// coarse level
	ierr = 	DMRestoreNamedGlobalVector(dmf,"Kphi",&Kphi_fine); 				CHKERRQ(ierr);	// fine level

						//ierr = 	DMRestoreNamedGlobalVector(dmf,"eta_creep",&EtaCreep_fine); 	CHKERRQ(ierr);	// fine level
						//ierr = 	DMRestoreNamedGlobalVector(dmc,"eta_creep",&EtaCreep_coarse); 	CHKERRQ(ierr);	// coarse level

	// Cleanup
	ierr = MatDestroy(&A);		CHKERRQ(ierr);
	ierr = VecDestroy(&Rscale);	CHKERRQ(ierr);
/*
	ierr =	VecDestroy(&Kphi0_coarse);		CHKERRQ(ierr);
	ierr =	VecDestroy(&Kphi0_fine);		CHKERRQ(ierr);
	ierr =	VecDestroy(&Phi_old_fine);		CHKERRQ(ierr);
	ierr =	VecDestroy(&Phi_old_coarse);	CHKERRQ(ierr);
	ierr =	VecDestroy(&rhol_fine);			CHKERRQ(ierr);
	ierr =	VecDestroy(&rhol_coarse);		CHKERRQ(ierr);
	ierr =	VecDestroy(&mul_fine);			CHKERRQ(ierr);
	ierr =	VecDestroy(&mul_coarse);		CHKERRQ(ierr);
	ierr =	VecDestroy(&EtaCreep_fine);		CHKERRQ(ierr);
	ierr =	VecDestroy(&EtaCreep_coarse);	CHKERRQ(ierr);
*/
	DMCoarsenHookAdd(dmc,DMCoarsenHook_DARCY,NULL,NULL);



	PetscFunctionReturn(0);
}


/*
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ExtractCoefficientsFromDA"
PetscErrorCode ExtractCoefficientsFromDA(DM da,const char *name,PetscScalar ***data)
{
	// This extracts material coefficients from named vectors that
	//	are attached to the DA, and returns a 3D array

	Vec 			local, global;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	ierr = 	DMGetLocalVector(da,&local); 								CHKERRQ(ierr);
	ierr = 	DMGetNamedGlobalVector(da,"Kphi0",&global); 				CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalBegin(da, global,INSERT_VALUES, local); 		CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalEnd  (da, global,INSERT_VALUES, local); 		CHKERRQ(ierr);
	ierr = 	DMDAVecGetArray(da, local, &data);  						CHKERRQ(ierr);		// 3D array with data that was extracted



	PetscFunctionReturn(0);
}
*/
