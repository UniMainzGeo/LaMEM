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
 **    filename:   JacRes.h
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
//... FDSTAG DARCY/2PhaseFlow JACOBIAN AND RESIDUAL ROUTINE FUNCTIONS  ......
//---------------------------------------------------------------------------
#ifndef __JacResDarcy_h__
#define __JacResDarcy_h__

//---------------------------------------------------------------------------
// * replace setting time parameters consistently in the entire code

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//......................   DARCY/LiquidPressure FUNCTIONS   ..................
//---------------------------------------------------------------------------

PetscErrorCode JacResGetDarcyParam(
		JacRes      *jr,
		PetscScalar *phRat,
		PetscScalar *Kphi_,     // Permeability
						//PetscScalar *rhol_, 	// Liquid density
		PetscScalar *mu_, 	// Liquid viscosity
		PetscScalar *Ss_);	// New: Specific storage

// check whether Darcy material parameters are properly defined
PetscErrorCode JacResCheckDarcyParam(JacRes *jr);

// setup Darcy/LiquidPressure parameters
PetscErrorCode JacResCreateDarcyParam(JacRes *jr);

// destroy Darcy/LiquidPressure parameters
PetscErrorCode JacResDestroyDarcyParam(JacRes *jr);

// initialize Darcy/LiquidPressure from markers
PetscErrorCode JacResInitDarcy(JacRes *jr);

// correct Darcy for diffusion (Newton update)
PetscErrorCode JacResUpdateDarcy(JacRes *jr);

// apply Darcy two-point constraints
PetscErrorCode JacResApplyDarcyBC(JacRes *jr);

// compute Darcy residual vector
PetscErrorCode JacResGetDarcyRes(SNES snes, Vec x, Vec f, JacRes *jr);

// assemble Darcy preconditioner matrix
//PetscErrorCode JacResGetDarcyMat(JacRes *jr);
PetscErrorCode JacResGetDarcyMat(SNES snes, JacRes *jr);

// BC for Darcy:
PetscErrorCode BCCreateDarcy(JacRes *jr, BCCtx *bc);


PetscErrorCode BCApplyBound_DARCY(BCCtx *bc, JacRes *jr);

PetscErrorCode FormFunction_DARCY(SNES snes,Vec x,Vec f, JacRes *jr);
PetscErrorCode FormJacobian_DARCY(SNES snes,Vec x, Mat P, Mat J, JacRes *jr);

PetscErrorCode UpdateDarcy_DA(JacRes *jr);
PetscErrorCode DMCoarsenHook_DARCY(DM dmf,DM dmc,void *ctx);

PetscErrorCode FormRHS_DARCY(SNES snes,JacRes *jr);

// New
PetscErrorCode JacResCopyDarcySol(SNES snes, JacRes *jr, Vec x);
PetscErrorCode JacResUpdateGhostPoints(SNES snes, Vec x, JacRes *jr);
PetscErrorCode BCDestroyDarcy(JacRes *jr, BCCtx *bc);
PetscErrorCode BCSetParamDarcy(JacRes *jr, BCCtx *bc, UserCtx *user);
PetscErrorCode IncreaseLiquidPressureBottom(JacRes *jr, BCCtx *bc, UserCtx *user);
PetscErrorCode DarcyPrintPl(JacRes *jr);

//PetscErrorCode ExtractCoefficientsFromDA(DM da,const char *name,PetscScalar ***data);	// extract names vectors from DA


//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

/*#define SET_TPC_DARCY(bc, a, k, j, i, pmdof) { \
	if(bc[k][j][i] == DBL_MAX) a[k][j][i].Pl = pmdof; \
	else                       a[k][j][i].Pl = 2.0*bc[k][j][i] - pmdof; }

#define SET_EDGE_CORNER_DARCY(n, a, K, J, I, k, j, i, pmdof) \
	a[K][J][I].Pl = a[k][j][I].Pl + a[k][J][i].Pl + a[K][j][i].Pl - 2.0*pmdof;*/

#define SET_TPC_DARCY(bc, a, k, j, i, pmdof) { \
	if(bc[k][j][i] == DBL_MAX) a[k][j][i] = pmdof; \
	else                       a[k][j][i] = 2.0*bc[k][j][i] - pmdof; }

#define SET_EDGE_CORNER_DARCY(n, a, K, J, I, k, j, i, pmdof) \
	a[K][J][I] = a[k][j][I] + a[k][J][i] + a[K][j][i] - 2.0*pmdof;


//---------------------------------------------------------------------------
#endif
