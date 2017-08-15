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
		PetscScalar *rhol_,      // Liquid density
		PetscScalar *mul_,       // Liquid viscosity
		PetscScalar *Kphi_,      // Permeability
		//PetscScalar *Ss_,      // Specific storage
		PetscScalar *betam_,     // Matrix compressibility
		PetscScalar *betal_,     // Liquid compressibility
		PetscScalar *Phi_,       // Effective porosity
		PetscScalar *Ts_,        // Tensile strength
		PetscScalar *nuu_,       // Undrained poisson's ratio
		PetscScalar *Kphiu_,     // Undrained permeability
		PetscScalar *Phiu_);      // Undrained porosity

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

PetscErrorCode JacResGetDarcyRes(JacRes *jr);
PetscErrorCode JacResUpdateDarcyPermeability(JacRes *jr);

// assemble Darcy preconditioner matrix
PetscErrorCode JacResGetDarcyMat(JacRes *jr);

PetscErrorCode UpdateDarcy_DA(JacRes *jr);
PetscErrorCode SolveDarcyKSP(JacRes *jr);

PetscErrorCode GetCellCoordinatesDarcySources(JacRes *jr);
PetscErrorCode DarcySourcePropInit(JacRes *jr, FILE *fp);
PetscErrorCode SourcePropGetStruct(FILE *fp,PetscInt numSources, DarcySourceParam *sources,PetscInt ils, PetscInt ile, UnitsType utype);


//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define SET_TPC_DARCY(bc, a, k, j, i, pmdof) { \
	if(bc[k][j][i] == DBL_MAX) a[k][j][i] = pmdof; \
	else                       a[k][j][i] = 2.0*bc[k][j][i] - pmdof; }

#define SET_EDGE_CORNER_DARCY(n, a, K, J, I, k, j, i, pmdof) \
	a[K][J][I] = a[k][j][I] + a[k][J][i] + a[K][j][i] - 2.0*pmdof;


//---------------------------------------------------------------------------
#endif
