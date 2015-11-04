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
 **    filename:   bc.h
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
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#ifndef __bc_h__
#define __bc_h__
//---------------------------------------------------------------------------

#define _max_periods_ 20
#define _max_path_points_ 30
#define _max_poly_points_ 15
#define _max_bc_blocks_ 3

//---------------------------------------------------------------------------
// index shift type
typedef enum
{
	_LOCAL_TO_GLOBAL_,
	_GLOBAL_TO_LOCAL_

} ShiftType;
//---------------------------------------------------------------------------
typedef struct
{
	// path description
	PetscInt    npath;                        // number of path points of Bezier curve
	PetscScalar theta[  _max_path_points_  ]; // orientation angles at path points
	PetscScalar time [  _max_path_points_  ]; // times at path points
	PetscScalar path [6*_max_path_points_-4]; // Bezier curve path & control points

	// block description
	PetscInt    npoly;                      // number of polygon vertices
	PetscScalar poly [2*_max_poly_points_]; // polygon coordinates
	PetscScalar bot, top;                   // bottom & top coordinates of the block

} BCBlock;
//---------------------------------------------------------------------------

PetscErrorCode BCBlockReadFromOptions(BCBlock *bcb, Scaling *scal);

PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *act, PetscScalar x[]);

PetscErrorCode BCBlockGetPolygon(BCBlock *bcb, PetscScalar Xb[], PetscScalar *cpoly);

//---------------------------------------------------------------------------
// boundary condition context
typedef struct
{
	//=====================================================================
	//
	// Boundary condition vectors contain prescribed DOF values:
	//
	//    *Internal points (marked with positive number in the index arrays)
	//        DBL_MAX   - active DOF flag
	//        otherwise - single-point constraint (SPC) value
	//
	//    *Boundary ghost point (marked with -1 in the index arrays)
	//        DBL_MAX   - free-slip (zero-flux) condition flag
	//        otherwise - two-point constraint (TPC) value
	//
	// Boundary ghost points require consistent setting of constraints
	// on the processor boundaries (since PETSc doesn't exchange boundary
	// ghost point values). Internal ghost points should be synchronized
	// after initializing the single-point constraints. Synchronization
	// can be skipped if all ghost points are initialized redundantly
	// on all the processes (DO THIS!).
	//
	// Single point constraints are additionally stored as lists
	// for constraining matrices and vectors. Matrices require global
	// index space, vectors require local index space.
	//
	// Global v-p index space can be either monolithic or block
	// Local v-p index space is ALWAYS coupled (since all solvers are coupled)
	//
	// NOTE! It may be worth storing TPC also as lists (for speedup).
	//=====================================================================

	FDSTAG   *fs;   // staggered grid
	TSSol    *ts;   // time stepping parameters
	Scaling  *scal; // scaling parameters

	// boundary conditions vectors (velocity, pressure, temperature)
	Vec bcvx, bcvy, bcvz, bcp, bcT; // local (ghosted)

	// single-point constraints
	ShiftType    stype;   // current index shift type
	PetscInt     numSPC;  // total number of constraints
	PetscInt    *SPCList; // local indices of SPC
	PetscScalar *SPCVals; // values of SPC

	// velocity
	PetscInt     vNumSPC;
	PetscInt    *vSPCList;
	PetscScalar *vSPCVals;

	// pressure
	PetscInt     pNumSPC;
	PetscInt    *pSPCList;
	PetscScalar *pSPCVals;

	// temperature
	PetscInt     tNumSPC;
	PetscInt    *tSPCList;
	PetscScalar *tSPCVals;

	// temperature on top and bottom boundaries
	PetscScalar  Tbot, Ttop;

	// horizontal background strain-rate parameters
	PetscBool    ExxAct;
	PetscInt     ExxNumPeriods;
	PetscScalar  ExxTimeDelims [_max_periods_-1];
	PetscScalar  ExxStrainRates[_max_periods_  ];

	PetscBool    EyyAct;
	PetscInt     EyyNumPeriods;
	PetscScalar  EyyTimeDelims [_max_periods_-1];
	PetscScalar  EyyStrainRates[_max_periods_  ];

	// Dirichlet pushing block constraints
	PetscBool     pbAct;  // flag for activating pushing
	PetscBool     pbApp;  // flag for applying pushing on a time step
	PetscScalar   theta;  // rotation angle
	PetscScalar   Vx, Vy; // Dirichlet values for Vx and Vy
	PushParams    *pb;    // major pushing block parameters

	// two-point constraints
//	PetscInt     numTPC;       // number of two-point constraints (TPC)
//	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
//	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
//	PetscScalar *TPCVals;      // values of TPC
//	PetscScalar *TPCLinComPar; // linear combination parameters

	BCBlock      blocks;  // BC block
//	PetscInt     bcbAct;  // BC block activation flag


} BCCtx;
//---------------------------------------------------------------------------

PetscErrorCode BCClear(BCCtx *bc);

// create boundary condition context
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs, TSSol *ts, Scaling *scal);

// set background strain-rates
PetscErrorCode BCSetParam(BCCtx *bc, UserCtx *user);

// set parameters from PETSc options
PetscErrorCode BCReadFromOptions(BCCtx *bc);

// get current background strain rates
PetscErrorCode BCGetBGStrainRates(
	BCCtx       *bc,
	PetscScalar *Exx_,
	PetscScalar *Eyy_,
	PetscScalar *Ezz_);

// destroy boundary condition context
PetscErrorCode BCDestroy(BCCtx *bc);

// apply boundary conditions
PetscErrorCode BCApply(BCCtx *bc);

// apply constraints on the boundaries
PetscErrorCode BCApplyBound(BCCtx *bc);

// shift indices of constrained nodes
PetscErrorCode BCShiftIndices(BCCtx *bc, ShiftType stype);

//---------------------------------------------------------------------------

// initialize pushing boundary conditions context
PetscErrorCode BCSetPush(BCCtx *bc, UserCtx *user);

// compute pushing parameters
PetscErrorCode BCCompPush(BCCtx *bc);

// apply pushing constraints
PetscErrorCode BCApplyPush(BCCtx *bc);

// advect the pushing block
PetscErrorCode BCAdvectPush(BCCtx *bc);

// stretch staggered grid if background strain rates are defined
PetscErrorCode BCStretchGrid(BCCtx *bc);

//---------------------------------------------------------------------------

PetscErrorCode BCApplyBezier(BCCtx *bc);

//---------------------------------------------------------------------------

#endif
