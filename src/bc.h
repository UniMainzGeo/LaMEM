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
#define _max_boxes_ 5
#define _max_path_points_ 25
#define _max_poly_points_ 50

//---------------------------------------------------------------------------

struct FB;
struct Scaling;
struct TSSol;
struct FDSTAG;
struct Marker;
struct DBMat;
struct JacRes;

//---------------------------------------------------------------------------
// index shift type
enum ShiftType
{
	_LOCAL_TO_GLOBAL_,
	_GLOBAL_TO_LOCAL_

};

//---------------------------------------------------------------------------
// Bezier block (rotating polygon moving along Bezier curve)
//---------------------------------------------------------------------------

struct BCBlock
{
	// path description
	PetscInt    npath;                        // number of path points of Bezier curve
	PetscScalar theta[  _max_path_points_  ]; // orientation angles at path points
	PetscScalar time [  _max_path_points_  ]; // times at path points
	PetscScalar path [6*_max_path_points_-4]; // Bezier curve path & control points (3*n-2)

	// block description
	PetscInt    npoly;                      // number of polygon vertices
	PetscScalar poly [2*_max_poly_points_]; // polygon coordinates
	PetscScalar bot, top;                   // bottom & top coordinates of the block

	// WARNING bottom coordinate should be advected (how? average?)
	// Top of the box can be assumed to be the free surface
	// sticky air nodes should never be constrained (this is easy to check)

};

//---------------------------------------------------------------------------

// setup data structures
PetscErrorCode BCBlockCreate(BCBlock *bcb, Scaling *scal, FB *fb);

// compute position along the path and rotation angle as a function of time
PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *f, PetscScalar x[]);

// compute current polygon coordinates
PetscErrorCode BCBlockGetPolygon(BCBlock *bcb, PetscScalar Xb[], PetscScalar *cpoly);

//---------------------------------------------------------------------------
// Dropping boxes (rectangular boxes moving with constant vertical velocity)
//---------------------------------------------------------------------------

struct DBox
{
	PetscInt    num;                   // number of boxes
	PetscScalar bounds[6*_max_boxes_]; // box bounds
	PetscScalar zvel;                  // vertical velocity

} ;

//---------------------------------------------------------------------------

PetscErrorCode DBoxReadCreate(DBox *dbox, Scaling *scal, FB *fb);

//---------------------------------------------------------------------------

// boundary condition context
struct BCCtx
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
	DBMat    *dbm;  // material database
	JacRes   *jr;   // Jacobian-residual context (CROSS-REFERENCE!)

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

	// two-point constraints
//	PetscInt     numTPC;       // number of two-point constraints (TPC)
//	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
//	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
//	PetscScalar *TPCVals;      // values of TPC
//	PetscScalar *TPCLinComPar; // linear combination parameters

	//=====================
	// VELOCITY CONSTRAINTS
	//=====================

	// horizontal background strain-rate parameters
	PetscInt     ExxNumPeriods;
	PetscScalar  ExxTimeDelims [_max_periods_-1];
	PetscScalar  ExxStrainRates[_max_periods_  ];

	PetscInt     EyyNumPeriods;
	PetscScalar  EyyTimeDelims [_max_periods_-1];
	PetscScalar  EyyStrainRates[_max_periods_  ];

	// Bezier block
	PetscInt 	 nblocks;             // number of Bezier blocks
	BCBlock      blocks[_max_boxes_]; // BC block

	// dropping boxes
	DBox         dbox;

	// velocity inflow & outflow boundary condition
	PetscInt     face, phase;   // face (1-left 2-right 3-front 4-back) & phase identifiers
	PetscScalar  bot, top;      // bottom & top coordinates of the plate
	PetscScalar  velin, velout; // inflow & outflow velocities

	// open boundary flag
	PetscInt     top_open;

	// no-slip boundary condition mask
	PetscInt     noslip[6];

	// fixed phase (no-flow condition)
	PetscInt     fixPhase;

	//========================
	// TEMPERATURE CONSTRAINTS
	//========================

	// temperature on top and bottom boundaries
	PetscScalar  Tbot, Ttop;

	//=====================
	// PRESSURE CONSTRAINTS
	//=====================

	// pressure on top and bottom boundaries
	PetscScalar  pbot, ptop;
};
//---------------------------------------------------------------------------

// create boundary condition context
PetscErrorCode BCCreate(BCCtx *bc, FB *fb);

// allocate internal vectors and arrays
PetscErrorCode BCCreateData(BCCtx *bc);

// destroy boundary condition context
PetscErrorCode BCDestroy(BCCtx *bc);

// apply ALL boundary conditions
PetscErrorCode BCApply(BCCtx *bc);

// apply SPC to global solution vector
PetscErrorCode BCApplySPC(BCCtx *bc);

// shift indices of constrained nodes
PetscErrorCode BCShiftIndices(BCCtx *bc, ShiftType stype);

//---------------------------------------------------------------------------
// Specific constraints
//---------------------------------------------------------------------------

// apply pressure constraints
PetscErrorCode BCApplyPres(BCCtx *bc);

// apply temperature constraints
PetscErrorCode BCApplyTemp(BCCtx *bc);

// apply default velocity constraints on the boundaries
PetscErrorCode BCApplyVelDefault(BCCtx *bc);

// apply Bezier blocks
PetscErrorCode BCApplyBezier(BCCtx *bc);

// apply inflow/outflow boundary velocities
PetscErrorCode BCApplyBoundVel(BCCtx *bc);

// apply dropping boxes
PetscErrorCode BCApplyDBox(BCCtx *bc);

// constraint all cells containing phase
PetscErrorCode BCApplyPhase(BCCtx *bc);

// create SPC constraint lists
PetscErrorCode BCListSPC(BCCtx *bc);

// apply two-point constraints on the boundaries
PetscErrorCode BCApplyVelTPC(BCCtx *bc);

//---------------------------------------------------------------------------
// Service functions
//---------------------------------------------------------------------------

// get current background strain rates
PetscErrorCode BCGetBGStrainRates(
	BCCtx       *bc,
	PetscScalar *Exx_,
	PetscScalar *Eyy_,
	PetscScalar *Ezz_);

// stretch staggered grid if background strain rates are defined
PetscErrorCode BCStretchGrid(BCCtx *bc);

// change phase of inflow markers if velocity boundary condition is defined
PetscErrorCode BCOverridePhase(BCCtx *bc, PetscInt cellID, Marker *P);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define LIST_SPC(bc, list, vals, cnt, iter)\
	if(bc[k][j][i] != DBL_MAX) { list[cnt] = iter; vals[cnt] = bc[k][j][i]; cnt++; }

//---------------------------------------------------------------------------
#endif
