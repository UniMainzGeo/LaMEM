//---------------------------------------------------------------------------
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#ifndef __bc_h__
#define __bc_h__
//---------------------------------------------------------------------------

#define _max_periods_ 20

//---------------------------------------------------------------------------
// index shift type
typedef enum
{
	_LOCAL_TO_GLOBAL_,
	_GLOBAL_TO_LOCAL_

} ShiftType;
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

#endif
