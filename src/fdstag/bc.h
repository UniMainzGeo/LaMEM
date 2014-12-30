//---------------------------------------------------------------------------
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#ifndef __bc_h__
#define __bc_h__
//---------------------------------------------------------------------------
// index shift type
typedef enum
{
	LOCAL_TO_GLOBAL,
	GLOBAL_TO_LOCAL

} ShiftType;
//---------------------------------------------------------------------------
// boundary condition context
typedef struct
{
	//=====================================================================
	// WARNING!
	//
	// Global v-p index space can be either coupled or uncoupled
	// (used to constrain rows & columns of preconditioner matrices).
	//
	// Local v-p index space is ALWAYS coupled, since all solvers are coupled
	// (used to constrain primary unknown & residual vectors).
	//=====================================================================

	// boundary conditions vectors (velocity, pressure, temperature)
	Vec bcvx,  bcvy, bcvz, bcp, bcT; // local (ghosted)

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

	// two-point constraints
//	PetscInt     numTPC;       // number of two-point constraints (TPC)
//	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
//	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
//	PetscScalar *TPCVals;      // values of TPC
//	PetscScalar *TPCLinComPar; // linear combination parameters

	PetscScalar  Tbot, Ttop; // temperature on top and bottom boundaries

	// background strain-rate parameters
	PetscBool    bgAct;    // flag for activating background strain-rates
	PetscScalar  Exx, Eyy; // horizontal background strain-rates

	// Dirichlet pushing block constraints
	PetscBool     pbAct;  // flag for activating pushing
	PetscBool     pbApp;  // flag for applying pushing on a time step
	PetscScalar   theta;  // rotation angle
	PetscScalar   Vx, Vy; // Dirichlet values for Vx and Vy
	PushingParams *pb;    // major pushing block parameters
	TSSol         *ts;    // time stepping parameters
	Scaling       *scal;  // scaling parameters


} BCCtx;
//---------------------------------------------------------------------------

PetscErrorCode BCClear(BCCtx *bc);

// create boundary condition context
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs, TSSol *ts, Scaling *scal);

// set background strain-rates
PetscErrorCode BCSetStretch(BCCtx *bc, UserContext *user);

// destroy boundary condition context
PetscErrorCode BCDestroy(BCCtx *bc);

// apply boundary conditions
PetscErrorCode BCApply(BCCtx *bc, FDSTAG *fs);

// apply constraints on the boundaries
PetscErrorCode BCApplyBound(BCCtx *bc, FDSTAG *fs);

// shift indices of constrained nodes
PetscErrorCode BCShiftIndices(BCCtx *bc, FDSTAG *fs, ShiftType stype);

//---------------------------------------------------------------------------

// initialize pushing boundary conditions context
PetscErrorCode BCSetPush(BCCtx *bc, UserContext *user);

// compute pushing parameters
PetscErrorCode BCCompPush(BCCtx *bc);

// apply pushing constraints
PetscErrorCode BCApplyPush(BCCtx *bc, FDSTAG *fs);

// advect the pushing block
PetscErrorCode BCAdvectPush(BCCtx *bc, TSSol *ts);

//---------------------------------------------------------------------------

#endif
