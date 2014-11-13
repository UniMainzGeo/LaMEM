//---------------------------------------------------------------------------
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#ifndef __bc_h__
#define __bc_h__
//---------------------------------------------------------------------------
// boundary condition context
typedef struct
{
	// boundary conditions vectors (velocity, pressure, temperature)
	Vec bcvx,  bcvy, bcvz, bcp, bcT; // local (ghosted)

	// single-point constraints
	PetscInt     numSPC;   // number of single point constraints (SPC)
	PetscInt    *SPCList;  // global indices of SPC (global layout)
	PetscScalar *SPCVals;  // values of SPC

	PetscInt     numSPCPres;   // number of pressure SPC
	PetscInt    *SPCListPres;  // global indices of pressure SPC (pressure layout)

	// two-point constraints
	PetscInt     numTPC;       // number of two-point constraints (TPC)
	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
	PetscScalar *TPCVals;      // values of TPC
	PetscScalar *TPCLinComPar; // linear combination parameters

	// background strain-rate parameters
	PetscBool    bgActive;     // flag for activating background strain-rates
	PetscScalar  Exx, Eyy;     // horizontal background strain-rates

	// Dirichlet pushing block constraints
	PetscBool     pActive;     // flag for activating pushing
	PetscBool     pApply;      // flag for applying pushing on a time step
	PetscScalar   theta;       // rotation angle
	PetscScalar   Vx, Vy;      // Dirichlet values for Vx and Vy
	PushingParams *pb;         // major pushing block parameters

} BCCtx;
//---------------------------------------------------------------------------

PetscErrorCode BCClear(BCCtx *bc);

// create boundary condition context
PetscErrorCode BCCreate(BCCtx *bc, FDSTAG *fs, idxtype idxmod);

// set background strain-rates
PetscErrorCode BCSetStretch(BCCtx *bc, UserContext *user);

// destroy boundary condition context
PetscErrorCode BCDestroy(BCCtx *bc);

// apply boundary conditions
PetscErrorCode BCApply(BCCtx *bc, FDSTAG *fs, TSSol *ts, Scaling *scal, idxtype idxmod);

// apply constraints on the boundaries
PetscErrorCode BCApplyBound(BCCtx *bc, FDSTAG *fs);

//---------------------------------------------------------------------------

// initialize pushing boundary conditions context
PetscErrorCode BCSetPush(BCCtx *bc, UserContext *user);

// compute pushing parameters
PetscErrorCode BCCompPush(BCCtx *bc, TSSol *ts, Scaling *scal);

// apply pushing constraints
PetscErrorCode BCApplyPush(BCCtx *bc, FDSTAG *fs);

// advect the pushing block
PetscErrorCode BCAdvectPush(BCCtx *bc, TSSol *ts);

//---------------------------------------------------------------------------

#endif
