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
	// NOTE: get rid of these vectors, by extending single- and two-point constraint specification
	Vec bcvx,  bcvy, bcvz, bcp, bcT; // local (ghosted)

	// single-point constraints
	PetscInt     numSPC;   // number of single point constraints (SPC)
	PetscInt    *SPCList;  // global indices of SPC (global layout)

	PetscInt     numSPCPres;   // number of pressure SPC
	PetscInt    *SPCListPres;  // global indices of pressure SPC (pressure layout)

	// two-point constraints
	PetscInt     numTPC;       // number of two-point constraints (TPC)
	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
	PetscScalar *TPCVals;      // values of TPC
	PetscScalar *TPCLinComPar; // linear combination parameters

} BCCtx;
//---------------------------------------------------------------------------

// create boundary condition context
PetscErrorCode FDSTAGCreateBCCtx(BCCtx *bc, FDSTAG *fs);

// destroy boundary condition context
PetscErrorCode FDSTAGDestroyBCCtx(BCCtx *bc);

// initialize boundary constraint vectors
PetscErrorCode FDSTAGInitBC(BCCtx *bc, FDSTAG *fs, idxtype idxmod);

//---------------------------------------------------------------------------
#endif
