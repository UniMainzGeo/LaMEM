//---------------------------------------------------------------------------
//..............   LaMEM - FDSTAG CANONICAL INTERFACE ROUTINES   ............
//---------------------------------------------------------------------------
#ifndef __interface_h__
#define __interface_h__
//---------------------------------------------------------------------------

// copy local coordinates from DMDA object to FDSTAG
//PetscErrorCode FDSTAGSetCoordDMDA(FDSTAG *fs, DM da);

// initialize material properties in the FDSTAG data structures
PetscErrorCode FDSTAGInitMaterialProps(JacResCtx *jrctx, UserContext *usr);

// initialize phase ratios in the FDSTAG data structures
// PetscErrorCode FDSTAGInitPhaseRatios(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr);

// project effective viscosity field from FDSTAG to LaMEM layout
PetscErrorCode FDSTAGProjectEffVisc(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr);

// copy residual from FDSTAG to LaMEM collocated layout
PetscErrorCode FDSTAGCopyResidual(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr, Vec f, Vec g);

// copy solution from LaMEM collocated layout to FDSTAG & set boundary values
PetscErrorCode FDSTAGCopySolution(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr, Vec v, Vec p);

//---------------------------------------------------------------------------
#endif
