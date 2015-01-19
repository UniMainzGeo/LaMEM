//---------------------------------------------------------------------------
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, UserCtx *user);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, UserCtx *user);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx, UserCtx *user);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx, UserCtx *user);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx, UserCtx *user);
PetscErrorCode ADVMarkInitBands        (AdvCtx *actx, UserCtx *user);

//---------------------------------------------------------------------------
#endif
