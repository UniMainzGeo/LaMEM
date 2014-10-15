//---------------------------------------------------------------------------
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, UserContext *user);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, UserContext *user);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx, UserContext *user);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx, UserContext *user);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx, UserContext *user);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx, UserContext *user);

//---------------------------------------------------------------------------
#endif
