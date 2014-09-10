//---------------------------------------------------------------------------
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx, FDSTAG *fs, UserContext *user);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx, FDSTAG *fs, UserContext *user);

// add random noise to marker coordinates if required
PetscErrorCode ADVMarkRandomNoise(AdvCtx *actx, FDSTAG *fs, UserContext *user);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx, FDSTAG *fs, UserContext *user);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckPhaseIDs(AdvCtx *actx, UserContext *user);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx, FDSTAG *fs, UserContext *user);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx, FDSTAG *fs, UserContext *user);

//---------------------------------------------------------------------------
#endif
