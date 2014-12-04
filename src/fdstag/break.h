//---------------------------------------------------------------------------
//........................  BREAKPOINT ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __break_h__
#define __break_h__

//---------------------------------------------------------------------------
// Routines for writing breakpoint files
//---------------------------------------------------------------------------
PetscErrorCode BreakWriteMain     (UserContext *user, AdvCtx *actx);
PetscErrorCode BreakWriteMark     (AdvCtx *actx);
PetscErrorCode BreakWriteSol      (JacRes *jr);
PetscErrorCode BreakWriteGrid     (UserContext *user, FDSTAG *fs, AdvCtx *actx);
PetscErrorCode BreakWriteInfo     (UserContext *user, AdvCtx *actx, JacRes *jr);
PetscErrorCode BreakWriteDiscret1D(int fid, Discret1D ds);

//---------------------------------------------------------------------------
// Routines for reading breakpoint files
//---------------------------------------------------------------------------
PetscErrorCode BreakReadMain     (UserContext *user, AdvCtx *actx, JacRes *jr);
PetscErrorCode BreakReadMark     (AdvCtx *actx);
PetscErrorCode BreakReadSol      (JacRes *jr);
PetscErrorCode BreakReadGrid     (UserContext *user, FDSTAG *fs);
PetscErrorCode BreakReadInfo     (UserContext *user, AdvCtx *actx, JacRes *jr);
PetscErrorCode BreakReadDiscret1D(int fid, Discret1D ds);

#endif
