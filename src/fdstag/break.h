//---------------------------------------------------------------------------
//........................  BREAKPOINT ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __break_h__
#define __break_h__

//---------------------------------------------------------------------------
// Routines for writing breakpoint files
//---------------------------------------------------------------------------
PetscErrorCode BreakWriteMain     (UserCtx *user, AdvCtx *actx, JacType jtype);
PetscErrorCode BreakWriteInfo     (UserCtx *user, AdvCtx *actx, JacType jtype);
PetscErrorCode BreakWriteSol      (JacRes *jr);
PetscErrorCode BreakWriteMark     (AdvCtx *actx);
PetscErrorCode BreakWriteGrid     (UserCtx *user, AdvCtx *actx);
PetscErrorCode BreakWriteDiscret1D(int fid, Discret1D ds);
PetscErrorCode BreakWriteStrech1D (int fid, MeshSeg1D ms);

//---------------------------------------------------------------------------
// Routines for reading breakpoint files
//---------------------------------------------------------------------------
PetscErrorCode BreakReadMain     (UserCtx *user, AdvCtx *actx, JacType *jtype);
PetscErrorCode BreakReadInfo     (UserCtx *user, AdvCtx *actx, JacType *jtype);
PetscErrorCode BreakReadSol      (JacRes *jr);
PetscErrorCode BreakReadMark     (AdvCtx *actx);
PetscErrorCode BreakReadGrid     (UserCtx *user, FDSTAG *fs);
PetscErrorCode BreakReadDiscret1D(int fid, Discret1D ds);
PetscErrorCode BreakReadStrech1D (int fid, MeshSeg1D ms);

#endif
