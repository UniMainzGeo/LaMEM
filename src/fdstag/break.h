//---------------------------------------------------------------------------
//........................  BREAKPOINT ROUTINES   ...........................
//---------------------------------------------------------------------------
#ifndef __break_h__
#define __break_h__

// check if breakpoints exist and restart new simulation if not
PetscErrorCode BreakCheck    (UserCtx *user);

//---------------------------------------------------------------------------
// Write breakpoint files
//---------------------------------------------------------------------------
// grid, mark, gsol, info are written together
PetscErrorCode BreakWrite    (UserCtx *user, AdvCtx *actx, JacType jtype);

//---------------------------------------------------------------------------
// Read breakpoint files
//---------------------------------------------------------------------------
// gsol, info are read together; grid, mark are read separately
PetscErrorCode BreakRead     (UserCtx *user, AdvCtx *actx, JacType *jtype);
PetscErrorCode BreakReadGrid (UserCtx *user, FDSTAG *fs);
PetscErrorCode BreakReadMark (AdvCtx *actx);

//---------------------------------------------------------------------------
// Read and write FDSTAG 1D structures
//---------------------------------------------------------------------------
void BreakWriteDiscret1D (FILE *fp, Discret1D ds, MeshSeg1D ms);
void BreakReadDiscret1D  (FILE *fp, Discret1D ds, MeshSeg1D ms);

#endif
