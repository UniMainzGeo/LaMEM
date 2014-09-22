//---------------------------------------------------------------------------
//..............   LaMEM - FDSTAG CANONICAL INTERFACE ROUTINES   ............
//---------------------------------------------------------------------------
#ifndef __interface_h__
#define __interface_h__
//---------------------------------------------------------------------------

// initialize material properties in the FDSTAG data structures
PetscErrorCode FDSTAGInitMaterialProps(JacRes *jr, UserContext *usr);

//---------------------------------------------------------------------------
#endif
