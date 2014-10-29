//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#ifndef __input_h__
#define __input_h__
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInitCode(UserContext *user);

PetscErrorCode FDSTAGReadInputFile(UserContext *user);

PetscErrorCode ReadMaterialProperties(UserContext *user);

//---------------------------------------------------------------------------

#endif
