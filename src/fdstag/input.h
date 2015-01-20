//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#ifndef __input_h__
#define __input_h__
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user);

PetscErrorCode FDSTAGSetDefaultValues(UserCtx *user);

PetscErrorCode FDSTAGReadInputFile(JacRes *jr, UserCtx *user);

PetscErrorCode FDSTAGReadCommLine(UserCtx *user );

PetscErrorCode ReadMaterialProperties(UserCtx *user);

PetscErrorCode FDSTAGInitMaterialProp(UserCtx *user );

//---------------------------------------------------------------------------
PetscErrorCode ReadMeshSegDir(
	FILE        *fp,
	const char  *name,
	PetscScalar  beg,
	PetscScalar  end,
	PetscInt    *tncels,
	MeshSegInp  *msi,
	PetscInt     dim,
	PetscScalar  charLength);

//---------------------------------------------------------------------------

#endif
