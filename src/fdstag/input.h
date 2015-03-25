//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#ifndef __input_h__
#define __input_h__
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user);

PetscErrorCode InputSetDefaultValues(JacRes *jr, UserCtx *user);

PetscErrorCode InputReadFile(JacRes *jr, UserCtx *user, FILE *fp);

PetscErrorCode InputReadCommLine(UserCtx *user);

//---------------------------------------------------------------------------
PetscErrorCode ReadMeshSegDir(
	FILE        *fp,
	const char  *name,
	PetscScalar  beg,
	PetscScalar  end,
	PetscInt    *tncels,
	MeshSegInp  *msi);

//---------------------------------------------------------------------------
#endif
