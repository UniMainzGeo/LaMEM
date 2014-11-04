//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#ifndef __input_h__
#define __input_h__
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInitCode(JacRes *jr, UserContext *user);

PetscErrorCode FDSTAGReadInputFile(JacRes *jr, UserContext *user);

PetscErrorCode ReadMaterialProperties(UserContext *user);

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
