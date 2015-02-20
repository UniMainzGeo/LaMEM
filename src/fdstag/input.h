//---------------------------------------------------------------------------
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#ifndef __input_h__
#define __input_h__
//---------------------------------------------------------------------------

PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user);

PetscErrorCode InputSetDefaultValues(UserCtx *user);

PetscErrorCode InputReadFile(JacRes *jr, UserCtx *user);

PetscErrorCode InputReadCommLine(UserCtx *user );

// old routines to input material properties
PetscErrorCode InitMaterialProp(UserCtx *user );
PetscErrorCode ReadMaterialProperties(UserCtx *user);

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
