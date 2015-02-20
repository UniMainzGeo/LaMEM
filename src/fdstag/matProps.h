//---------------------------------------------------------------------------
//.................Routines related to  Material Properties .................
//---------------------------------------------------------------------------
#ifndef __matProps_h__
#define __matProps_h__
//---------------------------------------------------------------------------
// allows some overhead when reading
#define _max_over_phases_ 5
//---------------------------------------------------------------------------

// new routine to input material properties
PetscErrorCode MatPropInit(JacRes *jr, UserCtx *usr);

PetscErrorCode SetMatSoftening(JacRes *jr, UserCtx *usr);

PetscErrorCode MatPropGetStruct( FILE *fp, Material_t *m, PetscInt ils, PetscInt ile);

PetscErrorCode MatPropErrorCheck(PetscInt id, PetscInt err);

void MatPropSet(Material_t *m, PetscInt dim);

void MatPropPrint(Material_t *m, PetscScalar eta0);

// routines to get info from file
void getLineStruct   (FILE *fp, PetscInt *ls, PetscInt *le, PetscInt *count, PetscInt *count1);
void getMatPropInt   (FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscInt *value, PetscInt *found );
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscScalar *value, PetscInt *found );

#endif
