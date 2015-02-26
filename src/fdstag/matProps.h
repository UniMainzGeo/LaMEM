//---------------------------------------------------------------------------
//................ Routines related to Material Properties ..................
//---------------------------------------------------------------------------
#ifndef __matProps_h__
#define __matProps_h__
//---------------------------------------------------------------------------
// allows some overhead when reading
#define _max_overhead_ 5
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//........................... MATERIAL PARAMETERS ...........................
//---------------------------------------------------------------------------
// main function
PetscErrorCode MatPropInit(JacRes *jr, UserCtx *usr);

// get material properties structure in file
PetscErrorCode MatPropGetStruct(FILE *fp, Material_t *m, PetscInt ils, PetscInt ile);

// error checking
PetscErrorCode MatPropErrorCheck(PetscInt id, PetscInt err);

// set default values for material parameters
void MatPropSet(Material_t *m, PetscInt dim);

// print info
void MatPropPrint(Material_t *m, PetscScalar eta);

//---------------------------------------------------------------------------
//............................ SOFTENING LAWS ...............................
//---------------------------------------------------------------------------
// main function
PetscErrorCode MatSoftInit(JacRes *jr, UserCtx *usr);

// get softening laws structure in file
PetscErrorCode MatSoftGetStruct( FILE *fp, Soft_t *s, PetscInt ils, PetscInt ile);

// points every phase to the correct softening laws
PetscErrorCode MatSoftSet(JacRes *jr);

//---------------------------------------------------------------------------
//................ Routines to get structure-info from file .................
//---------------------------------------------------------------------------
// gets the file positions of a structure
void getLineStruct( FILE *fp, PetscInt *ls, PetscInt *le, PetscInt *count, PetscInt *count1, const char key[], const char key_end[]);

// gets an integer within specified positions in file
void getMatPropInt   (FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscInt *value, PetscInt *found);

// gets a scalar within specified positions in file
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscScalar *value, PetscInt *found);

#endif
