//---------------------------------------------------------------------------
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#ifndef __matProps_h__
#define __matProps_h__
//---------------------------------------------------------------------------
//........................... MATERIAL PARAMETERS ...........................
//---------------------------------------------------------------------------

// read all phases
PetscErrorCode MatPropInit(JacRes *jr, UserCtx *usr);

// read single phase
PetscErrorCode MatPropGetStruct(FILE *fp,
	PetscInt numPhases, Material_t *phases,
	PetscInt numSoft,   Soft_t     *matSoft,
	PetscInt ils, PetscInt ile, UnitsType utype);

//---------------------------------------------------------------------------
//............................ SOFTENING LAWS ...............................
//---------------------------------------------------------------------------

// read all softening laws
PetscErrorCode MatSoftInit(JacRes *jr, UserCtx *usr);

// read single softening laws
PetscErrorCode MatSoftGetStruct(FILE *fp,
	PetscInt numSoft, Soft_t *matSoft,
	PetscInt ils, PetscInt ile);

//---------------------------------------------------------------------------
//................ Routines to get structure-info from file .................
//---------------------------------------------------------------------------

// gets the file positions of a structure
void getLineStruct(
	FILE *fp, PetscInt *ls, PetscInt *le, PetscInt mux_num,
	PetscInt *count_starts, PetscInt *count_ends,
	const char key[], const char key_end[]);

// gets an integer within specified positions in file
void getMatPropInt(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscInt *value, PetscInt *found);

// gets a scalar within specified positions in file
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscScalar *value, PetscInt *found);

//---------------------------------------------------------------------------

#endif
