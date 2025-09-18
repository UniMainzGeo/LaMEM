#ifndef LAMEMOBJFUNCT_C_H
#define LAMEMOBJFUNCT_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Observation type constants
#define LAMEM_VELX  0
#define LAMEM_VELY  1  
#define LAMEM_VELZ  2
#define LAMEM_TOPO  3
#define LAMEM_BOUG  4
#define LAMEM_ISA   5
#define LAMEM_SHMAX 6

// Thin wrappers around actual objective function computation functions
PetscErrorCode LaMEMObjFunct_ObjFunctDestroy(void *objf);
PetscErrorCode LaMEMObjFunct_ObjFunctCreate(void *objf, void *IOparam, void *surf, void *fb);
PetscErrorCode LaMEMObjFunct_ObjFunctReadFromOptions(void *objf, const char *on[], void *fb);
PetscErrorCode LaMEMObjFunct_VecErrSurf(Vec mod, void *objf, PetscInt field, PetscScalar scal);
PetscErrorCode LaMEMObjFunct_ObjFunctCompErr(void *objf);

// Helper function to set observation type names
void LaMEMObjFunct_SetObservationNames(const char **on);

#ifdef __cplusplus
}
#endif

#endif