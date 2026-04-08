#include "LaMEMObjFunct_C.h"
#include "../LaMEM.h"
#include "../objFunct.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../fdstag.h"
#include "../surf.h"
#include "../JacRes.h"
#include "../tools.h"

extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMObjFunct_ObjFunctDestroy(void *objf) {
    return ObjFunctDestroy((ObjFunct*)objf);
}

PetscErrorCode LaMEMObjFunct_ObjFunctCreate(void *objf, void *IOparam, void *surf, void *fb) {
    return ObjFunctCreate((ObjFunct*)objf, (ModParam*)IOparam, (FreeSurf*)surf, (FB*)fb);
}

PetscErrorCode LaMEMObjFunct_ObjFunctReadFromOptions(void *objf, const char *on[], void *fb) {
    return ObjFunctReadFromOptions((ObjFunct*)objf, on, (FB*)fb);
}

PetscErrorCode LaMEMObjFunct_VecErrSurf(Vec mod, void *objf, PetscInt field, PetscScalar scal) {
    return VecErrSurf(mod, (ObjFunct*)objf, field, scal);
}

PetscErrorCode LaMEMObjFunct_ObjFunctCompErr(void *objf) {
    return ObjFunctCompErr((ObjFunct*)objf);
}

// Helper function to set observation type names
void LaMEMObjFunct_SetObservationNames(const char **on) {
    static const char *velx_name = "velx";
    static const char *vely_name = "vely";
    static const char *velz_name = "velz";
    static const char *topo_name = "topo";
    static const char *boug_name = "boug";
    static const char *isa_name = "isa";
    static const char *shmax_name = "shmax";
    
    on[LAMEM_VELX] = velx_name;
    on[LAMEM_VELY] = vely_name;
    on[LAMEM_VELZ] = velz_name;
    on[LAMEM_TOPO] = topo_name;
    on[LAMEM_BOUG] = boug_name;
    on[LAMEM_ISA] = isa_name;
    on[LAMEM_SHMAX] = shmax_name;
}

}