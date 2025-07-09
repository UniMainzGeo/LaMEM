#include "LaMEMSurf_C.h"
#include "../LaMEM.h"
#include "../surf.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../phase.h"
#include "../JacRes.h"
#include "../interpolate.h"
#include "../tools.h"

extern "C" {

// Free surface setup and initialization
PetscErrorCode LaMEMSurf_FreeSurfCreate(void *surf, void *fb) {
    return FreeSurfCreate((FreeSurf*)surf, (FB*)fb);
}

PetscErrorCode LaMEMSurf_FreeSurfCreateData(void *surf) {
    return FreeSurfCreateData((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfDestroy(void *surf) {
    return FreeSurfDestroy((FreeSurf*)surf);
}

// Free surface evolution and advection
PetscErrorCode LaMEMSurf_FreeSurfAdvect(void *surf) {
    return FreeSurfAdvect((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfAdvectTopo(void *surf) {
    return FreeSurfAdvectTopo((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfGetVelComp(void *surf, 
                                            PetscErrorCode (*interp)(void*, Vec, Vec, void*),
                                            Vec vcomp_grid, Vec vcomp_surf) {
    return FreeSurfGetVelComp((FreeSurf*)surf, 
                              (PetscErrorCode (*)(FDSTAG*, Vec, Vec, InterpFlags))interp,
                              vcomp_grid, vcomp_surf);
}

// Free surface smoothing and processing
PetscErrorCode LaMEMSurf_FreeSurfSmoothMaxAngle(void *surf) {
    return FreeSurfSmoothMaxAngle((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfGetAirPhaseRatio(void *surf) {
    return FreeSurfGetAirPhaseRatio((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfGetAvgTopo(void *surf) {
    return FreeSurfGetAvgTopo((FreeSurf*)surf);
}

// Erosion and sedimentation
PetscErrorCode LaMEMSurf_FreeSurfAppErosion(void *surf) {
    return FreeSurfAppErosion((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfAppSedimentation(void *surf) {
    return FreeSurfAppSedimentation((FreeSurf*)surf);
}

// Initial conditions and file I/O
PetscErrorCode LaMEMSurf_FreeSurfSetInitialPerturbation(void *surf) {
    return FreeSurfSetInitialPerturbation((FreeSurf*)surf);
}

PetscErrorCode LaMEMSurf_FreeSurfSetTopoFromFile(void *surf, void *fb) {
    return FreeSurfSetTopoFromFile((FreeSurf*)surf, (FB*)fb);
}

// Restart I/O
PetscErrorCode LaMEMSurf_FreeSurfReadRestart(void *surf, FILE *fp) {
    return FreeSurfReadRestart((FreeSurf*)surf, fp);
}

PetscErrorCode LaMEMSurf_FreeSurfWriteRestart(void *surf, FILE *fp) {
    return FreeSurfWriteRestart((FreeSurf*)surf, fp);
}

// Utility functions
PetscInt LaMEMSurf_InterpolateTriangle(PetscScalar *x, PetscScalar *y, PetscScalar *f, 
                                       PetscInt *i, PetscScalar xp, PetscScalar yp, 
                                       PetscScalar tol, PetscScalar *fp) {
    return InterpolateTriangle(x, y, f, i, xp, yp, tol, fp);
}

PetscScalar LaMEMSurf_IntersectTriangularPrism(PetscScalar *x, PetscScalar *y, PetscScalar *z, 
                                               PetscInt *i, PetscScalar vcell, 
                                               PetscScalar bot, PetscScalar top, PetscScalar tol) {
    return IntersectTriangularPrism(x, y, z, i, vcell, bot, top, tol);
}

}