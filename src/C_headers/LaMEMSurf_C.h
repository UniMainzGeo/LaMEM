#ifndef LAMEMSURF_C_H
#define LAMEMSURF_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Free surface setup and initialization
PetscErrorCode LaMEMSurf_FreeSurfCreate(void *surf, void *fb);
PetscErrorCode LaMEMSurf_FreeSurfCreateData(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfDestroy(void *surf);

// Free surface evolution and advection
PetscErrorCode LaMEMSurf_FreeSurfAdvect(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfAdvectTopo(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfGetVelComp(void *surf, 
                                            PetscErrorCode (*interp)(void*, Vec, Vec, void*),
                                            Vec vcomp_grid, Vec vcomp_surf);

// Free surface smoothing and processing
PetscErrorCode LaMEMSurf_FreeSurfSmoothMaxAngle(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfGetAirPhaseRatio(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfGetAvgTopo(void *surf);

// Erosion and sedimentation
PetscErrorCode LaMEMSurf_FreeSurfAppErosion(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfAppSedimentation(void *surf);

// Initial conditions and file I/O
PetscErrorCode LaMEMSurf_FreeSurfSetInitialPerturbation(void *surf);
PetscErrorCode LaMEMSurf_FreeSurfSetTopoFromFile(void *surf, void *fb);

// Restart I/O
PetscErrorCode LaMEMSurf_FreeSurfReadRestart(void *surf, FILE *fp);
PetscErrorCode LaMEMSurf_FreeSurfWriteRestart(void *surf, FILE *fp);

// Utility functions
PetscInt LaMEMSurf_InterpolateTriangle(PetscScalar *x, PetscScalar *y, PetscScalar *f, 
                                       PetscInt *i, PetscScalar xp, PetscScalar yp, 
                                       PetscScalar tol, PetscScalar *fp);
PetscScalar LaMEMSurf_IntersectTriangularPrism(PetscScalar *x, PetscScalar *y, PetscScalar *z, 
                                               PetscInt *i, PetscScalar vcell, 
                                               PetscScalar bot, PetscScalar top, PetscScalar tol);

#ifdef __cplusplus
}
#endif

#endif