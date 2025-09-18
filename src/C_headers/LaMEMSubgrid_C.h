#ifndef LAMEMSUBGRID_C_H
#define LAMEMSUBGRID_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Main subgrid marker management functions
PetscErrorCode LaMEMSubgrid_ADVMarkSubGrid(void *actx);
PetscErrorCode LaMEMSubgrid_ADVMarkClone(void *actx, PetscInt icell, PetscInt isubcell, 
                                         PetscScalar s[3], PetscScalar h[3], 
                                         void *dist, void *iclone);
PetscErrorCode LaMEMSubgrid_ADVMarkCheckMerge(void *actx, PetscInt ib, PetscInt ie, 
                                              PetscInt *nmerge, void *mark, void *cell, 
                                              void *iclone, void *imerge);
PetscErrorCode LaMEMSubgrid_ADVMarkMerge(void *mark, PetscInt nmark, PetscInt npmax, 
                                         PetscInt *sz);
PetscErrorCode LaMEMSubgrid_ADVCollectGarbageVec(void *actx, void *recvbuf, void *idel);

// Free surface interaction functions
PetscErrorCode LaMEMSubgrid_ADVMarkCrossFreeSurf(void *actx);
PetscErrorCode LaMEMSubgrid_ADVGetSedPhase(void *actx, Vec vphase);

#ifdef __cplusplus
}
#endif

#endif