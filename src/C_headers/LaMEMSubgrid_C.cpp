#include "LaMEMSubgrid_C.h"
#include "../LaMEM.h"
#include "../subgrid.h"
#include "../advect.h"
#include "../phase.h"
#include "../parsing.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../fdstag.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../surf.h"
#include "../marker.h"
#include "../AVD.h"
#include "../cvi.h"
#include "../Tensor.h"
#include "../tools.h"
#include "../phase_transition.h"


extern "C" {

// Main subgrid marker management functions
PetscErrorCode LaMEMSubgrid_ADVMarkSubGrid(void *actx) {
    return ADVMarkSubGrid((AdvCtx*)actx);
}

PetscErrorCode LaMEMSubgrid_ADVMarkClone(void *actx, PetscInt icell, PetscInt isubcell, 
                                         PetscScalar s[3], PetscScalar h[3], 
                                         void *dist, void *iclone) {
    return ADVMarkClone((AdvCtx*)actx, icell, isubcell, s, h, 
                        *(std::vector<spair>*)dist, *(std::vector<Marker>*)iclone);
}

PetscErrorCode LaMEMSubgrid_ADVMarkCheckMerge(void *actx, PetscInt ib, PetscInt ie, 
                                              PetscInt *nmerge, void *mark, void *cell, 
                                              void *iclone, void *imerge) {
    return ADVMarkCheckMerge((AdvCtx*)actx, ib, ie, *nmerge, 
                             *(std::vector<Marker>*)mark, *(std::vector<ipair>*)cell,
                             *(std::vector<Marker>*)iclone, *(std::vector<PetscInt>*)imerge);
}

PetscErrorCode LaMEMSubgrid_ADVMarkMerge(void *mark, PetscInt nmark, PetscInt npmax, 
                                         PetscInt *sz) {
    return ADVMarkMerge(*(std::vector<Marker>*)mark, nmark, npmax, *sz);
}

PetscErrorCode LaMEMSubgrid_ADVCollectGarbageVec(void *actx, void *recvbuf, void *idel) {
    return ADVCollectGarbageVec((AdvCtx*)actx, *(std::vector<Marker>*)recvbuf, 
                                *(std::vector<PetscInt>*)idel);
}

// Free surface interaction functions
PetscErrorCode LaMEMSubgrid_ADVMarkCrossFreeSurf(void *actx) {
    return ADVMarkCrossFreeSurf((AdvCtx*)actx);
}

PetscErrorCode LaMEMSubgrid_ADVGetSedPhase(void *actx, Vec vphase) {
    return ADVGetSedPhase((AdvCtx*)actx, vphase);
}

}