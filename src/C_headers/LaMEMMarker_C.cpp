#include "LaMEMMarker_C.h"
#include "../LaMEM.h"
#include "../marker.h"
#include "../parsing.h"
#include "../advect.h"
#include "../fdstag.h"
#include "../scaling.h"
#include "../JacRes.h"
#include "../phase.h"
#include "../tools.h"
#include "../bc.h"
#include "../surf.h"
#include "../phase_transition.h"


extern "C" {

// Thin wrappers - just cast and call existing functions
PetscErrorCode LaMEMMarker_ADVMarkInit(void *actx, void *fb) {
    return ADVMarkInit((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMMarker_ADVMarkInitCoord(void *actx) {
    return ADVMarkInitCoord((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkPerturb(void *actx) {
    return ADVMarkPerturb((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkSave(void *actx) {
    return ADVMarkSave((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkCheckMarkers(void *actx) {
    return ADVMarkCheckMarkers((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkSetTempGrad(void *actx) {
    return ADVMarkSetTempGrad((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkSetTempPhase(void *actx) {
    return ADVMarkSetTempPhase((AdvCtx*)actx);
}

PetscErrorCode LaMEMMarker_ADVMarkSetTempFile(void *actx, void *fb) {
    return ADVMarkSetTempFile((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMMarker_ADVMarkSetTempVector(void *actx) {
    return ADVMarkSetTempVector((AdvCtx*)actx);
}

// Specific initialization routines
PetscErrorCode LaMEMMarker_ADVMarkInitFiles(void *actx, void *fb) {
    return ADVMarkInitFiles((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMMarker_ADVMarkInitGeom(void *actx, void *fb) {
    return ADVMarkInitGeom((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMMarker_ADVMarkInitPolygons(void *actx, void *fb) {
    return ADVMarkInitPolygons((AdvCtx*)actx, (FB*)fb);
}

PetscErrorCode LaMEMMarker_ADVMarkReadCtrlPoly(void *fb, void *CtrlPoly, PetscInt *VolID, PetscInt *nCP) {
    return ADVMarkReadCtrlPoly((FB*)fb, (CtrlP*)CtrlPoly, *VolID, *nCP);
}

// Phase diagram functions
PetscErrorCode LaMEMMarker_LoadPhaseDiagram(void *actx, void *phases, PetscInt i) {
    return LoadPhaseDiagram((AdvCtx*)actx, (Material_t*)phases, i);
}

// Geometric primitives functions
void LaMEMMarker_SetPhaseSphere(void *sphere, void *P) {
    setPhaseSphere((GeomPrim*)sphere, (Marker*)P);
}

void LaMEMMarker_SetPhaseEllipsoid(void *ellipsoid, void *P) {
    setPhaseEllipsoid((GeomPrim*)ellipsoid, (Marker*)P);
}

void LaMEMMarker_SetPhaseBox(void *box, void *P) {
    setPhaseBox((GeomPrim*)box, (Marker*)P);
}

void LaMEMMarker_SetPhaseRidge(void *ridge, void *P) {
    setPhaseRidge((GeomPrim*)ridge, (Marker*)P);
}

void LaMEMMarker_SetPhaseLayer(void *layer, void *P) {
    setPhaseLayer((GeomPrim*)layer, (Marker*)P);
}

void LaMEMMarker_SetPhaseHex(void *hex, void *P) {
    setPhaseHex((GeomPrim*)hex, (Marker*)P);
}

void LaMEMMarker_SetPhaseCylinder(void *cylinder, void *P) {
    setPhaseCylinder((GeomPrim*)cylinder, (Marker*)P);
}

// Temperature computation
void LaMEMMarker_ComputeTemperature(void *geom, void *P, PetscScalar *T) {
    computeTemperature((GeomPrim*)geom, (Marker*)P, T);
}

// Utility functions
void LaMEMMarker_ADVMarkSecIdx(void *actx, PetscInt dir, PetscInt Islice, PetscInt *idx) {
    ADVMarkSecIdx((AdvCtx*)actx, dir, Islice, idx);
}

void LaMEMMarker_HexGetBoundingBox(PetscScalar *coord, PetscScalar *bounds) {
    HexGetBoundingBox(coord, bounds);
}

PetscInt LaMEMMarker_TetPointTest(PetscScalar *coord, PetscInt *ii, PetscScalar *xp, PetscScalar tol) {
    return TetPointTest(coord, ii, xp, tol);
}

}