#include "LaMEMAdjoint_C.h"
#include "../LaMEM.h"
#include "../adjoint.h"
#include "../phase.h"
#include "../tools.h"
#include "../fdstag.h"
#include "../scaling.h"
#include "../tssolve.h"
#include "../bc.h"
#include "../JacRes.h"
#include "../interpolate.h"
#include "../surf.h"
#include "../multigrid.h"
#include "../matrix.h"
#include "../lsolve.h"
#include "../nlsolve.h"
#include "../objFunct.h"
#include "../constEq.h"
#include "../parsing.h"
#include "../gravity.h"

extern "C" {

// Main adjoint functions
PetscErrorCode LaMEMAdjoint_Main(void *IOparam) {
    return LaMEMAdjointMain((ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_OptimizationTAO(void *tao, void *P, PetscReal *F, void *grad, void *ctx) {
    return AdjointOptimisationTAO((Tao)tao, (Vec)P, F, (Vec)grad, ctx);
}

PetscErrorCode LaMEMAdjoint_Optimization(void *P, PetscScalar F, void *grad, void *ctx) {
    return AdjointOptimisation((Vec)P, F, (Vec)grad, ctx);
}

// Adjoint gradient object creation/destruction
PetscErrorCode LaMEMAdjoint_Create(void *aop, void *jr, void *IOparam) {
    return AdjointCreate((AdjGrad*)aop, (JacRes*)jr, (ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_Destroy(void *aop, void *IOparam) {
    return AdjointDestroy((AdjGrad*)aop, (ModParam*)IOparam);
}

// Adjoint vectors creation/destruction
PetscErrorCode LaMEMAdjoint_VectorsCreate(void *Adjoint_vectors, void *IOparam) {
    return AdjointVectorsCreate((Adjoint_Vecs*)Adjoint_vectors, (ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_VectorsDestroy(void *Adjoint_vectors, void *IOparam) {
    return AdjointVectorsDestroy((Adjoint_Vecs*)Adjoint_vectors, (ModParam*)IOparam);
}

// Core gradient computation functions
PetscErrorCode LaMEMAdjoint_ObjectiveAndGradientFunction(void *aop, void *jr, void *nl, void *IOparam, void *snes, void *surf) {
    return AdjointObjectiveAndGradientFunction((AdjGrad*)aop, (JacRes*)jr, (NLSol*)nl, (ModParam*)IOparam, (SNES)snes, (FreeSurf*)surf);
}

PetscErrorCode LaMEMAdjoint_ComputeGradientsAndObjectiveFunction(void *Parameters, PetscScalar *ObjectiveValue, void *Gradient, void *IOparam) {
    return ComputeGradientsAndObjectiveFunction((Vec)Parameters, ObjectiveValue, (Vec)Gradient, (ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_ComputeGradients(void *jr, void *aop, void *nl, void *snes, void *IOparam) {
    return AdjointComputeGradients((JacRes*)jr, (AdjGrad*)aop, (NLSol*)nl, (SNES)snes, (ModParam*)IOparam);
}

// Objective function computation
PetscErrorCode LaMEMAdjoint_ObjectiveFunction(void *aop, void *jr, void *IOparam, void *surf) {
    return AdjointObjectiveFunction((AdjGrad*)aop, (JacRes*)jr, (ModParam*)IOparam, (FreeSurf*)surf);
}

// Finite difference gradients
PetscErrorCode LaMEMAdjoint_FiniteDifferenceGradients(void *IOparam) {
    return AdjointFiniteDifferenceGradients((ModParam*)IOparam);
}

// Input/Output functions
PetscErrorCode LaMEMAdjoint_ReadInputSetDefaults(void *IOparam, void *Adjoint_Vectors) {
    return LaMEMAdjointReadInputSetDefaults((ModParam*)IOparam, (Adjoint_Vecs*)Adjoint_Vectors);
}

PetscErrorCode LaMEMAdjoint_PrintGradientsAndObservationPoints(void *IOparam) {
    return PrintGradientsAndObservationPoints((ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_PrintCostFunction(void *IOparam) {
    return PrintCostFunction((ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_PrintScalingLaws(void *IO_param) {
    return PrintScalingLaws((ModParam*)IO_param);
}

// Material parameter scanning and modification
PetscErrorCode LaMEMAdjoint_ScanForMaterialParameters(void *fb, void *scaling, PetscInt *numParam, char parameterType, PetscInt *ParamID, PetscScalar *ParamValue, PetscInt *PhaseID, PetscScalar *Bounds) {
    return Adjoint_ScanForMaterialParameters((FB*)fb, (Scaling*)scaling, numParam, parameterType, ParamID, ParamValue, PhaseID, Bounds);
}

PetscErrorCode LaMEMAdjoint_AddMaterialParameterToCommandLineOptions(char *name, PetscInt ID, PetscScalar val) {
    return AddMaterialParameterToCommandLineOptions(name, ID, val);
}

PetscErrorCode LaMEMAdjoint_DeleteMaterialParameterFromCommandLineOptions(char *name, PetscInt ID) {
    return DeleteMaterialParameterFromCommandLineOptions(name, ID);
}

PetscErrorCode LaMEMAdjoint_CreateModifiedMaterialDatabase(void *IOparam) {
    return CreateModifiedMaterialDatabase((ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_CopyParameterToLaMEMCommandLine(void *IOparam, PetscScalar CurVal, PetscInt j) {
    return CopyParameterToLaMEMCommandLine((ModParam*)IOparam, CurVal, j);
}

// Boundary conditions
PetscErrorCode LaMEMAdjoint_ApplyBCs(void *dF, void *bc) {
    return Adjoint_ApplyBCs((Vec)dF, (BCCtx*)bc);
}

// Projection and interpolation
PetscErrorCode LaMEMAdjoint_PointInPro(void *jr, void *aop, void *IOparam, void *surf) {
    return AdjointPointInPro((JacRes*)jr, (AdjGrad*)aop, (ModParam*)IOparam, (FreeSurf*)surf);
}

// Finite difference specific functions
PetscErrorCode LaMEMAdjoint_Get_F_dFdu_Center(void *jr, void *aop, void *IOparam) {
    return AdjointGet_F_dFdu_Center((JacRes*)jr, (AdjGrad*)aop, (ModParam*)IOparam);
}

PetscErrorCode LaMEMAdjoint_GradientResetParameter(void *nl, PetscInt CurPar, PetscInt CurPhase, void *aop) {
    return AdjointGradientResetParameter((NLSol*)nl, CurPar, CurPhase, (AdjGrad*)aop);
}

PetscErrorCode LaMEMAdjoint_FormResidualFieldFD(void *snes, void *x, void *psi, void *nl, void *aop, void *IOparam) {
    return AdjointFormResidualFieldFD((SNES)snes, (Vec)x, (Vec)psi, (NLSol*)nl, (AdjGrad*)aop, (ModParam*)IOparam);
}

// Constitutive equation finite difference functions
PetscErrorCode LaMEMAdjoint_devConstEqFD(void *ctx, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk) {
    return devConstEqFD((ConstEqCtx*)ctx, (AdjGrad*)aop, (ModParam*)IOparam, ii, jj, k, ik, jk, kk);
}

PetscErrorCode LaMEMAdjoint_cellConstEqFD(void *ctx, void *svCell, PetscScalar dxx, PetscScalar dyy, PetscScalar dzz, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz, PetscScalar *gres, PetscScalar *rho, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk) {
    return cellConstEqFD((ConstEqCtx*)ctx, (SolVarCell*)svCell, dxx, dyy, dzz, *sxx, *syy, *szz, *gres, *rho, (AdjGrad*)aop, (ModParam*)IOparam, ii, jj, k, ik, jk, kk);
}

PetscErrorCode LaMEMAdjoint_setUpPhaseFD(void *ctx, PetscInt ID, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk) {
    return setUpPhaseFD((ConstEqCtx*)ctx, ID, (AdjGrad*)aop, (ModParam*)IOparam, ii, jj, k, ik, jk, kk);
}

PetscErrorCode LaMEMAdjoint_edgeConstEqFD(void *ctx, void *svEdge, PetscScalar d, PetscScalar *s, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk) {
    return edgeConstEqFD((ConstEqCtx*)ctx, (SolVarEdge*)svEdge, d, *s, (AdjGrad*)aop, (ModParam*)IOparam, ii, jj, k, ik, jk, kk);
}

// Utility functions
PetscErrorCode LaMEMAdjoint_Parameter_SetFDgrad_Option(PetscInt *FD_grad, char *name) {
    return Parameter_SetFDgrad_Option(FD_grad, name);
}

PetscErrorCode LaMEMAdjoint_swapStruct(void *A, void *B) {
    return swapStruct((struct Material_t*)A, (struct Material_t*)B);
}

}