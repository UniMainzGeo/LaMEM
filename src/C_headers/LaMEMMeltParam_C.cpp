#include "LaMEMMeltParam_C.h"
#include "../LaMEM.h"
#include "../meltParam.h"

extern "C" {

// Global variable access
extern PetscScalar Pc;

// Convert C struct to C++ struct
static meltPar_Katz convertMeltPar(LaMEMMeltPar_Katz c_mp) {
    meltPar_Katz cpp_mp;
    cpp_mp.A1 = c_mp.A1; cpp_mp.A2 = c_mp.A2; cpp_mp.A3 = c_mp.A3;
    cpp_mp.B1 = c_mp.B1; cpp_mp.B2 = c_mp.B2; cpp_mp.B3 = c_mp.B3;
    cpp_mp.C1 = c_mp.C1; cpp_mp.C2 = c_mp.C2; cpp_mp.C3 = c_mp.C3;
    cpp_mp.r1 = c_mp.r1; cpp_mp.r2 = c_mp.r2;
    cpp_mp.beta1 = c_mp.beta1; cpp_mp.beta2 = c_mp.beta2;
    cpp_mp.K = c_mp.K; cpp_mp.gamma = c_mp.gamma;
    cpp_mp.D_water = c_mp.D_water; cpp_mp.lambda = c_mp.lambda;
    cpp_mp.chi1 = c_mp.chi1; cpp_mp.chi2 = c_mp.chi2;
    cpp_mp.Cp = c_mp.Cp; cpp_mp.DS = c_mp.DS;
    return cpp_mp;
}

// Convert C++ struct to C struct
static LaMEMMeltPar_Katz convertMeltParBack(meltPar_Katz cpp_mp) {
    LaMEMMeltPar_Katz c_mp;
    c_mp.A1 = cpp_mp.A1; c_mp.A2 = cpp_mp.A2; c_mp.A3 = cpp_mp.A3;
    c_mp.B1 = cpp_mp.B1; c_mp.B2 = cpp_mp.B2; c_mp.B3 = cpp_mp.B3;
    c_mp.C1 = cpp_mp.C1; c_mp.C2 = cpp_mp.C2; c_mp.C3 = cpp_mp.C3;
    c_mp.r1 = cpp_mp.r1; c_mp.r2 = cpp_mp.r2;
    c_mp.beta1 = cpp_mp.beta1; c_mp.beta2 = cpp_mp.beta2;
    c_mp.K = cpp_mp.K; c_mp.gamma = cpp_mp.gamma;
    c_mp.D_water = cpp_mp.D_water; c_mp.lambda = cpp_mp.lambda;
    c_mp.chi1 = cpp_mp.chi1; c_mp.chi2 = cpp_mp.chi2;
    c_mp.Cp = cpp_mp.Cp; c_mp.DS = cpp_mp.DS;
    return c_mp;
}

// Thin wrappers - convert structs and call existing functions
void LaMEMMeltParam_SetMeltParamsToDefault_Katz(LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp;
    setMeltParamsToDefault_Katz(&cpp_mp);
    *mp = convertMeltParBack(cpp_mp);
}

PetscScalar LaMEMMeltParam_MPgetFReactive(PetscScalar P, PetscScalar T, PetscScalar Cf, PetscScalar M, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return MPgetFReactive(P, T, Cf, M, &cpp_mp);
}

PetscScalar LaMEMMeltParam_MPgetFEquilib(PetscScalar P, PetscScalar T, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return MPgetFEquilib(P, T, X, M, &cpp_mp);
}

PetscScalar LaMEMMeltParam_MPgetTSolidus(PetscScalar P, PetscScalar X, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return MPgetTSolidus(P, X, &cpp_mp);
}

PetscScalar LaMEMMeltParam_MPgetTEquilib(PetscScalar P, PetscScalar F, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return MPgetTEquilib(P, F, X, M, &cpp_mp);
}

PetscScalar LaMEMMeltParam_MPgetFconsH(PetscScalar P, PetscScalar Ti, PetscScalar X, PetscScalar M, PetscScalar *Tf, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return MPgetFconsH(P, Ti, X, M, Tf, &cpp_mp);
}

// Helper functions
PetscScalar LaMEMMeltParam_CalcDT(PetscScalar P, PetscScalar X, PetscScalar F, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return calcDT(P, X, F, &cpp_mp);
}

PetscScalar LaMEMMeltParam_CalcF(PetscScalar T, PetscScalar dT, PetscScalar P, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return calcF(T, dT, P, Fcpx, &cpp_mp);
}

PetscScalar LaMEMMeltParam_FZero(PetscScalar F, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return FZero(F, T, P, X, Fcpx, &cpp_mp);
}

PetscScalar LaMEMMeltParam_HZero(PetscScalar F, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return HZero(F, T, P, X, M, &cpp_mp);
}

PetscScalar LaMEMMeltParam_FX_bal(PetscScalar x1, PetscScalar x2, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return FX_bal(x1, x2, T, P, X, Fcpx, &cpp_mp);
}

PetscScalar LaMEMMeltParam_FT_bal(PetscScalar x1, PetscScalar x2, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp) {
    meltPar_Katz cpp_mp = convertMeltPar(*mp);
    return FT_bal(x1, x2, T, P, X, M, &cpp_mp);
}

// Global parameter access
void LaMEMMeltParam_SetPc(PetscScalar pc_value) {
    Pc = pc_value;
}

PetscScalar LaMEMMeltParam_GetPc(void) {
    return Pc;
}

}