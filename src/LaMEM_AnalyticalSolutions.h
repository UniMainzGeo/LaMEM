
#ifndef __LaMEM_AnalyticalSolutions_h__
#define __LaMEM_AnalyticalSolutions_h__

void solCx(
		double pos[],
		double _eta_A, double _eta_B,
		double _x_c, PetscInt _n,
		double vel[], double* presssure,
		double total_stress[], double strain_rate[] );


PetscErrorCode LaMEM_Initialize_AnalyticalSolution( UserContext *user,  LaMEMVelPressureDA C);


PetscErrorCode LaMEM_CompareNumerics_vs_AnalyticalSolution( UserContext *user,  LaMEMVelPressureDA C, Vec sol, Vec Pressure);


#endif
