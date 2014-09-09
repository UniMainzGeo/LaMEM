
#ifndef __Material_h__
#define __Material_h__

PetscErrorCode SetMaterialProperties( LaMEMVelPressureDA C, UserContext *user, DM DA_MG);


PetscErrorCode Average_EffectiveDensity_FDSTAG_Cell( const double VolumeFraction_phases[], UserContext *user,
		double T, double *EffectiveDensity_withAir, double *EffectiveDensity_withoutAir);

PetscErrorCode Average_EffectiveViscosity_FDSTAG_Cell( const double VolumeFraction_phases[], UserContext *user,
		double T, double E2nd, double T2nd, double P, double P_lithos,
		double PlasticStrain,
		double *EffectiveViscosity_withAir, double *EffectiveViscosity_withoutAir);

PetscErrorCode Average_EffectiveThermalProperties_FDSTAG_Cell(const double VolumeFraction_phases[], UserContext *user,
		PetscScalar *Effective_HeatCapacity_withAir, 		PetscScalar *Effective_HeatCapacity_withoutAir,
		PetscScalar *Effective_ThermalConductivity_withAir, PetscScalar *Effective_ThermalConductivity_withoutAir,
		PetscScalar *Effective_RadioactiveHeat_withAir, 	PetscScalar *Effective_RadioactiveHeat_withoutAir);


PetscErrorCode ComputeEffectiveViscosity( PetscInt phase, UserContext *user, double T, double E2nd, double T2nd, double P, double P_lithos,
		double PlasticStrain, PetscInt PlasticityCutoff, double *mu_viscous, double *mu_plastic, double *mu_eff);

PetscErrorCode ComputeEffectiveDensity( PetscInt phase, UserContext *user, double T, double *rho_eff);

PetscErrorCode ComputePointWiseProperties( LaMEMVelPressureDA C, PointWiseInformation *PointInformation, double Point[3], double V_element[],
		double P_element[], double T_element[], double T_old_element[], double DevStress[6], DMDACoor3d coord_elem[], double mu, double mu_viscous,
		PetscInt ComputeFull);

PetscErrorCode AverageViscosityDensityPerElement( double mu[MAX_ngp_vel], double mu_v[MAX_ngp_vel], double mu_p[MAX_ngp_vel], double *mu_average, double rho[MAX_ngp_vel], UserContext *user, PetscInt ngp_vel);

PetscErrorCode FDSTAG_DominantPhaseAtCorners(UserContext *user, PetscInt ielx, PetscInt iely, PetscInt ielz, PetscInt *DominantPhase_WithAir, PetscInt *DominantPhase_WithoutAir, PetscInt *DominantPhase_WithoutAir_TopCell);

PetscErrorCode LaMEMSetMaterialDataMemoryFromArray(
		MaterialsElementDynamic *data,
		const PetscInt i, const PetscInt j, const PetscInt k,
		const PetscInt ngp,
		PetscScalar ***list );

PetscErrorCode ResetStressesBeforeNonlinearIterations( LaMEMVelPressureDA C, UserContext *user);

PetscErrorCode UpdateMaterialProperties_FDSTAG(UserContext *user);




#endif
