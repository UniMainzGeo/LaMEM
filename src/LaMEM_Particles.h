
#ifndef __LaMEM_Particles_h__
#define __LaMEM_Particles_h__

PetscErrorCode InitializeParticles( DM da, UserContext *user );

PetscErrorCode InitializeParticles_ElementWise(LaMEMVelPressureDA C, DM da, UserContext *user );

PetscErrorCode ElementInfo( DMDACoor3d coord_elem[], PetscScalar *minx, PetscScalar *maxx, PetscScalar *miny,
		PetscScalar *maxy, PetscScalar *minz, PetscScalar *maxz, PetscScalar *max_leftx, PetscScalar *min_rightx,
		PetscScalar *max_fronty, PetscScalar *min_backy, PetscScalar *max_botz, PetscScalar *min_topz );

PetscErrorCode Inside_8NodeElement( PetscScalar x_real, PetscScalar y_real, PetscScalar z_real,
		PetscScalar *eta_out, PetscScalar *zetha_out, PetscScalar *phi_out, DMDACoor3d coord_elem[],	PetscInt *inside, PetscInt DisplayIterations);

PetscErrorCode NaturalCoordinates( const PetscInt _nnel, PetscScalar x_real, PetscScalar y_real, PetscScalar z_real,
		PetscScalar *eta_out, PetscScalar *zetha_out, PetscScalar *phi_out, DMDACoor3d coord_elem[], PetscInt *inside );

PetscErrorCode GetParticleNodes( LaMEMVelPressureDA C, DM da, UserContext *user);

PetscErrorCode AdvectParticles( LaMEMVelPressureDA C, DM da, UserContext *user, Vec Velocity, PetscScalar dt );

PetscErrorCode InterpolateParticles( LaMEMVelPressureDA C, DM da, UserContext *user);

PetscErrorCode SetInitialTracerProperties( DM da, UserContext *user );

PetscErrorCode SetInitialTracerPhasesFromFile( DM da, UserContext *user );

PetscErrorCode MaterialPropertiesFromTracers( LaMEMVelPressureDA C, UserContext *user);

PetscErrorCode OutputVTK_ViscosityDensityFields_FDSTAG( UserContext *user, Vec Temp);

PetscErrorCode WriteInitialParticlesToDisc( UserContext *user);
PetscErrorCode LoadInitialParticlesFromDisc( UserContext *user);

PetscErrorCode WriteParticlesToDisc( UserContext *user, PetscInt itime );

PetscErrorCode ParticlePhaseTransitions( UserContext *user );

PetscErrorCode ComputePropertiesAtParticles(LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Pressure,
		DM da_temp, Vec Temp, Vec Temp_old , UserContext *user, PetscScalar dt );

PetscErrorCode ComputeHistoryPropertiesFromParticles_FDSTAG(LaMEMVelPressureDA C, UserContext *user);
PetscErrorCode AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(Vec GlobalVector, UserContext *user, const char PropertyToBeAveraged[], const char NodeAtWhichAveragingIsPerformed[]);

PetscErrorCode GetParticleElementAndLocalCoordinates(
		PetscInt *ix_in, 	PetscInt *iy_in, 	PetscInt *iz_in,		PetscInt xs_g, 		PetscInt ys_g, 	PetscInt zs_g,
		PetscInt  xm_g, 	PetscInt ym_g, 		PetscInt zm_g, 			PetscInt xs, 		PetscInt ys, 	PetscInt zs,
		PetscInt  xm,		PetscInt ym, 		PetscInt zm,			DMDACoor3d ***coords, Particles 		*ParticleLocal_in,
		UserContext *user, 				DAVPElementType vpt_element_type,	PetscInt nx, 		PetscInt ny,
		PetscInt 	*TracerDeactivate_in, 		PetscInt *ElementFound_Outside,
		PetscInt *cpu_x_in, PetscInt *cpu_y_in, PetscInt *cpu_z_in,  PetscInt *DisplayOutput);

PetscErrorCode SendParticleToCorrectCPU(PetscInt NumParticlesMove[MaxNumCPU],  PetscInt StorageSpaceNeighborCPU[3][3][3], PetscInt NeighborCPU[3][3][3], Particles **ParticlesSendAway,
		UserContext *user, PetscInt DisplayOutput);

PetscErrorCode CorrectPhasesForInternalFreeSurface(UserContext *user );

PetscErrorCode LoadInitialParticlesFromDisc_FDSTAG( UserContext *user);

#endif
