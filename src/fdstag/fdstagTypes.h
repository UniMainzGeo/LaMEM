// LaMEM user data structures definitions

#ifndef __fdstagTypes_h__
#define __fdstagTypes_h__
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar    Vy_front, Vy_back, Vx_left, Vx_right, Vz_bot, Vz_top, Exx, Eyy;
	PetscInt       UpperBound, LowerBound, LeftBound, RightBound, FrontBound, BackBound, InternalBound;
	DMBoundaryType BCType_x;
	DMBoundaryType BCType_y;
	DMBoundaryType BCType_z;
} SBC;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Defines physical properties for each of the phases
typedef struct {
	PetscInt	ViscosityLaw[max_num_phases],		DensityLaw[max_num_phases],	  	  PlasticityLaw[max_num_phases];
	PetscScalar mu[max_num_phases], 	  			rho[max_num_phases], 			  n_exponent[max_num_phases];
	PetscScalar A[max_num_phases], E[max_num_phases];
	PetscScalar	ElasticShearModule[max_num_phases],	ElasticBulkModule[max_num_phases];
	PetscScalar Cohesion[max_num_phases],		  FrictionAngle[max_num_phases];
	PetscScalar	T_Conductivity[max_num_phases],		HeatCapacity[max_num_phases],	  RadioactiveHeat[max_num_phases];
	PetscScalar	ThermalExpansivity[max_num_phases], FrankKamenetskii[max_num_phases], Density_T0[max_num_phases];
	PetscScalar Powerlaw_e0[max_num_phases];
	PetscScalar	CohesionAfterWeakening[max_num_phases];
	PetscScalar FrictionAngleAfterWeakening[max_num_phases];
	PetscScalar	Weakening_PlasticStrain_Begin[max_num_phases];
	PetscScalar Weakening_PlasticStrain_End[max_num_phases];
	PetscScalar	Ra[max_num_phases];
} PhProps;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Structure that holds characteristic values for nondimensionalisation
typedef struct {
	PetscScalar Length, Time, Stress, Velocity, Temperature, Viscosity;
	PetscScalar Density, kg, Strainrate, ThermalExpansivity, km, SecYear, cmYear, Myrs, MPa, Force, Watt;
	PetscScalar T_conductivity, RadioactiveHeat, Joule, HeatCapacity, Jmol;
} nonDimUnits;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Structure that holds gravity parameters
typedef struct {
	PetscInt		GetIt;
	PetscInt		SaveDebug,SaveVTK,SaveRef;
	PetscBool		UseNumerics, UseAnalytics;
	PetscInt		survey_nx, survey_ny;
	PetscScalar		survey_xs, survey_xm;
	PetscScalar		survey_ys, survey_ym;
	PetscScalar		survey_z ;
	PetscScalar     ReferenceDensity;
	PetscScalar     StdDev;
	PetscScalar		LithColDens[9],LithColDepth[8];
	PetscInt		num_intp,LithColNum;
	char			RefDatFile2load[PETSC_MAX_PATH_LEN];
} gravityParams;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Structure that holds pushing parameters
typedef struct {
	PetscScalar		L_block, W_block, H_block;
	PetscScalar		x_center_block, y_center_block, z_center_block;
	PetscScalar     V_push[10], omega[10];
	PetscInt		num_changes, reset_pushing_coord, ind_change;
	PetscScalar		time[11];
	PetscInt	    coord_advect[10], dir[10];
	PetscScalar		theta;
	Vec 			PV_rhs;
} PushParams;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Mesh segments input data structures
typedef struct
{
	PetscInt    nsegs;                    // number of segments
	PetscScalar delims[MaxNumMeshSegs-1]; // coordinates of the delimiters
	PetscInt    ncells[MaxNumMeshSegs  ]; // number of cells for each segment
	PetscScalar biases[MaxNumMeshSegs  ]; // biases for each segment
} MeshSegInp;
//-----------------------------------------------------------------------------
/* variables that need to be stored in the breakpoint file:
 * HorizontalFreeSurfaceHeight
 */
//-----------------------------------------------------------------------------
// marker initialization type enumeration
typedef enum
{
	PARALLEL,    // read coordinates, phase and temperature from files in parallel
	REDUNDANT,   // read phase and temperature from file redundantly (uniform coordinates)
	DIAPIR,      // diapir setup
	BLOCK,       // falling block
	SUBDUCTION,  // subduction setup with air
	FOLDING,     // multilayer folding setup (Zagros)
	DETACHMENT,  // 1-layer over detachment (Grasemann & Schmalholz 2012)
	SLAB,        // slab detachment (Thieulot et al. 2014)
	SPHERES,     // multiple falling spheres
	BANDS,       // shear band formation 3D
	RESTART      // restart of simulation
	// ... add more
} SetupType;
//-----------------------------------------------------------------------------
typedef struct {

	// mesh segments
	MeshSegInp       mseg_x;
	MeshSegInp       mseg_y;
	MeshSegInp       mseg_z;

	// marker initialization type
	SetupType        msetup;

	// domain info
	PetscScalar      W, L, H;
	PetscScalar      x_left, y_front, z_bot;
	PetscInt         nel_x, nel_y, nel_z;
	PetscInt         NumPartX, NumPartY, NumPartZ;
	PetscScalar      Setup_Diapir_Hi; //migrated from ModelSetup - for 'DIAPIR' setup
	PetscInt         nnode_x, nnode_y, nnode_z; // NOT NECESSARY
	PetscInt         cpu_x, cpu_y, cpu_z; // NOT NECESSARY - only used in FDSTAGCreate - not even initialized
	//PetscInt         finest_nnode_x, finest_nnode_y, finest_nnode_z, finest_nelx, finest_nely, finest_nelz;
	//PetscInt         remesh;
	//PetscInt         refinex, refiney, refinez;
	//PetscScalar      ampl2D,ampl3D,amplNoise,mumax, Hinterface, amp; // perturbations to grid

	// material properties
	PetscInt         num_phases;
	Material         PhaseMaterialProperties;
	PhProps          PhaseProperties;
	//PetscInt         num_particle_local;
	//PetscInt         baselevelx0, baselevely0, baselevelx1, baselevely1;

	// boundary conditions
	SBC              BC;
	//PetscInt         internalBC_frontel, internalBC_backel, internalBC_node, zdepth_BC_el, zdepth_BC_node, internalBC, internalBC_coord;
	//PetscScalar      Vx_Front, Vy_Front, Vy_Back, Vx_Back, Vy_Partx, Vy_Partz;
	//PetscScalar      MaximumSurfaceAngle;

	// time-stepping
	PetscInt         save_timesteps;
	PetscScalar      CFL;
	PetscInt         time_end;
	PetscScalar      dt_max;
	PetscScalar      dt;

	//PetscScalar      time;
	//PetscInt         itime, time_start;
	//PetscInt         time_end_temp;
	//PetscInt         EulerianAfterTimestep, temp_initialize;

	// temperature - not active
	PetscScalar      Temp_bottom, Temp_top;
	PetscScalar      GasConstant;
	//PetscScalar      Xi;
	//PetscScalar      CriticalDiagonalRatio, dt_temp;

	// optimization
	PetscInt         mpi_group_id; //migrated from OptimiseParams
	PetscScalar      LowerViscosityCutoff, UpperViscosityCutoff; // JacRes
	//PetscInt         MaxNonlinearIterations;
	//PetscScalar      NonlinearIterationsAccuracy;
	//PetscInt         StokesSolver;   // 1 - Powell-Hesteness iterations; 2 - Schur Complement Reduction; 3 - Fully Coupled Solver; 4 - MatVec Test;
	//PetscInt         VelocitySolver; // 0 - User defined; 1 - Direct (MUMPS); 2 - Galerkin geometric multigrid; 3 - Fieldsplit + Algebraic Multigrid (ML)
	//PetscBool        VelocityTest;   // Request to perform single velocity solve for test purposes
	//PetscBool        ScaleSystem;    // Request to scale linear system before solution
	PetscBool        use_fdstag_canonical; // request native staggered grid discretization ?? is it necessary

	// restart
	PetscInt         save_breakpoints, break_point_number;
	PetscInt         restart;
	//PetscInt         incr_breakpoints, fileno;

	//markers
	char             ParticleFilename[PETSC_MAX_PATH_LEN];
	char             LoadInitialParticlesDirectory[PETSC_MAX_PATH_LEN];
	char             SaveInitialParticlesDirectory[PETSC_MAX_PATH_LEN];
	PetscInt         SaveParticles;
	PetscInt         ParticleInput;

	// input/output
	char             OutputFile[PETSC_MAX_PATH_LEN];
	char             ParamFile[PETSC_MAX_PATH_LEN];
	PetscBool        InputParamFile;

	// flags
	PetscBool        SkipStokesSolver;
	PetscBool        SavePartitioning;
	//PetscInt         PlasticityCutoff;
	//PetscBool        ArtTemp;
	//PetscInt         PlasticityModel;
	//PetscInt         InitialMantleLevel;
	//PetscInt         GridAdvectionMethod, NumSurfaceNodes, num_subdt, num_phase_transitions;
	//PetscInt         LoadInitialParticlesFromDisc;
	//PetscInt         MuMeanMethod, ApplyErosion;
	//char             InitialMeshFileName[PETSC_MAX_PATH_LEN];
	//PetscBool        AnalyticalBenchmark;
	//PetscInt         *NodesDistributionCPUsFineGrid;
	//PetscInt         MaxNumLocalParticles;
	//PetscInt         NumParticlesToStartInjection;
	//PetscInt         ParticleInjectionPhase;
	//PetscInt         NonlinearIterations;
	//PetscInt         MatlabOutputFiles, VTKOutputFiles, AVDPhaseViewer;
	//PetscInt		   InitialErosionSurfaceFromFile,InitialMeshFromFile;

	// gravity
	gravityParams    GravityField;
	PetscScalar      Gravity;
	PetscScalar      GravityAngle;

	// free surface
	DM               DA_SurfaceTopography;
	Vec              SurfaceTopography;
	Vec              SurfaceTopography_Vx, SurfaceTopography_Vy, SurfaceTopography_Vz;
	Vec              BottomTopography;
	//PetscScalar      FSSA;
	//PetscScalar      FactorSurfaceLayer;

	// pushing
	PetscInt         AddPushing;
	PushParams       Pushing;

	// scaling
	PetscInt         DimensionalUnits;
	nonDimUnits      Characteristic;

	// solution vectors (part of fdstag canonical implementation, this will be abandoned)
	Vec              sol, sol_advect;
	Vec              Pressure;
	Vec              ViscosityScaling;

	// other - not bothered to sort
//PetscScalar 				SurfaceAngle, SurfaceNoiseAmplitude;
//PetscScalar 				fluvial_erosion, diffusion_erosion;
//DM						DA_Materials, DA_Processors, DA_Quadrature;		// DA that contains material properties, and that is used to distribute grid
//DM 						DA_Vel, DA_Pres;								// DA's that store Pressure-Velocity Stokes solutions
//DM 						DA_Temp;										// DA for temperature solution
//DM 						DA_BottomTopography;		// DA that stores the surface topography
//Mat 						VV_MAT,  VP_MAT;								// Matrices that contain P-V coefficients for Stokes
//Mat 						PV_MAT,  PP_MAT;								// Matrices that contain P-V coefficients for Stokes
//Mat	 					approx_S;										// Preconditioning matrix for Stokes
//Mat						TEMP_MAT;										// Matrix that contains coefficients for Temperature
//Vec						Materials;				                        // Contains material properties


} UserCtx;
//-----------------------------------------------------------------------------
#endif
