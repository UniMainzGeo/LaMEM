//---------------------------------------------------------------------------
//..........................   FDSTAG DATA TYPES  ...........................
//---------------------------------------------------------------------------
#ifndef __fdstagTypes_h__
#define __fdstagTypes_h__
//-----------------------------------------------------------------------------
// Structure that holds boundary conditions info
typedef struct
{
	PetscScalar    Vy_front, Vy_back, Vx_left, Vx_right, Vz_bot, Vz_top, Exx, Eyy;
	PetscInt       UpperBound, LowerBound, LeftBound, RightBound, FrontBound, BackBound, InternalBound;
	DMBoundaryType BCType_x;
	DMBoundaryType BCType_y;
	DMBoundaryType BCType_z;
} SBC;
//-----------------------------------------------------------------------------
// Structure that holds gravity parameters - not yet used
typedef struct
{
	PetscInt     GetIt;
	PetscInt     SaveDebug,SaveVTK,SaveRef;
	PetscBool    UseNumerics, UseAnalytics;
	PetscInt     survey_nx, survey_ny;
	PetscScalar  survey_xs, survey_xm;
	PetscScalar  survey_ys, survey_ym;
	PetscScalar  survey_z ;
	PetscScalar  ReferenceDensity;
	PetscScalar  StdDev;
	PetscScalar  LithColDens[9],LithColDepth[8];
	PetscInt     num_intp,LithColNum;
	char         RefDatFile2load[MAX_PATH_LEN];
} gravityParams;
//-----------------------------------------------------------------------------
// Structure that holds pushing parameters
typedef struct
{
	PetscScalar  L_block, W_block, H_block;
	PetscScalar  x_center_block, y_center_block, z_center_block;
	PetscScalar  V_push[10], omega[10], time[11];
	PetscScalar  theta;
	PetscInt     num_changes, reset_pushing_coord, ind_change;
	PetscInt     coord_advect[10], dir[10];
} PushParams;
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
// marker initialization type enumeration
typedef enum
{
	PARALLEL,    // read coordinates, phase and temperature from files in parallel
	REDUNDANT,   // read phase and temperature from file redundantly (uniform coordinates)
	POLYGONS,    // read polygons from file redundantly
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
// Structure that holds user input data
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
	PetscScalar      Setup_Diapir_Hi; // for 'DIAPIR' setup
	PetscInt         nnode_x, nnode_y, nnode_z; // NOT NECESSARY if -nel is specified
	//PetscInt         cpu_x, cpu_y, cpu_z; // not necessary - only used in FDSTAGCreate - not even initialized
	//PetscInt         finest_nnode_x, finest_nnode_y, finest_nnode_z, finest_nelx, finest_nely, finest_nelz;
	//PetscInt         remesh;
	//PetscInt         refinex, refiney, refinez;
	//PetscScalar      ampl2D,ampl3D,amplNoise,mumax, Hinterface, amp; // perturbations to grid

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
	PetscScalar      LowerViscosityCutoff, UpperViscosityCutoff, InitViscosity, PlastViscosity; // JacRes

	// initial guess
	PetscScalar      DII_ref;

	//PetscInt         MaxNonlinearIterations;
	//PetscScalar      NonlinearIterationsAccuracy;
	//PetscInt         StokesSolver;   // 1 - Powell-Hesteness iterations; 2 - Schur Complement Reduction; 3 - Fully Coupled Solver; 4 - MatVec Test;
	//PetscInt         VelocitySolver; // 0 - User defined; 1 - Direct (MUMPS); 2 - Galerkin geometric multigrid; 3 - Fieldsplit + Algebraic Multigrid (ML)
	//PetscBool        VelocityTest;   // Request to perform single velocity solve for test purposes
	//PetscBool        ScaleSystem;    // Request to scale linear system before solution
	PetscBool        use_fdstag_canonical; // request native staggered grid discretization

	// restart
	PetscInt         save_breakpoints, break_point_number;
	PetscInt         restart;
	//PetscInt         incr_breakpoints, fileno;

	//markers
	char             ParticleFilename[MAX_PATH_LEN];
	char             LoadInitialParticlesDirectory[MAX_PATH_LEN];
	char             SaveInitialParticlesDirectory[MAX_PATH_LEN];
	PetscInt         SaveParticles;
	PetscInt         ParticleInput; // this needs to be connected in relation to marker setups

	// input/output
	char             OutputFile[MAX_PATH_LEN];
//	char             ParamFile[MAX_PATH_LEN];
//	PetscBool        InputParamFile;

	// flags
	PetscBool        SkipStokesSolver;
	PetscBool        SavePartitioning;
	//PetscInt         PlasticityCutoff;
	//PetscBool        ArtTemp;
	//PetscInt         PlasticityModel;
	//PetscInt         InitialMantleLevel;
	//PetscInt         GridAdvectionMethod, NumSurfaceNodes, num_subdt, num_phase_transitions;
	//PetscInt         LoadInitialParticlesFromDisc;
	//PetscInt         MuMeanMethod;
	//char             InitialMeshFileName[MAX_PATH_LEN];
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
	PetscScalar      FSSA;
	//PetscScalar      FactorSurfaceLayer;

	// pushing
	PetscInt         AddPushing;
	PushParams       Pushing;

	// solution vectors (part of fdstag canonical implementation, this will be abandoned)
	Vec              sol, sol_advect;
	Vec              Pressure;
	Vec              ViscosityScaling;

	// other - erosion
	//PetscInt         ApplyErosion;
	//PetscScalar      SurfaceAngle, SurfaceNoiseAmplitude;
	//PetscScalar      fluvial_erosion, diffusion_erosion;

} UserCtx;
//-----------------------------------------------------------------------------
#endif
