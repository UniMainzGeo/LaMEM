
// LaMEM data structures definitions
/* $Id: LaMEM.h 2911 2012-06-08 15:40:29Z apopov $ */

#ifndef __LaMEM_Types_h__
#define __LaMEM_Types_h__
//-----------------------------------------------------------------------------
// NOTE! This header is explained far too little. More comments are necessary!
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar Vx, Vy, Vz;
} Field;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar p;
} P_array;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar P[MAX_npres];
} PressureElem;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar *P;
} PressureElemDynamic;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar    Vy_front, Vy_back, Vx_left, Vx_right, Vz_bot, Vz_top, Exx, Eyy;
	PetscInt       UpperBound, LowerBound, LeftBound, RightBound, FrontBound, BackBound, InternalBound;
	DMBoundaryType BCType_x;
	DMBoundaryType BCType_y;
	DMBoundaryType BCType_z;
} sBC;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar x,y,z,num,phase,cpu,ix,iy,iz,eta,zetha,phi,T,P,PlasticStrain,Strain, Txx,Tyy,Tzz,Txy,Txz,Tyz,E2nd,T2nd,E2nd_pl, Mu_eff, Plastic, Mu_viscous;
} Particles;
//-----------------------------------------------------------------------------
// Defines the materials @ each integration point
typedef struct {
	PetscScalar Phases[MAX_ngp_vel], 			Coord[3][MAX_ngp_vel], 			 	Pressure[MAX_ngp_vel], 			Viscosity[MAX_ngp_vel];
	PetscScalar NumParticles[MAX_ngp_vel],  	ElasticShearModule[MAX_ngp_vel],	Density[MAX_ngp_vel],			DevStress[6][MAX_ngp_vel];
	PetscScalar	DevStrainrate[6][MAX_ngp_vel], 	SecondInvariantDevStress[MAX_ngp_vel], SecondInvariantDevStrainrate[MAX_ngp_vel];
	PetscScalar Temperature[MAX_ngp_vel],		Strain[MAX_ngp_vel],			PlasticStrain[MAX_ngp_vel],			ShearHeat[MAX_ngp_vel];
	PetscScalar PreExpFactor[MAX_ngp_vel],		ActivationEnergy[MAX_ngp_vel];
	PetscScalar	PlasticViscosity[MAX_ngp_vel],	Plastic[MAX_ngp_vel],			TrueViscosity[MAX_ngp_vel],			ElementVolume[MAX_ngp_vel];
} MaterialsElement;
//-----------------------------------------------------------------------------
typedef struct {
	PetscScalar *Phases, 					*Coord[3], 			 				*Pressure, 				*Viscosity;
	PetscScalar *NumParticles,  			*ElasticShearModule,				*Density,				*DevStress[6];
	PetscScalar	*DevStrainrate[6],			*SecondInvariantDevStress, 			*SecondInvariantDevStrainrate;
	PetscScalar	*Temperature,				*Strain,							*PlasticStrain,	 		*ShearHeat;
	PetscScalar *PreExpFactor,				*ActivationEnergy;
	PetscScalar *PlasticViscosity,			*Plastic,							*TrueViscosity, 		*ElementVolume;
} MaterialsElementDynamic;
//-----------------------------------------------------------------------------
// Defines physical properties for each of the phases
typedef struct {
	PetscInt	ViscosityLaw[max_num_phases],		DensityLaw[max_num_phases],	  	  PlasticityLaw[max_num_phases];
	PetscScalar mu[max_num_phases], 	  			rho[max_num_phases], 			  n_exponent[max_num_phases];
	PetscScalar A[max_num_phases], E[max_num_phases];
	PetscScalar	ElasticShearModule[max_num_phases],	Cohesion[max_num_phases],		  FrictionAngle[max_num_phases];
	PetscScalar	T_Conductivity[max_num_phases],		HeatCapacity[max_num_phases],	  RadioactiveHeat[max_num_phases];
	PetscScalar	ThermalExpansivity[max_num_phases], FrankKamenetskii[max_num_phases], Density_T0[max_num_phases];
	PetscScalar Powerlaw_e0[max_num_phases];
	PetscScalar	CohesionAfterWeakening[max_num_phases];
	PetscScalar FrictionAngleAfterWeakening[max_num_phases];
	PetscScalar	Weakening_PlasticStrain_Begin[max_num_phases];
	PetscScalar Weakening_PlasticStrain_End[max_num_phases];
	PetscScalar	Ra[max_num_phases];
} PhaseProps;
//-----------------------------------------------------------------------------
// Define structure that holds information for particle based phase transitions
typedef struct {
	PetscInt 	TransitionType, 	TransitionBelow, 	InitialPhase, 		TransformedPhase;
	PetscScalar TransitionDepth, 	TransitionP0,   	TransitionAlpha;
} PhaseTransitionProps;
//-----------------------------------------------------------------------------
typedef struct{
	PetscInt    Model, ind_fold_bot, ind_fold_top, ind_Hi_diapir;
	PetscScalar Diapir_Hi, SingleFold_H, amp_1, amp_2, Qnum, Qana, relerr;
} ModelSetup;
//-----------------------------------------------------------------------------
// Structure that holds characteristic values for nondimensionalisation
typedef struct {
	PetscScalar Length, Time, Stress, Velocity, Temperature, Viscosity;
	PetscScalar Density, kg, Strainrate, ThermalExpansivity, km, SecYear, cmYear, Myrs, MPa, Force, Watt;
	PetscScalar T_conductivity, RadioactiveHeat, Joule, HeatCapacity, Jmol;
} NonDimUnits;
//-----------------------------------------------------------------------------
// Structure that holds things like pressure, stress, strainrates at given points
typedef struct {
	PetscScalar DeviatoricStrainRate[6], DeviatoricStress[6], Pressure, Temperature, x,y,z;
	PetscScalar SecondInvariantDeviatoricStress, SecondInvariantDeviatoricStrainrate, ShearHeat;
	PetscScalar PlasticStrainrate2ndInvariant, ElementVolume, Temperature_diff;
} PointWiseInformation;
//-----------------------------------------------------------------------------
// Structure that holds time-dependent data
typedef struct {
	PetscScalar Time, Vrms, Vx_max, Vx_min, Vy_max, Vy_min, Vz_max, Vz_min;
	PetscScalar MinXCoordPhase[max_num_phases], MaxXCoordPhase[max_num_phases];
	PetscScalar MinYCoordPhase[max_num_phases], MaxYCoordPhase[max_num_phases];
	PetscScalar MinZCoordPhase[max_num_phases], MaxZCoordPhase[max_num_phases];
	PetscScalar	DevStress[6],					DevStrainrate[6];
	PetscScalar MaxTopography, 		MinTopography, MeanTopography;
} GlobalTimeDependentData;
//-----------------------------------------------------------------------------
// Structure that holds data related to the serial Finite Difference Erosion Code ('external' package)
typedef struct {
	PetscInt	ResolutionFactorX;			// how much denser is the erosion grid compared to the numerical grid we use for the Stokes solver?
	PetscInt	ResolutionFactorY;			// how much denser is the erosion grid compared to the numerical grid we use for the Stokes solver?
	DM			DA_FE_ErosionCode;	 		// DMDA that contains the erosion surface
	Vec 		ErosionSurface;				// vector with the surface for use in the erosion code
	PetscScalar	dt;							// approximate timestep used for erosion code
	PetscScalar InitialRandomNoise_m;		// initial random noise (in m) on the eroded surface
	PetscScalar InitialUpliftedSide_m;		// by how much meters is the right side of the surface topography uplifted?
	PetscScalar	rain_m_year, rain;			// how much meter of rain fell last year [uniformly over the model]/ [and in m/s]
	PetscScalar	k0,n,c;						// erodability constants in k=k0 + c*q^n
    PetscInt    BC, fill_lake, nbre_river, mode_river;
    PetscScalar location_river[100];
    PetscScalar rain_river_year;
} sFE_ErosionCode;

//-----------------------------------------------------------------------------
// Structure that holds erosion and sedimentation parameters
typedef struct {
	PetscInt 		ErosionModel;			//	0 - no erosion [default], 1 - infinitely fast,  2 - finite difference erosion model
	PetscInt 		ApplyErosion, UseInternalFreeSurface;
	PetscScalar 	fluvial_erosion, diffusion_erosion;
	PetscScalar 	SurfaceAngle, SurfaceNoiseAmplitude;
	sFE_ErosionCode	FE_ErosionCode;			// parameters used in the finite-difference erosion code
	PetscInt		SedimentationModel; 	// 0 - none [default], 1 - constant rate
	PetscScalar 	InitialFreeSurfaceHeight,SedimentationRate_cmYr, SedimentationRate;
	PetscScalar		HorizontalFreeSurfaceHeight, SedimentLayerThicknessYears;
	PetscInt		PhaseFirstSedimentedLayer, PhaseLastSedimentedLayer, StickyAirPhase, PhaseSedimented;
	PetscInt 		baselevelx0, baselevely0, baselevelx1, baselevely1;
} ErosionParams;
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
} GravityParams;
//-----------------------------------------------------------------------------
// Structure that holds isostasy parameters
typedef struct {
	PetscInt		GetIt;
	PetscInt        SaveRef;
	PetscInt		ref_xi,ref_yi;
	PetscScalar		corr_topo;
	PetscScalar		ref_rho;
	PetscScalar     TisoStdDev;
	char			RefDatFile2load[PETSC_MAX_PATH_LEN];
} IsostasyParams;
//-----------------------------------------------------------------------------
// Structure that holds surface velocity parameters
typedef struct {
	PetscInt		GetIt;
	PetscInt		SaveRef;
	PetscScalar     VxStdDev,VyStdDev,VzStdDev;
	char			RefDatFile2load[PETSC_MAX_PATH_LEN];
} SurfVelParams;
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
} PushingParams;
//-----------------------------------------------------------------------------
// Structure that holds output parameters
typedef struct {
	PetscInt		velocity,temperature,surface_topography,bottom_topography,quadrature;

} OutputParams;
//-----------------------------------------------------------------------------
// Structure that holds optimisation parameters
typedef struct {
	PetscInt		GetIt;
	PetscScalar		MisfitSurfVel, MisfitGravity, MisfitTiso;
	PetscScalar     NSurfVel,NGrav,NTiso ;
	PetscScalar     SumAbsVel,SumSqrsVel,SumAbsGrav,SumSqrsGrav,SumAbsTiso,SumSqrsTiso;
	PetscInt		mpi_group_id;
 } OptimiseParams;
//-----------------------------------------------------------------------------
// Structure that holds all information required to build the FDSTAG stencil at all points necessary
typedef struct {
	// DA's that handle FDSTAG discretization
	DM 	DA_CENTER, DA_CORNER, DA_XY_POINTS, DA_XZ_POINTS, DA_YZ_POINTS;
	// History-dependent variables:
	Vec		Center_PhaseProportions[max_num_phases], 	Center_Temperature, 	Center_Strain, 		Center_PlasticStrain, 	Center_NumParticles;
	Vec		Corner_PhaseProportions[max_num_phases], 	Corner_Temperature,		Corner_PhaseProportions_local[max_num_phases];
	Vec		XYPoints_PhaseProportions[max_num_phases],	XYPoints_Temperature, 	XYPoints_Strain, 	XYPoints_PlasticStrain;
	Vec		XZPoints_PhaseProportions[max_num_phases],	XZPoints_Temperature, 	XZPoints_Strain, 	XZPoints_PlasticStrain;
	Vec		YZPoints_PhaseProportions[max_num_phases], 	YZPoints_Temperature, 	YZPoints_Strain, 	YZPoints_PlasticStrain;
	// Material parameters which are NOT history dependent:
	Vec 	Center_T2nd, 		Center_E2nd, 		Center_Pressure, 		Center_EffectiveViscosity, 		Center_Density;
	Vec 	XYPoints_T2nd, 		XYPoints_E2nd, 		XYPoints_Pressure, 		XYPoints_EffectiveViscosity, 	XYPoints_Density;
	Vec 	XZPoints_T2nd, 		XZPoints_E2nd, 		XZPoints_Pressure, 		XZPoints_EffectiveViscosity, 	XZPoints_Density;
	Vec 	YZPoints_T2nd, 		YZPoints_E2nd, 		YZPoints_Pressure, 		YZPoints_EffectiveViscosity, 	YZPoints_Density;
	Vec 	Corner_Pressure, 	Corner_Density, 	Corner_HeatCapacity, 	Corner_Conductivity, 			Corner_RadioactiveHeat;
} FDSTAG_MaterialProperties;
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
 *
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
	SPHERES      // multiple falling spheres
	// ... add more
} SetupType;
//-----------------------------------------------------------------------------
typedef struct {
  PetscScalar				W,L,H,x_left,y_front,z_bot,ampl2D,ampl3D,amplNoise,mumax, Hinterface, amp;
  PetscScalar				LowerViscosityCutoff, UpperViscosityCutoff;
  PetscInt					StokesSolver;   // 1 - Powell-Hesteness iterations; 2 - Schur Complement Reduction; 3 - Fully Coupled Solver; 4 - MatVec Test;
  PetscInt					VelocitySolver; // 0 - User defined; 1 - Direct (MUMPS); 2 - Galerkin geometric multigrid; 3 - Fieldsplit + Algebraic Multigrid (ML)
  PetscBool					VelocityTest;   // Request to perform single velocity solve for test purposes
  PetscBool					ScaleSystem;    // Request to scale linear system before solution
  PetscBool                 use_fdstag_canonical; // request native staggered grid discretization
  PetscInt					nnode_x, nnode_y, nnode_z, refinex, refiney, refinez, itime, time_start, time_end, time_end_temp, save_timesteps;
  PetscInt 					cpu_x, cpu_y, cpu_z, LoadInitialParticlesFromDisc;
  PetscInt					MatlabOutputFiles, VTKOutputFiles, AVDPhaseViewer;
  PetscInt					InitialErosionSurfaceFromFile,InitialMeshFromFile, save_breakpoints, break_point_number, incr_breakpoints, fileno, MaxNonlinearIterations, internalBC_frontel, internalBC_backel, internalBC_node, zdepth_BC_el, zdepth_BC_node, internalBC, internalBC_coord;
  PetscInt					EulerianAfterTimestep, temp_initialize;
  PetscScalar				time, dt, CFL, FactorSurfaceLayer;
  PetscScalar				Gravity, Temp_bottom, Temp_top, Xi, dt_max, CriticalDiagonalRatio, GasConstant, dt_temp;
  PetscScalar           	Vx_Front, Vy_Front, Vy_Back, Vx_Back, MaximumSurfaceAngle,Vy_Partx, Vy_Partz, NonlinearIterationsAccuracy;
  PetscInt					num_particle_local, num_phases;
  PetscInt					MuMeanMethod, restart, ApplyErosion;
  PetscInt 					baselevelx0, baselevely0, baselevelx1, baselevely1, PlasticityCutoff;
  PetscScalar 				SurfaceAngle, SurfaceNoiseAmplitude;
  PetscScalar 				fluvial_erosion, diffusion_erosion, GravityAngle, FSSA;
  DM						DA_Materials, DA_Processors, DA_Quadrature;		// DA that contains material properties, and that is used to distribute grid
  DM 						DA_Vel, DA_Pres;								// DA's that store Pressure-Velocity Stokes solutions
  DM 						DA_Temp;										// DA for temperature solution
  DM 						DA_SurfaceTopography, DA_BottomTopography;		// DA that stores the surface topography
  Mat 						VV_MAT,  VP_MAT;								// Matrices that contain P-V coefficients for Stokes
  Mat 						PV_MAT,  PP_MAT;								// Matrices that contain P-V coefficients for Stokes
  Mat	 					approx_S;										// Preconditioning matrix for Stokes
  Mat						TEMP_MAT;										// Matrix that contains coefficients for Temperature
  Vec						Materials;				                        // Contains material properties
  FDSTAG_MaterialProperties FDSTAG;                                         // FDSTAG-specific data
  Vec						SurfaceTopography, BottomTopography, SurfaceTopography_Vx, SurfaceTopography_Vy, SurfaceTopography_Vz;
  ModelSetup				Setup;
  PhaseProps    			PhaseProperties;
  sBC						BC;
  char						OutputFile[PETSC_MAX_PATH_LEN], ParamFile[PETSC_MAX_PATH_LEN], ParticleFilename[PETSC_MAX_PATH_LEN];
  char						InitialMeshFileName[PETSC_MAX_PATH_LEN], LoadInitialParticlesDirectory[PETSC_MAX_PATH_LEN], SaveInitialParticlesDirectory[PETSC_MAX_PATH_LEN];
  PetscBool				    InputParamFile, AnalyticalBenchmark,SkipStokesSolver, SavePartitioning;
  Particles					*ParticlesLocal;
  PetscInt      			*NodesDistributionCPUsFineGrid, MaxNumLocalParticles, NumParticlesToStartInjection, ParticleInjectionPhase;
  PetscInt					ParticleInput, SaveParticles, DimensionalUnits;
  PetscInt					PlasticityModel, InitialMantleLevel;
  PetscInt					GridAdvectionMethod, NumSurfaceNodes, num_subdt, num_phase_transitions;
  PetscInt					NumPartX, NumPartY, NumPartZ, NonlinearIterations;
  PetscInt 					nel_x, nel_y, nel_z, finest_nnode_x, finest_nnode_y, finest_nnode_z, finest_nelx, finest_nely, finest_nelz;
  PetscInt					remesh;
  NonDimUnits				Characteristic;
  PhaseTransitionProps		PhaseTransitions[100];
  GlobalTimeDependentData	*TimeDependentData;
  ErosionParams 			ErosionParameters;
  Material 					PhaseMaterialProperties;
  OptimiseParams			Optimisation;
  SurfVelParams				SurfVelField;
  GravityParams				GravityField;
  PetscInt					AddPushing;
  PushingParams				Pushing;
  OutputParams				Output;
  IsostasyParams            Isostasy;
  PetscBool                 ArtTemp;
  // mesh segments
  MeshSegInp                mseg_x;
  MeshSegInp                mseg_y;
  MeshSegInp                mseg_z;
  // marker initialization type
  SetupType                 msetup;
  // solution vectors (part of fdstag canonical implementation, this will be abandoned)
  Vec                       sol, sol_advect;
  Vec                       Pressure;
  Vec                       ViscosityScaling;
} UserContext;
//-----------------------------------------------------------------------------
#endif
