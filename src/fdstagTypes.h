/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   fdstagTypes.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

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
	PARALLEL,   // read coordinates, phase and temperature from files in parallel
	REDUNDANT,  // read phase and temperature from file redundantly (uniform coordinates)
	POLYGONS,   // read polygons from file redundantly
	DIAPIR,     // diapir setup
	HOMO,       // homogebeous model
	BLOCK,      // falling block
	SUBDUCTION, // subduction setup with air
	FOLDING,    // multilayer folding setup (Zagros)
	DETACHMENT, // 1-layer over detachment (Grasemann & Schmalholz 2012)
	SLAB,       // slab detachment (Thieulot et al. 2014)
	SPHERES,    // multiple falling spheres
	BANDS,      // shear band formation 3D
	DOMES,      // salt domes 2D
	ROTATION,   // rotation benchmark 2D
	RESTART     // restart of simulation
	// ... add more
} SetupType;

//-----------------------------------------------------------------------------
// model parameter type enumeration
#define _max_num_ModParam_type_ 30
#define _max_num_MatParam_type_ 30
typedef enum // List of model parameter types (30)
{
	// -- material model parameter types --
	_RHO0_, _RHON_, _RHOC_,                             // density
	_ETA_, _BD_, _ED_, _VD_,                            // Newtonian linear diffusion creep
	_ETA0_,	_E0_, _BN_, _N_, _EN_, _VN_,                // power-law (dislocation) creep
	_BP_, _TAUP_, _GAMMA_, _Q_, _EP_, _VP_,             // Peierls creep
	_SHEAR_, _BULK_, _KP_,                              // elasticity
	_COHESION_, _FRICTION_, _CHSOFTID_, _FRSOFTID_,     // plasticity (Drucker-Prager)
	_ALPHA_, _CP_, _K_, _A_                             // energy
	// -- others --
	// ... (geometry, pushing box , etc. ...)

} PTypes;

/*
const char *PTypesName[] ={
		// -- material model parameter types --
		"rho0","rho_n","rho_c",                        // density
		"eta","Bd","Ed","Vd",                          // newtonian linear diffiusion
		"eta0","e0","Bn","n","En","Vn",                // power-law (dislocation) creep
		"Bp","taup","gamma","q","Ep","Vp",             // Peierls creep
		"shear","bulk","Kp",                           // elasticity
		"cohesion","friction","chSoftID","frSoftID",   // plasticity (Drucker Prager)
		"alpha","cp","k","A"                           // energy

		// -- others --
		// ... (geometry, pushing box , etc. ...)
};
*/
//-----------------------------------------------------------------------------
// Structure that holds inversion parameters
typedef struct
{
	PetscInt         use;  // use inersion parameters to redefine model parameters
	PetscInt         mdN;  // number of model parameters
	PetscInt         mID;  // current model number
	PetscInt        *phs;  // model phase number
	PetscInt        *typ;  // model parameter type 
	PetscScalar     *val;  // model value
	PetscScalar     *grd;  // gradient value
	PetscScalar      mfit; // misfit value for current model parameters
} ModParam;
//-----------------------------------------------------------------------------


// Source types
typedef enum
{
	POINT,
	PLANE,
	// ... add more
} SourceType;


// Structure that holds source parameters
typedef struct // Improve and put more options
{
	SourceType  source_type;
	PetscScalar x;
	PetscScalar y;
	PetscScalar z;
	PetscInt i;
	PetscInt j;
	PetscInt k;

} SourceParam;

// Structure that holds seismic stations
typedef struct
{
	PetscScalar x;
	PetscScalar y;
	PetscScalar z;
	PetscInt i;
	PetscInt j;
	PetscInt k;
	FILE 	*output_file;
	// ...
} Station;




//-----------------------------------------------------------------------------
// Structure that holds user input data
typedef struct
{
	// mesh segments
	MeshSegInp       mseg_x;
	MeshSegInp       mseg_y;
	MeshSegInp       mseg_z;

	// marker initialization type
	SetupType        msetup;
	PetscBool		 SeismicSource;

	// domain info
	PetscScalar      W, L, H;
	PetscScalar      x_left, y_front, z_bot;
	PetscInt         nel_x, nel_y, nel_z;
	PetscInt         NumPartX, NumPartY, NumPartZ;
	PetscScalar      Setup_Diapir_Hi; // for 'DIAPIR' setup
	PetscInt         nnode_x, nnode_y, nnode_z; // NOT NECESSARY if -nel is specified

	//PetscScalar      ampl2D,ampl3D,amplNoise,mumax, Hinterface, amp; // perturbations to grid

	//PetscInt         num_particle_local;
	//PetscInt         baselevelx0, baselevely0, baselevelx1, baselevely1;

	// boundary conditions
	SBC              BC;
	//PetscInt         internalBC_frontel, internalBC_backel, internalBC_node, zdepth_BC_el, zdepth_BC_node, internalBC, internalBC_coord;
	//PetscScalar      Vx_Front, Vy_Front, Vy_Back, Vx_Back, Vy_Partx, Vy_Partz;

	// time-stepping
	PetscInt         save_timesteps;
	PetscScalar      CFL;
	PetscInt         time_end;
	PetscScalar      dt_max;
	PetscScalar      dt;

	// temperature - not active
	PetscScalar      Temp_bottom, Temp_top;
	PetscScalar      GasConstant;
	char             TemperatureFilename[MAX_PATH_LEN];

	// optimization
	PetscInt         mpi_group_id; //migrated from OptimiseParams
	PetscScalar      LowerViscosityCutoff, UpperViscosityCutoff, InitViscosity, PlastViscosity; // JacRes

	// initial guess
	PetscScalar      DII_ref;

	//PetscBool        VelocityTest;   // Request to perform single velocity solve for test purposes
	//PetscBool        ScaleSystem;    // Request to scale linear system before solution

	// restart
	PetscInt         save_breakpoints, break_point_number;
	PetscInt         restart;

	//markers
	char             ParticleFilename[MAX_PATH_LEN];
	char             LoadInitialParticlesDirectory[MAX_PATH_LEN];
	char             SaveInitialParticlesDirectory[MAX_PATH_LEN];
	PetscInt         SaveParticles;
	PetscInt         ParticleInput; // this needs to be connected in relation to marker setups

	// input/output
	char             OutputFile[MAX_PATH_LEN];
	PetscInt         PolyInVolSkip[30];
	// flags
	PetscBool        SkipStokesSolver;
	PetscBool        SavePartitioning;


	PetscBool		 ExplicitSolver; //  True => for the moment, wave propagation


	// gravity
	gravityParams    GravityField;
	PetscScalar      Gravity;
	PetscScalar      GravityAngle;

	PetscScalar      FSSA;

	// pushing
	PetscInt         AddPushing;
	PushParams       Pushing;
	
	// topography
	char             TopoFilename[MAX_PATH_LEN];
	
	// source
	SourceParam SourceParams;

	// Seismic station coordinates (in meters)
	Station Station;



} UserCtx;
//-----------------------------------------------------------------------------

#endif
