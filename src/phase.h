/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#ifndef __phase_h__
#define __phase_h__
//---------------------------------------------------------------------------

struct Scaling;
struct FB;
struct FDSTAG;
struct JacRes;
struct ModParam;

//---------------------------------------------------------------------------
//.....................   Rheology experiment type  .........................
//---------------------------------------------------------------------------

enum ExpType
{
	_UniAxial_,     // Uni-axial experiment
	_SimpleShear_,  // Simple shear experiment
	_None_
};

//---------------------------------------------------------------------------
//.......................   Softening Law Parameters  .......................
//---------------------------------------------------------------------------


struct Soft_t
{
public:

	PetscInt    ID;   // softening law ID
	PetscScalar APS1; // begin of softening APS
	PetscScalar APS2; // end of softening APS
	PetscScalar APSheal2; //APS when healTau2 is activated
	PetscScalar A;    // reduction ratio
	PetscScalar Lm;   // material length scale
  	PetscScalar healTau;   // material healing timescale parameter [Myr]  NEW FOR HEALING IN SOFTENING
	PetscScalar healTau2; // timescale parameter after APS2
};

//---------------------------------------------------------------------------
//.......................   Phase Transition Law Parameters  ................
//---------------------------------------------------------------------------

enum type
{
	_Constant_,
	_Clapeyron_,
	_Box_,
	_NotInAirBox_
};

enum Parameter
{
	_T_,
	_Pressure_,
	_Depth_,
	_X_,
	_Y_,
	_PlasticStrain_,
	_MeltFraction_,
	_Time_
};

// Structure that contains info for phase transitions
struct Ph_trans_t
{
public:

	PetscInt    ID ;				                // Phase Transition ID
	type        Type ; 					            // Type Constant or Clapeyron
	Parameter   Parameter_transition; 	            // Parameter in Constant
	char        Name_clapeyron[_str_len_] ;         // Type [Constant or Clapeyron]
	PetscInt    PhaseDirection;                     // Direction in which PT goes [0-both; 1-below2above; 2-above2below]
	PetscScalar ConstantValue ;                     // Value (if Constant)
	PetscInt    Reset;								// Rset parameters

	// Clapeyron slope
	PetscInt    neq ;                               // number of equation
	PetscScalar P0_clapeyron[_max_num_eq_] ;        // for clapeyron
	PetscScalar T0_clapeyron[_max_num_eq_] ;
	PetscScalar clapeyron_slope[_max_num_eq_] ;

	// Box-like condition
	PetscScalar     bounds[6];                      //  left, right etc. of box
	PetscInt        TempType;                       //  Temp condition [0=none, 1=constant; 2=linear; 3=halfspace]
	PetscInt 		  BoxVicinity;					  //  0-check all particles; 1-only apply PT to particles in the vicinity of the box (*2 of bounds)

	PetscInt        number_phases;
	PetscInt        PhaseBelow[_max_tr_];
	PetscInt        PhaseAbove[_max_tr_];
	PetscInt        PhaseInside[_max_tr_];
	PetscInt        PhaseOutside[_max_tr_];
	PetscScalar     dT_within;
	PetscScalar     DensityAbove[_max_tr_];
	PetscScalar     DensityBelow[_max_tr_];

	PetscScalar     topTemp;
	PetscScalar     botTemp;
	PetscScalar     cstTemp;
	PetscScalar     thermalAge;

	//Segmented NotInAirBox
	PetscInt        nsegs;                    // number of segments
	PetscScalar     xbounds[2*( _max_NotInAir_segs_ +1)]; // number of bounds in x
	PetscScalar     ybounds[2*( _max_NotInAir_segs_ +1)]; // number of bounds in y
	PetscScalar     zbounds[2*( _max_NotInAir_segs_ +1)]; // number of bounds in z
	PetscScalar    *celly_xboundL;  //left boundary of segment evaluated at ycoord of cell
	PetscScalar    *celly_xboundR;  //right boundary of segment evaluated at ycoord of cell
	PetscScalar    *cbuffL;    // memory buffer for celly_xboundL
	PetscScalar    *cbuffR;    // memory buffer for celly_xboundR

	// for moving NotInAirBox
	PetscScalar     t0_box;
	PetscScalar     t1_box;
	PetscScalar     v_box;

	// for linking NotInAirBoxes
	PetscInt      phtr_link_left;
	PetscInt      phtr_link_right;
  
};

//---------------------------------------------------------------------------
//......................   Material parameter table   .......................
//---------------------------------------------------------------------------

struct Material_t
{
public:

	PetscInt     ID;                // material ID
	PetscInt     visID;             // visualization ID
    char         Name[_str_len_];   // name (description) of the phase

	// density parameters
	PetscScalar  rho;               // reference density                          [kg/m^3]
	PetscScalar  rho_n;             // depth-dependent density model parameter    [ ]
	PetscScalar  rho_c;             // depth-dependent density model parameter    [1/m]
	PetscScalar  beta;              // pressure-dependent density model parameter [1/Pa]
	// elasticity parameters
	PetscScalar  Kb;                // bulk modulus                               [Pa]
	PetscScalar  Kp;                // pressure dependence parameter              [ ]
	PetscScalar  G;                 // shear modulus                              [Pa]
	// diffusion creep parameters
	PetscScalar  Bd;                // pre-exponential constant                   [1/Pa/s]
	PetscScalar  Ed;                // activation energy                          [J/mol]
	PetscScalar  Vd;                // activation volume                          [m^3/mol]
	// dislocation creep parameters
	PetscScalar  Bn;                // pre-exponential constant                   [1/Pa^n/s]
	PetscScalar  n;                 // power law exponent                         [ ]
	PetscScalar  En;                // activation energy                          [J/mol]
	PetscScalar  Vn;                // activation volume                          [m^3/mol]
	// Peierls creep parameters
	PetscScalar  Bp;                // pre-exponential constant                   [1/s]
	PetscScalar  Ep;                // activation energy                          [J/mol]
	PetscScalar  Vp;                // activation volume                          [m^3/mol]
	PetscScalar  taup;              // scaling stress                             [Pa]
	PetscScalar  gamma;             // approximation parameter                    [ ]
	PetscScalar  q;                 // stress-dependence parameter                [ ]
	// Frank-Kamenetzky parameters
    PetscScalar  gamma_fk;          // parameter in Frank-Kamenetzky approximation [1/K]
	PetscScalar  TRef_fk;           // Frank-Kamenetzky reference Temperature [K]
	PetscScalar  eta_fk;            // reference viscosity for Frank-Kamenetzky [Pas]
	// plasticity parameters
	PetscScalar  fr;                // friction angle                             [deg]
	PetscScalar  ch;                // cohesion
	PetscScalar  eta_vp;            // visco-plastic regularization viscosity
	PetscScalar  rp;                // ratio of pore pressure to overburden stress
	PetscInt     frSoftID;          // friction softening law ID (-1 if not defined)
	PetscInt     chSoftID;          // cohesion softening law ID (-1 if not defined)
	PetscInt     healID;            // healing ID (-1 if not defined)   
	// thermal parameters
	PetscScalar  alpha;             // thermal expansivity                        [1/K]
	PetscScalar  Cp;                // cpecific heat (capacity)                   [J/kg/K]
	PetscScalar  k;                 // thermal conductivity                       [W/m/k]
	PetscScalar  A;                 // radiogenic heat production                 [W/kg]
	PetscScalar  T;                 // optional temperature to set within the phase
	PetscScalar  nu_k;              // optional multiplication factor that is used to compute the higher conductivtiy below the conductivity boundary temperature
	PetscScalar  T_Nu;              // optional definition of conductivity boundary temperature  gito
	PetscScalar  T_liq;             // optional magma liquidus temperature for Behn & Ito [2005] dike heating model
	PetscScalar  T_sol;             // optional magma solidus temperature for Behn & Ito [2005] dike heating model
	PetscScalar  Latent_hx;          // optional magma latent heat of crystalization for Behn & Ito [2005] dike heating model
	// phase diagram
	char         pdn[_pd_name_sz_]; // Unique phase diagram number
	char         pdf[_pd_name_sz_]; // Unique phase diagram number
	PetscInt     pdAct;             // phase diagram activity flag
	PetscScalar  mfc;               // melt fraction viscosity correction
	PetscScalar  rho_melt;          // rho melt
	PetscInt     Phase_Diagram_melt;// flag that allows only to consider the melt quantity from a phase diagram
};

//---------------------------------------------------------------------------
//............   Phase diagram data   .......................................
//---------------------------------------------------------------------------

struct PData
{
	// Stores data related to Phase Diagrams

	// Size of the phase diagram in P-T space
	PetscScalar  minT[_max_num_pd_];                      // minimum temperature of diagram
	PetscScalar  maxT[_max_num_pd_];                      // maximum temperature of diagram
	PetscScalar  dT[_max_num_pd_];                        // temperature increment
	PetscInt     nT[_max_num_pd_];                        // number of temperature points

	PetscScalar  minP[_max_num_pd_];                      // minimum pressure of diagram
	PetscScalar  maxP[_max_num_pd_];                      // maximum pressure of diagram
	PetscScalar  dP[_max_num_pd_];                        // pressure increment
	PetscInt     nP[_max_num_pd_];                        // number of pressure points
	PetscInt     numProps[_max_num_pd_];                  // number of collumns (or stored properties) in phase diagram

	char         rho_pdns[_pd_name_sz_][_max_num_pd_];    // loaded phase diagram numbers
	PetscScalar  rho_v[_max_pd_sz_][_max_num_pd_];        // Array containing the actual density data (= bulk density, including that of partial melt)
	PetscScalar  rho;

	// Melt content data
	PetscScalar  Me_v[_max_pd_sz_][_max_num_pd_];          // Array containing the actual melt content data
	PetscScalar  mf;					
	
	// Rho fluid data
	PetscScalar rho_f_v[_max_pd_sz_][_max_num_pd_];
	PetscScalar rho_f;
};

//---------------------------------------------------------------------------

struct DBMat
{
	Scaling *scal;

	// phase parameters
	PetscInt     numPhases;                // number phases
	Material_t   phases[_max_num_phases_]; // phase parameters
	PetscInt     numSoft;                  // number material softening laws
	Soft_t       matSoft[_max_num_soft_];  // material softening law parameters
	Ph_trans_t   matPhtr[_max_num_tr_];    // phase transition properties
	PetscInt     numPhtr;                  // number phase transitions
};

// read material database
PetscErrorCode DBMatCreate(DBMat *dbm, FB *fb, PetscBool PrintOutput);

// read single softening law
PetscErrorCode DBMatReadSoft(DBMat *dbm, FB *fb, PetscBool PrintOutput);

// read single material phase
PetscErrorCode DBMatReadPhase(DBMat *dbm, FB *fb, PetscBool PrintOutput);

// print single material parameter
void MatPrintScalParam(
		PetscScalar par,  const char key[],   const char label[],
		Scaling    *scal, const char title[], PetscInt   *print_title);

//---------------------------------------------------------------------------
//............ PREDEFINED RHEOLOGICAL PROFILES (from literature) ............
//---------------------------------------------------------------------------

// read profile name from file
PetscErrorCode GetProfileName(FB *fb, Scaling *scal, char name[], const char key[]);

// diffusion creep profiles
PetscErrorCode SetDiffProfile(Material_t *m, char name[]);

// dislocation creep profiles
PetscErrorCode SetDislProfile(Material_t *m, char name[]);

// Peierls creep profiles
PetscErrorCode SetPeirProfile(Material_t *m, char name[]);

// correct experimental creep prefactor to tensor units
PetscErrorCode CorrExpPreFactor(PetscScalar &B, PetscScalar  n, ExpType type, PetscInt MPa);

// correct experimental stress and strain rate parameters to tensor units
PetscErrorCode CorrExpStressStrainRate(PetscScalar &D, PetscScalar &S, ExpType type, PetscInt MPa);

//---------------------------------------------------------------------------

// read phases from command line [Note: this is now directly possible]
// PetscErrorCode MatPropSetFromCL(JacRes *jr);

// Print overview (for debugging purposes only)
PetscErrorCode PrintMatProp(Material_t *MatProp);

//---------------------------------------------------------------------------
#endif
