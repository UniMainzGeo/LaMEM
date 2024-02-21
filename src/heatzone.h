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
#ifndef __heatzone_h__
#define __heatzone_h__
//---------------------------------------------------------------------------

struct FB;
struct ConstEqCtx;
struct DBMat;
struct FDSTAG;
struct TSSol;
struct JacRes;
struct Controls;
struct AdvCtx;

//---------------------------------------------------------------------------
//.......................   HeatZone Parameters  .......................
//---------------------------------------------------------------------------

struct HeatZone
{
public:
  PetscInt ID;               // heatzone ID
  PetscInt HeatFunction;     // heating function to use [0=q_hotspot, 1=q_ridge] from Mittelstaedt et. al., 2008
  PetscInt FunctType;        // heat zone dimensionality [1d_x-gauss, 2d_xy-gauss, 3d_xyz-gauss], gaussian dependent on width in x-direction and heat zone center coordinate
  PetscScalar bounds[6];     //  left, right, front, back, bottom, top of heatzone
  PetscScalar rho;           // density of heating material
  PetscScalar Cp;            // cpecific heat of heating material
  PetscScalar asthenoTemp;   // required: asthenospheric temperature
  PetscScalar heatRate;      // required for q_hotspot: heating rate of the hotspot
  PetscScalar spreadingRate; // optional parameter for q_ridge that indicates the spreading velocity of the plate; if not defined it uses bvel_velin specified

/*   
  PetscInt PhaseID, PhaseTransID, nPtr; // associated material phase and phase transition IDs

  PetscInt istep_count, nD, j1, j2;

  PetscInt dynheatzone_start;  //starting timestep for dynamic diking if 0 then no dynamic diking
  PetscScalar A; // Smoothing parameter for variable M calculation
  PetscScalar B; // Value to prevent NaNs
  PetscScalar knee; // Determines the transition from min to max M in the M_val equation
	PetscScalar Ts; // Tensile strength of rock for variable M calculation (Pa)
  PetscScalar zeta_0; // Initial bulk viscosity for variable M calculation (Pa*s) *revisit [local initial bulk viscosity]
  PetscScalar Mf;        // amount of magma-accomodated extension in front of box 
  PetscScalar Mb;        // amount of magma-accommodated extension in back of box
  PetscScalar Mc;        // amount of magma-accommodated extension in center of box
  PetscInt istep_nave;      // number of timesteps for time averaging
  PetscInt nstep_locate;    // Locate heatzone every nstep_locate timestep to allow elastic stresses to settle down between relocations
  PetscInt out_stress;      // option to output mean stresses to std out
  PetscInt out_heatzoneloc; // option to output heatzone location to std out
  PetscScalar y_Mc;         // location in y direction of Mc, if in x-direction x_Mc needs to be given or in z-direction z_Mc
  PetscScalar x_Mc;
  PetscScalar z_Mc;
  PetscScalar Tsol;
  PetscScalar filtx;
  PetscScalar filty;
  PetscScalar drhomagma;
  PetscScalar zmax_magma;
  PetscScalar magPfac;
  PetscScalar magPwidth;
  // PetscScalar ymindyn;
  // PetscScalar ymaxdyn;
  Vec sxx_eff_ave;
  Vec magPressure;
  Vec sxx_eff_ave_hist; */
};

struct DBPropHeatZone
{
  PetscInt numHeatZone;                     // number of heatzones
  HeatZone matHeatZone[_max_num_heatzone_]; // heatzone properties per heatzone ID
};

// create the heatzone strutures for read-in
PetscErrorCode DBHeatZoneCreate(DBPropHeatZone *dbheatzone, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// read in heatzone parameters
PetscErrorCode DBReadHeatZone(DBPropHeatZone *dbheatzone, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// compute heatzone heat after Mittelstaedt et. al., 2008
PetscErrorCode GetHeatZoneSource(JacRes *jr,
                                     Material_t *phases,
                                     PetscScalar &Tc,
                                     PetscScalar *phRat, // phase ratios in the control volume
                                     PetscScalar &rho_A,
                                     PetscScalar &y_c,
                                     PetscScalar &x_c,
                                     PetscScalar &z_c,
                                     PetscInt J,
                                     PetscScalar sxx_eff_ave_cell);

/* // compute and subtract dike heat for overlapping heat zone contributions
PetscErrorCode SubtractDikeHeatSource(JacRes *jr,
                                Material_t *phases,
                                PetscScalar &Tc,
                                PetscScalar *phRat,          // phase ratios in the control volume
                                PetscScalar &k,
                                PetscScalar &rho_A,
                                PetscScalar &y_c,
                                PetscInt J); */

//---------------------------------------------------------------------------
#endif
