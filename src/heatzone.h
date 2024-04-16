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
  PetscScalar asthenoTemp;   // required: asthenospheric temperature
  PetscScalar rho;           // density of heating material
  PetscScalar Cp;            // cpecific heat of heating material
  PetscScalar timeStart;     // optional parameter to delay heat zone activation
  PetscScalar tempStart;     // optional parameter to delay heat zone activation
  PetscScalar heatRate;      // required for q_hotspot: heating rate of the hotspot
  PetscScalar spreadingRate; // optional parameter for q_ridge that indicates the spreading velocity of the plate; if not defined it uses bvel_velin specified
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

// compensate for dike heat contributions if within heat zone
PetscErrorCode SubtractDikeHeatSource(JacRes *jr,
                                      Material_t *phases,
                                      PetscScalar &Tc,
                                      PetscScalar *phRat,
                                      PetscScalar &rho_A,
                                      PetscScalar &y_c,
                                      PetscInt J,
                                      PetscScalar sxx_eff_ave_cell);

//---------------------------------------------------------------------------
#endif
