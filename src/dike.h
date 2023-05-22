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
#ifndef __dike_h__
#define __dike_h__
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
//.......................   Dike Parameters  .......................                                                                                                      
//---------------------------------------------------------------------------

struct Dike
{
public:
  PetscInt ID;        // dike ID
  PetscInt dyndike_start;  //starting timestep for dynamic diking if 0 then no dynamic diking
  PetscInt PhaseID, PhaseTransID, nPtr;      // associated material phase and phase transition IDs
  PetscInt istep_count, nD, j1, j2;
  PetscInt istep_nave;       //number of timesteps for time averaging
  PetscScalar Mf;        // amount of magma-accomodated extension in front of box 
  PetscScalar Mb;        // amount of magma-accommodated extension in back of box
  PetscScalar Mc;        // amount of magma-accommodated extension in center of box
  PetscScalar y_Mc;      // location in y direction of Mc, if in x-direction x_Mc needs to be given or in z-direction z_Mc
  PetscScalar x_Mc;
  PetscScalar z_Mc;
  PetscScalar Tsol;
  PetscScalar filtx; 
  PetscScalar filty;
  PetscScalar drhomagma;
  PetscScalar zmax_magma;
  Vec sxx_eff_ave;
  Vec sxx_eff_ave_hist;
};
      
struct DBPropDike
{
  PetscInt numDike;                   // number of dikes
  Dike     matDike[_max_num_dike_];   // dike properties per dike ID
};

// create the dike strutures for read-in 
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// read in dike parameters
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// compute the added RHS of the dike for the continuity equation
PetscErrorCode GetDikeContr(ConstEqCtx *ctx, PetscScalar *phRat, PetscInt &Airphase, PetscScalar &dikeRHS, PetscScalar &y_c, PetscInt J);

// compute dike heat after Behn & Ito, 2008
PetscErrorCode Dike_k_heatsource(JacRes *jr,
                                Material_t *phases,
                                PetscScalar &Tc,
                                PetscScalar *phRat,          // phase ratios in the control volume
                                PetscScalar &k,
                                PetscScalar &rho_A,
                                PetscScalar &y_c,
                                PetscInt J); 

PetscErrorCode Compute_sxx_eff(JacRes *jr, PetscInt nD);
PetscErrorCode Smooth_sxx_eff(JacRes *jr, PetscInt nD, PetscInt  j1, PetscInt j2);
PetscErrorCode Set_dike_zones(JacRes *jr, PetscInt nD, PetscInt nPtr, PetscInt  j1, PetscInt j2);
PetscErrorCode Locate_Dike_Zones(AdvCtx *actx);
PetscErrorCode DynamicDike_ReadRestart(DBPropDike *dbdike, DBMat *dbm, JacRes *jr, FB *fb, FILE *fp);
PetscErrorCode DynamicDike_WriteRestart(JacRes *jr, FILE *fp);
PetscErrorCode DynamicDike_Destroy(JacRes *jr);

//---------------------------------------------------------------------------
#endif
