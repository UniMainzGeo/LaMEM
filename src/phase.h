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
 **    filename:   phase.h
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
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#ifndef __phase_h__
#define __phase_h__
//---------------------------------------------------------------------------

// max number of phases
#define max_num_phases 32

// max number of soft laws
#define max_num_soft   10

// maximums for Pds
#define max_num_pd    8     // max no of phase diagrams
#define max_num_ro    40100  // max grid size of Pd
#define max_name      54     // Length of the unique face diagram name

//---------------------------------------------------------------------------

struct Scaling;
struct FB;
struct JacRes;
struct ModParam;

//---------------------------------------------------------------------------
//.......................   Softening Law Parameters  .......................
//---------------------------------------------------------------------------

struct Soft_t
{
public:

	PetscInt    ID;   // softening law ID
	PetscScalar APS1; // begin of softening APS
	PetscScalar APS2; // end of softening APS
	PetscScalar A;    // reduction ratio

};

//---------------------------------------------------------------------------
//......................   Material parameter table   .......................
//---------------------------------------------------------------------------

struct Material_t
{
public:

	PetscInt     ID;       // material ID
	// density parameters
	PetscScalar  rho;      // reference density                          [kg/m^3]
	PetscScalar  rho_n;    // depth-dependent density model parameter    [ ]
	PetscScalar  rho_c;    // depth-dependent density model parameter    [1/m]
	PetscScalar  beta;     // pressure-dependent density model parameter [1/Pa]
	// elasticity parameters
	PetscScalar  K;        // bulk modulus                               [Pa]
	PetscScalar  Kp;       // pressure dependence parameter              [ ]
	PetscScalar  G;        // shear modulus                              [Pa]
	// diffusion creep parameters
	PetscScalar  Bd;       // pre-exponential constant                   [1/Pa/s]
	PetscScalar  Ed;       // activation energy                          [J/mol]
	PetscScalar  Vd;       // activation volume                          [m^3/mol]
	// dislocation creep parameters
	PetscScalar  Bn;       // pre-exponential constant                   [1/Pa^n/s]
	PetscScalar  n;        // power law exponent                         [ ]
	PetscScalar  En;       // activation energy                          [J/mol]
	PetscScalar  Vn;       // activation volume                          [m^3/mol]
	// Peierls creep parameters
	PetscScalar  Bp;       // pre-exponential constant                   [1/s]
	PetscScalar  Ep;       // activation energy                          [J/mol]
	PetscScalar  Vp;       // activation volume                          [m^3/mol]
	PetscScalar  taup;     // scaling stress                             [Pa]
	PetscScalar  gamma;    // approximation parameter                    [ ]
	PetscScalar  q;        // stress-dependence parameter                [ ]
	// plasticity parameters
	PetscScalar  fr;       // friction angle                             [deg]
	PetscScalar  ch;       // cohesion
	PetscScalar  rp;       // ratio of pore pressure to overburden stress
	PetscInt     frSoftID; // friction softening law ID (-1 if not defined)
	PetscInt     chSoftID; // cohesion softening law ID (-1 if not defined)
	// thermal parameters
	PetscScalar  alpha;    // thermal expansivity                        [1/K]
	PetscScalar  Cp;       // cpecific heat (capacity)                   [J/kg/K]
	PetscScalar  k;        // thermal conductivity                       [W/m/k]
	PetscScalar  A;        // radiogenic heat production                 [W/kg]
	PetscScalar  T;        // optional temperature to set within the phase
	// Phase diagram
	char         pdn[max_name];   // Unique phase diagram number
	char         pdf[max_name];   // Unique phase diagram number
	PetscInt     Pd_rho;          // density from phase diagram?
	// MeltExtraction Parameter
	PetscScalar Mtrs;          // Threshold melt extraction                   []
	PetscScalar Mleft;         // Minimum amount of melt left in the source  []
	PetscScalar Mmax;          // Maximum Melt extractable from phase        []
	PetscScalar RelInt;        // Relative amount of intrusion               []
	PetscScalar TInt;          // Temperature of the intrusion               [Deg C]
	PetscScalar TExt;          // Temperature of extrusion                   [Deg C]
	PetscInt 	PhInt;  // Phase Id of the intrusion                  []
	PetscInt	PhExt;  // Phase Id of the effusion                   []
	PetscInt	PhNext;  // Phase Id of the effusion                   []
	PetscScalar DInt;          // Depth of intrusion                         [m]
	PetscScalar DExt;          // Depth of intrusion                         [m]
	PetscScalar pMant;         // Specify if a phase is mantle or not [0 or 1]
	PetscScalar S;             // Random Source
	PetscInt    MeltE;         // Control value that states if
	};

//---------------------------------------------------------------------------
//............   Phase diagram data   .......................................
//---------------------------------------------------------------------------

struct PData
{
	// Stores data related to Phase Diagrams

	// Size of the phase diagram in P-T space
	PetscScalar  minT[max_num_pd];                      // minimum temperature of diagram
	PetscScalar  maxT[max_num_pd];                      // maximum temperature of diagram
	PetscScalar  dT[max_num_pd];                        // temperature increment
	PetscInt     nT[max_num_pd];                        // number of temperature points

	PetscScalar  minP[max_num_pd];                      // minimum pressure of diagram
	PetscScalar  maxP[max_num_pd];                      // maximum pressure of diagram
	PetscScalar  dP[max_num_pd];                        // pressure increment
	PetscInt     nP[max_num_pd];                        // number of pressure points
	PetscInt     numProps[max_num_pd];                  // number of collumns (or stored properties) in phase diagram

	char         rho_pdns[max_name][max_num_pd];        // loaded phase diagram numbers
	PetscScalar  rho_v[max_num_ro][max_num_pd];         // Array containing the actual density data (= bulk density, including that of partial melt)
	PetscScalar  rho;

	// Melt content data
	PetscScalar  Me_v[max_num_ro][max_num_pd];          // Array containing the actual melt content data
	PetscScalar  mf;
	PetscScalar  mfext;
	PetscScalar  mfextot;
	
	// Rho fluid data
	PetscScalar rho_f_v[max_num_ro][max_num_pd];
	PetscScalar rho_f;
};

//---------------------------------------------------------------------------

struct DBMat
{
public:

	Scaling *scal;

	// phase parameters
	PetscInt     numPhases;              // number phases
	Material_t   phases[max_num_phases]; // phase parameters
	PetscInt     numSoft;                // number material softening laws
	Soft_t       matSoft[max_num_soft];  // material softening law parameters

};

// read material database
PetscErrorCode DBMatCreate(DBMat *dbm, FB *fb);

// read single softening law
PetscErrorCode DBMatReadSoft(DBMat *dbm, FB *fb);

// read single material phase
PetscErrorCode DBMatReadPhase(DBMat *dbm, FB *fb);

// print single material parameter
void MatPrintScalParam(
		PetscScalar par,  const char key[],   const char label[],
		Scaling    *scal, const char title[], PetscInt   *print_title);

//---------------------------------------------------------------------------
//............ PREDEFINED RHEOLOGICAL PROFILES (from literature) ............
//---------------------------------------------------------------------------
enum TensorCorrection
{
	_UniAxial_,      // Uni-axial experiment
	_SimpleShear_,   // Simple shear experiment
	_None_           // geological-scale units

};

// read profile name from file
PetscErrorCode GetProfileName(FB *fb, Scaling *scal, char name[], const char key[]);

// diffusion creep profiles
PetscErrorCode SetDiffProfile(Material_t *m, char name[]);

// dislocation creep profiles
PetscErrorCode SetDislProfile(Material_t *m, char name[]);

// Peierls creep profiles
PetscErrorCode SetPeirProfile(Material_t *m, char name[]);

// units and tensor correction
PetscErrorCode SetProfileCorrection(PetscScalar *B, PetscScalar n, TensorCorrection tensorCorrection, PetscInt MPa);

//---------------------------------------------------------------------------

// read phases from command line
// PetscErrorCode MatPropSetFromCL(JacRes *jr);

// assign phases from calling function
PetscErrorCode MatPropSetFromLibCall(JacRes *jr, ModParam *mod, FB *fb);

//---------------------------------------------------------------------------
#endif
