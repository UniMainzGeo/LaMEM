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

//---------------------------------------------------------------------------

struct Scaling;
struct FB;

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
void MatPrintScalParam(PetscScalar par, const char key[], const char label[], Scaling *scal);

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
//PetscErrorCode MatPropSetFromLibCall(JacRes *jr, ModParam *mod);

//---------------------------------------------------------------------------
#endif
