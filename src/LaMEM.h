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
//.........................   Main include file   ...........................
//---------------------------------------------------------------------------

#ifndef __LaMEM_h__
#define __LaMEM_h__

//-----------------------------------------------------------------------------
// EXTERNAL INCLUDES
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <petsc.h>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#ifdef _WIN32
#include "asprintf.h"       // required for some windows compilers
#endif

using namespace std;

//-----------------------------------------------------------------------------
//   PREFERABLE VARIABLES
//
//   PetscInt    - for indices                    (int or long long int)
//   PetscScalar - for floating point variables   (float or double)
//   float       - for reduced size output
//   size_t      - for variable sizes
//   uint64_t    - for offsets & counters
//   PetscMPIInt - for passing service integer parameters to MPI functions (int)
//   MPIU_SCALAR - appropriate MPI Data Type for sending/receiving PetsScalar
//   MPIU_INT    - appropriate MPI Data Type for sending/receiving PetscInt
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SIZE LIMITS
//-----------------------------------------------------------------------------

// number of neighbor domains in 3D lattice (including self)
const PetscInt _num_neighb_ = 27;

// string length (two null characters are reserved in the end, i.e. 128)
const PetscInt _str_len_ = 260;

// marker storage capacity overhead
//const PetscScalar _cap_overhead_ = 1.61803398875;

const PetscScalar _cap_overhead_ = (1.0 + PetscSqrtReal(5.0))/2.0;

// maximum marker per cell per direction
const PetscInt _max_nmark_ = 5;

// minimum marker per cell per direction
const PetscInt _min_nmark_ = 2;

// cell marker buffer size
const PetscInt _mark_buff_sz_ = 256;

// local marker buffer size as percentage of local number of markers
const PetscInt _mark_buff_ratio_ = 5;

// maximum number of strain rate application periods
const PetscInt _max_periods_ = 20;

// maximum number of time steps
const PetscInt _max_num_steps_ = 2000;

// maximum number of Bezier blocks
const PetscInt _max_boxes_ = 5;

// maximum number of Bezier path points
const PetscInt _max_path_points_ = 25;

// maximum number of polygon points of Bezier block
const PetscInt _max_poly_points_ = 50;

// maximum number of mesh segments in every direction
const PetscInt _max_num_segs_ = 10;

// maximum number of cells per mesh segment
const PetscInt _max_num_cells_ = 4096;

// maximum number of processes in every direction
const PetscInt _max_num_procs_ = 1024;

// FDSTAG near null space size
const PetscInt _max_nullsp_sz_ = 4;

// maximum number of geometry primitives
const PetscInt _max_geom_ = 100;

// maximum number of polygons on the same level
const PetscInt _max_polygons_ = 10;

// maximum number of observation types
const PetscInt _max_num_obs_ = 7;

// maximum number of components in the output vector (3D)
const PetscInt _max_num_comp_ = 9;

// maximum number of components in the output vector (surface)
const PetscInt _max_num_comp_surf_ = 3;

// maximum number of phase aggregates for output
const PetscInt _max_num_phase_agg_ = 5;

// maximum number of phases
const PetscInt _max_num_phases_ = 32;

// maximum number of softening laws
const PetscInt _max_num_soft_ = 10;

// max char length
const PetscInt char_ph_tr = 10;

// maximum number of phase transition law
const PetscInt _max_num_tr_ = 20;

// maximum number of segments of NotInAirBoxes
const PetscInt _max_NotInAir_segs_ = 6;

// maximum number of dikes
const PetscInt _max_num_dike_ = 12;

// maximum number of equation parameter
const PetscInt _max_num_eq_ = 2;

// maximum number of phase transition per each phase
const PetscInt _max_tr_  = 8;

// maximum number of phase diagrams
const PetscInt _max_num_pd_  = 8;

// maximum grid size of phase diagram
const PetscInt _max_pd_sz_ = 40100;

// length of unique phase diagram name
const PetscInt _pd_name_sz_ = 260;

// length of scaling unit label
const PetscInt _lbl_sz_ = 23;

// maximum number of sedimentary layers (free surface sedimentation)
const PetscInt _max_sed_layers_ = 50;

// maximum number of erosion phases (free surface erosion)
const PetscInt _max_er_phases_ = 50;

// maximum number of adjoint parameters
const PetscInt _max_adj_par_ = 50;

// maximum number of adjoint points
const PetscInt _max_adj_point_ = 100;

// maximum number of passive tracers
const PetscInt _max_passive_tracer = 100000;

// maximum number of control polygons
const PetscInt _max_ctrl_poly_ = 20;

// maximum number of multigrid levels
const PetscInt _max_num_mg_levels_ = 16;

// maximum number of matrix-free levels
const PetscInt _max_num_mat_free_levels_ = 8;

// adjoint parameter limits
const PetscInt _MAX_PAR_ = 100;
const PetscInt _MAX_OBS_ = 100;

//-----------------------------------------------------------------------------
// TYPE DEFINITIONS
//-----------------------------------------------------------------------------

typedef pair <PetscScalar, PetscInt> spair;
typedef pair <PetscInt,    PetscInt> ipair;

//-----------------------------------------------------------------------------
// UNUSED PARAMETERS MACRO
//-----------------------------------------------------------------------------

#define UNUSED(x) (void)(x)

//-----------------------------------------------------------------------------
// PRINT FORMAT MACRO
//-----------------------------------------------------------------------------

#define PetscMPIInt_FMT "d"

//-----------------------------------------------------------------------------
// PROTOTYPES
//-----------------------------------------------------------------------------

struct FB;

// LaMEM library main function

PetscErrorCode LaMEMLibMain(void *param, FB *fb);

//-----------------------------------------------------------------------------
#endif
