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
 **    filename:   LaMEM.h
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
//.........................   Main include file   ...........................
//---------------------------------------------------------------------------

#ifndef __LaMEM_h__
#define __LaMEM_h__

//-----------------------------------------------------------------------------
//   PREFERABLE VARIABLES
//
//   PetscInt    - for all indices                    (can be int or long long int)
//   PetscScalar - for all floating point variables   (can be float or double)
//   float       - for reduced size output
//   size_t      - for all sizes offsets & counters   (unsigned long long int)
//   PetscMPIInt - for passing service integer parameters to MPI functions (int)
//   MPIU_SCALAR - appropriate MPI Data Type for sending/receiving PetsScalar
//   MPIU_INT    - appropriate MPI Data Type for sending/receiving PetscInt
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SIZE LIMITS
//-----------------------------------------------------------------------------

// number of neighbor domains in 3D lattice (including self)
#define _num_neighb_ 27

// string length (two null characters are reserved in the end, i.e. 128)
#define _str_len_ 130

// marker storage capacity overhead
#define _cap_overhead_ 1.61803398875

// maximum marker per cell per direction
#define _max_nmark_ 5

// minimum marker per cell per direction
#define _min_nmark_ 2

// cell marker buffer size
#define _mark_buff_sz_ 256

// local marker buffer size as percentage of local number of markers
#define _mark_buff_ratio_ 5

// maximum number of strain rate application periods
#define _max_periods_ 20

// maximum number of Bezier blocks
#define _max_boxes_ 5

// maximum number of Bezier path points
#define _max_path_points_ 25

// maximum number of polygon points of Bezier block
#define _max_poly_points_ 50

// maximum number of mesh segments in every direction
#define _max_num_segs_ 10

// maximum number of cells per mesh segment
#define _max_num_cells_ 4096

// maximum number of processes in every direction
#define _max_num_procs_ 1024

// FDSTAG near null space size
#define _max_nullsp_sz_ 4

// maximum number of geometry primitives
#define _max_geom_ 100

// maximum number of observation types
#define _max_num_obs_ 7

// maximum number of components in the output vector (3D)
#define _max_num_comp_ 9

// maximum number of components in the output vector (surface)
#define _max_num_comp_surf_ 3

// maximum number of phase aggregates for output
#define _max_num_phase_agg_ 5

// maximum number of phases
#define _max_num_phases_ 32

// maximum number of softening laws
#define _max_num_soft_ 10

// max char length
#define char_ph_tr 10

// maximum number of phase transition law
#define _max_num_tr_ 20

// maximum number of segments of NotInAirBoxes
#define _max_NotInAir_segs_ 6

// maximum number of dikes
#define _max_num_dike_ 6

// maximum number of equation parameter
#define _max_num_eq_ 2

// maximum number of phase transition per each phase
#define _max_tr_ 5

// maximum number of phase diagrams
#define _max_num_pd_ 8

// maximum grid size of phase diagram
#define _max_pd_sz_ 40100

// length of unique phase diagram name
#define _pd_name_sz_ 54

// length of scaling unit label
#define _lbl_sz_ 23

// maximum number of sedimentary layers (free surface sedimentation)
#define _max_sed_layers_ 50

// maximum number of erosion phases (free surface erosion)
#define _max_er_phases_ 50

// maximum number of adjoint parameters
#define _max_adj_par_ 50

// maximum number of adjoint points
#define _max_adj_point_ 100

// maximum number of passive tracers
#define _max_passive_tracer 100000

// maximum number of control polygons
#define _max_ctrl_poly_ 20

// cast macros
#define LLD long long int

//-----------------------------------------------------------------------------
// EXTERNAL INCLUDES
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <sys/stat.h>
#include <errno.h>
#include <petsc.h>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;

//-----------------------------------------------------------------------------
// TYPE DEFINITIONS
//-----------------------------------------------------------------------------

typedef pair <PetscScalar, PetscInt> spair;
typedef pair <PetscInt,    PetscInt> ipair;

//-----------------------------------------------------------------------------
// PROTOTYPES
//-----------------------------------------------------------------------------

// LaMEM library main function

PetscErrorCode LaMEMLibMain(void *param);

//-----------------------------------------------------------------------------
#endif
