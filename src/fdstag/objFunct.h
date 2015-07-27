/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This sofware was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   objFunct.h
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
//.....................   OBJECTIVE FUNCTION ROUTINES   .....................
//---------------------------------------------------------------------------
#ifndef __objFunct_h__
#define __objFunct_h__
//---------------------------------------------------------------------------
// maximum number of obervational types
#define _max_num_obs_ 7
#define _max_len_name_ 8


typedef enum 	// observation type
{
	_VELX_,            // 1: horizontal velocity (x) at the surface
	_VELY_,            // 2: horizontal velocity (y) at the surface
	_VELZ_,            // 3: vertical velocity (z) at the surface
	_TOPO_,            // 4: surface topography
	_BOUG_,            // 5: bouguer anomaly
	_ISA_,             // 6: orientation of isa (<-> sks-seismic anisotropy)
	_SHMAX_            // 7: orientation of SHmax
} ObsType;

//---------------------------------------------------------------------------
//........................ Objective function object ........................
//---------------------------------------------------------------------------
typedef struct
{
	FreeSurf     *surf;                 // free surface object
	char         *infile;               // input file name
	PetscBool    CompMfit;              // Compute misfit?
	PetscInt     otUse[_max_num_obs_+1];// array of boolean USED flags
	PetscInt     otN;                   // number of USED observation types
	PetscInt     ocN;                   // total number of observational constraints
	PetscScalar  err[_max_num_obs_];    // array containing individual sums of errors
	PetscScalar  errtot;                // total error
	Vec          obs[_max_num_obs_];    // vectors containing the observations
	Vec          qul[_max_num_obs_];    // vectors containing quality info (quality/sigma)^2, where quality (0..1)

	// missing ...
	// (data) covariance matrix

} ObjFunct;
//---------------------------------------------------------------------------

// destroy object
PetscErrorCode ObjFunctDestroy(ObjFunct *objf);

// create objective function object
PetscErrorCode ObjFunctCreate(ObjFunct *objf, FreeSurf *surf);

// read command line options
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[]);

// compute error
PetscErrorCode ObjFunctCompErr(ObjFunct *objf);

// compute weighted Least square error for surface vectors
PetscErrorCode VecErrSurf(Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal);

//---------------------------------------------------------------------------

#endif
