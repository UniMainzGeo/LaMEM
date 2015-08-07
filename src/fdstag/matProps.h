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
 **    filename:   matProps.h
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
#ifndef __matProps_h__
#define __matProps_h__
//---------------------------------------------------------------------------
//........................... MATERIAL PARAMETERS ...........................
//---------------------------------------------------------------------------
// read all phases
PetscErrorCode MatPropInit(JacRes *jr, FILE *fp);

// read single phase
PetscErrorCode MatPropGetStruct(FILE *fp,
	PetscInt numPhases, Material_t *phases,
	PetscInt numSoft,   Soft_t     *matSoft,
	PetscInt ils, PetscInt ile, UnitsType utype);

// read phases from command line
PetscErrorCode MatPropSetFromCL(JacRes *jr);

// assign phases from calling function
PetscErrorCode MatPropSetFromLibCall(JacRes *jr, ModParam *mod);

//---------------------------------------------------------------------------
//............................ SOFTENING LAWS ...............................
//---------------------------------------------------------------------------

// read all softening laws
PetscErrorCode MatSoftInit(JacRes *jr, FILE *fp);

// read single softening laws
PetscErrorCode MatSoftGetStruct(FILE *fp,
	PetscInt numSoft, Soft_t *matSoft,
	PetscInt ils, PetscInt ile);

//---------------------------------------------------------------------------
//............ PREDEFINED RHEOLOGICAL PROFILES (from literature) ............
//---------------------------------------------------------------------------
typedef enum
{
	_UniAxial_,      // Uni-axial experiment
	_SimpleShear_,   // Simple shear experiment
	_None_           // geological-scale units
} TensorCorrection;

// diffusion creep profiles
PetscErrorCode SetDiffProfile(Material_t *m, char name[]);

// dislocation creep profiles
PetscErrorCode SetDislProfile(Material_t *m, char name[]);

// Peierls creep profiles
PetscErrorCode SetPeirProfile(Material_t *m, char name[]);

// units and tensor correction
PetscErrorCode SetProfileCorrection(PetscScalar *B, PetscScalar n, TensorCorrection tensorCorrection, PetscInt MPa);

//---------------------------------------------------------------------------
//................ Routines to get structure-info from file .................
//---------------------------------------------------------------------------

// gets the file positions of a structure
void getLineStruct(
	FILE *fp, PetscInt *ls, PetscInt *le, PetscInt mux_num,
	PetscInt *count_starts, PetscInt *count_ends,
	const char key[], const char key_end[]);

// gets an integer within specified positions in file
void getMatPropInt(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscInt *value, PetscInt *found);

// gets a scalar within specified positions in file
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscScalar *value, PetscInt *found);

// gets a string within specified positions in file
void getMatPropString(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], char value[], PetscInt max_L, PetscInt *found );

//---------------------------------------------------------------------------

#endif
