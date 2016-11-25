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

// read material parameter limits
PetscErrorCode MatParLimRead(
		FB        *fb,
		Scaling   *scal,
		MatParLim *matLim);

// read all material phases and softening laws from file
PetscErrorCode MatPropsReadAll(
		FB         *fb,
		Scaling    *scal,
		PetscInt   *numPhases,
		Material_t *phases,
		PetscInt   *numSoft,
		Soft_t     *matSoft);

// read single softening law
PetscErrorCode MatSoftRead(
		FB       *fb,
		PetscInt  numSoft,
		Soft_t   *matSoft);

// read single material phase
PetscErrorCode MatPhaseRead(
		FB         *fb,
		Scaling    *scal,
		PetscInt    numPhases,
		Material_t *phases,
		PetscInt    numSoft,
		Soft_t     *matSoft);

static inline void MatPrintScalParam(PetscScalar par, const char key[], const char label[], Scaling *scal)
{
	if(par == 0.0) return;

	if(scal->utype == _NONE_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s = %g [ ]  ", key, par);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "%s = %g [] %s  ", key, par, label);
	}
}

//---------------------------------------------------------------------------
//............ PREDEFINED RHEOLOGICAL PROFILES (from literature) ............
//---------------------------------------------------------------------------
typedef enum
{
	_UniAxial_,      // Uni-axial experiment
	_SimpleShear_,   // Simple shear experiment
	_None_           // geological-scale units

} TensorCorrection;

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
