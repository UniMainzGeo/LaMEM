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
 **    filename:   outFunct.h
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
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __outFunct_h__
#define __outFunct_h__
//---------------------------------------------------------------------------

struct OutBuf;
struct JacRes;

//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------

typedef struct OutVec OutVec;

struct OutVec
{
	JacRes   *jr;
	OutBuf   *outbuf;
	char     *name;    // output vector name
	PetscInt  ncomp;   // number of components

//	PetscInt  phase_mask[_max_num_phases_];

	PetscErrorCode (*OutVecWrite)(OutVec*); // output function pointer

};

void OutVecCreate(
	OutVec         *outvec,
	JacRes         *jr,
	OutBuf         *outbuf,
	const char     *name,
	const char     *label,
	PetscErrorCode (*OutVecWrite)(OutVec*),
	PetscInt        ncomp);

void OutVecDestroy(OutVec *outvec);

//---------------------------------------------------------------------------

PetscErrorCode PVOutWritePhase       (OutVec*);
PetscErrorCode PVOutWriteDensity     (OutVec*);
PetscErrorCode PVOutWriteViscTotal   (OutVec*);
PetscErrorCode PVOutWriteViscCreep   (OutVec*);
PetscErrorCode PVOutWriteViscoPlastic(OutVec*);
PetscErrorCode PVOutWriteVelocity    (OutVec*);
PetscErrorCode PVOutWritePressure    (OutVec*);
PetscErrorCode PVOutWriteTotalPress  (OutVec*);
PetscErrorCode PVOutWriteEffPress    (OutVec*);
PetscErrorCode PVOutWriteOverPress   (OutVec*);
PetscErrorCode PVOutWriteLithoPress  (OutVec*);
PetscErrorCode PVOutWritePorePress   (OutVec*);
PetscErrorCode PVOutWriteTemperature (OutVec*);
PetscErrorCode PVOutWriteDevStress   (OutVec*);
PetscErrorCode PVOutWriteJ2DevStress (OutVec*);
PetscErrorCode PVOutWriteStrainRate  (OutVec*);
PetscErrorCode PVOutWriteJ2StrainRate(OutVec*);
PetscErrorCode PVOutWriteMeltFraction(OutVec*);
PetscErrorCode PVOutWriteFluidDensity(OutVec*);
PetscErrorCode PVOutWriteVolRate     (OutVec*);
PetscErrorCode PVOutWriteVorticity   (OutVec*);
PetscErrorCode PVOutWriteAngVelMag   (OutVec*);
PetscErrorCode PVOutWriteTotStrain   (OutVec*);
PetscErrorCode PVOutWritePlastStrain (OutVec*);
PetscErrorCode PVOutWritePlastDissip (OutVec*);
PetscErrorCode PVOutWriteTotDispl    (OutVec*);
PetscErrorCode PVOutWriteSHmax       (OutVec*);
PetscErrorCode PVOutWriteEHmax       (OutVec*);
PetscErrorCode PVOutWriteISA         (OutVec*);
PetscErrorCode PVOutWriteGOL         (OutVec*);
PetscErrorCode PVOutWriteYield       (OutVec*);
// === debug vectors ===============================================
PetscErrorCode PVOutWriteMomentRes   (OutVec*);
PetscErrorCode PVOutWriteContRes     (OutVec*);
PetscErrorCode PVOutWritEnergRes     (OutVec*);
PetscErrorCode PVOutWriteMeltFraction(OutVec*);
PetscErrorCode PVOutWriteFluidDensity(OutVec*);

// ... add more output functions here

//---------------------------------------------------------------------------
#endif
