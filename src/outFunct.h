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

PetscErrorCode PVOutWritePhase         (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDensity       (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteViscTotal     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteViscCreep     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteViscoPlastic  (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVelocity      (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePressure      (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteOverPressure  (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteLithosPressure(JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTemperature   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDevStress     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteJ2DevStress   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteStrainRate    (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteJ2StrainRate  (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVolRate       (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVorticity     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteAngVelMag     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTotStrain     (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePlastStrain   (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePlastDissip   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTotDispl      (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteSHmax         (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteEHmax         (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteISA           (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteGOL           (JacRes*, OutBuf*);
// === debug vectors ===============================================
PetscErrorCode PVOutWriteJacTest     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteMomentRes   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteContRes     (JacRes*, OutBuf*);
PetscErrorCode PVOutWritEnergRes     (JacRes*, OutBuf*);

// ... add more output functions here

//---------------------------------------------------------------------------

#endif
