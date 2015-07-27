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
 **    filename:   dfzero.h
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
//..................... DFZERO ROOT FINDING ALGORITHM .......................
//---------------------------------------------------------------------------
#ifndef __dfzero_h__
#define __dfzero_h__
//---------------------------------------------------------------------------
static inline PetscScalar FDMIN(PetscScalar a, PetscScalar b)
{
    if(a <= b) return a;
    else       return b;
}
//---------------------------------------------------------------------------
static inline PetscScalar FDMAX(PetscScalar a, PetscScalar b)
{
    if(a >= b) return a;
    else       return b;
}
//---------------------------------------------------------------------------
static inline PetscScalar FDSIGN(PetscScalar x, PetscScalar y)
{
    if(y >= 0.0) return  PetscAbsScalar(x);
    else         return -PetscAbsScalar(x);
}
//---------------------------------------------------------------------------
void DFZERO(
	PetscScalar (*F)(PetscScalar, void *), // nonlinear function with parameter and context
	void         *FCTX,                    // pointer to a function evaluation context
	PetscScalar  *_B,                      // left bound of root interval (output)
	PetscScalar  *_C,                      // right bound of root interval (output)
	PetscScalar    R,                      // initial guess
	PetscScalar    RE,                     // relative tolerance
	PetscScalar    AE,                     // absolute tolerance
	PetscInt     *_IFLAG);                 // error code (output)
//---------------------------------------------------------------------------
#endif
