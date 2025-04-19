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
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------
#ifndef __bfbt_h__
#define __bfbt_h__
//---------------------------------------------------------------------------
struct Pmat;
//---------------------------------------------------------------------------

PetscErrorCode PMatBFBTCreate(PMat pm);

PetscErrorCode PMatBFBTAssemble(PMat pm);

PetscErrorCode PMatBFBTDestroy(PMat pm);

PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------
#endif
