/*
 * BFBT.h
 *
 *  Created on: 04.03.2020
 *      Author: daniel
 */

//---------------------------------------------------------------------------
#ifndef __bfbt_h__
#define __bfbt_h__
//---------------------------------------------------------------------------
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------

PetscErrorCode PMatBFBTCreate(PMat pm);

PetscErrorCode PMatBFBTAssemble(PMat pm);

PetscErrorCode PMatBFBTDestroy(PMat pm);

PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------
#endif
