/*
 * BFBT.h
 *
 *  Created on: 04.03.2020
 *      Author: daniel
 */

//---------------------------------------------------------------------------
#ifndef __bfbt_h__
#define __bfbt_h__



struct _p_PMat;
typedef struct _p_PMat *PMat;
struct JacRes;
struct GeomPrim;



//---------------------------------------------------------------------------
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------

PetscErrorCode PMatBFBTCreate(PMat pm);

PetscErrorCode PMatBFBTAssemble(PMat pm);

PetscErrorCode PMatBFBTDestroy(PMat pm);

PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y);

PetscErrorCode BFBTGaussianSmoothing(JacRes *jr, PetscInt nSpheres, GeomPrim *geom);

//---------------------------------------------------------------------------
// compute |a-b|, a,b are vectors
#define VEC_ABS(vecabs, a, b) {vecabs = sqrt( (a[2]-b[2])*(a[2]-b[2]) + (a[1]-b[1])*(a[1]-b[1]) + (a[0]-b[0])*(a[0]-b[0]) );}

//---------------------------------------------------------------------------
#endif
