/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

#ifndef __bfbt_h__
#define __bfbt_h__



struct _p_PMat;
typedef struct _p_PMat *PMat;
struct AdvCtx;
struct FB;
struct GeomPrim;
struct JacRes;

struct SphereData{
	PetscScalar centerx;
	PetscScalar centery;
	PetscScalar centerz;
	PetscScalar radius;
	PetscInt	phase;
};



//---------------------------------------------------------------------------
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------

PetscErrorCode PMatBFBTCreate(PMat pm);

PetscErrorCode PMatBFBTAssemble(PMat pm);

PetscErrorCode PMatBFBTDestroy(PMat pm);

PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y);

//PetscErrorCode BFBTGaussianSmoothing(AdvCtx *actx, FB *fb, GeomPrim *geom, AllGeom1 blubb);
PetscErrorCode BFBTGaussianSmoothing(JacRes *jr);

//---------------------------------------------------------------------------
// compute |a-b|, a,b are vectors
#define VEC_ABS(vecabs, a, b) {vecabs = sqrt( (a[2]-b[2])*(a[2]-b[2]) + (a[1]-b[1])*(a[1]-b[1]) + (a[0]-b[0])*(a[0]-b[0]) );}

//---------------------------------------------------------------------------
#endif
