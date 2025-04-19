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
//..................... MATRIX EVALUATION CONTEXT  .........................
//---------------------------------------------------------------------------
#ifndef __matData_h__
#define __matData_h__
//---------------------------------------------------------------------------

struct JacRes;
struct FDSTAG;

//---------------------------------------------------------------------------

// preconditioning matrix storage format
enum PMatType
{
	_MONOLITHIC_,
	_BLOCK_
};

//---------------------------------------------------------------------------

struct MatData
{
	PMatType    type;                             // matrix type
	FDSTAG     *fs;                               // staggered grid
	Vec         ivx, ivy, ivz, ip;                // index vectors
	Vec         bcvx, bcvy, bcvz, bcp;            // boundary condition vectors
	PetscInt    numSPC,  *SPCList;                // single points constraints (SPC)
	PetscInt    vNumSPC, *vSPCList;               // velocity SPC
	PetscInt    pNumSPC, *pSPCList;               // pressure SPC
	Vec         K, rho, eta, etaxy, etaxz, etayz; // parameter vectors
	PetscScalar dt;                               // time step
	PetscScalar fssa;                             // density gradient penalty parameter
	PetscScalar grav[3];                          // global gravity components
	PetscScalar pgamma;                           // penalty parameter
	PetscInt    coarsened;                        // coarsening flag
};

//---------------------------------------------------------------------------

PetscErrorCode MatDataSetFromOptions(MatData *md);

PetscErrorCode MatDataCreate(MatData *md, JacRes *jr);

PetscErrorCode MatDataInitParam(MatData *md, JacRes *jr);

PetscErrorCode MatDataDestroy(MatData *md);

PetscErrorCode MatDataCreateIndex(MatData *md);

PetscErrorCode MatDataCoarsen(MatData *coarse, MatData *fine);

PetscErrorCode MatDataRestricParam(MatData *coarse, MatData *fine);

PetscErrorCode MatDataRestricBC(MatData *coarse, MatData *fine);





//---------------------------------------------------------------------------

#endif
