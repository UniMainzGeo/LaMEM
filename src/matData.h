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

// matrix type
enum PMatType
{
	_MONOLITHIC_,
	_BLOCK_
};

//---------------------------------------------------------------------------

struct MatData
{
	FDSTAG     *fs;                                // staggered grid
	Vec         ivx, ivy, ivz, ip;                 // index vectors
	Vec         bcvx, bcvy, bcvz, bcp;             // boundary condition vectors
	PetscInt    numSPC,  *SPCList;                 // single points constraints (SPC)
	PetscInt    vNumSPC, *vSPCList;                // velocity SPC
	PetscInt    pNumSPC, *pSPCList;                // pressure SPC
	Vec         Kb, rho, eta, etaxy, etaxz, etayz; // parameter vectors
	PetscScalar dt;                                // time step
	PetscScalar fssa;                              // density gradient penalty parameter
	PetscScalar grav[3];                           // global gravity components
	PetscInt    coarse;                            // coarsening flag
	PMatType    type;                              // matrix type
	PetscScalar pgamma;                            // penalty parameter
};

//---------------------------------------------------------------------------

PetscErrorCode MatDataSetFromOptions(MatData *md);

PetscErrorCode MatDataCreate(MatData *md, JacRes *jr);

PetscErrorCode MatDataCreateData(MatData *md);

PetscErrorCode MatDataDestroy(MatData *md);

PetscErrorCode MatDataCoarsen(MatData *coarse, MatData *fine);

PetscErrorCode MatDataComputeIndex(MatData *md);

PetscErrorCode MatDataInitParam(MatData *md, JacRes *jr);

PetscErrorCode MatDataRestricParam(MatData *coarse, MatData *fine);

PetscErrorCode MatDataRestricBC(MatData *coarse, MatData *fine);

PetscErrorCode MatDataListSPC(MatData *md);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define LIST_SPC_IND(bc, list, cnt, iter) \
	if(bc[k][j][i] != DBL_MAX) { list[cnt] = iter; cnt++; }

//---------------------------------------------------------------------------


#endif
