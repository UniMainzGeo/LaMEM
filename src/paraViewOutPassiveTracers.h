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
//..............   MARKER PARAVIEW XML OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __paraViewOutPassiveTracers_h__
#define __paraViewOutPassiveTracers_h__
//---------------------------------------------------------------------------
//................ ParaView marker output driver object .....................
//---------------------------------------------------------------------------

struct FB;
struct AdvCtx;

//---------------------------------------------------------------------------

struct PVPtr
{
	AdvCtx    *actx;              // advection context
	char      outfile[_str_len_+20]; // output file name
	long int  offset;             // pvd file offset
	PetscInt  outptr;             // marker output flag
	PetscInt  outpvd;             // pvd file output flag
	PetscInt  Temperature;
	PetscInt  Pressure;
	PetscInt  Phase;
	PetscInt  MeltFraction;
	PetscInt  ID;
	PetscInt  Active;
	PetscInt  Grid_mf;

};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVPtrCreate(PVPtr *pvptr, FB *fb);

// write all time-step output files to disk (PVD, PVTU, VTU)
PetscErrorCode PVPtrWriteTimeStep(PVPtr *pvptr, const char *dirName, PetscScalar ttime);

// .vtu marker output
PetscErrorCode PVPtrWriteVTU(PVPtr *pvptr, const char *dirName);

// .pvtu marker output
PetscErrorCode PVPtrWritePVTU(PVPtr *pvptr, const char *dirName);


#endif
