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
#ifndef __paraViewOutMark_h__
#define __paraViewOutMark_h__
//---------------------------------------------------------------------------
//................ ParaView marker output driver object .....................
//---------------------------------------------------------------------------

struct FB;
struct AdvCtx;

//---------------------------------------------------------------------------

struct PVMark
{
	AdvCtx    *actx;              // advection context
	char      outfile[_str_len_+20]; // output file name
	long int  offset;             // pvd file offset
	PetscInt  outmark;            // marker output flag
	PetscInt  outpvd;             // pvd file output flag

};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVMarkCreate(PVMark *pvmark, FB *fb);

// write all time-step output files to disk (PVD, PVTU, VTU)
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime);

// .vtu marker output
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName);

// .pvtu marker output
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName);

//---------------------------------------------------------------------------

#endif
