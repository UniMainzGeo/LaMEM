//---------------------------------------------------------------------------
//..............   MARKER PARAVIEW XML OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __paraViewOutMark_h__
#define __paraViewOutMark_h__
//---------------------------------------------------------------------------
//................ ParaView marker output driver object .....................
//---------------------------------------------------------------------------
typedef struct
{
	AdvCtx     *actx;       // free surface object
	char       *outfile;    // output file name
	long int    offset;     // pvd file offset
	PetscInt    outmark;    // pvd file output flag
	PetscInt    outpvd;     // pvd file output flag
} PVMark;
//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVMarkCreate(PVMark *pvmark, AdvCtx *actx, const char *filename);

// read options
PetscErrorCode PVMarkReadFromOptions(PVMark *pvmark);

// write all time-step output files to disk (PVD, PVTS, VTS)
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime, PetscInt tindx);

// vtu marker output
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName);

// pvtu marker output
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName);

//---------------------------------------------------------------------------

#endif
