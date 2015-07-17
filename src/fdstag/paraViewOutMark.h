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
	AdvCtx     *actx;       // advection context
	char       *outfile;    // output file name
	long int    offset;     // pvd file offset
	PetscInt    outmark;    // marker output flag
	PetscInt    outpvd;     // pvd file output flag
} PVMark;
//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVMarkCreate(PVMark *pvmark, AdvCtx *actx, const char *filename);

// free memory
PetscErrorCode PVMarkDestroy(PVMark *pvmark);

// read options
PetscErrorCode PVMarkReadFromOptions(PVMark *pvmark);

// write all time-step output files to disk (PVD, PVTU, VTU)
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime, PetscInt tindx);

// .vtu marker output
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName);

// .pvtu marker output
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName);

//---------------------------------------------------------------------------

#endif
