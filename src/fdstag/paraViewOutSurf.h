//---------------------------------------------------------------------------
//..............   FREE SRUFACE PARAVIEW XML OUTPUT ROUTINES   ..............
//---------------------------------------------------------------------------
#ifndef __paraViewOutSurf_h__
#define __paraViewOutSurf_h__
//---------------------------------------------------------------------------
// maximum number of components in the output vector
#define _max_num_comp_surf_ 3
//---------------------------------------------------------------------------
//................ ParaView free surface output driver object ...............
//---------------------------------------------------------------------------
typedef struct
{
	FreeSurf *surf;       // free surface object
	char     *outfile;    // output file name
	float    *buff;       // direct output buffer
	long int  offset;     // pvd file offset
	PetscInt  outpvd;     // pvd file output flag
	PetscInt  velocity;   // velocity output flag
	PetscInt  topography; // surface topography output flag
	PetscInt  amplitude;  // topography amplitude output flag

} PVSurf;
//---------------------------------------------------------------------------

// clear object
PetscErrorCode PVSurfClear(PVSurf *pvout);

// create ParaView output driver
PetscErrorCode PVSurfCreate(PVSurf *pvout, FreeSurf *surf, const char *filename);

// read options
PetscErrorCode PVSurfReadFromOptions(PVSurf *pvout);

// destroy ParaView output driver
PetscErrorCode PVSurfDestroy(PVSurf *pvout);

// write all time-step output files to disk (PVD, PVTS, VTS)
PetscErrorCode PVSurfWriteTimeStep(PVSurf *pvout, const char *dirName, PetscScalar ttime, PetscInt tindx);

// parallel output file .pvts
PetscErrorCode PVSurfWritePVTS(PVSurf *pvout, const char *dirName);

// sequential output file .vts
PetscErrorCode PVSurfWriteVTS(PVSurf *pvout, const char *dirName);

//---------------------------------------------------------------------------
#endif
