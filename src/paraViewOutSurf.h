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
//..............   FREE SRUFACE PARAVIEW XML OUTPUT ROUTINES   ..............
//---------------------------------------------------------------------------
#ifndef __paraViewOutSurf_h__
#define __paraViewOutSurf_h__

//---------------------------------------------------------------------------

struct FB;
struct FreeSurf;

//---------------------------------------------------------------------------
//................ ParaView free surface output driver object ...............
//---------------------------------------------------------------------------

struct PVSurf
{
	FreeSurf  *surf;               // free surface object
	char       outfile[_str_len_+20]; // output file name
	char       outfile_refine[_str_len_+20]; // output file name
	float     *buff;               // direct output buffer
	float     *buff_refine;        // direct output buffer
	long int   offset;             // pvd file offset
	long int   offset_refine;             // pvd file offset
	PetscInt   outsurf;            // free surface output flag
	PetscInt   outsurf_refine;            // free surface output flag
	PetscInt   outpvd;             // pvd file output flag
	PetscInt   velocity;           // velocity output flag
	PetscInt   topography;         // surface topography output flag
	PetscInt   amplitude;          // topography amplitude output flag
	PetscInt   topoRefine;         // refined topography output flag
};

//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVSurfCreate(PVSurf *pvsurf, FB *fb);

// create buffer array
PetscErrorCode PVSurfCreateData(PVSurf *pvsurf);

// destroy ParaView output driver
PetscErrorCode PVSurfDestroy(PVSurf *pvsurf);

// write all time-step output files to disk (PVD, PVTS, VTS)
PetscErrorCode PVSurfWriteTimeStep(PVSurf *pvsurf, const char *dirName, PetscScalar ttime);

// parallel output file .pvts
PetscErrorCode PVSurfWritePVTS(PVSurf *pvsurf, const char *dirName);

// sequential output file .vts
PetscErrorCode PVSurfWriteVTS(PVSurf *pvsurf, const char *dirName);

// sequential output file .vts
PetscErrorCode PVSurfWriteVTSRefine(PVSurf *pvsurf, const char *dirName, PetscInt nx_refine, PetscInt ny_refine, 
    PetscScalar rangeX, PetscScalar rangeY, PetscScalar *topo);

//---------------------------------------------------------------------------

void OutputBufferWrite(
	FILE     *fp,
	float    *buff,
	PetscInt  cn);

PetscErrorCode PVSurfWriteCoord(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteCoordRefine(PVSurf *pvsurf, FILE *fp, PetscInt nx_refine, PetscInt ny_refine, 
    PetscScalar rangeX, PetscScalar rangeY, PetscScalar *topo);

PetscErrorCode PVSurfWriteVel(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteTopo(PVSurf *pvsurf, FILE *fp);

PetscErrorCode PVSurfWriteTopoRefine(PVSurf *pvsurf, FILE *fp, PetscInt nx_refine, PetscInt ny_refine, PetscScalar *topo);

PetscErrorCode PVSurfWriteAmplitude(PVSurf *pvsurf, FILE *fp);

//---------------------------------------------------------------------------
#endif
