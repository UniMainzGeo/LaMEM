//---------------------------------------------------------------------------
//......................   FDSTAG GMT OUTPUT ROUTINES   .....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "scaling.h"
#include "gmtOut.h"
#include "outFunct.h"
#include "Utils.h"
//---------------------------------------------------------------------------
// * comments go into here
// * ...
//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "gmtInterface"
PetscErrorCode gmtInterface()
{
	// GMT
	void *API = NULL;
	void ID;
	void fileID;

	// LaMEM
	FDSTAG   *fs;
	PetscScalar *data;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// Start a GMT session
	API = GMT_Create_Session ("Session name", 2, 0, NULL);
    if(API) {
        PetscPrintf(PETSC_COMM_WORLD,"Created GMT API!\n");
    }

    // Register in or output resources (Here: output entire 2D grid file)
    ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE,  GMT_IS_SURFACE, GMT_OUT,  NULL, &fileID);

    // Object ID encoding
    GMT_Encode_ID (API, "2dGrid.grd", ID);

    // Read additional command line arguments ??
    //GMT_Init_IO (void *API, unsigned int family, unsigned int geometry, unsigned int direction, unsigned int mode, unsigned int n_args, void *args);

	// Use data stored in memory
	data = GMT_Get_Data (API, ID, GMT_GRID_DATA_ONLY, NULL);

	// write data (to parallelize: GMT_WRITE_SEGMENT)
	GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, &fileID, data);

	// Destroy allocated data ?
	GMT_Destroy_Data (API, data);

	// Terminate GMT session
	GMT_Destroy_Session (API);

	PetscFunctionReturn(0);
}
//---------
