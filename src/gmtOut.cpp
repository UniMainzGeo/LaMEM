/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   gmtOut.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//......................   FDSTAG GMT OUTPUT ROUTINES   .....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "gmtOut.h"


//---------------------------------------------------------------------------
// * comments go into here
// * ...
//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
/*
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
	PetscFunctionBeginUser;


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
*/
//---------------------------------------------------------------------------
