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
 **    filename:   input.c
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
//....................... FILE INPUT / INITIALIZATION .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "input.h"
#include "matProps.h"
#include "fdstagTypes.h"
//---------------------------------------------------------------------------
// set default code parameters and read input file, if required
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitCode"
PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user, ModParam *iop)
{
	FILE      *fp;
	char      *all_options;
	char      ParamFile[PETSC_MAX_PATH_LEN];
	PetscBool InputParamFile;

	PetscMPIInt  size;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED
	PetscOptionsGetAll( &all_options ); // copy all command line args

	// set default values for parameters
	ierr = InputSetDefaultValues(jr, user); CHKERRQ(ierr);

	// check whether input file is specified
	ierr = PetscOptionsGetString(PETSC_NULL, "-ParamFile", ParamFile, PETSC_MAX_PATH_LEN, &InputParamFile); CHKERRQ(ierr);

	// read additional PETSc options from input file if required
	if(InputParamFile == PETSC_TRUE)
	{
		// open file
		fp = fopen(ParamFile, "r");

		if(!fp)
		{
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot open input file %s", ParamFile);
		}
		else
		{
			// read an input file using the parser
			ierr = InputReadFile(jr, user, fp); CHKERRQ(ierr);

			// read softening laws from file
			ierr = MatSoftInit(jr, fp); CHKERRQ(ierr);

			// initialize and read material properties from file
			ierr = MatPropInit(jr, fp); CHKERRQ(ierr);

			// close file
			fclose(fp);
		}
	}

	// Interpret command line options and overwrite material properties
	ierr = MatPropSetFromCL(jr); CHKERRQ(ierr);

	// Interpret IO parameter structure to set material properties
	ierr = MatPropSetFromLibCall(jr, iop); CHKERRQ(ierr);

	// Interpret command line options
	ierr = InputReadCommLine(user); CHKERRQ(ierr);

	// error check and info print:
	// we need at least 2 cells in each direction (not 1) and dt > 0 ! - IS THIS NECESSARY?
	if (user->nel_x==1){ user->nel_x=2; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the x-direction! I have increased this. \n");}
	if (user->nel_y==1){ user->nel_y=2; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the y-direction! I have increased this. \n");}
	if (user->dt==0.0) { user->dt=1e-3; PetscPrintf(PETSC_COMM_WORLD," With FDSTAG dt cannot be zero. Value has been changed to dt = 1e-3. \n");}

	// print information about simulation
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Total # of cpu's               : %lld \n",(LLD)size);

	// info about particles if used
	if (user->ParticleInput==1)
	{
		// check whether a sufficient amount of particles are specified; if not increase it.
		if(user->NumPartX < 2) { user->NumPartX = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in x-direction; Increasing\n"); }
		if(user->NumPartY < 2) { user->NumPartY = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in y-direction; Increasing\n"); }
		if(user->NumPartZ < 2) { user->NumPartZ = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in z-direction; Increasing\n"); }
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of tracers/cell         : [%lld,%lld,%lld] \n",(LLD)(user->NumPartX),(LLD)(user->NumPartY),(LLD)(user->NumPartZ));

	// show initial timestep
	if(jr->scal.utype == _GEO_)
	{
		PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g [MYr] \n",user->dt);
	}
	else if(jr->scal.utype == _SI_)
	{
		PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g  [s] \n",user->dt);
	}
	else if(jr->scal.utype == _SI_)
	{
		PetscPrintf(PETSC_COMM_WORLD," Initial time step              : %g  [ ] \n",user->dt);
	}

	// boundary conditions
	PetscPrintf(PETSC_COMM_WORLD," BC employed                    : BC.[LeftBound=%lld RightBound=%lld; FrontBound=%lld BackBound=%lld; LowerBound=%lld UpperBound=%lld] \n",
			(LLD)(user->BC.LeftBound), (LLD)(user->BC.RightBound), (LLD)(user->BC.FrontBound), (LLD)(user->BC.BackBound), (LLD)(user->BC.LowerBound), (LLD)(user->BC.UpperBound) );

	ierr = PetscFree(all_options); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// set default values for parameters
#undef __FUNCT__
#define __FUNCT__ "InputSetDefaultValues"
PetscErrorCode InputSetDefaultValues(JacRes *jr, UserCtx *user)
{
	// domain
	user->W                = 1.0;
	user->L                = 1.0;
	user->H                = 1.0;
	user->x_left           = -0.5*(user->W)*0.0;
	user->y_front          = -0.5*(user->L)*0.0;
	user->z_bot            = -0.5*(user->H)*0.0;

	// set default grid sizes
	user->nel_x            = 8;
	user->nel_y            = 8;
	user->nel_z            = 8;

	user->NumPartX         = 2; // Number of particles per cell in x-direction
	user->NumPartY         = 2; // Number of particles per cell in y-direction
	user->NumPartZ         = 2; // Number of particles per cell in z-direction

	// scaling
	jr->scal.utype = _NONE_;

	// time-stepping
	user->time_end         = 1.0;
	user->save_timesteps   = 1;
	user->CFL              = 0.5;
	user->dt_max           = 1e6; // maximum timestep
	user->dt               = 1e-3;// initial timestep

	// temperature
	user->Temp_top         = 0.0;
	user->Temp_bottom      = 1.0;
	//user->GasConstant      = 1.0;            // REMOVE

	user->Gravity          = 1.0;
	user->GravityAngle     = 0.0; // angle of gravity with z-axis (can be changed in the x-z plane)

	// input/output
	sprintf(user->OutputFile,                    "FallingBlock3D");
	sprintf(user->ParticleFilename,              "ParticlesInput.dat");
	sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");
//	user->PolyInVolSkip[0] = 0;

	// FDSTAG Canonical Default Model Setup
	user->msetup            = BLOCK;

	// write partitioning
	user->SavePartitioning  = PETSC_FALSE;

	// linear solver settings
	user->SkipStokesSolver     = PETSC_FALSE;

	// boundary conditions
	user->BC.Vy_front      = 0;
	user->BC.Vy_back       = 0;
	user->BC.Vz_bot        = 0;
	user->BC.Vz_top        = 0;
	user->BC.Vx_left       = 0;
	user->BC.Vx_right      = 0;
	user->BC.Exx           = 0;
	user->BC.Eyy           = 0;
	user->BC.InternalBound = 0; // 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.UpperBound    = 1; // 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.LowerBound    = 1; // 0-free surface, 1-free slip					, 2-no-slip
	user->BC.LeftBound     = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no slip
	user->BC.RightBound    = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.FrontBound    = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.BackBound     = 1; // 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip

	user->ParticleInput      = 1;  // 0-do not use particles to track phases; 1-do use particles to track phases ???
	user->SaveParticles      = 0;  // Save particles or not?
	user->restart            = 0;  // Are we restarting a simulation?
	user->save_breakpoints   = 10; // After how many steps do we make a breakpoint?
	user->break_point_number = 0;  // The number of the breakpoint file

	// optimization
	user->LowerViscosityCutoff  = 1.0;  // lowermost viscosity cutoff in model
	user->UpperViscosityCutoff  = 1e30; // uppermost viscosity cutoff in model
	user->InitViscosity         = 0.0;  // initial viscosity
	user->PlastViscosity        = 0.0;  // plasticity regularization viscosity
	user->DII_ref               = 0.0;  // initial guess strain rate

	user->mpi_group_id = 0;

	// phases
	jr->numPhases = 2;  // the default # of phases.
	// In case we set props from the command line, and we use a folding
	// setup combined with FDSTAG we might have to do something smarter
	// (check the setup and increase this)

	// Default Pushing BC parameters
	user->AddPushing                    = 0;
	user->Pushing.reset_pushing_coord   = 0;
	user->Pushing.theta                 = 0.0;

	user->FSSA                          =	0.0;

	// set this option to monitor actual option usage
	PetscOptionsInsertString("-options_left");
    
    // Add a few default options
    PetscOptionsInsertString("-options_left");
    
    

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InputReadFile"
PetscErrorCode InputReadFile(JacRes *jr, UserCtx *user, FILE *fp)
{
	 // parse the input file

	PetscInt found;
	double d_values[1000];
	PetscInt i_values[1000];
	PetscInt nv, i;
	char setup_name[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read number of cells
	parse_GetInt(fp, "nel_x",   &user->nel_x, &found);
	parse_GetInt(fp, "nel_y",   &user->nel_y, &found);
	parse_GetInt(fp, "nel_z",   &user->nel_z, &found);

	// only number of cells is a relevant parameter for FDSTAG
	user->nnode_x = user->nel_x+1;
	user->nnode_y = user->nel_y+1;
	user->nnode_z = user->nel_z+1;

	// Non-dimensionalization parameters
	ierr = ScalingReadFromFile(&jr->scal, fp); CHKERRQ(ierr);

	// domain parameters
	parse_GetDouble( fp, "L", &user->L, &found );
	parse_GetDouble( fp, "W", &user->W, &found );
	parse_GetDouble( fp, "H", &user->H, &found );
	parse_GetDouble( fp, "x_left",  &user->x_left,  &found );
	parse_GetDouble( fp, "y_front", &user->y_front, &found );
	parse_GetDouble( fp, "z_bot",   &user->z_bot,   &found );

	//==============
	// MESH SEGMENTS
	//==============
	ierr = ReadMeshSegDir(fp, "seg_x", user->x_left,  user->x_left  + user->W, &user->nel_x, &user->mseg_x); CHKERRQ(ierr);
	ierr = ReadMeshSegDir(fp, "seg_y", user->y_front, user->y_front + user->L, &user->nel_y, &user->mseg_y); CHKERRQ(ierr);
	ierr = ReadMeshSegDir(fp, "seg_z", user->z_bot,   user->z_bot   + user->H, &user->nel_z, &user->mseg_z); CHKERRQ(ierr);

	// read model setup
	parse_GetString(fp, "msetup", setup_name, MAX_NAME_LEN, &found);
	if(found)
	{
		if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "polygons"))   user->msetup = POLYGONS;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else if(!strcmp(setup_name, "domes"))      user->msetup = DOMES;
		else if(!strcmp(setup_name, "rotation"))   user->msetup = ROTATION;
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#ERROR! Incorrect model setup: %s", setup_name);
	}

	parse_GetString( fp, "OutputFile", user->OutputFile, PETSC_MAX_PATH_LEN, &found );
	parse_GetInt( fp,    "save_breakpoints", &user->save_breakpoints, &found );
	parse_GetInt( fp,    "restart", &user->restart, &found );
	parse_GetInt( fp,    "save_timesteps", &user->save_timesteps, &found );
	parse_GetInt( fp,    "time_end", &user->time_end, &found );
	parse_GetDouble( fp, "CFL", &user->CFL, &found );
	parse_GetDouble( fp, "dt_max", &user->dt_max, &found );
	parse_GetDouble( fp, "FSSA", &user->FSSA, &found );	// FSSA parameter

	// Particle related variables
	parse_GetInt( fp,    "ParticleInput", &user->ParticleInput, &found );
	parse_GetString( fp, "ParticleFilename", user->ParticleFilename, PETSC_MAX_PATH_LEN, &found );
	parse_GetString( fp, "TemperatureFilename", user->TemperatureFilename, PETSC_MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->TemperatureFilename, "noTemperatureFileName");
	}
	parse_GetString( fp, "TopoFilename", user->TopoFilename, PETSC_MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->TopoFilename, "noTopoFileName");
	}
	parse_GetString( fp, "LoadInitialParticlesDirectory", user->LoadInitialParticlesDirectory, PETSC_MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetString( fp, "SaveInitialParticlesDirectory", user->SaveInitialParticlesDirectory, PETSC_MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->SaveInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetInt( fp,    "SaveParticles", &user->SaveParticles, &found );

	parse_GetInt( fp,    "NumPartX", &user->NumPartX, &found );
	parse_GetInt( fp,    "NumPartY", &user->NumPartY, &found );
	parse_GetInt( fp,    "NumPartZ", &user->NumPartZ, &found );

	// boundary conditions
	parse_GetDouble( fp, "BC.Eyy", &user->BC.Eyy, &found );
	parse_GetDouble( fp, "BC.Exx", &user->BC.Exx, &found );

	parse_GetInt( fp, "BC.InternalBound", &user->BC.InternalBound, &found );
	parse_GetInt( fp, "BC.LeftBound", &user->BC.LeftBound, &found );
	parse_GetInt( fp, "BC.RightBound", &user->BC.RightBound, &found );
	parse_GetInt( fp, "BC.FrontBound", &user->BC.FrontBound, &found );
	parse_GetInt( fp, "BC.BackBound", &user->BC.BackBound, &found );
	parse_GetInt( fp, "BC.LowerBound", &user->BC.LowerBound, &found );
	parse_GetInt( fp, "BC.UpperBound", &user->BC.UpperBound, &found );
	parse_GetDouble( fp, "Temp_top", &user->Temp_top, &found );
	parse_GetDouble( fp, "Temp_bottom", &user->Temp_bottom, &found );

	// optimization
	parse_GetDouble( fp, "LowerViscosityCutoff", &user->LowerViscosityCutoff, &found );
	parse_GetDouble( fp, "UpperViscosityCutoff", &user->UpperViscosityCutoff, &found );
	parse_GetDouble( fp, "InitViscosity",        &user->InitViscosity,        &found );
	parse_GetDouble( fp, "PlastViscosity",       &user->PlastViscosity,       &found );

	// initial guess strain rate
	parse_GetDouble( fp, "DII_ref",              &user->DII_ref,        &found );

	parse_GetDouble( fp, "Gravity", &user->Gravity, &found );

	// Pushing Parameters
	parse_GetInt( fp,    "AddPushing", &user->AddPushing, &found );
	parse_GetInt( fp,    "Pushing.num_changes", &user->Pushing.num_changes, &found );
	parse_GetInt( fp,    "Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, &found );
	parse_GetDouble( fp, "Pushing.L_block", &user->Pushing.L_block, &found );
	parse_GetDouble( fp, "Pushing.W_block", &user->Pushing.W_block, &found );
	parse_GetDouble( fp, "Pushing.H_block", &user->Pushing.H_block, &found );
	parse_GetDouble( fp, "Pushing.theta", &user->Pushing.theta, &found );
	parse_GetDouble( fp, "Pushing.x_center_block", &user->Pushing.x_center_block, &found );
	parse_GetDouble( fp, "Pushing.y_center_block", &user->Pushing.y_center_block, &found );
	parse_GetDouble( fp, "Pushing.z_center_block", &user->Pushing.z_center_block, &found );
	parse_GetDoubleArray( fp, "Pushing.time",   &nv, d_values, &found );  for( i=0; i<user->Pushing.num_changes+1; i++ ) { user->Pushing.time[i]   = d_values[i];}
	parse_GetDoubleArray( fp, "Pushing.V_push", &nv, d_values, &found );  for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.V_push[i] = d_values[i];}
	parse_GetDoubleArray( fp, "Pushing.omega",  &nv, d_values, &found  ); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.omega[i]  = d_values[i];}
	parse_GetIntArray( fp, "Pushing.coord_advect",&nv, i_values, &found); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.coord_advect[i] = i_values[i];}
	parse_GetIntArray( fp, "Pushing.dir",         &nv, i_values, &found); for( i=0; i<user->Pushing.num_changes;   i++ ) { user->Pushing.dir[i]          = i_values[i];}

/*
	// Marker setting: skip certain volumes that are defined in input file
	parse_GetIntArray( fp, "PolyInVolSkip",         &nv, i_values, &found);
	if (found==1)
	{
		for( i=0; i<=i_values[0];   i++ )
		{
			user->PolyInVolSkip[i] = i_values[i];
		}
	}
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// change values from the command line
#undef __FUNCT__
#define __FUNCT__ "InputReadCommLine"
PetscErrorCode InputReadCommLine(UserCtx *user )
{
	PetscInt       i, nel_array[3], nel_input_max;
	PetscBool      found,flg;
	char           setup_name[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscOptionsGetReal(PETSC_NULL,"-W"      , &user->W      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-L"      , &user->L      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-H"      , &user->H      , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-y_front", &user->y_front, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-nel_x",   &user->nel_x,   PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_y",   &user->nel_y,   PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_z",   &user->nel_z,   PETSC_NULL);

	// alternative: specify the # of elements as -nel 8,16,32  which gives nel_x=8, nel_y=16, nel_z=32
	nel_input_max=3;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-nel", nel_array, &nel_input_max, &found); CHKERRQ(ierr);

	if (found==PETSC_TRUE) {
		user->nel_x = nel_array[0];
		user->nel_y = nel_array[1];
		user->nel_z = nel_array[2];
	}
	user->nnode_x = user->nel_x + 1;
	user->nnode_y = user->nel_y + 1;
	user->nnode_z = user->nel_z + 1;

	PetscOptionsGetInt(PETSC_NULL ,"-time_end",       &user->time_end,       PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-save_timesteps", &user->save_timesteps, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-CFL",            &user->CFL,            PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt",             &user->dt,             PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt_max",         &user->dt_max,         PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-FSSA",           &user->FSSA,           PETSC_NULL); // FSSA parameter [should be between 0-1]

	// FDSTAG Canonical Model Setup
	PetscOptionsGetString(PETSC_NULL,"-msetup", setup_name, MAX_NAME_LEN, &found);
	if(found == PETSC_TRUE)
	{	if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "polygons"))   user->msetup = POLYGONS;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else if(!strcmp(setup_name, "domes"))      user->msetup = DOMES;
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"ERROR! Incorrect model setup: %s", setup_name);
	}

	// boundary conditions
	PetscOptionsGetInt(PETSC_NULL, "-BC.InternalBound", &user->BC.InternalBound, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.UpperBound",    &user->BC.UpperBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.LowerBound",    &user->BC.LowerBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.LeftBound",     &user->BC.LeftBound,     PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.RightBound",    &user->BC.RightBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.FrontBound",    &user->BC.FrontBound,    PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-BC.BackBound",     &user->BC.BackBound,     PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vy_front",      &user->BC.Vy_front,      PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vy_back",       &user->BC.Vy_back,       PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vz_top",        &user->BC.Vz_top,        PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vz_bot",        &user->BC.Vz_bot,        PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vx_left",       &user->BC.Vx_left,       PETSC_NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL,"-BC.Vx_right",      &user->BC.Vx_right,      PETSC_NULL); // y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL,"-BC.Exx",           &user->BC.Exx,           PETSC_NULL); // Exx background strain-rate
	PetscOptionsGetReal(PETSC_NULL,"-BC.Eyy",           &user->BC.Eyy,           PETSC_NULL); // Eyy background strain-rate

	// number of markers
	PetscOptionsGetInt(PETSC_NULL,"-NumPartX",          &user->NumPartX,         PETSC_NULL); // # of tracers per cell in x-direction
	PetscOptionsGetInt(PETSC_NULL,"-NumPartY",          &user->NumPartY,         PETSC_NULL); // # of tracers per cell in y-direction
	PetscOptionsGetInt(PETSC_NULL,"-NumPartZ",          &user->NumPartZ,         PETSC_NULL); // # of tracers per cell in z-direction

	// flags
	PetscOptionsGetInt(PETSC_NULL,"-restart",          &user->restart,          PETSC_NULL); // # restart a simulation if possible?
	PetscOptionsGetInt(PETSC_NULL,"-save_breakpoints", &user->save_breakpoints, PETSC_NULL); // after how many steps do we create a breakpoint file?
	PetscOptionsGetInt(PETSC_NULL,"-SaveParticles",    &user->SaveParticles,    PETSC_NULL); // save particles to disk?

	PetscOptionsGetBool(PETSC_NULL,"-SavePartitioning",&user->SavePartitioning, PETSC_NULL);
	PetscOptionsGetBool(PETSC_NULL,"-SkipStokesSolver",&user->SkipStokesSolver, PETSC_NULL);

	// optimization
	PetscOptionsGetReal(PETSC_NULL,"-LowerViscosityCutoff", &user->LowerViscosityCutoff, PETSC_NULL); // lower viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL,"-UpperViscosityCutoff", &user->UpperViscosityCutoff, PETSC_NULL); // upper viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL,"-InitViscosity",        &user->InitViscosity,        PETSC_NULL); // upper viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL,"-PlastViscosity",       &user->PlastViscosity,       PETSC_NULL); // upper viscosity cutoff

	// initial guess strain rate
    PetscOptionsGetReal(PETSC_NULL,"-DII_ref",              &user->DII_ref,        PETSC_NULL);

	// gravity
	PetscOptionsGetReal(PETSC_NULL ,"-GravityAngle",   &user->GravityAngle,     PETSC_NULL); // Gravity angle in x-z plane

	// pushing boundary conditions related parameters
	PetscOptionsGetInt(PETSC_NULL,"-AddPushing"                 , &user->AddPushing                 , PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.num_changes"        , &user->Pushing.num_changes        , PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.L_block"           , &user->Pushing.L_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.W_block"           , &user->Pushing.W_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.H_block"           , &user->Pushing.H_block            , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.x_center_block"    , &user->Pushing.x_center_block     , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.y_center_block"    , &user->Pushing.y_center_block     , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.z_center_block"    , &user->Pushing.z_center_block     , PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.theta"             , &user->Pushing.theta              , PETSC_NULL);

	// pushing - array variables
	char matprop_opt[PETSC_MAX_PATH_LEN];
	flg = PETSC_FALSE;

	for(i=0;i<user->Pushing.num_changes;i++){
		// v_push in cm/year
		sprintf(matprop_opt,"-Pushing.V_push_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.V_push[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    V_push[%lld] = %g \n",(LLD)i,user->Pushing.V_push[i]);

		// rate of rotation deg/yr
		sprintf(matprop_opt,"-Pushing.omega_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.omega[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Omega[%lld] = %g \n",(LLD)i,user->Pushing.omega[i]);

		// advect block or not for each time segment
		sprintf(matprop_opt,"-Pushing.coord_advect_%lld",(LLD)i);
		ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.coord_advect[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Coord_advect[%lld] = %g \n",(LLD)i,user->Pushing.coord_advect[i]);

		// direction of pushing
				sprintf(matprop_opt,"-Pushing.dir_%lld",(LLD)i);
				ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.dir[i], &flg); CHKERRQ(ierr);
				if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Direction[%lld] = %g \n",(LLD)i,user->Pushing.dir[i]);
	}

	for(i=0;i<user->Pushing.num_changes+1;i++){
		// time intervals is a num_changes+1 array
		sprintf(matprop_opt,"-Pushing.time_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.time[i], &flg); CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Time[%lld] = %g \n",(LLD)i,user->Pushing.time[i]);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ReadMeshSegDir"
PetscErrorCode ReadMeshSegDir(
	FILE        *fp,
	const char  *name,
	PetscScalar  beg,
	PetscScalar  end,
	PetscInt    *tncels,
	MeshSegInp  *msi)
{
	// read mesh refinement data for a direction from the input file
	// NOTE: parameter "tncels" passes negative number of segments & returns total number of cells

	PetscInt    i, jj, arsz, found;
	PetscScalar buff[3*MaxNumMeshSegs-1];

	PetscFunctionBegin;

	// check whether segments are not specified
	if((*tncels) > 0)
	{
		msi->nsegs = 0;
		PetscFunctionReturn(0);
	}

	// set number of segments
	msi->nsegs = -(*tncels);

	// read segments
	parse_GetDoubleArray(fp, name, &arsz, buff, &found);

	if(!found) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Mesh refinement segments are not specified\n");

	if(arsz != 3*msi->nsegs-1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect number entries in the mesh refinement array\n");

	// load the data ... delimiters
	for(i = 0, jj = 0; i < msi->nsegs-1; i++, jj++) msi->delims[i] = buff[jj];

	// ... number of cells
	for(i = 0; i < msi->nsegs; i++, jj++) msi->ncells[i] = (PetscInt)buff[jj];

	// ... biases
	for(i = 0; i < msi->nsegs; i++, jj++) msi->biases[i] = buff[jj];

	// check the data ... delimiter sequence
	for(i = 1; i < msi->nsegs-1; i++)
	{
		if(msi->delims[i] <= msi->delims[i-1])
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Refinement segments are unordered/overlapping\n");
		}
	}

	// ... delimiter bounds
	for(i = 0; i < msi->nsegs-1; i++)
	{
		if(msi->delims[i] <= beg || msi->delims[i] >= end)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Refinement segments out of bound\n");
		}
	}

	// ... number of cells
	for(i = 0; i < msi->nsegs; i++)
	{
		if(msi->ncells[i] <= 0)
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Number of cells must be non-negative\n");
		}
	}

	// ... biases
	for(i = 0; i < msi->nsegs; i++)
	{
		if(!msi->biases[i])
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Bias factors must be non-zero\n");
		}
	}

	// ... compute total number of cells
	for(i = 0, (*tncels) = 0; i < msi->nsegs; i++) (*tncels) += msi->ncells[i];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
