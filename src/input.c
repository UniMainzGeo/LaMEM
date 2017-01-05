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
#include "scaling.h"
#include "tssolve.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "JacRes.h"
#include "input.h"
#include "matProps.h"

/*
//---------------------------------------------------------------------------
// set default code parameters and read input file, if required
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitCode"
PetscErrorCode FDSTAGInitCode(JacRes *jr, UserCtx *user, ModParam *iop)
{
	FILE      *fp;
	char      *all_options;
	char      ParamFile[MAX_PATH_LEN];
	PetscBool InputParamFile;

	PetscMPIInt  size;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED
	PetscOptionsGetAll(NULL, &all_options ); // copy all command line args

	// set default values for parameters
	ierr = InputSetDefaultValues(jr, user); CHKERRQ(ierr);

	// check whether input file is specified
	ierr = PetscOptionsGetString(NULL, NULL, "-ParamFile", ParamFile, MAX_PATH_LEN, &InputParamFile); CHKERRQ(ierr);

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

			// read pushing parameter from file
			ierr = PushInputReadFile(user,fp); CHKERRQ(ierr);

			// read bezier parameter from file
			ierr = BezierInputReadFile(user,fp); CHKERRQ(ierr);

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

	// set pressure limit for plasticity
	jr->matLim.presLimAct  = PETSC_TRUE;

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
//	user->Pushing.reset_pushing_coord   = 0;
//	user->Pushing.theta                 = 0.0;

	// default bezier parameters
	user->AddBezier 					=	0;

	user->FSSA                          =	0.0;

	// set this option to monitor actual option usage
	PetscOptionsInsertString(NULL, "-options_left");

	// Add a few default options
	PetscOptionsInsertString(NULL, "-options_left");

	// Resolve SuperLU_DIST repetitive factorization issue (temporary ad hoc solution)
	PetscOptionsInsertString(NULL, "-mat_superlu_dist_fact SamePattern_SameRowPerm");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InputReadFile"
PetscErrorCode InputReadFile(JacRes *jr, UserCtx *user, FILE *fp)
{
	// parse the input file

	PetscInt found;
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

	parse_GetString( fp, "OutputFile", user->OutputFile, MAX_PATH_LEN, &found );
	parse_GetInt( fp,    "save_breakpoints", &user->save_breakpoints, &found );
	parse_GetInt( fp,    "restart", &user->restart, &found );
	parse_GetInt( fp,    "save_timesteps", &user->save_timesteps, &found );
	parse_GetInt( fp,    "time_end", &user->time_end, &found );
	parse_GetDouble( fp, "CFL", &user->CFL, &found );
	parse_GetDouble( fp, "dt_max", &user->dt_max, &found );
	parse_GetDouble( fp, "FSSA", &user->FSSA, &found );	// FSSA parameter

	// Particle related variables
	parse_GetInt( fp,    "ParticleInput", &user->ParticleInput, &found );
	parse_GetString( fp, "ParticleFilename", user->ParticleFilename, MAX_PATH_LEN, &found );
	parse_GetString( fp, "TemperatureFilename", user->TemperatureFilename, MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->TemperatureFilename, "noTemperatureFileName");
	}
	parse_GetString( fp, "TopoFilename", user->TopoFilename, MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->TopoFilename, "noTopoFileName");
	}
	parse_GetString( fp, "LoadInitialParticlesDirectory", user->LoadInitialParticlesDirectory, MAX_PATH_LEN, &found );
	if (!found){
		sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetString( fp, "SaveInitialParticlesDirectory", user->SaveInitialParticlesDirectory, MAX_PATH_LEN, &found );
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

	// bezier flag
	parse_GetInt( fp,    "AddBezier",  &user->AddBezier,  &found );

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// change values from the command line
#undef __FUNCT__
#define __FUNCT__ "InputReadCommLine"
PetscErrorCode InputReadCommLine(UserCtx *user )
{
	PetscInt       nel_array[3], nel_input_max;
	PetscBool      found;
	char           setup_name[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscOptionsGetReal(NULL, NULL,"-W"      , &user->W      , NULL);
	PetscOptionsGetReal(NULL, NULL,"-L"      , &user->L      , NULL);
	PetscOptionsGetReal(NULL, NULL,"-H"      , &user->H      , NULL);
	PetscOptionsGetReal(NULL, NULL,"-y_front", &user->y_front, NULL);

	PetscOptionsGetInt(NULL, NULL ,"-nel_x",   &user->nel_x,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-nel_y",   &user->nel_y,   NULL);
	PetscOptionsGetInt(NULL, NULL ,"-nel_z",   &user->nel_z,   NULL);

	// alternative: specify the # of elements as -nel 8,16,32  which gives nel_x=8, nel_y=16, nel_z=32
	nel_input_max=3;
	ierr = PetscOptionsGetIntArray(NULL, NULL,"-nel", nel_array, &nel_input_max, &found); CHKERRQ(ierr);

	if (found==PETSC_TRUE) {
		user->nel_x = nel_array[0];
		user->nel_y = nel_array[1];
		user->nel_z = nel_array[2];
	}
	user->nnode_x = user->nel_x + 1;
	user->nnode_y = user->nel_y + 1;
	user->nnode_z = user->nel_z + 1;

	PetscOptionsGetInt (NULL, NULL ,"-time_end",       &user->time_end,       NULL);
	PetscOptionsGetInt (NULL, NULL ,"-save_timesteps", &user->save_timesteps, NULL);
	PetscOptionsGetReal(NULL, NULL,"-CFL",            &user->CFL,            NULL);
	PetscOptionsGetReal(NULL, NULL,"-dt",             &user->dt,             NULL);
	PetscOptionsGetReal(NULL, NULL,"-dt_max",         &user->dt_max,         NULL);
	PetscOptionsGetReal(NULL, NULL,"-FSSA",           &user->FSSA,           NULL); // FSSA parameter [should be between 0-1]

	// FDSTAG Canonical Model Setup
	PetscOptionsGetString(NULL, NULL,"-msetup", setup_name, MAX_NAME_LEN, &found);
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
	PetscOptionsGetInt (NULL, NULL, "-BC.InternalBound", &user->BC.InternalBound, NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.UpperBound",    &user->BC.UpperBound,    NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.LowerBound",    &user->BC.LowerBound,    NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.LeftBound",     &user->BC.LeftBound,     NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.RightBound",    &user->BC.RightBound,    NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.FrontBound",    &user->BC.FrontBound,    NULL);
	PetscOptionsGetInt (NULL, NULL, "-BC.BackBound",     &user->BC.BackBound,     NULL);
	PetscOptionsGetReal(NULL, NULL,"-BC.Vy_front",      &user->BC.Vy_front,      NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(NULL, NULL,"-BC.Vy_back",       &user->BC.Vy_back,       NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(NULL, NULL,"-BC.Vz_top",        &user->BC.Vz_top,        NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(NULL, NULL,"-BC.Vz_bot",        &user->BC.Vz_bot,        NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(NULL, NULL,"-BC.Vx_left",       &user->BC.Vx_left,       NULL); // y-velocity @ front boundary
	PetscOptionsGetReal(NULL, NULL,"-BC.Vx_right",      &user->BC.Vx_right,      NULL); // y-velocity @ front boundary

	PetscOptionsGetReal(NULL, NULL,"-BC.Exx",           &user->BC.Exx,           NULL); // Exx background strain-rate
	PetscOptionsGetReal(NULL, NULL,"-BC.Eyy",           &user->BC.Eyy,           NULL); // Eyy background strain-rate

	// number of markers
	PetscOptionsGetInt(NULL, NULL,"-NumPartX",          &user->NumPartX,         NULL); // # of tracers per cell in x-direction
	PetscOptionsGetInt(NULL, NULL,"-NumPartY",          &user->NumPartY,         NULL); // # of tracers per cell in y-direction
	PetscOptionsGetInt(NULL, NULL,"-NumPartZ",          &user->NumPartZ,         NULL); // # of tracers per cell in z-direction

	// flags
	PetscOptionsGetInt(NULL, NULL,"-restart",          &user->restart,          NULL); // # restart a simulation if possible?
	PetscOptionsGetInt(NULL, NULL,"-save_breakpoints", &user->save_breakpoints, NULL); // after how many steps do we create a breakpoint file?
	PetscOptionsGetInt(NULL, NULL,"-SaveParticles",    &user->SaveParticles,    NULL); // save particles to disk?

	PetscOptionsGetBool(NULL, NULL,"-SavePartitioning",&user->SavePartitioning, NULL);
	PetscOptionsGetBool(NULL, NULL,"-SkipStokesSolver",&user->SkipStokesSolver, NULL);

	// optimization
	PetscOptionsGetReal(NULL, NULL,"-LowerViscosityCutoff", &user->LowerViscosityCutoff, NULL); // lower viscosity cutoff
	PetscOptionsGetReal(NULL, NULL,"-UpperViscosityCutoff", &user->UpperViscosityCutoff, NULL); // upper viscosity cutoff
	PetscOptionsGetReal(NULL, NULL,"-InitViscosity",        &user->InitViscosity,        NULL); // upper viscosity cutoff
	PetscOptionsGetReal(NULL, NULL,"-PlastViscosity",       &user->PlastViscosity,       NULL); // upper viscosity cutoff

	// initial guess strain rate
    PetscOptionsGetReal(NULL, NULL,"-DII_ref",              &user->DII_ref,        NULL);

	// gravity
	PetscOptionsGetReal(NULL, NULL ,"-GravityAngle",   &user->GravityAngle,     NULL); // Gravity angle in x-z plane

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
*/
