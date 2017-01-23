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
#include "nlsolveExplicit.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PushInputReadFile"
PetscErrorCode PushInputReadFile(UserCtx *user, FILE *fp)
{
	// initialize Push parameters fromfile
	PetscInt *ls, *le;
	PetscInt i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of push parameters from file
	PetscPrintf(PETSC_COMM_WORLD, "Reading pushing parameters: \n\n");

	// clear memory
	ierr = PetscMemzero(user->Pushing, sizeof(PushParams)*(size_t)MAX_PUSH_BOX); CHKERRQ(ierr);

	// initialize ID for consistency check
	for(i = 0; i < MAX_PUSH_BOX; i++) user->Pushing[i].ID = -1;

	// allocate memory for array to store line info
	ierr = makeIntArray(&ls, NULL, MAX_PUSH_BOX); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, MAX_PUSH_BOX); CHKERRQ(ierr);

	// read number of entries
	getLineStruct(fp, ls, le, MAX_PUSH_BOX, &count_starts, &count_ends, "<PushingBlockStart>", "<PushingBlockEnd>");

	// error checking
	if(count_starts > MAX_PUSH_BOX || count_ends > MAX_PUSH_BOX)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)MAX_PUSH_BOX);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incomplete material structures! <PushingBlockStart> & <PushingBlockEnd> don't match");
	}

	// store actual number of pushing block
	user->nPush = count_starts;

	// read each individual pushing block
	for(i = 0; i < user->nPush; i++)
	{
		ierr = PushingBlockGetStruct(fp, user->nPush, user->Pushing, ls[i], le[i]); CHKERRQ(ierr);
	}

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PushingBlockGetStruct"
PetscErrorCode PushingBlockGetStruct(FILE *fp, PetscInt nPush, PushParams *Pushing, PetscInt ils, PetscInt ile)
{
	// read pushing box from file
	PushParams *p;
	PetscInt   found, ID, nv;

	PetscFunctionBegin;

	// pushing box ID
	getMatPropInt(fp, ils, ile, "PushID", &ID, &found);

	// error checking
	if(!found) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "No pushing block ID specified! ");
	if(ID > nPush-1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect pushing block numbering!");

	// get pointer to specified pushing block
	p = Pushing + ID;

	// check ID
	if(p->ID != -1)  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect Pushing block numbering!");

	// set ID
	p->ID = ID;

	// read and store pushing box parameters
	getMatPropInt      (fp, ils, ile, "num_changes",         &p->num_changes, 		  &found);
	getMatPropInt      (fp, ils, ile, "reset_pushing_coord", &p->reset_pushing_coord, &found);
	getMatPropInt      (fp, ils, ile, "ind_change", 		 &p->ind_change, 		  &found);
	getMatPropScalar   (fp, ils, ile, "theta", 				 &p->theta, 			  &found);
	getMatPropScalar   (fp, ils, ile, "L_block", 			 &p->L_block, 			  &found);
	getMatPropScalar   (fp, ils, ile, "W_block", 			 &p->W_block, 			  &found);
	getMatPropScalar   (fp, ils, ile, "H_block", 			 &p->H_block, 			  &found);
	getMatPropScalar   (fp, ils, ile, "x_center_block", 	 &p->x_center_block, 	  &found);
	getMatPropScalar   (fp, ils, ile, "y_center_block", 	 &p->y_center_block, 	  &found);
	getMatPropScalar   (fp, ils, ile, "z_center_block", 	 &p->z_center_block, 	  &found);
	getMatPropScalArray(fp, ils, ile, "V_push",         &nv,  p->V_push,			  &found);
	getMatPropScalArray(fp, ils, ile, "omega",          &nv,  p->omega,				  &found);
	getMatPropScalArray(fp, ils, ile, "time",           &nv,  p->time,  			  &found);
	getMatPropIntArray (fp, ils, ile, "coord_advect",   &nv,  p->coord_advect,		  &found);
	getMatPropIntArray (fp, ils, ile, "dir",            &nv,  p->dir,				  &found);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BezierInputReadFile"
PetscErrorCode BezierInputReadFile(UserCtx *user, FILE *fp)
{
	// initialize bezier block parameters from file
	PetscInt *ls, *le;
	PetscInt i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear memory
	ierr = PetscMemzero(user->blocks, sizeof(BCBlock)*(size_t)_max_bc_blocks_); CHKERRQ(ierr);

	// initialize ID for consistency check
	for(i = 0; i < _max_bc_blocks_; i++) user->blocks[i].ID = -1;

	// allocate memory for array to store line info
	ierr = makeIntArray(&ls, NULL, _max_bc_blocks_); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, _max_bc_blocks_); CHKERRQ(ierr);

	// read number of entries
	getLineStruct(fp, ls, le, _max_bc_blocks_, &count_starts, &count_ends, "<BezierBlockStart>", "<BezierBlockEnd>");

	// error checking
	if(count_starts > MAX_PUSH_BOX || count_ends > _max_bc_blocks_)
	{
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Too many bezier structures specified! Max allowed: %lld", (LLD)_max_bc_blocks_);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incomplete bezier structures! <BezierBlockStart> & <BezierBlockEnd> don't match");
	}

	// store actual number of bezier block
	user->nblo = count_starts;

	// read each individual bezier block
	for(i = 0; i< user->nblo; i++)
	{
		ierr = BezierBlockGetStruct(fp, user->nblo, user->blocks, ls[i], le[i]); CHKERRQ(ierr);
	}

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BezierBlockGetStruct"
PetscErrorCode BezierBlockGetStruct(FILE *fp, PetscInt nblo, BCBlock *blocks, PetscInt ils, PetscInt ile)
{
	// read bezier block from file
	BCBlock  *b;
	PetscInt found, ID, nv;

	PetscFunctionBegin;

	// bezier box ID
	getMatPropInt(fp, ils, ile, "BezierID", &ID, &found);

	// error checking
	if(!found) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "No bezier block ID specified! ");
	if(ID > nblo-1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect bezier block numbering!");

	// get pointer to specified bezier block
	b = blocks + ID;

	// check ID
	if(b->ID != -1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect Bezier block numbering!");

	// set ID
	b->ID = ID;

	// read and store bezier block parameters & error checking
	getMatPropInt(fp, ils, ile, "npath", &b->npath, &found);
	if(!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "No npath specified! ");
	if(b->npath < 1 || b->npath > _max_path_points_) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER, "npath should be range: [%d - %d]\n", 1, _max_path_points_);

	getMatPropScalArray(fp,ils,ile,"theta",&nv, b->theta, &found);
	if(nv != b->npath) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "wrong number of theta specified! ");

	getMatPropScalArray(fp,ils,ile,"time",&nv, b->time, &found);
	if(nv != b->npath) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "wrong number of time specified! ");

	getMatPropScalArray(fp,ils,ile,"path",&nv, b->path, &found);
	if(nv != 6*b->npath -4) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "wrong number of path specified! ");

	getMatPropInt(fp,ils,ile,"npoly",&b->npoly, &found);
	if(!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "number of polygon not specified! ");
	if(b->npoly < 1 || b->npoly > _max_poly_points_) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER, "npoly should be range: [%d - %d]\n", 1, _max_poly_points_);

	getMatPropScalArray(fp,ils,ile,"poly",&nv, b->poly, &found);
	if(nv != 2*b->npoly) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "wrong number of poly specified! ");

	getMatPropScalar(fp, ils, ile, "bot", &b->bot, &found);
	getMatPropScalar(fp, ils, ile, "top", &b->top, &found);

	PetscFunctionReturn(0);
}
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

	// Explicit solver (wave propagation)
	jr->ExplicitSolver = user->ExplicitSolver;

	// Force damping factor used to absorb oscillation energy
	jr->ForceDampingFactor  = user->ForceDampingFactor;

	jr->DensityFactor  = user->DensityFactor;
	//jr->OutputStressFile  = user->OutputStressFile;

	jr->SeismicSource  = user->SeismicSource;
	if (jr->SeismicSource == PETSC_TRUE) {
		jr->SourceParams = user->SourceParams;
	}

	jr->SeismicStation= user->SeismicStation;
	if (jr->SeismicStation == PETSC_TRUE) {
		jr->StationParams = user->StationParams;
	}

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

	PetscOptionsInsertString(NULL, "-options_left");

	// Add a few default options
	PetscOptionsInsertString(NULL, "-options_left");

	// Resolve SuperLU_DIST repetitive factorization issue (temporary ad hoc solution)
	PetscOptionsInsertString(NULL, "-mat_superlu_dist_fact SamePattern_SameRowPerm");

    // Explicit solver (wave propagation)
    user->ExplicitSolver = PETSC_FALSE;
    user->SeismicSource	 = PETSC_FALSE;
    user->DensityFactor = 1.0;

    user->SourceParams.x = -1;
    user->SourceParams.y = -1;
    user->SourceParams.z = -1;
	user->SourceParams.i = -1;
	user->SourceParams.j = -1;
	user->SourceParams.k = -1;

    user->StationParams.x = -1;
    user->StationParams.y = -1;
    user->StationParams.z = -1;
    user->StationParams.i = -1;
	user->StationParams.j = -1;
	user->StationParams.k = -1;

    user->AbsBoundaries = PETSC_TRUE;
    user->AB.NxL = 20;
    user->AB.NxR = 20;
    user->AB.NyL = 20;
    user->AB.NyR = 20;
    user->AB.NzL = 20;
    user->AB.NzR = 20;

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
	PetscInt nv, i, ab, source, station;
	char setup_name[MAX_NAME_LEN], source_type_name[MAX_NAME_LEN];

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
		else if(!strcmp(setup_name, "homo"))       user->msetup = HOMO;		// Created to try wave propagation
		else if(!strcmp(setup_name, "layer"))      user->msetup = LAYER;	// Created to try wave propagation
		else if(!strcmp(setup_name, "heterogeneous"))       user->msetup = HETEROGENEOUS;	// Created to try wave propagation
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

	// Force damping factor used to absorb oscillation energy
	parse_GetDouble(fp, "ForceDampingFactor",&user->ForceDampingFactor, &found);
	if (found!=PETSC_TRUE ) user->ForceDampingFactor=0.8;

	// Scaling density factor
	parse_GetDouble(fp, "density_factor",&user->DensityFactor, &found);
	if (found!=PETSC_TRUE ) user->DensityFactor=1.0;
	parse_GetString(fp,"OutputStressFile", &user->OutputStressFile, MAX_PATH_LEN, &found);
	if (found !=PETSC_TRUE) sprintf(&user->OutputStressFile, "stress_file");

	// Absorbing boundaries
	//if (user->ExplicitSolver == PETSC_TRUE)
	//{
		parse_GetInt(fp, "abs_boundaries",&ab, &found);
		if (found==PETSC_TRUE && ab==0)
		{
			user->AbsBoundaries=PETSC_FALSE;
			user->AB.NxL = 0;
			user->AB.NxR = 0;
			user->AB.NyL = 0;
			user->AB.NyR = 0;
			user->AB.NzL = 0;
			user->AB.NzR = 0;
		}else{	// Default value is true
			// Number of absorbing boundaries
			parse_GetInt( fp, "AB.NxL", &user->AB.NxL, &found );
			if ( user->AB.NxL >= user->nel_x/2) user->AB.NxL = 0;
			parse_GetInt( fp, "AB.NxR", &user->AB.NxR, &found );
			if ( user->AB.NxR >= user->nel_x/2) user->AB.NxR = 0;
			parse_GetInt( fp, "AB.NyL", &user->AB.NyL, &found );
			if ( user->AB.NyL >= user->nel_y/2) user->AB.NyL = 0;
			parse_GetInt( fp, "AB.NyR", &user->AB.NyR, &found );
			if ( user->AB.NyR >= user->nel_y/2) user->AB.NyR = 0;
			parse_GetInt( fp, "AB.NzL", &user->AB.NzL, &found );
			if ( user->AB.NzL >= user->nel_z/2) user->AB.NzL = 0;
			parse_GetInt( fp, "AB.NzR", &user->AB.NzR, &found );
			if ( user->AB.NzR >= user->nel_z/2) user->AB.NzR = 0;
		}
	//}


	// Seismic source
	parse_GetInt(fp, "seismic_source",&source, &found);
	if (found==PETSC_TRUE && source==1)
	{
		user->SeismicSource=PETSC_TRUE;
		parse_GetString(fp,"source_type", source_type_name, MAX_NAME_LEN, &found);
		if(found == PETSC_TRUE)
		{
			if      (!strcmp(source_type_name, "point"))   user->SourceParams.source_type = POINT;
			else if (!strcmp(source_type_name, "plane"))   user->SourceParams.source_type = PLANE;
			else if (!strcmp(source_type_name, "uniaxial_compression"))   user->SourceParams.source_type = COMPRES;
			else if (!strcmp(source_type_name, "moment_tensor"))   user->SourceParams.source_type = MOMENT;
			else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"ERROR! Incorrect source type: %s", source_type_name);

			// Source coordinates
			if (user->SourceParams.source_type == POINT || user->SourceParams.source_type == MOMENT) {
				parse_GetDouble( fp, "x_source", &user->SourceParams.x, &found );
				parse_GetDouble( fp, "y_source", &user->SourceParams.y, &found );
				parse_GetDouble( fp, "z_source", &user->SourceParams.z, &found );

				if ( (user->SourceParams.x < 0.0) || (user->SourceParams.x > user->W) || (user->SourceParams.y < 0.0) || (user->SourceParams.y > user->L) || (user->SourceParams.z < 0.0) || (user->SourceParams.z > user->H))
				{
					SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect source coordinates\n");
				}
			}
			// Source function
			parse_GetDouble( fp, "amplitude", &user->SourceParams.amplitude, &found );
			parse_GetDouble( fp, "alfa", &user->SourceParams.alfa, &found );
			parse_GetDouble( fp, "t0", &user->SourceParams.t0, &found );
			parse_GetDouble( fp, "frequency", &user->SourceParams.frequency, &found );

			// Moment tensor source
			if (user->SourceParams.source_type == MOMENT) {
				parse_GetDouble( fp, "M0", &user->SourceParams.moment_tensor.M0, &found );
				parse_GetDouble( fp, "Mxx", &user->SourceParams.moment_tensor.Mxx, &found );
				parse_GetDouble( fp, "Myy", &user->SourceParams.moment_tensor.Myy, &found );
				parse_GetDouble( fp, "Mzz", &user->SourceParams.moment_tensor.Mzz, &found );
				parse_GetDouble( fp, "Mxy", &user->SourceParams.moment_tensor.Mxy, &found );
				parse_GetDouble( fp, "Mxz", &user->SourceParams.moment_tensor.Mxz, &found );
				parse_GetDouble( fp, "Myz", &user->SourceParams.moment_tensor.Myz, &found );
			}

			/*// Initialize other fields
			user->SourceParams.i = -1;
			user->SourceParams.j = -1;
			user->SourceParams.k = -1;
			user->SourceParams.xrank = -1;
			user->SourceParams.yrank = -1;
			user->SourceParams.zrank = -1;*/
		}
	}else
	{
		user->SeismicSource=PETSC_FALSE;
	}

	// Seismic station
	parse_GetInt(fp, "seismic_station",&station, &found);
	if (found==PETSC_TRUE && station==1)
	{
		user->SeismicStation=PETSC_TRUE;
		parse_GetDouble( fp, "x_rec", &user->StationParams.x, &found );
		parse_GetDouble( fp, "y_rec", &user->StationParams.y, &found );
		parse_GetDouble( fp, "z_rec", &user->StationParams.z, &found );
		if ( user->StationParams.x >= user->W) user->StationParams.x=user->W/2;
		if ( user->StationParams.y >= user->L) user->StationParams.y=user->L/2;
		if ( user->StationParams.z >= user->H) user->StationParams.z=user->H;
		//parse_GetString(fp,"output_file_name", user->StationParams.output_file_name, MAX_PATH_LEN, &found);
		//if (found ==PETSC_FALSE) sprintf(user->StationParams.output_file_name, "seismogram");
	}else
	{
		user->SeismicStation=PETSC_FALSE;
	}



	// Seismic station coordinates
	//parse_GetDouble( fp, "x_rec", &user->Station.x, &found );
	//parse_GetDouble( fp, "y_rec", &user->Station.y, &found );
	//parse_GetDouble( fp, "z_rec", &user->Station.z, &found );
	//if (&user->Station.x < 0 || &user->Station.x > &user->W || &user->Station.y < 0 || &user->Station.y > &user->L || &user->Station.z < 0 || &user->Station.z > &user->H)
	//{
		//SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect station coordinates\n");
	//}


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
		else if(!strcmp(setup_name, "homo"))       user->msetup = HOMO;				// Created to try wave propagation
		else if(!strcmp(setup_name, "layer"))      user->msetup = LAYER;			// Created to try wave propagation
		else if(!strcmp(setup_name, "heterogeneous")) user->msetup = HETEROGENEOUS;	// Created to try wave propagation
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

	PetscOptionsGetBool(NULL, NULL,"-ExplicitSolver",&user->ExplicitSolver, NULL);

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
