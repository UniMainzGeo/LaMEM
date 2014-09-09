/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

Breakpoint.c: Contains the following subroutines:

SaveBreakPoint					-	Save a breakpoint file, for later restarting of the code.
LoadBreakPoint					-	Loads a breakpoint file to restart the code.

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
 */


#include "LaMEM.h"
#include "Breakpoint.h"
#include "Mesh.h"


/*==========================================================================================================*/
/* Create a breakpoint file, which can be employed to restart the simulation 								*/
#undef __FUNCT__
#define __FUNCT__ "SaveBreakPoint"
PetscErrorCode SaveBreakPoint( UserContext *user, DM da_nodes, DM da_pres, Vec Velocity,
		Vec Temp, Vec Pressure, PetscInt itime, PetscInt FileNumber )
{
	PetscMPIInt         rank, size;
	PetscErrorCode 		ierr;
	PetscInt			mx,my,mz;
	PetscInt			xs,ys,zs,xm,ym,zm;
	PetscInt			xs_p,ys_p,zs_p,xm_p,ym_p,zm_p;
	Vec					information, Particle_Vec, coord, TimeDependent_Data;
	char				SaveFileName_VCoord[PETSC_MAX_PATH_LEN], SaveInformationFileName[PETSC_MAX_PATH_LEN], SaveFileNameParticles[PETSC_MAX_PATH_LEN];
	char				SaveFileName_P[PETSC_MAX_PATH_LEN], SaveFileName_T[PETSC_MAX_PATH_LEN], SaveFileName_Surface[PETSC_MAX_PATH_LEN], SaveFileName_ErosionSurface[PETSC_MAX_PATH_LEN];
	PetscViewer			view_out;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	/* Two types of files are created here:
	 * 	(1) One file that contains global data such as velocity, pressure and material properties
	 *  (2) One file with particle data - each PROC saves it's own file here.
	 */


	/* (1) Save velocity, pressure and material properties at integration points all in a single global file*/
	sprintf(SaveFileName_VCoord,	"%s%lld.breakpoint",		"BreakPoint_Global_VCoord_"		, (LLD)FileNumber	);  // construct the filename with global vectors of velocity & coordinates [both with blocksize 3]
	sprintf(SaveFileName_P,			"%s%lld.breakpoint",		"BreakPoint_Global_P_"		, (LLD)FileNumber		);  // construct the filename with global vectors of pressure [can have blocksize 1 or 4]
	sprintf(SaveFileName_T,			"%s%lld.breakpoint",		"BreakPoint_Global_T_"		, (LLD)FileNumber		);  // construct the filename with global vectors of temperature [blocksize 1]
	sprintf(SaveInformationFileName,"%s%lld.breakpoint",		"BreakPoint_Global_Information_"	, (LLD)FileNumber		);  // construct the filename that contains info (local for every proc)
	sprintf(SaveFileNameParticles  ,"%s%lld.%lld.breakpoint",	"BreakPoint_Particles_"				, (LLD)FileNumber,(LLD)rank	); 	// construct the filename with Particles
	sprintf(SaveFileName_Surface,	"%s%lld.breakpoint",		"BreakPoint_Global_InternalFreeSurface_"	, (LLD)FileNumber		);  // construct the filename with that contains the internal free surface
	sprintf(SaveFileName_ErosionSurface,	"%s%lld.breakpoint",		"BreakPoint_Global_ErosionSurface_"	, (LLD)FileNumber		);  // construct the filename with that contains the internal free surface


	ierr = DMDAGetInfo(da_nodes, 0, &mx,    &my, &mz, 0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da_nodes,  	&xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_pres,      	&xs_p, &ys_p, &zs_p, &xm_p, &ym_p, &zm_p);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_nodes,		&coord); CHKERRQ(ierr);		/* Get global coordinate vector */


	/* Create a vector with various types of info */
	ierr = VecCreate(PETSC_COMM_SELF,&information); CHKERRQ(ierr);
	ierr = VecSetSizes(information,PETSC_DECIDE,40); CHKERRQ(ierr);
	ierr = VecSetFromOptions(information); CHKERRQ(ierr);
	ierr = VecSetValue(information,0 ,	(PetscScalar)(mx), 						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,1 ,	(PetscScalar)(my), 						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,2 ,	(PetscScalar)(mz), 						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,3 ,	(PetscScalar)(size), 						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,4,	(PetscScalar)(xs)	,						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,5,	(PetscScalar)(ys)	,						INSERT_VALUES); CHKERRQ(ierr);

	ierr = VecSetValue(information,6,	(PetscScalar)(zs)	,						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,7,	(PetscScalar)(xm)	,						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,8,	(PetscScalar)(ym)	,						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,9,	(PetscScalar)(zm)	,						INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,10,	(PetscScalar)(itime),						INSERT_VALUES);	 CHKERRQ(ierr);	// step that is being saved
	ierr = VecSetValue(information,11,	(PetscScalar)(user->time),				INSERT_VALUES);	 CHKERRQ(ierr);	// time
	ierr = VecSetValue(information,12,	(PetscScalar)(user->dt),					INSERT_VALUES);	 CHKERRQ(ierr);	// current dt
	ierr = VecSetValue(information,13,	(PetscScalar)(user->break_point_number),	INSERT_VALUES);	 CHKERRQ(ierr);	// the # of the breakpoint files that were written already
	ierr = VecSetValue(information,14,	(PetscScalar)(user->ErosionParameters.HorizontalFreeSurfaceHeight),	INSERT_VALUES);	 CHKERRQ(ierr);	// the # of the breakpoint files that were written already

	/* Store W,L,H,x_left,y_front,z_bot */
	ierr = VecSetValue(information,15,	user->x_left,	INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,16,	user->y_front,	INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,17,	user->z_bot,	INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,18,	user->W,		INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,19,	user->L,		INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,20,	user->H,		INSERT_VALUES);	 CHKERRQ(ierr);

	/* Store Pushing coordinates */
	ierr = VecSetValue(information,21,	user->Pushing.x_center_block,	INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,22,	user->Pushing.y_center_block,	INSERT_VALUES);	 CHKERRQ(ierr);
	ierr = VecSetValue(information,23,	user->Pushing.z_center_block,	INSERT_VALUES);	 CHKERRQ(ierr);


	ierr = VecAssemblyBegin(information);  CHKERRQ(ierr);
	ierr = VecAssemblyEnd(information); CHKERRQ(ierr);

	// ad-hoc solution; should be replaced
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,
			1,
			NumTimeDepData*user->time_end,
			(PetscScalar*)(user->TimeDependentData),
			&TimeDependent_Data); CHKERRQ(ierr);


	ierr = VecAssemblyBegin(TimeDependent_Data); CHKERRQ(ierr); 	
	ierr = VecAssemblyEnd(TimeDependent_Data); CHKERRQ(ierr);


	/* Write the actual file */
	if (rank==0){

		/* Write file with information */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveInformationFileName,FILE_MODE_WRITE, &view_out); CHKERRQ(ierr);
		ierr = VecView(information, 		view_out); 				CHKERRQ(ierr);
		ierr = VecView(TimeDependent_Data, 	view_out); 				CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out); 						CHKERRQ(ierr);

		if(user->ErosionParameters.ErosionModel == 2){
			// Save erosion surface information
			ierr = PetscViewerCreate(PETSC_COMM_SELF,&view_out);				CHKERRQ(ierr);
			ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
			ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
			ierr = PetscViewerFileSetName(view_out,SaveFileName_ErosionSurface);		CHKERRQ(ierr);
			ierr = VecView(user->ErosionParameters.FE_ErosionCode.ErosionSurface, view_out); 	CHKERRQ(ierr);		// erosion surface

			ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
		}


	}


	/* Write the file with global vectors */
	ierr = PetscViewerCreate(PETSC_COMM_WORLD,&view_out);				CHKERRQ(ierr);
	ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
	//ierr = PetscViewerBinarySkipInfo(view_out);							CHKERRQ(ierr);			// important, since we add Vec with 3 dof's as well as Vec's with 1 dof to the file!
	ierr = PetscViewerFileSetName(view_out,SaveFileName_VCoord); 		CHKERRQ(ierr);


	ierr = VecView(Velocity, 		view_out); CHKERRQ(ierr);	// blocksize 3
	ierr = VecView(coord, 			view_out); CHKERRQ(ierr);	// blocksize 3

	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);




	// Save pressure
	ierr = PetscViewerCreate(PETSC_COMM_WORLD,&view_out);				CHKERRQ(ierr);
	ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(view_out,SaveFileName_P);				CHKERRQ(ierr);
	ierr = VecView(Pressure, 			view_out); CHKERRQ(ierr);	// blocksize 4 or 1
	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);


	// Save temperature & upper/lower surface of model
	ierr = PetscViewerCreate(PETSC_COMM_WORLD,&view_out);				CHKERRQ(ierr);
	ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(view_out,SaveFileName_T);				CHKERRQ(ierr);

#ifdef TEMPERATURE
	ierr = VecView(Temp, 			view_out); CHKERRQ(ierr);
#endif


	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);



	// Save internal free surface information
	ierr = PetscViewerCreate(PETSC_COMM_WORLD,&view_out);				CHKERRQ(ierr);
	ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(view_out,SaveFileName_Surface);		CHKERRQ(ierr);

	ierr = VecView(user->BottomTopography, view_out); 	CHKERRQ(ierr);		// surface & bottom topography
	ierr = VecView(user->SurfaceTopography, view_out); 	CHKERRQ(ierr);

	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);


	/* (2) Save particle data of each PROC in separate files */
#ifdef PARTICLES
	/* Put all particles (and corresponding data such as T, P, location etc.) in a single vector */

	// ad-hoc solution should be removed
	// ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,user->num_particle_local*particle_props,(PetscScalar*)(user->ParticlesLocal),&Particle_Vec); CHKERRQ(ierr);

	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,particle_props*user->num_particle_local,(PetscScalar*)(user->ParticlesLocal),&Particle_Vec); CHKERRQ(ierr);


	ierr = VecAssemblyBegin(Particle_Vec); 	VecAssemblyEnd(Particle_Vec); CHKERRQ(ierr);


	/* write file */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileNameParticles,FILE_MODE_WRITE, &view_out); CHKERRQ(ierr);
	ierr = VecView(Particle_Vec, 		view_out); 				CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out); 						CHKERRQ(ierr);
#endif


	/* Print info to the screen */
	PetscPrintf(PETSC_COMM_WORLD,"  \n");
	PetscPrintf(PETSC_COMM_WORLD,"**********************  \n");
	PetscPrintf(PETSC_COMM_WORLD," Created a global breakpoint file %s \n",SaveFileName_VCoord);
	PetscPrintf(PETSC_COMM_WORLD,"  # particles = %lld | rank[%lld] \n", (LLD)(user->num_particle_local), (LLD)rank);

	PetscPrintf(PETSC_COMM_WORLD,"**********************  \n");
	PetscPrintf(PETSC_COMM_WORLD,"  \n");

	/* Clean up */
	ierr = VecDestroy(&information); 	 	CHKERRQ(ierr);
	ierr = VecDestroy(&TimeDependent_Data); 	CHKERRQ(ierr);
#ifdef PARTICLES
	ierr = VecDestroy(&Particle_Vec); 		CHKERRQ(ierr);
#endif


	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Load a breakpoint file, in order to restart a simulation */
#undef __FUNCT__
#define __FUNCT__ "LoadBreakPoint"
PetscErrorCode LoadBreakPoint(LaMEMVelPressureDA C, UserContext *user, DM da_nodes, Vec Velocity,
		Vec Temp, Vec Pressure, PetscInt FileNumber )
{
	PetscMPIInt					rank, size;
	PetscErrorCode				ierr;
	Vec							information,  coord, Particle_Vec, TimeDependent_Data;
	PetscInt					nx,ny,nz, ind[1];
	PetscInt					itime,i;
	PetscInt					ipart, numParticles;
	PetscScalar					n_read[1], *information_vec;
	PetscViewer					view_in,view_in1;
	char						LoadFileName_VCoord[PETSC_MAX_PATH_LEN], LoadInformationFileName[PETSC_MAX_PATH_LEN], LoadFileNameParticles[PETSC_MAX_PATH_LEN];
	char						LoadFileName_P[PETSC_MAX_PATH_LEN], LoadFileName_T[PETSC_MAX_PATH_LEN], LoadFileName_Surface[PETSC_MAX_PATH_LEN], LoadFileName_ErosionSurface[PETSC_MAX_PATH_LEN];
	Particles					*Particle_array;
	GlobalTimeDependentData		*TimeDependent_array;
	DAVPElementType 			element_type;
	FILE *fd;
	PetscBool create_break_point;

	ierr = DMDAGetInfo(da_nodes,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions

	element_type = C->type;

	/* Two types of files need to be loaded here:
	 * 	(1) One file that contains global data such as velocity, pressure and material properties
	 *  (2) One file with particle data - each PROC saves it's own file here.
	 */

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	/* Read the input coordinates on cpu 0 */

	/* Load the file */
	FileNumber = user->fileno;
	sprintf(LoadFileName_VCoord,	"%s%lld.breakpoint",		"BreakPoint_Global_VCoord_"		, (LLD)FileNumber);  		// construct the filename
	sprintf(LoadFileName_P,			"%s%lld.breakpoint",		"BreakPoint_Global_P_"		, (LLD)FileNumber);  		// construct the filename
	sprintf(LoadFileName_T,			"%s%lld.breakpoint",		"BreakPoint_Global_T_"		, (LLD)FileNumber);  		// construct the filename

	sprintf(LoadInformationFileName,"%s%lld.breakpoint",		"BreakPoint_Global_Information_"	, (LLD)FileNumber);  		// construct the filename that contains info
	sprintf(LoadFileNameParticles  ,"%s%lld.%lld.breakpoint",	"BreakPoint_Particles_"				, (LLD)FileNumber,(LLD)rank	); 	// construct the filename with Particles
	sprintf(LoadFileName_Surface,	"%s%lld.breakpoint",		"BreakPoint_Global_InternalFreeSurface_"	, (LLD)FileNumber		);  // construct the filename with that contains the internal free surface
	sprintf(LoadFileName_ErosionSurface,	"%s%lld.breakpoint",		"BreakPoint_Global_ErosionSurface_"	, (LLD)FileNumber		);  // construct the filename with that contains the internal free surface


	create_break_point = PETSC_FALSE;
	fd = fopen(LoadInformationFileName,"r");						 // check whether file exists
	if( fd != NULL ) {
		create_break_point = PETSC_TRUE;
		fclose(fd);
	}

	if (create_break_point==PETSC_FALSE){
		/* No restart file detected */
		PetscPrintf(PETSC_COMM_WORLD," No breakpoint detected; starting new simulation \n");
	}
	else {
		/* Breakpoint file detected - load it*/

		PetscPrintf(PETSC_COMM_WORLD,"  \n");
		PetscPrintf(PETSC_COMM_WORLD,"**********************  \n");
		PetscPrintf(PETSC_COMM_WORLD," Restarting simulation from breakpoint file %s \n",LoadFileName_VCoord);


		/* (1) Load data from the single global breakpoint file */

		/* Load information vector */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadInformationFileName,FILE_MODE_READ, &view_in1); CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&information);CHKERRQ(ierr);
		ierr = VecLoad(information,view_in1); 		CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&TimeDependent_Data);CHKERRQ(ierr);
		ierr = VecLoad(TimeDependent_Data,view_in1); 	CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in1); 						CHKERRQ(ierr);


		/* Update some info */
		VecGetArray(information,&information_vec);

		user->x_left 	=  information_vec[15];
		user->y_front 	=  information_vec[16];
		user->z_bot 	=  information_vec[17];
		user->W 		=  information_vec[18];
		user->L 		=  information_vec[19];
		user->H 		=  information_vec[20];

		/* Upload pushing block coordinates if required */
		if (user->Pushing.reset_pushing_coord == 0)
		{
			user->Pushing.x_center_block =  information_vec[21];
			user->Pushing.y_center_block =  information_vec[22];
			user->Pushing.z_center_block =  information_vec[23];
		}

		VecRestoreArray(information,&information_vec);

		/*Load erosion surface if FE Erosion Model was used*/
		if(rank==0){
			if (user->ErosionParameters.ErosionModel == 2){
				/* Load vector that contains erosion surface */
				ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileName_ErosionSurface,FILE_MODE_READ, &view_in1); CHKERRQ(ierr);
				ierr = VecLoad(user->ErosionParameters.FE_ErosionCode.ErosionSurface,view_in1); 		CHKERRQ(ierr);
				ierr = PetscViewerDestroy(&view_in1); 						CHKERRQ(ierr);
				/*Set coordinates of the Erosion_Code DM*/
				ierr = DMDASetUniformCoordinates(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,0,0); CHKERRQ(ierr);
			}
		}

		/* Load global data */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,LoadFileName_VCoord,FILE_MODE_READ, &view_in); CHKERRQ(ierr);

		ierr = VecLoad(Velocity, view_in); 			CHKERRQ(ierr);

		/* Read coordinates and add them in a consistent manner */
		ierr = DMGetCoordinates(da_nodes,			&coord); 	CHKERRQ(ierr);		/* Get global coordinate vector */
		ierr = VecLoad(coord, view_in); 				CHKERRQ(ierr);
		ierr = DAUpdatedGhostedCoordinates( da_nodes ); 		CHKERRQ(ierr);

		ierr = PetscViewerDestroy(&view_in); 					CHKERRQ(ierr);




		// Load Pressure	
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,LoadFileName_P,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
		ierr = VecLoad(Pressure, view_in); 			CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in); 					CHKERRQ(ierr);


		// Load temperature & coordinates [blocksize 1]
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,LoadFileName_T,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
#ifdef TEMPERATURE
		ierr = VecLoad(Temp, view_in); 				CHKERRQ(ierr);
#endif
		ierr = PetscViewerDestroy(&view_in); 					CHKERRQ(ierr);


		// Load internal free surface
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,LoadFileName_Surface,FILE_MODE_READ, &view_in); CHKERRQ(ierr);

		ierr = VecLoad(user->BottomTopography, view_in); 	CHKERRQ(ierr);		// surface & bottom topography
		ierr = VecLoad(user->SurfaceTopography, view_in); 	CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in); 					CHKERRQ(ierr);








		/* Check that input grid size is conform with current grid size */
		ind[0] = 0; ierr = VecGetValues(information,1, ind, n_read);CHKERRQ(ierr);
		if (  ((PetscInt) n_read[0]) != nx  ){
			PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in x-direction !! %lld %lld \n",(LLD)n_read[0], (LLD)nx);
			MPI_Abort(PETSC_COMM_WORLD,1);
		}
		ind[0] = 1; ierr = VecGetValues(information,1, ind, n_read);CHKERRQ(ierr);
		if (  ((PetscInt) n_read[0]) != ny  ){
			PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in y-direction !! %lld %lld \n",(LLD)n_read[0], (LLD)ny);
			MPI_Abort(PETSC_COMM_WORLD,1);
		}
		ind[0] = 2; ierr = VecGetValues(information,1, ind, n_read);CHKERRQ(ierr);
		if (  ((PetscInt) n_read[0]) != nz  ){
			PetscPrintf(PETSC_COMM_WORLD,"The mesh you want to read from file has a different size in z-direction !! %lld %lld \n",(LLD)n_read[0], (LLD)nz);
			MPI_Abort(PETSC_COMM_WORLD,1);
		}
		ind[0] = 3; ierr = VecGetValues(information,1, ind, n_read);CHKERRQ(ierr);
		if (  ((PetscInt) n_read[0]) != size  ){
			PetscPrintf(PETSC_COMM_WORLD,"The breakpoint file was created on %lld processors and you are running with %lld processors !! \n",(LLD)n_read[0], (LLD)size);
			MPI_Abort(PETSC_COMM_WORLD,1);
		}

		//
		ind[0] = 10; ierr = VecGetValues(information,1, ind, n_read); CHKERRQ(ierr);	itime 						= ((PetscInt)  n_read[0]	);
		ind[0] = 11; ierr = VecGetValues(information,1, ind, n_read); CHKERRQ(ierr);	user->time 					= ( n_read[0] 		);
		ind[0] = 12; ierr = VecGetValues(information,1, ind, n_read); CHKERRQ(ierr);	user->dt   					= ( n_read[0]		);
		ind[0] = 13; ierr = VecGetValues(information,1, ind, n_read); CHKERRQ(ierr);	user->break_point_number   	= ( (PetscInt) n_read[0] );
		ind[0] = 14; ierr = VecGetValues(information,1, ind, n_read); CHKERRQ(ierr);	user->ErosionParameters.HorizontalFreeSurfaceHeight   	= ( n_read[0] );




		user->time_start 		 = itime+1;
		user->break_point_number = user->break_point_number+1;

		/* Add time-dependent data to correct place */
		ierr = VecGetArray(TimeDependent_Data, (PetscScalar**)&TimeDependent_array); CHKERRQ(ierr);
		for (i=0; i<itime; i++){
			user->TimeDependentData[i] = TimeDependent_array[i];
		}
		ierr = VecRestoreArray(TimeDependent_Data,(PetscScalar**)&TimeDependent_array); CHKERRQ(ierr);



		/* (2) Load Particle data - each PROC has its own file */

#ifdef PARTICLES
		/* load particles */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileNameParticles,FILE_MODE_READ, &view_in1); CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&Particle_Vec);CHKERRQ(ierr);
		ierr = VecLoad(Particle_Vec,view_in1); 		CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in1); CHKERRQ(ierr);

		/* The number of particles loaded on the current PROC is equal to the length of the vector */
		ierr = VecGetSize(Particle_Vec, &numParticles);				CHKERRQ(ierr);
		user->num_particle_local = numParticles/particle_props;


		// Put particles at the correct location
		ierr = VecGetArray(Particle_Vec,(PetscScalar**)&Particle_array); CHKERRQ(ierr);
		for (ipart=0; ipart<user->num_particle_local; ipart++){
			user->ParticlesLocal[ipart] = Particle_array[ipart];
		}
		ierr = VecRestoreArray(Particle_Vec,(PetscScalar**)&Particle_array); CHKERRQ(ierr);

		ierr = VecDestroy(&Particle_Vec); 		CHKERRQ(ierr);

#endif





		/* Update coordinates of DMDA's based on the potentially new size of the computational domain */
		ierr = DMDASetUniformCoordinates(user->DA_SurfaceTopography,   user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L,	user->z_bot, user->z_bot+ user->H);	CHKERRQ(ierr);
		ierr = DMDASetUniformCoordinates(user->DA_BottomTopography,  	user->x_left, user->x_left + user->W, user->y_front, user->y_front + user->L,	user->z_bot, user->z_bot+ user->H);	CHKERRQ(ierr);

		if (element_type == DAVP_FDSTAG){
			ierr = SetUniformCoordinates_FDSTAG( user ); CHKERRQ(ierr);
		}



		PetscPrintf(PETSC_COMM_SELF,"  Reading timestep %lld | time=%g | dt=%g  | # particles = %lld | rank[%lld] \n", (LLD)itime, user->time, user->dt, (LLD)(user->num_particle_local), (LLD)rank);


		/* Load timestep */
		ierr = PetscOptionsGetReal(PETSC_NULL ,"-dt_start",	&user->dt		, PETSC_NULL); CHKERRQ(ierr);	// use a different dt
		PetscPrintf(PETSC_COMM_WORLD,"  Timestep used =%g  \n", user->dt);

		PetscPrintf(PETSC_COMM_WORLD,"**********************  \n");
		PetscPrintf(PETSC_COMM_WORLD,"  \n");


		/* Tidy up */
		ierr = VecDestroy(&information); 		CHKERRQ(ierr);
		ierr = VecDestroy(&TimeDependent_Data); 	CHKERRQ(ierr);


	}
	PetscFunctionReturn(0);
}
/*==========================================================================================================*/






