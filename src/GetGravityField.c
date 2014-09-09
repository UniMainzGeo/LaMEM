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

GetGravityField.c

Created on: 06.02.2012
Author: Tobias Baumann

Major revisions July 2012

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "GetGravityField.h"
#include "Utils.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "GetSurfaceVelocity.h"




#undef __FUNCT__
#define __FUNCT__ "GetAiryIsostasy"
PetscErrorCode GetAiryIsostasy(UserContext *user, const PetscInt ngp_vel,grid3D_struct	grid3D, PetscInt itime)
{
	// This is an ADHOC implementation to compute the lithostatic pressure
	// it will soon be adapted to LaMEM canonical
	// Isostatic topography as a function of density of the asthenosphere and correction topography (T. Becker at al. 2014)

	PetscScalar              ***materials_array;
	MaterialsElementDynamic 	material_data;
	char 				       *FileName,*DirectoryName;
	PetscViewer					view_out;
	PetscMPIInt                 rank;
	DMDALocalInfo 		        linfo;
	PetscScalar                 Pref, Pref_send=0.0, Pcol;
	PetscInt                    ix,iy,iz,xs,ys,zs,xm,ym,zm;
	Vec                         gvec_Tiso, lvec_Tiso;
	PetscScalar                *array_Tiso;
	PetscInt                    cpu_z;
	PetscScalar                 global_SumAbsTiso,global_SumSqrsTiso;

	PetscErrorCode              ierr;
	PetscFunctionBegin;


	user->Optimisation.SumAbsVel  = 0.0;
	user->Optimisation.SumSqrsVel = 0.0;


	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

	// ADD_HOC assumption : 1 proc in z direction !!
	ierr = DMDAGetInfo(user->DA_Processors,0,0,0,0, 0,0,&cpu_z,0,0,0,0,0,0); CHKERRQ(ierr);
	if (cpu_z != 1)
	{
		SETERRQ2(PETSC_COMM_WORLD,1," Sorry, the adhoc implementation of isostasy doesn't allow more than 1 cpu in z-dir ",0,0);


	}
	// get dimensions of the DA
	ierr = DMDAGetLocalInfo(user->DA_Processors,&linfo); CHKERRQ(ierr);


	// get indices of density and coordinates
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	// create local vectors to store the ENTIRE isostatic topography
	ierr = VecCreateSeq(PETSC_COMM_SELF,linfo.mx*linfo.my,&lvec_Tiso); CHKERRQ(ierr);
	ierr = VecSet(lvec_Tiso,0.0); CHKERRQ(ierr);


	// access the pressure vector
	ierr = VecGetArray(lvec_Tiso,&array_Tiso); CHKERRQ(ierr);

	// access the materials
	ierr = DMDAVecGetArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);

	// compute the reference pressure
	if(user->Isostasy.ref_xi >= xs &&  user->Isostasy.ref_xi < xs+xm && user->Isostasy.ref_yi >= ys &&  user->Isostasy.ref_yi < ys+ym)
	{
		Pcol = 0.0;
		for (iz=zs; iz<zs+zm; iz++)
		{
			LaMEMSetMaterialDataMemoryFromArray( &material_data , user->Isostasy.ref_xi-xs, user->Isostasy.ref_yi-ys,iz-zs, ngp_vel, materials_array );
			Pcol  += material_data.Density[0]*user->Characteristic.Density *grid3D.dzh;
		}
		Pref_send = Pcol;
	}

	// send the reference value to all
	ierr = MPI_Allreduce( &Pref_send, &Pref, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     Reference pressure Pref = %g\n",Pref);	CHKERRQ(ierr);


	for (iy=ys; iy<ys+ym; iy++)
	{
		for (ix=xs; ix<xs+xm; ix++)
		{
			// sum over the entire column and compute the lithostatic pressure
			Pcol = 0.0;
			for (iz=zs; iz<zs+zm; iz++)
			{
				LaMEMSetMaterialDataMemoryFromArray( &material_data , ix-xs,iy-ys, iz-zs, ngp_vel, materials_array );
				Pcol  += material_data.Density[0] * user->Characteristic.Density * grid3D.dzh;
			}

			// Isostatic topography
			array_Tiso[iy*(linfo.mx-1)+ix+iy] = (Pref-Pcol) / user->Isostasy.ref_rho - user->Isostasy.corr_topo;

		}
	}

	// restore the materials
	ierr = DMDAVecRestoreArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);

	// restore the array
	ierr = VecRestoreArray(lvec_Tiso,&array_Tiso); CHKERRQ(ierr);


	// store all in all local vectors (we can safely sum because each cpu only contributed its local part)
	ierr = Sum_Vectors(PETSC_COMM_WORLD ,&lvec_Tiso,linfo.mx*linfo.my);					CHKERRQ(ierr);

	// scatter sequential vector lvec_x -> distributed vector gvec_x
	ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,linfo.mx*linfo.my,&gvec_Tiso);		CHKERRQ(ierr);
	ierr = VecSet(gvec_Tiso,0.0); CHKERRQ(ierr);
	ierr = VecSeq2VecMpi(rank,lvec_Tiso,&gvec_Tiso);									CHKERRQ(ierr);

	// destroy sequential vectors
	ierr = VecDestroy(&lvec_Tiso); 														CHKERRQ(ierr);


	// Save binaries
	if(user->Isostasy.SaveRef==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save reference data -- \n");					CHKERRQ(ierr);
		asprintf(&DirectoryName, "ReferenceData_%1.6lld",(LLD)itime);
		ierr = LaMEM_CreateOutputDirectory(DirectoryName); 										CHKERRQ(ierr);
		asprintf( &FileName,"%s/REF_IsostaticTopography.bin",DirectoryName);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     save file: %s \n",FileName);					CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileName, FILE_MODE_WRITE, &view_out);			CHKERRQ(ierr);
		ierr = VecView(gvec_Tiso,view_out);														CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out);													CHKERRQ(ierr);

		free(FileName);
		free(DirectoryName);
	}

	// Get misfit
	if(user->Optimisation.GetIt==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");	CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#-- get isostatic topography misfit -------\n");	CHKERRQ(ierr);

		ierr = GetMisfit_IsostaticTopography(user,gvec_Tiso,PETSC_COMM_WORLD); CHKERRQ(ierr);


		if (rank==0){
			global_SumAbsTiso = user->Optimisation.SumAbsTiso;
			global_SumSqrsTiso = user->Optimisation.SumSqrsTiso;
		}
		else
		{
			global_SumAbsTiso = 0.0;
			global_SumSqrsTiso = 0.0;
		}


		ierr = MPI_Allreduce(&global_SumAbsTiso,&user->Optimisation.SumAbsTiso,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
		ierr = MPI_Allreduce(&global_SumSqrsTiso,&user->Optimisation.SumSqrsTiso,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
	}

	// Destroy vectors
	ierr = VecDestroy(&gvec_Tiso); 														CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetMisfit_IsostaticTopography"
PetscErrorCode GetMisfit_IsostaticTopography(UserContext *user,Vec gvec_Tiso, MPI_Comm PLANE_COMM)
{
//	--- Declarations --------------------------------------------------------------------------------------------------
	// --- general ---
	PetscErrorCode		ierr;

	// --- file i/o ---
	PetscViewer			view_in;
	char				filename_in[PETSC_MAX_PATH_LEN];

	//  --- misfit calculation ---
	Vec 				gvec_TisoRef;
	PetscScalar 		err_Tiso;
	PetscScalar         l1_Tiso;

	PetscInt			size_TisoRef,i,num,global_num, NumOfPoints;
	PetscScalar			*garray_TisoRef,*garray_Tiso;



//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

//	0.) --- initialize variables ---
	err_Tiso 				= 0.0;
	l1_Tiso                 = 0.0;

// 1.) --- load reference vectors ---
	// file name for reference data
	sprintf(filename_in,"%s",user->Isostasy.RefDatFile2load);

	// initialize & zero vectors to store data
	ierr = VecDuplicate(gvec_Tiso,&gvec_TisoRef);										CHKERRQ(ierr);
	ierr = VecZeroEntries(gvec_TisoRef);												CHKERRQ(ierr);

	// load data
	ierr = PetscViewerBinaryOpen(PLANE_COMM,filename_in,FILE_MODE_READ,&view_in);		CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#  -- load reference data --\n");					CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     load file: %s\n",filename_in);					CHKERRQ(ierr);
	ierr = VecLoad(gvec_TisoRef,view_in);												CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in);												CHKERRQ(ierr);

// 2.) --- check whether reference velocities contain NANs ---
	/* There are a couple of possibilities why a vel-node contain NAN:
	 *  - it's a boundary node which was set to NAN to prevent the total error being influenced by the boundary conditions
	 *  - the vel-node was set to NAN on purpose since there is no vel data available
	 *  - it's NAN because interpolation failed because it lacked for data while creating a the ref velocities
	 */

	ierr = VecGetLocalSize(gvec_TisoRef,&size_TisoRef);									CHKERRQ(ierr);
	ierr = VecGetSize(gvec_TisoRef,&NumOfPoints);										CHKERRQ(ierr);

	num  = 0;
	ierr = VecGetArray(gvec_TisoRef	,&garray_TisoRef);									CHKERRQ(ierr);
	ierr = VecGetArray(gvec_Tiso		,&garray_Tiso);									CHKERRQ(ierr);
	for (i=0; i<size_TisoRef; i++) {
		if(isnan(garray_TisoRef[i])){

			// set NAN-components to zero so that Petsc vector routines work fine
			// also set the corresponding grid velocity to zero so that total error isn't effected
			garray_TisoRef[i] 	= 0.0;
			garray_Tiso[i] 		= 0.0;


			// if velocity at this node is NAN, also the total number of points needs to be corrected.
			num = num + 1;
		}
	}
	ierr = VecRestoreArray(gvec_Tiso		,&garray_Tiso);								CHKERRQ(ierr);
	ierr = VecRestoreArray(gvec_TisoRef	,&garray_TisoRef);								CHKERRQ(ierr);

	ierr = MPI_Allreduce(&num,&global_num,1,MPIU_INT,MPI_SUM,PLANE_COMM);				CHKERRQ(ierr);
	NumOfPoints = NumOfPoints - global_num;
	user->Optimisation.NSurfVel = (PetscScalar) NumOfPoints;

// 3.) --- Normalization --
	// Normalize with the maximum of the reference "true" signal
	ierr = VecScale(gvec_TisoRef,1/(user->Isostasy.TisoStdDev));						CHKERRQ(ierr);
	ierr = VecScale(gvec_Tiso   ,1/(user->Isostasy.TisoStdDev));						CHKERRQ(ierr);

// 4.) --- calculate error ---
	//  lvec_TisoRef = (TisoRef-Tiso),
	ierr = VecAXPBY(gvec_TisoRef,1.0,-1.0	,gvec_Tiso);								CHKERRQ(ierr);

	// l1 norm of lvec_VRef
	ierr = VecNorm(gvec_TisoRef,NORM_1,&l1_Tiso);									    CHKERRQ(ierr);

	user->Optimisation.SumAbsTiso = l1_Tiso;

	//lvec_TisoRef = Tiso.^2,
	ierr = VecPointwiseMult(gvec_TisoRef,gvec_TisoRef,gvec_TisoRef);					CHKERRQ(ierr);

	// lvec_TisoRef = Tiso.^2 * quality,
	//ierr = VecPointwiseMult(gvec_TisoRef,gvec_TisoRef,gvec_Vqual);						CHKERRQ(ierr);

	// lvec_TisoRef = sum(Tiso.^2 *quality),
	ierr = VecSum(gvec_TisoRef,&err_Tiso);												CHKERRQ(ierr);

	user->Optimisation.SumSqrsTiso = err_Tiso;

	// RMS = sqrt(1/N * sum(x_i.^2)  )
	//*global_VelMisfit = sqrt( (err_Tiso + err_Vy + err_Vz)/(3*NumOfPoints) );

	ierr = PetscPrintf(PLANE_COMM,"#     Global err_Tiso: %g \n",err_Tiso);				CHKERRQ(ierr);
	ierr = PetscPrintf(PLANE_COMM,"#     NumberofPoints: %lld \n",(LLD)NumOfPoints);	CHKERRQ(ierr);

// 5.) --- clean reference vectors ---
	ierr = VecDestroy(&gvec_TisoRef); 													CHKERRQ(ierr);


	PetscFunctionReturn(0);
}




















#undef __FUNCT__
#define __FUNCT__ "GetGravityField"
PetscErrorCode GetGravityField(UserContext *user, const PetscInt ngp_vel,PetscInt itime)
{
// --- Declarations ---------------------------------------------------------------------------------------------------

	// --- common ------------
	PetscErrorCode				ierr;
	PetscInt 					i,j,n,num,nondim;
	PetscViewer					view_out;
	PetscMPIInt					rank;
	char						filename_out[PETSC_MAX_PATH_LEN],*FileName,*DirectoryName;
	PetscScalar					result,rho_max;
	const PetscScalar 			G = 6.67384 * pow(10,-11);


	// --- 3d density grid ---
	PetscInt					xs,ys,zs,xm,ym,zm,ix,iy,iz,intp;
	grid3D_struct				grid3D;
	PetscScalar					**cornervec;
	PetscScalar 				***materials_array;
	MaterialsElementDynamic 	material_data;
	Vec 						lvec_rho, lvec_intpx, lvec_intpy, lvec_intpz;
	PetscScalar					*larray_rho,*larray_intpx,*larray_intpy,*larray_intpz, ifsheight;
	PetscInt					*larray_phase;

	// --- 2d survey ---------
	survey_struct				survey;
	Vec							gvec_survey_dg,lvec_survey_dg,lvec_survey_dg2save,lvec_survey_coords;
	PetscScalar					*larray_survey_coords,*larray_survey_dg;


//	Vec                         lvec_kernel;
//	PetscScalar                 *larray_kernel;
	PetscInt		            zeros;
	PetscLogDouble	            cputime_start, cputime_end_single,cputime_end;
	PetscScalar                 xq[2],yq[2],zq[2];

// --- Code -----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	if(ngp_vel !=1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"--- ERROR ngp_vel != 1 \n");								CHKERRQ(ierr);
	}

	// --- initialize variables ---
	if(  ((user->GravityField.UseNumerics == PETSC_FALSE) && (user->GravityField.UseAnalytics == PETSC_FALSE)) || ((user->GravityField.UseNumerics == PETSC_TRUE) && (user->GravityField.UseAnalytics == PETSC_TRUE))){
		user->GravityField.UseNumerics  = PETSC_TRUE;
		user->GravityField.UseAnalytics = PETSC_FALSE;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     No integration option available \n");				CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     Integration performed numerically (default) \n");	CHKERRQ(ierr);
	}
	else if((user->GravityField.UseNumerics == PETSC_TRUE) && (user->GravityField.UseAnalytics == PETSC_FALSE)){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     Integration performed numerically \n");				CHKERRQ(ierr);
	}
	else if((user->GravityField.UseNumerics == PETSC_FALSE) && (user->GravityField.UseAnalytics == PETSC_TRUE)){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     Integration performed analytically \n");				CHKERRQ(ierr);
	}

	// --- non-dimensionalize ? ---
	nondim = 0;
	result		= 0.0;
	//rho_max 	= 0.0;
	rho_max     = user->GravityField.ReferenceDensity;

	grid3D.xs 	= 0.0; 	grid3D.xm = 0.0;
	grid3D.ys 	= 0.0; 	grid3D.ym = 0.0;
	grid3D.zs 	= 0.0; 	grid3D.zm = 0.0;


	// --- define 3D grid spacing + grid properties (dimensional units) ---
	if (nondim == 1){
		grid3D.xs = (user->x_left);		grid3D.xm = ((user->x_left) +(user->W));
		grid3D.ys = (user->y_front);	grid3D.ym = ((user->y_front)+(user->L));
		grid3D.zs = (user->z_bot);		grid3D.zm = ((user->z_bot)	+(user->H));
	}
	else{
		grid3D.xs = (user->x_left) *(user->Characteristic.Length);	grid3D.xm = ((user->x_left) +(user->W))*(user->Characteristic.Length);
		grid3D.ys = (user->y_front)*(user->Characteristic.Length); 	grid3D.ym = ((user->y_front)+(user->L))*(user->Characteristic.Length);
		grid3D.zs = (user->z_bot)  *(user->Characteristic.Length); 	grid3D.zm = ((user->z_bot)	+(user->H))*(user->Characteristic.Length);

	}


	// --- define 3D grid spacing + grid properties ---
	grid3D.dx		= (grid3D.xm-grid3D.xs)/(PetscScalar)(user->nel_x);	grid3D.dxh = grid3D.dx/2.0;
	grid3D.dy		= (grid3D.ym-grid3D.ys)/(PetscScalar)(user->nel_y);	grid3D.dyh = grid3D.dy/2.0;
	grid3D.dz		= (grid3D.zm-grid3D.zs)/(PetscScalar)(user->nel_z);	grid3D.dzh = grid3D.dz/2.0;
	grid3D.cellvol	= grid3D.dx*grid3D.dy*grid3D.dz;




	// Compute isostasy
	if (user->Isostasy.GetIt)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,">>>>>>>>>>>>>> ISO 1\n");	CHKERRQ(ierr);
		ierr = GetAiryIsostasy(user,ngp_vel,grid3D,itime); CHKERRQ(ierr);
	}


// --- initialize survey ----------------------------------------------------------------------------------------------

	// --- create survey + its coordinates (dimensionalized) ---
	survey.nx = (user->GravityField.survey_nx);
	survey.ny = (user->GravityField.survey_ny);

	survey.xs = (user->GravityField.survey_xs);
	survey.xm = (user->GravityField.survey_xm);
	survey.ys = (user->GravityField.survey_ys);
	survey.ym = (user->GravityField.survey_ym);
	survey.z  = (user->GravityField.survey_z) ;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     gravity profile [m]: \n");								CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     xs: %g, xm: %g, nx: %lld\n",survey.xs, survey.xm,(LLD)survey.nx);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     ys: %g, ym: %g, ny: %lld\n",survey.ys, survey.ym,(LLD)survey.ny);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     z : %g\n",survey.z);										CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     rhoref : %g\n",rho_max);										CHKERRQ(ierr);
	if (nondim == 1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     dimensionless calculation \n");						CHKERRQ(ierr);
	}
	else{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     calculation with SI-dimensional units \n");			CHKERRQ(ierr);
	}

	for(i=0;i<=user->GravityField.LithColNum;i++)
	{
		if (i<user->GravityField.LithColNum)
		{
			ierr = PetscPrintf(PETSC_COMM_WORLD,"#     LithColDepth[%lld]: %g \n",(LLD)i,user->GravityField.LithColDepth[i]);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     LithColDens[%lld] : %g \n",(LLD)i,user->GravityField.LithColDens[i]);CHKERRQ(ierr);
	}

	// --- non dimensionalize ? ---
	if (nondim == 1){
		survey.xs = survey.xs/(user->Characteristic.Length);
		survey.xm = survey.xm/(user->Characteristic.Length);
		survey.ys = survey.ys/(user->Characteristic.Length);
		survey.ym = survey.ym/(user->Characteristic.Length);
		survey.z  = survey.z /(user->Characteristic.Length);

	}



	// create survey + its coordinates
	ierr = VecCreateSeq(PETSC_COMM_SELF,2*(survey.nx*survey.ny),&lvec_survey_coords);					CHKERRQ(ierr);
	ierr = VecSet(lvec_survey_coords,0.0);																CHKERRQ(ierr);

	ierr = VecCreateSeq(PETSC_COMM_SELF,(survey.nx*survey.ny),&lvec_survey_dg);							CHKERRQ(ierr);
	ierr = VecSet(lvec_survey_dg,0.0);																	CHKERRQ(ierr);


	// define 2D survey spacing
	survey.dx = (survey.xm-survey.xs)/(survey.nx-1);
	survey.dy = (survey.ym-survey.ys)/(survey.ny-1);

	// create a coordinate array and store it to the vector
	ierr = VecGetArray(lvec_survey_coords,&larray_survey_coords);										CHKERRQ(ierr);
	n	 = 0;
	for (i=0; i<survey.nx; i++) {
		for (j=0; j<((survey.ny*2)); j+=2) {
			n 							=	i*survey.ny*2+j;
			larray_survey_coords[n]		=	survey.xs+i*survey.dx;
			larray_survey_coords[n+1]	=	survey.ys+j/2*survey.dy;
		}
	}
	ierr = VecRestoreArray(lvec_survey_coords,&larray_survey_coords);									CHKERRQ(ierr);





// --- initialize 3d grid ---------------------------------------------------------------------------------------------

	// --- get indices of density and coordinates ---
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); 								CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscInt),&larray_phase);								CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&larray_rho); 								CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&larray_intpx); 							CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&larray_intpy); 							CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(xm*ym*zm*ngp_vel)*sizeof(PetscScalar),&larray_intpz); 							CHKERRQ(ierr);


	// --- get maximum density - locally ---
	ierr = DMDAVecGetArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);
	num  = 0;
	ifsheight = user->ErosionParameters.InitialFreeSurfaceHeight*user->Characteristic.Length;


	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for (ix=xs; ix<xs+xm; ix++){

				// does this needs to be destroyed again??
				LaMEMSetMaterialDataMemoryFromArray( &material_data, ix-xs,iy-ys,iz-zs, ngp_vel, materials_array );

				intp = 0;

				if (nondim==1){
					larray_intpx[num] = material_data.Coord[0][intp];
					larray_intpy[num] = material_data.Coord[1][intp];
					larray_intpz[num] = material_data.Coord[2][intp];
				}
				else
				{
					larray_intpx[num] = material_data.Coord[0][intp]*user->Characteristic.Length;
					larray_intpy[num] = material_data.Coord[1][intp]*user->Characteristic.Length;
					larray_intpz[num] = material_data.Coord[2][intp]*user->Characteristic.Length;
				}

				larray_rho[num]	= material_data.Density[intp]*user->Characteristic.Density;

				// What is the difference between material_data.Density[intp] and user.FDSTAG.Center_Density

				//PetscPrintf(PETSC_COMM_,SELF"> Rho: %g \n",larray_rho[num]);




				if (larray_intpz[num] >= (ifsheight-user->GravityField.LithColDepth[0]))
				{
					larray_rho[num]	= larray_rho[num] - user->GravityField.LithColDens[0];
				}
				else
				{

					for(i=0;i<(user->GravityField.LithColNum -1);i++)
					{

						if((larray_intpz[num] <  (ifsheight-user->GravityField.LithColDepth[i] )) && (larray_intpz[num] >=  (ifsheight-user->GravityField.LithColDepth[i+1])))
						{
							larray_rho[num]	= larray_rho[num] - user->GravityField.LithColDens[i+1];
						}

					}
					if(larray_intpz[num] <  (ifsheight-user->GravityField.LithColDepth[user->GravityField.LithColNum-1]))
					{
						larray_rho[num]	= larray_rho[num] - user->GravityField.LithColDens[user->GravityField.LithColNum];
					}

				}


/*
				if (larray_intpz[num] >= (ifsheight-user->GravityField.LithColDens[0]))
				{
					larray_rho[num]	= larray_rho[num] - user->GravityField.LithColDens[0];
				}
				else if ((larray_intpz[num] <  (ifsheight-20.0e3)) && (larray_intpz[num] >=  (ifsheight-30.0e3)))
				{
					larray_rho[num]	= larray_rho[num] - 2900;
				}
				else if (larray_intpz[num] <  (ifsheight-30.0e3) && (larray_intpz[num] >=  (ifsheight-90.0e3)))
				{
					larray_rho[num]	= larray_rho[num] - 3350;
				}
				else if (larray_intpz[num] <  (ifsheight-user->GravityField.LithColDepth[user->GravityField.LithColNum-1]) )
				{
					larray_rho[num]	= larray_rho[num] - user->GravityField.LithColDens[user->GravityField.LithColNum];
				}
*/

				// set sticky air density to zero, we don't want to have the density normalized to sticky air density
				larray_phase[num] =	(PetscInt) material_data.Phases[intp];
				if((user->ErosionParameters.UseInternalFreeSurface == 1) && (larray_phase[num] == user->ErosionParameters.StickyAirPhase))
				{
					larray_rho[num]	= 0.0;
				}

//				PetscPrintf(PETSC_COMM_WORLD,"> Rho_corrected: %g \n",larray_rho[num]);

				num = num+1;
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials,user->Materials, &materials_array); CHKERRQ(ierr);

	// --- get global maximum density ---
//	ierr = MPI_Allreduce(&rho_max_local,&rho_max,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);				CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     rho_max: %g \n",rho_max);								CHKERRQ(ierr);

// --- basic calculations ---------------------------------------------------------------------------------------------

/*
	// --- kernel vector -------------------------------------------------------
	ierr = VecCreateSeq(PETSC_COMM_SELF,(survey.nx*survey.ny*xm*ym*zm),&lvec_kernel);                   CHKERRQ(ierr);
	ierr = VecGetArray(lvec_kernel,&larray_kernel);			                                            CHKERRQ(ierr);
	//nkernel = 0;
*/

	// get arrays
	ierr = VecGetArray(lvec_survey_coords,&larray_survey_coords);										CHKERRQ(ierr);
	ierr = VecGetArray(lvec_survey_dg,&larray_survey_dg);												CHKERRQ(ierr);
	ierr = LaMEMCreate2dArray(3,8,&cornervec,PETSC_NULL);												CHKERRQ(ierr);

	zeros = 0;
	num  = 0;
	n	 = 0;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);														CHKERRQ(ierr);
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscTime(&cputime_start);

	for (i=0; i<survey.nx; i++)                 // --- SURVEY LOOP ---
	{
		for (j=0; j<survey.ny*2; j+=2)          // --- SURVEY LOOP ---
		{
			for (iz=zs; iz<zs+zm; iz++)         // --- GRID LOOP ---
			{
				for (iy=ys; iy<ys+ym; iy++)     // --- GRID LOOP ---
				{
					for (ix=xs; ix<xs+xm; ix++) // --- GRID LOOP ---
					{


						survey.x 		= 	larray_survey_coords[i*survey.ny*2+j];
						survey.y 		= 	larray_survey_coords[i*survey.ny*2+j+1];

						if ((user->ErosionParameters.UseInternalFreeSurface == 1) && (larray_phase[num] == user->ErosionParameters.StickyAirPhase)){
							//ierr = PetscPrintf(PETSC_COMM_WORLD,"--- DEBUG sticky air phase \n");CHKERRQ(ierr);
						}
						else if (larray_rho[num] == 0.0){
							//ierr = PetscPrintf(PETSC_COMM_WORLD,"--- DEBUG not sticky air but still zero \n");CHKERRQ(ierr);
							zeros++;
						}
						else
						{
							if (user->GravityField.UseNumerics  == PETSC_TRUE )
							{
								cornervec[0][0] =larray_intpx[num]-grid3D.dxh; cornervec[1][0] =larray_intpy[num]-grid3D.dyh; cornervec[2][0] =larray_intpz[num]-grid3D.dzh; //1
								cornervec[0][1] =larray_intpx[num]+grid3D.dxh; cornervec[1][1] =larray_intpy[num]-grid3D.dyh; cornervec[2][1] =larray_intpz[num]-grid3D.dzh; //4
								cornervec[0][2] =larray_intpx[num]+grid3D.dxh; cornervec[1][2] =larray_intpy[num]+grid3D.dyh; cornervec[2][2] =larray_intpz[num]-grid3D.dzh; //6
								cornervec[0][3] =larray_intpx[num]-grid3D.dxh; cornervec[1][3] =larray_intpy[num]+grid3D.dyh; cornervec[2][3] =larray_intpz[num]-grid3D.dzh; //7
								cornervec[0][4] =larray_intpx[num]-grid3D.dxh; cornervec[1][4] =larray_intpy[num]-grid3D.dyh; cornervec[2][4] =larray_intpz[num]+grid3D.dzh; //2
								cornervec[0][5] =larray_intpx[num]+grid3D.dxh; cornervec[1][5] =larray_intpy[num]-grid3D.dyh; cornervec[2][5] =larray_intpz[num]+grid3D.dzh; //3
								cornervec[0][6] =larray_intpx[num]+grid3D.dxh; cornervec[1][6] =larray_intpy[num]+grid3D.dyh; cornervec[2][6] =larray_intpz[num]+grid3D.dzh; //5
								cornervec[0][7] =larray_intpx[num]-grid3D.dxh; cornervec[1][7] =larray_intpy[num]+grid3D.dyh; cornervec[2][7] =larray_intpz[num]+grid3D.dzh; //8

								ierr = GetGravityEffectNumerical(grid3D.cellvol,(user->GravityField.num_intp),survey.x,survey.y,survey.z,cornervec,&result);CHKERRQ(ierr);
//								larray_kernel[nkernel] = result;
								larray_survey_dg[n] = larray_survey_dg[n] + larray_rho[num]*result;
								result              = 0.0;

							}
							else if (user->GravityField.UseAnalytics  == PETSC_TRUE )
							{
								xq[0] = survey.x - (larray_intpx[num]-grid3D.dxh);
								xq[1] = survey.x - (larray_intpx[num]+grid3D.dxh);
								yq[0] = survey.y - (larray_intpy[num]-grid3D.dyh);
								yq[1] = survey.y - (larray_intpy[num]+grid3D.dyh);
								zq[0] = survey.z - (larray_intpz[num]-grid3D.dzh);
								zq[1] = survey.z - (larray_intpz[num]+grid3D.dzh);

								ierr = GetGravityEffectAnalytical2(xq,yq,zq,&result);CHKERRQ(ierr);
//								larray_kernel[nkernel] = result;
								larray_survey_dg[n] = larray_survey_dg[n] + larray_rho[num]*result;
								result = 0.0;

					}
							else
							{
								PetscPrintf(PETSC_COMM_WORLD,"#    WARNING: method not implemented \n");CHKERRQ(ierr);
							}
						}

//						nkernel          = nkernel+1;
						num			  	 = num+1;
					}								// --- GRID LOOP ---
				}									// --- GRID LOOP ---
			}										// --- GRID LOOP ---
			// set counter back to zero
			num = 0;
			//next survey point
			n=n+1;
		}											// --- SURVEY LOOP ---
	}												// --- SURVEY LOOP ---

	PetscTime(&cputime_end_single);
	PetscPrintf(PETSC_COMM_SELF,"# [%lld] Gravity loop took [%f s]for %lld cells and %lld zerocells\n",rank,cputime_end_single-cputime_start,(LLD) (xm*ym*zm),(LLD) zeros);
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"# Gravity loop took [%f s]\n",cputime_end-cputime_start);
//	PetscPrintf(PETSC_COMM_SELF,"rank %lld : nkernel  %lld \n",(LLD)rank,(LLD)nkernel);


	// --- restore arrays of local vectors ---
//	ierr = VecRestoreArray(lvec_kernel,&larray_kernel);     			                                CHKERRQ(ierr);
	ierr = VecRestoreArray(lvec_survey_dg,&larray_survey_dg);			                                CHKERRQ(ierr);
	ierr = VecRestoreArray(lvec_survey_coords,&larray_survey_coords);                                   CHKERRQ(ierr);
	ierr = LaMEMDestroy2dArray(&cornervec,PETSC_NULL);					                                CHKERRQ(ierr);



// --- gathering process ----------------------------------------------------------------------------------------------

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);														CHKERRQ(ierr);
	ierr = Sum_Vectors(PETSC_COMM_WORLD ,&lvec_survey_dg,(survey.nx*survey.ny));						CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,(survey.nx*survey.ny),&gvec_survey_dg);			CHKERRQ(ierr);
	ierr = VecSet(gvec_survey_dg,0.0);																	CHKERRQ(ierr);
	ierr = VecSeq2VecMpi(rank,lvec_survey_dg,&gvec_survey_dg);											CHKERRQ(ierr);


// --- Unnormalized data ----------------------------------------------------------------------------------------------

// --- save binaries --------------------------------------------------------------------------------------------------

/*
	if(1)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save kernel -- \n");								CHKERRQ(ierr);

		// --- create directory ---
		asprintf(&DirectoryName, "GravityKernels");
		ierr = LaMEM_CreateOutputDirectory(DirectoryName); 												CHKERRQ(ierr);

		// --- create filename ---
		sprintf(filename_out,"%s/gk_%lld.bin",DirectoryName,(LLD)rank);
		ierr = PetscPrintf(PETSC_COMM_SELF,"#     save file: %s (%g ) \n",filename_out,(survey.nx*survey.ny*xm*ym*zm));					CHKERRQ(ierr);


		// --- save binary file ---
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename_out,FILE_MODE_WRITE,&view_out); 			CHKERRQ(ierr);
		ierr = VecView(lvec_kernel, 				view_out); 												CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out);			 												CHKERRQ(ierr);

		// --- destroy output vectors ---
		ierr = VecDestroy(&lvec_kernel 	);																CHKERRQ(ierr);
		free(DirectoryName);

	}
*/




	if(user->GravityField.SaveDebug==1){

		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save debug data -- \n");								CHKERRQ(ierr);
		// --- create directory ---
		asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
		ierr = LaMEM_CreateOutputDirectory(DirectoryName); 												CHKERRQ(ierr);

		// --- create filename ---
		if(user->GravityField.UseAnalytics==PETSC_TRUE){
			//sprintf(filename_out,"GravityField_%s_Timestep_%1.6lld_analytical_%lld.dat",user->OutputFile,(LLD)itime,(LLD)rank);
			sprintf(filename_out,"%s/%s_GravityFieldAnalytical.%lld.%lld_ID%lld.out",DirectoryName, user->OutputFile,(LLD)rank,(LLD)itime+1000000LL,(LLD) user->Optimisation.mpi_group_id);
		}
		else{
			//sprintf(filename_out,"GravityField_%s_Timestep_%1.6lld_numerical_%lld.dat",user->OutputFile,(LLD)itime,(LLD)rank);
			sprintf(filename_out,"%s/%s_GravityFieldNumerical.%lld.%lld_ID%lld.out",DirectoryName, user->OutputFile,(LLD)rank,(LLD)itime+1000000LL,(LLD) user->Optimisation.mpi_group_id);
		}

		ierr = PetscPrintf(PETSC_COMM_SELF,"#     save file: %s \n",filename_out);						CHKERRQ(ierr);

		// --- create output vectors ---
		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,larray_rho  ,&lvec_rho);  		CHKERRQ(ierr);	ierr = VecAssemblyBegin(lvec_rho);    CHKERRQ(ierr);			ierr = VecAssemblyEnd(lvec_rho);  CHKERRQ(ierr);
		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,larray_intpx,&lvec_intpx);  	CHKERRQ(ierr);	ierr = VecAssemblyBegin(lvec_intpx);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(lvec_intpx);  CHKERRQ(ierr);
		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,larray_intpy,&lvec_intpy);  	CHKERRQ(ierr);	ierr = VecAssemblyBegin(lvec_intpy);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(lvec_intpy);  CHKERRQ(ierr);
		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,xm*ym*zm*ngp_vel,larray_intpz,&lvec_intpz); 		CHKERRQ(ierr);	ierr = VecAssemblyBegin(lvec_intpz);  CHKERRQ(ierr);			ierr = VecAssemblyEnd(lvec_intpz);  CHKERRQ(ierr);

		// --- multiplying with gravitational constant ---
		ierr = VecScale(lvec_survey_dg,G);																CHKERRQ(ierr);

		// --- save binary file ---
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename_out,FILE_MODE_WRITE,&view_out); 			CHKERRQ(ierr);
		ierr = VecView(lvec_survey_coords,		view_out); 												CHKERRQ(ierr);
		ierr = VecView(lvec_survey_dg,			view_out); 												CHKERRQ(ierr);
		ierr = VecView(lvec_rho, 				view_out); 												CHKERRQ(ierr);
		ierr = VecView(lvec_intpx, 				view_out); 												CHKERRQ(ierr);
		ierr = VecView(lvec_intpy, 				view_out); 												CHKERRQ(ierr);
		ierr = VecView(lvec_intpz, 				view_out); 												CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out);			 												CHKERRQ(ierr);

		// --- destroy output vectors ---
		ierr = VecDestroy(&lvec_rho 	);																CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_intpx 	); 																CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_intpy 	); 																CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_intpz 	); 																CHKERRQ(ierr);

		free(DirectoryName);
	}

	if(user->GravityField.SaveRef==1){


		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save reference data -- \n");							CHKERRQ(ierr);
		// --- create directory ---
		asprintf(&DirectoryName, "ReferenceData_%1.6lld",(LLD)itime);
		ierr = LaMEM_CreateOutputDirectory(DirectoryName); 												CHKERRQ(ierr);
		// --- create filename ---
		asprintf(&FileName,"%s/REF_Gravity.bin",DirectoryName);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     save file: %s \n",FileName);							CHKERRQ(ierr);
		// --- save binary file ---
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileName , FILE_MODE_WRITE, &view_out);			CHKERRQ(ierr);
		ierr = VecView(gvec_survey_dg,view_out); 														CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out); 															CHKERRQ(ierr);

		free(FileName);
		free(DirectoryName);

	}

// --- save VTKs ------------------------------------------------------------------------------------------------------
	if(user->GravityField.SaveVTK==1){

		// --- multiplying with gravitational constant ---
		ierr = VecDuplicate(lvec_survey_dg,&lvec_survey_dg2save);										CHKERRQ(ierr);
		ierr = VecCopy(lvec_survey_dg,lvec_survey_dg2save);												CHKERRQ(ierr);
		ierr = VecScale(lvec_survey_dg2save,G);															CHKERRQ(ierr);

		ierr = SaveGravityField2VTK(user,lvec_survey_dg2save, lvec_survey_coords,itime);				CHKERRQ(ierr);

		ierr = VecDestroy(&lvec_survey_dg2save);														CHKERRQ(ierr);
	}

// --- get misfit -----------------------------------------------------------------------------------------------------
	if(user->Optimisation.GetIt==1){
		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- get gravity field misfit --------------\n");

		ierr = GetMisfit_GravityField(user,gvec_survey_dg);												CHKERRQ(ierr);
	}

// --- clean up -------------------------------------------------------------------------------------------------------

	ierr = VecDestroy(&lvec_survey_coords);																CHKERRQ(ierr);
	ierr = VecDestroy(&lvec_survey_dg	);																CHKERRQ(ierr);
	ierr = VecDestroy(&gvec_survey_dg	);																CHKERRQ(ierr);

	ierr = PetscFree(larray_phase);																		CHKERRQ(ierr);
	ierr = PetscFree(larray_rho); 																		CHKERRQ(ierr);
	ierr = PetscFree(larray_intpx); 																	CHKERRQ(ierr);
	ierr = PetscFree(larray_intpy); 																	CHKERRQ(ierr);
	ierr = PetscFree(larray_intpz);																		CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetGravityEffectNumerical"
PetscErrorCode GetGravityEffectNumerical(PetscScalar dV,PetscInt num_int,PetscScalar survey_x,PetscScalar survey_y ,PetscScalar survey_z,PetscScalar **C ,PetscScalar *result)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	PetscInt		i,j;
	PetscScalar		dummy,x,y,z,r1p5;
	PetscScalar		detJacobian;
	PetscScalar		exp_pow;

//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;


	// --- Initialization ---
	detJacobian = 0.125*dV;
	i			= 0;
	j			= 0;
	dummy 		= 0.0;
	x 			= 0.0;
	y 			= 0.0;
	z 			= 0.0;
	r1p5 		= 0.0;
	exp_pow 	= 1.5;



	// 2 integration points per dimension
	if(num_int == 2)
	{

		PetscScalar		x_grid[8]={0.0},y_grid[8]={0.0},z_grid[8]={0.0};
		PetscScalar		A_int2[8][8];

		// hardcoded shapefunctions (trilinear,quadrilateral,"brickshaped"/cuboid):8 nodes (>) for 8 integrationpoints(v)
		// data was produced using following matlab routines: get_GLpointsweights_2intp.m / ShapeFunction_3D.m
		A_int2[0][0]=0.490562612162344; A_int2[0][1]=0.131445855765802; A_int2[0][2]=0.035220810900864; A_int2[0][3]=0.131445855765802;
		A_int2[1][0]=0.131445855765802; A_int2[1][1]=0.035220810900864; A_int2[1][2]=0.009437387837656; A_int2[1][3]=0.035220810900864;
		A_int2[2][0]=0.131445855765802; A_int2[2][1]=0.035220810900864; A_int2[2][2]=0.131445855765802; A_int2[2][3]=0.490562612162344;
		A_int2[3][0]=0.035220810900864; A_int2[3][1]=0.009437387837656; A_int2[3][2]=0.035220810900864; A_int2[3][3]=0.131445855765802;
		A_int2[4][0]=0.131445855765802; A_int2[4][1]=0.490562612162344; A_int2[4][2]=0.131445855765802; A_int2[4][3]=0.035220810900864;
		A_int2[5][0]=0.035220810900864; A_int2[5][1]=0.131445855765802; A_int2[5][2]=0.035220810900864; A_int2[5][3]=0.009437387837656;
		A_int2[6][0]=0.035220810900864; A_int2[6][1]=0.131445855765802; A_int2[6][2]=0.490562612162344; A_int2[6][3]=0.131445855765802;
		A_int2[7][0]=0.009437387837656; A_int2[7][1]=0.035220810900864; A_int2[7][2]=0.131445855765802; A_int2[7][3]=0.035220810900864;

		A_int2[0][4]=0.131445855765802; A_int2[0][5]=0.035220810900864; A_int2[0][6]=0.009437387837656; A_int2[0][7]=0.035220810900864;
		A_int2[1][4]=0.490562612162344; A_int2[1][5]=0.131445855765802; A_int2[1][6]=0.035220810900864; A_int2[1][7]=0.131445855765802;
		A_int2[2][4]=0.035220810900864; A_int2[2][5]=0.009437387837656; A_int2[2][6]=0.035220810900864; A_int2[2][7]=0.131445855765802;
		A_int2[3][4]=0.131445855765802; A_int2[3][5]=0.035220810900864; A_int2[3][6]=0.131445855765802; A_int2[3][7]=0.490562612162344;
		A_int2[4][4]=0.035220810900864; A_int2[4][5]=0.131445855765802; A_int2[4][6]=0.035220810900864; A_int2[4][7]=0.009437387837656;
		A_int2[5][4]=0.131445855765802; A_int2[5][5]=0.490562612162344; A_int2[5][6]=0.131445855765802; A_int2[5][7]=0.035220810900864;
		A_int2[6][4]=0.009437387837656; A_int2[6][5]=0.035220810900864; A_int2[6][6]=0.131445855765802; A_int2[6][7]=0.035220810900864;
		A_int2[7][4]=0.035220810900864; A_int2[7][5]=0.131445855765802; A_int2[7][6]=0.490562612162344; A_int2[7][7]=0.131445855765802;

		// perform A(8x8)*C(8x3) => results in x-,y-,z-coordinates for all 8 integration points
		for (i=0; i<8; i++) {
			for (j=0; j<8; j++) {
				x_grid[i] = A_int2[i][j]*C[0][j] + x_grid[i];
				y_grid[i] = A_int2[i][j]*C[1][j] + y_grid[i];
				z_grid[i] = A_int2[i][j]*C[2][j] + z_grid[i];
			}
		}

		// perform integral calculation for integration point coordinates; sum up and multiply with det(jacobian)
		// notice: this is only allowed because det(jacobian) is constant

		for (i=0; i<8; i++) {
			x		=	(double)survey_x	- (double)x_grid[i];
			y		=	(double)survey_y	- (double)y_grid[i];
			z		=	(double)survey_z	- (double)z_grid[i];

			r1p5	=   pow((double)(x*x+y*y+z*z),(double)exp_pow);
			dummy 	=  	z/r1p5 + dummy;
		}
		*result		= 	dummy * detJacobian;

	}



	// 1 integration point per dimension
	else if(num_int == 1) {
		PetscScalar A_int1[8] = {0.125000000000000,0.125000000000000,0.125000000000000,0.125000000000000,0.125000000000000,0.125000000000000,0.125000000000000,0.125000000000000};
		PetscScalar x_grid=0.0,y_grid=0.0,z_grid=0.0;

		// perform A(1x8)*C(8x3) => results in x-,y-,z-coordinates for a single integration point
		for (j=0; j<8; j++) {
			x_grid = A_int1[j]*C[0][j] + x_grid;
			y_grid = A_int1[j]*C[1][j] + y_grid;
			z_grid = A_int1[j]*C[2][j] + z_grid;
		}

		// weighting factor: 8 = 2*2*2
		x		=	survey_x	- x_grid;
		y		=	survey_y	- y_grid;
		z		=	survey_z	- z_grid;
		r1p5	=   pow((x*x+y*y+z*z),1.5);
		dummy 	= 	z/r1p5 *8* detJacobian;
		*result =	dummy;
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD,"WARNING ! (funct GetGravityEffectNumerical): your choice for num_int is not implemented yet !");
	}

	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetGravityEffectAnalytical"
PetscErrorCode GetGravityEffectAnalytical(PetscScalar survey_x,PetscScalar survey_y ,PetscScalar survey_z,PetscScalar **cornervec,PetscScalar *result)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	PetscScalar		x,y,z,r;
	PetscInt		ind;
	PetscScalar		dummy = 0.0;

//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;


	// get gravity effect of CURRENT cell
	for (ind=0; ind<4; ind++) {
		x	=	survey_x	- cornervec[0][ind];
		y	=	survey_y	- cornervec[1][ind];
		z	=	survey_z	- cornervec[2][ind];
		r	=   pow((x*x+y*y+z*z),0.5);
		dummy = dummy - (y*log((x+r)) + x*log((y+r)) - z* atan2((x*y),(z*r)));
	}

	for (ind=4; ind<8; ind++) {
		x	=	survey_x	- cornervec[0][ind];
		y	=	survey_y	- cornervec[1][ind];
		z	=	survey_z	- cornervec[2][ind];
		r	=   pow((x*x+y*y+z*z),0.5);

		dummy = dummy + (y*log((x+r)) + x*log((y+r)) - z* atan2((x*y),(z*r)));
	}



	*result = dummy;


	PetscFunctionReturn(0);
}



//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetGravityEffectAnalytical2"
PetscErrorCode GetGravityEffectAnalytical2(PetscScalar *x,PetscScalar *y, PetscScalar *z,PetscScalar *gsum)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	PetscInt		ind,i,j,k;
	PetscScalar		rijk=0.0,arg1,arg2,arg3;
	PetscScalar     ijk[] ={-1.0,1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0};
	PetscScalar     lsum=0.0;
//  --- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	ind  = 0;
	for (i=0;i<2;i++)
	{
	    for (j=0;j<2;j++)
	    {
	        for (k=0;k<2;k++)
	        {
	            //rijk = sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
	            rijk = pow(x[i]*x[i]+y[j]*y[j]+z[k]*z[k],0.5); // factor 10 faster on mogon
	            arg1 = atan2( x[i]*y[j], z[k]*rijk );
	            if(arg1 < 0 )
	            {
	            	arg1 =arg1 + 2*M_PI;
	            }
	            arg2 = log(rijk +y[j]);
	            arg3 = log(rijk +x[i]);
	            lsum = lsum + ijk[ind] *( z[k]*arg1  - x[i]*arg2 -y[j]*arg3 );
	            ind =ind+1;
	        }
	    }
	}

	*gsum = lsum;

	PetscFunctionReturn(0);
}

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "GetMisfit_GravityField"
PetscErrorCode GetMisfit_GravityField(UserContext *user, Vec lvec_dg)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	// --- general ---
	PetscErrorCode		ierr;

	// --- file i/o ---
	PetscViewer			view_in;
	char				filename_in[PETSC_MAX_PATH_LEN];

	// --- misfit calculation ---
	Vec 				lvec_dgRef;
	PetscScalar 		max_dgRefAbs,l1_norm,sum_squares,err_l1,err_l2;


//  --- Code ----------------------------------------------------------------------------------------------------------

	PetscFunctionBegin;

// 	0.) --- initialize variables ---
	max_dgRefAbs 			= 0.0;


// 	1.) --- load reference vector ---

	// filename
	//sprintf(filename_in,"GravityField_%s_REF.bin",user->OutputFile);
	sprintf(filename_in,"%s",user->GravityField.RefDatFile2load);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- load reference data -- \n");								CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     load file: %s \n",filename_in);							CHKERRQ(ierr);

	// initialize & zero vectors to store data
	ierr = VecDuplicate(lvec_dg,&lvec_dgRef);															CHKERRQ(ierr);
	ierr = VecZeroEntries(lvec_dgRef);																	CHKERRQ(ierr);

	// load data
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename_in,FILE_MODE_READ,&view_in);					CHKERRQ(ierr);
	ierr = VecLoad(lvec_dgRef, view_in);										        				CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in);																CHKERRQ(ierr);


// 	2.) --- Normalization + Error ---

	// global max(|dgRef|),

	// --- get error ---
	//max(|dgRef|)
	//ierr = VecAbsMax(lvec_dgRef,&max_dgRefAbs);															CHKERRQ(ierr);
	// normalize with the maximum of the reference gravity
	//ierr = VecScale(lvec_dgRef,1/max_dgRefAbs);															CHKERRQ(ierr);
	//ierr = VecScale(lvec_dg   ,1/max_dgRefAbs);															CHKERRQ(ierr);

	// normalize with std deviation
	ierr = VecScale(lvec_dgRef,1/(user->GravityField.StdDev));											CHKERRQ(ierr);
	ierr = VecScale(lvec_dg   ,1/(user->GravityField.StdDev));											CHKERRQ(ierr);



	// dgRef = dgRef - dg
	ierr = VecAXPBY(lvec_dgRef,1.0,-1.0,lvec_dg);														CHKERRQ(ierr);

	// sum of abs dgRef values
	ierr = VecNorm(lvec_dgRef,NORM_1,&l1_norm);															CHKERRQ(ierr);

	// dgRef = (dgRef).^2
	ierr = VecPointwiseMult(lvec_dgRef,lvec_dgRef,lvec_dgRef);							                CHKERRQ(ierr);
	// sum or squares
	ierr = VecSum(lvec_dgRef,&sum_squares);                                                             CHKERRQ(ierr);


	user->Optimisation.SumAbsGrav = l1_norm;
	user->Optimisation.SumSqrsGrav = sum_squares;
	user->Optimisation.NGrav = (PetscScalar) (user->GravityField.survey_nx*user->GravityField.survey_ny);


	// L2-norm, sqrt(sum(|dgRef -dg|^2))
	//ierr = VecNorm(lvec_dgRef,NORM_2,&err_dg);															CHKERRQ(ierr);
	// L1-norm, sum(|dgRef -dg|)
	//ierr = VecNorm(lvec_dgRef,NORM_1,&err_dg);															CHKERRQ(ierr);
	// rms = L2-norm/sqrt(length(dg))
	//user->Optimisation.MisfitGravity= err_dg / sqrt(user->GravityField.survey_nx * user->GravityField.survey_ny);
	// arithmetic mean = L1-norm/(length(dg))
	//user->Optimisation.MisfitGravity= err_dg / (user->GravityField.survey_nx * user->GravityField.survey_ny);

	// error based on l2 and l1 norm

	err_l2 = sqrt(user->Optimisation.SumSqrsGrav/(user->Optimisation.NGrav));
	err_l1 =      user->Optimisation.SumAbsGrav /(user->Optimisation.NGrav);

	PetscPrintf(PETSC_COMM_WORLD,"#     gravity max value:  %g\n",max_dgRefAbs);
	PetscPrintf(PETSC_COMM_WORLD,"#     gravity L2/sqrt(N): %g\n",err_l2);
	PetscPrintf(PETSC_COMM_WORLD,"#     gravity L1/(N):     %g\n",err_l1);



// 	3.) --- clean vectors ---
	ierr = VecDestroy(&lvec_dgRef); 																	CHKERRQ(ierr);



	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "SaveGravityField2VTK"
PetscErrorCode SaveGravityField2VTK(UserContext *user, Vec lvec_dg, Vec lvec_coords,PetscInt itime)
{
//	--- Declarations --------------------------------------------------------------------------------------------------

	// --- general ---
	PetscErrorCode		ierr;
	PetscInt 			i,j,n;
	PetscInt			survey_nx,survey_ny;
	PetscScalar			*larray_coords,*larray_dg;

	// --- file i/o ---
	char 				*vtk_filename,*DirectoryName;
	FILE 				*vtk_fp;


//  --- Code ----------------------------------------------------------------------------------------------------------


		PetscFunctionBegin;

		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save VTK file -- \n");								CHKERRQ(ierr);

		// --- create directory ---
		asprintf(&DirectoryName, "Timestep_%1.6lld",(LLD)itime);
		ierr = LaMEM_CreateOutputDirectory(DirectoryName); 												CHKERRQ(ierr);

		// --- create file name ---
		if(user->GravityField.UseAnalytics==PETSC_TRUE){
			asprintf( &vtk_filename,"%s/%s_GravityFieldAnalytical.0.%lld.vts",DirectoryName, user->OutputFile,(LLD)itime+1000000LL);
			//asprintf( &vtk_filename,"%s_GravityFieldAnalytical_%lld.vts", user->OutputFile,(LLD) user->mpi_group_id);
		}
		else{
			asprintf( &vtk_filename,"%s/%s_GravityFieldNumerical.0.%lld.vts",DirectoryName, user->OutputFile,(LLD)itime+1000000LL);
			//asprintf( &vtk_filename, "%s_GravityFieldNumerical_%lld.vts", user->OutputFile,(LLD) user->mpi_group_id);
		}


		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     save file: %s \n",vtk_filename);						CHKERRQ(ierr);


		// open file and write header
		vtk_fp = fopen( vtk_filename, "w" );
		if( vtk_fp == NULL ) {
			printf("ERROR(DAView_2DVTK_StructuredGrid): Cannot open file = %s \n", vtk_filename );
			exit(0);
		}

		// coords
		survey_nx = (user->GravityField.survey_nx);
		survey_ny = (user->GravityField.survey_ny);

		//	header
		fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
		fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
		fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n",(LLD)0,(LLD)(survey_ny-1),(LLD)0,(LLD)(survey_nx-1), 0LL,0LL);
		fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");
		fprintf( vtk_fp, "\t\t\t<Points>\n");

		//	coordinates
		ierr = VecGetArray(lvec_coords,&larray_coords);CHKERRQ(ierr);
		n=0;
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for( j=0; j<survey_nx; j++ ) {
			for( i=0; i<survey_ny; i++ ) {

				fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n",(double)larray_coords[n]*0.001,(double)larray_coords[n+1]*0.001,(double)(user->GravityField.survey_z)*0.001 );
				n=n+2;
			}
		}
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
		fprintf( vtk_fp, "\t\t\t</Points>\n");
		ierr = VecRestoreArray(lvec_coords,&larray_coords);CHKERRQ(ierr);


		// field names
		fprintf( vtk_fp, "\t\t\t<PointData Scalars=\" ");
		fprintf( vtk_fp, "%s ", "dg_z" );
		fprintf( vtk_fp, "\">\n");

		// field data
		ierr = VecGetArray(lvec_dg,&larray_dg);CHKERRQ(ierr);
		n=0;
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "Gravity in mGal");
		for( j=0; j<(survey_nx); j++ ) {
			for( i=0; i<(survey_ny); i++ ) {
				// convert m/s^2 in mGal* 10 ^5
				fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", (double)(larray_dg[n]*100000.0) );
				//fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", (double)(larray_dg[n]) );
				n=n+1;
			}
		}
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
		ierr = VecRestoreArray(lvec_dg,&larray_dg);CHKERRQ(ierr);

		fprintf( vtk_fp, "\t\t\t</PointData>\n");
		fprintf( vtk_fp, "\t</StructuredGrid>\n");
		fprintf( vtk_fp, "</VTKFile>\n");


		if( vtk_fp!= NULL ) {
			fclose( vtk_fp );
			vtk_fp = NULL;
		}


		free(vtk_filename);
		free(DirectoryName);

	PetscFunctionReturn(0);
}
