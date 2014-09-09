#include "LaMEM.h"
#include "Addition_FDSTAG.h"
#include "Assembly_FDSTAG.h"


/*==========================================================================================================*/
/* Modify and get information from PV_MAT - this needs to be done before the MatAssembly in ComputeStiffnessMatrix*/
#undef __FUNCT__
#define __FUNCT__ "AddPushingToModelPV_MAT"
PetscErrorCode AddPushingToModelPV_MAT(UserContext *user,Mat PV_MAT,Mat PV_MAT_push,Vec pv_rhs_push)

{
	PetscErrorCode  ierr;
	PetscInt	    *rowidx_array = PETSC_NULL;
	PetscInt	    *rowidy_array = PETSC_NULL;
	PetscInt        numRowsx,numRowsy;
	PetscBool		flg;
	PetscInt		i, ind_change;

	// Check if pushing option is activated
	if (user->AddPushing == 1)
	{
		//Add pushing BC only within the specified time interval
		if ((user->time>=user->Pushing.time[0]) && (user->time<=user->Pushing.time[user->Pushing.num_changes]))
		{
			// check which pushing stage (i.e. between 0-1 Myr: omega=1 deg/yr, 1-2 Myr: omega=0 deg/yr etc.)
			ind_change = 0;
			for (i=0;i<user->Pushing.num_changes;i++){
				if (user->time>=user->Pushing.time[i] && user->time<=user->Pushing.time[i+1]){
					ind_change = i;
				}}

			// Count number of indices rowid
			flg  = PETSC_FALSE;
			ierr = GetGlobalIndicesForPushing(user, PETSC_NULL, &numRowsx, PETSC_NULL, &numRowsy, flg); CHKERRQ(ierr);

			// Allocate memory to variable rowid_array (for Vx and Vy)
			if(numRowsx > 0)
			{	ierr = PetscMalloc((size_t)numRowsx*sizeof(PetscInt), &rowidx_array); CHKERRQ(ierr);
			}

			if(numRowsy > 0)
			{	ierr = PetscMalloc((size_t)numRowsy*sizeof(PetscInt), &rowidy_array); CHKERRQ(ierr);
			}

			// Get global indices for dofs in variable rowid_array
			flg  = PETSC_TRUE;
			ierr = GetGlobalIndicesForPushing(user, rowidx_array, &numRowsx, rowidy_array, &numRowsy, flg); CHKERRQ(ierr);
			ierr = PetscSortInt(numRowsx,rowidx_array); CHKERRQ(ierr);
			ierr = PetscSortInt(numRowsy,rowidy_array); CHKERRQ(ierr);

			// Modify PV_MAT and save the vector values to be added to the rhs
			ierr = ZeroColumnsPV_MAT(user,pv_rhs_push,PV_MAT,PV_MAT_push,rowidx_array,numRowsx,rowidy_array,numRowsy, ind_change);  CHKERRQ(ierr);

			// Free allocated memory
			ierr = PetscFree(rowidx_array); CHKERRQ(ierr);
			ierr = PetscFree(rowidy_array); CHKERRQ(ierr);
		}
	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Option for pushing */
#undef __FUNCT__
#define __FUNCT__ "AddPushingToModel"
PetscErrorCode AddPushingToModel(UserContext *user, Mat VV_MAT,Mat VP_MAT,Vec rhs,Vec rhs_p, Vec pv_rhs_push)
{ 
	PetscErrorCode  ierr;
	PetscInt	    *rowidx_array = PETSC_NULL;
	PetscInt	    *rowidy_array = PETSC_NULL;
	PetscInt        numRowsx,numRowsy;
	PetscLogDouble	cputime_start,cputime_end;
	PetscBool		flg;
	PetscInt		i, ind_change;
	PetscScalar		time0,time1,vpush,omega;
	PetscInt		c_adv,dir;

	// Check if pushing option is activated
	if (user->AddPushing == 1)
	{
		//Add pushing BC only within the specified time interval
		if ((user->time>=user->Pushing.time[0]) && (user->time<=user->Pushing.time[user->Pushing.num_changes]))
		{
			PetscTime(&cputime_start);

			// check which pushing stage (i.e. between 0-1 Myr: omega=1 deg/yr, 1-2 Myr: omega=0 deg/yr etc.)
			ind_change = 0;
			for (i=0;i<user->Pushing.num_changes;i++){
				if (user->time>=user->Pushing.time[i] && user->time<=user->Pushing.time[i+1]){
					ind_change = i;
				}
			}

			time0 = user->Pushing.time[ind_change]/1e6/user->Characteristic.SecYear*user->Characteristic.Time;
			time1 = user->Pushing.time[ind_change+1]/1e6/user->Characteristic.SecYear*user->Characteristic.Time;
			vpush = user->Pushing.V_push[ind_change]/1e-2*user->Characteristic.SecYear*user->Characteristic.Velocity;
			omega = user->Pushing.omega[ind_change]/(M_PI/180.0)*(user->Characteristic.SecYear/user->Characteristic.Time);
			c_adv = user->Pushing.coord_advect[ind_change];
			dir   = user->Pushing.dir[ind_change];

			PetscPrintf(PETSC_COMM_WORLD,"#  Pushing BC: Time interval = [%g - %g] Ma, V_push = %g cm/yr, Omega = %g deg/yr, mobile_block = %d, direction = %d\n",time0,time1,vpush,omega,c_adv,dir);
			PetscPrintf(PETSC_COMM_WORLD,"#  Pushing BC: Block center coordinates x = [%g], y = [%g], z = [%g]\n",user->Pushing.x_center_block*1000,user->Pushing.y_center_block*1000,user->Pushing.z_center_block*1000);

			// Count number of indices rowid
			flg  = PETSC_FALSE;
			ierr = GetGlobalIndicesForPushing(user, PETSC_NULL, &numRowsx, PETSC_NULL, &numRowsy, flg); CHKERRQ(ierr);

			//PetscPrintf(PETSC_COMM_SELF,"#  Adding pushing bc affects %lld Vx entries and  %lld Vy entries \n", (LLD)numRowsx,(LLD)numRowsy);

			// Allocate memory to variable rowid_array (for Vx, Vy)
			if(numRowsx > 0)
			{	ierr = PetscMalloc((size_t)numRowsx*sizeof(PetscInt), &rowidx_array); CHKERRQ(ierr);
			}

			if(numRowsy > 0)
			{	ierr = PetscMalloc((size_t)numRowsy*sizeof(PetscInt), &rowidy_array); CHKERRQ(ierr);
			}

			// Get global indices for dofs in variable rowid_array
			flg  = PETSC_TRUE;
			ierr = GetGlobalIndicesForPushing(user, rowidx_array, &numRowsx, rowidy_array, &numRowsy, flg); CHKERRQ(ierr);
			ierr = PetscSortInt(numRowsx,rowidx_array); CHKERRQ(ierr);
			ierr = PetscSortInt(numRowsy,rowidy_array); CHKERRQ(ierr);

			// Modify stiffness matrix in such way that modifying the rhs is not necessary
			ierr = ModifyStiffnessMatrixForPushing(user,VV_MAT,VP_MAT,rhs,rowidx_array,numRowsx,rowidy_array,numRowsy,ind_change); CHKERRQ(ierr);


			// Add the vector from PV_MAT to the rhs
			ierr = VecAXPY(rhs_p,1.0,pv_rhs_push); CHKERRQ(ierr);

			// Advect the Pushing block coordinates
			ierr = AdvectThePushingBlockCoordinates(user,ind_change); CHKERRQ(ierr);

			// Free allocated memory
			ierr = PetscFree(rowidx_array); CHKERRQ(ierr);
			ierr = PetscFree(rowidy_array); CHKERRQ(ierr);

			PetscTime(&cputime_end);
			PetscPrintf(PETSC_COMM_WORLD,"#  Adding pushing boundary conditions took %g s \n",cputime_end-cputime_start);
		}
	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
/* Get global indices for the velocity points in the pushing block */
#undef __FUNCT__
#define __FUNCT__ "GetGlobalIndicesForPushing"
PetscErrorCode GetGlobalIndicesForPushing(UserContext *user, PetscInt *rowidx_array, PetscInt *numRowsx, PetscInt *rowidy_array, PetscInt *numRowsy, PetscBool flg)
{
	//=======================================================
	// WARNING! CODE WILL ONLY WORK FOR CONSTANT MESH SPACING
	//=======================================================

	PetscErrorCode  ierr;
	DM              cda;
	Vec             gc;
	PetscInt        i,j,k,xm,ym,zm,xs,ys,zs;
	PetscInt        idx,idy,rowidx, rowidy;
	PetscScalar		xx, yx, zx, xy, yy, zy;
	PetscScalar		rxx, ryx, rzx, rxy, ryy, rzy;
	PetscScalar     xb_left, xb_right, yb_front, yb_back, zb_bottom, zb_top;
	DMDACoor3d      ***coord;
	PetscScalar     dx, dy, dz;
	PetscScalar		theta;

	theta = user->Pushing.theta; //total rotation angle

	ierr = DMGetCoordinates(user->DA_Vel, &gc);      CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(user->DA_Vel, &cda);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, gc, &coord); CHKERRQ(ierr);

	// loop strictly over the local part, not to list points multiple times
	ierr = DMDAGetCorners(user->DA_Vel, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

	// get half of mesh spacing
	dx = (user->W/(PetscScalar)user->nnode_x)*0.5;
	dy = (user->L/(PetscScalar)user->nnode_y)*0.5;
	dz = (user->H/(PetscScalar)user->nnode_z)*0.5;

	// Block corner coordinates - non-dimensional
	xb_left   = user->Pushing.x_center_block - user->Pushing.L_block*0.5;
	xb_right  = user->Pushing.x_center_block + user->Pushing.L_block*0.5;
	yb_front  = user->Pushing.y_center_block - user->Pushing.W_block*0.5;
	yb_back   = user->Pushing.y_center_block + user->Pushing.W_block*0.5;
	zb_bottom = user->Pushing.z_center_block - user->Pushing.H_block*0.5;
	zb_top    = user->Pushing.z_center_block + user->Pushing.H_block*0.5;

	idx = 0;
	idy = 0;

	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for (i=xs; i<xs+xm; i++){

				// prepare point coordinates *** implicit assumption constant mesh size!!! ***

				xx = coord[k][j][i].x;			//Vx
				yx = coord[k][j][i].y + dy;
				zx = coord[k][j][i].z + dz;

				xy = coord[k][j][i].x + dx;		//Vy
				yy = coord[k][j][i].y;
				zy = coord[k][j][i].z + dz;

				// Clockwise ROTATION since block rotation is taken positive when anti-clockwise
				// Rotation matrix R=[ cos() sin()
				//					  -sin() cos() ]
				rxx = cos(theta)*xx + sin(theta)*yx;
				ryx =-sin(theta)*xx + cos(theta)*yx;
				rzx = zx;

				rxy = cos(theta)*xy + sin(theta)*yy;
				ryy =-sin(theta)*xy + cos(theta)*yy;
				rzy = zy;

				// Check for Vx points inside the block
				if(rxx > xb_left   && rxx < xb_right	&& ryx > yb_front  && ryx < yb_back	&& rzx > zb_bottom && rzx < zb_top)
				{
					// PetscPrintf(PETSC_COMM_WORLD,"# C [%d %d %d] \n",i,j,k);
					//PetscPrintf(PETSC_COMM_WORLD,"# Row Index X = %d \n",idx);

					if (flg == PETSC_TRUE)
					{
						// get global index and save rowidx if inside the box
						ierr = DAGetGlobalIndex(user->DA_Vel, i, j, k, 0, &rowidx); CHKERRQ(ierr); // 0 - Vx; 1 - Vy; 2 - Vz
						rowidx_array[idx] = rowidx;
					}

					// update counter
					idx++;
				}
				// Check for Vy points inside the block
				if(rxy > xb_left   && rxy < xb_right	&& ryy > yb_front  && ryy < yb_back	&& rzy > zb_bottom && rzy < zb_top)
				{
					// PetscPrintf(PETSC_COMM_WORLD,"# C [%d %d %d] \n",i,j,k);
					//PetscPrintf(PETSC_COMM_WORLD,"# Row Index Y = %d \n",idy);

					if (flg == PETSC_TRUE)
					{
						// get global index and save rowidy if inside the box
						ierr = DAGetGlobalIndex(user->DA_Vel, i, j, k, 1, &rowidy); CHKERRQ(ierr); // 0 - Vx; 1 - Vy; 2 - Vz
						rowidy_array[idy] = rowidy;

					}

					// update counter
					idy++;
				}
			}}}

	*numRowsx = idx;
	*numRowsy = idy;

	ierr = DMDAVecRestoreArray(cda,gc,&coord); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/*Modify Stiffness matrix for pushing */
#undef __FUNCT__
#define __FUNCT__ "ModifyStiffnessMatrixForPushing"
PetscErrorCode ModifyStiffnessMatrixForPushing(UserContext *user,Mat VV_MAT,Mat VP_MAT,Vec rhs,PetscInt *rowidx_array, PetscInt numRowsx,PetscInt *rowidy_array, PetscInt numRowsy,PetscInt ind_change)
{
	PetscErrorCode  ierr;
	PetscScalar     v_vv, v_vp;
	PetscMPIInt     rank,size;
	IS				isx,isy;
	Vec				x_push;
//	PetscViewer		viewer;
	/*
	PetscInt 		x,y,i;
	FILE          *fp;
	char          *fname;

		//for debugging
	asprintf(&fname,"Pushing_dofs.dat");
	fp = fopen(fname,"w");

	fprintf(fp,"%lld\n", (LLD)numRowsx);
	for (i=0; i<numRowsx; i++){
		x = rowidx_array[i];
		fprintf(fp,"%lld\n", (LLD)x);
	}

	fprintf(fp,"%lld\n", (LLD)numRowsy);
		for (i=0; i<numRowsy; i++){
		y = rowidy_array[i];
		fprintf(fp,"%lld\n", (LLD)y);
	}

	fclose(fp);
	free(fname);
		//end output
	 */

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

	v_vv = 1.0;
	v_vp = 0.0;

	// Create vector to include the modified velocity values
	ierr = VecDuplicate(rhs,&x_push); CHKERRQ(ierr);
	ierr = VecZeroEntries(x_push); CHKERRQ(ierr);

	// Create a vector x_push that contains all the [V_push] elements that need to be replaced in the rhs
	ierr = ModifySolutionsVectorForPushing2(user,x_push,rowidx_array,numRowsx,rowidy_array,numRowsy,ind_change); CHKERRQ(ierr);

	// Create an IS required by MatZeroRows()
	ierr = ISCreateGeneral(PETSC_COMM_WORLD,numRowsx,rowidx_array,PETSC_COPY_VALUES,&isx);  CHKERRQ(ierr);
	ierr = ISCreateGeneral(PETSC_COMM_WORLD,numRowsy,rowidy_array,PETSC_COPY_VALUES,&isy);  CHKERRQ(ierr);

	// MatZeroRowsColumns() seems to produce parallel artefacts
//	ierr = MatZeroRowsColumnsIS(VV_MAT,isx,v_vv,x_push,rhs);  CHKERRQ(ierr);
//	ierr = MatZeroRowsColumnsIS(VV_MAT,isy,v_vv,x_push,rhs);  CHKERRQ(ierr);

	ierr = MatZeroRowsIS(VV_MAT,isx,v_vv,x_push,rhs);  CHKERRQ(ierr);
	ierr = MatZeroRowsIS(VV_MAT,isy,v_vv,x_push,rhs);  CHKERRQ(ierr);

	ierr = MatZeroRowsIS(VP_MAT,isx,v_vp,PETSC_NULL,PETSC_NULL);  CHKERRQ(ierr);
	ierr = MatZeroRowsIS(VP_MAT,isy,v_vp,PETSC_NULL,PETSC_NULL);  CHKERRQ(ierr);

	ierr = ISDestroy(&isx); CHKERRQ(ierr);
	ierr = ISDestroy(&isy); CHKERRQ(ierr);

	ierr = VecDestroy(&x_push); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


/*==========================================================================================================*/
/* Modify the vector of solutions for zeroed rows for pushing */
#undef __FUNCT__
#define __FUNCT__ "ModifySolutionsVectorForPushing2"
PetscErrorCode ModifySolutionsVectorForPushing2(UserContext *user,Vec x_push,PetscInt *rowidx_array,PetscInt numRowsx,PetscInt *rowidy_array,PetscInt numRowsy,PetscInt ind_change)
{
	PetscErrorCode  ierr;
	PetscScalar     Vx,Vy,theta;
	PetscInt 		ix, iy,i,rowidx,rowidy;
	PetscScalar 	*x_push_local;
	PetscInt        low, high;

	// access the x_push vector
	ierr = VecGetOwnershipRange(x_push, &low, &high); CHKERRQ(ierr);
	ierr = VecGetArray(x_push,&x_push_local); CHKERRQ(ierr);

	Vx = 0.0;
	Vy = 0.0;

	if (user->Pushing.dir[ind_change] == 0)
	{
		theta = user->Pushing.theta;
		Vx    = cos(theta)*user->Pushing.V_push[ind_change];
		Vy    = sin(theta)*user->Pushing.V_push[ind_change];
	}
	if (user->Pushing.dir[ind_change] == 1)
	{
		Vx    = user->Pushing.V_push[ind_change];
		Vy    = 0.0;
	}
	if (user->Pushing.dir[ind_change] == 2)
	{
		Vx    = 0.0;
		Vy    = user->Pushing.V_push[ind_change];
	}

	// Fill in Vx values
	ix=0;
	for (i=ix; i<numRowsx; i++)
	{
		rowidx 			 	 = rowidx_array[i] - low;
		x_push_local[rowidx] = Vx;
	}
	// Fill in Vy values
	iy=0;
	for (i=iy; i<numRowsy; i++)
	{
		rowidy 			 	 = rowidy_array[i] - low;
		x_push_local[rowidy] = Vy;
	}

	ierr = VecRestoreArray(x_push,&x_push_local); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Modify the PV_MAT by zeroing the columns */
/* PV part is taken care separately before the MatAssembly */
#undef __FUNCT__
#define __FUNCT__ "ZeroColumnsPV_MAT"
PetscErrorCode ZeroColumnsPV_MAT(UserContext *user,Vec pv_rhs_push,Mat PV_MAT,Mat PV_MAT_push,PetscInt *rowidx_array,PetscInt numRowsx,PetscInt *rowidy_array,PetscInt numRowsy,PetscInt ind_change)
{
	PetscErrorCode  ierr;
	PetscInt 		i,j,pv_start,pv_end,pvrow,pv_ncols,column_id;
	PetscScalar     Vx,Vy, theta,rhs_value,sum_pv_x,sum_pv_y;
	const PetscInt	*pv_cols;
	const PetscScalar *pv_vals;
	PetscInt        low, high;
	const PetscScalar zero=0.0;

	ierr = MatGetOwnershipRange(PV_MAT_push,&pv_start,&pv_end); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(pv_rhs_push,&low,&high); CHKERRQ(ierr);

	Vx = 0.0;
	Vy = 0.0;

	if (user->Pushing.dir[ind_change] == 0)
	{
		theta = user->Pushing.theta;
		Vx    = cos(theta)*user->Pushing.V_push[ind_change];
		Vy    = sin(theta)*user->Pushing.V_push[ind_change];
	}
	if (user->Pushing.dir[ind_change] == 1)
	{
		Vx    = user->Pushing.V_push[ind_change];
		Vy    = 0.0;
	}
	if (user->Pushing.dir[ind_change] == 2)
	{
		Vx    = 0.0;
		Vy    = user->Pushing.V_push[ind_change];
	}

	for (pvrow=pv_start; pvrow<pv_end; pvrow++){
		sum_pv_x = 0;
		sum_pv_y = 0;

		ierr = MatGetRow(PV_MAT_push,pvrow,&pv_ncols,&pv_cols,&pv_vals);CHKERRQ(ierr);
		for (j=0; j<pv_ncols; j++){
			for (i=0; i<numRowsx; i++){
				if (pv_cols[j]==rowidx_array[i]){
					column_id = rowidx_array[i];
					sum_pv_x += pv_vals[j];
					// Zero the corresponding column values
					ierr = MatSetValues(PV_MAT,1,&pvrow,1,&column_id,&zero,INSERT_VALUES);CHKERRQ(ierr);
				}}
			for (i=0; i<numRowsy; i++){
				if (pv_cols[j]==rowidy_array[i]){
					column_id = rowidy_array[i];
					sum_pv_y += pv_vals[j];
					// Zero the corresponding column values
					ierr = MatSetValues(PV_MAT,1,&pvrow,1,&column_id,&zero,INSERT_VALUES);CHKERRQ(ierr);
				}}
		}
		ierr = MatRestoreRow(PV_MAT_push,pvrow,&pv_ncols,&pv_cols,&pv_vals);CHKERRQ(ierr);
		//VecSetValues
		rhs_value = -Vx*sum_pv_x - Vy*sum_pv_y;
		ierr = VecSetValue(pv_rhs_push,pvrow,rhs_value,INSERT_VALUES);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
/* Modify RHS for pushing*/
#undef __FUNCT__
#define __FUNCT__ "AdvectThePushingBlockCoordinates"
PetscErrorCode AdvectThePushingBlockCoordinates(UserContext *user, PetscInt ind_change)
{
	PetscScalar 	x_c,y_c,z_c,Vx,Vy,dx,dy,dt,omega,theta,dtheta;
	PetscInt 		advect_coord;

	dy    = 0.0;
	dx    = 0.0;
	dt 	  = user->dt;
	advect_coord = user->Pushing.coord_advect[ind_change];

	x_c = user->Pushing.x_center_block;
	y_c = user->Pushing.y_center_block;
	z_c = user->Pushing.z_center_block;

	Vx = 0.0;
	Vy = 0.0;

	if (user->Pushing.dir[ind_change] == 0)
	{
		omega = user->Pushing.omega[ind_change];
		theta = user->Pushing.theta;
		Vx    = cos(theta)*user->Pushing.V_push[ind_change];
		Vy    = sin(theta)*user->Pushing.V_push[ind_change];

		dtheta = omega*dt;
		user->Pushing.theta	= theta + dtheta;
	}
	if (user->Pushing.dir[ind_change] == 1)
	{
		Vx    = user->Pushing.V_push[ind_change];
		Vy    = 0.0;
	}
	if (user->Pushing.dir[ind_change] == 2)
	{
		Vx    = 0.0;
		Vy    = user->Pushing.V_push[ind_change];
	}

	dx     = dt * Vx;
	dy     = dt * Vy;

	// Advect the coordinates
	if (advect_coord == 0){ 									//stationary block
		user->Pushing.x_center_block = x_c 	   ;
		user->Pushing.y_center_block = y_c 	   ;
		user->Pushing.z_center_block = z_c     ;
	}
	else if (advect_coord == 1){									// advect the block
		user->Pushing.x_center_block = x_c + dx;
		user->Pushing.y_center_block = y_c + dy;
		user->Pushing.z_center_block = z_c;
	}


	PetscFunctionReturn(0);
}
