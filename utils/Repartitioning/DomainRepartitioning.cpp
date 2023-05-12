
#include "DomainRepartitioning.h"


static char help[] = "***********************************************************************\n"
		"Run as ./Repartition -timestep X \n Repartition of particle input files .\n"
		"***********************************************************************\n";



/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "LoadPartitioning"
PetscErrorCode LoadPartitioning(Partition &part,const char *filename)
{
/*
 * 	Reads processor partitioning file
 * 	INPUT:
 * 	- pointer to partitioning structure
 *
 * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
 */

	PetscErrorCode		ierr;
	PetscInt            fileID;

	PetscFunctionBeginUser;

    ierr = PetscBinaryOpen(filename,FILE_MODE_READ,&fileID);      CHKERRQ(ierr);

    ierr = PetscBinaryRead(fileID,&part.Nx,1,PETSC_INT);          CHKERRQ(ierr);
    ierr = PetscBinaryRead(fileID,&part.Ny,1,PETSC_INT);          CHKERRQ(ierr);
    ierr = PetscBinaryRead(fileID,&part.Nz,1,PETSC_INT);          CHKERRQ(ierr);

    // Allocate memory for three vectors x,y,z
    ierr = PetscMalloc((part.Nx+1)*sizeof(PetscScalar),&part.px); CHKERRQ(ierr);
    ierr = PetscMalloc((part.Ny+1)*sizeof(PetscScalar),&part.py); CHKERRQ(ierr);
    ierr = PetscMalloc((part.Nz+1)*sizeof(PetscScalar),&part.pz); CHKERRQ(ierr);

    ierr = PetscBinaryRead(fileID,part.px,part.Nx+1,PETSC_SCALAR);CHKERRQ(ierr);
    ierr = PetscBinaryRead(fileID,part.py,part.Ny+1,PETSC_SCALAR);CHKERRQ(ierr);
    ierr = PetscBinaryRead(fileID,part.pz,part.Nz+1,PETSC_SCALAR);CHKERRQ(ierr);

    ierr = PetscBinaryClose(fileID);                              CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}
/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "DestroyPartitioning"
PetscErrorCode DestroyPartitioning(Partition &part)
{
/*
 * 	Deallocates partitioning structure
 * 	INPUT:
 * 	- pointer to partitioning structure
 *
 * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
 */
	PetscErrorCode		ierr;

	PetscFunctionBeginUser;

	ierr = PetscFree(part.px);                                    CHKERRQ(ierr);
	ierr = PetscFree(part.py);                                    CHKERRQ(ierr);
	ierr = PetscFree(part.pz);                                    CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}
/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "LoadParticles"
PetscErrorCode LoadParticles(PetscInt size,PetscInt itime,Particles &particles)
{
/*  Loads particles to structure
 *
 * 	INPUT:
 * 	- pointer to Particles structure
 *
 * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
 */
	PetscErrorCode		ierr=0;
	char                fname[50];
	PetscViewer         view_in;
	Vec                 lvec_info,lvec_prtcls;

	PetscScalar        *larray_info,*larray_prtcls;
	PetscInt            rank,Nprops=0,Nprocs=0,Nprtcls=0,Ntot_prtcls=0,n=0,k;

	PetscFunctionBeginUser;

	// Read info vectors only
	for(rank=0;rank<size;rank++){
		sprintf(fname,"Particles.%d.%d.out",rank,itime+1000000);

		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fname,FILE_MODE_READ,&view_in);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_info);             CHKERRQ(ierr);
		ierr = VecLoad(lvec_info,view_in);                        CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in);                      CHKERRQ(ierr);

		ierr = VecGetArray(lvec_info,&larray_info);               CHKERRQ(ierr);
		Ntot_prtcls = Ntot_prtcls + (PetscInt) larray_info[0];
		ierr = VecRestoreArray(lvec_info,&larray_info);           CHKERRQ(ierr);

		ierr = VecDestroy(&lvec_info);                            CHKERRQ(ierr);
	}

	// Allocate structure, which holds all particles
	particles.N  = Ntot_prtcls;
	particles.x  = new PetscScalar[particles.N];
	particles.y  = new PetscScalar[particles.N];
	particles.z  = new PetscScalar[particles.N];
	particles.phs= new unsigned char [particles.N];

	// Read all
	for(rank=0;rank<size;rank++){

		sprintf(fname,"Particles.%d.%d.out",rank,itime+1000000);

		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fname,FILE_MODE_READ,&view_in);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_info);             CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_prtcls);           CHKERRQ(ierr);
		ierr = VecLoad(lvec_info,view_in);                        CHKERRQ(ierr);
		ierr = VecLoad(lvec_prtcls,view_in);                      CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in);                      CHKERRQ(ierr);

		ierr = VecGetArray(lvec_info,&larray_info);               CHKERRQ(ierr);
		Nprtcls = (PetscInt) larray_info[0];
		Nprocs  = (PetscInt) larray_info[1];
		if(size != Nprocs)
			PetscPrintf(PETSC_COMM_WORLD,"[Error] ... while loading particles: size != Nprocs");
		Nprops  = (PetscInt) larray_info[2];
		ierr = VecRestoreArray(lvec_info,&larray_info);           CHKERRQ(ierr);

		ierr = VecGetArray(lvec_prtcls,&larray_prtcls);           CHKERRQ(ierr);
		for(k=0;k<Nprtcls;k++){
			// x,y,z,num,phase,cpu, ... (28 in total)
			particles.x[n]   = larray_prtcls[(k*Nprops)+0];
			particles.y[n]   = larray_prtcls[(k*Nprops)+1];
			particles.z[n]   = larray_prtcls[(k*Nprops)+2];
			particles.phs[n] = (char unsigned)larray_prtcls[(k*Nprops)+4];
			n++;
		}
		ierr = VecRestoreArray(lvec_prtcls,&larray_prtcls);       CHKERRQ(ierr);

		ierr = VecDestroy(&lvec_info);                            CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_prtcls);                          CHKERRQ(ierr);
	}


	PetscFunctionReturn(ierr);
}
/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "SaveParticles"
PetscErrorCode SaveParticles(Partition target,Particles *particles)
{
/*  Loads particles to structure
 *
 * 	INPUT:
 * 	- pointer to Particles structure
 *
 * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
 */
	PetscErrorCode		ierr = 0;
	char                fname[PETSC_MAX_PATH_LEN];
	PetscViewer         view_out;
	Vec                 lvec_info,lvec_prtcls;
	PetscScalar        *larray_prtcls;
	PetscInt            rank,Nprops=28,k,n=0,num;

	PetscFunctionBeginUser;

	PetscPrintf(PETSC_COMM_WORLD,"Write %d InitialParticles files \n",target.Ntot);
	for(rank=0;rank<target.Ntot;rank++){

		// Create info vector
		ierr = VecCreateSeq(PETSC_COMM_SELF,10,&lvec_info);       CHKERRQ(ierr);
		ierr = VecSetValue(lvec_info,0,target.Ntot,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecSetValue(lvec_info,1,particles[rank].N,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecSetValue(lvec_info,2,Nprops,INSERT_VALUES);     CHKERRQ(ierr);
		ierr = VecAssemblyBegin(lvec_info);                       CHKERRQ(ierr);
		ierr = VecAssemblyEnd(lvec_info);                         CHKERRQ(ierr);

		// Create vector that holds particles
	    ierr = VecCreateSeq(PETSC_COMM_SELF,particles[rank].N*Nprops,&lvec_prtcls);CHKERRQ(ierr);
	    ierr = VecGetArray(lvec_prtcls,&larray_prtcls);           CHKERRQ(ierr);
		n = 0;
	    for(k=0;k<particles[rank].N;k++){
			// x,y,z,num,phase,cpu, ... (28 in total)
			larray_prtcls[(n*Nprops)+0] = particles[rank].x[k];
			larray_prtcls[(n*Nprops)+1] = particles[rank].y[k];
			larray_prtcls[(n*Nprops)+2] = particles[rank].z[k];
			larray_prtcls[(n*Nprops)+3] = k;
			larray_prtcls[(n*Nprops)+4] = (PetscScalar) particles[rank].phs[k];
			n++;
		}
		ierr = VecRestoreArray(lvec_prtcls,&larray_prtcls);       CHKERRQ(ierr);

	    // Save file
		num = mkdir("InitialParticles", S_IRWXU);
		if(num)
			PetscPrintf(PETSC_COMM_WORLD,"Write output to existing directory: InitialParticles, %d particles \n",n);
		else
			PetscPrintf(PETSC_COMM_WORLD,"Write output to created  directory: InitialParticles, %d particles \n",n);

		sprintf(fname,"InitialParticles/Initial_Particles.%d.out",rank);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fname,FILE_MODE_WRITE,&view_out);CHKERRQ(ierr);
		ierr = VecView(lvec_info,view_out);                       CHKERRQ(ierr);
		ierr = VecView(lvec_prtcls,view_out);                     CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out);                     CHKERRQ(ierr);

		// clean up vectors
		VecDestroy(&lvec_prtcls);                                 CHKERRQ(ierr);
		VecDestroy(&lvec_info);                                   CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}
/*==========================================================================================================*/
/* Main routine */
#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc,char **argv )
{/*
  * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
  */
	PetscErrorCode 	ierr;
    PetscInt        n,k,itime;
    Partition       init,target;
    Particles       all_prtcls, *sub_prtcls;


	// --- Initialize PETSC ---
	ierr = PetscInitialize(&argc,&argv,(char *)0,help);           CHKERRQ(ierr);
	
	itime = 0;
	PetscOptionsGetInt(PETSC_NULL ,"-timestep",	&itime	, PETSC_NULL);

    // Load initial partitioning file
    PetscPrintf(PETSC_COMM_WORLD,"> Load initial partitioning file \n");
    ierr = LoadPartitioning(init,"ProcessorPartitioning_init.bin");CHKERRQ(ierr);
    init.Ntot = init.Nx*init.Ny*init.Nz;

    // Load target partitioning file
    PetscPrintf(PETSC_COMM_WORLD,"> Load target partitioning file \n");
    ierr = LoadPartitioning(target,"ProcessorPartitioning_target.bin");CHKERRQ(ierr);
    target.Ntot = target.Nx*target.Ny*target.Nz;

    // Load initial particles
    PetscPrintf(PETSC_COMM_WORLD,"> Load particles\n");
    ierr = LoadParticles(init.Ntot,itime,all_prtcls);           CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"> Loaded %d particles\n",all_prtcls.N);

    // Allocate memory for target files
    sub_prtcls = new Particles[target.Ntot];
    for (k=0;k<target.Ntot;k++){
    	sub_prtcls[k].N = 0;
    	sub_prtcls[k].m = 0;
    }

    // Sorting to target partitioning using counter
    for(k=0;k<all_prtcls.N;k++){
		n = getGlobalRank(Bisection(target.px,target.Nx,all_prtcls.x[k]),
	                      Bisection(target.py,target.Ny,all_prtcls.y[k]),
	                      Bisection(target.pz,target.Nz,all_prtcls.z[k]),
	                      target.Nx,target.Ny,target.Nz);

		sub_prtcls[n].N++;
	}

	// Allocate memory for target file arrays
    for (k=0;k<target.Ntot;k++){
    	sub_prtcls[k].phs= new unsigned char[sub_prtcls[k].N];
    	sub_prtcls[k].x  = new PetscScalar[sub_prtcls[k].N];
    	sub_prtcls[k].y  = new PetscScalar[sub_prtcls[k].N];
    	sub_prtcls[k].z  = new PetscScalar[sub_prtcls[k].N];
    }

    // Transferring particles from initial partition to target partition
    for(k=0;k<all_prtcls.N;k++){
		n = getGlobalRank(Bisection(target.px,target.Nx,all_prtcls.x[k]),
	                      Bisection(target.py,target.Ny,all_prtcls.y[k]),
	                      Bisection(target.pz,target.Nz,all_prtcls.z[k]),
	                      target.Nx,target.Ny,target.Nz);

		sub_prtcls[n].x[sub_prtcls[n].m]   = all_prtcls.x[k];
		sub_prtcls[n].y[sub_prtcls[n].m]   = all_prtcls.y[k];
		sub_prtcls[n].z[sub_prtcls[n].m]   = all_prtcls.z[k];
		sub_prtcls[n].phs[sub_prtcls[n].m] = all_prtcls.phs[k];

		sub_prtcls[n].m++;
	}

    // save target partion
   	ierr = SaveParticles(target,sub_prtcls);                      CHKERRQ(ierr);


    // Deallocate particles structures
    for (k=0;k<target.Ntot;k++){
    	delete [] sub_prtcls[k].phs;
    	delete [] sub_prtcls[k].x;
    	delete [] sub_prtcls[k].y;
    	delete [] sub_prtcls[k].z;
    }
    delete [] sub_prtcls;

    // Deallocate partitioning structures
    ierr = DestroyPartitioning(init);                             CHKERRQ(ierr);
    ierr = DestroyPartitioning(target);                           CHKERRQ(ierr);

	// --- cleanup PETSC ---
	ierr = PetscFinalize();                                       CHKERRQ(ierr);

	return 0;
}
/* End of code */
/*==========================================================================================================*/

