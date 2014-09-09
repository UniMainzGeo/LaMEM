

#include "petsc.h"
#include "petscda.h"


PetscErrorCode DAVerifyDecomposition( 
		const PetscInt M, const PetscInt N, const PetscInt P,
		const PetscInt size,
		const PetscInt stencil_width,
		DAPeriodicType wrap,DAStencilType stencil_type, 
		PetscInt *_m, PetscInt *_n, PetscInt *_p,
		PetscTruth *pass )
{
	PetscInt m,n,p;
	PetscInt pm;
	
	*pass = PETSC_TRUE;
	
	/* try for squarish distribution */
	n = (PetscInt)(0.5 + pow(((PetscReal)N*N)*((PetscReal)size)/((PetscReal)P*M),(PetscReal)(1./3.)));
	if (!n) n = 1;
	while (n > 0) {
		pm = size/n;
		if (n*pm == size) break;
		n--;
	}   
	if (!n) n = 1; 
	m = (PetscInt)(0.5 + sqrt(((PetscReal)M)*((PetscReal)size)/((PetscReal)P*n)));
	if (!m) m = 1;
	while (m > 0) {
		p = size/(m*n);
		if (m*n*p == size) break;
		m--;
	}
	if (M > P && m < p) {PetscInt _m = m; m = p; p = _m;}
	
	if (m*n*p != size) {
		PetscPrintf(PETSC_COMM_WORLD, "\t\tCould not find good partition");  
		*pass = PETSC_FALSE;
	}
	if (M < m) {
		PetscPrintf(PETSC_COMM_WORLD,"\t\tPartition in x direction is too fine! %D %D\n",M,m);
		*pass = PETSC_FALSE;
	}
	if (N < n) {
		PetscPrintf(PETSC_COMM_WORLD,"\t\tPartition in y direction is too fine! %D %D\n",N,n);
		*pass = PETSC_FALSE;
	}
	if (P < p) {
		PetscPrintf(PETSC_COMM_WORLD,"\t\tPartition in z direction is too fine! %D %D\n",P,p);
		*pass = PETSC_FALSE;
	}
	
	*_m = m;
	*_n = n;
	*_p = p;
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
	PetscErrorCode ierr;
	PetscInt nlevels, nrx,nry,nrz;
	PetscTruth coarsen;
	PetscInt k;
	PetscInt Mi,Ni,Pi, Mk,Nk,Pk;
	PetscInt nprocs;
	DAPeriodicType wrap;
	PetscInt m,n,p;
	
	PetscInitialize(&argc,&args,(char *)0,0);
	
	
	
	nlevels = 1;
	PetscOptionsGetInt( PETSC_NULL, "-nlevels", &nlevels, 0 );
	
	nrx = nry = nrz = 2;
	PetscOptionsGetInt( PETSC_NULL, "-nrx", &nrx, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-nry", &nry, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-nrz", &nrz, 0 );
	
	Mi = Ni = Pi = 2;
	PetscOptionsGetInt( PETSC_NULL, "-M", &Mi, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-N", &Ni, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-P", &Pi, 0 );
	
	nprocs = 1;
	PetscOptionsGetInt( PETSC_NULL, "-nprocs", &nprocs, 0 );
	
	coarsen = PETSC_FALSE;
	PetscOptionsGetTruth( PETSC_NULL, "-coarsen", &coarsen, 0 );
	
	wrap = DA_NONPERIODIC;
	Mk = Mi;
	Nk = Ni;
	Pk = Pi;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "Verifying DA decomposition\n");
	PetscPrintf( PETSC_COMM_WORLD, "  npus=%d    (-nprocs)\n", nprocs );
	PetscPrintf( PETSC_COMM_WORLD, "  refinement_x=%d    (-nrx)\n", nrx );
	PetscPrintf( PETSC_COMM_WORLD, "  refinement_y=%d    (-nry)\n", nry );
	PetscPrintf( PETSC_COMM_WORLD, "  refinement_z=%d    (-nrz)\n", nrz );
	PetscPrintf( PETSC_COMM_WORLD, "  nlevels=%d    (-nlevels)\n", nlevels );
	PetscPrintf( PETSC_COMM_WORLD, "  M_initial=%d    (-M)\n", Mi );
	PetscPrintf( PETSC_COMM_WORLD, "  N_initial=%d    (-N)\n", Ni );
	PetscPrintf( PETSC_COMM_WORLD, "  P_initial=%d    (-P)\n", Pi );
	
	if( coarsen == PETSC_TRUE ){
		PetscPrintf( PETSC_COMM_WORLD, "Applying coarsening to initial grid\n" );
		for( k=0; k<nlevels; k++ ) {
			PetscTruth pass = PETSC_FALSE;
			
			PetscPrintf( PETSC_COMM_WORLD, "Level[%d]: \n", k );
			PetscPrintf( PETSC_COMM_WORLD, "\tM=%d : N=%d : P=%d \n", Mk,Nk,Pk );
			DAVerifyDecomposition(Mk,Nk,Pk, nprocs,0, wrap, DA_STENCIL_STAR, &m,&n,&p, &pass ); 
			PetscPrintf( PETSC_COMM_WORLD, "\tcpu_x=%d : cpu_y=%d : cpu_z=%d \n", m,n,p );
			if(pass==PETSC_FALSE) {
				PetscPrintf( PETSC_COMM_WORLD, "\tDecomposition is not valid \n\n");
			}
			else {
				PetscPrintf( PETSC_COMM_WORLD, "\tDecomposition is valid \n\n" );
			}
			
			
			if (DAXPeriodic(wrap) ){
				if(nrx)
						Mk = Mk / nrx;
				else
						Mk = Mk;
			} else {
				if(nrx)
						Mk = 1 + (Mk - 1) / nrx;
				else
						Mk = Mk;
			}
			if (DAYPeriodic(wrap) ){
				if(nry)
						Nk = Nk / nry;
				else
						Nk = Nk;
			} else {
				if(nry)
						Nk = 1 + (Nk - 1) / nry;
				else
						Nk = Nk;
			}
			if (DAZPeriodic(wrap) ){
				if(nrz)
						Pk = Pk / nrz;
				else
						Pk = Pk;
			} else {
				if(nrz)
						Pk = 1 + (Pk - 1) / nrz;
				else
						Pk = Pk;
			}
		}
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "Refining initial grid\n" );
		for( k=0; k<nlevels; k++ ) {
			PetscTruth pass = PETSC_FALSE;
			
			PetscPrintf( PETSC_COMM_WORLD, "Level[%d]: \n", k );
			PetscPrintf( PETSC_COMM_WORLD, "\tM=%d : N=%d : P=%d \n", Mk,Nk,Pk );
			DAVerifyDecomposition(Mk,Nk,Pk, nprocs,0, wrap, DA_STENCIL_STAR, &m,&n,&p, &pass ); 
			PetscPrintf( PETSC_COMM_WORLD, "\tcpu_x=%d : cpu_y=%d : cpu_z=%d \n", m,n,p );
			if(pass==PETSC_FALSE) {
				PetscPrintf( PETSC_COMM_WORLD, "\tDecomposition is not valid \n\n");
			}
			else {
				PetscPrintf( PETSC_COMM_WORLD, "\tDecomposition is valid \n\n" );
			}
			
			if (DAXPeriodic(wrap) ){
				Mk = nrx*Mk;
			} else {
				Mk = 1 + nrx*(Mk - 1);
			}
			if (DAYPeriodic(wrap) ){
				Nk = nry*Nk;
			} else {
				Nk = 1 + nry*(Nk - 1);
			}
			if (DAZPeriodic(wrap) ){
				Pk = nrz*Pk;
			} else {
				Pk = 1 + nrz*(Pk - 1);
			}
		}
	}
	
	
	
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}



