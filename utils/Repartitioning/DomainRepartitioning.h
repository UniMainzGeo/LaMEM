#ifndef __DomainRepartitioning_h__
#define __DomainRepartitioning_h__

//-----------------------------------------------------------------------------
// EXTERNAL INCLUDES
//-----------------------------------------------------------------------------

#define _GNU_SOURCE
#include "stdio.h"
//#include "stdlib.h"
//#include "string.h"
#include "petsc.h"
#include "sys/stat.h"
//---------------------------------------------------------------------------
// processor partitioning structure
typedef struct
{
	PetscInt       Nx,Ny,Nz,Ntot;
	PetscScalar    *px,*py,*pz;

} Partition;
//---------------------------------------------------------------------------
// structure holding the particles
typedef struct
{
	PetscInt       N,m;
	PetscScalar    *x,*y,*z;
	unsigned char  *phs;

} Particles;

//---------------------------------------------------------------------------
// subfunctions functions
// this function returns global rank of processor in DMDA
static inline PetscMPIInt getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p)
{
	if (i < 0 || i >= m || j < 0 || j >= n || k < 0 || k >= p) return -1;
	return (PetscMPIInt)(i + j*m + k*m*n);
}

//-----------------------------------------------------------------------------
template <class T>
PetscInt Bisection(T *px,PetscInt Nx,PetscScalar x){

	PetscInt      L,R,M;

	L = 0;
	R = Nx;

	while((R-L)>1){
		M = (L+R)/2;
		if(px[M]<= x)
			L=M;
		if(px[M]>= x)
			R=M;
	}

	return(L);
}
//---------------------------------------------------------------------------

#define LLD long long int


#endif
