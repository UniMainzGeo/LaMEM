//---------------------------------------------------------------------------
//..........................  FDSTAG UTILITIES ..............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "tools.h"
//---------------------------------------------------------------------------


// this are a couple of basic statistic functions
//---------------------------------------------------------------------------
PetscScalar getArthMean(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar sum = 0.0;

    for (k=0; k<n; k++)
        sum += data[k];
    return (sum/(PetscScalar)n);
}
//---------------------------------------------------------------------------
PetscScalar getVar(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar mean = getArthMean(data,n);
    PetscScalar temp = 0.0;
	PetscPrintf(PETSC_COMM_WORLD,"mean=%g \n",mean);
    for (k=0; k<n; k++)
        temp += (mean-data[k])*(mean-data[k]);
    return (temp/(PetscScalar)n);
}
//---------------------------------------------------------------------------
PetscScalar getStdv(PetscScalar *data, PetscInt n)
{

	PetscPrintf(PETSC_COMM_WORLD,"data[0..%lld]=%g,%g ,%g ,%g \n",(LLD) n, data[0],data[1],data[2],data[3], data[n-1]);
    return sqrt(getVar(data,n));
}

