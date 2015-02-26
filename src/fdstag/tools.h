/* utilites for fdstag */
/* $Id: tools.h 5681 2015-02-20 20:57:42Z ltbaumann $ */

#ifndef __tools_h__
#define __tools_h__


//-----------------------------------------------------------------------------
// INLINE FUNCTIONS
//-----------------------------------------------------------------------------

// this are a couple of basic statistic functions
//---------------------------------------------------------------------------
static inline PetscScalar getArthMean(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar sum = 0.0;

    for (k=0; k<n; k++)
        sum += data[k];
    return sum/(PetscScalar)n;
}
//---------------------------------------------------------------------------
static inline PetscScalar getVar(PetscScalar *data, PetscInt n)
{
	PetscInt    k;
    PetscScalar mean = getArthMean(data,n);
    PetscScalar temp = 0.0;

    for (k=0; k<n; k++)
        temp += (mean-data[k])*(mean-data[k]);
    return temp/(PetscScalar)n;
}
//---------------------------------------------------------------------------
static inline PetscScalar getStdv(PetscScalar *data, PetscInt n)
{
    return sqrt(getVar(data,n));
}

#endif /* __tools_h__ */
