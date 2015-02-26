/* utilites for fdstag */
/* $Id: tools.h 5681 2015-02-20 20:57:42Z ltbaumann $ */

#ifndef __tools_h__
#define __tools_h__

// this are a couple of basic statistic functions
//---------------------------------------------------------------------------
PetscScalar getArthMean(PetscScalar *data, PetscInt n);
PetscScalar getVar(PetscScalar *data, PetscInt n);
PetscScalar getStdv(PetscScalar *data, PetscInt n);

#endif /* __tools_h__ */
