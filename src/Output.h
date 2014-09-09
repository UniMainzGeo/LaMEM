
#ifndef __Output_h__
#define __Output_h__

PetscErrorCode WriteOutputFileMatlab( UserContext *user, DM da_nodes, DM da_temp, Vec Velocity, Vec Temp, PetscInt itime, const PetscInt _ElementType, const PetscInt _ngp_vel, const char DirectoryName[]);

#endif
