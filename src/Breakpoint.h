

#ifndef __Breakpoint_h__
#define __Breakpoint_h__


PetscErrorCode SaveBreakPoint( UserContext *user, DM da_nodes, DM da_pres, Vec Velocity,
		Vec Temp, Vec Pressure, PetscInt itime, PetscInt FileNumber );


PetscErrorCode LoadBreakPoint(LaMEMVelPressureDA C, UserContext *user, DM da_nodes, Vec Velocity,
		Vec Temp, Vec Pressure, PetscInt FileNumber );

PetscErrorCode LoadBreakPoint_old( UserContext *user, DM da_nodes, DM da_pres, DM da_temp, Vec Velocity,
		Vec Temp, Vec Pressure, PetscInt FileNumber );

#endif
