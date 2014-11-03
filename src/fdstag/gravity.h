//---------------------------------------------------------------------------
//.........................   GRAVITY FIELD    ..............................
//---------------------------------------------------------------------------
#ifndef __gravity_h__
#define __gravity_h__

//---------------------------------------------------------------------------
// survey context
typedef struct
{
	PetscInt     i,j,nx,ny;
	PetscInt     xs,xm,ys,ym;
	PetscInt     iter;
	PetscScalar  x,y,z,dx,dy;
	Vec          lvec_dg,lvec_dg2save,gvec_dg;
	PetscScalar *coord,*dg;
	PetscMPIInt  rank;
} GravitySurvey;
//---------------------------------------------------------------------------


#endif
