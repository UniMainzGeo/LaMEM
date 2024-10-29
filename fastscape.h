#ifndef _FASTSCAPELIB_H_
#define _FASTSCAPELIB_H_

PetscErrorCode fastscape(FreeSurf*, AdvCtx*);

#ifdef __cplusplus
extern "C"
{
#endif
    
    double *fastscapeFortran(int*,int*,double*,double*,double*,double*,int*,double*,double*);
    void clearArray();
#ifdef __cplusplus
}
#endif
/*
#define START_PLANE_LOOP_FS \
	for(j = 0; j < ny_fs; j++) \
	{	for(i = 0; i < nx_fs; i++) \
		{

// finalize plane access loop
#define END_PLANE_LOOP_FS \
		} \
	}
*/
#endif

