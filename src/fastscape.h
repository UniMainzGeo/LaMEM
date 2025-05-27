#ifndef _FASTSCAPELIB_H_
#define _FASTSCAPELIB_H_

#include "paraViewOutSurf.h"

struct FastScapeLib
{
	FreeSurf *surf;
	PVSurf *pvsurf;
	JacRes *jr;             // global residual context
	DBMat  *dbm;    // material database

	Vec vz_fs, vz_collect;
	Vec gtopo_fs;
	VecScatter ctx;

	PetscScalar kf;		// the bedrock river incision (SPL) rate parameter (or Kf) in meters (to the power 1-2m) per year
	PetscScalar kfsed;  // sediment river incision (SPL) rate parameter (or Kf) in meters (to the power 1-2m) per year; note that when kfsed < 0, 
						// its value is not used, i.e., kf for sediment and bedrock have the same value, regardless of sediment thickness
	PetscScalar m;		// drainage area exponent in the SPL
	PetscScalar n;		// slope exponent in the SPL
	PetscScalar kd;		// the bedrock transport coefficient (or diffusivity) for hillslope processes in meter squared per year
	PetscScalar kdsed;  // sediment transport coefficient (or diffusivity) for hillslope processes in meter squared per year, note that when kdsed < 0, 
						// its value is not used, i.e., kd for sediment and bedrock have the same value, regardless of sediment thickness
	PetscScalar g;		// bedrock dimensionless deposition/transport coefficient for the enriched SPL 
	PetscScalar gsed;   //sediment dimensionless deposition/transport coefficient for the enriched SPL, note that when gsed < 0, 
						//its value is not used, i.e., g for sediment and bedrock have the same value, regardless of sediment thickness   
	PetscScalar p;		// slope exponent for multi-direction flow; the distribution of flow among potential receivers 
						//(defined as the neighbouring nodes that define a negative slope)is proportional to local slope to power p

	PetscInt setMarine; // flag of using marine process
	PetscScalar sealevel; // sea level in meters
	PetscScalar poro_silt; // reference/surface porosity for silt
	PetscScalar poro_sand; // reference/surface porosity for sand
	PetscScalar zporo_silt; // e-folding depth for exponential porosity law for silt 
	PetscScalar zporo_sand; // e-folding depth for exponential porosity law for sand
	PetscScalar ratio; // silt fraction for material leaving the continent
	PetscScalar Lsolve; // averaging depth/thickness needed to solve the silt-sand equation in meters
	PetscScalar kds_silt; // marine transport coefficient (diffusivity) for silt in meters squared per year
	PetscScalar kds_sand; // marine transport coefficient (diffusivity) for sand in meters squared per year

	PetscScalar rangeX; // range in x-direction
	PetscScalar rangeY; // range in y-direction
	PetscInt    refine;	// whether refine the grid in FastScape
	PetscInt    nx_refine; // nodes in x-direction after refinement
	PetscInt    ny_refine; // nodes in y-direction after refinement

	PetscScalar Max_dt; // max dt used in FastScape
	PetscInt    FS_BC;  // boundary condition in FastScape
	PetscInt    sedPhases;   // sediment layers phase numbers

	PetscInt	bgphase; // phase of background
};

PetscErrorCode FastScapeCreate(FastScapeLib*, FB*);

PetscErrorCode fastscape(FastScapeLib*);
PetscErrorCode bilinearInterpolate(FastScapeLib*, PetscScalar, PetscScalar, Scaling*, PetscInt);
PetscErrorCode savePvtsFS(PVSurf*, PetscInt, PetscInt, PetscScalar, PetscScalar, PetscScalar, PetscInt, const char*, PetscScalar*);

#ifdef __cplusplus
extern "C"
{
#endif
    
    double *fastscapeFortran(int*,int*,double*,double*,double*,double*,int*,double*,double*,
    	double*,double*,double*,double*,double*,double*,double*, double*, double*, int*, 
    	double*,double*,double*,double*,double*,double*,double*, double*, double*, int*);
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


