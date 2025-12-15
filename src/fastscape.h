#ifndef _FASTSCAPELIB_H_
#define _FASTSCAPELIB_H_

#define NO_NEED -1

#include "paraViewOutSurf.h"

struct JacRes;
// for gather information 
typedef struct 
{
    PetscInt sx, sy, sz;
    PetscInt nx, ny, nz;
    PetscMPIInt rankZ_id;
} ProcInfo;
// for interploation
struct GridIndex 
{
  PetscInt m;
  PetscInt mm;
  PetscInt n;
  PetscInt nn;
  bool found;
};
struct MeshSeg1DFS
{
	PetscInt    nsegs;                    // number of segments
	PetscInt    istart[_max_num_segs_+1]; // indices of the first nodes plus last index
	PetscScalar xstart[_max_num_segs_+1]; // coordinates of the first nodes plus total size
	PetscScalar biases[_max_num_segs_  ]; // biases for each segment
	PetscInt    tcels;                    // total number of cells
	PetscInt    bias;                  // uniform grid flag
	PetscInt    periodic;                 // periodic topology flag
  PetscScalar grid_spacing_min[_max_num_segs_  ];
  PetscScalar grid_spacing_max[_max_num_segs_  ];  
  PetscInt    nnodes_nug;
  PetscScalar *ncoor_nug;  
  PetscScalar min_spacing;
  PetscScalar max_spacing;

};
struct FSGrid
{
  PetscInt nodes;
  PetscInt nodes_extend;
  PetscInt nodes_refine;
  PetscScalar dx;
  PetscScalar dx_refine;
  PetscScalar dx_extend;
  PetscScalar *ncoor;
  PetscScalar *ncoor_extend;
  PetscScalar *ncoor_refine;
};

struct FastScapeLib
{
	FreeSurf *surf;
	PVSurf   *pvsurf;
	JacRes   *jr;         
  Scaling  *scal;
  ProcInfo *proc_info;

  // grid information
  // non uniform grid
  PetscInt    non_uniform_grid;
  MeshSeg1DFS   msx_fs;
  MeshSeg1DFS   msy_fs;
  Vec gtopo_nug, vx_nug, vy_nug, vz_nug;
  FSGrid  fsX;
  FSGrid  fsY;
  PetscScalar *ncoor_ori_x, *ncoor_ori_y;

  // LaMEM grid
  PetscInt    nx_LaMEM;
  PetscInt    ny_LaMEM;
  PetscScalar rangeX_begin;
  PetscScalar rangeX_end;
  PetscScalar rangeY_begin;
  PetscScalar rangeY_end; 
	PetscScalar rangeX; // range in x-direction
	PetscScalar rangeY; // range in y-direction
	PetscScalar rangeZ_begin;
	PetscScalar rangeZ_end;
	Vec vx_LaMEM, vy_LaMEM, vz_LaMEM, gtopo_fs, gtopo_collect;
  Vec vx_collect, vy_collect, vz_collect;
  VecScatter ctx;

  // refined LaMEM grid
  PetscInt    refine; // whether refine the grid in FastScape
	PetscInt    nx_refine; // nodes in x-direction after refinement
	PetscInt    ny_refine; // nodes in y-direction after refinement
	Vec gtopo_refine, vx_refine, vy_refine, vz_refine;

  // 2D grid
	PetscInt 	  fs2D;	 // flag of 2D model
	PetscScalar extendedRange; // extended range in x or y direction
	PetscInt    extendedNodes; // nodes in x or y direction after extending
	PetscInt    extendedXNodes;
	PetscScalar extendedXRange;
	PetscInt    extendedYNodes; 
	PetscScalar extendedYRange;
	
  PetscInt    extendedX;
	PetscInt    extendedY;
	Vec gtopo_extend, vx_extend, vy_extend, vz_extend;

  // refined 2D grid
	PetscInt    etRefineXNodes;
	PetscInt    etRefineYNodes;
	Vec gtopo_et_refine, vx_et_refine, vy_et_refine, vz_et_refine;

  // max grid
  PetscInt    nx_solve;
  PetscInt    ny_solve;
  PetscInt    nodes_solve;

  // surface parameter
  PetscScalar Max_dt; // max dt used in FastScape
  char    FS_BC[_str_len_];  // topography boundary condition in FastScape
  char    FS_VELBC[_str_len_]; // velocity boundary condition in FastScape
  PetscInt    random_noise; // random noise flag for topography

  // erosion
  PetscScalar kf;   // the bedrock river incision (SPL) rate parameter (or Kf) in meters (to the power 1-2m) per year
  PetscScalar kfsed;  // sediment river incision (SPL) rate parameter (or Kf) in meters (to the power 1-2m) per year; note that when kfsed < 0, 
            // its value is not used, i.e., kf for sediment and bedrock have the same value, regardless of sediment thickness
  PetscScalar m;    // drainage area exponent in the SPL
  PetscScalar n;    // slope exponent in the SPL
  PetscScalar kd;   // the bedrock transport coefficient (or diffusivity) for hillslope processes in meter squared per year
  PetscScalar kdsed;  // sediment transport coefficient (or diffusivity) for hillslope processes in meter squared per year, note that when kdsed < 0, 
            // its value is not used, i.e., kd for sediment and bedrock have the same value, regardless of sediment thickness
  PetscScalar g;    // bedrock dimensionless deposition/transport coefficient for the enriched SPL 
  PetscScalar gsed;   //sediment dimensionless deposition/transport coefficient for the enriched SPL, note that when gsed < 0, 
            //its value is not used, i.e., g for sediment and bedrock have the same value, regardless of sediment thickness   
  PetscScalar p;    // slope exponent for multi-direction flow; the distribution of flow among potential receivers 
            //(defined as the neighbouring nodes that define a negative slope)is proportional to local slope to power p
  // sedimentaion
  PetscInt    setMarine; // flag of using marine process
  PetscInt    sedPhases;   // sediment layers phase numbers
  PetscScalar sealevel; // sea level in meters
  PetscScalar poro_silt; // reference/surface porosity for silt
  PetscScalar poro_sand; // reference/surface porosity for sand
  PetscScalar zporo_silt; // e-folding depth for exponential porosity law for silt 
  PetscScalar zporo_sand; // e-folding depth for exponential porosity law for sand
  PetscScalar ratio; // silt fraction for material leaving the continent
  PetscScalar Lsolve; // averaging depth/thickness needed to solve the silt-sand equation in meters
  PetscScalar kds_silt; // marine transport coefficient (diffusivity) for silt in meters squared per year
  PetscScalar kds_sand; // marine transport coefficient (diffusivity) for sand in meters squared per year


  // output parameter
  PetscInt    surf_out_nstep;
  PetscScalar vec_times;
  Vec silt_fraction, basement, total_erosion, drainage_area;
  Vec erosion_rate, slope, curvature, chi, catchment, lake_depth;
  // output
  char       outfile_fs[_str_len_+20]; // output file name
  char       outfile_refine[_str_len_+20]; // output file name
  char       outfile_extend[_str_len_+20]; // output file name
  char       outfile_et_refine[_str_len_+20]; 

  float     *buff_fs;            // direct output buffer
  long int   offset_fs;             // pvd file offset
  PetscInt   outsurf_fs;            // free surface output flag
  PetscInt   outpvd_fs;             // pvd file output flag
  PetscInt   out_topofs;         // refined topography output flag
  PetscInt   out_silt_fraction;
  PetscInt   out_basement;
  PetscInt   out_total_erosion;
  PetscInt   out_drainage_area;
  PetscInt   out_erosion_rate;
  PetscInt   out_slope;
  PetscInt   out_curvature;
  PetscInt   out_chi;
  PetscInt   out_catchment;
  PetscInt   out_lake_depth;
};

PetscErrorCode FastScapeCreate(FastScapeLib*, FB*);
PetscErrorCode FastScapeCreateData(FastScapeLib*);
PetscErrorCode FastScapeLoadGridInf(FastScapeLib*);
PetscErrorCode FSLoadNonUniformGrid(MeshSeg1DFS*, PetscScalar, Scaling*);
void           GenerateGridCoordinates(PetscScalar*, PetscInt, PetscScalar, PetscScalar);
PetscErrorCode FastScapeCreateSurfaceGrid(FastScapeLib*, PetscInt);
PetscErrorCode ScalingFastScapeCreate(Scaling*);
PetscErrorCode PVSurfFastScapeCreate(FastScapeLib*, FB*);
PetscErrorCode FastScapeCopyVelocity(FastScapeLib*);
PetscErrorCode FastScapeCreateGlobalGrid(PetscScalar*, MeshSeg1DFS, PetscInt, Scaling*);
GridIndex      FindIndexForInterpolation(FastScapeLib*, PetscScalar, PetscScalar, PetscScalar, PetscScalar);
PetscScalar    ReturnBiInterFunction(PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscInt, PetscInt, PetscInt, PetscInt);
PetscErrorCode InterpolationFor3DNonUniformGrid(FastScapeLib*, PetscScalar*, PetscInt);
PetscErrorCode InterpolationFor2DNonUniformGrid(FastScapeLib*, PetscScalar*, PetscScalar*);
PetscErrorCode GatherVariableFromLaMEM(FastScapeLib*, PetscScalar*, PetscScalar*, PetscScalar*, PetscScalar*, PetscInt);
PetscErrorCode FastScapeStretchGrid(FastScapeLib*);
PetscErrorCode BilinearInterpolate(FastScapeLib*, PetscScalar*, PetscScalar*, PetscScalar*, Scaling*, PetscInt, PetscInt, PetscInt);
PetscErrorCode Extended2D(FastScapeLib*, PetscScalar*, PetscScalar*, PetscScalar*, Scaling*, PetscInt, PetscInt, PetscInt);
PetscErrorCode FastScapeRun(FastScapeLib*);
PetscErrorCode PassValue2D(FastScapeLib*, PetscScalar*, PetscScalar*);
PetscErrorCode PassValue3D(FastScapeLib*, PetscScalar*, PetscScalar*);
PetscErrorCode PVSurfWriteVTSFS(FastScapeLib*, const char*, PetscScalar*, PetscInt);
PetscErrorCode PVSurfWriteCoordFS(FastScapeLib*, FILE*, PetscScalar*, PetscInt);
PetscErrorCode PVSurfWriteInfFS(FastScapeLib*, FILE*, PetscScalar*, PetscInt);
PetscErrorCode UpdatePVDFileFS(const char*, const char*, const char*,long int*, PetscScalar, PetscInt, PetscInt);
PetscErrorCode SavePvtsFS(FastScapeLib*, PetscScalar, PetscInt, const char*, PetscScalar*);
PetscErrorCode FastScapeSave(FastScapeLib*, PetscInt, PetscScalar);
PetscErrorCode FastScapeFortranCppAdvc(FastScapeLib*, PetscScalar, PetscScalar, PetscInt, PetscScalar, PetscScalar*, PetscScalar*, PetscScalar*, PetscScalar*);
PetscErrorCode FastScapeReadRestart(FastScapeLib*, FILE*);
PetscErrorCode FastScapeWriteRestart(FastScapeLib*, FILE*);
PetscErrorCode FastScapeDestroy(FastScapeLib*);

#ifdef __cplusplus
extern "C"
{
#endif
    /*
    Define FastScape functions as C functions. Must use the exact same function/variable name
    and type as used in FastScape. All function names must be made lowercase, and an
    underscore added at the end. Types must be defined as pointers, and sent to
    FastScape as a reference. Additional functions are available within FastScape,
    see https://fastscape.org/fastscapelib-fortran/ for a list of all functions and
    their input parameters. These functions must be defined at the top here before
    they are used.
    */

    // Function to initialize FastScape.
    void fastscape_init_();

    // Set number of grid points in x (nx) and y (ny)
    void fastscape_set_nx_ny_(const int *nnx, const int *nny);

    // Allocate memory, must be called after set nx/ny.
    void fastscape_setup_();

    // Set the x and y extent of the FastScape model.
    void fastscape_set_xl_yl_(const double *xxl, const double *yyl);

    // Set FastScape timestep. This will vary based on the LaMEM timestep.
    void fastscape_set_dt_(const double *dtt);

    // Initialize FastScape topography.
    void fastscape_init_h_(double *hp);
    
    // Initialize FastScape silt fraction during a restart.
    void fastscape_init_f_(double *Fmixp);

    // Set FastScape erosional parameters on land. These parameters will apply to the stream power law (SPL)
    // and hillslope diffusion for basement and sediment. This can be set between timesteps.
    void fastscape_set_erosional_parameters_(double *kkf,
                                             const double *kkfsed,
                                             const double *mm,
                                             const double *nnn,
                                             double *kkd,
                                             const double *kkdsed,
                                             const double *gg1,
                                             const double *gg2,
                                             const double *pp);

    // Set FastScape marine erosional parameters. This can be set between timesteps.
    void fastscape_set_marine_parameters_(const double *sl,
                                          const double *p1,
                                          const double *p2,
                                          const double *z1,
                                          const double *z2,
                                          const double *r,
                                          const double *l,
                                          const double *kds1,
                                          const double *kds2);

    // Set FastScape boundary conditions.
    void fastscape_set_bc_(const int *jbc);

    // Set FastScape uplift rate. This can be set between timesteps.
    void fastscape_set_u_(double *up);

    // Set advection velocities for FastScape. This can be set between timesteps.
    void fastscape_set_v_(double *ux, double *uy);

    // Set the precipitation rate in meters per year. This can be set between timesteps.
    void fastscape_set_precip_(double *precipp);

    // Run FastScape for a single FastScape timestep.
    void fastscape_execute_step_();

    // Extract from the model the current time step
    void fastscape_get_step_(int *sstep);

    // Set FastScape topography. This can be set between timesteps.
    void fastscape_set_h_(double *hp);

    /*
    Set FastScape basement. This can be set between timesteps. Sediment within FastScape
    is considered as the difference between the topography and basement, though this may differ
    from sediment as seen in LaMEM because the FastScape basement only takes the surface
    velocities into consideration.
    */
    void fastscape_set_basement_(double *b);

    // increment (or uplift) the topography h, the basement height b and the stratigraphic horizons
    void fastscape_set_all_layers_(double *dhp);

    // set the convergence parameters for the Gauss-Seidel iterations performed while 
    // numerically solving the Stream Power law.
    void fastscape_set_tolerance_(const double *tol_relp, 
                                  const double *tol_absp, 
                                  const int *nGSStreamPowerLawMaxp);

    // get the actual number of Gauss-Seidel iterations performed while numerically solving 
    // the Stream Power law during the last time step.
    void fastscape_get_gssiterations_(int *nGSSp); 

    // Copy the current FastScape topography.
    void fastscape_copy_h_(double *hp);

    // Copy the current FastScape silt fraction.
    void fastscape_copy_f_(double *sf);

    // Copy the current FastScape basement.
    void fastscape_copy_basement_(double *b);

    // Copy the current total erosion in meters.
    void fastscape_copy_total_erosion_(double *etotp);
    
    // Reset the total erosion to zero
    void fastscape_reset_cumulative_erosion_();

    // Copy the current drainage area in meters squared.
    void fastscape_copy_drainage_area_(double *ap);

    // Copy the current erosion rate in meters per year.
    void fastscape_copy_erosion_rate_(double *eratep);

    // Copy the current FastScape slopes.
    void fastscape_copy_slope_(double *slopep);

    // Copy the current curvature.
    void fastscape_copy_curvature_(double *curvaturep);

    // Copy the current chi parameter.
    void fastscape_copy_chi_(double *chip);

    // Copy the current catchment area in meter squared.
    void fastscape_copy_catchment_(double *catchp);
    
    // Copy the the geometry and depth of lakes.
    void fastscape_copy_lake_depth_(double *Lp);

    // Extract from the model the model dimensions.
    void fastscape_get_sizes_(int *nnx, int *nny);

    /*
    Extract three fluxes from the model at the current time step: the tectonic flux which is the integral
    over the model of the uplift/subsidence function, the erosion flux which is the integral over the model 
    of the erosion/deposition rate and the boundary flux which is the integral of sedimentary flux across 
    the four boundaries (all in m3/yr)Extract from the model the model dimensions.
    */
    void fastscape_get_fluxes_(double *ttectonic_flux, 
                               double *eerosion_flux, 
                               double *bboundary_flux);

    // Display on the screen basic information about the model
    void fastscape_view_();

    // Display debug information and routine timing
    void fastscape_debug_();

    // Destroy FastScape.
    void fastscape_destroy_();

    // Creates a .vtk file
    void fastscape_vtk_(const double *fp, const double *vexp);

    // Creates a set of .vtk files containing stratigraphic information
    void fastscape_strati_(const int *nstepp, 
                           const int *nreflectorp,
                           const int *nfreqp, 
                           const double *vexp);

#ifdef __cplusplus
}
#endif

#endif