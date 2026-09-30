/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//..............   USER-DEFINED DYLIB PLUGINS (Julia)   .....................
//---------------------------------------------------------------------------
// Loads a juliac-built shared library (dylib_plugin = <path> in the .dat
// file, or -dylib_plugin <path> on the command line, which overrides the
// .dat value). This is the plugin ABI: the LaMEMPluginScaling struct it
// shares with every hook below, plus the hooks themselves. Current hooks:
//   lamem_phase_transition() - called once per time step, per MPI rank.
//     ABI v2: all arrays are LaMEM's internal (non-dimensional) units; the
//     plugin dimensionalises using the LaMEMPluginScaling struct passed in:
//       dim        = internal*scale        (length/time/stress/viscosity/
//                                            strain_rate/velocity/density)
//       T_dim      = T*temperature - Tshift
//       p_rheology = (p_raw + pShift)*stress
//     Outputs: phase_out, T_out (pre-filled with phase_in/T, internal units).
//     Returns markers with phase or T changed, or <0 on failure (-1 generic,
//     -2 bad scaling struct, -3 ABI mismatch).
// Optional symbol lamem_plugin_abi_version(): if present, must return
// DYLIB_PLUGIN_ABI_VERSION or loading fails.
#ifndef dylib_plugins_h_
#define dylib_plugins_h_
//---------------------------------------------------------------------------

#define DYLIB_PLUGIN_ABI_VERSION 2

// fixed-width types (not PetscScalar/PetscInt): the ABI must not depend on
// this PETSc build's configuration; field order/types must match the
// plugin's Julia-side struct exactly
struct LaMEMPluginScaling
{
	int32_t abi_version, utype;
	double  length, time, stress, temperature, viscosity, strain_rate, velocity, density;
	double  Tshift, pShift, dt;
	int64_t step;
};

typedef int (*DylibPluginFn)(
	size_t n,
	double *x,  double *y,  double *z,
	double *T,  double *p,
	double  time,
	double *sxx, double *syy, double *szz,
	double *sxy, double *sxz, double *syz,
	double *j2_stress_cell, double *j2_strainrate_cell,
	double *eta_cell,       double *aps_cell,
	int    *phase_in, int *phase_out,
	double *T_out,
	const LaMEMPluginScaling *scaling);

struct AdvCtx;
struct FB;
struct DBMat;

PetscErrorCode DylibPluginLoad(AdvCtx *actx, FB *fb);
PetscErrorCode DylibPluginPhaseTransition(AdvCtx *actx);
PetscBool      DylibPluginIsActive(void);
PetscErrorCode DylibPluginCheckPhaseTr(DBMat *dbm);

//---------------------------------------------------------------------------
#endif
