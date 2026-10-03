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
// .dat value). This header is the plugin ABI (v3): the plain C structs below,
// mirrored field-by-field in scripts/dylib_plugins/LaMEMPlugin.jl, plus the
// symbols a plugin exports:
//
//   int32_t lamem_plugin_abi_version(void)               - mandatory, must
//       return DYLIB_PLUGIN_ABI_VERSION
//   int32_t lamem_plugin_struct_sizes(int64_t *s, int32_t n) - mandatory,
//       writes sizeof() of the n <= DYLIB_PLUGIN_NUM_STRUCTS structs, in the
//       order Markers, Cells, Step, Scaling, and returns
//       DYLIB_PLUGIN_NUM_STRUCTS; the loader compares them with its own
//       sizeof()s, so a layout mismatch fails at load time
//   int32_t lamem_phase_transition(const LaMEMPluginMarkers*,
//       const LaMEMPluginCells*, const LaMEMPluginStep*,
//       const LaMEMPluginScaling*) - called once per time step, per MPI rank
//       (also with n == 0), after the built-in phase transitions.
//       Returns the number of markers it changed, or <0 on failure
//       (-1 generic, -2 bad scaling struct, -3 ABI mismatch).
//
// All values are LaMEM's internal (non-dimensional) units; the plugin
// dimensionalises them using LaMEMPluginScaling:
//   dim        = internal*scale  (scale: length, time, stress, ... below)
//   T_dim      = T*temperature - Tshift
//   p_rheology = (p_raw + pShift)*stress     (marker p and cell pn)
//
// Markers: one entry per local marker. *_in arrays are read-only copies of
// the marker fields; *_out arrays are pre-filled with the same values and
// are what LaMEM writes back. Writable are phase, T, p, APS, ATS, the
// deviatoric stress S and the displacement U; only the position is
// read-only (moving markers is advection's job). The marker pressure is a
// history variable (incremented by the grid pressure change in advection,
// projected to the cells as svBulk.pn = p_old): changing it imposes a
// pressure-history jump that acts through the volumetric elastic term
// IKdt*(p - pn) of the next solve, i.e. only for compressible phases.
// LaMEM copies an *_out value back only if it differs from *_in, so a
// plugin must leave unchanged entries bit-for-bit as pre-filled.
//
// Cells: structure-of-arrays over the local cells of this rank (ncells), a
// read-only copy of SolVarCell taken once per step, plus the cell-centred J2
// invariants and the phase ratios (phRat[c*numPhases + ph]). A marker's cell
// is cell_index[i] (0-based).
//
// The ABI part of this header is plain C: define LAMEM_PLUGIN_ABI_ONLY
// before including it to get only the structs, without LaMEM's own
// declarations (used by test/t40_PhaseTransitionPlugin/scaling_guard_test.c).
#ifndef dylib_plugins_h_
#define dylib_plugins_h_
//---------------------------------------------------------------------------
#include <stddef.h>
#include <stdint.h>

#define DYLIB_PLUGIN_ABI_VERSION 3
#define DYLIB_PLUGIN_NUM_STRUCTS 4

// fixed-width types only (not PetscScalar/PetscInt): the ABI must not depend
// on this PETSc build's configuration. Field order/types must match the
// plugin's Julia-side structs exactly (checked via lamem_plugin_struct_sizes)
typedef struct LaMEMPluginMarkers
{
	size_t         n;          // number of local markers
	const int32_t *cell_index; // local cell of each marker, 0-based

	// read-only
	const double  *x, *y, *z;  // position

	// writable fields: input values (p: raw, without pShift)
	const int32_t *phase_in;
	const double  *T_in, *p_in, *aps_in, *ats_in;
	const double  *sxx_in, *syy_in, *szz_in, *sxy_in, *sxz_in, *syz_in;
	const double  *ux_in, *uy_in, *uz_in;

	// writable fields: output values (pre-filled with the input values)
	int32_t       *phase_out;
	double        *T_out, *p_out, *aps_out, *ats_out;
	double        *sxx_out, *syy_out, *szz_out, *sxy_out, *sxz_out, *syz_out;
	double        *ux_out, *uy_out, *uz_out;
} LaMEMPluginMarkers;

typedef struct LaMEMPluginCells
{
	size_t         ncells;     // number of local cells
	int32_t        numPhases;  // phases in the material table
	int32_t        reserved;   // explicit padding, always 0

	// SolVarCell::svDev
	const double  *eta, *eta_st, *I2Gdt, *Hr, *aps, *psr;

	// SolVarCell::svBulk
	const double  *theta, *rho, *IKdt, *alpha, *Tn, *pn, *rho_pf, *mf, *phi, *Ha, *cond;

	// remaining SolVarCell fields
	const double  *sxx, *syy, *szz;  // deviatoric stress
	const double  *hxx, *hyy, *hzz;  // history stress (elastic)
	const double  *dxx, *dyy, *dzz;  // total deviatoric strain rate
	const int32_t *free_surf;        // SolVarCell::FreeSurf
	const double  *ux, *uy, *uz;     // total displacement
	const double  *ats;              // accumulated total strain
	const double  *eta_cr;           // creep viscosity
	const double  *DIIdif, *DIIdis, *DIIprl, *DIIfk, *DIIpl; // relative strain rates
	const double  *yield;            // average yield stress

	// derived, cell-centred second invariants of the deviatoric
	// stress and strain rate (see DylibPluginComputeJ2)
	const double  *j2_stress, *j2_strainrate;

	// phase ratios, ncells*numPhases, cell-major: phRat[c*numPhases + ph]
	const double  *phRat;
} LaMEMPluginCells;

typedef struct LaMEMPluginStep
{
	double  time, dt;
	int64_t step;
} LaMEMPluginStep;

typedef struct LaMEMPluginScaling
{
	int32_t abi_version, utype;
	double  length, time, stress, temperature, viscosity, strain_rate, velocity, density;
	double  conductivity, expansivity, dissipation_rate;
	double  Tshift, pShift;
} LaMEMPluginScaling;

typedef int32_t (*DylibPluginFn)(
    const LaMEMPluginMarkers *markers,
    const LaMEMPluginCells   *cells,
    const LaMEMPluginStep    *step,
    const LaMEMPluginScaling *scaling);

#ifndef LAMEM_PLUGIN_ABI_ONLY
struct AdvCtx;
struct FB;
struct DBMat;

PetscErrorCode DylibPluginLoad(AdvCtx *actx, FB *fb);
PetscErrorCode DylibPluginPhaseTransition(AdvCtx *actx);
PetscBool      DylibPluginIsActive(void);
PetscErrorCode DylibPluginCheckPhaseTr(DBMat *dbm);
#endif

//---------------------------------------------------------------------------
#endif
