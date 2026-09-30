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
//..............   USER-DEFINED PHASE TRANSITION PLUGIN (Julia)   ...........
//---------------------------------------------------------------------------
//
// Phase 1 proof-of-concept: load a shared library compiled from Julia with
// `juliac` (a Julia "bundle") and call a hardcoded C-callable function once
// per time step, once per MPI rank, to let the plugin decide phase changes
// for the local markers.
//
// The plugin library must export a symbol with EXACTLY this C ABI
// (this is ABI v2; see "ABI VERSIONING" below):
//
//   int lamem_phase_transition(
//       size_t        n,                    // number of local markers
//       double       *x,                    // marker x [INTERNAL, non-dimensional]
//       double       *y,                    // marker y [INTERNAL, non-dimensional]
//       double       *z,                    // marker z [INTERNAL, non-dimensional]
//       double       *T,                    // marker temperature [INTERNAL, non-dimensional
//                                            //   -- i.e. P->T itself, unconverted]
//       double       *p,                    // marker pressure [INTERNAL, non-dimensional,
//                                            //   RAW solution pressure -- i.e. P->p itself,
//                                            //   WITHOUT jr->ctrl.pShift folded in; the
//                                            //   pressure the rheology/plasticity actually
//                                            //   sees is (P->p + pShift), see `scaling` below]
//       double        time,                 // current simulation time [INTERNAL,
//                                            //   non-dimensional -- jr->bc->ts->time itself]
//       double       *sxx,                  // marker deviatoric stress xx [INTERNAL]
//       double       *syy,                  // marker deviatoric stress yy [INTERNAL]
//       double       *szz,                  // marker deviatoric stress zz [INTERNAL]
//       double       *sxy,                  // marker deviatoric stress xy [INTERNAL]
//       double       *sxz,                  // marker deviatoric stress xz [INTERNAL]
//       double       *syz,                  // marker deviatoric stress yz [INTERNAL]
//       double       *j2_stress_cell,       // host-cell J2 (deviatoric stress) [INTERNAL],
//       double       *j2_strainrate_cell,   // host-cell J2 (strain rate) [INTERNAL]
//                                            // Both use the same cell/edge-averaging geometry
//                                            // as LaMEM's own ParaView j2_dev_stress /
//                                            // j2_strain_rate fields (diagonal cell components
//                                            // + the 4 surrounding edge components per
//                                            // off-diagonal direction, averaged, see
//                                            // PhTrPluginComputeJ2() in
//                                            // phase_transition_plugin.cpp) -- NOT a
//                                            // diagonal-only approximation -- but are left
//                                            // in LaMEM's internal (non-dimensional) units
//                                            // here; multiply by scaling->stress /
//                                            // scaling->strain_rate to dimensionalise.
//       double       *eta_cell,             // host-cell effective viscosity [INTERNAL,
//                                            //   svCell->svDev.eta itself]
//       double       *aps_cell,             // host-cell accumulated plastic strain [-]
//                                            //   (already dimensionless in LaMEM itself)
//       int          *phase_in,             // marker phase, current
//       int          *phase_out,            // marker phase, requested by plugin
//       double       *T_out,                // marker temperature, requested by plugin
//                                            //   [INTERNAL, non-dimensional, same units as
//                                            //   the T array above]. Pre-filled by LaMEM with
//                                            //   a copy of the T array (i.e. "no change" is
//                                            //   the default); the built-in Box-type
//                                            //   transitions can also reset a marker's
//                                            //   temperature (e.g. a linear/halfspace T
//                                            //   profile inside the box), so this gives the
//                                            //   plugin the same capability. Every value MUST
//                                            //   be finite (checked with isfinite() before any
//                                            //   write-back; NaN/Inf is treated as a plugin
//                                            //   failure, collectively across all MPI ranks --
//                                            //   see PhTrPluginApply()).
//       const LaMEMPluginScaling *scaling);  // LaMEM's characteristic scales (see below);
//                                            // read-only for the duration of the call, valid
//                                            // only until lamem_phase_transition returns (it
//                                            // points at a function-local stack struct in
//                                            // PhTrPluginApply, NOT retained across calls)
//                                            // returns: number of markers changed (phase
//                                            // and/or T; see PhTrPluginApply for exactly how
//                                            // this is counted and reported), or a NEGATIVE
//                                            // value to signal a plugin-side failure (LaMEM
//                                            // aborts the run with SETERRQ, collectively
//                                            // across all MPI ranks, if a negative value is
//                                            // returned on ANY rank). The Julia side must wrap
//                                            // its body in try/catch and return -1 on any
//                                            // caught exception: an uncaught exception inside
//                                            // a @ccallable function aborts the whole process
//                                            // outside PETSc's error handling.
//
// UNIT HANDLING IS DELIBERATELY ALL ON THE JULIA SIDE. The C side (LaMEM
// itself) does NOT dimensionalise or non-dimensionalise anything in the
// arrays above: every array (in and out) carries LaMEM's raw internal,
// non-dimensional values, exactly as stored on the Marker/SolVarCell
// structs. The ONLY thing LaMEM adds is the read-only `scaling` struct, so
// that unit conversion -- which is where a Phase 1 plugin author most needs
// to reason in familiar units (km, Myr, MPa, C, ...) -- happens entirely in
// Julia, where it is easy to read, test and get exactly right, rather than
// being baked into (and hidden inside) the C plugin loader. The conversion
// formulas a plugin needs (and that
// test/t40_PhaseTransitionPlugin/LaMEMPlugin.jl implements once, for reuse
// by any plugin that includes it) are:
//
//   dimensional_length      = internal_length      * scaling->length
//   dimensional_time        = internal_time         * scaling->time
//   dimensional_stress      = internal_stress       * scaling->stress
//   dimensional_strain_rate = internal_strain_rate  * scaling->strain_rate
//   dimensional_viscosity   = internal_viscosity    * scaling->viscosity
//   dimensional_velocity    = internal_velocity     * scaling->velocity
//   dimensional_density     = internal_density      * scaling->density
//   dimensional_T           = internal_T * scaling->temperature - scaling->Tshift
//   pressure the rheology/plasticity sees, dimensional:
//                             = (internal_p + scaling->pShift) * scaling->stress
//
// and the exact inverses (needed to fill T_out, which must be handed back
// in internal units):
//
//   internal_T = (dimensional_T + scaling->Tshift) / scaling->temperature
//   internal_length      = dimensional_length      / scaling->length
//   internal_time        = dimensional_time        / scaling->time
//   internal_stress      = dimensional_stress       / scaling->stress
//   internal_strain_rate = dimensional_strain_rate  / scaling->strain_rate
//
// These are exactly the formulas src/scaling.h/.cpp and
// src/phase_transition.cpp's Set_Constant_Phase_Transition (which
// non-dimensionalises its own ConstantValue threshold as
// `(ConstantValue + Tshift)/temperature` before comparing it against a
// marker's raw internal P->T) already use internally -- see "ABI v2:
// internal units + scaling" in doc/phase_transition_plugin_PHASE1_REPORT.md
// for where each formula was verified against LaMEM's own source, and for
// the exact Julia code that reproduces the Constant transition's
// comparison bit-for-bit using them.
//
// LaMEMPluginScaling: a plain, packed, read-only struct, isbits on the
// Julia side (test/t40_PhaseTransitionPlugin/LaMEMPlugin.jl defines the
// matching Julia struct). Every "scal->FIELD" / "jr->FIELD" comment below
// names the exact source struct field this value is copied from (verified
// against src/scaling.h, src/JacRes.h and src/tssolve.h, not assumed):
//
//   typedef struct
//   {
//       int32_t  abi_version;   // = 2 for this signature (see "ABI VERSIONING")
//       int32_t  utype;         // scal->utype (enum UnitsType, src/scaling.h):
//                                //   0 = _NONE_ (non-dimensional in/out)
//                                //   1 = _SI_   (SI units in/out)
//                                //   2 = _GEO_  (geological units in/out: the
//                                //       units named below for each field)
//       double   length;        // scal->length      [km  in geo, m   in SI]
//       double   time;          // scal->time        [Myr in geo, s   in SI]
//       double   stress;        // scal->stress      [MPa in geo, Pa  in SI]
//       double   temperature;   // scal->temperature [K in both geo AND si --
//                                //   scal->temperature is never rescaled by
//                                //   utype, only Tshift differs; the "C" you
//                                //   see in geo-mode dimensional T comes
//                                //   entirely from Tshift below, not from a
//                                //   Celsius-scaled `temperature` factor]
//       double   viscosity;     // scal->viscosity   [Pa*s in both geo and SI]
//       double   strain_rate;   // scal->strain_rate [1/s in both geo and SI
//                                //   -- confirmed from scaling.cpp: computed
//                                //   as 1/time_SI, i.e. NOT 1/Myr even in geo
//                                //   mode, despite scal->time itself being
//                                //   in Myr there]
//       double   velocity;      // scal->velocity    [cm/yr in geo, m/s in SI]
//       double   density;       // scal->density     [kg/m^3 in both geo and SI]
//       double   Tshift;        // scal->Tshift: dimensional T = internal*temperature - Tshift
//                                //   (Tshift = 273.15 in geo mode -> Celsius output;
//                                //    Tshift = 0 in SI/none mode -> Kelvin output)
//       double   pShift;        // jr->ctrl.pShift: dimensional pressure seen by the
//                                //   rheology/plasticity = (internal_p + pShift)*stress
//                                //   (the `p` array above is the RAW internal_p, i.e.
//                                //   WITHOUT pShift folded in -- add it yourself)
//       double   dt;            // jr->bc->ts->dt [INTERNAL, non-dimensional; multiply
//                                //   by `time` above to dimensionalise]
//       int64_t  step;          // jr->bc->ts->istep: number of time steps ALREADY
//                                //   COMPLETED before this call (0 at the very first
//                                //   step; the step about to be solved is step+1,
//                                //   matching the "STEP N" banner LaMEM itself prints
//                                //   via PrintStep(ts->istep + 1) in tssolve.cpp)
//       const char *lbl_length, *lbl_time, *lbl_stress, *lbl_temperature,
//                  *lbl_viscosity, *lbl_strain_rate, *lbl_velocity, *lbl_density;
//                                // scal->lbl_FIELD: the exact unit-label strings LaMEM's
//                                //   own screen output uses for each field above (e.g.
//                                //   "[km]", "[Myr]", "[MPa]", "[C]", "[Pa*s]", "[1/s]",
//                                //   "[cm/yr]", "[kg/m^3]" in geo mode) -- NUL-terminated,
//                                //   valid only for the duration of the call (they point
//                                //   into jr->scal itself, not a copy)
//   } LaMEMPluginScaling;
//
// ABI VERSIONING: this struct, the trailing `scaling` argument, and the new
// `T_out` argument were all introduced together as ABI v2; there was no
// shipped v1 (the original Phase 1 signature was changed before ever being
// released, so there is exactly one C ABI in this repository's history that
// matters going forward: this one). A plugin built against a different
// argument list will crash (wrong argument count/stack layout) if loaded.
// To fail loudly instead of crashing on a future ABI break: if the plugin
// ALSO exports an optional symbol `int lamem_phase_transition_abi_version(void)`
// (a @ccallable Julia function is enough, no struct needed), PhTrPluginLoad()
// calls it and SETERRQs with a clear message if the returned value does not
// equal 2. A plugin that does not export this optional symbol at all is NOT
// detected this way -- exporting the version symbol is therefore effectively
// mandatory in practice for a plugin that wants a safe failure mode instead
// of a hard crash on a future ABI break; this is documented, not silently
// worked around. test/t40_PhaseTransitionPlugin/LaMEMPlugin.jl exports it
// for any plugin that `include`s it.
//
// enabled by the runtime option:
//
//   -phase_transition_lib <path-to-shared-library>
//
// If the option is absent, this module does nothing and LaMEM's built-in
// phase transitions run exactly as before.
//
// IMPORTANT (BLAS/LAPACK): a juliac-built plugin bundle carries its own
// copy of Julia's base runtime, which includes LinearAlgebra/OpenBLAS_jll
// even if the plugin code itself never uses linear algebra. During
// jl_init_with_image_handle, Julia's LinearAlgebra.__init__ calls
// libblastrampoline's lbt_forward(..., clear=1, ...) to register its OWN
// bundled OpenBLAS as the default BLAS/LAPACK target (and, incidentally,
// resets the BLAS thread count to its own default). Because LaMEM and the
// plugin share the SAME process-wide libblastrampoline instance -- this is
// ordinary dyld install-name de-duplication: both LaMEM and the plugin
// bundle reference the library as "@rpath/libblastrampoline.5.dylib", so
// dyld resolves the plugin's copy to the identical already-loaded image;
// it has nothing to do with two-level-namespace symbol *resolution*, just
// the fact that both processes' load commands name the same install name
// -- this OVERWRITES the forwarding table entries PETSc's own startup had
// set up via the LBT_DEFAULT_LIBS environment variable, and silently
// changes PETSc's own BLAS thread count, breaking every subsequent PETSc
// BLAS/LAPACK call ("no BLAS/LAPACK library loaded for idamax_()" etc.)
// and/or silently multithreading PETSc's BLAS under MPI, unless corrected.
//
// PhTrPluginLoad() fixes both automatically: it snapshots libblastrampoline's
// forwarding table (which libraries are registered, and their per-library
// suffix) and the current BLAS thread count BEFORE jl_init_with_image_handle
// runs (the snapshot strings must be copied out, not just pointer-saved,
// because Julia's clear=1 call frees them), then re-forwards each
// snapshotted library (additively, clear=0, which does not remove Julia's
// own registrations, only re-points the specific symbols PETSc needs) and
// restores the thread count immediately after Julia init returns. This
// requires no options and no environment variable beyond whatever PETSc
// itself already needed (e.g. LBT_DEFAULT_LIBS for PETSc_jll>=3.25 -- see
// the project's own petsc-jll-325-lbt-default-libs memory note). Two
// options remain as an explicit override / fallback, only used if
// libblastrampoline's lbt_get_config() is unavailable or the snapshot
// found nothing (in which case the fallback additionally tries parsing
// the LBT_DEFAULT_LIBS environment variable on ';' before giving up):
//
//   -phase_transition_lbt_ilp64 <path-to-ILP64-openblas>
//   -phase_transition_lbt_lp64  <path-to-LP64-openblas>
//
// See doc/phase_transition_plugin_PHASE1_REPORT.md for how this was
// diagnosed and verified fixed.
//
// IMPORTANT CAVEATS (Phase 1 scope; see doc/phase_transition_plugin_PHASE1_REPORT.md):
//  - The plugin (and thus the embedded Julia runtime) is loaded and
//    initialised AT MOST ONCE PER PROCESS, on the first call to
//    PhTrPluginLoad() that finds -phase_transition_lib set. LaMEMLibSolve()
//    can run multiple times in one process (adjoint/inversion drivers call
//    it repeatedly -- see src/adjoint.cpp). Julia's runtime cannot be
//    re-initialised after jl_atexit_hook(), and dlclose()-ing a library that
//    holds an initialised Julia runtime is not supported by Julia at all --
//    the plugin handle is deliberately NEVER closed (leaked for the life of
//    the process); only jl_atexit_hook is called, at most once, from a
//    PetscRegisterFinalize() callback (not from LaMEMLibSolve() itself), so
//    it fires exactly once at PetscFinalize(), after the last
//    LaMEMLibSolve() call. A second PhTrPluginLoad() call in the same
//    process that names a DIFFERENT -phase_transition_lib than the one
//    already loaded is treated as an error (SETERRQ), not silently ignored.
//  - This design is for the STANDALONE LaMEM executable only. It CANNOT be
//    used from a Julia host process (e.g. LaMEM.jl / LaMEM_jll calling into
//    a LaMEM shared library that in turn tries to dlopen this kind of
//    plugin): on macOS, the plugin bundle's own libjulia would be a SECOND,
//    distinct libjulia image (a different install name/path than the
//    host's own already-loaded libjulia), which Julia does not support; on
//    Linux, jl_init_with_image_handle would be called on an
//    already-initialised Julia runtime, which is also unsupported. This is
//    a fundamental limitation of embedding a second Julia runtime inside a
//    process that is itself already a Julia runtime, not something a
//    different loading strategy in this file could work around.
//  - The Julia side must stay single-threaded: jl_parse_opts() is called
//    with --threads=1 --gcthreads=1 in addition to --handle-signals=no, so
//    that JULIA_NUM_THREADS/JULIA_NUM_GC_THREADS in the calling user's
//    environment cannot start additional Julia threads whose GC safepoint
//    mechanism depends on the very signal handling --handle-signals=no just
//    disabled. Combined effects of --handle-signals=no to keep in mind:
//    deep/runaway recursion in the plugin becomes a hard SEGV (reported by
//    PETSc's own handler) rather than a catchable Julia StackOverflowError;
//    Julia's own Ctrl-C/SIGINT handling is not installed either; and
//    jl_parse_opts() itself calls the C library's exit() directly on an
//    unparseable option, which would terminate the whole LaMEM process, not
//    just the plugin, so the argv passed to it must stay exactly as tested.
//  - Cell-level quantities (j2_stress_cell, j2_strainrate_cell, eta_cell,
//    aps_cell) reflect the PREVIOUS converged nonlinear-solver state (they
//    are read at the top of the time step, before the current step's
//    SNESSolve), i.e. they are lagged by one step, exactly like the
//    built-in Phase_Transition()'s use of svCell.
//  - The plugin must be idempotent within a single time step:
//    ADVSelectTimeStep() can force the step to be redone (cutting dt) via
//    `continue`, in which case Phase_Transition() and this plugin both run
//    again for the same nominal step on markers that may already carry the
//    plugin's previous verdict.
//  - Julia's BLAS calls go through libblastrampoline's own dlsym-based
//    dispatch (see the BLAS/LAPACK note above), not through RTLD_GLOBAL
//    symbol binding directly; there is no separate RTLD_GLOBAL-specific
//    BLAS-binding concern beyond the libblastrampoline forwarding-table
//    interaction already described above.
//---------------------------------------------------------------------------
#ifndef phase_transition_plugin_h_
#define phase_transition_plugin_h_
//---------------------------------------------------------------------------

// ABI v2 scaling struct passed by pointer as the last argument to
// lamem_phase_transition (see the header comment above for the exact
// meaning/units of every field). Plain, packed-by-natural-alignment (no
// #pragma pack: field order was chosen doubles-before-int64-before-pointers
// to avoid any padding surprises, but this is not itself part of the ABI
// contract - what matters is that the Julia-side struct in
// test/t40_PhaseTransitionPlugin/LaMEMPlugin.jl mirrors this EXACT field
// order and EXACT field types, which is verified by the "loud failure on a
// wrong struct" test in the t40 testset).
// NOTE: every numeric field here is a fixed-width/fixed-layout C type
// (int32_t/double/int64_t), deliberately NOT PetscScalar/PetscInt: the ABI
// must be identical regardless of how THIS PARTICULAR PETSc build defines
// those (e.g. PetscInt is 64-bit in the Int64 configuration this was built
// and tested against, but is not guaranteed 64-bit in general, and
// PetscScalar could in principle be complex or single-precision in another
// build) - the plugin's Julia-side struct is fixed at Cdouble/Cint/Clong,
// so the C side must match that exactly, not whatever this build's PETSc
// happens to use.
struct LaMEMPluginScaling
{
	int32_t     abi_version;
	int32_t     utype;
	double      length;
	double      time;
	double      stress;
	double      temperature;
	double      viscosity;
	double      strain_rate;
	double      velocity;
	double      density;
	double      Tshift;
	double      pShift;
	double      dt;
	int64_t     step;
	const char *lbl_length;
	const char *lbl_time;
	const char *lbl_stress;
	const char *lbl_temperature;
	const char *lbl_viscosity;
	const char *lbl_strain_rate;
	const char *lbl_velocity;
	const char *lbl_density;
};

#define PHASE_TRANSITION_PLUGIN_ABI_VERSION 2

struct AdvCtx;

// load the plugin library (no-op if -phase_transition_lib is not given, or
// if a plugin has already been loaded earlier in this process)
PetscErrorCode PhTrPluginLoad(AdvCtx *actx);

// call the plugin once for all local markers (no-op if no plugin is loaded)
PetscErrorCode PhTrPluginApply(AdvCtx *actx);

// true if a plugin library has been successfully loaded
PetscBool PhTrPluginIsActive(void);

//---------------------------------------------------------------------------
#endif
