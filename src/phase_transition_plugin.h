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
// The plugin library must export a symbol with EXACTLY this C ABI:
//
//   int lamem_phase_transition(
//       size_t        n,                    // number of local markers
//       double       *x,                    // marker x [dimensional]
//       double       *y,                    // marker y [dimensional]
//       double       *z,                    // marker z [dimensional]
//       double       *T,                    // marker temperature [scal units: Celsius
//                                            //   in "geo" mode, Kelvin in "SI"/"none" mode
//                                            //   -- i.e. P->T*scal->temperature - scal->Tshift,
//                                            //   the SAME quantity/units LaMEM's own marker
//                                            //   I/O and ParaView temperature output use]
//       double       *p,                    // marker pressure, INCLUDING the same pressure
//                                            // shift the rheology/plasticity and the built-in
//                                            // "Pressure" Constant transition see:
//                                            // (P->p + jr->ctrl.pShift)*scal->stress
//                                            // [dimensional, MPa in "geo" mode]
//       double        time,                 // current simulation time [dimensional]
//       double       *sxx,                  // marker deviatoric stress xx [dimensional]
//       double       *syy,                  // marker deviatoric stress yy [dimensional]
//       double       *szz,                  // marker deviatoric stress zz [dimensional]
//       double       *sxy,                  // marker deviatoric stress xy [dimensional]
//       double       *sxz,                  // marker deviatoric stress xz [dimensional]
//       double       *syz,                  // marker deviatoric stress yz [dimensional]
//       double       *j2_stress_cell,       // host-cell J2 (deviatoric stress) [dimensional],
//       double       *j2_strainrate_cell,   // host-cell J2 (strain rate) [dimensional]
//                                            // Both are computed the same way LaMEM's own
//                                            // ParaView j2_dev_stress / j2_strain_rate fields
//                                            // are (diagonal cell components + the 4
//                                            // surrounding edge components per off-diagonal
//                                            // direction, averaged, see PhTrPluginComputeJ2()
//                                            // in phase_transition_plugin.cpp) -- NOT a
//                                            // diagonal-only approximation.
//       double       *eta_cell,             // host-cell effective viscosity [dimensional,
//                                            // linear Pa*s -- NOT log10, unlike the ParaView
//                                            // visc_total field which IS log10-scaled]
//       double       *aps_cell,             // host-cell accumulated plastic strain [-]
//       int          *phase_in,             // marker phase, current
//       int          *phase_out);           // marker phase, requested by plugin
//                                            // returns: number of markers changed,
//                                            // or a NEGATIVE value to signal a plugin-side
//                                            // failure (LaMEM aborts the run with SETERRQ
//                                            // if a negative value is returned). The Julia
//                                            // side must wrap its body in try/catch and
//                                            // return -1 on any caught exception: an
//                                            // uncaught exception inside a @ccallable
//                                            // function aborts the whole process outside
//                                            // PETSc's error handling.
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
