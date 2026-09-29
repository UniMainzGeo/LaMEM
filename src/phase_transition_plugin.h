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
// libblastrampoline's lbt_forward() to register its OWN bundled OpenBLAS as
// the default BLAS/LAPACK target. Because LaMEM and the plugin share the
// SAME process-wide libblastrampoline instance (only one copy is ever
// dlopen'd, since it is referenced by the same install name), this
// OVERWRITES the forwarding table entries PETSc's own startup had set up
// via the LBT_DEFAULT_LIBS environment variable -- even though
// LBT_DEFAULT_LIBS is still set, Julia's init does not consult it and
// unconditionally registers its own default, breaking every subsequent
// PETSc BLAS/LAPACK call ("no BLAS/LAPACK library loaded for idamax_()"
// etc.) unless this is corrected. PhTrPluginLoad() re-registers PETSc's own
// BLAS libraries (additively, via lbt_forward(..., clear=0, ...), which
// does not remove Julia's registrations, only re-points the specific
// symbols PETSc needs) immediately after Julia init, using:
//
//   -phase_transition_lbt_ilp64 <path-to-ILP64-openblas>
//   -phase_transition_lbt_lp64  <path-to-LP64-openblas>
//
// i.e. the SAME two libraries named in LBT_DEFAULT_LIBS (order: ILP64
// first, then LP64 -- matches PETSc_jll's own convention for PETSc built
// with 64-bit BLAS indices). If these options are not given, a plugin is
// still loaded, but PETSc's own linear solves will likely fail on the very
// next SNES/KSP solve with a libblastrampoline "no BLAS/LAPACK library
// loaded" error; a warning is printed in that case. See
// doc/phase_transition_plugin_PHASE1_REPORT.md for how this was diagnosed
// and verified fixed.
//
// IMPORTANT CAVEATS (Phase 1 scope; see doc/phase_transition_plugin_PHASE1_REPORT.md):
//  - The plugin (and thus the embedded Julia runtime) is loaded and
//    initialised AT MOST ONCE PER PROCESS, on the first call to
//    PhTrPluginLoad() that finds -phase_transition_lib set. LaMEMLibSolve()
//    can run multiple times in one process (adjoint/inversion drivers call
//    it repeatedly -- see src/adjoint.cpp). Julia's runtime cannot be
//    re-initialised after jl_atexit_hook(), and dlclose()-ing a library that
//    holds an initialised Julia runtime is not supported, so subsequent
//    PhTrPluginLoad() calls in the same process are no-ops that reuse the
//    already-loaded plugin, and PhTrPluginDestroy() must be called AT MOST
//    ONCE per process, from a PetscRegisterFinalize() callback -- not from
//    LaMEMLibSolve() -- so it fires exactly once at PetscFinalize(), after
//    the last LaMEMLibSolve() call.
//  - The Julia side must stay single-threaded (no @spawn/Threads.@threads in
//    the plugin) and Julia's own stack-overflow guard page is disabled by
//    --handle-signals=no (see PhTrPluginLoad() for why), so a runaway
//    recursive plugin function will segfault instead of raising a catchable
//    StackOverflowError.
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
//  - On Linux, RTLD_GLOBAL (used to load the plugin, see PhTrPluginLoad())
//    means a BLAS/LAPACK symbol referenced by the Julia side could bind to
//    PETSc's already-loaded BLAS instead of the plugin's own -- harmless for
//    a plugin that does no linear algebra (as in Phase 1), but relevant if a
//    Phase 2 plugin uses LinearAlgebra.
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
