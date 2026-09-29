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
//       double       *T,                    // marker temperature [dimensional]
//       double       *p,                    // marker pressure [dimensional]
//       double        time,                 // current simulation time [dimensional]
//       double       *sxx,                  // marker deviatoric stress xx [dimensional]
//       double       *syy,                  // marker deviatoric stress yy [dimensional]
//       double       *szz,                  // marker deviatoric stress zz [dimensional]
//       double       *sxy,                  // marker deviatoric stress xy [dimensional]
//       double       *sxz,                  // marker deviatoric stress xz [dimensional]
//       double       *syz,                  // marker deviatoric stress yz [dimensional]
//       double       *j2_stress_cell,       // host-cell J2 (deviatoric stress) [dimensional]
//       double       *j2_strainrate_cell,   // host-cell J2 (strain rate) [dimensional]
//       double       *eta_cell,             // host-cell effective viscosity [dimensional]
//       double       *aps_cell,             // host-cell accumulated plastic strain [-]
//       int          *phase_in,             // marker phase, current
//       int          *phase_out);           // marker phase, requested by plugin
//                                            // returns: number of markers changed
//
// enabled by the runtime option:
//
//   -phase_transition_lib <path-to-shared-library>
//
// If the option is absent, this module does nothing and LaMEM's built-in
// phase transitions run exactly as before.
//---------------------------------------------------------------------------
#ifndef phase_transition_plugin_h_
#define phase_transition_plugin_h_
//---------------------------------------------------------------------------

struct AdvCtx;

// load the plugin library (no-op if -phase_transition_lib is not given)
PetscErrorCode PhTrPluginLoad(AdvCtx *actx);

// call the plugin once for all local markers (no-op if no plugin is loaded)
PetscErrorCode PhTrPluginApply(AdvCtx *actx);

// release plugin resources (no-op if no plugin is loaded)
PetscErrorCode PhTrPluginDestroy(void);

// true if a plugin library has been successfully loaded
PetscBool PhTrPluginIsActive(void);

//---------------------------------------------------------------------------
#endif
