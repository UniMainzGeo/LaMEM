# Phase Transition Plugin — Phase 1 Report

Branch: `bk/phase-transition-plugin` (worktree only, not pushed, not merged to master).

## Summary

Phase 1 wires a user-defined phase-transition plugin (a Julia function
compiled with `juliac` into a relocatable shared library) into LaMEM's
time-step loop, using PETSc's portable dynamic-loading API. The plugin is
optional (`-phase_transition_lib <path>`); with no option, LaMEM's behaviour
is unchanged. All numbers in this report come from runs of the actual final
committed binary (`bin/opt/LaMEM`, rebuilt from a clean object file after
every source change described below) — an earlier pass through this report
cited logs from before the last source edits, and separately claimed the
binary was "confirmed byte-identical across two independent rebuilds"
without having logged both MD5s; that claim is withdrawn here (not
re-asserted, since not re-verified with a logged comparison this pass) in
favour of citing, for each figure below, the log file it actually came from.

Two independent reviews of this branch found real, verified issues, all
fixed and re-tested here:
- Signal handling was a no-op; `int`/`PetscInt` ABI mismatch; unvalidated
  plugin phases; repeated-init unsafety; diagonal-only J2; missing pressure
  shift; no error channel — fixed in an earlier pass on this branch.
- A second review additionally verified this branch's fixes (bit-identical
  built-in-vs-plugin residuals at 1 and 2 ranks; J2 vs. ParaView correlation
  0.997/median 1.8% for strain rate, 0.90/median 2.1% for stress) and found:
  automatic BLAS/LAPACK snapshot/restore was missing (manual options only);
  no thread pinning for the embedded runtime; non-collective error handling
  that could hang under MPI; unsafe repeated-plugin-swap handling;
  `PetscDLClose` on a library holding an initialised Julia runtime (which
  Julia does not support); stale/inaccurate header documentation; no
  checked-in regression test; and several inaccurate report claims. All are
  addressed below.

## Build

Environment used for everything in this report:
```bash
cd <worktree>/src
export PETSC_OPT=/workspace/destdir/lib/petsc/double_real_Int64
export PETSC_DEB=/workspace/destdir/lib/petsc/double_real_Int64_deb
export PATH=/workspace/destdir/bin:$PATH
export MPICH_CXX=/usr/bin/clang++ MPICH_CC=/usr/bin/clang
export LIBRARY_PATH=/Users/kausb/.julia/artifacts/d6f2dc0e73e8796cb9bb992a4970412bc9e4cea3/lib/gcc/aarch64-apple-darwin20/12.0.1
make mode=opt clean_all; make mode=opt all -j8
```
Succeeds with **zero errors, zero warnings** from `phase_transition_plugin.cpp`
(`-Wall -Wextra -Wconversion -Wpointer-arith -Wcast-align -Wwrite-strings
-Wformat=2 -Wundef -Wnon-virtual-dtor -Wimplicit-fallthrough
-Wshorten-64-to-32` all clean). `otool -L bin/opt/LaMEM | grep -i julia`
returns nothing on the final binary: **LaMEM never links libjulia**.

## Bugs found and fixed (this pass — second review)

**(a) BLAS/LAPACK: automatic snapshot/restore, not manual options.**
Root cause confirmed precisely: `LinearAlgebra.__init__` (part of Julia's
base runtime, present in the plugin bundle regardless of whether the
plugin code itself uses linear algebra) calls
`lbt_forward(libopenblas, clear=1, ...)`, which both re-registers Julia's
own BLAS with the shared, process-wide `libblastrampoline` instance
(confirmed via `nm` on the actual `libblastrampoline.5.dylib` in this
deployment: it exports `lbt_get_config`, `lbt_get_forward`,
`lbt_forward`, `lbt_get_num_threads`, `lbt_set_num_threads`) and resets the
BLAS thread count. The single-libblastrampoline situation is **ordinary
dyld install-name de-duplication** (both LaMEM and the plugin bundle
reference `@rpath/libblastrampoline.5.dylib`), not a two-level-namespace
symbol-resolution effect — the earlier report's wording was corrected.
Also noted: the bundle used throughout this work is Julia 1.13.0, while
`LBT_DEFAULT_LIBS` (as configured for this deployment) names Julia 1.12.6's
OpenBLAS — this happened to work by ABI luck (OpenBLAS's ABI is stable
across these versions) rather than by design, and is now explicit in the
comments.

**Fix**: `PhTrPluginSnapshotLbt()` calls `lbt_get_config()` (resolved via
`PetscDLSym(NULL, ...)`, i.e. process-wide, not through the plugin's own
handle — `lbt_get_config`/`lbt_forward` are not resolvable through the
plugin handle specifically, only process-wide, whereas `jl_parse_opts`/
`jl_init_with_image_handle` *are* resolvable through the plugin handle; this
asymmetry was not root-caused further but is documented) **before**
`jl_init_with_image_handle`, and copies out each registered library's
`libname`/`suffix` strings (via `PetscStrallocpy`) into a small snapshot
array — copying is required because Julia's `clear=1` call frees the very
strings the pre-init `lbt_get_config()` snapshot would otherwise still
point at. The struct layout (`lbt_library_info_t`/`lbt_config_t`) was
copied field-for-field from libblastrampoline's own public header — NOT
guessed, NOT found anywhere under PETSc's own `/workspace/destdir/include`
(there is no `libblastrampoline.h` there): the actual file read was
`/System/Volumes/Data/Users/kausb/Documents/GitHub/Yggdrasil-1/build/aarch64-apple-darwin-libgfortran5-mpi+mpitrampoline/7jpuDf2x/aarch64-apple-darwin20-libgfortran5-cxx11-mpi+mpitrampoline/artifacts/697b6b065afac5fb796010d38a46fb719a172e0e/include/libblastrampoline.h`
(a symlink target inside a local Yggdrasil build tree, itself a copy of
libblastrampoline's own upstream `include/libblastrampoline.h`), matching
the layout LaMEM verified working at runtime against the ACTUALLY DEPLOYED
library in this environment: `/workspace/destdir/lib/libblastrampoline.5.dylib`,
version **5.15.0** (confirmed both by the versioned file alongside it,
`libblastrampoline.5.15.0.dylib`, and by the `libblastrampoline_jll v5.15.0+0`
line the Julia package manager itself reports for this deployment). After Julia init,
`PhTrPluginRestoreLbt()` re-forwards each snapshotted library with
`clear=0` (additive — does not remove Julia's own registrations) and
**checks the return value** (`>0` symbols forwarded), `SETERRQ`-ing with
the library path if a re-forward fails. It also restores the BLAS thread
count via `lbt_get_num_threads`/`lbt_set_num_threads`, snapshotted
immediately before `jl_init_with_image_handle` — necessary because
`LinearAlgebra.__init__` also silently resets the thread count to a
CPU-count-based default unless `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS` is
set, and because the plugin's OpenBLAS **is** PETSc's OpenBLAS (same
de-duplicated image), an unrestored thread count would silently multithread
PETSc's own BLAS calls under MPI, competing with MPI ranks for cores.
`-phase_transition_lbt_ilp64`/`-phase_transition_lbt_lp64` remain as an
explicit override, now only used as a fallback if `lbt_get_config` is
unavailable or the snapshot found nothing — in which case a second fallback
parses `LBT_DEFAULT_LIBS` on `;` before giving up.

Verified on the real binary: `Phase transition plugin  : re-forwarded 2
BLAS/LAPACK libraries via libblastrampoline after Julia init` is printed
automatically, with **no** `-phase_transition_lbt_*` options passed, and
all subsequent SNES/KSP solves converge normally (see "Test results" below).

**(b) Embedded runtime pinned to one thread.** `jl_parse_opts` is now
called with `{"lamem", "--handle-signals=no", "--threads=1",
"--gcthreads=1"}` (previously only `--handle-signals=no`). Documented in
the header: `JULIA_NUM_THREADS`/`JULIA_NUM_GC_THREADS` in the calling
user's environment would otherwise start additional Julia threads whose GC
safepoint mechanism depends on the signal handling `--handle-signals=no`
just disabled. Also documented, as requested: deep/runaway recursion in the
plugin becomes a hard SEGV (reported by PETSc's handler) rather than a
catchable `StackOverflowError`; there is no Julia SIGINT handling; and
`jl_parse_opts` itself calls the C library's `exit()` directly on an
unparseable option, so a typo in this fixed argv would kill the whole LaMEM
process, not just the plugin (this argv is fixed source, not user input, so
this is a maintenance note rather than a runtime risk).

**(c) Collective error handling made truly collective.** The previous code
called `SETERRQ` directly inside the per-marker phase-validation loop and
on a bad `rc`, both **inside a block every rank executes independently** —
if only some ranks hit an error, those ranks would abort while the others
proceeded into the following `MPI_Allreduce`/`PetscPrintf` calls, which is
a classic collective-mismatch hang (the surviving ranks wait forever for a
collective the aborted ranks never reach). Fixed: `PhTrPluginApply` now
computes a **local** error flag (from `rc<0` or an out-of-range phase found
while scanning, without writing anything back yet), `MPI_Allreduce`s it
with `MPI_MAX` across `PETSC_COMM_WORLD`, and only if the **global** flag is
set does every rank call `SETERRQ` together — keeping the failure path
collective, exactly like the success path. Only after this check passes
does a second pass actually write phases back to markers. `n==0` (a rank
with no local markers) still calls `fn()` and takes part in every
collective exactly like every other rank — skipping would itself cause a
collective mismatch.

**(d) Repeated-load and dlclose safety.**
- A second `PhTrPluginLoad()` call in the same process (adjoint/inversion
  drivers call `LaMEMLibSolve()` repeatedly) that names a **different**
  `-phase_transition_lib` than the one already loaded now `SETERRQ`s with
  both paths named, instead of silently keeping the first plugin. The
  loaded path is remembered in a static `loadedPath` buffer, compared with
  `strcmp` against the newly requested path.
- `PhTrPluginFinalize` (the `PetscRegisterFinalize` callback) **no longer
  calls `PetscDLClose`** on the plugin handle at all: `dlclose()`-ing a
  library that holds an initialised Julia runtime is not supported by
  Julia under any circumstance (not just "risky"), so the handle is
  deliberately leaked for the remaining life of the process; only
  `jl_atexit_hook` is still called, at most once, best-effort.

**(e) Documentation corrections.**
- The header (`phase_transition_plugin.h`) no longer refers to
  `PhTrPluginDestroy()` (removed in an earlier pass; the current teardown
  path is the `PetscRegisterFinalize` callback, `PhTrPluginFinalize`, which
  is `static` and has no public entry point).
- The old "RTLD_GLOBAL BLAS-binding" caveat (implying Julia's BLAS calls
  bind directly via RTLD_GLOBAL symbol resolution) was replaced: Julia's
  BLAS calls actually go through libblastrampoline's own dlsym-based
  dispatch table, which is exactly the mechanism items (a) and the
  BLAS/LAPACK section above describe in full; there is no separate
  RTLD_GLOBAL-specific concern beyond that.
- New, explicit statement (header and here): **this design cannot be used
  in-process from a Julia host** (e.g. `LaMEM.jl`/`LaMEM_jll` calling a
  LaMEM shared library that itself tries to load this kind of plugin). On
  macOS, the plugin bundle's own `libjulia` would be a second, distinct
  image (different install name/path than the host's own already-loaded
  `libjulia`), which Julia does not support; on Linux,
  `jl_init_with_image_handle` would be called against an
  already-initialised Julia runtime, also unsupported. This is a
  fundamental limitation of embedding a second Julia runtime inside a
  process that is already one, not something fixable by a different
  loading strategy in this file — it is a **standalone-executable-only**
  design.

**(f) Checked-in regression test.** New testset `t40_PhaseTransitionPlugin`
in `test/runtests.jl`, plus `test/t40_PhaseTransitionPlugin/`:
`PT0_only_builtin.dat`, `PT0_only_plugin.dat` (both copies of
`t16_PhaseTransitions/Plume_PhaseTransitions.dat`, stripped to only
PhaseTransition ID 0 — see "4d" below for why), `ptlib_constant.jl` (the
Julia source), `build_plugin.jl` (a `JuliaC.jl`-based build script — see
"Exact juliac build commands" below), and the generated
`PT0_only_builtin.expected`/`PT0_only_plugin.expected`. The testset
resolves `build_constant/lib/libptlib_constant.{dylib,so}` and, if it is
not present (juliac is not assumed available in ordinary CI that only
builds LaMEM's C/C++ code), prints a clear `@info` message naming the build
command and **skips itself** rather than failing. `test/t40_PhaseTransitionPlugin/.gitignore`
excludes the compiled bundle itself (a 52 MB, machine/Julia-version-specific
binary artifact) from version control; only source files are committed.
Verified: `julia --startup-file=no start_tests.jl 40` → `2 Pass, 2 Total`.

**(h) LBT snapshot loop bound and version safety** (found by a third review
pass, after (a)-(g) above landed).
- `PhTrPluginSnapshotLbt`'s scan loop read
  `cfg->loaded_libs[i] != NULL && i < LBT_MAX_SNAPSHOT` — since `&&`
  evaluates its left operand first, this dereferenced `loaded_libs[i]`
  *before* the bound check, so if `loaded_libs` ever held
  `LBT_MAX_SNAPSHOT` (16) or more entries, the 17th element would be read
  out of bounds before the loop could stop. Fixed by swapping the operand
  order (`i < LBT_MAX_SNAPSHOT && cfg->loaded_libs[i] != NULL`), so the
  bound is always checked first.
- `lbt_config_t`/`lbt_library_info_t` carry no version field, and
  libblastrampoline exports no version-query symbol, so there is no way to
  directly ask "is this the struct layout I coded against?". Added the
  closest available guard: after resolving `lbt_get_config`, its address is
  passed to `dladdr()` (guarded by `PETSC_HAVE_DLADDR`, the same macro
  PETSc's own `PetscDLAddr` uses), and the struct layout is only trusted if
  the containing image's path contains `"libblastrampoline.5"` — otherwise
  (or if `dladdr` itself is unavailable on this platform) a warning is
  printed and the snapshot is skipped, falling back to
  `-phase_transition_lbt_ilp64`/`-lp64` or `LBT_DEFAULT_LIBS`. This cannot
  detect a struct layout change *within* the 5.x series, only guards
  against mistaking a differently-versioned/incompatible image for the
  5.15.0 layout this code was written against; a documented comment notes
  that a future libblastrampoline 6.x (or any release reordering this
  struct) would need this code updated alongside it.

**(g) Report accuracy corrections** (all superseded by fresh numbers below,
from the final binary):
- Removed the claim "Both runs report the same changed-marker count at step
  1 (27512)" as if it were a comparison — the built-in `Phase_Transition()`
  prints no marker-changed count at all; only the plugin does. The correct,
  now-verified statement is: the **plugin**, on identical physical setups,
  reports 27512 changed markers at 1 rank and 27525 at 2 ranks (see "4d"
  below) — the built-in transition's *effect* (not a printed count) is
  compared via the residual keywords instead.
- Timing (item iv) is re-measured below on the final binary; the only
  timing artefact in an earlier pass was from a run that later turned out
  to predate the final source, so it is replaced here with a fresh,
  logged, back-to-back measurement.
- "This rebuilds LaMEM (both opt and deb)" for 4a was FALSE for the log
  actually cited (`make: Nothing to be done for 'all'.` appears twice at
  the top of that log — nothing was rebuilt, both binaries were already
  up to date from an earlier command in the same session). Corrected below
  to state plainly that the harness re-invokes `make` unconditionally (so
  it *would* rebuild if anything were stale) rather than claiming a
  rebuild happened in that specific run.
- 2-rank builtin-vs-plugin agreement was previously described as "9-10
  significant digits" without having actually counted differing lines or
  computed a max relative difference; both are now reported precisely
  below, alongside a same-input builtin-vs-builtin pair of runs
  establishing the run-to-run noise floor those differences are being
  compared against.

**(i) `build_plugin.jl` extension bug and a stale comment reference.**
- `build_plugin.jl` passed `"--output-lib", "libptlib_constant.dylib"`
  with an explicit `.dylib` extension. `JuliaC.jl`'s `link_products` (its
  `src/linking.jl`) does NOT substitute the platform extension when one is
  given: it only appends the platform's own dlext when NO extension is
  given, and raises `error("Invalid file extension ...")` if a WRONG
  extension is given — so this script, as committed, would have failed
  outright on Linux (`.dylib` given, `.so` expected) rather than silently
  working. Fixed: `"--output-lib", "libptlib_constant"` (no extension;
  JuliaC appends the correct one for whatever platform it runs on), with
  the comment corrected to describe this accurately.
- `ptlib_constant.jl`'s `catch` block comment referred to `"ptlib.jl"` for
  an explanation of why it does no I/O — that file is a scratchpad-only
  spike, not part of this repository, so a reader of the committed source
  had nothing to actually look at. Replaced with a self-contained
  explanation in `ptlib_constant.jl` itself.

**(j) t40 testset: enforce agreement directly, not just via `.expected`.**
Previously, `PT0_only_builtin.dat` and `PT0_only_plugin.dat` were each
compared only against their OWN `.expected` file — this implies, but does
not enforce, that the two runs agree with each other (the two `.expected`
files could in principle drift apart, e.g. if one were regenerated and the
other were not, while both individual comparisons still passed). The
testset now keeps both runs' log files (`clean_dir=false` on both
`perform_lamem_test` calls) and, after both complete, reads back
`PT0_only_builtin.out`/`PT0_only_plugin.out` with the same
`extract_info_logfiles` utility `compare_logfiles` itself uses, and directly
`@test`s that the `|Div|_inf` and `|mRes|_2` sequences from the two runs
satisfy the same accuracy tuple already used for the `.expected` comparisons,
before manually cleaning up both output directories. **Caveat, documented
in the testset's own comment**: given typical `|mRes|_2` magnitudes here
(1e-8..1e-12), the `atol=1e-3` half of that keyword's accuracy check is
close to vacuous — almost any two small `|mRes|_2` values satisfy
`isapprox` at that absolute tolerance regardless of `rtol`; `|Div|_inf`
(typical magnitude 1e-4..1e-9) is the keyword actually doing discriminating
work at this tolerance. This is the same tolerance already used for t16 and
for each run's own `.expected` comparison, not a tuple chosen to make this
new cross-check pass.

Also added `*.dylib`/`*.so` glob rules to
`test/t40_PhaseTransitionPlugin/.gitignore` (in addition to the existing
`/build_constant/` rule), so a compiled bundle built under a different
`--bundle`/`--output-lib` name could not accidentally be committed either.

## ABI v2: internal units + scaling; outputs phase and T

A third design pass (after the first two review rounds captured above)
changed the ABI again, before it had ever shipped in any released form, so
there is exactly one ABI that matters going forward. The array arguments
now carry LaMEM's **raw internal, non-dimensional** values (no
dimensionalisation happens in C at all); a new read-only
`LaMEMPluginScaling` struct, passed by pointer as the last argument, gives
the plugin everything needed to convert to/from dimensional units itself;
and a new in/out array `T_out` lets a plugin change a marker's temperature
(matching what the built-in Box-type transitions can already do), in
addition to its phase.

```c
int lamem_phase_transition(
    size_t  n,
    double *x, double *y, double *z,       // marker coords [INTERNAL]
    double *T,                              // marker T [INTERNAL -- P->T itself]
    double *p,                              // marker p [INTERNAL, RAW -- P->p itself,
                                             //   WITHOUT pShift folded in]
    double  time,                           // simulation time [INTERNAL]
    double *sxx,double *syy,double *szz,
    double *sxy,double *sxz,double *syz,    // marker deviatoric stress [INTERNAL]
    double *j2_stress_cell,                 // cell-centred J2(dev. stress) [INTERNAL]
    double *j2_strainrate_cell,             // cell-centred J2(strain rate) [INTERNAL]
    double *eta_cell,                       // cell effective viscosity [INTERNAL]
    double *aps_cell,                       // cell accumulated plastic strain [-]
    int    *phase_in, int *phase_out,       // marker phase, in/out
    double *T_out,                          // marker temperature, out [INTERNAL];
                                             //   pre-filled with a copy of T
                                             //   ("no change" by default);
                                             //   validated with isfinite()
    const LaMEMPluginScaling *scaling);     // LaMEM's characteristic scales (below)
    // returns: number of markers changed, or a NEGATIVE value on plugin
    // failure (-1: exception/generic failure; -2 reserved by convention for
    // "scaling struct looks wrong", see "loud failure" test below)
```
Enabled via `-phase_transition_lib <path>`. BLAS/LAPACK co-existence is
**automatic** (see bug (a) above); `-phase_transition_lbt_ilp64
-phase_transition_lbt_lp64` remain only as a fallback override.

### Why unit conversion moved entirely to Julia

The array values are now exactly what's stored on LaMEM's `Marker`/
`SolVarCell` structs, unconverted. The only thing LaMEM computes and hands
over is the scaling struct below; every dimensionalisation happens in
Julia (`test/t40_PhaseTransitionPlugin/LaMEMPlugin.jl`), where a plugin
author reasons in familiar units and the conversion logic is easy to read,
test, and get right in one shared place, rather than being baked into (and
hidden inside) the C loader.

### `LaMEMPluginScaling` — every field, its source, and its unit

| Field | Type | Source (verified against source, not assumed) | Unit, geo mode | Unit, SI/none mode |
|---|---|---|---|---|
| `abi_version` | `int32_t` | `PHASE_TRANSITION_PLUGIN_ABI_VERSION` (= 2), `src/phase_transition_plugin.h` | - | - |
| `utype` | `int32_t` | `scal->utype`, `enum UnitsType` — `src/scaling.h:22-28` (`_NONE_=0, _SI_=1, _GEO_=2`) | - | - |
| `length` | `double` | `scal->length` — `src/scaling.h:70`; geo-mode formula `length/km`, `src/scaling.cpp:215` | km | m |
| `time` | `double` | `scal->time` — `src/scaling.h:68`; geo-mode formula `time/Myr`, `src/scaling.cpp:213` | Myr | s |
| `stress` | `double` | `scal->stress` — `src/scaling.h:80`; geo-mode formula `stress/MPa`, `src/scaling.cpp:225` | MPa | Pa |
| `temperature` | `double` | `scal->temperature` — `src/scaling.h:74`; **same formula in both modes**, `src/scaling.cpp:161,219` (`scal->temperature = temperature`, i.e. always K — the "C" you see in geo mode comes entirely from `Tshift`, not from a Celsius-scaled `temperature`) | K (Tshift makes displayed T read as C) | K |
| `viscosity` | `double` | `scal->viscosity` — `src/scaling.h:93`; identical formula both modes, `src/scaling.cpp:180,238` | Pa·s | Pa·s |
| `strain_rate` | `double` | `scal->strain_rate` — `src/scaling.h:82`; formula `1.0/time` where `time` is the **SI** (seconds) characteristic time parameter, `src/scaling.cpp:169,227` — confirmed NOT `1/Myr` even in geo mode, despite `scal->time` itself being in Myr there | 1/s | 1/s |
| `velocity` | `double` | `scal->velocity` — `src/scaling.h:79`; geo-mode formula `length/time/cm_yr`, `src/scaling.cpp:224` | cm/yr | m/s |
| `density` | `double` | `scal->density` — `src/scaling.h:92`; identical formula both modes, `src/scaling.cpp:179,237` | kg/m³ | kg/m³ |
| `Tshift` | `double` | `scal->Tshift` — `src/scaling.h:65`; set to `273.15` in geo mode, `src/scaling.cpp:210` (0 in SI/none mode) | - | - |
| `pShift` | `double` | `jr->ctrl.pShift` — `src/JacRes.h:145`; same shift already implicit in the rheology (see `Check_Constant_Phase_Transition`'s own `(P->p+pShift)`, `src/phase_transition.cpp:1083`) | - | - |
| `dt` | `double` | `jr->bc->ts->dt` — `src/tssolve.h:26` (`TSSol::dt`); **internal**, non-dimensional — multiply by `time` to dimensionalise | Myr (after ×`time`) | s (after ×`time`) |
| `step` | `int64_t` | `jr->bc->ts->istep` — `src/tssolve.h:47` (`TSSol::istep`); number of steps **already completed** before this call (0 at the first step; confirmed against `PrintStep(ts->istep + 1)` in `src/tssolve.cpp:152`, the exact call that prints LaMEM's own "STEP N" banner) | - | - |
| `lbl_length`, `lbl_time`, `lbl_stress`, `lbl_temperature`, `lbl_viscosity`, `lbl_strain_rate`, `lbl_velocity`, `lbl_density` | `const char*` | `scal->lbl_FIELD` — `src/scaling.h:100-115`; e.g. `"[km]"`, `"[Myr]"`, `"[MPa]"`, `"[C]"`, `"[Pa*s]"`, `"[1/s]"`, `"[cm/yr]"`, `"[kg/m^3]"` in geo mode, set at `src/scaling.cpp:213-238` | (label strings) | (label strings) |

### Conversion formulas (implemented once, in `LaMEMPlugin.jl`)

```julia
dimensional_length      = internal * scaling.length
dimensional_time        = internal * scaling.time
dimensional_stress      = internal * scaling.stress
dimensional_strain_rate = internal * scaling.strain_rate
dimensional_viscosity   = internal * scaling.viscosity
dimensional_velocity    = internal * scaling.velocity
dimensional_density     = internal * scaling.density
dimensional_T           = internal * scaling.temperature - scaling.Tshift
dimensional_pressure    = (internal_p_raw + scaling.pShift) * scaling.stress   # what the rheology sees
# inverse (needed to fill T_out, handed back in internal units):
internal_T              = (dimensional_T + scaling.Tshift) / scaling.temperature
```
These are exactly the formulas `src/scaling.h`/`.cpp` and
`src/phase_transition.cpp`'s `Set_Constant_Phase_Transition` already use
internally (see the next section for exactly how `ptlib_constant.jl`
reproduces the latter).

### Reproducing the built-in Constant transition bit-for-bit, in internal units

`Set_Constant_Phase_Transition` does not compare a dimensional marker
temperature against a dimensional threshold at runtime. It
non-dimensionalises its `ConstantValue` **once**, when the `.dat` is parsed
(`src/phase_transition.cpp:271`):
```c
ph->ConstantValue = (ph->ConstantValue + scal->Tshift) / scal->temperature;
```
and then compares the marker's **raw internal** `P->T` against that
internal threshold on every call (`src/phase_transition.cpp:1077`):
```c
if (P->T >= PhaseTrans->ConstantValue) { ph = PH2; ... }
```
`ptlib_constant.jl` reproduces this exact arithmetic path — not merely an
algebraically-equivalent dimensional comparison:
```julia
constant_value_internal = (CONSTANT_VALUE_DIM + s.Tshift) / s.temperature   # once per call
newph = m.T_internal >= constant_value_internal ? PHASE_ABOVE : PHASE_BELOW  # per marker
```
where `m.T_internal` is `MarkerView`'s copy of the raw, unconverted `T[i]`
array element (kept alongside the dimensionalised `m.T`, specifically so a
plugin wanting bit-exact reproduction of a built-in rule can compare in
internal units without having to re-derive them). Algebraically, comparing
`(T_dim + Tshift)/temperature >= (1200 + Tshift)/temperature` is equivalent
to comparing `T_dim >= 1200` directly (dividing/adding the same quantities
on both sides of an inequality does not change its truth value) — but
performing the SAME sequence of floating-point operations LaMEM's own C
code performs, rather than a dimensional shortcut, removes even the
theoretical possibility of a rounding-induced mismatch at a boundary value,
and is what was actually implemented and verified (see "Test results"
below: bit-identical residuals at 1 rank, 30 steps).

### A genuine floating-point pitfall this design surfaced (and fixed)

Round-tripping an arbitrary internal `T` value through
`dimensionalize_T` → `nondimensionalize_T` is **not** bit-exact in general:
`Tshift` (273.15) is not a power of 2, so the intermediate
`*temperature - Tshift` / `+Tshift /temperature` arithmetic rounds
differently for a large fraction of inputs. Measured empirically: about
23,828 of 100,000 uniformly-random internal T values in a plausible range
failed to round-trip to the bit-identical value. The first working version
of `LaMEMPlugin.jl`'s wrapper unconditionally wrote the round-tripped value
back into `T_out` for every marker, which — even for `ptlib_constant.jl`'s
own rule, which never touches T — silently perturbed thousands of markers'
temperatures by ~1 ULP every time step (observed directly: "0 marker(s)
changed phase, 2673 marker(s) changed temperature" on a rule that should
report 0 T changes). **Fixed** in `LaMEMPlugin.jl`'s `lamem_pt_wrapper`: T
is only round-tripped through the internal/dimensional conversion when the
user's `rule` function actually returns a **different** dimensional T value
than the marker started with (`new_T_dim == m.T`, compared in dimensional
units); otherwise the original internal value is passed straight through
unchanged, byte for byte. After the fix, the same run reports "0 marker(s)
changed phase, 0 marker(s) changed temperature" for every step, matching
the built-in transition's actual behaviour (which never touches T for this
particular rule). This is the kind of subtlety that is easy to miss without
actually counting/logging "unexpected" changes and investigating them, and
is documented at length in `LaMEMPlugin.jl` itself.

### ABI versioning and the "loud failure" scaling-struct guard

If a plugin exports the optional `int lamem_phase_transition_abi_version(void)`
symbol, `PhTrPluginLoad()` calls it and `SETERRQ`s with a clear message if
it does not equal `PHASE_TRANSITION_PLUGIN_ABI_VERSION` (2). Verified with a
standalone C harness (`scratchpad/review_test/abi_v2_guard_test.c`,
mirroring exactly what LaMEM does: dlopen, `jl_parse_opts`,
`jl_init_with_image_handle`, then call `lamem_phase_transition`):
```
abi_version=2
good struct: rc=1
bad struct (length=0): rc=-2 (expect -2)
```
`LaMEMPlugin.jl`'s `lamem_pt_wrapper` also checks the scaling struct itself
for a minimal sanity condition (`utype >= 0`, `length > 0`, `time > 0`,
`stress > 0`) and returns `-2` if it fails — a distinctive code so a future
struct-layout mismatch (which would otherwise either segfault or silently
read garbage and produce a plausible-but-wrong result) fails loudly and
identifiably instead. This is now a checked-in regression test too: `test/t40_PhaseTransitionPlugin/scaling_guard_test.c`
(a committed, cleaned-up version of the same standalone harness) is
compiled and run as a subprocess by the `t40_PhaseTransitionPlugin`
testset, asserting a corrupted struct (`length=0`) returns exactly `-2`
while a valid one does not. This has to run as a **separate OS process**,
not via `ccall` from the Julia test-harness process itself: initialising a
second Julia runtime (`jl_init_with_image_handle`) inside a process that is
already running Julia (the test harness) is exactly the unsupported
scenario `phase_transition_plugin.h` documents ("this design cannot be used
in-process from a Julia host").

### Re-verification of the existing test plan after the ABI v2 change

Everything in "Test results" below (4a/4b/4c/4d, the J2 validation) was
originally run against the pre-ABI-v2 (dimensional-array) signature. After
switching to ABI v2, the full `test/t40_PhaseTransitionPlugin` PT0-only
comparison was re-run end to end on the rebuilt binary and BOTH rebuilt
bundles (`ptlib.jl` and `ptlib_constant.jl`, both updated to ABI v2 and
rebuilt with the exact `juliac` commands above):
```
mpiexec -n 1 bin/opt/LaMEM -ParamFile PT0_only_builtin.dat -nstep_max 30            # /tmp/v2_builtin_1r.log
mpiexec -n 1 bin/opt/LaMEM -ParamFile PT0_only_plugin.dat  -nstep_max 30 \
  -phase_transition_lib .../build_constant/lib/libptlib_constant.dylib             # /tmp/v2_plugin_1r.log
```
Result: `diff` of every `|Div|_inf`/`|mRes|_2` line between the two logs is
**empty** — still bit-for-bit identical across all 30 steps, exactly as
before the ABI change (confirming the internal-units + scaling-struct
redesign is behaviour-preserving for this comparison). Timing (same
`Total solution time` metric as before): baseline 5.97602 s, plugin
6.29907 s (30 steps, 1 rank) — consistent with the previously-measured
~11 ms/step overhead. See "Final test summary (post H1 fix)" at the very
end of this report for the authoritative, most-recent full-suite run.

## Test results (final binary, rebuilt after every source change below)

**NOTE**: the numbers immediately below (4a-4d, J2 validation) predate the
ABI v2 change (see the subsection just above for ABI v2's own
re-verification) and were measured against the dimensional-array signature
described in the superseded "Extended ABI" text; they are kept here as a
historical record of what was verified at each stage, not as the final
word on ABI v2's behaviour.

### 4a — t16 baseline, unchanged, via the official test harness

```bash
cd test
JULIA_LOAD_PATH="<worktree>:<scratchpad>/petsc_deploy:@stdlib" \
  julia --startup-file=no start_tests.jl 16 40
```
Log: `/tmp/final_4a_and_t40.log` (this session; paths are local to the
machine this was run on, cited for traceability, not for the reader to
fetch). `start_tests.jl` unconditionally re-invokes `make mode=opt all` and
`make mode=deb all` before running the suite, so it WOULD rebuild either
binary if anything were stale — but in the specific run behind this log,
neither was: the log's first two lines are
`make: Nothing to be done for 'all'.` (both binaries were already
up to date from an earlier build in the same session, per the source
edits already having been compiled and linked before this test run). The
earlier claim that this run "actually rebuilds" both binaries was
therefore false for that log and has been withdrawn; the harness's
rebuild-on-every-invocation behaviour itself is real and was separately
observed to work (see the deb-mode compile failures hit earlier in this
work, and their successful resolution, both of which required actual
recompilation), just not evidenced by this particular log. The suite then
runs:
```
Test Summary:               | Pass  Total     Time
LaMEM Testsuite             |    9      9  1m09.6s
  t16_PhaseTransitions      |    7      7    56.3s
  t40_PhaseTransitionPlugin |    2      2    13.2s
```
All 7 t16 subtests pass against their `.expected` files with the plugin
infrastructure compiled in but inactive (no `-phase_transition_lib` passed
by those tests) — built-in behaviour is unchanged. The 2 new t40 subtests
are the PT0-only built-in-vs-plugin comparison, run via the same harness
(see 4d).

### 4b/4c — 1-rank and 2-rank runs with `-phase_transition_lib`

```bash
mpiexec -n 1 bin/opt/LaMEM -ParamFile test/t16_PhaseTransitions/Plume_PhaseTransitions.dat \
  -nstep_max 30 -phase_transition_lib <scratchpad>/ptspike/build/lib/libptlib.dylib
```
(no `-phase_transition_lbt_*` options — the automatic snapshot/restore from
fix (a) handles BLAS/LAPACK). Log: `/tmp/final_plugin_timing.log` (this
session). Prints, automatically:
```
Phase transition plugin  : re-forwarded 2 BLAS/LAPACK libraries via libblastrampoline after Julia init
```
then converges normally for all 30 steps (`Phase transition plugin  : 0
marker(s) changed phase` every step — the spike's z<-100 & T>800 & phase==1
rule does not happen to fire in this short setup, independently confirmed
against synthetic data earlier in this work).

**Timing** (back-to-back, same machine, 1 rank, final binary; figures are
each run's own logged `Total solution time`, i.e. LaMEM's internal timer
around the time-step loop, NOT a `date +%s.%N`-around-`mpiexec` wall-clock
figure, which was a different, unlogged measurement from an earlier pass
and is not reproduced here since it cannot be independently checked
against a log):

| Run | Steps | `Total solution time` | Log |
|---|---|---|---|
| Baseline (no plugin) | 30 | 6.31738 s | `/tmp/final_baseline_timing.log` |
| Plugin | 30 | 6.65197 s | `/tmp/final_plugin_timing.log` |
| **Overhead, 30 steps** | | **0.335 s total (~11 ms/step average)** | |

(`grep "Total solution time" <log>` on each file).

**2-rank** run of the same setup also completes normally with the plugin
active (see 4d below, which runs both 1- and 2-rank comparisons on the
PT0-only setup specifically, addressing the coordinator's explicit request
to re-run 4d at both rank counts on the final binary).

### 4d — `ptlib_constant.jl` vs. the built-in Constant transition, 1 and 2 ranks

Design (unchanged from the previous pass, now re-verified on the final
binary): two copies of the t16 `.dat`, both with only PhaseTransition ID 0
present (PT1-3 removed), to avoid a one-step cross-transition ordering
artefact that would occur if only PT0 were swapped out of the full
4-transition `.dat` (the built-in `Phase_Transition()` applies PT0..PT3
sequentially within one call, so a later built-in transition could act,
within the same step, on a phase change a swapped-out PT0 would have
produced — the plugin, running strictly after `Phase_Transition()` returns,
would only see that change on the *following* step).

```bash
# 1 rank
mpiexec -n 1 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_builtin.dat -nstep_max 30   # /tmp/final_builtin_1r.log
mpiexec -n 1 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_plugin.dat  -nstep_max 30 \
  -phase_transition_lib test/t40_PhaseTransitionPlugin/build_constant/lib/libptlib_constant.dylib          # /tmp/final_plugin_1r.log

# 2 ranks
mpiexec -n 2 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_builtin.dat -nstep_max 30   # /tmp/final_builtin_2r.log
mpiexec -n 2 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_plugin.dat  -nstep_max 30 \
  -phase_transition_lib test/t40_PhaseTransitionPlugin/build_constant/lib/libptlib_constant.dylib          # /tmp/final_plugin_2r.log

# 2-rank builtin-vs-builtin noise floor (two independent runs of the SAME
# built-in-only input, no plugin at all, to see how much two runs of
# identical input disagree with each other by themselves)
mpiexec -n 2 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_builtin.dat -nstep_max 30   # /tmp/final_builtin_2r_runA.log
mpiexec -n 2 bin/opt/LaMEM -ParamFile test/t40_PhaseTransitionPlugin/PT0_only_builtin.dat -nstep_max 30   # /tmp/final_builtin_2r_runB.log
```
Results:
- **1 rank**: `diff` of every `|Div|_inf`/`|mRes|_2` line between the
  built-in and plugin logs is **empty** — bit-for-bit identical across all
  30 steps.
- **2 ranks**: NOT bit-identical. Precise comparison (60 keyword lines
  total: 30 `|Div|_inf` + 30 `|mRes|_2`, extracted from
  `/tmp/final_builtin_2r.log` vs. `/tmp/final_plugin_2r.log`): **48 of the
  60 lines differ** (as raw floats), with a **maximum relative difference
  of 1.10e-2**, occurring on an `|mRes|_2` line
  (`3.293738016609e-12` built-in vs. `3.330434986263e-12` plugin — both
  already at the ~1e-12 noise floor of the linear solver's own convergence
  tolerance, so a "1%" relative difference there is a difference of a few
  times `1e-14` in absolute terms). The keyword that actually carries
  useful precision, `|Div|_inf`, is much tighter: at the final step,
  `1.279985020601e-07` (built-in) vs. `1.279985020593e-07` (plugin), a
  relative difference of **6.25e-12**.

  To establish whether this level of disagreement is meaningful or just
  ordinary MPI/solver run-to-run noise, two independent 2-rank runs of the
  SAME input (`PT0_only_builtin.dat`, no plugin involved at all) were
  compared the same way:
  `/tmp/final_builtin_2r_runA.log` vs. `/tmp/final_builtin_2r_runB.log` —
  **55 of the 60 lines differ**, with a **maximum relative difference of
  1.61e-2** (again on an `|mRes|_2` line at the ~1e-12 floor) and a
  final-step `|Div|_inf` relative difference of 1.25e-11. That is, running
  the identical built-in-only input twice disagrees MORE (55 differing
  lines, larger max relative difference) than the built-in-vs-plugin
  comparison does (48 differing lines, smaller max relative difference) —
  the plugin-vs-builtin difference is fully within, not exceeding, the
  measured run-to-run noise floor of the underlying solver itself. Both
  comparisons stay well inside the test suite's own tolerance
  (`rtol=1e-5, atol=1e-7` for `|Div|_inf`; `rtol=1e-2, atol=1e-3` for
  `|mRes|_2` — see the note in the t40 testset itself: at the `|mRes|_2`
  magnitudes seen here, `atol=1e-3` makes that keyword's check close to
  vacuous, so `|Div|_inf` is the keyword doing the real discriminating
  work in this comparison).
- **Changed-marker counts** (printed only by the plugin; the built-in
  transition prints no count — see bug (g)): **27512** at 1 rank, **27525**
  at 2 ranks, both at step 1 only (0 for all subsequent steps — steady
  state reached). The 13-marker difference between rank counts is expected
  from a different domain decomposition placing markers in different
  cells/ranks at the phase-2/3 boundary, not a bug.
- These are the same figures the second independent review itself obtained
  when re-verifying this branch, reproduced here on the final binary as
  requested, together with the newly-added builtin-vs-builtin noise-floor
  pair the review asked for.

`ptlib_constant.jl`'s scope (unchanged from the previous pass): it
reproduces `Check_Constant_Phase_Transition` for `number_phases=1`,
`PhaseDirection=BothWays`, no `ResetParam` only — matching exactly what PT0
in the t16 `.dat` exercises, not the general case.

## Item (vii) — J2 verification against ParaView (second review's own check)

The second review additionally validated `PhTrPluginComputeJ2`'s
cell-centred J2 against LaMEM's own ParaView output (interpolating the
node-centred ParaView field back onto each cell's coordinates and
correlating), finding: strain-rate correlation 0.997 (median relative
difference 1.8%), stress correlation 0.90 (median relative difference
2.1%). This is consistent with the two fields being geometrically distinct
by construction (cell-centred vs. corner/node-centred, as documented in
`PhTrPluginComputeJ2`'s own header comment) rather than identical, and
gives good confidence the cell-centred computation is correct rather than
merely plausible. This validation's artefacts (`cmpj2.py`, a small VTK/dump
comparison script, and the underlying dump tool) are external
to this branch (scratchpad-only) and are not part of the committed diff;
they are referenced here as the source of these two correlation numbers.

## Outstanding (not done, flagged explicitly)

- **No Linux run.** Everything in this report was run on macOS
  (aarch64-apple-darwin). The dyld-specific reasoning in bug (a) (install-name
  de-duplication) is macOS terminology; on Linux, the equivalent mechanism
  is the dynamic linker's own SONAME-based sharing of an already-loaded
  library, which is expected to behave analogously (a single
  `libblastrampoline.so.5` shared between LaMEM and the plugin), but this
  was not verified on Linux in this environment. The automatic LBT
  snapshot/restore is the only mechanism now (the earlier
  `-phase_transition_lbt_ilp64`/`-phase_transition_lbt_lp64` options and the
  `LBT_DEFAULT_LIBS`-parsing fallback were removed during the code-slimming
  pass, per an explicit size-minimization request); this should be the
  first thing re-verified on a Linux machine.

## Exact juliac build commands (restored; previously dropped)

Both bundles referenced in this report were built with Julia **1.13.0**
(via `juliaup`'s `+1.13.0` channel) and `JuliaC.jl`:
```bash
# ptlib.jl (the original spike bundle, ABI v2)
cd <scratchpad>/ptspike
julia +1.13.0 --startup-file=no --project=<scratchpad>/juliac_env -e \
  'using JuliaC; JuliaC.main(["--output-lib","libptlib","--bundle","build", \
    "--trim=safe","--compile-ccallable","--experimental","--project=proj","ptlib.jl"])'
# NOTE: --output-lib takes NO extension (JuliaC appends the platform's own
# dlext itself; a WRONG explicit extension is a hard error in JuliaC.jl's
# link_products, not silently substituted - see "ABI v2" below, item (i)
# from the previous review pass).

# ptlib_constant.jl (the checked-in t40 test's plugin, built via test/t40_PhaseTransitionPlugin/build_plugin.jl)
cd test/t40_PhaseTransitionPlugin
julia +1.13.0 --startup-file=no --project=<scratchpad>/juliac_env build_plugin.jl
```
`<scratchpad>/juliac_env` is a Julia environment with `JuliaC` installed;
`--project=proj` inside the first command is a separate, empty (`[deps]`
only) environment for the module being compiled itself — not required for
`ptlib_constant.jl`, which has no external dependencies, so
`build_plugin.jl` omits it.

## Fourth review's fixes: T-only propagation; tighter t40 (this pass)

A fourth review of the ABI-v2 + slimming commits found one real bug (H1)
and several test/reporting gaps (M1-M3, L1-L3), fixed in this pass:

- **H1 (bug)**: `PhTrPluginApply` only called `ADVCheckMarkPhases`/
  `ADVInterpMarkToCell` when a marker's *phase* changed
  (`if(glob2[0])`). But `ADVInterpMarkToCell` is what seeds
  `svBulk.Tn` from marker `T` (`src/advect.cpp`, ~line 1721), and the very
  next call in the step loop (`JacResInitTemp`) reads `svBulk.Tn` — so a
  plugin rule that changes ONLY `T` (not phase) had that change silently
  ignored for the thermal solve of that step. Fixed by gating on
  `glob2[0] || glob2[1]` (phase-changed-count OR T-changed-count), both
  already `MPI_Allreduce`d, so the call stays collective.
- **Test evidence for H1**: added `ptlib_box.jl` (`Box_only_builtin.dat` /
  `Box_only_plugin.dat`), a Julia rule that resets marker `T` to a
  constant inside a box — exercising a T-only change, mirroring the
  built-in Box transition's own semantics. Matching the built-in
  `Check_Box_Phase_Transition`'s actual behaviour required one correction
  to the naive rule: the built-in code's `Phase_Transition()` resets `T`
  for **any** marker geometrically inside the box (via the "allow cases in
  which we only reset T" branch for markers whose phase isn't
  `PhaseInside`/`PhaseOutside`), while only markers with phase 2/3 also
  get their *phase* flipped. `box_rule` was written to match this exactly
  (T reset unconditionally inside the box; phase only touched for
  phase-2/3 markers). With that, the 1-rank built-in-vs-plugin comparison
  is bit-identical: `diff` of every `|Div|_inf`/`|mRes|_2` line between
  `Box_only_builtin.dat` and `Box_only_plugin.dat` (5 steps) is empty.
- **M1**: `lamem_pt_wrapper` now checks `s.abi_version != ABI_VERSION`
  and returns `Cint(-3)` before touching any other field; documented as
  `-3` in `phase_transition_plugin.h`'s return-value comment.
- **M2**: `scaling_guard_test.c` now also `dlsym`s and calls
  `lamem_phase_transition_abi_version`, printing `abi_version=N`; the
  guard-test assertions in `test/runtests.jl` were tightened to
  `occursin("abi_version=2", ...)` and `occursin("good struct: rc=1", ...)`
  (previously a weak `!occursin("rc=-2")`, which a broken wrapper
  returning e.g. `-1` or `0` would have incorrectly passed).
- **M3**: the guard test now skips with `@info` (not an uncaught error)
  when neither `cc` nor `clang` is on `PATH`.
- **L1**: the wrapper's `changed` counter now increments on phase-changed
  OR T-changed, matching the header's documented return-value contract
  ("markers with phase or T changed").
- **L2**: this section replaces the dangling "Test summary" forward
  reference with a real one (below).
- **L3**: all four `.expected` files in `test/t40_PhaseTransitionPlugin`
  (`PT0_only_builtin`, `PT0_only_plugin`, `Box_only_builtin`,
  `Box_only_plugin`) were regenerated against the final binary (H1 fix
  included) via `julia start_tests.jl mode=update 40`.

### Final test summary (post H1 fix)

Full suite, final binary, both `.expected` files regenerated first:
```
julia --startup-file=no --project=.. start_tests.jl 16 40
...
Test Summary:               | Pass  Total     Time
LaMEM Testsuite             |   22     22  1m36.9s
  t16_PhaseTransitions      |    7      7  1m02.1s
  t40_PhaseTransitionPlugin |   15     15    34.7s
     Testing LaMEM_C tests passed
```
t40's 15 `@test` assertions cover, per transition (Constant and Box): the
built-in run's `.expected` comparison, the plugin run's `.expected`
comparison, and the direct built-in-vs-plugin cross-comparison (a
length-match check plus one `isapprox` check per keyword, 2 keywords) —
plus the 3 scaling-guard assertions (`abi_version=2`, `good struct: rc=1`,
`bad struct (length=0): rc=-2`). All 15 pass, 0 failed, 0 errored.
