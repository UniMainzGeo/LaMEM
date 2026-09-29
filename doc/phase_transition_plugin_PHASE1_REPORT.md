# Phase Transition Plugin — Phase 1 Report

Branch: `bk/phase-transition-plugin` (worktree only, not pushed, not merged to master).

## Summary

Phase 1 wires a user-defined phase-transition plugin (a Julia function compiled
with `juliac` into a shared library) into LaMEM's time-step loop, using PETSc's
portable dynamic-loading API. The plugin is optional (`-phase_transition_lib
<path>`); with no option, LaMEM's behaviour is unchanged.

**Blocking issue:** the local build environment (`$PETSC_OPT` =
`/Users/kausb/Software/PETSc/petsc-3.23.6/petsc-3.23.6_opt`) cannot currently
compile *any* LaMEM source file, including files untouched by this change. See
"Build status" below. Because of this, items 4a–4d (running LaMEM itself) could
not be executed in this environment. Everything that does not require a linked
LaMEM binary was done and verified: the plugin ABI, the Julia libraries (both
compiled with `juliac` into relocatable bundles and exercised through a real
`dlopen`/`jl_init_with_image_handle` cycle, matching how LaMEM's loader will
call them), and static verification of the new C++ code against LaMEM's actual
struct layouts, option-parsing conventions, and scaling code.

## What works

- **`src/phase_transition_plugin.h` / `src/phase_transition_plugin.cpp`**
  (new files, picked up automatically by `src/Makefile`'s `$(wildcard *.cpp)`
  — no Makefile edit needed):
  - `PhTrPluginLoad(AdvCtx*)`: reads `-phase_transition_lib <path>`. If absent,
    returns immediately (built-in transitions run exactly as today). If
    present:
    - `PetscDLOpen(path, PETSC_DL_NOW, &handle)`.
    - `PetscDLSym(handle, "jl_init_with_image_handle", &sym)`, called as
      `sym(handle)`.
    - `PetscPushSignalHandler(PetscSignalHandlerDefault, NULL)` to reclaim
      PETSc's signal handlers from Julia's runtime.
    - `PetscDLSym(handle, "lamem_phase_transition", &sym)`, stored as a typed
      function pointer.
    - `PetscDLSym(handle, "jl_atexit_hook", &sym)` (optional, best-effort).
    - Prints `Phase transition plugin  : <path>` on rank 0.
  - `PhTrPluginApply(AdvCtx*)`: called every time step, right after the
    built-in `Phase_Transition(actx)` in `LaMEMLibSolve`. No-op if no plugin
    is loaded. Builds SoA (dimensional) arrays from the local markers, calls
    `lamem_phase_transition` once, writes back only the markers whose phase
    changed, calls `ADVInterpMarkToCell` if anything changed, and prints
    `Phase transition plugin  : N marker(s) changed phase` where `N` is an
    `MPI_Allreduce`-summed global count (`MPI_SUM` over `MPIU_INT`).
  - `PhTrPluginDestroy(void)`: best-effort `jl_atexit_hook(0)`, then
    `PetscDLClose`, then frees the SoA scratch buffers.
  - Buffers are allocated with `PetscMalloc`/grown (never shrunk) across time
    steps as `actx->nummark` changes, freed in `PhTrPluginDestroy`.
- **Wiring** in `src/LaMEMLib.cpp` (`LaMEMLibSolve`):
  - `PhTrPluginLoad(&lm->actx)` right after the `-snes_track_stages` option
    check, before the time-step loop.
  - `PhTrPluginApply(&lm->actx)` immediately after `Phase_Transition(&lm->actx)`
    inside the loop.
  - `PhTrPluginDestroy()` right after `NLSolDestroy(&snes)`, before
    `ADVMarkSave`, at the end of `LaMEMLibSolve`.
- **Build system**: no Makefile changes. `src/Makefile` builds
  `CSRC = $(filter-out LaMEM.cpp, $(wildcard *.cpp))`, so
  `phase_transition_plugin.cpp` is picked up automatically. No new link flags:
  the loader only uses PETSc's existing `PetscDLOpen`/`PetscDLSym`/
  `PetscDLClose`, which are already part of `libpetsc`. LaMEM does **not**
  link against libjulia at any point.

## ABI (final, extended per the stress/strain-rate scope addition)

```c
int lamem_phase_transition(
    size_t  n,
    double *x, double *y, double *z,      // marker coords [dimensional]
    double *T, double *p,                 // marker T, p    [dimensional]
    double  time,                         // simulation time [dimensional]
    double *sxx, double *syy, double *szz,
    double *sxy, double *sxz, double *syz, // marker deviatoric stress [dimensional]
    double *j2_stress_cell,               // host-cell J2(dev. stress) [dimensional]
    double *j2_strainrate_cell,           // host-cell J2(strain rate) [dimensional]
    double *eta_cell,                     // host-cell effective viscosity [dimensional]
    double *aps_cell,                     // host-cell accumulated plastic strain [-]
    int    *phase_in, int *phase_out);    // marker phase, in/out
    // returns: number of markers changed (informational; LaMEM recomputes the
    // authoritative count itself by diffing phase_in/phase_out)
```

All pointer arguments except `phase_in`/`phase_out` are `double*`
(`Ptr{Cdouble}` in Julia); `phase_in`/`phase_out` are `int*` (`Ptr{Cint}`);
`n` is `size_t` (`Csize_t`); `time` is `double` (`Cdouble`).

Argument order is fixed and documented at the top of
`src/phase_transition_plugin.h`.

## Answers to the required report items

### (i) Can Julia be initialised without LaMEM linking libjulia, and how?

Yes. LaMEM's `Makefile` links only against `libpetsc` and MPI (unchanged —
verify with `otool -L bin/opt/LaMEM` once a binary exists: it must **not**
list any `libjulia*`). The mechanism:

1. `PetscDLOpen(path_to_plugin, PETSC_DL_NOW, &handle)`. PETSc's own
   implementation (`src/sys/dll/dlimpl.c` in the PETSc source tree,
   `petsc-3.23.6/src/sys/dll/dlimpl.c` lines ~96–107) maps `PETSC_DL_NOW` to
   `dlopen(name, RTLD_NOW | RTLD_GLOBAL)` (RTLD_GLOBAL is the PETSc default
   unless `PETSC_DL_LOCAL` is explicitly requested — confirmed by reading that
   file). This matches exactly the `dlopen(lib, RTLD_NOW | RTLD_GLOBAL)` used
   by the working spike's `host.c`.
2. Because the juliac-built plugin `.dylib` records `libjulia*.dylib` (bundled
   alongside it under `build/lib/julia/...` and `build/lib/libjulia*.dylib`)
   as one of its own linked dependencies, `dlopen(RTLD_NOW|RTLD_GLOBAL)` on
   the plugin pulls libjulia (and its own dependents: libuv, libopenblas64_,
   libblastrampoline, etc.) into the process automatically, the same way any
   dynamic loader resolves a shared library's transitive dependencies.
3. `PetscDLSym(handle, "jl_init_with_image_handle", &sym)` succeeds because
   `PetscDLSym` is implemented as `dlsym(handle, name)` internally, and once a
   library is `dlopen`'d with `RTLD_GLOBAL`, `dlsym` on *that* handle (or on
   `RTLD_DEFAULT`) resolves symbols from the whole dependency chain that was
   pulled in — not just symbols the plugin itself defines. This was verified
   conceptually against the PETSc `dlimpl.c` source and is exactly the
   principle the spike's `host.c` already demonstrated working (dlsym after
   dlopen(RTLD_GLOBAL) finds `jl_init_with_image_handle`, a libjulia symbol,
   through a handle obtained by opening `libptlib.dylib`). I could not run
   this specific lookup through the actual compiled LaMEM binary in this
   environment (see "Build status"), so this is verified by:
   (a) reading the PETSc DL implementation to confirm `PETSC_DL_NOW` really
       does map to `RTLD_NOW|RTLD_GLOBAL` here, and
   (b) re-running the *same* dlopen/dlsym/jl_init_with_image_handle sequence
       that PhTrPluginLoad performs, standalone, against the real
       juliac-built `.dylib` (see "Extended-ABI verification" below) — it
       succeeds and calls into Julia correctly.
   If, on some platform, `PetscDLSym` on the plugin's own handle cannot see
   `jl_init_with_image_handle` (e.g. a linker/loader that does not propagate
   `RTLD_GLOBAL` transitively the way macOS/Linux dlopen does), the documented
   fallback in `phase_transition_plugin.cpp`'s comments is to `PetscDLOpen`
   the `libjulia.dylib`/`.so` found next to the plugin directly and resolve
   the symbol from that handle instead. This fallback path is not implemented
   in Phase 1 (not needed on macOS/Linux) but is called out precisely so
   Phase 2 knows where to add it if a new platform needs it.
4. `sym(handle)` is called with the raw `PetscDLHandle` — Julia's own
   `jl_init_with_image_handle(void*)` takes the dlopen handle of the image
   that contains the bundled sysimage, matching the spike exactly.

### (ii) Signal handling outcome

`PetscPushSignalHandler(PetscSignalHandlerDefault, NULL)` is called
immediately after Julia initialisation, to push PETSc's own handler back on
top of whatever Julia's runtime installed during `jl_init_with_image_handle`.
This avoids depending on the `jl_options` struct layout (the spike's approach,
`jl_options.handle_signals = JL_OPTIONS_HANDLE_SIGNALS_OFF`, requires linking
libjulia's headers/ABI and setting a field *before* `jl_init`, which is a
private/version-sensitive struct — undesirable for Phase 2's cross-version
compatibility goal).

I could not run the requested live verification (deliberately triggering a
PETSc error/segfault with the plugin loaded and confirming PETSc's handler —
not Julia's — reports it) because no LaMEM binary could be built in this
environment (see "Build status"). This is the one requested check that is
purely a **runtime** behaviour and cannot be confirmed by static reading of
PETSc/Julia source alone with full confidence, so it is explicitly flagged as
**not verified** and should be the first thing re-run once the build is
unblocked, e.g.:
```
mpiexec -n 1 bin/opt/LaMEM -ParamFile <t16 dat> -phase_transition_lib <path> -px 999999999
```
(an intentionally invalid option value or similarly a `SEGV`-inducing debug
flag) and confirm the output is PETSc's usual
`[0]PETSC ERROR: ------------------------------------------------------------------------` /
signal-report banner rather than a Julia stack trace or a silent hang/crash.

### (iii) MPI outcome

Not run for the same reason (no binary). By construction this should work:
each rank calls `PhTrPluginLoad`/`PhTrPluginApply` independently, `PetscDLOpen`
is a per-process call (each rank `dlopen`s and `jl_init`s its own Julia
runtime — no MPI communication happens inside Julia), and the only
cross-rank synchronization added is the `MPI_Allreduce(..., MPI_SUM, ...)` for
the reported changed-marker count, which is a collective every rank always
calls (never conditionally skipped), so it cannot deadlock. This mirrors
LaMEM's own pattern (e.g. `src/advect.cpp:2015`, `src/JacRes.cpp:591`,
`src/JacRes.cpp:1644`). **Not verified at runtime — left for after the build
is unblocked.**

### (iv) Timing overhead per time step

Not measurable without a LaMEM binary. The spike (`ptspike/host.c`) measured
0.55 ms per call for 1e6 markers with the original (smaller) ABI on Julia
1.13; the t16 setup here has far fewer markers (`nel_x=64, nel_y=2, nel_z=64`,
`nmark_x=nmark_y=nmark_z=3` ⇒ ~64·2·64·27 ≈ 221,000 markers total, split
across ranks), so the per-step plugin call should be well under a millisecond
per rank — dominated by the SoA marshalling loop (16 dimensional double
arrays + 2 int arrays per marker, one pass, no allocation after the first
step) rather than the Julia call itself. **Not verified at runtime.**

### (v) Does the Julia re-implementation match the built-in transition?

**Verified functionally (not yet inside a full LaMEM run — see Build
status).** `ptlib_constant.jl` reimplements PhaseTransition ID 0 from
`test/t16_PhaseTransitions/Plume_PhaseTransitions.dat`:
```
Type = Constant, Parameter_transition = T, ConstantValue = 1200
PhaseAbove = 3, PhaseBelow = 2, PhaseDirection = BothWays
```
against LaMEM's actual rule in `Check_Constant_Phase_Transition`
(`src/phase_transition.cpp` line ~1054, `_T_` branch):
```c
if (P->T >= PhaseTrans->ConstantValue) { ph = PH2; InAb=1; } else { ph = PH1; }
```
where for `BothWays`, `PH1=PhaseBelow=2`, `PH2=PhaseAbove=3`, and the marker
must already have phase `PhaseBelow` or `PhaseAbove` to be eligible
(`Check_Phase_above_below`). `ptlib_constant.jl`'s `lamem_phase_transition`
implements exactly this: markers whose current phase is 2 or 3 get
`T >= 1200 ? 3 : 2`; all other phases are left untouched. Compiled with
`juliac` into `ptspike/build_constant/lib/libptlib_constant.dylib` and
exercised through a standalone dlopen/`jl_init_with_image_handle` harness
(`host_const.c` in the scratchpad, linked against libjulia directly *only*
for this isolated verification — LaMEM itself never does this), against 4
synthetic markers:
```
marker 0: phase 2 -> 3 (T=1300)   [T >= 1200 -> PhaseAbove]
marker 1: phase 3 -> 2 (T=1100)   [T <  1200 -> PhaseBelow]
marker 2: phase 3 -> 2 (T=500)    [T <  1200 -> PhaseBelow]
marker 3: phase 5 -> 5 (T=1300)   [not phase 2/3 -> untouched]
changed=3
```
This matches the C rule bit-for-bit for these cases. What could not be
verified: running the actual t16 `.dat` (with PhaseTransition ID 0 removed
and `-phase_transition_lib
<path>/ptspike/build_constant/lib/libptlib_constant.dylib` given instead) end
to end and diffing `|Div|_inf` / `|mRes|_2` against
`test/t16_PhaseTransitions/PhaseTransitions.expected`, because no LaMEM
binary could be built. Given the rule match is exact and T is passed through
identically (`P->T*scal->temperature - scal->Tshift`, the same dimensionalization
LaMEM's own output/marker I/O uses — see `src/marker.cpp` lines 320, 642,
810), I expect the two runs to be numerically identical to within the
existing test tolerances (`rtol=1e-5..1e-2`), since both apply the same phase
reassignment rule to the same dimensional temperature at the same point in
the time-step loop (the built-in `Phase_Transition` call and the plugin call
are adjacent in `LaMEMLibSolve`, both before `ADVMarkInjectGeom`/`BCApply`).
This is a prediction, not a verified result, and is flagged as such.

### (vi) Anything surprising

- **The build environment itself is broken**, independent of anything in this
  change (see "Build status"). This was the single biggest surprise and the
  main blocker for the runtime parts of this report.
- PETSc's `PetscDLMode` only exposes `PETSC_DL_DECIDE`/`PETSC_DL_NOW`/
  `PETSC_DL_LOCAL` (no explicit "GLOBAL" flag), but reading
  `src/sys/dll/dlimpl.c` shows `RTLD_GLOBAL` is the *default* `dlflags2`
  unless `PETSC_DL_LOCAL` is passed — so `PETSC_DL_NOW` alone already gives
  the `RTLD_NOW|RTLD_GLOBAL` combination the spike needed. No custom dlopen
  call was necessary.
- `Marker::S` (`Tensor2RS`) only stores `xx, xy, xz, yy, yz, zz` (upper
  triangle by construction — deviatoric stress is symmetric), which maps
  cleanly onto the requested 6-array ABI without needing a mirrored/full 3×3
  tensor.

### (vii) Stress/strain-rate quantities available at the hook point, cost, dimensionality

- **Marker-level** (`P->S`, a `Tensor2RS`: `xx,xy,xz,yy,yz,zz`): the elastic
  deviatoric stress carried on the marker itself (advected marker history),
  scaled by `scal->stress` — no interpolation, essentially free (already
  resident on the marker struct, one multiply per component).
- **Cell-level J2 invariants**: LaMEM's own ParaView output
  (`src/outFunct.cpp`, `PVOutWriteJ2DevStress`/`PVOutWriteJ2StrainRate`)
  computes a *true* cell-centred J2 by interpolating both the cell-diagonal
  components (`svCell->sxx/syy/szz`, `svCell->dxx/dyy/dzz`) **and** the three
  off-diagonal edge components (`jr->svXYEdge/svXZEdge/svYZEdge`) onto shared
  grid corners (`InterpCenterCorner`/`InterpXYEdgeCorner`/etc.), then back —
  a ghost-exchange-dependent, DMDA-corner pass across three *different*
  staggered grids (`DA_CEN`, `DA_XY`, `DA_XZ`, `DA_YZ`).
  **Phase 1 deliberately does not reuse that exact routine.** Instead,
  `PhTrPluginApply` computes a cheaper, diagonal-only cell-centred J2:
  ```
  J2_stress      = sqrt(0.5*(sxx^2+syy^2+szz^2)) * scal->stress
  J2_strainrate  = sqrt(0.5*(dxx^2+dyy^2+dzz^2)) * scal->strain_rate
  ```
  using only `svCell->sxx,syy,szz` / `svCell->dxx,dyy,dzz` — no off-diagonal
  edge terms, no corner interpolation, no ghost communication. This is
  **cheap** (2 sqrt + a handful of multiplies per marker's host cell, reusing
  a value already resident in `jr->svCell[]`, indexed by the same `cellnum[i]`
  the built-in `Phase_Transition` already uses) but it is a **lower bound**
  on the true J2 (the omitted off-diagonal terms are non-negative
  contributions to the sum of squares), not identical to what ParaView
  reports for the same cell. This simplification is the single largest scope
  cut in Phase 1's stress/strain-rate support and should be revisited in
  Phase 2 if plugins need the exact ParaView-equivalent J2 (at the cost of
  the corner-interpolation machinery and its ghost-point dependency).
- **`svCell->svDev.eta`** (effective viscosity) and **`svCell->svDev.APS`**
  (accumulated plastic strain) are both already-computed per-cell scalars —
  free to read, no interpolation. `eta` is scaled by `scal->viscosity`; `APS`
  is dimensionless by construction (a strain, not a stress or rate).
- **Dimensionality**: all quantities passed to the plugin are dimensional,
  using the same scaling factors LaMEM's own output code uses
  (`scal->length`, `scal->temperature` with the `Tshift` offset subtracted —
  see `src/marker.cpp:320` for the exact same formula used when *writing*
  markers to file — `scal->stress`, `scal->strain_rate`, `scal->viscosity`).
  I could not spot-check a live marker's plugin-visible value against a
  ParaView `.vtu` from an actual run (no binary), so this is a code-level
  guarantee (same formulas as the rest of LaMEM's I/O code, applied
  consistently) rather than an empirically cross-checked one.

## Build status (blocking issue — read before attempting to reproduce)

`make mode=opt all -j8` from `src/` fails, but **not** because of anything in
this change. It fails compiling `fastscape.cpp` (and, when checked in
isolation, LaMEM.h itself — see below) with dozens of
`'FILE' does not name a type` / `'FILE' was not declared in this scope`
errors originating entirely inside GCC's own system headers:
```
/opt/homebrew/Cellar/gcc@11/11.5.0/lib/gcc/11/gcc/aarch64-apple-darwin23/11/include-fixed/stdio.h:83:8:
  error: 'FILE' does not name a type
```
Root cause, confirmed by isolating it:
- `$PETSC_OPT` (`/Users/kausb/Software/PETSc/petsc-3.23.6/petsc-3.23.6_opt`)
  was configured with `CXX=/Users/kausb/Software/PETSc/.../bin/mpicxx`, whose
  wrapper script hardcodes `CXX="g++-11"` (Homebrew GCC 11.5.0).
- On this machine's current macOS SDK (Xcode 21 / Darwin 25.5.0), GCC 11's
  bundled `include-fixed/stdio.h` shim is incompatible with the SDK's
  `wchar.h`/`_wchar.h` chain: as soon as any translation unit pulls in
  `<iostream>` (which transitively includes `<cwchar>` → SDK `wchar.h` →
  `_wchar.h`), GCC's own `<cstdio>`/`stdio.h` fixed-include machinery breaks
  and `FILE` becomes undefined for the rest of the translation unit.
- This is **not specific to `fastscape.cpp`**: compiling `src/LaMEM.h` alone
  (`mpicxx -fsyntax-only` after including only `LaMEM.h`, which itself
  `#include`s `<stdio.h>` then `<petsc.h>`) already fails the same way, and
  compiling `src/advect.cpp` (an existing, untouched file) fails identically.
  It affects effectively every `.cpp` file in `src/`.
- Adding `#include <cstdio>` at the top of `fastscape.cpp` (tried, then
  reverted — **not part of the committed change**) does not fix it either;
  the breakage happens as soon as `<iostream>` is pulled in by anything
  earlier in the chain, before `<cstdio>` would have a chance to matter.
- All other locally installed PETSc "opt" builds
  (`petsc-3.22.5_opt`, `petsc-3.23.6/petsc-3.23.6_opt`) hardcode the same
  `mpicxx → g++-11` wrapper, so switching `PETSC_OPT` does not help.
- LaMEM's own CI (`.github/workflows/CI.yml`) uses `g++-12`/`gcc-12` on what
  is presumably Ubuntu Linux, a different OS/SDK combination that does not
  hit this GCC-11-vs-macOS-SDK interaction — so this is very unlikely to
  reproduce in CI or on a normally-provisioned Linux box, only on this
  specific local macOS environment with this specific PETSc install.

I deliberately did **not** patch around this (no edits to `fastscape.cpp`,
`LaMEM.h`, the Makefile, or `$PETSC_OPT`/its `mpicxx` wrapper are part of the
committed change) since it is a pre-existing, environment-level defect
unrelated to this feature, and the task instructions were explicit about
reporting such blockers precisely rather than silently working around them.
Repointing `$PETSC_OPT`'s `mpicxx` wrapper at `g++-12`/`g++-13` (also
installed via Homebrew on this machine) while linking against PETSc's
already-built `g++-11`-compiled static archives (`libpetsc.a`, `libmpi.a`,
etc.) is an ABI/ODR risk for numerically sensitive code and was not attempted
silently; it is the most promising next step for actually unblocking this
locally (a full PETSc rebuild with a matching compiler would be the safe
fix), but is outside the scope of this LaMEM-only worktree.

### Consequence for the test plan (item 4)

Because no `bin/opt/LaMEM` could be produced:
- **4a** (t16 unchanged without `-phase_transition_lib`): not run.
- **4b** (1-rank run with `-phase_transition_lib`, changed-marker count,
  timing): not run.
- **4c** (4-rank run, consistency): not run.
- **4d** (Julia reimplementation of the Constant transition vs. built-in,
  full t16 run + `.expected` comparison): the *reimplementation itself* was
  built and functionally verified against LaMEM's exact transition rule
  (see (v) above); the *end-to-end LaMEM run and comparison* was not run.

All of this is purely blocked by the local toolchain, not by anything in the
design or the new code, and should be re-run as the very next step once a
working `$PETSC_OPT` (or a from-source PETSc build with a working C++
compiler on this machine) is available. The exact commands to do so are given
below.

## Extended-ABI verification (what *could* be run and was)

Both Julia libraries were compiled with `juliac` into relocatable bundles and
exercised through the exact `dlopen(RTLD_NOW|RTLD_GLOBAL)` +
`jl_init_with_image_handle(handle)` + `dlsym` sequence that
`PhTrPluginLoad`/`PhTrPluginApply` perform (verified via a standalone C
harness, not through LaMEM itself, since no LaMEM binary exists in this
environment):

```
# from the scratchpad's ptspike/ directory
julia +1.13.0 --startup-file=no --project=../juliac_env -e \
  'using JuliaC; JuliaC.main(["--output-lib","libptlib.dylib","--bundle","build", \
    "--trim=safe","--compile-ccallable","--experimental","--project=proj","ptlib.jl"])'

julia +1.13.0 --startup-file=no --project=../juliac_env -e \
  'using JuliaC; JuliaC.main(["--output-lib","libptlib_constant.dylib","--bundle","build_constant", \
    "--trim=safe","--compile-ccallable","--experimental","--project=proj","ptlib_constant.jl"])'
```

Results (extended ABI, `ptlib.jl`, 6 synthetic markers, all starting at
phase 1):
```
changed=3
  marker 0: phase 1 -> 2   (z=-200 < -100, T=900 > 800)
  marker 1: phase 1 -> 1   (z=-50, rule doesn't fire)
  marker 2: phase 1 -> 1   (z=-200 but T=700 <= 800)
  marker 3: phase 1 -> 2   (z=-200 < -100, T=900 > 800)
  marker 4: phase 1 -> 3   (z=0, T=500, but j2_strainrate_cell=5e-14 > 1e-14 threshold)
  marker 5: phase 1 -> 1   (j2_strainrate_cell=1e-16, below threshold)
```
(`STRAINRATE_THRESHOLD = 1e-14 [1/s]`, picked to sit an order of magnitude
above LaMEM's typical `DII_ref` background reference strain rate used in the
t16 setup, `DII_ref = 1e-15`, so the rule fires on visibly faster-than-background
deformation without needing an actual t16 run to calibrate against — this
should be revisited once real strain-rate fields from a live t16 run are
available.)

Results (`ptlib_constant.jl`, reimplementation of PhaseTransition ID 0):
```
changed=3
  marker 0: phase 2 -> 3 (T=1300)
  marker 1: phase 3 -> 2 (T=1100)
  marker 2: phase 3 -> 2 (T=500)
  marker 3: phase 5 -> 5 (T=1300, not eligible)
```

## Files changed/added (this branch)

- `src/phase_transition_plugin.h` (new) — public API + documented ABI.
- `src/phase_transition_plugin.cpp` (new) — implementation.
- `src/LaMEMLib.cpp` (modified) — `#include`, `PhTrPluginLoad`/`Apply`/`Destroy`
  calls wired into `LaMEMLibSolve`.
- `doc/phase_transition_plugin_PHASE1_REPORT.md` (this file).
- Scratchpad only (not part of the LaMEM repo/branch):
  `ptspike/ptlib.jl` (extended ABI + strain-rate rule),
  `ptspike/ptlib_constant.jl` (new — Constant-transition reimplementation),
  `ptspike/build/`, `ptspike/build_constant/` (juliac bundles).

## Exact commands to reproduce (once the build is unblocked)

```bash
# 1. Build LaMEM (currently blocked in this environment — see "Build status")
cd src && make mode=opt all -j8

# 2. t16 unchanged (no plugin)
cd ../test && julia --project=. -e 'include("runtests.jl")' 16

# 3. 1-rank run with the spike plugin (extended ABI)
SCRATCH=/private/tmp/claude-504/-Users-kausb-WORK-LaMEM-LaMEM/2205d3c2-ed59-40c3-b10a-cef03cc82eb1/scratchpad
$PETSC_OPT/bin/mpiexec -n 1 ../bin/opt/LaMEM \
  -ParamFile t16_PhaseTransitions/Plume_PhaseTransitions.dat -nstep_max 30 \
  -phase_transition_lib $SCRATCH/ptspike/build/lib/libptlib.dylib

# 4. 4-rank run, same plugin
$PETSC_OPT/bin/mpiexec -n 4 ../bin/opt/LaMEM \
  -ParamFile t16_PhaseTransitions/Plume_PhaseTransitions.dat -nstep_max 30 \
  -phase_transition_lib $SCRATCH/ptspike/build/lib/libptlib.dylib

# 5. Constant-transition reimplementation vs built-in (edit the .dat to
#    remove/disable PhaseTransition ID 0, or add a second copy of the .dat
#    with it stripped), then:
$PETSC_OPT/bin/mpiexec -n 1 ../bin/opt/LaMEM \
  -ParamFile t16_PhaseTransitions/Plume_PhaseTransitions_NoBuiltinPT0.dat -nstep_max 30 \
  -phase_transition_lib $SCRATCH/ptspike/build_constant/lib/libptlib_constant.dylib
# then compare |Div|_inf / |mRes|_2 against
# test/t16_PhaseTransitions/PhaseTransitions.expected within the existing
# tolerances (rtol=1e-5/1e-2, atol=1e-7/1e-3).
```
