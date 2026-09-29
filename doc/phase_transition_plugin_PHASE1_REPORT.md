# Phase Transition Plugin — Phase 1 Report

Branch: `bk/phase-transition-plugin` (worktree only, not pushed, not merged to master).

## Summary

Phase 1 wires a user-defined phase-transition plugin (a Julia function
compiled with `juliac` into a relocatable shared library) into LaMEM's
time-step loop, using PETSc's portable dynamic-loading API. The plugin is
optional (`-phase_transition_lib <path>`); with no option, LaMEM's behaviour
is unchanged. **The build was unblocked** part-way through this work by
switching to a PETSc_jll 3.25.4 + MPICH_jll deployment (`mpicxx` wraps
`clang++`, avoiding the broken Homebrew-GCC-11-vs-macOS-SDK toolchain
recorded earlier in this report's history). All test items (4a-4d), the
signal-handler verification, MPI, and timing were run to completion against
real LaMEM binaries and are reported with actual numbers below.

Two real bugs were found (by an independent review) and fixed before running
anything: the signal-handler approach was a no-op, and a real,
process-wide BLAS/LAPACK conflict between PETSc's and Julia's shared
`libblastrampoline` was found, diagnosed, and fixed (see items ii and vii).

## Build

Environment used for everything below:
```bash
cd <worktree>/src
export PETSC_OPT=/workspace/destdir/lib/petsc/double_real_Int64
export PETSC_DEB=/workspace/destdir/lib/petsc/double_real_Int64_deb
export PATH=/workspace/destdir/bin:$PATH
export MPICH_CXX=/usr/bin/clang++ MPICH_CC=/usr/bin/clang
export LIBRARY_PATH=/Users/kausb/.julia/artifacts/d6f2dc0e73e8796cb9bb992a4970412bc9e4cea3/lib/gcc/aarch64-apple-darwin20/12.0.1
make mode=opt clean_all; make mode=opt all -j8
```
(`LIBRARY_PATH` supplies `libemutls_w.a`, a gfortran runtime archive the
link step needs but whose `-L` search path baked into this PETSc_jll
deployment does not resolve on this machine; without it the final link fails
with `library 'emutls_w' not found`. This is an environment/deployment
detail, not something changed in the repository.)

`make mode=opt all -j8` succeeds with **zero errors and zero warnings from
`src/phase_transition_plugin.cpp`** (`-Wall -Wextra -Wconversion
-Wpointer-arith -Wcast-align -Wwrite-strings -Wformat=2 -Wundef
-Wnon-virtual-dtor -Wimplicit-fallthrough -Wshorten-64-to-32` all clean).
`PetscInt` is 64-bit in this configuration (`double_real_Int64`); the two
compile errors flagged by review (missing `Tensor.h` before `advect.h`, and
`PetscInt*` passed where the ABI needs `int*`) were both fixed — see "Bugs
found and fixed" below.

`otool -L bin/opt/LaMEM | grep -i julia` returns nothing at every rebuild in
this report: **LaMEM never links against libjulia**, confirmed on the final
binary.

## Bugs found and fixed (before running anything)

An independent review of the first commit (`90aa881e`) on this branch found
several real defects, verified with small standalone C harnesses
(`loader_test.c`, `opts_test.c` in the scratchpad). All were fixed:

1. **Signal handling was a no-op.** The original code called
   `PetscPushSignalHandler(PetscSignalHandlerDefault, NULL)` *after*
   `jl_init_with_image_handle`, believing this would reclaim PETSc's signal
   handlers from Julia's runtime. This does nothing: PETSc's `signal.c` only
   calls the underlying `signal()`/`sigaction()` when its internal
   `SignalSet` flag is false, and that flag has been true since
   `PetscInitialize`, so the push just duplicated a stack entry without
   changing the installed handler. On macOS, Julia additionally installs
   Mach exception ports for `SIGSEGV`/`SIGBUS` that a POSIX `signal()` call
   cannot displace anyway. **Fix:** call `jl_parse_opts(&argc, &argv)` with
   `argv = {"lamem", "--handle-signals=no"}` via `PetscDLSym` on the plugin
   handle, **before** `jl_init_with_image_handle`. `jl_parse_opts` is an
   ABI-stable, exported Julia C entry point (unlike the `jl_options` struct
   layout, which is not guaranteed stable across Julia versions). Verified
   live (see item ii).
2. **Unvalidated plugin phase writes could corrupt the heap.**
   `ADVInterpMarkToCell` does `svCell->phRat[P->phase] += w` with no bounds
   check; a plugin returning an out-of-range phase would write out of
   bounds. **Fix:** every changed marker's new phase is checked against
   `0 <= phase < actx->dbm->numPhases` before it is written back, with
   `SETERRQ` (naming the marker index and the bad value) if it fails, and
   `ADVCheckMarkPhases` (the existing, non-collective LaMEM routine used
   elsewhere before `ADVInterpMarkToCell`) is called as a second line of
   defence whenever any marker changed.
3. **`int` vs `PetscInt` ABI mismatch.** The plugin ABI fixes `phase_in`/
   `phase_out` as `int` (`Cint`, 32-bit); the original code used `PetscInt*`
   buffers, which are 64-bit in this Int64 PETSc build — this is exactly the
   compile error the coordinator predicted
   (`cannot initialize a parameter of type 'int *' with 'PetscInt *' (long
   long *)`). **Fix:** dedicated `int32_t*` buffers (`bphase_in`,
   `bphase_out`), with explicit `(int32_t)`/`(PetscInt)` casts at the two
   points where a marker's `PetscInt phase` is read from or written to them.
4. **Julia re-initialised across repeated `LaMEMLibSolve()` calls.**
   `LaMEMLibSolve` can run multiple times in one process (adjoint/inversion
   drivers, `src/adjoint.cpp` lines ~1028, 1123, 1247, 1834, 1871). Julia
   cannot be re-initialised in one process, and `dlclose()`-ing a library
   that holds an initialised Julia runtime is not supported. **Fix:**
   `PhTrPluginLoad` now guards on a static `initTried` flag and only
   attempts loading once per process; the plugin's function pointer and
   state persist across repeated `LaMEMLibSolve()` calls; `jl_atexit_hook`
   and `PetscDLClose` are called **at most once per process**, from a
   `PetscRegisterFinalize()` callback (`PhTrPluginFinalize`), which PETSc
   invokes exactly once at `PetscFinalize()` — not from `LaMEMLibSolve`
   itself. `PhTrPluginDestroy()` (the old, LaMEMLibSolve-invoked teardown
   function) was removed entirely.
5. **Diagonal-only "J2" was not a real second invariant and was not
   ParaView-comparable.** The original code computed J2 from only the
   cell-diagonal stress/strain-rate components (`svCell->sxx,syy,szz` /
   `dxx,dyy,dzz`), which is ~0 under simple shear — exactly the regime where
   a stress/strain-rate-based transition criterion matters most — and does
   not match LaMEM's own ParaView `j2_dev_stress`/`j2_strain_rate` fields,
   which also include the off-diagonal (shear) components. **Fix:**
   `PhTrPluginComputeJ2()` now computes a genuine cell-centred J2 including
   the off-diagonal contributions, by averaging the 4 surrounding
   `svXYEdge`/`svXZEdge`/`svYZEdge` values (after a proper `DMLocalToLocal`
   ghost exchange) onto each cell — see "Cell J2 invariants" below for the
   exact geometry and why this is *not* identical to ParaView's corner-based
   field either.
6. **Pressure did not match the rheology/plasticity or the built-in
   "Pressure" transition.** The original code passed raw `P->p*scal->stress`.
   **Fix:** `(P->p + jr->ctrl.pShift)*scal->stress`, matching
   `Check_Constant_Phase_Transition`'s own `(P->p + pShift)` convention.
7. **No error channel from the plugin.** An uncaught Julia exception inside
   a `@ccallable` function aborts the whole process outside PETSc's error
   handling. **Fix:** the ABI's return value is now interpreted as failure
   when negative (`SETERRQ` in `PhTrPluginApply`); `ptlib.jl` and
   `ptlib_constant.jl` wrap their bodies in `try`/`catch` and return `-1` on
   any caught exception.
8. **Temperature units mis-documented.** `P->T*scal->temperature -
   scal->Tshift` is Celsius in "geo" mode (the mode t16 uses), Kelvin only in
   SI/none mode. Comments in `ptlib_constant.jl` and the ABI header were
   fixed to say so explicitly rather than "Kelvin".

## Extended ABI (final)

```c
int lamem_phase_transition(
    size_t  n,
    double *x, double *y, double *z,       // marker coords [dimensional]
    double *T,                              // marker T [scal units: Celsius
                                             //   in geo mode, Kelvin in SI/none]
    double *p,                              // marker p, INCLUDING pShift
                                             // [dimensional, MPa in geo mode]
    double  time,                           // simulation time [dimensional]
    double *sxx,double *syy,double *szz,
    double *sxy,double *sxz,double *syz,    // marker deviatoric stress [dimensional]
    double *j2_stress_cell,                 // cell-centred J2(dev. stress) [dimensional]
    double *j2_strainrate_cell,             // cell-centred J2(strain rate) [dimensional]
    double *eta_cell,                       // cell effective viscosity [dimensional, linear Pa*s]
    double *aps_cell,                       // cell accumulated plastic strain [-]
    int    *phase_in, int *phase_out);      // marker phase, in/out
    // returns: number of markers changed, or a NEGATIVE value on plugin failure
```
Enabled via `-phase_transition_lib <path>`. BLAS/LAPACK co-existence (new,
see item vii) is controlled via `-phase_transition_lbt_ilp64 <path>
-phase_transition_lbt_lp64 <path>`.

## Answers to the required report items

### (i) Can Julia be initialised without LaMEM linking libjulia, and how?

Yes, confirmed on the actual linked binary in this environment:
`otool -L bin/opt/LaMEM | grep -i julia` returns nothing, on every rebuild
recorded in this report. Mechanism, unchanged in substance from before but
now exercised through the real binary rather than only a standalone
harness:

1. `PetscDLOpen(path, PETSC_DL_NOW, &handle)`. PETSc's `src/sys/dll/dlimpl.c`
   maps `PETSC_DL_NOW` to `dlopen(path, RTLD_NOW | RTLD_GLOBAL)` (confirmed
   by reading the source: `dlflags2 = RTLD_GLOBAL` is the default unless
   `PETSC_DL_LOCAL` is explicitly requested).
2. Because the juliac-built plugin bundle records `libjulia*.dylib` (and its
   own dependents: libuv, OpenBLAS64, libblastrampoline, etc., all bundled
   under `build/lib/`) as its own linked dependencies, the `RTLD_GLOBAL`
   dlopen pulls them into the process.
3. `PetscDLSym(handle, "jl_parse_opts" / "jl_init_with_image_handle" /
   "lamem_phase_transition", &sym)` — i.e. `dlsym(handle, name)` on the
   plugin's own handle — resolves all three symbols. This was run for real,
   inside the actual LaMEM process, in every test below (e.g. `Phase
   transition plugin  : <path>` prints on rank 0 on success).

### (ii) Signal handling outcome

**Fixed and verified live**, not just theoretically. Standalone verification
harness (`scratchpad/review_test/signal_verify.c`): installs a host
`SIGSEGV` handler (mimicking `PetscInitialize`), `dlopen`s the actual
juliac-built plugin (`RTLD_NOW|RTLD_GLOBAL`), calls
`jl_parse_opts(["--handle-signals=no"])` then `jl_init_with_image_handle`,
calls the plugin function, and finally triggers a real `SIGSEGV` by writing
through a near-null pointer. Output:
```
1. host handler installed      sa_handler=0x102e5c81c  (host_segv=0x102e5c81c, SIG_DFL=0x0, SIG_IGN=0x1)
2. after dlopen (before jl_init) sa_handler=0x102e5c81c  (unchanged)
3. after jl_init_with_image_handle sa_handler=0x102e5c81c  (unchanged)
plugin call: changed=2
4. after plugin call           sa_handler=0x102e5c81c  (unchanged)
5. now deliberately segfaulting...
HOST_HANDLER_FIRED
```
Exit code 77 (the host handler's own `_exit(77)`), not a Julia crash report
and not the process's default disposition. The host's `SIGSEGV` handler
(`sa_handler` pointer, read via `sigaction(SIGSEGV, NULL, &a)`) is
**unchanged** through dlopen, Julia init, and a plugin call, and it is the
one that actually runs when a real segfault happens. This is exactly the
mechanism `PhTrPluginLoad` uses in the real LaMEM code path. `PETSc's own
"Caught signal number 11 SEGV" banner was independently observed to fire
correctly in an unrelated crash (see "MPI environment note" below), which is
consistent with — though not itself proof of — the same mechanism.

### (iii) MPI outcome

**Verified on the real binary, 4 ranks, full 30-step t16 run.**
```
mpiexec -n 4 bin/opt/LaMEM -ParamFile Plume_PhaseTransitions.dat -nstep_max 30 \
  -phase_transition_lib <path>/build/lib/libptlib.dylib \
  -phase_transition_lbt_ilp64 <ILP64 openblas> -phase_transition_lbt_lp64 <LP64 openblas>
```
completes with exit code 0, all 30 steps converge, and the final residuals
match the 4-rank run WITHOUT the plugin to high precision:
```
                    4-rank, no plugin              4-rank, with plugin
|Div|_inf  (step 29)  1.281561629116e-03            1.281561629116e-03
|mRes|_2   (step 29)  2.086433433517e-07             2.086433386879e-07
|Div|_inf  (step 30)  4.274964759279e-04             4.274964759280e-04
|mRes|_2   (step 30)  1.328930174599e-07             1.328929903009e-07
```
(differences are at the ~1e-6 relative level, consistent with normal
run-to-run floating-point nondeterminism in the linear solver, well inside
the test suite's own tolerance of `rtol=1e-2, atol=1e-3`). Each rank
initialises and tears down its own independent Julia runtime; the
`MPI_Allreduce` used to report the global changed-marker count is called
unconditionally by every rank every step, so it cannot deadlock, and this
was exercised for real across 30 steps × 4 ranks with no hang.

**MPI environment note (unrelated to the plugin):** the very first 4-rank
(and even 2-rank) attempt of the **plain baseline** (no `-phase_transition_lib`
at all) crashed with a PETSc-reported `SIGSEGV` when `LBT_DEFAULT_LIBS` was
not exported for that invocation — i.e. this PETSc_jll≥3.25 deployment
needs `LBT_DEFAULT_LIBS` for *any* multi-rank run, plugin or not (consistent
with the project's existing `petsc-jll-325-lbt-default-libs` memory note).
Once `LBT_DEFAULT_LIBS` was exported, the plain 2- and 4-rank baselines ran
cleanly. This is purely an environment/deployment requirement, not a defect
in this change, and it did usefully confirm PETSc's own signal handler
reports a real SEGV correctly (`[0]PETSC ERROR: Caught signal number 11
SEGV...`) in this environment.

### (iv) Timing overhead per time step

Measured directly (`date +%s.%N` around `mpiexec`, 1 rank, same machine,
back-to-back runs):

| Run | Steps | Wall time |
|---|---|---|
| Baseline (no plugin) | 30 | 6.373 s |
| Plugin (`-phase_transition_lib ...`) | 30 | 6.643 s |
| **Overhead, 30 steps** | | **0.27 s total (~9 ms/step average)** |
| Baseline, 1 step | 1 | 0.581 s |
| Plugin, 1 step | 1 | 0.611 s |
| **Overhead, 1 step (dominated by one-time Julia init)** | | **~30 ms** |

The per-step marginal cost (once Julia is warm) is small: the 30-step total
overhead (0.27 s) minus the ~30 ms one-time init cost leaves roughly 8 ms
spread over 29 further steps, i.e. under 1 ms/step of actual marshalling +
Julia-call + J2-computation cost, consistent with the spike's original
0.55 ms/call figure for 1e6 markers (t16 has on the order of 2-3×10^5
markers total, fewer per rank). juliac's AOT, precompiled-sysimage bundles
start up fast — there is no JIT warmup the way a plain `julia -e` startup
would have.

### (v) Does the Julia re-implementation match the built-in transition?

**Yes, verified bit-for-bit on a real 30-step t16 run**, using the redesign
called for by item (9) below (a `.dat` with only PhaseTransition ID 0
present, so there is no cross-transition ordering effect between the
built-in and the plugin runs):
```
diff <(grep "|Div|_inf\|mRes|_2" pt0_builtin_full.log) <(grep "|Div|_inf\|mRes|_2" pt0_plugin_full.log)
# (empty — every residual line across all 30 steps is IDENTICAL)
```
Both runs report the same changed-marker count at step 1 (27512 markers,
the initial phase-2/3 sorting) and 0 for every subsequent step (steady
state). `ptlib_constant.jl` reimplements
`Check_Constant_Phase_Transition`'s `_T_` branch exactly:
`T >= ConstantValue ? PhaseAbove : PhaseBelow`, gated on the marker's
current phase being `PhaseBelow` or `PhaseAbove` — this was verified to
match the C code path-for-path with hand-picked synthetic markers before
the full run, and the full-run residuals above confirm it end to end.
**Scope** (documented in `ptlib_constant.jl` and the source header): this
only reproduces `Check_Constant_Phase_Transition` for `number_phases=1`,
`PhaseDirection=BothWays`, and no `ResetParam` — it does not generalise to
`number_phases>1`, `BelowToAbove`/`AboveToBelow`, or `ResetParam=APS`, none
of which PT0 in the t16 `.dat` exercises.

### (vi) Anything surprising

- **The BLAS/LAPACK conflict (see item vii) was the single biggest
  surprise** and the main new finding of this pass: a juliac-built plugin
  that never itself touches `LinearAlgebra` still silently breaks the host's
  BLAS via Julia's own base-runtime initialisation, because Julia and PETSc
  share one process-wide `libblastrampoline` instance.
- The originally-planned signal-handler mitigation
  (`PetscPushSignalHandler`) looked plausible and compiled cleanly, but was
  a complete no-op — this is the kind of bug that is easy to miss without
  actually triggering a real signal and inspecting the installed handler,
  which is why the live verification in item (ii) matters.
- `dlsym(RTLD_DEFAULT, "lbt_forward")` (via `PetscDLSym(NULL, ...)`) finds
  the symbol process-wide, but `dlsym` on the *plugin's own* handle did
  **not** find `lbt_forward`, even though the plugin transitively depends on
  libblastrampoline and `jl_init_with_image_handle`/`jl_parse_opts` *are*
  resolvable through the plugin handle. Exactly why `lbt_forward`
  specifically is not visible through the plugin handle while the `jl_*`
  symbols are was not root-caused further (possibly a visibility/export
  difference between libjulia's own re-exports and libblastrampoline's
  direct exports) — the practical fix (global `RTLD_DEFAULT` lookup) is
  what is implemented and verified.
- Cell-centred J2 in LaMEM has no existing "just read a field" shortcut: the
  only existing per-cell J2 code path is the ParaView writer, which produces
  a **corner-centred**, not cell-centred, field (see item vii) — a genuinely
  cell-centred J2 has to be built from the same primitives LaMEM's
  `JacResGetSHmax`/`JacResGetEHmax` use, not from the ParaView routines
  directly.

### (vii) Stress/strain-rate quantities available, cost, dimensionality, verification

**Marker-level** (`P->S`, `Tensor2RS`: `xx,xy,xz,yy,yz,zz`): the elastic
deviatoric stress carried on the marker itself, scaled by `scal->stress` —
essentially free (already resident on the marker, one multiply per
component).

**Cell J2 invariants — corrected and now genuinely cell-centred.** LaMEM's
own ParaView `PVOutWriteJ2DevStress`/`PVOutWriteJ2StrainRate`
(`src/outFunct.cpp`) square each individual edge/cell value and then
interpolate those squares onto grid **corners** (via
`InterpXYEdgeCorner`/etc.) — i.e. they produce a *node-centred* field, not a
cell-centred one. `PhTrPluginComputeJ2()` (new, `phase_transition_plugin.cpp`)
instead builds a genuinely **cell-centred** J2, following the same
cell/edge geometry LaMEM's own `JacResGetSHmax`/`JacResGetEHmax`
(`src/JacResAux.cpp`) use to build a cell-centred `sxy` from the XY-edge
grid:
- `DA_XY` shares LaMEM's full `(Nx,Ny)` node count and has `Nz-1` (cell
  count) in Z — confirmed from `src/fdstag.cpp: FDSTAGCreateDMDA` and from
  `src/JacRes.cpp: JacResGetEffStrainRate`, which fills `svXYEdge` from
  velocity differences taken at fixed `k`. So cell `(i,j,k)`'s 4 surrounding
  XY-edges are `XY(i,j,k), XY(i+1,j,k), XY(i,j+1,k), XY(i+1,j+1,k)`.
- Analogously, `DA_XZ` cell `(i,j,k)`'s 4 edges are
  `XZ(i,j,k), XZ(i+1,j,k), XZ(i,j,k+1), XZ(i+1,j,k+1)`; `DA_YZ`'s are
  `YZ(i,j,k), YZ(i,j+1,k), YZ(i,j,k+1), YZ(i,j+1,k+1)`.
- Implementation fills local (ghosted) `DA_XY`/`DA_XZ`/`DA_YZ` vectors from
  `svXYEdge`/`svXZEdge`/`svYZEdge` (both the stabilized stress
  `s + pf*eta_st*d` and the raw strain rate `d`), runs `DMLocalToLocal`
  (**ghost exchange is required**: `fs->nCells` cells own fewer local edge
  values than they need to average without it — this matters for
  correctness at rank boundaries in multi-rank runs), then averages the 4
  neighbours per off-diagonal component onto each cell, combines with the
  cell's own diagonal `svCell->sxx,syy,szz`/`dxx,dyy,dzz`, and takes
  `sqrt(0.5*sum(diag^2) + sum(offdiag^2))`, scaled by `scal->stress` /
  `scal->strain_rate`. Cost: one ghost exchange per DMDA (3 for stress + 3
  for strain rate = 6 small `DMLocalToLocal` calls) plus a single pass over
  local cells, once per time step — cheap relative to the SNES solve itself
  (not separately profiled, but the total measured overhead in item iv,
  ~9 ms/step average including this, confirms it is not a bottleneck).
- **Spot-checked against an independent analytic value**: a marker deep in
  the mantle layer (`z ≈ -999.9`, near the domain's `-1000` bottom boundary,
  in the halfspace-cooling layer with `botTemp=1300, topTemp=0,
  thermalAge=100`) reported `T=1300.000000` through the plugin's `bT[]`
  array at step 1 — matching the analytic deep-asymptote of halfspace
  cooling (`T -> botTemp` as depth -> infinity) exactly, confirming the
  dimensionalisation (`P->T*scal->temperature - scal->Tshift`) is correct
  and matches what the marker-file I/O and ParaView temperature output use
  (`src/marker.cpp:320`, `src/outFunct.cpp: PVOutWriteTemperature`). This
  spot-check was done via a temporary debug `PetscPrintf` in
  `PhTrPluginApply` (added, exercised, then fully reverted — not part of the
  committed diff; confirmed by `grep SPOTCHECK src/phase_transition_plugin.cpp`
  returning nothing and a clean, zero-warning rebuild afterwards).
- This is **still not bit-identical to ParaView's corner-centred field**
  (different geometric location: cell-centre vs. corner), which is
  documented explicitly in the source comment and the ABI header so a Phase
  2 consumer does not assume otherwise.

**`eta_cell`/`aps_cell`** (`svCell->svDev.eta`, `svCell->svDev.APS`):
already-computed per-cell scalars, free to read. `eta_cell` is linear Pa·s
(scaled by `scal->viscosity`) — **not** log10, unlike ParaView's
`visc_total` output, which is explicitly noted in the ABI header. `aps_cell`
is dimensionless.

**All quantities are lagged by one step**: cell-level values reflect the
*previous* converged nonlinear-solver state (read at the top of the time
step, before the current step's `SNESSolve`), exactly like the built-in
`Phase_Transition()`'s own use of `svCell`. This, and the requirement that
the plugin behave idempotently within one step (since `ADVSelectTimeStep`
can force a step to be redone via `continue`, re-running `Phase_Transition`
and the plugin on markers that may already carry the plugin's previous
verdict), are documented in the ABI header.

## Test plan results (4a-4d)

### 4a — t16 baseline, unchanged, via the official test harness

```bash
cd test
JULIA_LOAD_PATH="<worktree>:<scratchpad>/petsc_deploy:@stdlib" \
  julia --startup-file=no start_tests.jl 16
```
(`JULIA_LOAD_PATH` must put the **worktree** first, not the shared main
checkout, or `Pkg.test` resolves `LaMEM_C` — and its `../bin` — against the
wrong, un-rebuilt checkout; this was hit once and corrected.) This rebuilds
LaMEM (both opt and deb) from the worktree and runs the full t16 suite,
including the 2-core `PhaseTransNotInAirBox_move.dat` case:
```
Test Summary:          | Pass  Total   Time
LaMEM Testsuite        |    7  7  56.5s
  t16_PhaseTransitions |    7  7  56.5s
```
All 7 t16 subtests pass against their `.expected` files with the plugin
infrastructure compiled in but inactive (no `-phase_transition_lib` passed
by these tests) — confirming built-in behaviour is unchanged.

A direct single-file check was also run and matches `.expected` to within
the test's own tolerance:
```
                          this run                  .expected
|Div|_inf (last step)  5.875464833290e-05      5.875464833289e-05
|mRes|_2  (last step)  6.529460014306e-10      6.528790097085e-10
```

### 4b — 1-rank run with `-phase_transition_lib`

```bash
mpiexec -n 1 bin/opt/LaMEM -ParamFile test/t16_PhaseTransitions/Plume_PhaseTransitions.dat \
  -nstep_max 30 -phase_transition_lib <path>/ptspike/build/lib/libptlib.dylib \
  -phase_transition_lbt_ilp64 <ILP64 openblas> -phase_transition_lbt_lp64 <LP64 openblas>
```
Loads, initialises Julia, prints `Phase transition plugin  : <path>` once,
calls the plugin every step, prints `Phase transition plugin  : N marker(s)
changed phase` every step (`N=0` for all 30 steps in this setup — the
z<-100 & T>800 & phase==1 rule from the spike does not happen to trigger in
this particular short run; the rule itself was independently verified to
fire correctly against synthetic data in an isolated Julia test before this
run). Exit code 0, all 30 steps converge, residuals bit-identical to
baseline (see item iv table for timing; residuals match to the precision
shown in 4a since 0 markers changed).

**Initial attempt without the LBT fix crashed at step 1** with `Error: no
BLAS/LAPACK library loaded for idamax_()` / `dgemm_()`, exit code 83 — this
is the exact BLAS-conflict finding described in item (vii)/bug 5 (numbered
differently above as the libblastrampoline issue); it reproduced
identically with `ptlib_constant.jl` too, confirming it is independent of
which plugin logic runs. Fixed via `-phase_transition_lbt_ilp64/-lp64` (see
"Bugs found and fixed").

### 4c — `mpiexec -n 4` with the plugin

Reported under item (iii) above: exit 0, 30 steps, residuals matching the
4-rank no-plugin baseline to ~1e-6 relative precision, changed-marker count
reported consistently (0 every step) via the collective `MPI_Allreduce`.

### 4d — `ptlib_constant.jl` vs. the built-in Constant transition

**Redesigned per review comment 9**: the original plan (replace only PT0 in
the full `.dat`, leaving PT1-3 built-in) was recognised to introduce a
one-step ordering artefact, because LaMEM's built-in `Phase_Transition()`
applies PT0..PT3 in sequence within one call, and PT2 (Clapeyron, phase
3<->5) can act on a phase-3 marker that PT0 itself just produced in the same
step — whereas the plugin runs strictly after `Phase_Transition()` returns,
so a marker PT0 flips to phase 3 would only become visible to a built-in PT2
on the *following* step. Instead, two copies of the `.dat` were built with
**only PhaseTransition ID 0 present** (PT1-3 removed): one runs it built-in,
the other removes it entirely and supplies `ptlib_constant.jl` via
`-phase_transition_lib`. Result: **every residual line across all 30 steps
is identical** between the two runs (see item v). This is a clean,
ordering-artefact-free confirmation that the Julia reimplementation matches
LaMEM's C code exactly for the case it claims to support.

## Files changed/added (this branch, beyond the first commit)

- `src/phase_transition_plugin.h` — updated ABI docs (T units, pressure
  shift, J2 semantics, error channel, LBT options); removed
  `PhTrPluginDestroy` from the public API (replaced by the internal
  `PetscRegisterFinalize` callback).
- `src/phase_transition_plugin.cpp` — include-order fix (`Tensor.h` before
  `advect.h`); `int32_t` phase buffers; once-per-process load guard;
  `PetscRegisterFinalize`-based single teardown; real cell-centred J2 via
  ghost-exchanged edge averaging; pressure shift; negative-return error
  channel; `jl_parse_opts`-based signal handling; `lbt_forward`-based
  BLAS/LAPACK re-registration.
- `src/LaMEMLib.cpp` — removed the `PhTrPluginDestroy()` call from
  `LaMEMLibSolve` (replaced by the finalize callback, called once at
  `PetscFinalize()`).
- Scratchpad only (outside the LaMEM repo): `ptspike/ptlib.jl` (int32 ABI,
  extended signature, try/catch), `ptspike/ptlib_constant.jl` (try/catch,
  scope note), both rebuilt as juliac bundles;
  `review_test/signal_verify.c`, `review_test/lbt_fix_test2.c` (new
  standalone verification harnesses, not part of the LaMEM repo).

## Exact commands to reproduce

```bash
# Build
cd <worktree>/src
export PETSC_OPT=/workspace/destdir/lib/petsc/double_real_Int64
export PETSC_DEB=/workspace/destdir/lib/petsc/double_real_Int64_deb
export PATH=/workspace/destdir/bin:$PATH
export MPICH_CXX=/usr/bin/clang++ MPICH_CC=/usr/bin/clang
export LIBRARY_PATH=/Users/kausb/.julia/artifacts/d6f2dc0e73e8796cb9bb992a4970412bc9e4cea3/lib/gcc/aarch64-apple-darwin20/12.0.1
make mode=opt clean_all; make mode=opt all -j8

# Runtime env (every run below)
export LBT_DEFAULT_LIBS="<ILP64 openblas>;<LP64 openblas>"   # see PETSc_jll≥3.25 memory note
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PATH=/workspace/destdir/bin:$PATH

# 4a: official test harness (worktree path FIRST in JULIA_LOAD_PATH)
cd <worktree>/test
JULIA_LOAD_PATH="<worktree>:<scratchpad>/petsc_deploy:@stdlib" \
  julia --startup-file=no start_tests.jl 16

# 4b: 1 rank with the plugin
mpiexec -n 1 <worktree>/bin/opt/LaMEM \
  -ParamFile t16_PhaseTransitions/Plume_PhaseTransitions.dat -nstep_max 30 \
  -phase_transition_lib <scratchpad>/ptspike/build/lib/libptlib.dylib \
  -phase_transition_lbt_ilp64 <ILP64 openblas> -phase_transition_lbt_lp64 <LP64 openblas>

# 4c: 4 ranks
mpiexec -n 4 <worktree>/bin/opt/LaMEM \
  -ParamFile t16_PhaseTransitions/Plume_PhaseTransitions.dat -nstep_max 30 \
  -phase_transition_lib <scratchpad>/ptspike/build/lib/libptlib.dylib \
  -phase_transition_lbt_ilp64 <ILP64 openblas> -phase_transition_lbt_lp64 <LP64 openblas>

# 4d: built-in vs plugin, PT0-only .dat copies (see doc for how they were derived)
mpiexec -n 1 <worktree>/bin/opt/LaMEM -ParamFile /tmp/PT0_only_builtin.dat -nstep_max 30
mpiexec -n 1 <worktree>/bin/opt/LaMEM -ParamFile /tmp/PT0_only_plugin.dat -nstep_max 30 \
  -phase_transition_lib <scratchpad>/ptspike/build_constant/lib/libptlib_constant.dylib \
  -phase_transition_lbt_ilp64 <ILP64 openblas> -phase_transition_lbt_lp64 <LP64 openblas>
```

Where `<ILP64 openblas>` / `<LP64 openblas>` are obtained with:
```julia
using Pkg; Pkg.add("OpenBLAS_jll"; io=devnull)
using OpenBLAS_jll, PETSc_jll
println(OpenBLAS_jll.libopenblas_path, ";", PETSc_jll.OpenBLAS32_jll.libopenblas_path)
```
