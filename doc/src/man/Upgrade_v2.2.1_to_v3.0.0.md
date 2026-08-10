# Upgrading from v2.2.1 to v3.0.0

A practical guide for updating your input files, solver setups, tests, and any local
source patches when moving from **LaMEM v2.2.1** (released 2026‑04‑08) to
**LaMEM v3.0.0** (initial release 2026‑05‑06, upstream PR
[#65](https://github.com/UniMainzGeo/LaMEM/pull/65)).

> **Targets the current v3.0.0 development line, not just the release tag.** This guide
> reflects HEAD commit `eaf75d4b` (2026‑06‑15), which includes the post‑release fixes
> merged after the tag (PRs #66, #68, #69, #70 — PBC hardening, integer‑type cleanup,
> passive‑tracer output, and the restored `t36` erosion test). Changes made since the
> `v3.0.0` tag are flagged inline with **_(post‑release)_**.
>
> This guide was produced by diffing the two source trees directly and by running
> unmodified v2.2.1 input files under the v3.0.0 binary to observe real behavior.
> File/line references point into the current v3.0.0 source unless noted.

---

## 0. TL;DR — the migration checklist

If you only read one section, read this.

1. **Rewrite your solver block.** The old top-level `SolverType`, `DirectSolver`,
   `DirectPenalty`, and all `MG*` parameters are **gone**. Replace them (and most of
   your hand-written `<PetscOptionsStart>` solver flags) with a new
   `<SolverOptionsStart>…<SolverOptionsEnd>` block. **→ [§4](@ref "4. Solver options: the headline change")**
   ⚠️ These removed parameters are **silently ignored** — your model will still run,
   but with v3.0.0 *default* solver settings, not yours. See the pitfall in [§8](@ref "8. Pitfalls & troubleshooting").
2. **2D setups: change `nel_y = 1` → `nel_y = 2`.** Fewer than 2 cells in any
   direction is now a hard error. **→ [§5.1](@ref "5.1 Minimum 2 cells per direction (breaking)")**
3. **Graded meshes: drop non-unit `bias_*`.** Bias ratios ≠ 1.0 are now rejected. **→ [§5.2](@ref "5.2 Mesh bias ratios deprecated")**
4. **Out-of-plane background shear removed.** `exz_strain_rates` / `eyz_strain_rates`
   (and their `*_time_delims`) no longer exist. **→ [§5.3](@ref "5.3 Background strain rates reduced")**
5. **Gravity-survey feature removed** entirely (`gravity.cpp/h` deleted). **→ [§5.4](@ref "5.4 Gravity-survey feature removed")**
6. **Paths changed:** `input_models/` → `examples/`; `matlab/`, `utils/`, top-level
   `docs/` removed; new `info/` reference dir. **→ [§3](@ref "3. Repository layout changes")**
7. **Tests renamed:** `t1_…` → `t01_…`, `.expected` files lost their `-pN` suffix.
   Update any downstream CI. **→ [§6](@ref "6. Tests & examples")**
8. **If you patched the C++ source:** every `CHKERRQ(ierr)` is now `PetscCall(...)`,
   the build enforces `-std=c++17` with strict warnings, and needs PETSc 3.19–3.25.
   **→ [§7](@ref "7. For people who patched the source")**

### Quick self-check when migrating a file

The fastest way to migrate one `.dat` file: work top to bottom. Items 1–7 cover every
breaking change and both pitfalls; item 8 is the runtime confirmation.

1. **Solver block** — search for `SolverType`, `DirectSolver`, `DirectPenalty`, `MG` — if
   present, rewrite the solver block (§4). ⚠️ Otherwise silently ignored (Pitfall #1) —
   this is the dangerous one.
2. **Grid resolution** — grep for `nel_x`/`nel_y`/`nel_z = 1` — bump any direction to ≥ 2
   (§5.1). Otherwise: `Less than two cells are specified…`.
3. **Mesh grading** — grep for `bias_` — drop non-unit ratios; use uniform `nseg_*`
   segments (§5.2). Otherwise: `Non-unit bias ratios are deprecated…`.
4. **Background shear** — grep for `exz_`, `eyz_` — the out-of-plane components are gone;
   remove them (§5.3). Otherwise: required-parameter-not-found on `exz_strain_rates`.
5. **Linear + thermal** — if you used `-snes_type ksponly` together with thermal diffusion
   (`act_temp_diff = 1`), replace it with `set_linear_problem = 1` (or `-snes_max_it 1`,
   `src/nlsolve.cpp:117`) — the old combination now aborts (§8, Pitfall #2).
6. **Periodic BCs** — if you set `periodic`, check it against the compatibility guards:
   no-slip bottom, *no* no-slip left/right, no `xx`/`xy` background strain, no `bvel_face`,
   no `fix_phase`/`fix_cell` (§5.5). Otherwise: `Periodic condition is incompatible with …`.
7. **Gravity survey** — grep for any gravity-survey options; the feature was removed with
   **no replacement** (§5.4).
8. **Verify at runtime** — run once with `monitor_solvers = 1` and confirm the
   penalty/tolerances/solver in the PETSc option dump match what you expect.

---

## 1. Requirements & build

| | v2.2.1 | v3.0.0 |
|---|---|---|
| PETSc | 3.19-ish | **3.19 – 3.25** (some features need 3.24/3.25) |
| C++ standard | not pinned in Makefile | **`-std=c++17`** enforced |
| Makefile | `Makefile` + `Makefile.in` | **single `Makefile`** (`Makefile.in` removed) |
| Warnings | minimal | **strict** (`-Wall -Wextra -Wconversion -Wshorten-64-to-32 …`) |

Build commands are unchanged in spirit:

```bash
export PETSC_OPT=/path/to/petsc/opt
export PETSC_DEB=/path/to/petsc/deb
cd src
make mode=opt all        # → bin/opt/LaMEM
make mode=deb all        # → bin/deb/LaMEM
make mode=opt dylib      # → lib/opt/LaMEMLib.* (Julia integration)
```

The `coupled_direct` Stokes solver and a few other paths are guarded on PETSc
3.24/3.25 (`src/options.cpp:65`, `src/lsolve.cpp:367`, etc.). On older PETSc in the
supported range, use a `block_*` solver instead.

---

## 2. What's new in v3.0.0 (the "why")

From the upstream release notes (PR #65):

- **Hybrid multigrid framework** — assembled *and* matrix-free levels in one solve.
- **Matrix-free Jacobian operator** (`src/matFree.cpp/h`, `src/matData.cpp/h`).
- **Block direct solver** and improved **2D multigrid coarsening** for quasi-2D models.
- **Periodic boundary conditions** (mechanical *and* thermal), end-to-end.
- **wBFBT preconditioner** (Rudi et al., 2017) — `src/matBFBT.cpp`.
- **Automated solver options & stop tolerances** — a new `<SolverOptionsStart>` DSL
  (`src/options.cpp/h`) that picks sensible PETSc settings for you.

The last point is the single biggest day-to-day change and is why most v2.2.1 input
files need editing.

---

## 3. Repository layout changes

| v2.2.1 | v3.0.0 | Notes |
|---|---|---|
| `input_models/` | `examples/` | renamed |
| `input_models/BuildInSetups/` | `examples/BuiltInSetups/` | typo fixed in dir name |
| `input_models/MCC/`, `ScalingTests/`, `input/` | — | removed |
| — | `examples/PeriodicFreeSurface/` | new worked PBC example |
| `matlab/`, `utils/` | — | removed |
| `docs/`, `doc/Manual`, `doc/Paraview`, `doc/installation` | `doc/src/` + `doc/make.jl` | now a Documenter.jl site |
| — | `info/` | **new** quick reference: `info/options/input_file.dat`, `info/options/solver_options.txt` |
| `scripts/*.py, *.sh` (loose) | `scripts/bash/`, `scripts/julia/` | reorganized |
| `src/gravity.cpp/h` | — | feature removed (see §5.4) |
| — | `src/options.*`, `src/matFree.*`, `src/matData.*`, `src/matBFBT.cpp`, `src/matAux.cpp` | new solver machinery |

**Two files worth bookmarking in v3.0.0:** `info/options/input_file.dat` (every input
parameter) and `info/options/solver_options.txt` (every solver prefix/option, with
worked examples). They are the canonical reference this guide condenses.

---

## 4. Solver options: the headline change

### 4.1 The old way (v2.2.1)

Solver configuration was split between **top-level parameters** and a hand-written
**`<PetscOptionsStart>`** block of raw PETSc flags:

```text
# Solver options
    SolverType     = direct       # direct or multigrid
    DirectSolver   = mumps        # mumps / superlu_dist / pastix
    DirectPenalty  = 1e5

<PetscOptionsStart>
    -snes_type ksponly            # linear problem
    -js_ksp_type gmres
    -js_ksp_max_it 25
    -js_ksp_monitor
    -js_ksp_rtol 1e-4
    -js_ksp_atol 1e-10
<PetscOptionsEnd>
```

### 4.2 The new way (v3.0.0)

A structured `<SolverOptionsStart>` block. LaMEM expands it into the correct PETSc
options for you (`src/options.cpp`):

```text
<SolverOptionsStart>
    set_linear_problem = 1
    monitor_solvers    = 1
    linear_tolerances  = 1e-4 1e-10 25   # rtol, atol, maxit
    stokes_solver      = block_direct
    direct_solver_type = mumps
    penalty            = 1e5
<SolverOptionsEnd>
```

You can still drop to raw PETSc flags: a `<PetscOptionsStart>` block is honored, and
`skip_defaults = 1` disables the auto-presets entirely so you control everything.

### 4.3 Removed → new parameter map

| v2.2.1 parameter | v3.0.0 replacement |
|---|---|
| `SolverType = direct` | `stokes_solver = block_direct` (or `coupled_direct`) |
| `SolverType = multigrid` | `stokes_solver = coupled_mg` (or `block_mg`) |
| `DirectSolver = mumps/superlu_dist` | `direct_solver_type = mumps/superlu_dist/default` |
| `DirectPenalty = P` | `penalty = P` |
| `MGLevels = N` | `num_mg_levels = N` |
| `MGSmoother = jacobi` | `smoother_type = light` (richardson+jacobi) |
| `MGSmoother = chebyshev` | `smoother_type = intermediate` (chebyshev+sor) |
| `MGJacobiDamp = d` | `smoother_damping = d` |
| `MGSweeps = n` | `smoother_num_sweeps = n` |
| `MGCoarseSolver = redundant` | `coarse_solver = direct` + `coarse_num_cpu = 1` |
| `MGCoarseSolver = direct/mumps` | `coarse_solver = direct` + `direct_solver_type = …` |
| `MGRedundantNum` / `MGRedundantSolver` | *(removed — redundant coarse solver deprecated; telescope is used internally)* |
| `-snes_type ksponly` | `set_linear_problem = 1` |
| `-js_ksp_monitor` / `-snes_monitor` | `monitor_solvers = 1` |
| `-js_ksp_rtol / atol / max_it` | `linear_tolerances = rtol atol maxit` |
| `-snes_rtol / atol / max_it` | `nonlinear_tolerances = rtol atol maxit` |

> An `atol` of `-1` means "choose automatically". Tolerance triplets are
> `rtol atol maxit`; the block/coarse tolerances are pairs `rtol maxit`.

### 4.4 `stokes_solver` presets — what each expands to

Defined in `src/options.cpp:170-225`:

| `stokes_solver` | Expands to (PETSc) | Use when |
|---|---|---|
| `coupled_direct` | `-jp_type user` (+ penalty) | small problems; ⚠️ upstream recommends a `block_*` solver instead |
| `block_direct` | `-jp_type bf -bf_vs_type user` | direct velocity solve with penalty (Powell–Hestenes) |
| `coupled_mg` | `-jp_type mg` | coupled Galerkin geometric multigrid |
| `block_mg` | `-jp_type bf -bf_vs_type mg -bf_schur_type inv_eta` | block factorization, MG velocity, inverse-viscosity Schur |
| `wbfbt` | `-jp_type bf -bf_vs_type mg -bf_schur_type wbfbt` | **new** scaled BFBT Schur (Rudi et al. 2017) |

`smoother_type` (`src/options.cpp:602-614`): `light` = richardson + jacobi,
`intermediate` = chebyshev + sor, `heavy` = gmres + bjacobi (**default**).

> **`penalty` applies to _all_ direct solvers** _(post‑release)_ — both `coupled_direct`
> and `block_direct` (`src/options.h:41`). Earlier v3.0.0 docs said "only for
> block_direct"; this was corrected after the tag, so map your old `DirectPenalty` to
> `penalty` regardless of which direct solver you pick.
>
> **`direct_solver_type = default`** uses the PETSc built-in LU factorization, which is
> **sequential only** — fine on 1 core, but for parallel runs choose `mumps` or
> `superlu_dist`.

### 4.5 New capabilities with no v2.2.1 equivalent

- **Matrix-free Jacobian:** `use_mat_free_jac` in the block (or `-js_mat_free`).
- **Matrix-free MG levels:** `num_mat_free_levels` (or `-gmg_mat_free_levels N`) —
  N top levels run matrix-free, followed by one assembled level. *(Not supported with
  the block preconditioner or 2D coarsening.)*
- **wBFBT** Schur preconditioner: `stokes_solver = wbfbt`.
- **Zero-config solver:** omit the solver block entirely and `solverOptionsSetDefaults()`
  picks a working solver for your problem type.

### 4.6 Worked example — full before/after

**v2.2.1 — `FallingBlock_mono_PenaltyDirect.dat` (excerpt):**

```text
    SolverType     = direct
    DirectSolver   = mumps
    DirectPenalty  = 1e5
<PetscOptionsStart>
    -snes_type ksponly
    -js_ksp_type gmres
    -js_ksp_max_it 25
    -js_ksp_monitor
    -js_ksp_rtol 1e-4
    -js_ksp_atol 1e-10
<PetscOptionsEnd>
```

**v3.0.0 — `test/t01_FB1_Direct/FallingBlock_Direct_Default.dat` (the migrated form):**

```text
<SolverOptionsStart>
    set_linear_problem = 1
    monitor_solvers    = 1
    linear_tolerances  = 1e-4 1e-10 25   # rtol, atol, maxit
    stokes_solver      = block_direct
    direct_solver_type = default
    penalty            = 1e5             # keep your old DirectPenalty here!
<SolverOptionsEnd>
```

---

## 5. Input-file (`.dat`) breaking changes

### 5.1 Minimum 2 cells per direction (breaking)

FDSTAG now aborts if any direction has fewer than 2 cells (`src/fdstag.cpp:69`).
Classic quasi-2D setups that used `nel_y = 1` must use `nel_y = 2`.

**Observed error (verified by running):**

```
[0]PETSC ERROR: Less than two cells are specified in the y - direction
```

The run aborts immediately at grid setup (exit code 83).

### 5.2 Mesh bias ratios deprecated

Non-unit `bias_x/bias_y/bias_z` (graded meshes) are now rejected (`src/fdstag.cpp:83`):

```
Non-unit bias ratios are deprecated (bias_x)
```

Use uniform segments (or split a direction into multiple `nseg_*` uniform segments).

### 5.3 Background strain rates reduced

The out-of-plane background shear components were removed. In v2.2.1 `src/bc.cpp` read
`exz_strain_rates` and `eyz_strain_rates` (plus their `*_time_delims`); v3.0.0 keeps
only the in-plane / normal set:

- **Still supported:** `exx_strain_rates`, `eyy_strain_rates`, `exy_strain_rates`
  (+ matching `*_time_delims`).
- **Removed:** `exz_strain_rates`, `exz_time_delims`, `eyz_strain_rates`, `eyz_time_delims`.

If your model imposed background `exz`/`eyz`, that loading path no longer exists.

### 5.4 Gravity-survey feature removed

`src/gravity.cpp` and `src/gravity.h` were deleted with **no replacement** _(post‑release:
the files still existed at the `v3.0.0` tag and were removed in the following dev commits)_.
Any input relying on the gravity-survey / synthetic-gravity output is unsupported.

### 5.5 Periodic boundary conditions (new)

New optional `periodic` flag (`src/fdstag.cpp:914`). It is checked for compatibility at
startup (`src/bc.cpp:755-780`) and will abort with a clear message if combined with an
incompatible option. The guards (verbatim):

- `Periodic condition requires no-slip on bottom boundary`
- `Periodic condition is incompatible with no-slip on left and right boundaries`
- `Periodic condition is incompatible with xx and xy background strain rates`
- `Periodic condition is incompatible with boundary velocity (bvel_face)`
- `Periodic condition is incompatible with fixed phase (fix_phase)` / `fixed cells (fix_cell)`

See `examples/PeriodicFreeSurface/` for a working setup.

> **Not a breaking change**, but note: `act_temp_diff`, `act_steady_temp`,
> `Plume_Phase_Mantle`, and `Adjoint_FieldSensitivity` exist in *both* versions —
> they are **not** new in v3.0.0.

---

## 6. Tests & examples

### 6.1 Directory renaming

Test directories are now zero-padded and a few were renumbered:

| v2.2.1 | v3.0.0 |
|---|---|
| `t1_FB1_Direct` … `t9_PhaseDiagrams` | `t01_FB1_Direct` … `t09_PhaseDiagrams` |
| `t34_TopoDiffusion` | `t35_TopoDiffusion` |
| `t34_spatially_limited_erosion` | `t36_spatially_limited_erosion` |

### 6.2 Expected-output files renamed

The core-count / build suffix was dropped from `.expected` filenames:

```
Passive_tracer-2D_p1.expected        →  Passive_tracer-2D.expected
t14_1D_VP_Direct_opt-p1.expected     →  1D_VP_Direct_opt.expected
PhaseTransitions-FreeSlip_p1.expected →  PhaseTransitions-FreeSlip.expected
```

If you have downstream CI or scripts that reference the old `-pN` names, update them.

### 6.3 Running the suite

Unchanged entry point, but `test/start_tests.jl` now builds **both** `opt` and `deb`
before running, so `make test` compiles for you:

```bash
cd test
make test                                     # compiles opt+deb, runs full suite
julia --project=../. start_tests.jl           # equivalent
julia --project=../. start_tests.jl is64bit   # 64-bit integer build
julia --project=../. start_tests.jl use_dynamic_lib
```

The v2.2.1 `no_superlu` argument is gone; test-set names inside `runtests.jl` use the
new zero-padded directory names.

### 6.4 Passive-tracer output _(post‑release)_

Passive-tracer output now writes **one `.vtu` file per time step** (PR
[#70](https://github.com/UniMainzGeo/LaMEM/pull/70)) instead of the older single-series
layout. If you have post-processing that globs passive-tracer output, expect one file
per step now. Some direct-solver types in the test suite were also retuned
(`mumps`↔`superlu_dist`) because a few tests are sensitive to the factorization backend
— relevant if you run the suite on a non-x86 platform.

---

## 7. For people who patched the source

If you carry local C++ modifications, expect to rewrite them:

- **Error handling:** every `ierr = Foo(); CHKERRQ(ierr);` became `PetscCall(Foo());`.
  v3.0.0 contains **zero** `CHKERRQ` (v2.2.1 had ~4000). Port your patches to the
  `PetscCall(...)` form; MPI calls use `PetscCallMPI(...)`.
- **Integer types cleaned up** _(post‑release)_: `PetscMPIInt` was removed from LaMEM
  data structures and integer types/formatting were made consistent so the tree compiles
  cleanly in **64‑bit integer mode** (`start_tests.jl … is64bit`). If your patch stored
  ranks/counts as `PetscMPIInt` in a LaMEM struct, switch to `PetscInt` and cast only at
  the MPI call site.
- **C++17 + strict warnings:** the build adds `-Wconversion`, `-Wshorten-64-to-32`,
  `-Wnon-virtual-dtor`, `-Wformat=2`, and more. Implicit narrowing that compiled
  silently before may now warn (and CI treats the tree as warning-clean). Watch
  `PetscInt`↔`int`↔`size_t` conversions especially.
- **New solver files:** `options.cpp/h`, `matFree.cpp/h`, `matData.cpp/h`,
  `matBFBT.cpp`, `matAux.cpp`. If you touched residual/Jacobian assembly and use
  matrix-free operators, mirror changes into `matFree.cpp`.
- **Removed files:** `gravity.cpp/h`, `Makefile.in`. Remove references.
- **PETSc version guards:** use `PETSC_VERSION_LT(3, 24, 0)` / `(3, 25, 0)` if you
  need to support the full 3.19–3.25 range (see `src/options.cpp:65`).

---

## 8. Pitfalls & troubleshooting

There are **two failure modes** when you run a v2.2.1 file under v3.0.0, and one is
much more dangerous than the other.

### ⚠️ Pitfall #1 (SILENT — the important one): removed solver params are ignored

Running an **unmodified** v2.2.1 falling-block input under the v3.0.0 binary
**succeeds with exit code 0** — no error, no warning. But:

- `SolverType`, `DirectSolver`, `DirectPenalty`, and the `MG*` params are silently
  dropped (they are no longer parsed).
- LaMEM substitutes its **default** solver preset. In the tested case the file asked
  for `DirectPenalty = 1e5` but the run actually used the default
  `penalty = 1.0e3` (`-jp_pgamma 1000`).
- The raw flags inside `<PetscOptionsStart>` *are* still applied.

**Consequence:** your model appears to work, but may run with a different solver /
penalty / tolerances than you intended, changing convergence and possibly results.
**There is no diagnostic.** You must audit each `.dat` solver block by hand and port it
to `<SolverOptionsStart>` per [§4.3](@ref "4.3 Removed → new parameter map").

### Pitfall #2 (LOUD): grid & required-parameter errors

These abort immediately with a clear PETSc error — easy to spot and fix:

| Symptom | Cause | Fix |
|---|---|---|
| `Less than two cells are specified in the y - direction` | `nel_y = 1` (or any dir < 2) | set to `≥ 2` (§5.1) |
| `Non-unit bias ratios are deprecated (bias_x)` | graded mesh bias | use uniform segments (§5.2) |
| `act_temp_diff = 1 and -snes_type ksponly are incompatible` | old linear-solve flag + thermal diffusion | use `set_linear_problem = 1`, or `-snes_max_it 1` (`src/nlsolve.cpp:117`) |
| `Periodic condition is incompatible with …` | `periodic` + unsupported option | remove the conflicting option (§5.5) |
| Required-parameter-not-found on `exz_strain_rates` | removed background shear | drop it (§5.3) |

> **Migrating a file?** Work through the [Quick self-check](@ref "Quick self-check when migrating a file")
> at the top of this guide ([§0](@ref "0. TL;DR — the migration checklist")) — it turns every pitfall
> above into a concrete grep-and-fix step.

---

## Appendix: parameter change reference

**Removed in v3.0.0** (verified absent from `src/`):
`SolverType`, `DirectSolver`, `DirectPenalty`, `MGLevels`, `MGCoarseSolver`,
`MGSmoother`, `MGSweeps`, `MGJacobiDamp`, `MGRedundantNum`, `MGRedundantSolver`,
`exz_strain_rates`, `exz_time_delims`, `eyz_strain_rates`, `eyz_time_delims`.

**Added in v3.0.0** (verified present only in v3.0.0 `src/`):
`stokes_solver`, `direct_solver_type`, `penalty`, `coarse_solver`, `coarse_tolerances`,
`init_thermal_solver`, `skip_defaults`, `linear_tolerances`, `nonlinear_tolerances`,
`block_tolerances`, `thermal_tolerances`, `smoother_type`, `smoother_ksp`,
`smoother_pc`, `smoother_omega`, `smoother_damping`, `picard_to_newton`, `periodic`,
plus internal knobs (`num_mg_levels`, `num_mat_free_levels`, `coarse_num_cpu`,
`subdomain_*`, …) documented in `info/options/solver_options.txt`.

*Guide generated by comparing the local `LaMEM_v221` and `LaMEM_v300` trees (v3.0.0 at
HEAD `eaf75d4b`, 2026‑06‑15) and by executing v2.2.1 inputs under the v3.0.0
`bin/opt/LaMEM` binary. Post‑release deltas (`v3.0.0` tag → HEAD) were reconciled from
`git log v3.0.0..HEAD` and are flagged **_(post‑release)_** throughout.*
