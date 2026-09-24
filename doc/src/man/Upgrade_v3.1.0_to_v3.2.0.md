# Upgrading from v3.1.0 to v3.2.0

LaMEM v3.2.0 is a **feature release for surface processes**. It adds an optional coupling to the
[FastScape](https://fastscape.org) landscape-evolution library, a new slope-dependent erosion law,
and a `surf_mode` switch between the built-in surface processes and FastScape. Almost all of it is
opt-in. An unmodified v3.1.0 `.dat` file runs under v3.2.0 and gives the same numbers, with one
narrow exception: `topo_diff = 1` combined with `units = none`. The other breaking changes are in
the **Julia test harness**, and in **source patches that touch the time loop or the free surface**.

> **What this guide reflects.** The `v3.2.0` tag, which is also upstream `master` at commit
> `406b6444` (2026‑09‑22). No commits have landed after the tag. The baseline is `7e7a012e`
> (2026‑08‑12), the `v3.1.0` tag, which is the exact commit that
> *Upgrading from v3.0.0 to v3.1.0* was written against. The two guides therefore cover consecutive
> ranges with no gap and no overlap. The range spans 89 commits and 7 merged PRs: #76 (FastScape
> coupling), #80 (slope-dependent erosion), #81 (the previous upgrade guide), #84 (PETSc_jll 3.25.4
> in CI), #85 (README), #86 (single-threaded BLAS/OpenMP in tests) and #87 (version bump).
>
> Every claim below was verified by direct `grep` in both trees. The pitfalls in
> [§8](@ref "8. Pitfalls when upgrading to v3.2.0") were confirmed on three binaries built for this
> guide against PETSc 3.22.5 on aarch64 Linux: a v3.1.0 build, a default v3.2.0 build, and a v3.2.0
> build with `surf=scape`. Unmodified v3.1.0 input files were run on all three.

---

## 0. TL;DR — the v3.2.0 upgrade checklist

**If you run LaMEM from `.dat` files and never set `topo_diff` under `units = none`, your existing
inputs need no changes.** Items 1–4 are about input files; the rest concern building, testing and
patching.

1. **`topo_diff = 1` now requires `units = geo` or `units = si`.** Under `units = none`, v3.1.0 ran
   anyway. It used the value as a raw non-dimensional diffusivity while printing the unit label
   `[m^2/s]`. v3.2.0 stops with an error. The new `slope_dependent_erosion` has the same rule.
   **→ [§5.1](@ref "5.1 topo_diff now requires dimensional units")**
2. ⚠️ **Never set `surf_mode = 0`.** It is accepted, the run exits 0, and the free surface **never
   moves**: advection, erosion and sedimentation are all skipped. The only valid values are 1 (the
   default, the built-in processes) and 2 (FastScape).
   **→ [§4.1](@ref "4.1 surf_mode (v3.2.0)")**
3. ⚠️ **With `surf_mode = 2`, the built-in surface keys are silently ignored.** These are
   `erosion_model`, `sediment_model`, `slope_dependent_erosion`, `topo_diff` and all their
   parameters. They are not even checked for missing required companions.
   **→ [§4.3](@ref "4.3 FastScape coupling (surf_mode = 2)")**
4. ⚠️ **`vel_boundary` is documented backwards.** `FastScape.md` and `info/options/input_file.dat`
   say `0` = zero velocity. In the code, `1` sets the boundary velocity to zero and `0` keeps the
   velocity from LaMEM. **→ [§4.4](@ref "4.4 vel_boundary: the docs are inverted")**
5. **FastScape has to be enabled at compile time**, with `make mode=opt surf=scape all` and
   `FASTSCAPE_LIB` set. A default binary given `surf_mode = 2` stops at startup. Check any binary
   with `LaMEM -fastscape_info`. **→ [§1.2](@ref "1.2 Optional FastScape build")**
6. **Test harness: the `opt=` keyword is gone** from `perform_lamem_test`, `run_lamem_local_test`
   and `CreatePartitioningFile_local`, so passing it raises a `MethodError`. **`deb=true` now really
   runs the debug binary**; in v3.1.0 it silently ran the optimized one.
   **→ [§6.1](@ref "6.1 The opt keyword was removed")**
7. **`make test` in `test/` builds with FastScape whenever `FASTSCAPE_LIB` is exported** in your
   shell. **→ [§6.3](@ref "6.3 FASTSCAPE_LIB switches the test build")**
8. **`make format` and `make check` in `src/` now require exactly astyle 3.1** and refuse any other
   version. **→ [§1.3](@ref "1.3 astyle is pinned to 3.1")**
9. **PETSc 3.25.4 is now the recommended version** and the CI version. The code still compiles
   against 3.22–3.24. **→ [§1.1](@ref "1.1 PETSc and toolchain (v3.2.0)")**
10. **If you patched the time loop or the free surface**, the surface-process calls moved inside
    `if(lm->surf.SurfMode == 1)`, and `FreeSurf` and `LaMEMLib` gained members.
    **→ [§7](@ref "7. If you patched the v3.1.0 source")**

### Quick self-check for a v3.1.0 input file

Work through these checks from top to bottom on one `.dat` file:

```bash
grep -nE "^\s*units\s*=\s*none"  my_model.dat   # together with…
grep -nE "^\s*topo_diff\s*=\s*1" my_model.dat   # → both hit? now a hard error (§5.1)
grep -n  "surf_mode"             my_model.dat   # → "= 0"? free surface frozen (§4.1)
grep -rn "opt *= *\(true\|false\)" my_tests/    # → keyword removed from test_utils.jl (§6.1)
grep -rn "deb *= *true"          my_tests/      # → now really runs bin/deb (§6.1)
echo "$FASTSCAPE_LIB"                           # → non-empty? make test builds surf=scape (§6.3)
astyle --version                                # → not 3.1? make format / make check refuse (§1.3)
```

Nothing was removed from the `.dat` vocabulary. LaMEM v3.1.0 could parse 471 parameters, and all
471 still parse in v3.2.0. There are 46 new parameters, all opt-in; see the
[Appendix](@ref "Appendix: v3.2.0 parameter change reference").

---

## 1. Requirements & build changes in v3.2.0

### 1.1 PETSc and toolchain (v3.2.0)

**The code's PETSc version guards did not change.** The four `PETSC_VERSION_LT` checks are the same
in both trees (`src/LaMEMLib.cpp:1003`, `src/JacResTemp.cpp:232`, `src/options.cpp:70`,
`src/lsolve.cpp:358`). What changed is the recommendation:

| | v3.1.0 | v3.2.0 |
|---|---|---|
| Recommended PETSc (`Installation.md`) | 3.22.5 | **3.25.4** ("also still compiles against 3.22.x – 3.24.x") |
| CI `PETSc_jll` (`test/setup_packages.jl:5`) | 3.22.0 | 3.25.4 |
| CI `MPICH_jll` | 4.2.3 | 5.0.1 |
| CI Julia (`.github/workflows/CI.yml:43`) | 1.10 | **1.12.7** (exact patch: 1.12.6 ships SuiteSparse_jll 7.8.0 and does not resolve) |

The claim that older versions still work was checked for this guide: v3.2.0 builds warning-free
against PETSc 3.22.5, in both the default and the `surf=scape` configurations.

`Installation.md` has a new **§1.2 "Compiling LaMEM against the precompiled PETSc"**, which is the
same route CI uses. In brief:

- Deploy `PETSc_jll` into `/workspace/destdir`. The path is hard-coded in the packaged
  `petscvariables`.
- Point `PETSC_OPT` at `lib/petsc/double_real_Int32` (or `…Int64`).
- Build as usual.
- From PETSc_jll 3.25 on, a binary built this way has to be told where BLAS is before it runs
  outside Julia: `LBT_DEFAULT_LIBS` must list **both** the ILP64 and LP64 OpenBLAS.

*(This route was not reproduced on the machine used for this guide. It is summarised from
`Installation.md` and `test/test_utils.jl:621`.)*

### 1.2 Optional FastScape build

FastScape is **not** part of the default build. `src/Makefile` gained a `surf` variable, which
defaults to `none`:

```bash
export FASTSCAPE_LIB=/dir/containing/libfastscapelib_fortran   # new
cd src && make mode=opt surf=scape all                          # adds -DWITH_FASTSCAPE
./bin/opt/LaMEM -fastscape_info                                 # prints FASTSCAPE_ENABLED / FASTSCAPE_DISABLED
```

- `surf=scape` without `FASTSCAPE_LIB` stops immediately:
  `Makefile:46: *** Environmental variable FASTSCAPE_LIB must be set to installation directory.  Stop.`
- The link line embeds an rpath to `FASTSCAPE_LIB` (`src/Makefile:116–120`), so the binary finds
  the library at run time without `LD_LIBRARY_PATH`.
- **Switching an existing build to `surf=scape` does not recompile it.** The object files do not
  depend on the `surf` flag, so on a tree already built without FastScape,
  `make mode=opt surf=scape all` prints `Nothing to be done for 'all'`, and the binary stays
  `FASTSCAPE_DISABLED`. Run `make mode=opt clean_all` (and `mode=deb`) first, then rebuild.
  Switching back the other way needs the same step.
- `-fastscape_info` (`src/LaMEM.cpp:30–43`) prints its answer and exits before reading any input
  file. The test suite uses it to decide whether to skip FastScape tests.
- Upstream tests against Fastscapelib 2.8.4. The easiest source is
  `Pkg.add("Fastscapelib_jll")`, then `Fastscapelib_jll.LIBPATH_list`. For this guide, the
  `surf=scape` binary was linked against the `Fastscapelib_jll` 2.8.4 shared library.

### 1.3 astyle is pinned to 3.1

The new `checkastyle` target (`src/Makefile:169`) runs before `make format`, `make checkformat`
and therefore `make check`. It refuses any astyle other than **3.1**, the version shipped with
Ubuntu 22.04/24.04 and used by CI. The reason is that newer astyle releases reformat files that 3.1
considers correct. With Homebrew's astyle, for example, you get:

```
ERROR: wrong astyle version.
  required : 3.1
  found    : Artistic Style Version 3.6.9
```

Install 3.1 (from distro packages or built from source) if you contribute code. Nothing changes if
you only compile and run LaMEM.

---

## 2. What's new in v3.2.0 (the "why")

- **FastScape coupling (#76).** On each time step LaMEM passes the surface velocities to FastScape.
  FastScape advances the topography, in sub-steps of at most `max_fs_dt`, solving for:
  - fluvial incision (stream-power law);
  - hillslope diffusion;
  - sediment transport and deposition;
  - optionally, marine transport with compaction.

  LaMEM then takes the topography back and converts markers between air, rock and the sediment
  phase. The whole coupling (`src/fastscape.cpp`, about 3 100 lines) is compiled only with
  `-DWITH_FASTSCAPE`. It works in 2D and 3D, on uniform and non-uniform grids. The user
  documentation is the new page [Surface processes with FastScape](FastScape.md).
- **Slope-dependent erosion (#80).** `E [m/yr] = prefactor_slope · |∇h|^n_slope` is applied on
  the internal free surface. It sits on top of any `erosion_model`, after erosion and before
  sedimentation and topographic diffusion. It can mimic stream-power incision in profile models,
  where full FastScape coupling would be excessive.
- **Unit-safety checks.** Both `topo_diffusivity` (m²/s) and `prefactor_slope` (m/yr) are now
  rejected under `units = none`, because there they would be silently reinterpreted as raw
  non-dimensional numbers.
- **CI modernisation (#84, #86).** PETSc_jll 3.25.4 on Julia 1.12.7. The harness now pins BLAS and
  OpenMP to one thread. This fixes a SuperLU_DIST/OpenBLAS race that made tests fail
  intermittently, and sometimes return wrong residuals, about two runs in three.

---

## 3. Repository layout changes in v3.2.0

Added:

| Path | What |
|------|------|
| `src/fastscape.cpp`, `src/fastscape.h` | FastScape coupling (only active with `surf=scape`) |
| `doc/src/man/FastScape.md` | User guide: building, input block, output, troubleshooting |
| `doc/src/man/Upgrade_v3.0.0_to_v3.1.0.md` | The previous upgrade guide, now in the docs under Release Notes |
| `test/t37_Collision_FastScape/` | FastScape-coupled collision test (skipped unless the binary has FastScape) |
| `test/t38_slope_dependent_erosion/` | Slope-dependent erosion on a dome, with analytic validation in the PR |

Renamed: `test/t02_FB2_MG/FB2_a_CoupledMG_opt.expected` → `FB2_a_CoupledMG_deb.expected`. The name
now matches what the test actually runs; see [§6.1](@ref "6.1 The opt keyword was removed").

No directories were moved or removed. No files were deleted.

---

## 4. New surface-process options

All of these belong to the free-surface section and are read only when `surf_use = 1`.

### 4.1 surf_mode (v3.2.0)

```
    surf_mode = 1     # 1 - built-in erosion/sedimentation (default), 2 - FastScape
```

Read at `src/surf.cpp:55`, default 1 (`src/surf.cpp:38`), upper bound 2. The time loop
(`src/LaMEMLib.cpp:700`, `:723`) has branches only for `== 1` and `== 2`.

!!! warning "surf_mode = 0 freezes the free surface — silently"
    `src/surf.h:44` describes 0 as "none", and the parser accepts it. **But neither branch of the
    time loop runs**, so `FreeSurfAdvect` is skipped along with every surface process. Observed
    here: `test/t15_RTI/RTI.dat` with `surf_mode = 0` added exited 0. Its surface topography after
    2 steps was still exactly the initial 0.499–0.501, whereas the unmodified run had already
    evolved to 0.49891–0.50109. The free-surface summary also stops printing the erosion,
    sedimentation and diffusion lines. To switch surface processes off, use `surf_mode = 1` with
    `erosion_model = 0` and `sediment_model = 0` (the defaults), not `surf_mode = 0`.

`surf_mode = 3` fails loudly:
`Entry 1 in parameter "[-]surf_mode" is larger than allowed : val=3, max=2`.

### 4.2 Slope-dependent erosion

```
    slope_dependent_erosion = 1     # default 0
    prefactor_slope         = 0.05  # [m/yr] — NOTE m/yr, unlike er_rates/sed_rates in cm/yr (default 1.0)
    n_slope                 = 2.0   # [-] (default 1.0)
```

- The parameters are read at `src/surf.cpp:123–141` and the law is implemented in
  `FreeSurfAppSlopeErosion` (`src/surf.cpp:1103`).
- The slope is a centred difference, one-sided at the global boundaries.
- The update is explicit, with automatic sub-stepping: the drop per sub-step is capped at half the
  minimum horizontal spacing (`src/surf.cpp:1208`). This keeps it well defined for any `n_slope`,
  including `n_slope < 1`.
- Each step prints `Applying slope-dependent erosion (E = … [m/yr] * slope^…) in N sub-step(s).`
- It works only with `surf_mode = 1`. It requires `units = geo` or `si`
  ([§5.1](@ref "5.1 topo_diff now requires dimensional units")).

### 4.3 FastScape coupling (surf_mode = 2)

The requirements are checked in this order. Each one fails loudly if it is not met:

| Requirement | Error if violated |
|---|---|
| Binary built with `surf=scape` | `surf_mode = 2 (FastScape) requested, but this LaMEM binary was built without FastScape support. Recompile with WITH_FASTSCAPE enabled, or select a different surf_mode.` (`src/surf.cpp:67`) |
| A `<FastScapeStart>` … `<FastScapeEnd>` block | `<FastScapeStart> - <FastScapeEnd> blocks must be defined` (`src/fastscape.cpp:74`) |
| `units = geo` or `units = si` | `Incorrect unit type for FastScape` (`src/fastscape.cpp:339`) |
| Required block keys (`max_fs_dt`, `topo_boundary`, `vel_boundary`, `random_noise`, `sed_phases`, `kf`, `kfsed`, `m`, `n`, `kd`, `kdsed`, `g`, `gsed`, `p`; plus `extendedRange` and `extendedNodes` if `fs2D = 1`; plus the nine marine keys if `set_marine = 1`) | the usual missing-parameter error |

A minimal block that ran successfully here on a 2D geo-unit model (`fs2D = 1`):

```
<FastScapeStart>
    fs2D          = 1        # 2D LaMEM model → extend the profile in y
    extendedRange = 100      # [km]
    extendedNodes = 33
    max_fs_dt     = 0.001    # [Myr]
    topo_boundary = 1111
    vel_boundary  = 0000     # see §4.4 — 0 KEEPS the LaMEM velocity
    random_noise  = 1
    sed_phases    = 1        # phase ID given to new sediment (single integer here)
    kf = 1e-6   kfsed = -1.0   m = 0.4   n = 1.0
    kd = 1e-2   kdsed = -1.0   g = 0.0   gsed = 0.0   p = -2.0
<FastScapeEnd>
```

(Written one key per line in a real file.)

!!! warning "Built-in surface keys are ignored under surf_mode = 2 — silently"
    The whole built-in block (`src/surf.cpp:77–157`) is skipped. `erosion_model`,
    `sediment_model`, `slope_dependent_erosion`, `topo_diff` and their companions (`er_*`, `sed_*`,
    `topo_diffusivity`, …) are never read. Observed here: `erosion_model = 2`, `sediment_model = 1`,
    `topo_diff = 1` and `slope_dependent_erosion = 1` were added **without** any of their required
    parameters. The run exited 0, FastScape ran on both steps, and no message mentioned any of
    them. If you switch a model to FastScape, remove those lines, so that nobody reading the file
    later assumes they are in effect.

!!! note "Reused key names inside the block"
    `n` and `sed_phases` already exist elsewhere in LaMEM: `n` is the power-law exponent in
    material blocks, and `sed_phases` is the **array** of sediment phases for `surf_mode = 1`.
    Inside `<FastScapeStart>`, `n` is the stream-power slope exponent and `sed_phases` is a
    **single** phase ID. Do not copy a `sed_phases = 2 3 4` line into the FastScape block.

Other constraints. [Surface processes with FastScape](FastScape.md) was written for the PR and may
lag the code, so each point below says whether it was checked against the source:

- FastScape runs serially on rank 0. **Checked**: the Fortran solve is inside
  `if(ISRankZero(...))` (`src/fastscape.cpp:1542`).
- Its surface grid is built once and never regenerated, so models whose LaMEM grid stretches (a
  background strain rate) are **not supported**. **Documented only**: nothing in `fastscape.cpp`
  detects a background strain rate, so such a model is not rejected. Treat the limitation as real,
  but do not expect an error message to warn you.
- With `surf_mode = 2`, FastScape state is appended to restart files (`src/surf.cpp:314`,
  `:335`). Restart files from `surf_mode = 1` runs are unchanged. **Checked** in the source; a
  restart was not run for this guide.
- The page states `max_fs_dt` in `[Myr]`. **That is only true for `units = geo`.** The value is
  read in LaMEM time units (`src/fastscape.cpp:102`, `:306`, `:324`), so under `units = si` it is
  in seconds.

### 4.4 vel_boundary: the docs are inverted

`vel_boundary` is a four-digit string with one digit per boundary, in the order bottom, right, top,
left. The implementation (`src/fastscape.cpp:2826–2851`) **zeroes** the velocity on a boundary
whose digit is `'1'` and leaves it untouched otherwise:

```c
if (FSLib->FS_VELBC[0] == '1' && j == 0) { vx_pass[ind] = 0; vy_pass[ind] = 0; vz_pass[ind] = 0; }
```

The comment at `src/fastscape.cpp:109` agrees ("1 -- boundary velocity == 0; 0 -- boundary velocity
from LaMEM"). The shipped documentation, including the published
[Surface processes with FastScape](FastScape.md) page, says the opposite:

- `doc/src/man/FastScape.md:102`: "`0` sets the boundary velocity to zero and `1` keeps the velocity
  transferred from LaMEM"
- `info/options/input_file.dat:206`: "0-zero velocity, 1-keep velocity transferred from LaMEM"

A correction of both files is proposed in UniMainzGeo/LaMEM#88. Until it is merged, **go by the code**: `0000` (as used in t37) keeps the LaMEM velocity on all four sides, and `1111`
pins all four to zero. `topo_boundary` is documented correctly (`0` reflective, `1` fixed height).

### 4.5 FastScape output

A coupled run also writes `<out_file_name>_fs.pvd`. This is **on by default**: `out_surf_fs`,
`out_fs_pvd` and every field flag default to 1 (`src/fastscape.cpp:800–814`), so you set a flag to 0
to drop that field. The fields are selected
with the 12 other `out_surf_*` / `out_fs_pvd` flags listed in the
[Appendix](@ref "Appendix: v3.2.0 parameter change reference"). The write frequency is set with
`surf_out_nstep` in the FastScape block.

---

## 5. Input-file breaking changes in v3.2.0

There is one, and it is loud.

### 5.1 topo_diff now requires dimensional units

`src/surf.cpp:148–155` (commit `cabca8c7`). `topo_diffusivity` is specified in m²/s and scaled by
`length_si² / time_si`. Under `units = none` both of those are 1, so v3.1.0 used the value as a
raw non-dimensional diffusivity while printing it as `[m^2/s]`.

**Observed**: `test/t15_RTI/RTI.dat` (`units = none`) with `topo_diff = 1` and
`topo_diffusivity = 1e-5` appended. Under v3.1.0 it exits 0 and prints
`Topographic diffusion : active (K = 1e-05 [m^2/s])`. Under v3.2.0 it stops with:

```
topo_diff defines topo_diffusivity in [m^2/s] and requires a dimensional unit system. Set units = geo or units = si, or deactivate topo_diff.
```

**Fix**: switch the model to `units = geo`/`si`. If you really meant a non-dimensional
diffusivity, drop `topo_diff`: there is no longer a way to express that.

No `units = none` input shipped with v3.1.0 uses `topo_diff`. The only in-tree `topo_diff` test,
t35, uses `units = geo`, so no example or test needed changing.

---

## 6. Test framework & CI in v3.2.0

### 6.1 The opt keyword was removed

`opt` was dropped from `run_lamem_local_test` (`test/test_utils.jl:49`), `perform_lamem_test`
(`:685`) and `CreatePartitioningFile_local`. All 81 `opt=true` occurrences were removed from
`runtests.jl`. The binary is now chosen by `deb` alone: `deb=false` (the default) runs `bin/opt`,
and `deb=true` runs `bin/deb`.

**Before (v3.1.0)**:

```julia
perform_lamem_test(dir, ParamFile, "MyTest_opt"; opt=true, keywords=keywords, ...)
perform_lamem_test(dir, ParamFile, "MyTest_deb"; deb=true, opt=false, ...)
```

**After (v3.2.0)**:

```julia
perform_lamem_test(dir, ParamFile, "MyTest_opt"; keywords=keywords, ...)
perform_lamem_test(dir, ParamFile, "MyTest_deb"; deb=true, ...)
```

Passing `opt` now fails when the function is called:

```
MethodError: no method matching perform_lamem_test(::String, ::String, ::String; opt::Bool)
This method does not support all of the given keyword arguments (and may not support any).
  … got unsupported keyword argument "opt"
```

!!! warning "deb=true used to run the optimized binary"
    In v3.1.0, `opt` defaulted to `true` and was checked **before** `deb` (v3.1.0
    `test_utils.jl:57`). So `deb=true` without `opt=false` silently ran `bin/opt`. Two upstream
    tests did exactly that: `runtests.jl:145` and `:1175`, the t24 test. In v3.2.0 they genuinely
    run `bin/deb`. If your own tests pass `deb=true`, they now run a different binary, and an
    expected file generated from the optimized build may need regenerating (`make update <N>`).
    This is also why `FB2_a_CoupledMG_opt.expected` was renamed `_deb`: that test had passed
    `opt=false` and always ran the debug build.

### 6.2 FastScape-only tests

`LaMEM_has_fastscape(; bin_dir, deb)` (`test/test_utils.jl:165`, exported) runs the binary with
`-fastscape_info`. t37 (`runtests.jl:1492`) uses it to gate the opt and deb runs separately. When
the binary lacks FastScape, the testset records `@test_skip "FastScape is not installed"` instead
of failing. Use the same pattern for any test that needs an optional build feature.

### 6.3 FASTSCAPE_LIB switches the test build

`test/start_tests.jl:29` rebuilds LaMEM with `surf=scape`, both opt and deb, **whenever
`FASTSCAPE_LIB` is set in the environment**, and with the plain build otherwise. So an exported
`FASTSCAPE_LIB` in your `.bashrc` silently turns every `make test` into a FastScape build. If it
points at the wrong directory, the link fails. To test the default build, use:

```bash
env -u FASTSCAPE_LIB make test 38
```

The switch does **not** force a rebuild. If `bin/` already holds a non-FastScape build, `make`
reports `Nothing to be done` ([§1.2](@ref "1.2 Optional FastScape build")), t37 skips itself, and
the suite still reports success, with t37 counted as `Broken`. This was observed while preparing
this guide. Run `make mode=opt clean_all` and `make mode=deb clean_all` in `src/` before switching,
and confirm with `../bin/opt/LaMEM -fastscape_info`.

CI (`test/compile_lamem.jl:41`) always builds with `surf=scape`, against `Fastscapelib_jll`
staged in `/workspace/destdir/lib/fastscape`.

### 6.4 Runtime environment for test runs

Every LaMEM process started by the harness now goes through `add_dylibs`
(`test/test_utils.jl:621`), which does three things:

- **Pins `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` to 1**, even for
  a locally built PETSc. SuperLU_DIST's OpenMP threads were intermittently passing corrupted
  arguments into OpenBLAS (`** On entry to DGEMM parameter number 8 had an illegal value`). The
  corrupted solves produced wrong residuals, for example `|mRes|_2` off by 82 % in t02, according
  to PR #86. If you run a SuperLU_DIST model by hand with threaded BLAS and see that message, set
  `OMP_NUM_THREADS=1`.
- Sets the platform's real loader variable (`LD_LIBRARY_PATH` on Linux,
  `DYLD_FALLBACK_LIBRARY_PATH` on macOS, `PATH` on Windows). v3.1.0 always set the macOS variable.
- Sets `LBT_DEFAULT_LIBS` for PETSc_jll ≥ 3.25 ([§1.1](@ref "1.1 PETSc and toolchain (v3.2.0)")).

`t23_Permeable` now uses `direct_solver_type = mumps` (`test/t23_Permeable/Permeable.dat:267`),
because the `superlu_dist` solve diverged under PETSc 3.25 with Int64. The description of PR #86
says SuperLU_DIST was restored there, but the merged tree still uses MUMPS, and the `add_dylibs`
docstring confirms it. `coupled_direct` + SuperLU_DIST is therefore not covered by the v3.2.0
suite.

### 6.5 New tests and numbering

t37 (`t37_Collision_FastScape`) and t38 (`t38_slope_dependent_erosion`) were added. PR #80 used
the number t37 but was renumbered to t38 at merge. **The next free number is t39.**

For this guide the complete v3.2.0 suite was run locally (`make test`, PETSc 3.22.5, 32-bit
indices, aarch64 Linux): **103 passed, 0 failed**, in 23 min 49 s. The only non-passing entries
were t37's two `@test_skip`s: the binary in `bin/` turned out to be a default build (see
[§6.3](@ref "6.3 FASTSCAPE_LIB switches the test build")), so t37 is covered here only by CI's
`surf=scape` build.

---

## 7. If you patched the v3.1.0 source

Much smaller than the v3.1.0 upgrade. Only one prototype was added to the existing headers
(`FreeSurfAppSlopeErosion`), and none were changed or removed. The edits that can conflict with a
patch are these:

**The time loop.** In `LaMEMLibSolve` (`src/LaMEMLib.cpp:699–755`), the free-surface steps are
now selected by mode.

Before:

```c
PetscCall(FreeSurfAdvect(&lm->surf));
...
PetscCall(FreeSurfAppErosion(&lm->surf));
PetscCall(FreeSurfAppSedimentation(&lm->surf));
PetscCall(FreeSurfAppTopoDiffusion(&lm->surf));
```

After:

```c
if(lm->surf.SurfMode == 1) PetscCall(FreeSurfAdvect(&lm->surf));
#ifdef WITH_FASTSCAPE
if(lm->surf.SurfMode == 2) PetscCall(FastScapeCopyVelocity(&lm->FSLib));
#endif
...
if(lm->surf.SurfMode == 1)
{
    PetscCall(FreeSurfAppErosion(&lm->surf));
    PetscCall(FreeSurfAppSlopeErosion(&lm->surf));   // new
    PetscCall(FreeSurfAppSedimentation(&lm->surf));
    PetscCall(FreeSurfAppTopoDiffusion(&lm->surf));
}
#ifdef WITH_FASTSCAPE
if(lm->surf.SurfMode == 2) { FastScapeRun; FreeSurfSmoothMaxAngle; FreeSurfGetAvgTopo; }
#endif
```

A custom surface process belongs inside the `SurfMode == 1` block. Its parameters must be read
inside the `if(surf->SurfMode == 1)` block of `FreeSurfCreate` (`src/surf.cpp:77`).

**Structs.**
- `FreeSurf` (`src/surf.h`) gained `FastScapeLib *FSLib`, `SurfMode`, and
  `slope_dependent_erosion` / `prefactor_slope` / `n_slope`.
- `LaMEMLib` gained `FastScapeLib FSLib` (`src/LaMEMLib.h:45`).
- These members are **unconditional**: they exist in default builds too. Only the code that uses
  FastScape is behind `#ifdef WITH_FASTSCAPE`.
- `LaMEMLibSetLinks` wires `surf.FSLib` and `FSLib.{surf,pvsurf,jr,scal}`.

**Marker phase conversion.** In `ADVMarkCrossFreeSurf` (`src/subgrid.cpp:545`) and its
passive-tracer counterpart (`src/passive_tracer.cpp:885`), an air marker that ends up below the
surface under `surf_mode = 2` always becomes `surf->phase`, the FastScape `sed_phases`. The v3.1.0
`SedimentModel` logic now runs only in mode 1.

**Unit conversions.** If you add a parameter with physical units, look at the new `units = none`
guards (`src/surf.cpp:129`, `:150`) before relying on `getScalarParam`'s scaling argument alone.

**Formatting.** `make check` needs astyle 3.1 ([§1.3](@ref "1.3 astyle is pinned to 3.1")).

To see exactly what moved:

```bash
git diff 7e7a012e 406b6444 -- src/surf.h src/surf.cpp src/LaMEMLib.h src/LaMEMLib.cpp src/subgrid.cpp src/passive_tracer.cpp
```

---

## 8. Pitfalls when upgrading to v3.2.0

### Silent: surf_mode = 0 freezes the free surface

The run exits 0, the topography never changes, and no message says so. See
[§4.1](@ref "4.1 surf_mode (v3.2.0)").

### Silent: built-in surface keys under surf_mode = 2

`erosion_model`, `sediment_model`, `slope_dependent_erosion`, `topo_diff` and their parameters are
ignored without even a missing-parameter check. See
[§4.3](@ref "4.3 FastScape coupling (surf_mode = 2)").

### Silent: vel_boundary means the opposite of what the docs say

If you follow the docs, `1111` looks like "keep LaMEM velocities everywhere", but it actually pins
all four boundaries to zero. See [§4.4](@ref "4.4 vel_boundary: the docs are inverted").

### Silent: deb=true in your own tests now runs a different binary

A test that passed in v3.1.0 was comparing the **opt** binary's output. It may now drift outside
tolerance against an expected file generated by opt. See
[§6.1](@ref "6.1 The opt keyword was removed").

### Silent: an exported FASTSCAPE_LIB changes what make test builds

See [§6.3](@ref "6.3 FASTSCAPE_LIB switches the test build").

### Silent: a stale build ignores surf=scape

`make … surf=scape all` on an existing default build does nothing, so the binary has no FastScape
and FastScape tests are skipped rather than failed. Run `clean_all` first
([§1.2](@ref "1.2 Optional FastScape build")).

### Loud: topo_diff with units = none

```
topo_diff defines topo_diffusivity in [m^2/s] and requires a dimensional unit system. Set units = geo or units = si, or deactivate topo_diff.
```

The same applies to `slope_dependent_erosion`
(`slope_dependent_erosion defines prefactor_slope in [m/yr] and requires a dimensional unit system. …`).

### Loud: FastScape requested but not available or not configured

```
surf_mode = 2 (FastScape) requested, but this LaMEM binary was built without FastScape support. …
<FastScapeStart> - <FastScapeEnd> blocks must be defined
Incorrect unit type for FastScape
Makefile:46: *** Environmental variable FASTSCAPE_LIB must be set to installation directory.  Stop.
```

### Loud: opt= in a test

A `MethodError` naming the unsupported keyword `opt`. Delete the keyword
([§6.1](@ref "6.1 The opt keyword was removed")).

### Loud: wrong astyle

`ERROR: wrong astyle version.` Install astyle 3.1.

### Not a pitfall: your v3.1.0 input files

Seven unmodified v3.1.0 free-surface inputs were run for 2 steps under both a v3.1.0 and a v3.2.0
binary: t10 `Compressible1D_withSaltandBasement`, t14 `1D_VP`, t15 `RTI`, t16 `TimeTransition`,
t20 `RTI_FSSA`, t22 `ridge_geom_2D` and t32 `BC_velocity_2D_LR`. Every `|Div|_inf`, `|Div|_2` and
`|mRes|_2` line was **bit-identical**, so the `SurfMode` refactor does not change `surf_mode = 1`
results. The expected files of t24 (erosion + sedimentation), t35 (topographic diffusion) and t36
(spatially limited erosion) did not change between the tags either.

---

## Appendix: v3.2.0 parameter change reference

Extracted with the `get*Param` vocabulary diff. **Each of the 46 additions was verified by direct
`grep` to be absent from v3.1.0 `src/`.** Line numbers refer to v3.2.0. No parameter was removed.

**Free surface (`surf_mode = 1` or general):**

| Parameter | Location | Default | Notes |
|-----------|----------|---------|-------|
| `surf_mode` | `src/surf.cpp:55` | 1 | 1 built-in, 2 FastScape; 0 accepted but freezes the surface ⚠️ |
| `slope_dependent_erosion` | `src/surf.cpp:123` | 0 | mode 1 only; needs dimensional units |
| `prefactor_slope` | `src/surf.cpp:140` | 1.0 | **m/yr** |
| `n_slope` | `src/surf.cpp:141` | 1.0 | dimensionless |

**`<FastScapeStart>` block (`src/fastscape.cpp:82–156`), `surf_mode = 2` only:**

| Group | Parameters | Required? |
|-------|-----------|-----------|
| Grid | `non_uniform_grid`, `fs2D`, `fs_refine` | optional (defaults 0, 0, 1) |
| 2D extension | `extendedRange` [km], `extendedNodes` (> 2; a value of 2 is rejected with the message "must be ≥ 2") | required if `fs2D = 1` (`extendedNodes` only on a uniform grid) |
| Coupling & BCs | `max_fs_dt` [LaMEM time unit: Myr for geo, s for si], `topo_boundary`, `vel_boundary` (see §4.4), `random_noise`, `sed_phases`* | required |
| Stream power & hillslope | `kf`, `kfsed`, `m`, `n`*, `kd`, `kdsed`, `g`, `gsed`, `p` | required |
| Marine | `set_marine` (default 0); `sealevel`, `poroSilt`, `poroSand`, `zporoSilt`, `zporoSand`, `ratio`, `depth_siltsand_solve`, `kdsSilt`, `kdsSand` | the nine are required if `set_marine = 1` |
| Output cadence | `surf_out_nstep`, `vec_times` | optional (default 1) |

\* `n` and `sed_phases` are not new names (see [§4.3](@ref "4.3 FastScape coupling (surf_mode = 2)")),
so the bulk diff does not list them. Inside the block they are new keys with different meanings.
That brings the count to 42 FastScape keys in the diff (29 in the block plus 13 output flags), and
2 reused ones.

**FastScape output flags (`src/fastscape.cpp:818–838`), all optional:**
`out_surf_fs`, `out_fs_pvd`, `out_surf_topofs`, `out_surf_silt_fraction`, `out_surf_basement`,
`out_surf_total_erosion`, `out_surf_drainage_area`, `out_surf_erosion_rate`, `out_surf_slope`,
`out_surf_curvature`, `out_surf_chi`, `out_surf_catchment`, `out_surf_lake_depth`.

Totals: v3.1.0 parses 471 parameters and v3.2.0 parses 517. That is 471 unchanged, 0 removed and
46 added: 4 free-surface, 29 FastScape block, 13 FastScape output.

Other user-visible changes that are not `.dat` parameters:

| Change | Location | Kind |
|--------|----------|------|
| `-fastscape_info` command-line probe | `src/LaMEM.cpp:30` | Added |
| `make … surf=scape`, `FASTSCAPE_LIB` | `src/Makefile:23, 44–49, 97–100` | Added |
| `checkastyle` (astyle 3.1 pin) | `src/Makefile:169` | Added |
| Error for `topo_diff` under `units = none` | `src/surf.cpp:150` | Behaviour change |
| `opt` keyword in `test_utils.jl` | `test/test_utils.jl` | Removed |
| `LaMEM_has_fastscape`, `add_dylibs` | `test/test_utils.jl:165`, `:621` | Added |

---

*Guide generated by diffing upstream `7e7a012e` (v3.1.0, 2026‑08‑12) against `406b6444`
(v3.2.0, 2026‑09‑22). It was checked by running unmodified v3.1.0 inputs and targeted variants
under freshly built v3.1.0, v3.2.0 and v3.2.0 + FastScape binaries (PETSc 3.22.5, aarch64 Linux).
Parameter claims were confirmed by direct `grep` in both trees rather than taken from the bulk
vocabulary diff.*
