# Upgrading from v3.0.0 to v3.1.0

LaMEM v3.1.0 is a **maintenance and refactoring release**. Unlike the v2.2.1 → v3.0.0 jump, it
changes almost nothing about how you write a `.dat` file: existing v3.0.0 input files run
unmodified. The disruption is concentrated in two places — **code that patches LaMEM's source**,
and **test or CI harnesses** built on the Julia test framework.

> **What this guide reflects.** Upstream `master` at commit `7e7a012e` (2026‑08‑12), where the
> `Version : 3.1.0` string first appears on `master`. The baseline is `eaf75d4b` (2026‑06‑15), the
> exact commit the *Upgrading from v2.2.1 to v3.0.0* guide was written against, so the two guides
> tile with no gap and no overlap. 54 commits and 7 merged PRs (#71, #75, #77, #78, #79 plus two
> branch merges). The version string moved 3.0.0 → 3.0.1 → 3.0.2 → 3.1.0 inside this window; all
> of those intermediate bumps are covered here.
>
> Every claim below was verified against both trees by direct `grep`, and the pitfalls in
> [§8](@ref "8. Pitfalls when upgrading to v3.1.0") were confirmed by running unmodified v3.0.0 input
> files under a freshly built v3.1.0 binary.

---

## 0. TL;DR — the v3.1.0 upgrade checklist

**If you only run LaMEM from `.dat` files, you are almost certainly done already.** Items 1–3 are
the only user-facing breaks, and two of them are narrow.

1. ⚠️ **`out_gradient` was removed and is now silently ignored.** If your `.dat` requested gradient
   output, you will simply stop getting it — no error, no warning, identical residuals.
   **→ [§5.1](@ref "5.1 Removed: out_gradient (silent)")**
2. ⚠️ **`Adjoint_FieldSensitivity` was removed and is now silently ignored.** The adjoint
   sensitivity-kernel feature was purged wholesale.
   **→ [§5.2](@ref "5.2 Removed: Adjoint_FieldSensitivity (silent)")**
3. **The permeability report file is now `<name>.darcy.out`, not `<name>.out`.** Post-processing
   that globs for the old name finds nothing.
   **→ [§5.3](@ref "5.3 Permeability report renamed")**
4. **`make purge` no longer exists** — in `test/` *or* in any `examples/` Makefile. Use
   `make clean`. **→ [§6.1](@ref "6.1 Makefile targets replaced")**
5. **The test framework was rewritten.** `maintenance = true` / `update_expected = true` inside
   `runtests.jl` are gone, replaced by `make test` / `make work` / `make update`; every `@testset`
   now needs a `should_run_test(...)` guard. **→ [§6](@ref "6. Test framework & examples")**
6. **If you patched LaMEM's source, expect real work.** All global scratch buffers were removed and
   29 header prototypes changed. **→ [§7](@ref "7. If you patched the v3.0.0 source")**
7. **New developer tool dependencies:** `astyle` and `julia` are now needed for `make check` in
   `src/`. **→ [§1](@ref "1. Requirements & build changes")**

### Quick self-check for a v3.0.0 input file

Walk these top to bottom on one `.dat` file. Each item names what you would otherwise never notice:

```bash
grep -n "out_gradient"             my_model.dat   # → hit? output silently disappears (§5.1)
grep -n "Adjoint_FieldSensitivity" my_model.dat   # → hit? silently ignored (§5.2)
grep -rn "\.out"                   my_postproc/   # → globbing permeability reports? renamed (§5.3)
grep -rn "make purge"              my_scripts/    # → target no longer exists (§6.1)
grep -rn "update_expected\|maintenance" my_tests/ # → flags removed from runtests.jl (§6.2)
```

Nothing else in the `.dat` vocabulary changed. Of the 472 parameters LaMEM v3.0.0 could parse, 470
are unchanged, 2 were removed, and 1 was added.

---

## 1. Requirements & build changes

**Unchanged:** C++17, PETSc 3.19–3.25, MPI, and the same supported platforms (x86_64, arm64/aarch64,
armv6l/armv7l, i686, powerpc64le, amd64/MSYS). The compile commands are identical:

```bash
export PETSC_OPT=/path/to/petsc/optimized
export PETSC_DEB=/path/to/petsc/debug
cd src && make mode=opt all
```

**New — static-check targets** (`src/Makefile:138–177`):

| Target | Effect | Needs |
|--------|--------|-------|
| `make format` | Rewrite `*.cpp`/`*.h` in place per `src/.astylerc` | `astyle` |
| `make checkformat` | Fail if `make format` would change anything | `astyle` |
| `make checkgetrestore` | Verify every PETSc `*GetArray*` pairs with its `*RestoreArray*` inside one function | `julia` |
| `make check` | `checkformat` + `checkgetrestore` | both |

These are enforced upstream, so a PR that has not run `make check` will be sent back. `astyle` and
`julia` are therefore **new build-time dependencies for contributors** — they are not needed to
merely compile and run LaMEM.

`src/.astylerc` (new file) pins Allman braces, tab indentation at width 4, and `convert-tabs`
(leading indent in tabs, in-line alignment in spaces).

!!! tip "Gotcha when writing new code"
    astyle flattens hand-aligned continuation lines unless they sit inside parentheses. Rather than
    fighting it, hoist long expressions into named temporaries:

    ```c
    // astyle will re-indent this continuation:
    gx = (topo[J][I2] - topo[J][I1])
         / (COORD_NODE(I2, sx, fs->dsx) - COORD_NODE(I1, sx, fs->dsx));

    // stable under `make checkformat`:
    dx = COORD_NODE(I2, sx, fs->dsx) - COORD_NODE(I1, sx, fs->dsx);
    gx = (topo[J][I2] - topo[J][I1])/dx;
    ```

---

## 2. What's new in v3.1.0 (the "why")

v3.1.0 is dominated by one piece of engineering: **removing every global scratch buffer from the
code**. In v3.0.0, long-lived `Vec`s hung off `JacRes` and `FDSTAG` and were reused as working
memory by whatever routine needed them. That is fast and allocation-free, but it makes ownership
implicit — two routines can quietly stomp on the same buffer, and a leaked or unrestored array is
invisible.

v3.1.0 replaces them with explicit borrow-and-return, backed by an automated checker. The other
changes follow from the same impulse toward mechanical verification:

- **`scripts/petsc_getrestore_check.jl`** (new) — static checker pairing `*GetArray*` /
  `*RestoreArray*` calls within each function body, wired to `make checkgetrestore`.
- **Enforced formatting** — `src/.astylerc` + `make checkformat`, ending style drift in review.
- **Valgrind integration** — `test/check_valgrind_report.jl` (new) plus `make grind` /
  `make report`; the old `scripts/bash/ValgrindCheck.sh` was deleted.
- **PETSc object-balance checking** — `make check` in `test/` runs the suite under `-log_view` and
  fails on any object type created a different number of times than it was destroyed.
- **Test selectors and explicit modes** — run one test instead of all 36, and stop hand-editing
  flags.

Two small runtime features also landed: `continue_on_fail` and `-snes_track_stages`
([§4](@ref "4. New solver options")).

---

## 3. Repository layout changes in v3.1.0

Minimal this time. Added:

| Path | What |
|------|------|
| `src/.astylerc` | astyle configuration for `make format` / `checkformat` |
| `scripts/petsc_getrestore_check.jl` | Get/Restore pairing checker |
| `test/check_valgrind_report.jl` | Valgrind XML report checker for `make grind` / `make report` |
| `doc/src/man/Upgrade_v2.2.1_to_v3.0.0.md` | The previous migration guide, now shipped in-tree |

Removed:

| Path | Replacement |
|------|-------------|
| `scripts/bash/ValgrindCheck.sh` | `make grind` / `make report` in `test/` |
| `test/t08_AdjointGradients/AdjointGradients_SensitivityKernel*.dat`, `PSDKernel*` | none — feature purged |
| `test/t05_Perm/permea.darcy.dat` | none |

No directories were added, moved, or removed.

---

## 4. New solver options

Two additions to the `<SolverOptionsStart>` block. Both are optional and default to off, so they
cannot break an existing file.

### 4.1 `continue_on_fail`

```
<SolverOptionsStart>
    continue_on_fail = 1     # default 0
<SolverOptionsEnd>
```

Maps to PETSc's `-snes_continue_on_fail` (`src/options.cpp:317`): the simulation continues after a
diverged nonlinear solve rather than aborting, unless a severe reason is detected.

!!! warning "Use it deliberately"
    Upstream added a warning against enabling this in the test suite — a test that continues past
    divergence can report a pass while producing meaningless numbers. It is a production convenience
    for long runs you do not want to lose, not a way to make a failing model pass.

### 4.2 `-snes_track_stages`

A raw PETSc-style flag (`src/LaMEMLib.cpp:586`) enabling staged logging. Off by default and active
only in normal runs.

---

## 5. Input-file breaking changes

All three are narrow. The first two are **silent** — the model runs to completion and gives you a
result that is simply missing something you asked for.

### 5.1 Removed: out_gradient (silent)

Removed from `src/paraViewOutBin.cpp` by the sensitivity-kernel purge. In v3.0.0 it toggled gradient
output in the ParaView binary writer.

**Observed behaviour under v3.1.0** — an unmodified v3.0.0 input file with `out_gradient = 1`
appended:

```
exit code 0
no error, no warning naming the parameter
|mRes|_2 identical to the same run without the line
```

The only "unused options" warning printed is PETSc's routine complaint about `-ParamFile`, which
appears in a clean run too. **Nothing tells you the option was dropped.** If you relied on gradient
output, it is gone and you must notice by its absence in the output files.

### 5.2 Removed: Adjoint_FieldSensitivity (silent)

Removed from both `src/JacRes.cpp` and `src/adjoint.cpp`. This was the entry point to the adjoint
**field-sensitivity / sensitivity-kernel** feature, which was purged as a whole — the associated
test inputs (`AdjointGradients_SensitivityKernel.dat`,
`AdjointGradients_SensitivityKernel_PSD.dat`, `PSDKernelPaper.dat`, and the ParaView state/camera
files) were deleted with it.

Behaviour is the same as §5.1: silently ignored, exit code 0.

The rest of the adjoint machinery (`adjoint.cpp`, gradient inversion, scaling laws) is untouched —
only the field-sensitivity kernel path was removed. One related cosmetic change: the documented
default for `Adjoint_ScalingLawFilename` in `info/options/input_file.dat` changed from
`ScalingLaw_Test.dat` to `ScalingLaw_Test.out`.

### 5.3 Permeability report renamed

`src/JacResAux.cpp:502` — the Darcy/permeability report file extension changed:

```
<out_file_name>.out   →   <out_file_name>.darcy.out
```

Loud in the sense that the file is plainly there under a new name, but silent to any script that
globs the old one. Update post-processing that collects permeability reports.

---

## 6. Test framework & examples

This is where an existing workflow is most likely to break, because the commands themselves changed.

### 6.1 Makefile targets replaced

`make purge` **no longer exists** — neither in `test/` nor in `examples/BuiltInSetups`,
`examples/Localization`, `examples/PeriodicFreeSurface`, or `examples/SubductionWithParticles`.
Every example Makefile was reduced to a single `clean` target that now also sweeps `*.out`, `*.log`,
`*.png`, and `*.vts` recursively.

The `test/` Makefile gained a full set of modes:

| Target | Effect |
|--------|--------|
| `make test` | Run tests, delete generated output (`mode=test`) |
| `make work` | Run tests, **keep** generated files for inspection (`mode=work`) |
| `make update` | Run tests and **OVERWRITE** the `.expected` files (`mode=update`) |
| `make grind` | Run under Valgrind, then check the report for leaks / uninitialised values |
| `make report` | Re-check `.xml` from a previous `make grind` without re-running it |
| `make check` | Run under `-log_view`; fail on PETSc object creation/destruction mismatch |
| `make clean` | Remove leftover generated files |

### 6.2 Regenerating expected files

**Before (v3.0.0)** — hand-edit the top of `runtests.jl`, run, then remember to revert:

```julia
maintenance     = true
update_expected = true    # OVERWRITES all .expected files
```

**After (v3.1.0)** — the mode comes from the Makefile target, so there is no flag to forget:

```bash
make update 37     # regenerate only test 37's expected file
make update        # regenerate ALL of them — rarely what you want
```

`runtests.jl:104` parses `mode=test|work|update` from `ARGS` into `test_mode`, then derives
`update_expected = (test_mode == "update")` and `clean_files = (test_mode != "work")`.

### 6.3 Test selectors

New in v3.1.0 (`runtests.jl:60`, `runtests.jl:83`). Any bare number or hyphenated range after the
target selects a subset; other arguments (`is64bit`, `valgrind`, …) are ignored as selectors:

```bash
make test 01 05 32       # only t01, t05, t32
make test 03-07 11       # a range plus a single test
make update 37           # regenerate one expected file
make check 12-17         # object-balance check over a range
```

With no selector the whole suite runs, so existing `make test` invocations behave as before.

### 6.4 Registering a test now requires a guard

Every `@testset` sits inside an `if should_run_test("t<N>_<Name>")` block. **A testset without the
guard runs on every invocation**, including `make test 05`, silently defeating selection. Note the
two closing `end`s:

```julia
if should_run_test("t37_my_test")
@testset "t37_my_test" begin
    ...
end
end
```

`should_run_test` parses the number with `^t0*(\d+)_` and *fails open* — a name it cannot parse
always runs.

### 6.5 perform_lamem_test gained two parameters

`valgrind` and `memcheck`, defaulting to the `use_valgrind` / `use_memcheck` globals set by
`make grind` and `make check`. Existing call sites need no change.

Also worth knowing: `clean_dir=true` (the default) now deletes `*.bin` along with `Timestep*/`,
`*.pvd`, `*.out`, `*.log`, `markers*`, and `restart` (`test_utils.jl:501`). **A test that generates
its own topography/marker binary and then cleans it up by hand will fail with `ENOENT`** — the
framework has already removed it. Delete such manual cleanup.

---

## 7. If you patched the v3.0.0 source

This is the expensive part of the upgrade, and there is no way around it: **29 function prototypes
in `src/*.h` changed**, and every global scratch buffer was removed. A patch written against v3.0.0
internals will not apply cleanly.

### 7.1 Global buffers are gone

These members no longer exist on `JacRes` / `FDSTAG`:

```c
Vec gvx,  gvy,  gvz;                      // global velocity
Vec lvx,  lvy,  lvz;                      // local (ghosted) velocity
Vec gfx,  gfy,  gfz;                      // global forces
Vec lfx,  lfy,  lfz;                      // local (ghosted) forces
Vec ldxx, ldyy, ldzz, ldxy, ldxz, ldyz;   // local strain-rate components
Vec       gdxy, gdxz, gdyz;               // global strain-rate components
Vec gp, lp, lp_lith, lp_pore;             // pressure, lithostatic, pore
Vec gc, lbcor;                            // corner buffers
Vec lT;                                   // temperature box stencil
Vec lgradfield, phi;                      // gradient field / PSD context
```

Code that reached for `jr->lvx` must now borrow a vector and give it back
(`src/fdstag.cpp:1690`, `1704`, …):

```c
Vec vx, vy, vz;

PetscCall(FDSTAGGetLocalVectorFace(fs, &vx, &vy, &vz));
//  ... use them ...
PetscCall(FDSTAGRestoreLocalVectorFace(fs, &vx, &vy, &vz));
```

Full replacement API:

| Function | Returns |
|----------|---------|
| `FDSTAGGetLocalVectorFace` / `Restore…` | ghosted face-centred `vx, vy, vz` |
| `FDSTAGGetGlobalVectorFace` / `Restore…` | global face-centred `vx, vy, vz` |
| `FDSTAGGetLocalVectorEdge` / `Restore…` | ghosted edge-centred `vxy, vxz, vyz` |
| `FDSTAGGetGlobalVectorEdge` / `Restore…` | global edge-centred `vxy, vxz, vyz` |
| `FDSTAGGetLocalVectorCenter` / `Restore…` | ghosted cell-centred `vxx, vyy, vzz` |
| `DMGetLocalVectorClean(dm, &g)` | a zeroed local vector |
| `FDSTAGCombineVectors` / `FDSTAGSplitVectors` | block ↔ monolithic conversion |
| `FDSTAGSetEdgeCornerCenter` / `FDSTAGSetEdgeCornerFaces` | edge/corner fill-in |

**Every Get must be matched by its Restore inside the same function.** `make checkgetrestore`
enforces this and will reject a PR that leaks one. The checker is a lexical brace-scope scanner, not
a full C parser, so a Get/Restore split across `#ifdef` branches is reported as unmatched — review
flagged cases by hand rather than assuming a false positive.

### 7.2 CHKERRQ is gone

`CHKERRQ` no longer appears anywhere in `src/`. Use `PetscCall(...)` and `PetscCallMPI(...)`. If
your patch mixes the two, it will still compile, but `make format` and review will push back.

### 7.3 Removed helpers

Among the 29 changed prototypes, some functions disappeared outright rather than changing shape —
for example `JacResCopySol` and `JacResCopyVel`. Diff the headers you depend on before porting:

```bash
git diff eaf75d4b 7e7a012e -- src/JacRes.h src/fdstag.h src/tools.h
```

### 7.4 Format your patch

Run `make check` in `src/` before submitting. Both `checkformat` and `checkgetrestore` are gates,
and reformatting after the fact produces a noisy diff that is hard to review.

---

## 8. Pitfalls when upgrading to v3.1.0

### Silent: a removed parameter changes your results without saying so

`out_gradient` and `Adjoint_FieldSensitivity` are accepted, ignored, and never mentioned. This was
confirmed empirically: an unmodified v3.0.0 input with both parameters appended ran to exit code 0
and produced residuals identical to the same run without them, with no message naming either
parameter.

**There is no diagnostic for this.** LaMEM's "unused database option" warning covers PETSc
command-line options, not `.dat` keys. Grep your input files
(see [§0](@ref "Quick self-check for a v3.0.0 input file")).

### Silent: a testset without should_run_test ignores selectors

It will run during `make test 05` and every other selective invocation. The suite still passes, so
nothing draws attention to it; you only notice when a "quick" single-test run takes twenty minutes.

### Silent: post-processing that globs `*.out` for permeability reports

The file is now `<name>.darcy.out`. Your script finds nothing and, depending on how it is written,
may report zero results rather than an error.

### Loud: make purge no longer exists

```
make: *** No rule to make target 'purge'.  Stop.
```

Use `make clean`.

### Loud: manual cleanup of a generated test input

```
IOError: unlink("t37_my_test/t37_topo.bin"): no such file or directory (ENOENT)
```

`clean_dir=true` already removed it. Delete the manual `rm`.

### Loud: make check fails in src/ on a fresh patch

```
ERROR: source files are not properly formatted.
Formatted  /path/to/src/yourfile.cpp
Run 'make format' locally and commit the changes.
```

or

```
N potential mismatch(es) found.
```

Run `make format`, and pair up any unmatched `*GetArray*` call.

### Not a pitfall: your v3.0.0 input files

Worth stating plainly, because the v2.2.1 → v3.0.0 upgrade conditioned people to expect trouble: an
unmodified v3.0.0 input file runs correctly under v3.1.0. Verified with
`test/t01_FB1_Direct/FallingBlock_Direct_Default.dat` taken from `eaf75d4b` and run unchanged
against a v3.1.0 binary — exit code 0, no errors.

---

## Appendix: v3.1.0 parameter change reference

Extracted with the `get*Param` vocabulary diff and **verified individually by direct `grep` in both
trees** — bulk diffs produce false positives for parameters read through indirection.

| Parameter | v3.0.0 (`eaf75d4b`) | v3.1.0 (`7e7a012e`) | Status | Notes |
|-----------|---------------------|---------------------|--------|-------|
| `continue_on_fail` | — | `src/options.cpp:317` | **Added** | `<SolverOptionsStart>`; maps to `-snes_continue_on_fail`; default 0 |
| `out_gradient` | `src/paraViewOutBin.cpp:338` | — | **Removed** ⚠️ | Silently ignored |
| `Adjoint_FieldSensitivity` | `src/JacRes.cpp:323`, `src/adjoint.cpp:443` | — | **Removed** ⚠️ | Silently ignored; sensitivity-kernel feature purged |

Totals: 472 parameters parseable in v3.0.0, 471 in v3.1.0 — 470 unchanged, 2 removed, 1 added.

Non-`get*Param` user-visible changes:

| Change | Location | Kind |
|--------|----------|------|
| `-snes_track_stages` staged logging flag | `src/LaMEMLib.cpp:586` | Added |
| Permeability report `.out` → `.darcy.out` | `src/JacResAux.cpp:502` | Renamed |
| `Adjoint_ScalingLawFilename` documented default `.dat` → `.out` | `info/options/input_file.dat` | Doc |

---

*Guide generated by diffing upstream `eaf75d4b` (2026‑06‑15) against `7e7a012e` (2026‑08‑12) and by
executing unmodified v3.0.0 inputs under a v3.1.0 `bin/opt/LaMEM` binary. Parameter claims were
confirmed by direct `grep` in both trees rather than from the bulk vocabulary diff.*
