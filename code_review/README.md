# LaMEM source code review

A critical review of `src/`, looking for bugs, problems and misconceptions.
This folder is tracked on the `code-review` branch only; keep it off `master`.

## How to resume (read this first)

1. Find the next task marked `todo` or `in progress` in the status table.
2. Open its task file (`NN-*.md`). If it doesn't exist yet, create it from the template below.
3. Note the current commit (`git rev-parse --short HEAD`) in the table.
4. Add each issue to the task file first. Add it to `findings.md` once it has been checked: take the next
   free ID of its group (S = major, M = minor, N = note) and put it in that group's table and
   details. New IDs go at the end; the quick/medium/design order is only for existing IDs.
5. Update the status table before the session ends.

## Fixing findings: build and test

- Check each finding in the code before fixing it; some fix ideas in `findings.md` are wrong.
- Always build LaMEM in both modes after a change: `make mode=deb all` and `make mode=opt all`
  in `src/`. Tests use both binaries (most use opt, a few use deb). Hand runs use `bin/opt` or
  `bin/deb` directly, so a stale binary gives wrong conclusions.
- Use test selectors and run only the tests that touch the changed code: `make test 30 31` or
  `make test 03-07` in `test/`. The full suite takes more than 10 minutes.
- If a fix changes results, say so in the commit. `make update NN` rewrites all references in
  those test folders; keep them all (references of tests that passed before do the same job).
- If the build fails to link with undefined symbols, the objects are stale (built against
  another PETSc): run `make mode=deb clean_all` (or `mode=opt`) in `src/`, then rebuild.

## Status

| # | Area | Files | ~Lines | Status | Date | Commit |
|---|------|-------|-------:|--------|------|--------|
| 1 | Grid and parallel layout | fdstag, tools, Tensor | 3,600 | done | 2026-10-16 | 4835b6ba |
| 2 | Equations and assembly | JacRes, JacResAux, JacResTemp | 3,500 | done | 2026-10-16 | 4835b6ba |
| 3 | Material laws | constEq, phase, meltParam | 3,000 | done | 2026-10-16 | 4835b6ba |
| 4 | Solvers | lsolve, nlsolve, multigrid, matrix, matFree, matAux, matBFBT, matData | 7,800 | done | 2026-10-16 | 4835b6ba |
| 5 | Markers and advection | marker, advect, AVD, interpolate, cvi, subgrid | 9,100 | done | 2026-10-16 | 4835b6ba |
| 6 | Boundary conditions | bc | 3,100 | done | 2026-10-16 | 4835b6ba |
| 7 | Phase transitions and dikes | phase_transition, dike | 3,200 | done | 2026-10-16 | 4835b6ba |
| 8 | Free surface and erosion | surf, fastscape | 5,800 | done | 2026-10-16 | 4835b6ba |
| 9 | Adjoint and inversion | adjoint, objFunct | 4,600 | done | 2026-10-16 | 4835b6ba |
| 10 | Input, time stepping, driver, passive tracers | parsing, options, scaling, tssolve, LaMEMLib, LaMEM, passive_tracer | 5,700 | done | 2026-10-16 | 4835b6ba |
| 11 | Output | outFunct, paraViewOutBin, paraViewOutAVD, paraViewOutMark, paraViewOutPassiveTracers, paraViewOutSurf | 4,600 | done | 2026-10-16 | 4835b6ba |

Statuses: `todo`, `in progress`, `done`.

Line numbers in `findings.md` were updated to commit 6124aa7d on 2026-10-16 (task files still use 4835b6ba).

## What we look for

- **Bugs:** wrong indices, off-by-one errors, missing ghost-cell updates, memory leaks or
  double frees, values that are never set
- **Physics and maths:** wrong formulas, wrong units or scaling, sign errors, terms missing
  from the Jacobian
- **Parallel issues:** results that depend on the number of CPUs, missing reductions across
  processes
- **Weak spots:** unchecked inputs, silent fallbacks, magic numbers, dead code
- **Misconceptions:** comments or docs that don't match what the code does

## Severity

- **critical:** wrong results or crashes in normal use
- **major:** wrong results in some setups, or a real risk of failure
- **minor:** an edge case, unclear code, or a small inaccuracy
- **note:** style, cleanup or a question to check

Confidence: **confirmed** (checked in the code, or reproduced) or **possible** (needs checking).

## Tags

Each finding has a `- Tags:` line: severity, confidence, kind (one or two), effort.
Example: `- Tags: major, confirmed, logic, quick`. Find all quick fixes with
`grep -B1 "Tags:.*quick" findings.md`.

Kind:
- **physics:** wrong formula, units, scaling or sign
- **logic:** wrong index, off-by-one, overwritten or stale value
- **parallel:** results change with CPU count, hangs, missing exchange between CPUs
- **memory:** out-of-bounds access, leaks, use after free
- **input:** typos or bad values silently accepted, missing checks, docs that don't match input
- **restart:** wrong behaviour after a restart
- **cleanup:** dead code, wrong comments or messages

Effort:
- **quick:** a one-line or local fix
- **medium:** one function or file
- **design:** needs a decision or a bigger rewrite

## Finding format (used in `findings.md`)

```
### S-001 [major, confirmed] Short title
- Tags: major, confirmed, logic, quick
- Where: src/File.cpp:123 (commit xxxxxxx)
- What: ...
- Why it matters: ...
- Fix idea: ...
- Status: open | fixed | not a bug
```

## Task file template (`NN-name.md`)

```
# Task N: <Area>

Commit: xxxxxxx   Started: YYYY-MM-DD   Status: in progress

## Coverage
- [ ] file.cpp: functions checked / skipped

## Notes and open questions

## Issues (draft)
```
