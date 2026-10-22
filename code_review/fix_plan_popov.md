# Fix plan for findings assigned to Anton Popov

Plan for the 83 findings assigned to Anton in `findings_audit_request_feedback.md`
(80 by authorship, 3 by area). Made 2026-10-20 at commit 61bb7aa0. Each finding is checked in
the code before it is fixed. Fixes that change results are reported; reference outputs are not
updated silently.

Testing:
- Always build LaMEM in both modes after a change: `make mode=deb all` and `make mode=opt all`
  in `src/`. Tests use both binaries (most use opt, a few use deb). Hand runs use `bin/opt` or
  `bin/deb` directly, so a stale binary gives wrong conclusions.
- Use test selectors and run only the tests that touch the changed code: `make test 30 31` or
  `make test 03-07` in `test/`. Run the full suite (>10 min) only when asked.
- `make update NN` rewrites all references in those test folders. Keep them all (no need to
  revert the ones whose tests passed before; they do the same job).

## A. Fix right away (small, clear, local)

| Group | IDs | Fix |
|---|---|---|
| One-line logic | S-001, S-007, S-010, M-017, M-027, M-037, M-040 | `fixPhase == -1` check, `abs(v)`, `> 1`, `>=`, `\|\|`, move `free()` |
| Input data | S-003 | `Vn = 0` (Boris, GeoParams) |
| Grid | S-006 | stretch global box about `ref`; add a test |
| Inflow BC | S-011 | store base velocity in BCCtx, not a `static` |
| Markers | S-013, S-019, M-020, M-022, M-025, M-026, M-031 | sort `idel`, remap after perturbation, empty/NULL checks, free leak |
| Parallel | M-038 | total cell count (sum over ranks) |
| Polygons | S-018 + M-012 | same bug, fix together |
| Surface (Tobias) | S-004, M-005 | Boris: fix |
| Small | M-018, M-041, N-013 | missing SNES/KSP reasons, truncate .pvd, guard division |
| Shared | S-024 (Andrea), M-008 (Wenrong) | trivial; taken over, unassigned from them |

## B. Fix next (more work, clear fix)

- S-009: "user velout given" flag; default stays mass-conserving.
- S-029: fill boundary ghosts of lithostatic and pore pressure.
- S-035, M-021, M-058, M-059: AVD marker control, one batch.
- S-039: read solver-default inputs with command-line options applied. Done (no clear).
- M-016, M-030: reject typos and bad values.
- M-051, M-053, M-063: one helper so a one-rank error stops all ranks (no hang). Done in
  batch 4: no hang (PETSc aborts all ranks); unchecked reads fixed.
- M-055, M-056, M-061: periodic seam and coarse-level BCs in matrix-free / MG.
- M-060, M-066: buffer size checks.
- M-069: average doubled-angle form.
- N-024, N-025, N-026: the parts Boris approved.

## C. Discuss first

1. S-002 (`Lm`): code uses `Le/Lm`; fixed softening displacement suggests `Lm/Le`. Depends on the
   intended meaning of `Lm`.
2. S-046 (dt): fix changes all viscoelastic results; test references need updating.
3. S-044 (`fix_phase` on processor boundaries): ghosted flags vs ADD/MAX before scatter.
4. S-045 (`exy_num_periods`): really wrap markers, or document deletion.
5. S-036 (Clapeyron units in SI): with Andrea; input unit per unit system.
6. M-048 (T = 0 means unset): DBL_MAX sentinel may change existing inputs.
7. Intended? M-047, M-049, N-003, N-016, N-019.
8. With Garrett: M-046, M-067.
9. Low value: M-009, M-023, M-054, M-057, N-008, N-014 (code or docs?), N-017, N-020.
10. Delete dead code? N-009, N-010, N-015, dead parts of N-021, N-022.
11. M-063 (from B): "Incorrect sedimentation phase" stop in subgrid.cpp when the surface rises
    more than one cell per step or has overhangs. Relax, or keep the stop?
12. Nothing to do: M-019 (mitigated by `-snes_max_funcs unlimited`); N-028, N-029 (Boris: leave).
    Small clean-ups: N-002, N-007 (mostly Wenrong/Yun), N-011, N-018.

## Moved from A to C while fixing

- S-014: the local wall tests also protect the upper side. With one ghost layer, node `nx+1`
  does not exist at an internal boundary, so global tests alone would read out of range.
  Needs a 2-layer ghost vector or a different stencil.
- N-004: `eta && Bd` is accepted on purpose for restarts (phases are re-read into a struct
  that still holds the old Bd/Bn, see c8ec3931). A strict check needs a decision on restart.

## Status of A (2026-10-20)

Done, not committed: S-001, S-003, S-004, S-006 (test not yet added), S-007, S-010, S-011,
S-013, S-018/M-012, S-019 (perturbed coordinate kept in host cell), S-024, M-005, M-008,
M-017, M-018, M-020, M-022, M-025, M-026 (APS still not saved), M-027, M-031, M-037, M-038,
M-040, M-041, N-013. S-011 adds a field to BCCtx, so older restart files no longer load.
M-041: only the .pvd is truncated; stale `Timestep_*` folders after a restart from an earlier
checkpoint stay. Checkpointing at output steps was discussed (checkpoint ~20-150x one output
step, mostly markers) and dropped for now.

## Open items after group A (committed 14425854..876978e4)

- `FindPointInCellAdjoint` (same bug as M-022): fixed. Copy deleted, `FindPointInCell` declared
  in AVD.h, `< 0` checks added. The 3 topography loops in `AdjointPointInPro` tested `z` against
  the global box, so with more than one rank in z they stopped with "Non-local marker"; they now
  test the local box, as in `surf.cpp`.
- AVD batch from B (not committed):
  - S-035: delete path now counts claimed cells per point; smallest Voronoi volumes go first.
  - M-058: one helper `AVDGetLimitsMV`; `nmax` is reduced like `nmin`. Per-rank half-volumes
    kept on purpose (a volume cut by a rank boundary is half-size on each rank).
  - M-059: edge volumes are now true edges (shifted in the two directions other than `dir`);
    corner volumes get a quarter of the limits; log labels fixed. Changes results under `avd`.
  - M-021: fine as is (no parent marker for an empty volume; sparse ones fill over steps).
  - New (found while testing): `AVDCellInit` clamped only the upper border; a marker exactly on
    the lower bound rounded into the mask cell and stopped the run ("Inserting cells into
    boundary cells"). Hit by injected markers at cell-centre coordinates. Lower clamp added.
  - No test uses `mark_ctrl = avd`. Checked by hand on t15 RTI (32x2x32, 1 and 2 ranks,
    `-nmark_lim 20,60` for deletion and `200,800` for injection): runs clean. `make test` passed.
- Tests: `make test` passed after A (115 pass, t37 FastScape broken as before). No
  reference output changed, so the fixed paths (S-006, S-007, S-013, S-024, M-008) are not covered.

## Lowest priority (or ask Boris)

- S-006 test: pure shear with non-zero `bg_ref_point`, check the global box (Boris asked for it).

## Order of work

One commit per group of A, then B, then C after discussion.

## Group B, batch 1: inputs and BCs (2026-10-21, not committed)

- S-009: new flag `velout_user` in BCCtx (restart format changes again). `BCGetVelins` computes
  `velout` from mass balance in one place, only when `bvel_velout` is not given. The log prints
  the outflow velocity when it is given by the user.
- M-016: unknown `bvel_face` and `bvel_temperature_inflow` strings stop the run. `exy` with
  no-slip on front/back is an error; no-slip on bottom/top for vx is now applied with `exy`.
  Note: `info/options/input_file.dat` (option listing, not a test) sets both.
- M-030: marker phases are checked right after setup (negative phase was used as an array index
  in `ADVMarkSetTempPhase`); `nmark_sub < 1`, negative `nmark_lim`/`nmark_avd`, `nmin > nmax`
  are rejected; `stagp_a = 0` is used as given (was replaced by 2/3).
- `make test` passed (t37 FastScape broken as before).

## Group B, batch 2: buffer checks (2026-10-21, not committed)

- M-066: schedule size is checked against `_max_num_steps_` (also in `TSSolGetPeriodSteps`,
  which writes n+1 values); first `time_dt_periods` entry must be 0; `TSSolAdjustSchedule`
  no longer writes past the buffer (no insert when full). Decided: the schedule overrides basic
  time stepping, so the run still ends after the scheduled steps; the log now prints the
  schedule end and a warning if `time_end` is later. Note: the CFL-inserted steps at the end
  let t30 run to t = 3 (schedule ends at 2); kept as is, t30 unchanged.
- M-060: polygon file checked for truncation, `Vol.num <= Nmax`, polygon length `<= Lmax`,
  at most `_max_polygons_` per slice; VLAs `SxAll`/`SyAll` replaced by heap arrays.
- Tests: t30, t31 pass (selected tests only).

## Group B, batch 3: S-029, M-069 (2026-10-21, not committed)

- S-029: new `FDSTAGSetBoundCenter` (fdstag.cpp) copies boundary cell values to the outside
  ghosts (zero gradient, as for pressure; periodic x ghosts kept), then calls
  `FDSTAGSetEdgeCornerCenter`. Called after `lp_lith` and `lp_pore` are computed. Changes
  results where lithostatic pressure is used at walls (creep default, `p_litho_plast`, pore
  pressure): t19 MigratingRidge_2D, t26 dyndike_4core, t29, t32 (FB, LR). References updated
  for t19 MigratingRidge_2D, t29, t32 FB/LR (`make update`, unrelated references reverted).
- M-069: SHmax/EHmax directions averaged to corners in doubled-angle form (output only),
  shared helper `PVOutPutAxialDir` in outFunct.cpp.
- Full suite before reference update: 107 pass, 5 fail (above), t37 broken as before.

### MAJOR open issue: t26 dyndike_4core test is not reliable

- Before S-029 the reference already needed exactly 50 = maxit nonlinear iterations at step 4
  (just converged). With S-029 it hits maxit; `nonlinear_tolerances` maxit raised 50 -> 80 in
  `dyndike_4core.dat` (not committed), then 51 and 48 iterations.
- Steps 4-5 are not reproducible between runs on 4 ranks (steps 1-3 agree to 1e-12, steps 4-5
  differ by 5-15% in |Div|_inf), so the updated reference passes only sometimes.
- The run prints ~1000 "Number of diverged points : ~72000" warnings (local plasticity
  iterations fail), with and without S-029; the old reference has none.
- User: "this test was always a big mess". Needs a proper look (setup, solver, local
  iterations), not just a tolerance change. Stopped here on 2026-10-21.

## Group B, batch 4: M-051, M-053, M-063 (2026-10-22, not committed)

- The hangs are not real. An error that reaches `main` makes PETSc (3.24) call `MPI_Abort`, which
  stops all ranks. Checked: failed `mkdir` in `DirMake` on 2 ranks (read-only run folder) stops at
  once with exit code 83. So no "stop all ranks" helper is needed.
- M-051: new `FileRead` (tools.cpp) reads a block and stops with "file is truncated or damaged"
  on a short read (zero size is a no-op). All 11 `fread` calls use it (restart readers in
  LaMEMLib, fdstag, bc, advect, phase_transition, `VecReadRestart`; input file in parsing).
  Checked: t39 restart file cut by 100 bytes stops in `ADVReadRestart`.
- M-053: `BCReadFixCell` checks the `stat` result, closes the file on the size error, reads with
  `FileRead`; its rank-local errors use `PETSC_COMM_SELF`.
- M-063: hang part not real (as above). The "Incorrect sedimentation phase" stop on fast uplift
  or overhangs is model behaviour: moved to C.
- Tests: t01, t39 pass (t41 skipped, plugin not built). No results change.

## Group B, S-039 (2026-10-22, not committed)

- `setSolverOptions` no longer clears the options before `solverOptionsSetDefaults`. The defaults
  now see the command line (grid, CPUs, `act_temp_diff`, solver block keys). The saved options are
  still put back last, so the command line still wins over the defaults.
- Keys inside `<SolverOptionsStart>` need the block index on the command line, e.g.
  `-num_mg_levels[0] 2` (as for all block keys).
- Checked on t02 with `-nel_x 12 ... -num_mg_levels[0] 4`: now stops with "Max. possible: 3"
  (12 cells); `-act_temp_diff 1` now sets `ts_ksp_*`.
- Tests: t02, t03, t07, t08, t19, t20, t28, t30 pass; t26 only dyndike_4core fails (no command-line
  options, known issue). No results change.
