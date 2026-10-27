# Fix plan for findings assigned to Anton Popov

Plan for the 83 findings assigned to Anton in `findings_audit_request_feedback.md`
(80 by authorship, 3 by area). Made 2026-10-20 at commit 61bb7aa0. Each finding is checked in
the code before it is fixed. Fixes that change results are reported; reference outputs are not
updated silently.

Testing:
- Always build LaMEM in both modes after a change: `make mode=deb surf=scape all` and
  `make mode=opt surf=scape all` in `src/`. Tests use both binaries (most use opt, a few use
  deb). Hand runs use `bin/opt` or `bin/deb` directly, so a stale binary gives wrong conclusions.
- Use test selectors and run only the tests that touch the changed code: `make test 30 31` or
  `make test 03-07` in `test/`. Run the full suite (>10 min) only when needed.
- `make update NN` rewrites all references in those test folders. Keep them all (no need to
  revert the ones whose tests passed before; they do the same job).

## Status (2026-10-27, checked against commits and code)

- A: done, committed (14425854..876978e4).
- B: done, committed (see each batch below for the commits).
- C done: items 1 (S-002, 5fb3350a), 2 (S-046, b41cef9a), 4 (S-045, ea809634),
  7 (92b20ac3; M-049 changes t10, 2eb1e72b), 9, 10, 12; S-023 and the local-iteration /
  T-coupling work; t16 and t26 workarounds removed; t16, t19 sped up (e38ee480).
- C done 2026-10-26: M-063 (bff40e9f, stop kept), N-022 (d3a18147, 88e7fbf5),
  S-014 (b7a5c143, cell-center velocities precomputed, see below).
- C done 2026-10-27: S-044 (1f57d0ed, ghosted fixed-cell flag, owned faces only);
  M-048 not an issue (air gets Ttop); N-004 (2c3f19d3, fresh phase read on restart, strict checks).
- C done 2026-10-27: N-021 rest (1a67766e, option order in nlsolve/multigrid, MFFD type
  set once, matData comment, lsolve bf-only checks, matAux duplicates). No results change
  (t01-t08, t15, t20 pass).
- C done 2026-10-27: 13, rest of N-025 (ea15ddba: nstep_out = 0 is output by dt_out only, < 0 no
  output, dt_out on a fixed time grid; no trimming at time_end) and N-026 (b71b97cd:
  total_pressure edge values; fwrite per value won't fix).
- C done 2026-10-27: M-056 2D (MatDataRestrictBC2D copies boundary ghost flags like 3D; no effect
  today, 2D has no rediscretized levels; t03 t16 t19 t22 t26 t28 t30 pass).
- C done 2026-10-27: S-036 (Clapeyron transitions require units = geo; geo unchanged).
- C open (checked in code 2026-10-27):
  - 8. M-046, M-067 (with Garrett).
- Core-count drift (found 2026-10-26 while testing S-014), fixed 2026-10-26: with
  `rand_noise 0`, 1-core and 4-core runs of
  examples/BuiltInSetups/Subduction2D_FreeSlip_DirectSolver.dat drifted apart (dt 0.23798 vs
  0.23777 at step 20). Cause: the random shift of markers injected in empty cell corners
  (ADVCheckCorners), seeded per rank and drawn in rank-local order. Noise removed; markers go to
  the corner center. Now both runs agree to ~1e-10 over 21 steps. Changes results (t16, t22,
  t27, t32, t34 references updated). Known and accepted: the parent marker is searched in local
  cells only (no ghost markers), so it can still differ next to rank boundaries.
- Lowest priority: S-006 test; M-026 APS not saved by ADVMarkSave.
- Older notes below say "t37 FastScape broken" or "dyndike fails (known)": both are resolved
  (t37 needs `surf=scape`, see CLAUDE.md; dyndike fixed with the lmaxit removal).

## Random noise seeding (checked 2026-10-26, dropped)

All random sources start from a fixed seed (PETSc default 0x12345678, unseeded rand() =
seed 1, FreeSurf_NoiseSeed, std::mt19937 default), none uses the clock. So runs are the same
from run to run on the same number of cores, which is all the tests need. Noise still depends
on the core count (each rank draws its own sequence). No code change; a note to always use a
fixed seed was added to the development guide (doc/src/man/LaMEM_Development.md, 6.3.2).
The corner-injection noise (seed 1 + rank + nproc*step) was removed later the same day, see
"Core-count drift" in the status above (6df70c04).

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

All done and committed.

- S-009: "user velout given" flag; default stays mass-conserving. Done (b7aa6557).
- S-029: fill boundary ghosts of lithostatic and pore pressure. Done (c7ce5d27).
- S-035, M-021, M-058, M-059: AVD marker control, one batch. Done (329482cf).
- S-039: read solver-default inputs with command-line options applied. Done (90dfb957, no clear).
- M-016, M-030: reject typos and bad values. Done (b7aa6557).
- M-051, M-053, M-063: one helper so a one-rank error stops all ranks (no hang). Done in
  batch 4 (f7cc8f88): no hang (PETSc aborts all ranks); unchecked reads fixed. M-063 rest in C.
- M-055, M-056, M-061: periodic seam and coarse-level BCs in matrix-free / MG.
  Batch 5 (2026-10-22), committed 7ffd9332 (M-055, M-056), 9a8c9eab (M-061): t02 subtests
  b (2 mat-free MG levels, 4 levels) and c (same + periodic, no-slip bottom) added to
  runtests.jl.
  M-055 fixed (matFree.cpp x checks skip periodic; matBFBT.cpp cfv guard removed);
  -js_mat_free and wBFBT periodic converge.
  M-056 fixed in 3D (matData.cpp MatDataRestrictBC3D: boundary ghost flags copied from the fine
  grid, x skipped when periodic). 2D left as is: rebuilt coarse levels need matrix-free MG,
  which is rejected in 2D. Subtest c now converges (8, 9 its; Galerkin only: 7, 9).
  t02 references made (`make update 02`); a changes only at round-off. Subtests b, c now use
  their own input files; added d (wBFBT + periodic, settings from FallingSpheres_wBFBT.dat):
  13, 21 its.
  M-061 fixed: interpolate.cpp skips the x clamps when periodic (Y/Z face, center, YZ edge).
  Checked by hand (block at x = 0..0.5): density at x = 0 and x = L is now 1.5, was 2 at x = 0.
  t02 passes (no test checks output files).
  Fixed: bc.cpp printed the periodic line inside the no-slip mask line, without a newline.
  Also: MatDataRestrictParam3D clamps at the periodic seam (Im1++); a 4-rank SETERRQ in MG
  setup ended as "signal 9" with no message printed.
- M-060, M-066: buffer size checks. Done (d9c7978e).
- M-069: average doubled-angle form. Done (c7ce5d27).
- N-024, N-025, N-026: the parts Boris approved. Done (59bf2236, 9977c260).
  Done in src: parsing.cpp checks strtol (base 10) / strtod results, detects overlapping
  blocks, checks ftell, drops dead `&& fb`; options.cpp passes the whole PETSc option line.
  Relaxed (agreed 2026-10-23): error only if nothing parses, or an integer is given as a real
  (`16.0`, `1e3`); other trailing characters (`;`) allowed; extra values ignored (2D inputs give
  3D values; bc.cpp:689 reads 1 value of bvel_velin). Skipped: rank-0 hang (not real), paths
  with spaces / `=` in values. Real bug: t26 dike_heating_kfac.dat / dike_heating_rhoA.dat had
  `xbounds = [-15.0 15.0]`, old code read `[-15.0` as 0; brackets removed. With the right bounds kfac
  stalled (|F| ~6e-4 after 50 its), so its nonlinear atol was raised 1e-4 -> 1e-3. t26
  references updated (dyndike_4core.expected kept). Scan of all 119 test/example inputs: no
  parser errors. Full suite: only t26 dyndike (known) and t37 (binary without FastScape) fail.
  N-025 (2026-10-23): its_ steady-state norm seeding jr->ts_ksp_ref_norm is intended (Anton):
  comment added in LaMEMLib.cpp only. Fixed: heat production label W/kg, "WARNING" typo and
  newline (scaling.cpp), GettingStarted.md uses nstep_max, set_scalar_option uses %.17g.
  Left for C: fix_dt schedule one step late, dt_out without effect when nstep_out = 1, no dt
  trimming at time_end / output times. nstep_out <= 0 switching off output is intended (Boris).
  N-026 (2026-10-23): energ_res was stale only in dry-run and first output without initial
  guess (in the time loop JacResViewRes recomputes it): now recomputed there. Passive tracer
  pressure gets pShift. out_avd_ref < 1 rejected; unclaimed AVD cells written as phase 255
  (was points[-1]). Negative aggregate phaseID rejected; output name uses snprintf.
  vol_rate, vorticity, ang_vel_mag implemented (from JacResGetVelGrad) with options
  out_vol_rate, out_vorticity, out_ang_vel_mag; checked on t18 simple shear (wz = 2e-15,
  angular velocity 1.82 deg/Myr, expected 1.81; vol_rate ~1e-22). plast_dissip kept (Boris).
  Not done: one value per fwrite in marker/tracer/AVD files, total_pressure edge values.

## C. Discuss first

1. S-002 (`Lm`): code uses `Le/Lm`; fixed softening displacement suggests `Lm/Le`. Depends on the
   intended meaning of `Lm`.
   Done 2026-10-26: bug confirmed (Anton); now `Lm/Le`, documented in input_file.dat. Committed 5fb3350a.
2. S-046 (dt): fix changes all viscoelastic results; test references need updating.
   Done 2026-10-26: `fix_dt` also on with temperature (dt-dependent solve); temperature-only test results changed. Committed b41cef9a.
3. S-044 (`fix_phase` on processor boundaries): ghosted flags vs ADD/MAX before scatter.
   Done 2026-10-27: ghosted flag (BCApplyFixCells). ADD does not fit: the free marker is DBL_MAX
   and non-zero constraints would double; MIN would work but touches all BCApply* routines.
   Only fixed regions with a corner on a rank-boundary corner changed. Committed 1f57d0ed.
4. S-045 (`exy_num_periods`): really wrap markers, or document deletion.
   Done 2026-10-26: periodic marker advection != periodic mechanical BC (Anton). Marker
   neighbour list made periodic in x for Exy (`FDSTAGGetNeighbProc(fs, ptx)`), outflow markers
   wrap instead of being deleted. t18 unchanged. Committed ea809634.
5. S-036 (Clapeyron units in SI). Done 2026-10-27 (Anton, without Andrea): any Clapeyron
   transition now requires `units = geo` (slope [MPa/K], T0 [C], P0 [Pa]); si and none stop with
   an error. Nobody runs mantle-scale models in SI. Geo was correct before, results unchanged.
   Comments and log labels fixed ([MPa/K], [C]). t16 t21 t29 pass (t41 skipped, plugin not built).
6. M-048 (T = 0 means unset): DBL_MAX sentinel may change existing inputs.
   Not an issue (Anton, 2026-10-27): air gets `Ttop`, not phase `T`; 0 [C] blocks are rare. No change.
7. Intended? M-047, M-049, N-003, N-016, N-019. Done, committed 92b20ac3 (t10 refs 2eb1e72b).
   M-047 (2026-10-26): setting only the pressure history (cells, markers) is intended (Anton);
   zero pressure sum (no gravity) now counts as converged, no NaN / false warning.
   M-049 (2026-10-26): G averaging kept (history stress is averaged the same way, so a non-elastic
   phase drops out); Kb now averages compressibility 1/Kb (non-elastic = incompressible).
   t03, t04 unchanged (all their phases have the same Kb).
   I2Gdt left as is: it goes away with the planned switch to a material point method (Anton).
   N-003 (2026-10-26): won't fix, soft eta_min/eta_max limits are intended and obvious (Anton).
   N-016 (2026-10-26): flat water table at average topography is intended (Anton); comment added.
   N-019 (2026-10-26): -snes_newton_maxit is now the Newton budget per solve; once used up, Picard
   only (no Picard/Newton ping-pong). MMFD label fixed to MFFD. Divergence test vs first residual
   of the step kept. t03, t13, t14, t26, t32 pass.
8. With Garrett: M-046, M-067.
9. Low value: M-009, M-023, M-054, M-057, N-008, N-014 (code or docs?), N-017, N-020.
   Done 2026-10-23: M-009 surface z clamped to box on bottom/top ranks; M-023 RK2 midpoint clamped
   (not deleted); M-054 closed box with unbalanced CompensatingInflow stops with an error (Anton:
   imbalance not allowed); M-057 seed per rank and step (noise removed 2026-10-26, 6df70c04), `>=` on upper side,
   total count printed;
   N-008 area-weighted avg_topo; N-014 Clapeyron parameters override a predefined law.
   N-017, N-020: non-issues (Anton). References updated: t22, t27, t32, t34; t15 RTI error
   0.00148 -> 0.00041 (closer to analytical). Only dyndike fails (known).
   Committed 9125982e.
10. Delete dead code? N-009, N-010, N-015, dead parts of N-021, N-022.
    Done 2026-10-23 (3c04c606), see "Group C, dead code" below.
11. M-063 (from B): "Incorrect sedimentation phase" stop in subgrid.cpp when the surface rises
    more than one cell per step or has overhangs. Relax, or keep the stop?
    Done 2026-10-26 (Anton: keep the stop). Clearer message (reduce the time step) in subgrid.cpp
    and passive_tracer.cpp; the `-2` stop mark in `ADVGetSedPhase` no longer overwrites a rock
    cell above an air cell (overhang). No results change (t03, t15, t20, t21, t24, t35, t36, t38).
12. Nothing to do: M-019 (mitigated by `-snes_max_funcs unlimited`); N-028, N-029 (Boris: leave).
    Small clean-ups: N-002, N-007 (mostly Wenrong/Yun), N-011, N-018.
    Done 2026-10-23: M-019, N-028, N-029 marked won't fix. N-002, N-007, N-011 (input checks
    only; initial-guess skip of boxes/cylinders kept), N-018 fixed. Committed a868c954.
13. Left from B (N-025, N-026): dt_out without effect when nstep_out = 1, no dt trimming at
    time_end / output times; one value per fwrite in marker/tracer/AVD files, total_pressure
    edge values. fix_dt schedule one step late: already fixed (Anton, 2026-10-26).
    N-025 done 2026-10-27: nstep_out = 0 is output by dt_out only, < 0 is no output; dt_out
    output on a fixed time grid; dike.cpp guards istep % nstep_out. No trimming at time_end or
    output times (Anton: not needed). t12 no longer reads output steps 1 and 3 (unused, step 3
    not saved any more).
    N-026 done 2026-10-27: total_pressure uses boundary ghost values (use_bound = 1) like
    pressure; checked on t01 with -pres_top 0.3 (old edge difference 0.017, now 0). fwrite per
    value: won't fix (stdio buffering, small cost).

## Moved from A to C while fixing

- S-014: the local wall tests also protect the upper side. With one ghost layer, node `nx+1`
  does not exist at an internal boundary, so global tests alone would read out of range.
  Needs a 2-layer ghost vector or a different stencil.
  Done 2026-10-26 without extra ghost layers: `ADVelGetCellVel` (cvi.cpp) computes cell-center
  velocities (face averages) on owned cells, exchanges one ghost layer of DA_CEN, then sets wall
  ghosts (tangential: two-point constraint of the lower face; normal: extrapolation through the
  wall velocity). STAGP interpolates `v_p` from these. Serial results equal the old code to
  rounding, except at upper walls (x right, y back, z top), where the old code used a 2-cell
  interval (`nxe = 2*ccx[nx]-ccx[nx-1]`, `ccx[nx]` is already the ghost cell); this changes t16.
- N-004: `eta && Bd` is accepted on purpose for restarts (phases are re-read into a struct
  that still holds the old Bd/Bn, see c8ec3931). A strict check needs a decision on restart.
  Done 2026-10-27: fresh read (Anton: re-reading is meant to override parameters); each phase
  is cleared in DBMatReadPhase, strict checks, eta_st default in DBMatSetDefaultEtaSt. Committed 2c3f19d3.

## Status of A (2026-10-20)

Done, committed 14425854..876978e4: S-001, S-003, S-004, S-006 (test not yet added), S-007, S-010, S-011,
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
- `FindPointInCellAdjoint` fix committed f8fb5175.
- AVD batch from B (committed 329482cf):
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

## Group B, batch 1: inputs and BCs (2026-10-21, committed b7aa6557)

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

## Group B, batch 2: buffer checks (2026-10-21, committed d9c7978e)

- M-066: schedule size is checked against `_max_num_steps_` (also in `TSSolGetPeriodSteps`,
  which writes n+1 values); first `time_dt_periods` entry must be 0; `TSSolAdjustSchedule`
  no longer writes past the buffer (no insert when full). Decided: the schedule overrides basic
  time stepping, so the run still ends after the scheduled steps; the log now prints the
  schedule end and a warning if `time_end` is later. Note: the CFL-inserted steps at the end
  let t30 run to t = 3 (schedule ends at 2); kept as is, t30 unchanged.
- M-060: polygon file checked for truncation, `Vol.num <= Nmax`, polygon length `<= Lmax`,
  at most `_max_polygons_` per slice; VLAs `SxAll`/`SyAll` replaced by heap arrays.
- Tests: t30, t31 pass (selected tests only).

## Group B, batch 3: S-029, M-069 (2026-10-21, committed c7ce5d27)

- S-029: new `FDSTAGSetBoundCenter` (fdstag.cpp) copies boundary cell values to the outside
  ghosts (zero gradient, as for pressure; periodic x ghosts kept), then calls
  `FDSTAGSetEdgeCornerCenter`. Called after `lp_lith` and `lp_pore` are computed. Changes
  results where lithostatic pressure is used at walls (creep default, `p_litho_plast`, pore
  pressure): t19 MigratingRidge_2D, t26 dyndike_4core, t29, t32 (FB, LR). References updated
  for t19 MigratingRidge_2D, t29, t32 FB/LR (`make update`, unrelated references reverted).
- M-069: SHmax/EHmax directions averaged to corners in doubled-angle form (output only),
  shared helper `PVOutPutAxialDir` in outFunct.cpp.
- Full suite before reference update: 107 pass, 5 fail (above), t37 broken as before.

### RESOLVED (2026-10-23): t26 dyndike_4core test was not reliable

- Cause: `lmaxit = 8` in the input. Bisection needs ~20 steps for `lrtol = 1e-6`, so ~72000
  points "diverged" per iteration and Newton hit maxit. See "Local iterations" below.
- Now 17-27 nonlinear iterations per step, three 4-rank runs agree to 6 digits; passes.

## Group B, batch 4: M-051, M-053, M-063 (2026-10-22, committed f7cc8f88)

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

## Group B, S-039 (2026-10-22, committed 90dfb957)

- `setSolverOptions` no longer clears the options before `solverOptionsSetDefaults`. The defaults
  now see the command line (grid, CPUs, `act_temp_diff`, solver block keys). The saved options are
  still put back last, so the command line still wins over the defaults.
- Keys inside `<SolverOptionsStart>` need the block index on the command line, e.g.
  `-num_mg_levels[0] 2` (as for all block keys).
- Checked on t02 with `-nel_x 12 ... -num_mg_levels[0] 4`: now stops with "Max. possible: 3"
  (12 cells); `-act_temp_diff 1` now sets `ts_ksp_*`.
- Tests: t02, t03, t07, t08, t19, t20, t28, t30 pass; t26 only dyndike_4core fails (no command-line
  options, known issue). No results change.

## Group C, dead code: N-009, N-010, N-015, N-021, N-022 (2026-10-23, committed 3c04c606)

- N-009: `Tensor2RNEigen`, `getISA`, `Tensor2RSSpectral` removed.
- N-010 (user: remove fully): `bias_x/y/z` no longer read; `biases` arrays, the "deprecated"
  stop and the biased branches in fdstag.cpp and fastscape.cpp removed. Old inputs with
  `bias_*` now run on a uniform mesh without notice; the upgrade guide (§5.2) says so. `bias_*`
  lines dropped from t16, t19, t26 inputs and info/options/input_file.dat.
- N-015: objFunct.cpp deleted; `ObsType`, `ObjFunct`, `_max_num_obs_` removed (objFunct.h still
  holds `ModParam`). adjoint.cpp: `cellConstEqFD`, `edgeConstEqFD`, `devConstEqFD`,
  `setUpPhaseFD` removed; undefined `AdjointGradientResetParameter` declaration removed; fields
  `Perturb`, `CurScalst`, `gradfield`, `Adjoint_Vecs.val`, `DII_ref` (input `Adjoint_DII_ref`)
  removed. `mdI < _MAX_OBS_` tests dropped: they were false with exactly `_MAX_OBS_` points,
  which then skipped the observation points.
- N-021 dead part: `compareVecs`, `genRandVec`, `VecSetBC`, `TestInterp`, `TestInterpBC` removed;
  `MGLevelCreate` is used, only its header declaration had 3 instead of 4 arguments (fixed).
- N-022 dead part: `AVDLoadPoints`, `AVDInjectDeletePoints`, `AVDExecuteMarkerInjection`,
  `_DISP_`, commented test block in subgrid.cpp removed.
- Tests: t02, t06-t08, t11, t15, t16, t19, t37 pass; t26 only dyndike_4core fails (known issue).
  No results change.

## Local iterations and temperature coupling (2026-10-23, committed ac3eda06, 235f3424, a2f65f5b)

- Non-converged local rheology iterations are now an error (`checkConvConstEq`), not a warning.
- Cause of the warnings in tests: `lmaxit = 8` (t16 PhaseTransNotInAirBox_move, t25, t26 x5,
  t27); removed with all other `lmaxit`/`lrtol` lines (t37 set the defaults). Defaults are used.
- t16 PhaseTransNotInAirBox_move ("fails to converge" in runtests.jl): Newton
  stalled at an absolute floor ~0.027. Temperature is updated in every nonlinear iteration
  (`SNESLineSearchPreCheck`) and conductivity / heat source jumped in T (useTk at T_Nu, dike at
  T_sol), so T kept flipping. With rtol 1e-4 it converges only with both features off.
- Fix: `RAMP` (tools.h) smooths the jumps over +-10 C. useTk applied per phase with its own
  T_Nu (was averaged with 0 for phases without T_Nu). S-023 fixed (dike kfac, assigned to Garrett
  Ito, taken over). With rtol 1e-4 steps 1-8 now converge in 8-13 iterations, but from step 13
  Newton still stalls near |F| ~ 1e-2, so this test keeps rtol 2.5e-2 and its loose check
  (comment in runtests.jl says why). Open: find the remaining floor (found, see next item).
- Remaining floor found (2026-10-23): the temperature matrix (`JacResGetTempMat`) had no
  derivatives of k(T) and of the dike heat source. With the steep ramps (nu_k = 5 and kfac over
  20 C) the T update alone fell into a 2-cycle (|ge| 0.0129 / 0.0142), which kept Newton at
  ~9e-3. The matrix now has the dk/dT terms (face k is a two-cell average) and d(rho_A)/dT on the
  diagonal (both by central differences, dT = 1e-3 C); it is no longer symmetric (fgmres).
  t16 PhaseTransNotInAirBox_move back to its original settings (rtol 1e-3, maxit 20) and the
  usual t16 checks; linear maxit 20 -> 50 (needed ~23). All 20 steps converge in 6-14
  iterations, final |eRes| ~1e-7 (was ~2e-3). Results change: references of t16, t26, t27
  updated. Tests: t12, t16, t26, t27, t28 pass. Committed a2f65f5b.
- Ramp half-width (10 C) and derivative perturbation (1e-3 C) are now inputs `T_ramp` and
  `T_pert` (Controls, JacRes.cpp; listed in info/options/input_file.dat). Defaults unchanged,
  no results change.
- Workarounds from f37db09c / N-024 / S-029 removed: t26 heat_rhoA and variableM back to
  nonlinear rtol 1e-3 (were 5e-2, 1e-2) and the usual checks; kfac atol back to 1e-4;
  dyndike maxit back to 50. t26 references updated; t16, t26 pass.
- Tests: t16, t25, t26, t27, t37 pass. Results change: references of t16, t26, t27 updated.
