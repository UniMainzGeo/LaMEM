# Fix plan for findings assigned to Anton Popov

Plan for the 83 findings assigned to Anton in `findings_audit_request_feedback.md`
(80 by authorship, 3 by area). Made 2026-10-20 at commit 61bb7aa0. Each finding is checked in
the code before it is fixed. Tests: `make test` in `test/` (most tests use the opt build); compile check with `make mode=deb all`. Fixes that change
results are reported; reference outputs are not updated silently.

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
| Shared | S-024 (Andrea), M-008 (Wenrong) | trivial; tell them |

## B. Fix next (more work, clear fix)

- S-009: "user velout given" flag; default stays mass-conserving.
- S-029: fill boundary ghosts of lithostatic and pore pressure.
- S-035, M-021, M-058, M-059: AVD marker control, one batch.
- S-039: read solver-default inputs with command-line options applied.
- M-016, M-030: reject typos and bad values.
- M-051, M-053, M-063: one helper so a one-rank error stops all ranks (no hang).
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
11. Nothing to do: M-019 (mitigated by `-snes_max_funcs unlimited`); N-028, N-029 (Boris: leave).
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

- S-006: Boris asked for a test (pure shear with non-zero `bg_ref_point`, check the global box).
- Tell Andrea about S-024 and Wenrong about M-008 (both fixed).
- `FindPointInCellAdjoint` (adjoint.cpp:102) has the same SETERRQ-in-PetscInt bug as M-022.
  Proposed: declare `FindPointInCell` in AVD.h, delete the copy, check `< 0` at the 8 callers.
  Boris's code: fix or add to his list (not decided).
- Next batch: AVD marker control from B (S-035, M-021, M-058, M-059). S-013 + S-035 together
  likely explain random phase loss/duplication under `mctrl = avd`.
- Tests: `make test` passed after A (115 pass, t37 FastScape broken as before). No
  reference output changed, so the fixed paths (S-006, S-007, S-013, S-024, M-008) are not covered.

## Order of work

One commit per group of A, then B, then C after discussion.
