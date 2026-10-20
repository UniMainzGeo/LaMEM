# Fix plan for findings assigned to Anton Popov

Plan for the 83 findings assigned to Anton in `findings_audit_request_feedback.md`
(80 by authorship, 3 by area). Made 2026-10-20 at commit 61bb7aa0. Each finding is checked in
the code before it is fixed. Tests run in debug mode (`make mode=deb test`). Fixes that change
results are reported; reference outputs are not updated silently.

## A. Fix right away (small, clear, local)

| Group | IDs | Fix |
|---|---|---|
| One-line logic | S-001, S-007, S-010, M-017, M-027, M-037, M-040 | `fixPhase == -1` check, `abs(v)`, `> 1`, `>=`, `\|\|`, move `free()` |
| Input data | S-003 | `Vn = 0` (Boris, GeoParams) |
| Grid | S-006 | stretch global box about `ref`; add a test |
| Inflow BC | S-011 | store base velocity in BCCtx, not a `static` |
| Markers | S-013, S-019, M-020, M-022, M-025, M-026, M-031 | sort `idel`, remap after perturbation, empty/NULL checks, free leak |
| Parallel | S-014, M-038 | global indices; total cell count |
| Polygons | S-018 + M-012 | same bug, fix together |
| Surface (Tobias) | S-004, M-005 | Boris: fix |
| Small | M-018, M-041, N-004, N-013 | missing SNES/KSP reasons, truncate .pvd, real creep-input checks, guard division |
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

## Order of work

One commit per group of A, then B, then C after discussion.
