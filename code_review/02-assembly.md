# Task 2: Equations and assembly

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] JacRes.h
- [x] JacRes.cpp: all functions
- [x] JacResAux.cpp: SHmax/EHmax (skimmed), OverPressure, LithoStaticPressure, PorePressure,
      GetPermea, GetVelGrad
- [x] JacResTemp.cpp: all functions
- Not checked: the matching Jacobian code in matrix.cpp / matFree.cpp (task 4), constEq
  internals (task 3), dike.cpp helpers (task 7).

## Notes and open questions
- `pLithoVisc` is **on by default** (JacRes.cpp:57), so creep with an activation volume uses
  lithostatic pressure unless the user turns this off. This makes S-029 relevant in
  default setups.
- `AdiabHeat` defaults to 0 and `FSSA` defaults to 1.
- `init_pres` is not used in any test or example, so S-001 is untested.
- t37_Collision_FastScape uses `Adiabatic_Heat = 1`. Running it would show whether S-030
  produces visible heating at the boundaries.
- SHmax/EHmax pick the eigenvector by its position (v1/v2). Whether that's right depends on
  the eigenvalue order in Tensor2RS2DSpectral. Check this in task 1.
- Output of lithostatic pressure is fine at boundaries: outFunct uses `use_bound = 0`, which
  clamps the indices.

## Issues
Moved to findings.md as S-001, S-029, S-030, M-001, M-045, M-046, N-028, M-002, N-001, M-047, N-002, N-016.
