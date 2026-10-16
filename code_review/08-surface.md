# Task 8: Free surface and erosion

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] surf.h / surf.cpp: all functions
- [x] fastscape.h / fastscape.cpp: all functions (the FastScape Fortran library itself not checked)
- [x] Call order in LaMEMLib.cpp:690-760, restart read/write order
- Not checked: paraViewOutSurf.cpp (task 10), subgrid/advect use of `surf->phase` (task 5).

## Notes and open questions
- Call order per step: FreeSurfAdvect (or FastScapeCopyVelocity) -> ADVAdvect -> BCStretchGrid
  -> erosion / slope erosion / sedimentation / diffusion (or FastScapeRun + SmoothMaxAngle)
  -> ADVRemap -> GetAirPhaseRatio -> TSSolStepForward. Erosion uses the start-of-step time.
- Whole LaMEMLib struct is fwrite'd to the restart file, so margin positions (sed. model 2),
  `surf->phase` and `count_nsteps` survive a restart. Pointers are re-linked in LaMEMLib.cpp:470.
- `unit_m_yr` code is right; only its comment is inverted (N-007).
- FastScape velocity/time scale factors look odd (velocity_fs = 1/yr in SI, 100 in GEO) but
  are always used as divisors, so the conversion is right.
- GatherVariableFromLaMEM relies on PETSc rank order (x fastest, z last) to rebuild the
  surface on rank 0. Correct for the default DMDA layout.
- FreeSurfSmoothMaxAngle does one pass only, so the max angle is not guaranteed after it.
- Slope-dependent erosion lowers every sloping node, also below sea level, and isn't
  mass-conserving. That is a model choice, but worth saying in the docs.
- FreeSurfGetAirPhaseRatio leaves cells alone when markers say 100% air but the surface says
  otherwise (WARNING comment at surf.cpp:930). Known.
- PassValue2D's first loop uses a swapped index formula, but on both sides, so it's harmless.

## Issues
Moved to findings.md as S-032, S-033, S-004, M-005, S-005, S-034, S-043, M-006, M-050, M-007, M-008, M-009, M-010, M-011, N-007, N-008.
