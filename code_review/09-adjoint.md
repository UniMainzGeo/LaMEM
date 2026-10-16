# Task 9: Adjoint and inversion

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] adjoint.h / adjoint.cpp: all functions read. cellConstEqFD, edgeConstEqFD, devConstEqFD,
  setUpPhaseFD are never called (dead), so only skimmed.
- [x] objFunct.h / objFunct.cpp: all functions. ObjFunct is never used anywhere (dead).
- [x] Call sites: LaMEM.cpp:46-76, LaMEMLib.cpp:620, 672-683, 787-789.
- PSD derivative formulas (adjoint.cpp:3420-3450) not checked term by term.

## Notes and open questions
- Flow: LaMEMAdjointMain -> read input -> params copied to the PETSc options database ->
  LaMEMLibMain is run once per FD gradient plus once for the adjoint. The adjoint gradient is
  computed after every time step (LaMEMLib.cpp:683); only the last one is kept.
- The adjoint solve uses J, not J^T (documented in the file header). Only right for symmetric J.
- dr/dp is a finite difference of the residual: the material DB is rebuilt from the options
  database and the phases are swapped into jr->dbm. Only `phases[]` is swapped, so softening
  laws, phase transitions etc. can't be adjoint parameters.
- Test run (debug build, 1 core) on t08 AdjointGradients_CompareGradients.dat with
  `Adjoint_GradientCalculation = CostFunction`:
  - point (0.5,0.5,0.5): cost 8.6676e-4 = 1/2 (v - target)^2, FD and adjoint agree (0.5%).
  - point (0.31,0.43,0.37): cost 8.8026e-4, but 1/2 (v - target)^2 = 8.6841e-4.
    eta[1]: FD -4.42e-10, adjoint +4.23e-9 (wrong sign). eta[0]: FD 1.72e-5, adjoint 1.49e-5.
  - Files: scratchpad only (not kept). No test covers CostFunction mode with velocity points.
- To build the debug binary here: `make mode=deb surf=scape` (FASTSCAPE_LIB is set).

## Issues
Moved to findings.md as S-040, S-025, S-026, S-027, S-028, S-047, S-041, S-042, S-048, M-070, M-042, M-071, M-072, M-043, M-044, N-027, N-015.
