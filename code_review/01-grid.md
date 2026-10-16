# Task 1: Grid and parallel layout

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] fdstag.h / fdstag.cpp: all functions
- [x] tools.h / tools.cpp: all functions
- [x] Tensor.h / Tensor.cpp: all functions
- Callers checked: advect.cpp:1993-1998 (CFL), bc.cpp:2727-2741 (grid stretch),
  advect.cpp:793-805 + JacRes.cpp:797 (stress rotation), JacResAux.cpp:56-110, 140-175
  (SHmax/EHmax), marker.cpp:1415-1440 (polygon stretch).

## Notes and open questions
- Task 2 question settled: Tensor2RS2DSpectral sorts eigenvalues in descending order, so v1
  belongs to the largest one. SHmax uses v2 (most compressive) and EHmax uses v1 (most
  extensional). Both are right.
- Stress rotation is fine: vorticity = curl v, angle = |w|*dt/2, Rodrigues matrix and
  R*S*R^T are written out correctly.
- Jacobi rotations (3D and 2D) match the Numerical Recipes form.
- Discret1DCoarsenCoord ghost exchange and Discret1DGatherCoord counts are right.
- Ghost coordinates at boundaries copy the first/last cell size (periodic: the opposite
  end's size), so COORD_NODE(-1) etc. always give real, shifted positions.
- FDSTAGSetEdgeCorner* fills edge/corner ghosts by linear extrapolation (a+b-c). For periodic
  x this also overwrites x-ghosts at y/z edges. Probably harmless, not checked further.
- Discret1DStretch only updates the fine grid. Whether multigrid coarse grids are
  re-coarsened after stretching is a task 4 question.
- getISA / Tensor2RNEigen / Tensor2RSSpectral have no callers (N-009).

## Issues
Moved to findings.md as S-006, S-007, M-012, N-009, N-010, M-051, N-018.
