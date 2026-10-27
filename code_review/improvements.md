# Improvements

Ideas found during the review that improve the code but don't fix a bug. Do them after the fixes.

## Multigrid transfer of viscosity (from N-021, 2026-10-27)

Viscosity enters the coarse levels in two ways: Galerkin levels compute R*A*P with fixed
geometric P (viscosity jumps are smoothed over); matrix-free levels average fine viscosity
arithmetically over 8 cells (matData.cpp MatDataRestrictParam3D/2D, also edge values), which
favours the strong phase.

- a. Option for coarse viscosity averaging (arithmetic / geometric / harmonic, default
  arithmetic). Compare multigrid iteration counts on t02, t03 and a falling block with
  viscosity contrast 1e4-1e6.
- b. P weights from node coordinates instead of fixed 9/16, 3/16, ... (non-uniform grids).
  Expected small gain.
- Dropped: R = P^T (R and P stencils differ by design); viscosity-weighted P (much more
  expensive).
