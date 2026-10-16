# Task 11: Output

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] outFunct.h / outFunct.cpp: all functions
- [x] paraViewOutBin.h / paraViewOutBin.cpp: all functions
- [x] paraViewOutAVD.h / paraViewOutAVD.cpp: all functions
- [x] paraViewOutMark.cpp: all functions
- [x] paraViewOutPassiveTracers.cpp: all functions
- [x] paraViewOutSurf.cpp: all functions
- Checked helpers: InterpCenterCorner / edge / face (interpolate.cpp), JacResGetOverPressure,
  JacResGetSHmax/EHmax, JacResGetVelGrad (ghost updates), restart re-creation (LaMEMLib.cpp:276-279).

## Notes and open questions
- Output buffer is filled per component with stride `ncomp`; offsets in the XML header match
  what is dumped for all live vectors.
- Pressure shift: buffer = cf*p - shift, shift = -cf*pShift, so output is p + pShift. Consistent
  with constEq.cpp:133.
- J2 invariants: center part 0.5*(sxx²+syy²+szz²) plus edge parts s² are added at corners, then
  sqrt. Correct (average of squares, not square of average).
- Stress output adds 2*eta_st*d back to the stored history stress (pf = 0 for the initial guess).
- `jr->ge` lives on DA_T, scattered with DA_CEN. Same layout and stencil, so it works.
- PVTR/PVTS/AVD piece order uses the same rank -> (rx, ry, rz) mapping as fdstag.cpp:967;
  surface files use `dsz.color = rx + ry*Px`, which matches `iproc` in the PVTS file.
- AVD uses a uniform grid inside each rank's box (fine as a viewer), phase as UInt8 (max 32 phases).
- Mixed cells give non-integer `phase` (visID weighted by phase ratio). By design.

## Issues
Moved to findings.md as M-039, M-040, M-041, M-069, N-026.
