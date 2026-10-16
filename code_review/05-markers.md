# Task 5: Markers and advection

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] marker.h / marker.cpp: all functions
- [x] advect.h / advect.cpp: all functions
- [x] interpolate.h / interpolate.cpp: all functions
- [x] AVD.h / AVD.cpp: all functions
- [x] cvi.h / cvi.cpp: all functions
- [x] subgrid.h / subgrid.cpp: all functions
- Done by three review agents; the main claim of every major item was re-checked by hand.
  Items marked "agent" in the drafts were not re-checked in detail.

## Notes and open questions
- Checked and fine: hex 5-tet split; linear T formula; °C/K handling in SetTempFile/Save/InitFiles;
  exchange pairing with periodic x and 2 ranks; ADVCollectGarbage is fine for the (ascending)
  exchange and cvi lists; cell/edge weight normalisation; APS/ATS invariants; stress rotation
  (Rodrigues, R S R^T) after the increment; MINMOD correction terms; RK2 scheme and home-rank
  reset; subgrid bookkeeping; CheckMarkers reduction before error.
- Order per step: ADVSelectTimeStep -> FreeSurfAdvect -> ADVAdvect (exchange) -> passive tracers
  -> BCStretchGrid -> erosion -> ADVRemap (marker control, CrossFreeSurf, projection).
- Open: size of S-046 effect on viscoelastic runs; whether msetup=files restart from
  ADVMarkSave output fails CheckMarkers (M-026).

## Issues (draft)
Moved to findings.md as S-013, S-035, S-014, S-045, S-046, S-015, S-016, S-017, S-018, S-019, M-020, M-057, M-058, M-021, M-059, M-022, M-023, M-024, M-060, M-025, M-026, M-027, M-028, M-061, M-029, M-062, M-030, M-031, M-063, N-029, N-013, N-022.

### marker.cpp (checked by me unless noted)
- M1 [major, confirmed] marker.cpp:2100 oblique ridge line uses `y` not `y - y_ridgeFront`; divides by 0 if yF == yB.
- M2 [major, confirmed] marker.cpp:873 `rand_ampl` divided by maxPhaseID instead of chLen (inf with 1 phase). `cosine` checked against maxPhaseID (864). Layer noise uses rand() -> depends on CPU count.
- M3 [major, confirmed] marker.cpp:1738-1742,1832-1836 phase diagram: buffer full / unknown format only prints and returns -> rho = 0 silently. nT*nP never checked against _max_pd_sz_ (40100). maxT/maxP printed one step too wide.
- M4 [major, confirmed] tools.cpp:627 interpStretch: test `CtrlPoly[nCP-1] != nCP-1` should be `!= numPoly-1`; later slices get uninitialized stretch. marker.cpp:1429 PolyPos < 1 / ordering / unknown VolID not checked.
- M5 [major, confirmed] advect.cpp:219-231 ADVMarkPerturb (rand_noiseGP) runs after ADVMapMarkToCells, no remap; with rand_noise also on, markers cross cells -> stale cellnum, wrong projection, possible non-local abort.
- M6 [minor, confirmed, agent] polygon file: Polys[_max_polygons_] (10), Vol.num/len/Fcount not bounds-checked.
- M7 [minor, confirmed] marker.cpp:599-604 SetTempFile loops over grid count, not actx->nummark (msetup=files). Ix+1 can be out of range when xp == ex.
- M8 [minor, confirmed] marker.cpp:289 ADVMarkSave: strrchr NULL if no '/' -> segfault. Saved without APS.
- M9 [minor, confirmed, agent] negative phase IDs not rejected at parse/CheckMarkers (caught later in ADVCheckMarkPhases).
- M10 [minor, confirmed, agent] ridge v_spread falls back to bc->velin, can be 0 -> inf age.
- M11 [minor, confirmed, agent] hard-coded kappa = 1e-6 for halfspace T (902,1031,1091); wrong in non-dim mode.
- M12 [note] ellipsoid missing fb->ID=jj (955); unknown Temperature string ignored; comment sin vs -cos (862/1943); polygon noise amplitude differs; rand seed same per rank; DirMake rank-0-only SETERRQ (tools.cpp:293).

### advect.cpp / interpolate.cpp
- A1 [major, confirmed] advect.cpp:275,1076-1086 actx->periodic set for exy_num_periods, but bc.cpp:760 forbids periodic grid with Exy -> no periodic neighbours, grank = -1, outflow markers deleted after x-shift. Flag also forces BASIC_EULER for no gain.
- A2 [minor, confirmed, agent] advect.cpp:1079 `X[0] > ex` vs rank test `x >= ex` (fdstag.h:327): marker exactly at ex sent to wrap neighbour without shift.
- A3 [minor, confirmed] advect.cpp:1475-1529 ADVCheckCorners: n == 0 (cell and local neighbours empty) -> sind uninitialized, reads markers[sind] from zero-size buffer.
- A4 [minor, confirmed, agent] ADVCheckCorners: parent search only in local cells (CPU-count dependent), same seed each call/rank, xp == xc counted in no corner, inject count printed for rank 0 only.
- A5 [major, confirmed mechanism; size possible] LaMEMLib.cpp:699 + tssolve.cpp:308: with non-fixed dt, solve uses old dt, ADVSelectTimeStep then sets ts->dt = dt_next before ADVAdvect -> positions, APS/ATS, stress rotation, time use new dt, stress/T increments used old dt.
- A6 [note, confirmed, agent] subgrid.cpp:525-603 ADVMarkCrossFreeSurf changes phase only; stress/APS/T kept (has WARNING comment).
- A7 [minor, confirmed, agent] interpolate.cpp:102,158,220,377 corner output at periodic x seam uses one side only (output only).
- A8 [minor, confirmed] advect.cpp:2043 adiabatic gradient uses |z - Ztop|; rock above InitLevel is heated. Same at bc.cpp:2793.
- A9 [note, agent] SETERRQ on PETSC_COMM_SELF in collective code (advect.cpp:1896, tools.cpp:798); nbyte cast to PetscMPIInt (1125); nmark_lim unchecked; unused _DISP_; redundant FreeSurfGetAirPhaseRatio call; PetscMalloc per injected marker.
