# Task 7: Phase transitions and dikes

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] phase_transition.h / phase_transition.cpp: all functions
- [x] dike.h / dike.cpp: all functions
- [x] Call sites: LaMEMLib.cpp (init, restart read/write, destroy, Phase_Transition at step
  start), JacRes.cpp:1004 (GetDikeContr), JacResTemp.cpp:125 (Dike_k_heatsource)
- [x] Example input info/options/input_file.dat (phase transition and dike blocks), test/t26_Dike
- Not run: nothing was reproduced by running LaMEM; all findings are from reading the code.

## Notes and open questions
- Call order per step: Phase_Transition (-> Locate_Dike_Zones -> marker loop ->
  ADVInterpMarkToCell) runs at the start of each step, before BCApply and the solve.
  It uses `ts->time` and `ts->dt` of the step about to be done.
- `DBMat` lives in the zeroed `LaMEMLib` struct, so every field not read from the input is 0.
  This matters: Box types never read `PhaseAbove/Below` (stay 0), Constant/Clapeyron never read
  `PhaseInside/Outside` (stay 0). So `nphc` in Phase_Transition is always 1 (S-020).
- Box: `InAbove = 1` means *inside*, which maps to `PhaseInside` = PH1 = "below" phase.
  For Constant/Clapeyron `InAbove = 1` means PH2 = "above" phase. The `ResetParam` logic uses
  `InsideAbove` and so means different things per type (M-064).
- Dynamic diking needs cpu_x = 1 (checked only in Set_dike_zones) and assumes the same ny on
  all y-ranks (MPI messages of size nx*ny / ny+1 to neighbours). Smoothing only looks one
  rank up and down in y.
- On restart, DBMatCreate is called again after DynamicPhTr_ReadRestart, so phase-transition
  parameters come from the input file but the moving/relocated box bounds come from the
  restart file. Looks intended.
- t26 dynamic-dike test uses `istep_nave = 1`, so the time averaging (S-022) is not tested.
- `Overwrite_density` writes `mat[].rho` for both phases of a transition; if a phase is in
  several transitions the last one wins. Fine as long as it is documented.
- Linear box temperature profile `T = d*(top-bot)+top` with d in [-1,0] is correct.
- Zircon_Reidite has T0 = 25 C and slope 1.9 MPa/K: just copied from the paper, not checked.

## Issues
Moved to findings.md as S-020, S-021, S-036, S-037, M-064, M-065, N-014, S-022, S-023, M-032, M-033, M-073, M-034, M-035, N-023. Halfspace kappa = 1e-6 added to M-062.
